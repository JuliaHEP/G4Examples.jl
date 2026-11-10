#---Optical photon ray tracer for the RICH detector----------------------------------------------------
# A stand-alone propagator of optical photons for the geometry defined in RICHDetector.jl. It does
# not call Geant4: it only needs the detector parameters and the optical properties (refractive
# indices, absorption and Rayleigh lengths, mirror reflectivity and quantum efficiency). It
# implements, for this particular geometry, the same physics as the Geant4 optical processes:
# - absorption in the aerogel and in the gas (G4OpAbsorption)
# - Rayleigh scattering in the aerogel (G4OpRayleigh)
# - Fresnel reflection and refraction, including the polarization, at the boundaries between the
#   aerogel and the gas (G4OpBoundaryProcess, polished dielectric_dielectric)
# - specular reflection on the spherical mirrors with a given reflectivity (polished dielectric_metal)
# - absorption at the boundary with any material without refractive index (air outside the vessel,
#   beam pipe, back of the mirrors and photodetectors)
# The photons are independent of each other, so the propagation can be parallelized trivially.
# The data structures are plain bits types, so they can also be copied to a GPU.

using LinearAlgebra: dot, normalize, cross, norm
using GeometryBasics: Vec3, Point3

const V3 = Vec3{Float64}

"""
    OpticalPhoton

An optical photon as created by Geant4, collected by the stacking action before being killed.
"""
struct OpticalPhoton
    position::V3
    direction::V3
    polarization::V3
    energy::Float64
    time::Float64
    primary::Bool               ## produced by the primary track
end

"""
    RICHTracer(det::RICHDetector)

Pre-computed geometry and optical constants used to propagate the photons.
"""
struct RICHTracer
    vesselLo::V3                ## gas vessel
    vesselHi::V3
    aerogelLo::V3               ## aerogel tile...
    aerogelHi::V3
    holeRadius::Float64         ## ...with a hole for the beam pipe
    aerogelRIndex::Float64
    aerogelAbsLength::Float64
    gasAbsLength::Float64
    pipeRadius::Float64         ## outer radius of the beam pipe
    mirrorRadius::Float64
    mirrorReflectivity::Float64
    mirrorCentre::NTuple{2,V3}  ## top and bottom mirrors
    mirrorLo::NTuple{2,V3}      ## box that cuts each spherical mirror
    mirrorHi::NTuple{2,V3}
    detectorLo::NTuple{2,V3}    ## top and bottom photodetector planes
    detectorHi::NTuple{2,V3}
    zRadiator::Float64          ## photons produced before this z come from the aerogel
end

function RICHTracer(det::RICHDetector; aerogelAbsLength=1m, gasAbsLength=100m)
    (; vesselHalfX, vesselHalfY, vesselZ, beamPipeRadius, beamPipeThickness, aerogelZ, aerogelThickness,
       aerogelHalfSize, mirrorRadius, detectorZ, detectorY, detectorHalfX, detectorThickness) = det
    yGap = beamPipeRadius + beamPipeThickness + 5mm      ## same as in richConstruct
    zmin = (vesselZ[1] + vesselZ[2])/2
    mlo(side) = V3(-(vesselHalfX - 1mm), side > 0 ? yGap : -(vesselHalfY - 1mm), zmin)
    mhi(side) = V3(vesselHalfX - 1mm, side > 0 ? vesselHalfY - 1mm : -yGap, vesselZ[2] - 1mm)
    dlo(side) = V3(-detectorHalfX, side > 0 ? detectorY[1] : -detectorY[2], detectorZ - detectorThickness)
    dhi(side) = V3(detectorHalfX, side > 0 ? detectorY[2] : -detectorY[1], detectorZ)
    RICHTracer(V3(-vesselHalfX, -vesselHalfY, vesselZ[1]), V3(vesselHalfX, vesselHalfY, vesselZ[2]),
               V3(-aerogelHalfSize, -aerogelHalfSize, aerogelZ), V3(aerogelHalfSize, aerogelHalfSize, aerogelZ + aerogelThickness),
               beamPipeRadius + beamPipeThickness + 1mm, det.aerogelRIndex, aerogelAbsLength, gasAbsLength,
               beamPipeRadius + beamPipeThickness, mirrorRadius, det.mirrorReflectivity,
               (V3(mirrorCentre(det, 1)...), V3(mirrorCentre(det, -1)...)),
               (mlo(1), mlo(-1)), (mhi(1), mhi(-1)), (dlo(1), dlo(-1)), (dhi(1), dhi(-1)),
               aerogelZ + aerogelThickness + 1mm)
end

#---Geometry helpers---------------------------------------------------------------------------------
const tolerance = 1e-6mm
axisVector(i, s) = V3(s*(i == 1), s*(i == 2), s*(i == 3))
inside(p, lo, hi) = all(lo .<= p .<= hi)

# Slab method: distances along the ray to enter and leave a box, and the axes of the faces
function slab(p, d, lo, hi)
    tn, tf, an, af = -Inf, Inf, 0, 0
    for i in 1:3
        if d[i] == 0
            (p[i] < lo[i] || p[i] > hi[i]) && return (Inf, -Inf, 0, 0)
        else
            t1, t2 = (lo[i] - p[i])/d[i], (hi[i] - p[i])/d[i]
            t1 > t2 && ((t1, t2) = (t2, t1))
            t1 > tn && ((tn, an) = (t1, i))
            t2 < tf && ((tf, af) = (t2, i))
        end
    end
    return (tn, tf, an, af)
end

# Distance to enter a box from outside (Inf if the box is missed) and the axis of the face
function boxEntry(p, d, lo, hi)
    tn, tf, an, _ = slab(p, d, lo, hi)
    return (tn > tolerance && tn <= tf) ? (tn, an) : (Inf, 0)
end

# Distances along the ray to the cylinder of radius r around the z axis, when the radius is
# decreasing (entering) and increasing (leaving)
function cylinderRoots(p, d, r)
    a = d[1]^2 + d[2]^2
    a == 0 && return (Inf, Inf)
    b = p[1]*d[1] + p[2]*d[2]
    disc = b^2 - a*(p[1]^2 + p[2]^2 - r^2)
    disc < 0 && return (Inf, Inf)
    s = sqrt(disc)
    return ((-b - s)/a, (-b + s)/a)
end
positive(t) = t > tolerance ? t : Inf

# Distance to the spherical mirror of the given side, Inf if the ray does not hit it
function mirrorDistance(tr::RICHTracer, p, d, side)
    oc = p - tr.mirrorCentre[side]
    b = dot(d, oc)
    disc = b^2 - dot(oc, oc) + tr.mirrorRadius^2
    disc < 0 && return Inf
    s = sqrt(disc)
    for t in (-b - s, -b + s)
        t > tolerance && inside(p + t*d, tr.mirrorLo[side], tr.mirrorHi[side]) && return t
    end
    return Inf
end

#---Optical processes--------------------------------------------------------------------------------
# Rotate the vector v given in a frame where u is the z axis (CLHEP rotateUz)
function rotateUz(u, v)
    up = u[1]^2 + u[2]^2
    if up > 0
        up = sqrt(up)
        return V3((u[1]*u[3]*v[1] - u[2]*v[2])/up + u[1]*v[3], (u[2]*u[3]*v[1] + u[1]*v[2])/up + u[2]*v[3], -up*v[1] + u[3]*v[3])
    else
        return u[3] > 0 ? v : V3(-v[1], v[2], -v[3])
    end
end

# Rayleigh scattering of a polarized photon (same algorithm as G4OpRayleigh)
function rayleigh(d, pol)
    while true
        cosθ = rand()
        rand() < 0.5 && (cosθ = -cosθ)
        sinθ, ϕ = sqrt(1 - cosθ^2), 2π*rand()
        newd = rotateUz(d, V3(sinθ*cos(ϕ), sinθ*sin(ϕ), cosθ))
        newpol = pol - dot(newd, pol)*newd
        if norm(newpol) < 1e-12           ## new direction along the old polarization
            ϕ = 2π*rand()
            newpol = rotateUz(newd, V3(cos(ϕ), sin(ϕ), 0.))
        else
            newpol = normalize(newpol)
            rand() < 0.5 && (newpol = -newpol)
        end
        dot(newpol, pol)^2 >= rand() && return (newd, newpol)
    end
end

# Fresnel reflection or refraction at a polished boundary between two dielectrics with indices
# n1 (incident) and n2, following G4OpBoundaryProcess::DielectricDielectric. N is the normal of
# the surface pointing back into the incident medium. Returns the new direction and polarization,
# and whether the photon was transmitted.
function fresnel(d, pol, N, n1, n2)
    cost1 = -dot(d, N)
    sint1 = sqrt(max(0., 1 - cost1^2))
    sint2 = sint1 * n1/n2
    reflected = d + 2cost1*N
    sint2 >= 1 && return (reflected, -pol + 2dot(pol, N)*N, false)      ## total internal reflection
    cost2 = sqrt(1 - sint2^2)
    if sint1 > 0
        A_trans = normalize(cross(d, N))
        E1_perp = dot(pol, A_trans)
        E1_parl = norm(pol - E1_perp*A_trans)
    else                                                               ## normal incidence
        A_trans, E1_perp, E1_parl = pol, 0., 1.
    end
    s1 = n1*cost1
    E2_perp = 2s1*E1_perp/(n1*cost1 + n2*cost2)
    E2_parl = 2s1*E1_parl/(n2*cost1 + n1*cost2)
    E2_total = E2_perp^2 + E2_parl^2
    if rand() > n2*cost2*E2_total/s1                                   ## reflection
        sint1 > 0 || return (reflected, n2 > n1 ? -pol : pol, false)
        E2_parl = n2*E2_parl/n1 - E1_parl
        E2_perp = E2_perp - E1_perp
        A_paral = normalize(cross(reflected, A_trans))
        return (reflected, (E2_parl*A_paral + E2_perp*A_trans)/sqrt(E2_parl^2 + E2_perp^2), false)
    else                                                               ## refraction
        newd = normalize(d + (cost1 - cost2*n2/n1)*N)
        sint1 > 0 || return (newd, pol, true)
        A_paral = normalize(cross(newd, A_trans))
        return (newd, (E2_parl*A_paral + E2_perp*A_trans)/sqrt(E2_total), true)
    end
end

#---Propagation of one photon------------------------------------------------------------------------
"""
    propagate(tr::RICHTracer, photon::OpticalPhoton)

Propagate a photon until it is absorbed, leaves the detector or reaches a photodetector. Returns
the position of the photon on the photodetector, or `nothing`. The quantum efficiency is not applied.
"""
function propagate(tr::RICHTracer, photon::OpticalPhoton; maxsteps=1000)
    p, d, pol, E = photon.position, photon.direction, photon.polarization, photon.energy
    nAerogel, nGas = tr.aerogelRIndex, c4f10RIndex(E)
    lo, hi, rh = tr.aerogelLo, tr.aerogelHi, tr.holeRadius
    inAerogel = inside(p, lo, hi) && p[1]^2 + p[2]^2 > rh^2
    for _ in 1:maxsteps
        if inAerogel
            ##---distances to the boundary of the aerogel, to absorption and to Rayleigh scattering
            _, tbox, axis, _ = slab(p, d, lo, hi)
            thole = positive(cylinderRoots(p, d, rh)[1])
            tabs = -log(rand()) * tr.aerogelAbsLength
            tray = -log(rand()) * aerogelRayleigh(E)
            t = min(tbox, thole, tabs, tray)
            p = p + t*d
            if t == tabs
                return nothing
            elseif t == tray
                d, pol = rayleigh(d, pol)
            else
                t == tbox && axis == 3 && d[3] < 0 && return nothing     ## front face: into the air
                N = t == tbox ? axisVector(axis, -sign(d[axis])) : normalize(V3(p[1], p[2], 0.))
                d, pol, transmitted = fresnel(d, pol, N, nAerogel, nGas)
                inAerogel = !transmitted
            end
        else
            ##---distances to the vessel walls, the aerogel, the beam pipe, the mirrors and the photodetectors
            tvessel = slab(p, d, tr.vesselLo, tr.vesselHi)[2]
            tbox, axis = boxEntry(p, d, lo, hi)
            tbox < Inf && (q = p + tbox*d; q[1]^2 + q[2]^2 < rh^2) && (tbox = Inf)   ## into the hole
            thole = positive(cylinderRoots(p, d, rh)[2])
            thole < Inf && !inside(p + thole*d, lo, hi) && (thole = Inf)
            tpipe = positive(cylinderRoots(p, d, tr.pipeRadius)[1])
            tmirror1, tmirror2 = mirrorDistance(tr, p, d, 1), mirrorDistance(tr, p, d, 2)
            tdet1, tdet2 = boxEntry(p, d, tr.detectorLo[1], tr.detectorHi[1])[1], boxEntry(p, d, tr.detectorLo[2], tr.detectorHi[2])[1]
            tabs = -log(rand()) * tr.gasAbsLength
            t = min(tvessel, tbox, thole, tpipe, tmirror1, tmirror2, tdet1, tdet2, tabs)
            p = p + t*d
            if t == tdet1 || t == tdet2
                return p
            elseif t == tmirror1 || t == tmirror2
                n = normalize(p - tr.mirrorCentre[t == tmirror1 ? 1 : 2])
                (dot(d, n) > 0 && rand() < tr.mirrorReflectivity) || return nothing
                d, pol = d - 2dot(d, n)*n, -pol + 2dot(pol, n)*n
            elseif t == tbox || t == thole
                N = t == tbox ? axisVector(axis, -sign(d[axis])) : -normalize(V3(p[1], p[2], 0.))
                d, pol, transmitted = fresnel(d, pol, N, nGas, nAerogel)
                inAerogel = transmitted
            else                                                        ## absorbed, beam pipe or outside
                return nothing
            end
        end
    end
    return nothing
end
