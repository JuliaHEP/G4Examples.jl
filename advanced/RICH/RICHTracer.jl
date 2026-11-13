#---Optical photon propagator using the Geom4hep geometry-----------------------------------------------
# A stand-alone propagator of optical photons. It does not call Geant4: the geometry is exported by
# Geant4 to a GDML file and read by Geom4hep.jl, which provides the navigation (in which volume is a
# point, distance to the next boundary, normal of the surface). The propagator only needs, in addition,
# the optical properties of the materials and the volumes with special surfaces. It implements the
# same physics as the Geant4 optical processes:
# - absorption (G4OpAbsorption) and Rayleigh scattering (G4OpRayleigh) in the materials with the
#   properties ABSLENGTH and RAYLEIGH
# - Fresnel reflection and refraction, including the polarization, at the boundaries between two
#   materials with a refractive index RINDEX (G4OpBoundaryProcess, polished dielectric_dielectric)
# - specular reflection with a given reflectivity at the surface of the mirror volumes (polished
#   dielectric_metal skin surface)
# - detection when the photon enters a photodetector volume
# - absorption at the boundary with any other material without refractive index
# The material properties are tabulated on the energy grid used to build the Geant4 material
# properties tables, and interpolated linearly, as in Geant4.
#
# The photons are independent, so the propagation can be parallelized trivially. The same code runs
# on the CPU and in GPU kernels (see RICHTracerGPU.jl): the tracer is parametric in the floating point
# type, all its data is in arrays that can be moved to a GPU, and the random number generator is an
# argument of the propagation.

module RICHOptics

using Geom4hep: Geom4hep, locate, compute_step, in_world, material_index, read_gdml, close_geometry, shape_call
using LinearAlgebra: dot, cross, norm, normalize
using Random: AbstractRNG, default_rng
using GeometryBasics: Vec3
import Adapt

export OpticalPhoton, PhotonState, OpticalTracer, trace, propagate, uniform

const P3 = Geom4hep.Point3          # Geom4hep points and vectors (GeometryBasics has its own Point3)
const V3 = Geom4hep.Vector3

#---Photons----------------------------------------------------------------------------------------
"""
    OpticalPhoton

An optical photon as created by Geant4, collected by the stacking action before being killed.
"""
struct OpticalPhoton
    position::Vec3{Float64}
    direction::Vec3{Float64}
    polarization::Vec3{Float64}
    energy::Float64
    time::Float64
    primary::Bool               ## produced by the primary track
end

"""
    PhotonState{T}

What the propagator needs of a photon, in the floating point type `T` of the tracer.
"""
struct PhotonState{T<:AbstractFloat}
    position::P3{T}
    direction::V3{T}
    polarization::V3{T}
    energy::T
end
PhotonState{T}(p::OpticalPhoton) where T = PhotonState{T}(P3{T}(p.position...), V3{T}(normalize(p.direction)...),
                                                          V3{T}(p.polarization...), T(p.energy))

#---Tracer-----------------------------------------------------------------------------------------
const NOSURFACE, MIRROR, DETECTOR = Int8(0), Int8(1), Int8(2)
const RINDEX, ABSLENGTH, RAYLEIGH = 1, 2, 3
const PROPERTIES = ("RINDEX", "ABSLENGTH", "RAYLEIGH")
const BIG = 1e30                    # "infinite" length that can be interpolated (Inf * 0 is NaN)

"""
    OpticalTracer{T}

Geometry model, tables of optical properties and surfaces used to propagate the photons. Use
`Adapt.adapt(backend, tracer)` to move it to a GPU.
"""
struct OpticalTracer{T<:AbstractFloat, M, A<:AbstractArray{T,3}, S<:AbstractVector{Int8}, R<:AbstractVector{T}}
    model::M                    ## Geom4hep GeometryModel
    tables::A                   ## tables[i, property, material] on a uniform energy grid
    emin::T                     ## first energy of the grid
    invde::T                    ## 1 / energy step of the grid
    ne::Int32                   ## number of energies
    surface::S                  ## for each volume of the model: NOSURFACE, MIRROR or DETECTOR
    reflectivity::R             ## for each volume: reflectivity of its surface (mirrors)
end
Adapt.@adapt_structure OpticalTracer

"""
    OpticalTracer{T}(gdml; optics, energies, mirrors=Dict(), detectors=String[])

Build a tracer for the geometry in the GDML file `gdml`, in the floating point type `T`:
- `optics`: for each material name, the optical properties as pairs `"RINDEX" => f`, where `f(E)` is
  a function of the photon energy. The properties used are RINDEX, ABSLENGTH and RAYLEIGH.
- `energies`: the uniform grid (a range) on which the properties are tabulated
- `mirrors`: for each name of mirror volume, the reflectivity of its surface
- `detectors`: the names of the photodetector volumes
"""
function OpticalTracer{T}(gdml::AbstractString; optics, energies::AbstractRange, mirrors=Dict{String,Float64}(),
                          detectors=String[]) where T
    model, index = close_geometry(read_gdml(gdml, T); with_index=true)
    matnames = [m.name for m in index.materials]
    labels = [v.label for v in index.volumes]
    for name in keys(optics)
        name in matnames || error("Material $name not found in the geometry")
    end
    for name in Iterators.flatten((keys(mirrors), detectors))
        name in labels || error("Volume $name not found in the geometry")
    end
    tables = Array{T}(undef, length(energies), length(PROPERTIES), length(matnames))
    for (m, name) in enumerate(matnames), (k, prop) in enumerate(PROPERTIES)
        f = get(Dict(get(optics, name, ())), prop, nothing)
        tables[:, k, m] .= isnothing(f) ? (k == RINDEX ? 0 : BIG) : min.(f.(energies), BIG)
    end
    surface = Int8[l in keys(mirrors) ? MIRROR : l in detectors ? DETECTOR : NOSURFACE for l in labels]
    reflectivity = T[get(mirrors, l, 0) for l in labels]
    OpticalTracer(model, tables, T(first(energies)), T(1/step(energies)), Int32(length(energies)), surface, reflectivity)
end

floattype(::OpticalTracer{T}) where T = T

# Linear interpolation of the tables at the photon energy: the bin and the fraction within the bin
@inline function energybin(tr::OpticalTracer{T}, E::T) where T
    x = clamp((E - tr.emin) * tr.invde, zero(T), T(tr.ne - 1))
    i = min(unsafe_trunc(Int32, x), tr.ne - Int32(2))
    return i + Int32(1), x - T(i)
end
@inline property(tr::OpticalTracer, k, m, i, f) = @inbounds tr.tables[i, k, m] * (1 - f) + tr.tables[i + 1, k, m] * f

#---Random numbers---------------------------------------------------------------------------------
"""
    uniform(T, rng) -> (u, rng)

A uniform random number of type `T` and the new state of the generator. The propagation takes the
generator as an argument and returns its new state, so that it can use a stateless generator in
GPU kernels (see `Philox` in RICHTracerGPU.jl).
"""
@inline uniform(::Type{T}, rng::AbstractRNG) where T = (rand(rng, T), rng)

#---Optical processes------------------------------------------------------------------------------
# Rotate the vector v given in a frame where u is the z axis (CLHEP rotateUz)
@inline function rotateUz(u::V3{T}, v::V3{T}) where T
    up = u[1]^2 + u[2]^2
    if up > 0
        up = sqrt(up)
        return V3{T}((u[1]*u[3]*v[1] - u[2]*v[2])/up + u[1]*v[3], (u[2]*u[3]*v[1] + u[1]*v[2])/up + u[2]*v[3], -up*v[1] + u[3]*v[3])
    else
        return u[3] > 0 ? v : V3{T}(-v[1], v[2], -v[3])
    end
end

# Rayleigh scattering of a polarized photon (same algorithm as G4OpRayleigh). The loop of the
# rejection method is bounded, so that a kernel can never hang (the acceptance probability is ~1/2).
@inline function rayleigh(d::V3{T}, pol::V3{T}, rng) where T
    for _ in 1:1000
        cosθ, rng = uniform(T, rng)
        u, rng = uniform(T, rng)
        u < T(0.5) && (cosθ = -cosθ)
        u, rng = uniform(T, rng)
        sinθ, ϕ = sqrt(max(zero(T), 1 - cosθ^2)), T(2π)*u
        newd = rotateUz(d, V3{T}(sinθ*cos(ϕ), sinθ*sin(ϕ), cosθ))
        newpol = pol - dot(newd, pol)*newd
        if norm(newpol) < T(1e-6)                 ## new direction along the old polarization
            u, rng = uniform(T, rng)
            newpol = rotateUz(newd, V3{T}(cos(T(2π)*u), sin(T(2π)*u), zero(T)))
        else
            newpol = normalize(newpol)
            u, rng = uniform(T, rng)
            u < T(0.5) && (newpol = -newpol)
        end
        u, rng = uniform(T, rng)
        dot(newpol, pol)^2 >= u && return newd, newpol, rng
    end
    return d, pol, rng
end

# Fresnel reflection or refraction at a polished boundary between two dielectrics with indices
# n1 (incident) and n2, following G4OpBoundaryProcess::DielectricDielectric. N is the normal of
# the surface pointing back into the incident medium. Returns the new direction and polarization,
# whether the photon was transmitted, and the state of the generator.
@inline function fresnel(d::V3{T}, pol::V3{T}, N::V3{T}, n1::T, n2::T, rng) where T
    cost1 = -dot(d, N)
    sint1 = sqrt(max(zero(T), 1 - cost1^2))
    sint2 = sint1 * n1/n2
    reflected = d + 2cost1*N
    sint2 >= 1 && return reflected, -pol + 2dot(pol, N)*N, false, rng     ## total internal reflection
    cost2 = sqrt(1 - sint2^2)
    # Geant4 tests sint1 > 0. In finite precision sint1 of a normal incidence is of the order of
    # √eps (from 1 - cost1²), and the plane of incidence is then undefined (cross(d, N) ≈ 0).
    oblique = sint1 > 4sqrt(eps(T))
    oblique || (sint1 = zero(T))
    if oblique
        A_trans = normalize(cross(d, N))
        E1_perp = dot(pol, A_trans)
        E1_parl = norm(pol - E1_perp*A_trans)
    else                                                               ## normal incidence
        A_trans, E1_perp, E1_parl = pol, zero(T), one(T)
    end
    s1 = n1*cost1
    E2_perp = 2s1*E1_perp/(n1*cost1 + n2*cost2)
    E2_parl = 2s1*E1_parl/(n2*cost1 + n1*cost2)
    E2_total = E2_perp^2 + E2_parl^2
    u, rng = uniform(T, rng)
    if u > n2*cost2*E2_total/s1                                        ## reflection
        sint1 > 0 || return reflected, (n2 > n1 ? -pol : pol), false, rng
        E2_parl = n2*E2_parl/n1 - E1_parl
        E2_perp = E2_perp - E1_perp
        A_paral = normalize(cross(reflected, A_trans))
        return reflected, (E2_parl*A_paral + E2_perp*A_trans)/sqrt(E2_parl^2 + E2_perp^2), false, rng
    else                                                               ## refraction
        newd = normalize(d + (cost1 - cost2*n2/n1)*N)
        sint1 > 0 || return newd, pol, true, rng
        A_paral = normalize(cross(newd, A_trans))
        return newd, (E2_parl*A_paral + E2_perp*A_trans)/sqrt(E2_total), true, rng
    end
end

#---Geometry---------------------------------------------------------------------------------------
# Outward normal (global frame) of the volume of the state s at the global point p, zero if p is
# not on its surface
@inline function volumenormal(model, s, p::P3{T}) where T
    in_world(s) || return zero(V3{T})
    shape = @inbounds(model.volumes[s.volume]).shape
    return shape_call(Geom4hep.normal, model, shape, s.tolocal * p) * s.tolocal
end

# Normal of the boundary crossed at p, going from the state s to ns, pointing back into the volume
# of s. The surface is the one of the volume entered (if deeper) or left. If p is not within the
# tolerance of the surface (the navigator can push the point by a minimum step), the previous
# position of the photon `prev` is also tried.
@inline function boundarynormal(model, s, ns, p::P3{T}, prev::P3{T}, d::V3{T}) where T
    a, b = ns.depth > s.depth ? (ns, s) : (s, ns)
    n = volumenormal(model, a, p)
    iszero(n) && (n = volumenormal(model, b, p))
    iszero(n) && (n = volumenormal(model, a, prev))
    iszero(n) && (n = volumenormal(model, b, prev))
    iszero(n) && return -d                    ## not found: normal incidence
    return dot(n, d) > 0 ? -n : n
end
@inline samestate(a, b) = a.volume == b.volume && a.depth == b.depth && a.path == b.path

#---Propagation of one photon----------------------------------------------------------------------
"""
    trace(tr::OpticalTracer{T}, photon::PhotonState{T}, rng; maxsteps=1000) -> (detected, position, rng)

Propagate a photon until it is absorbed, leaves the detector or reaches a photodetector. Returns
whether it was detected, its last position (on the photodetector if detected) and the state of the
random number generator. The quantum efficiency is not applied.
"""
function trace(tr::OpticalTracer{T}, photon::PhotonState{T}, rng; maxsteps::Int32=Int32(1000)) where T
    model = tr.model
    p, d, pol = photon.position, photon.direction, photon.polarization
    i, f = energybin(tr, photon.energy)
    s = locate(model, p)
    for _ in Int32(1):maxsteps
        in_world(s) || break
        m = material_index(model, s)
        n1 = property(tr, RINDEX, m, i, f)
        n1 > 0 || break                          ## in a material without refractive index
        ##---distances to absorption and Rayleigh scattering, then to the next boundary
        u1, rng = uniform(T, rng)
        u2, rng = uniform(T, rng)
        tabs = -log(u1) * property(tr, ABSLENGTH, m, i, f)
        tray = -log(u2) * property(tr, RAYLEIGH, m, i, f)
        limit = min(tabs, tray)
        step, ns = compute_step(model, s, p, d, limit)
        prev, p = p, p + step * d
        if samestate(s, ns)                      ## no boundary crossed
            step < limit && continue             ## (the navigator only made a small move)
            tabs <= tray && break                ## absorbed
            d, pol, rng = rayleigh(d, pol, rng)
            continue
        end
        ##---at a boundary
        in_world(ns) || break                    ## leaving the world
        kind = @inbounds tr.surface[ns.volume]
        kind == DETECTOR && return true, p, rng
        N = boundarynormal(model, s, ns, p, prev, d)
        if kind == MIRROR
            u, rng = uniform(T, rng)
            u < @inbounds(tr.reflectivity[ns.volume]) || break
            d, pol = d - 2dot(d, N)*N, -pol + 2dot(pol, N)*N
        else
            n2 = property(tr, RINDEX, material_index(model, ns), i, f)
            n2 > 0 || break                      ## into a material without refractive index
            d, pol, transmitted, rng = fresnel(d, pol, N, n1, n2, rng)
            transmitted && (s = ns)
        end
        isfinite(d[1] + d[2] + d[3] + pol[1] + pol[2] + pol[3]) || break    ## numerical failure: lost
    end
    return false, p, rng
end

"""
    propagate(tr::OpticalTracer, photon::OpticalPhoton; rng=default_rng())

Propagate a photon on the CPU. Returns its position on the photodetector, or `nothing`.
"""
function propagate(tr::OpticalTracer{T}, photon::OpticalPhoton; rng=default_rng()) where T
    detected, p, _ = trace(tr, PhotonState{T}(photon), rng)
    return detected ? Vec3{Float64}(p...) : nothing
end

end # module RICHOptics

using .RICHOptics

#---The RICH detector-------------------------------------------------------------------------------
"""
    exportGDML(file=tempname() * ".gdml")

Write the geometry of the initialized Geant4 application to a GDML file and return its name.
"""
function exportGDML(file=tempname() * ".gdml")
    isfile(file) && rm(file)                 ## Geant4 does not overwrite an existing file
    Write(G4GDMLParser(), file, CxxPtr(GetWorldVolume()), false)
    return file
end

"""
    RICHTracer(det::RICHDetector, gdml; T=Float64)

Optical tracer for the RICH detector, with the geometry in the GDML file `gdml` (see `exportGDML`),
the optical properties of `richOptics(det)`, the two mirrors and the photodetectors.
"""
function RICHTracer(det::RICHDetector, gdml::AbstractString; T::Type=Float64)
    OpticalTracer{T}(gdml; optics = richOptics(det),
                     energies = range(first(photonEnergies), last(photonEnergies), length=length(photonEnergies)),
                     mirrors = Dict("MirrorTop" => det.mirrorReflectivity, "MirrorBottom" => det.mirrorReflectivity),
                     detectors = ["Photodetector"])
end
