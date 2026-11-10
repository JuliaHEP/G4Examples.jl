#---GPU optical photon ray tracer for the RICH detector-----------------------------------------------
# A port of the propagator in RICHTracer.jl to a KernelAbstractions.jl kernel, with one GPU thread
# per photon. The same kernel runs on any KernelAbstractions backend: `CPU()` (Julia threads),
# `MetalBackend()` (Apple GPUs), `CUDABackend()` (NVIDIA), `ROCBackend()` (AMD), ...
#
# The differences with the CPU version are:
# - all the structures are parametric in the floating point type `T`, because Apple GPUs do not
#   support Float64. With Float32 the geometrical tolerance must be larger (1 μm instead of 1 fm).
# - the energy dependent properties (refractive index of the gas, Rayleigh scattering length and
#   quantum efficiency) are evaluated on the host for each photon and copied with it, so the kernel
#   does not need lookup tables
# - the kernel returns a fixed size result for each photon (a status and a position), instead of
#   `nothing` or a position
# - the random numbers come from a counter-based generator (Philox2x32-10) whose key is a seed
#   for the batch and whose counter starts at the index of the photon. Each thread therefore has
#   its own independent stream of random numbers without any state in memory, and the results
#   do not depend on the backend or on the scheduling of the threads.

module RICHGPU

using KernelAbstractions
using GeometryBasics: Vec3

export GPUTracer, GPUPhoton, GPUHit, offload

#---Data structures----------------------------------------------------------------------------------
"""
    GPUPhoton{T}

The state of an optical photon at its creation, and its energy-dependent optical properties.
"""
struct GPUPhoton{T}
    position::Vec3{T}
    direction::Vec3{T}
    polarization::Vec3{T}
    nGas::T                     ## refractive index of the gas at the energy of the photon
    rayleighLength::T           ## Rayleigh scattering length in the aerogel
    qe::T                       ## quantum efficiency of the photodetector
end

"""
    GPUHit{T}

The result of the propagation of a photon: `detected` is true if the photon reached a
photodetector and passed the quantum efficiency test, and `position` is the hit position.
"""
struct GPUHit{T}
    detected::Bool
    position::Vec3{T}
end

"""
    GPUTracer{T}(tracer)

The constants of the geometry in the floating point type `T`, built from a `RICHTracer`.
"""
struct GPUTracer{T}
    vesselLo::Vec3{T}
    vesselHi::Vec3{T}
    aerogelLo::Vec3{T}
    aerogelHi::Vec3{T}
    holeRadius::T
    aerogelRIndex::T
    aerogelAbsLength::T
    gasAbsLength::T
    pipeRadius::T
    mirrorRadius::T
    mirrorReflectivity::T
    mirrorCentre::NTuple{2,Vec3{T}}
    mirrorLo::NTuple{2,Vec3{T}}
    mirrorHi::NTuple{2,Vec3{T}}
    detectorLo::NTuple{2,Vec3{T}}
    detectorHi::NTuple{2,Vec3{T}}
    tolerance::T
end

function GPUTracer{T}(tr; tolerance = T == Float32 ? 1e-3 : 1e-6) where T
    v(x) = Vec3{T}(x)
    GPUTracer{T}(v(tr.vesselLo), v(tr.vesselHi), v(tr.aerogelLo), v(tr.aerogelHi), tr.holeRadius,
                 tr.aerogelRIndex, tr.aerogelAbsLength, tr.gasAbsLength, tr.pipeRadius, tr.mirrorRadius,
                 tr.mirrorReflectivity, map(v, tr.mirrorCentre), map(v, tr.mirrorLo), map(v, tr.mirrorHi),
                 map(v, tr.detectorLo), map(v, tr.detectorHi), tolerance)
end

#---Counter-based random number generator (Philox2x32-10)--------------------------------------------
# J. K. Salmon, M. A. Moraes, R. O. Dror and D. E. Shaw, "Parallel random numbers: as easy as 1, 2, 3",
# Proceedings of SC '11, ACM (2011), https://doi.org/10.1145/2063384.2063405
# Same rounds and constants as the Philox2x32-10 generator of the Random123 library, but only the
# first of the two 32-bit output words is used.
struct Philox
    key::UInt32
    counter::UInt32
    stream::UInt32
end

@inline function philox2x32(key::UInt32, c0::UInt32, c1::UInt32)
    for _ in 1:10
        prod = UInt64(0xd256d193) * UInt64(c0)
        hi, lo = (prod >> 32) % UInt32, prod % UInt32
        c0, c1 = hi ⊻ key ⊻ c1, lo
        key += 0x9e3779b9
    end
    return c0
end

# A uniform random number in (0, 1] and the new state of the generator
@inline function uniform(::Type{T}, rng::Philox) where T
    x = philox2x32(rng.key, rng.counter, rng.stream)
    return T((x >> 8) + one(UInt32)) * T(5.9604644775390625e-8), Philox(rng.key, rng.counter + one(UInt32), rng.stream)
end

#---Vector helpers (written explicitly to keep the kernel simple)------------------------------------
@inline vdot(a, b) = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]
@inline vcross(a, b) = typeof(a)(a[2]*b[3] - a[3]*b[2], a[3]*b[1] - a[1]*b[3], a[1]*b[2] - a[2]*b[1])
@inline vnorm(a) = sqrt(vdot(a, a))
@inline vnormalize(a) = a / vnorm(a)
@inline inside(p, lo, hi) = lo[1] <= p[1] <= hi[1] && lo[2] <= p[2] <= hi[2] && lo[3] <= p[3] <= hi[3]
@inline axisVector(::Type{T}, i, s) where T = Vec3{T}(s*(i == 1), s*(i == 2), s*(i == 3))

#---Geometry (same algorithms as in RICHTracer.jl)---------------------------------------------------
@inline function slab(p::Vec3{T}, d, lo, hi) where T
    tn, tf, an, af = T(-Inf), T(Inf), 0, 0
    for i in 1:3
        if d[i] == 0
            (p[i] < lo[i] || p[i] > hi[i]) && return (T(Inf), T(-Inf), 0, 0)
        else
            t1, t2 = (lo[i] - p[i])/d[i], (hi[i] - p[i])/d[i]
            t1 > t2 && ((t1, t2) = (t2, t1))
            t1 > tn && ((tn, an) = (t1, i))
            t2 < tf && ((tf, af) = (t2, i))
        end
    end
    return (tn, tf, an, af)
end

@inline function boxEntry(p::Vec3{T}, d, lo, hi, tol) where T
    tn, tf, an, _ = slab(p, d, lo, hi)
    return (tn > tol && tn <= tf) ? (tn, an) : (T(Inf), 0)
end

@inline function cylinderRoots(p::Vec3{T}, d, r) where T
    a = d[1]^2 + d[2]^2
    a == 0 && return (T(Inf), T(Inf))
    b = p[1]*d[1] + p[2]*d[2]
    disc = b^2 - a*(p[1]^2 + p[2]^2 - r^2)
    disc < 0 && return (T(Inf), T(Inf))
    s = sqrt(disc)
    return ((-b - s)/a, (-b + s)/a)
end

@inline positive(t::T, tol) where T = t > tol ? t : T(Inf)

@inline function mirrorDistance(tr::GPUTracer{T}, p, d, C, lo, hi) where T
    oc = p - C
    b = vdot(d, oc)
    disc = b^2 - vdot(oc, oc) + tr.mirrorRadius^2
    disc < 0 && return T(Inf)
    s = sqrt(disc)
    t = -b - s
    t > tr.tolerance && inside(p + t*d, lo, hi) && return t
    t = -b + s
    t > tr.tolerance && inside(p + t*d, lo, hi) && return t
    return T(Inf)
end

#---Optical processes (same algorithms as in RICHTracer.jl)------------------------------------------
@inline function rotateUz(u::Vec3{T}, v) where T
    up = u[1]^2 + u[2]^2
    if up > 0
        up = sqrt(up)
        return Vec3{T}((u[1]*u[3]*v[1] - u[2]*v[2])/up + u[1]*v[3], (u[2]*u[3]*v[1] + u[1]*v[2])/up + u[2]*v[3], -up*v[1] + u[3]*v[3])
    else
        return u[3] > 0 ? v : Vec3{T}(-v[1], v[2], -v[3])
    end
end

@inline function rayleigh(d::Vec3{T}, pol, rng) where T
    while true
        cosθ, rng = uniform(T, rng)
        r, rng = uniform(T, rng)
        r < T(0.5) && (cosθ = -cosθ)
        u, rng = uniform(T, rng)
        sinθ, ϕ = sqrt(max(zero(T), 1 - cosθ^2)), T(2π)*u
        newd = rotateUz(d, Vec3{T}(sinθ*cos(ϕ), sinθ*sin(ϕ), cosθ))
        newpol = pol - vdot(newd, pol)*newd
        r, rng = uniform(T, rng)
        if vnorm(newpol) < T(1e-6)
            newpol = rotateUz(newd, Vec3{T}(cos(T(2π)*r), sin(T(2π)*r), zero(T)))
        else
            newpol = vnormalize(newpol)
            r < T(0.5) && (newpol = -newpol)
        end
        r, rng = uniform(T, rng)
        vdot(newpol, pol)^2 >= r && return (newd, newpol, rng)
    end
end

@inline function fresnel(d::Vec3{T}, pol, N, n1, n2, rng) where T
    cost1 = -vdot(d, N)
    sint1 = sqrt(max(zero(T), 1 - cost1^2))
    sint2 = sint1 * n1/n2
    reflected = d + 2cost1*N
    sint2 >= 1 && return (reflected, -pol + 2vdot(pol, N)*N, false, rng)     ## total internal reflection
    cost2 = sqrt(1 - sint2^2)
    if sint1 > 0
        A_trans = vnormalize(vcross(d, N))
        E1_perp = vdot(pol, A_trans)
        E1_parl = vnorm(pol - E1_perp*A_trans)
    else
        A_trans, E1_perp, E1_parl = pol, zero(T), one(T)
    end
    s1 = n1*cost1
    E2_perp = 2s1*E1_perp/(n1*cost1 + n2*cost2)
    E2_parl = 2s1*E1_parl/(n2*cost1 + n1*cost2)
    E2_total = E2_perp^2 + E2_parl^2
    r, rng = uniform(T, rng)
    if r > n2*cost2*E2_total/s1                                         ## reflection
        sint1 > 0 || return (reflected, n2 > n1 ? -pol : pol, false, rng)
        E2_parl = n2*E2_parl/n1 - E1_parl
        E2_perp = E2_perp - E1_perp
        A_paral = vnormalize(vcross(reflected, A_trans))
        return (reflected, (E2_parl*A_paral + E2_perp*A_trans)/sqrt(E2_parl^2 + E2_perp^2), false, rng)
    else                                                                ## refraction
        newd = vnormalize(d + (cost1 - cost2*n2/n1)*N)
        sint1 > 0 || return (newd, pol, true, rng)
        A_paral = vnormalize(vcross(newd, A_trans))
        return (newd, (E2_parl*A_paral + E2_perp*A_trans)/sqrt(E2_total), true, rng)
    end
end

#---Propagation of one photon (same logic as `propagate` in RICHTracer.jl)--------------------------
@inline function propagate(tr::GPUTracer{T}, photon::GPUPhoton{T}, rng::Philox) where T
    p, d, pol = photon.position, photon.direction, photon.polarization
    nAerogel, nGas, tol = tr.aerogelRIndex, photon.nGas, tr.tolerance
    lo, hi, rh = tr.aerogelLo, tr.aerogelHi, tr.holeRadius
    lost = GPUHit{T}(false, p)
    inAerogel = inside(p, lo, hi) && p[1]^2 + p[2]^2 > rh^2
    for _ in 1:1000
        if inAerogel
            _, tbox, axis, _ = slab(p, d, lo, hi)
            thole = positive(cylinderRoots(p, d, rh)[1], tol)
            u1, rng = uniform(T, rng)
            u2, rng = uniform(T, rng)
            tabs = -log(u1) * tr.aerogelAbsLength
            tray = -log(u2) * photon.rayleighLength
            t = min(tbox, thole, tabs, tray)
            p = p + t*d
            if t == tabs
                return lost
            elseif t == tray
                d, pol, rng = rayleigh(d, pol, rng)
            else
                t == tbox && axis == 3 && d[3] < 0 && return lost        ## front face: into the air
                N = t == tbox ? axisVector(T, axis, -sign(d[axis])) : vnormalize(Vec3{T}(p[1], p[2], zero(T)))
                d, pol, transmitted, rng = fresnel(d, pol, N, nAerogel, nGas, rng)
                inAerogel = !transmitted
            end
        else
            tvessel = slab(p, d, tr.vesselLo, tr.vesselHi)[2]
            tbox, axis = boxEntry(p, d, lo, hi, tol)
            if tbox < Inf
                q = p + tbox*d
                q[1]^2 + q[2]^2 < rh^2 && (tbox = T(Inf))                 ## into the hole
            end
            thole = positive(cylinderRoots(p, d, rh)[2], tol)
            thole < Inf && !inside(p + thole*d, lo, hi) && (thole = T(Inf))
            tpipe = positive(cylinderRoots(p, d, tr.pipeRadius)[1], tol)
            tmirror1 = mirrorDistance(tr, p, d, tr.mirrorCentre[1], tr.mirrorLo[1], tr.mirrorHi[1])
            tmirror2 = mirrorDistance(tr, p, d, tr.mirrorCentre[2], tr.mirrorLo[2], tr.mirrorHi[2])
            tdet1 = boxEntry(p, d, tr.detectorLo[1], tr.detectorHi[1], tol)[1]
            tdet2 = boxEntry(p, d, tr.detectorLo[2], tr.detectorHi[2], tol)[1]
            u, rng = uniform(T, rng)
            tabs = -log(u) * tr.gasAbsLength
            t = min(tvessel, tbox, thole, tpipe, tmirror1, tmirror2, tdet1, tdet2, tabs)
            p = p + t*d
            if t == tdet1 || t == tdet2
                u, rng = uniform(T, rng)
                return GPUHit{T}(u < photon.qe, p)
            elseif t == tmirror1 || t == tmirror2
                n = vnormalize(p - (t == tmirror1 ? tr.mirrorCentre[1] : tr.mirrorCentre[2]))
                u, rng = uniform(T, rng)
                (vdot(d, n) > 0 && u < tr.mirrorReflectivity) || return lost
                d, pol = d - 2vdot(d, n)*n, -pol + 2vdot(pol, n)*n
            elseif t == tbox || t == thole
                N = t == tbox ? axisVector(T, axis, -sign(d[axis])) : -vnormalize(Vec3{T}(p[1], p[2], zero(T)))
                d, pol, transmitted, rng = fresnel(d, pol, N, nGas, nAerogel, rng)
                inAerogel = transmitted
            else
                return lost
            end
        end
    end
    return lost
end

#---Kernel and host interface------------------------------------------------------------------------
@kernel function propagate_kernel!(hits, @Const(photons), tracer, key)
    i = @index(Global, Linear)
    # each photon has its own stream of random numbers, given by its index
    @inbounds hits[i] = propagate(tracer, photons[i], Philox(key, zero(UInt32), i % UInt32))
end

"""
    offload(backend, tracer::GPUTracer{T}, photons::Vector{GPUPhoton{T}}; seed=rand(UInt32), workgroupsize=256)

Copy the photons to the device of the `backend`, propagate them with one thread per photon, and
return a vector of `GPUHit{T}` on the host, in the same order as the photons.
"""
function offload(backend, tracer::GPUTracer{T}, photons::Vector{GPUPhoton{T}}; seed=rand(UInt32), workgroupsize=256) where T
    n = length(photons)
    n == 0 && return GPUHit{T}[]
    dphotons = KernelAbstractions.allocate(backend, GPUPhoton{T}, n)
    dhits = KernelAbstractions.allocate(backend, GPUHit{T}, n)
    copyto!(dphotons, photons)
    propagate_kernel!(backend, workgroupsize)(dhits, dphotons, tracer, seed % UInt32; ndrange=n)
    KernelAbstractions.synchronize(backend)
    return Array(dhits)
end

end # module RICHGPU
