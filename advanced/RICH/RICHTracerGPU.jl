#---Offloading the optical photons to a GPU------------------------------------------------------------
# Runs the propagator of RICHTracer.jl (`trace`) in a KernelAbstractions.jl kernel, with one GPU thread
# per photon. The same kernel runs on any KernelAbstractions GPU backend: `MetalBackend()` (Apple GPUs),
# `CUDABackend()` (NVIDIA), `ROCBackend()` (AMD), ... With the `CPU()` backend the photons are
# propagated with Julia threads instead (as Geom4hep does for its own kernels: the CPU backend of
# KernelAbstractions ≥ 0.10 inlines all the functions, which makes the navigation kernels very large).
#
# What is specific to the GPU:
# - the tracer and the geometry model are copied to the device with `Adapt.adapt(backend, tracer)`.
#   Apple GPUs do not support Float64, so the tracer is then built in Float32 (`RICHTracer(det, gdml;
#   T=Float32)`): Geom4hep uses a larger geometrical tolerance in single precision (1 μm).
# - the kernel returns a fixed size result for each photon (a status and a position), instead of
#   `nothing` or a position
# - the random numbers come from a counter-based generator (Philox2x32-10) whose key is a seed
#   for the batch and whose counter starts at the index of the photon. Each thread therefore has
#   its own independent stream of random numbers without any state in memory, and the results
#   do not depend on the backend or on the scheduling of the threads.

module RICHGPU

using KernelAbstractions
using Adapt: adapt
using ..RICHOptics: RICHOptics, OpticalTracer, PhotonState, trace
import ..RICHOptics: uniform
const P3 = RICHOptics.P3

export Philox, GPUHit, offload, todevice

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

#---Kernel and host interface------------------------------------------------------------------------
"""
    GPUHit{T}

Result of the propagation of a photon: whether it was detected, and where.
"""
struct GPUHit{T<:AbstractFloat}
    detected::Bool
    position::P3{T}
end

@inline function tracephoton(tracer, photon::PhotonState{T}, key, i) where T
    # each photon has its own stream of random numbers, given by its index
    detected, p, _ = trace(tracer, photon, Philox(key, zero(UInt32), i % UInt32))
    return GPUHit{T}(detected, p)
end

@kernel function propagate_kernel!(hits, @Const(photons), tracer, key, offset)
    i = @index(Global, Linear) + offset
    @inbounds hits[i] = tracephoton(tracer, photons[i], key, i)
end

"""
    todevice(backend, tracer::OpticalTracer)

Copy the tracer (geometry model and tables) to the device of the `backend`.
"""
todevice(backend, tracer::OpticalTracer) = adapt(backend, tracer)
todevice(::CPU, tracer::OpticalTracer) = tracer

"""
    offload(backend, tracer, photons::Vector{PhotonState{T}}; seed=rand(UInt32), workgroupsize=256, maxlaunch=200_000)

Copy the photons to the device of the `backend`, propagate them with one thread per photon, and
return a vector of `GPUHit{T}` on the host, in the same order as the photons. The `tracer` must be
on the device (see `todevice`).

Large batches are propagated in several kernel launches of at most `maxlaunch` photons: the
operating system can abort (or hang) a long-running kernel on a GPU that also drives a display
(macOS error "Impacting Interactivity"). The random numbers of a photon do not depend on how the
batch is split.
"""
function offload(backend, tracer::OpticalTracer{T}, photons::Vector{PhotonState{T}}; seed=rand(UInt32), workgroupsize=256,
                 maxlaunch=200_000) where T
    n = length(photons)
    n == 0 && return GPUHit{T}[]
    dphotons = KernelAbstractions.allocate(backend, PhotonState{T}, n)
    dhits = KernelAbstractions.allocate(backend, GPUHit{T}, n)
    copyto!(dphotons, photons)
    kernel! = propagate_kernel!(backend, workgroupsize)
    for offset in 0:maxlaunch:n-1
        kernel!(dhits, dphotons, tracer, seed % UInt32, offset; ndrange=min(maxlaunch, n - offset))
        KernelAbstractions.synchronize(backend)
    end
    return Array(dhits)
end

function offload(::CPU, tracer::OpticalTracer{T}, photons::Vector{PhotonState{T}}; seed=rand(UInt32), kwargs...) where T
    hits = Vector{GPUHit{T}}(undef, length(photons))
    key = seed % UInt32
    Threads.@threads for i in eachindex(photons)
        @inbounds hits[i] = tracephoton(tracer, photons[i], key, i)
    end
    return hits
end

end # module RICHGPU
