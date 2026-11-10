using Geant4
using Geant4.SystemOfUnits
using Geant4.PhysicalConstants
using FHist, GeometryBasics, Printf
using LinearAlgebra: dot, normalize          ## Geant4 also exports `dot`
using KernelAbstractions
if Sys.isapple()
    using Metal
else
    using CUDA
end
using GLMakie

backend, T = if @isdefined(CUDA) && CUDA.functional()
    CUDABackend(), Float32
elseif @isdefined(Metal) && Metal.functional()
    MetalBackend(), Float32
else
    CPU(), Float64
end
@printf "Backend: %s, precision: %s, Julia threads: %d\n" backend T Threads.nthreads()

include(joinpath(@__DIR__, "RICHDetector.jl"))
rich = RICHDetector()

struct RICHPhysics <: G4VUserPhysicsList
    function RICHPhysics(verbose)
        pl = FTFP_BERT(verbose)
        RegisterPhysics(pl, move!(G4OpticalPhysics(verbose)))
        return pl
    end
end

mutable struct RICHGunData <: G4JLGeneratorData
    gun::Union{Nothing, CxxPtr{G4ParticleGun}}
    particles::Vector{String}
    pmin::Float64
    pmax::Float64
    θmin::Float64
    θmax::Float64
    ##---the current track
    particle::String
    mass::Float64
    momentum::Float64
    direction::Vec3{Float64}
end

function RICHGunGenerator(; particles=["pi+", "kaon+", "proton"], pmin=1GeV, pmax=80GeV, θmin=50mrad, θmax=250mrad)
    data = RICHGunData(nothing, particles, pmin, pmax, θmin, θmax, "", 0., 0., Vec3(0.,0.,1.))
    function init(data::RICHGunData, ::Any)
        data.gun = move!(G4ParticleGun())
        SetParticlePosition(data.gun, G4ThreeVector(0, 0, 0))
    end
    function generate(evt::G4Event, data::RICHGunData)::Nothing
        gun = data.gun
        data.particle = rand(data.particles)
        SetParticleByName(gun, data.particle)
        data.mass = gun |> GetParticleDefinition |> GetPDGMass
        data.momentum = data.pmin * (data.pmax/data.pmin)^rand()
        cosθ = cos(data.θmin) - rand() * (cos(data.θmin) - cos(data.θmax))    ## uniform in solid angle
        sinθ, ϕ = sqrt(1 - cosθ^2), 2π * rand()
        data.direction = Vec3(sinθ * cos(ϕ), sinθ * sin(ϕ), cosθ)
        SetParticleEnergy(gun, sqrt(data.momentum^2 + data.mass^2) - data.mass)    ## kinetic energy
        SetParticleMomentumDirection(gun, G4ThreeVector(data.direction...))
        GeneratePrimaryVertex(gun, CxxPtr(evt))
    end
    G4JLPrimaryGenerator("RICHGun", data; init_method=init, generate_method=generate)
end
generator = RICHGunGenerator()

include(joinpath(@__DIR__, "RICHTracer.jl"))
include(joinpath(@__DIR__, "RICHTracerGPU.jl"))
using .RICHGPU
tracer = RICHTracer(rich)
gputracer = GPUTracer{T}(tracer)

toGPU(::Type{T}, p::OpticalPhoton) where T = GPUPhoton{T}(p.position, p.direction, p.polarization, c4f10RIndex(p.energy),
                                                          aerogelRayleigh(p.energy), quantumEfficiency(p.energy))
floattype(::GPUTracer{T}) where T = T

struct Hit
    position::Point3{Float64}
    energy::Float64
    radiator::Symbol            ## :aerogel or :gas
    primary::Bool               ## produced by the primary track
end
struct Photon
    event::Int
    θc::Float64                 ## reconstructed Cherenkov angle
    radiator::Symbol
    primary::Bool
end
struct Track
    particle::String
    mass::Float64
    momentum::Float64
    direction::Vec3{Float64}
end

mutable struct RICHSimData{B, GT} <: G4JLSimulationData
    ##---Configuration------------------------------------------------------------------------------
    mode::Symbol                                   ## :cpu, :gpu or :collect
    batchSize::Int                                 ## number of events per offload
    backend::B                                     ## KernelAbstractions backend
    tracer::RICHTracer                             ## CPU propagator
    gputracer::GT                                  ## GPU propagator
    opticalphoton::CxxPtr{G4ParticleDefinition}
    ##---Batch data---------------------------------------------------------------------------------
    pending::Vector{OpticalPhoton}                 ## photons waiting to be offloaded
    pendingEvent::Vector{Int}                      ## index of the event of each pending photon
    ##---Run data-----------------------------------------------------------------------------------
    tracks::Vector{Track}
    photons::Vector{Photon}
    sample::Vector{OpticalPhoton}                  ## photons kept in :collect mode
    nOffloads::Int                                 ## number of batches offloaded
    nOptPhotons::Int                               ## number of photons offloaded
    offloadTime::Float64                           ## time spent in the propagator (and the copies)
end
RICHSimData(backend, tracer, gputracer; mode=:gpu, batchSize=100) =
    RICHSimData(mode, batchSize, backend, tracer, gputracer, CxxPtr{G4ParticleDefinition}(C_NULL),
                OpticalPhoton[], Int[], Track[], Photon[], OpticalPhoton[], 0, 0, 0.)

function reflectionPoint(E, D, C, R)
    a, b = E - C, D - C
    â = normalize(a)
    ĉ = normalize(b - dot(b, â) * â)
    mirrorpoint(φ) = C + R * (cos(φ) * â + sin(φ) * ĉ)
    function g(φ)                                  ## difference of the cosines of the two angles
        M = mirrorpoint(φ)
        m = normalize(M - C)
        dot(normalize(E - M), m) - dot(normalize(D - M), m)
    end
    lo, hi = 0., atan(dot(b, ĉ), dot(b, â))
    for _ in 1:40
        mid = (lo + hi) / 2
        g(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return mirrorpoint((lo + hi) / 2)
end

function sphereIntersection(t, C, R)
    b = dot(t, C)
    return b + sqrt(b^2 - dot(C, C) + R^2)
end

const nGas = c4f10RIndex(3.2eV)

function cherenkovAngle(det::RICHDetector, track::Track, hit::Hit)
    t = track.direction
    side = hit.position[2] > 0 ? 1 : -1
    C = Point3{Float64}(mirrorCentre(det, side)...)
    zAerogel = det.aerogelZ + det.aerogelThickness
    if hit.radiator == :aerogel
        s = (det.aerogelZ + det.aerogelThickness/2) / t[3]
    else
        s = (zAerogel / t[3] + sphereIntersection(t, C, det.mirrorRadius)) / 2
    end
    E = Point3{Float64}(s * t)
    M = reflectionPoint(E, hit.position, C, det.mirrorRadius)
    d = normalize(M - E)
    if hit.radiator == :aerogel
        dt = Vec3(d[1], d[2], 0.) * nGas / det.aerogelRIndex
        d = Vec3(dt[1], dt[2], sqrt(1 - dot(dt, dt)))
    end
    return acos(clamp(dot(d, t), -1., 1.))
end

function expectedAngle(n, m, p)
    β = p / sqrt(p^2 + m^2)
    nβ = n * β
    return nβ > 1 ? acos(1/nβ) : NaN
end

radiatorOf(det, zvertex) = zvertex < det.aerogelZ + det.aerogelThickness + 1mm ? :aerogel : :gas

function flush!(data::RICHSimData, det::RICHDetector)
    isempty(data.pending) && return
    if data.mode == :collect
        append!(data.sample, data.pending)
    else
        detected = Tuple{Int, Point3{Float64}}[]
        data.offloadTime += @elapsed begin
            if data.mode == :gpu
                gpuphotons = [toGPU(floattype(data.gputracer), p) for p in data.pending]
                for (i, h) in enumerate(offload(data.backend, data.gputracer, gpuphotons))
                    h.detected && push!(detected, (i, Point3{Float64}(h.position)))
                end
            else
                for (i, p) in enumerate(data.pending)
                    pos = propagate(data.tracer, p)
                    isnothing(pos) || rand() < quantumEfficiency(p.energy) && push!(detected, (i, Point3{Float64}(pos)))
                end
            end
        end
        data.nOffloads += 1
        data.nOptPhotons += length(data.pending)
        for (i, pos) in detected
            photon, ievt = data.pending[i], data.pendingEvent[i]
            hit = Hit(pos, photon.energy, radiatorOf(det, photon.position[3]), photon.primary)
            push!(data.photons, Photon(ievt, cherenkovAngle(det, data.tracks[ievt], hit), hit.radiator, hit.primary))
        end
    end
    empty!(data.pending)
    empty!(data.pendingEvent)
    return
end

function beginrun(run::G4Run, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    data.opticalphoton = FindParticle("opticalphoton")
    empty!(data.tracks)
    empty!(data.photons)
    empty!(data.sample)
    data.nOffloads, data.nOptPhotons, data.offloadTime = 0, 0, 0.
    return
end

toVec3(v) = Vec3{Float64}(x(v), y(v), z(v))
function stackaction(track::G4Track, app::G4JLApplication)::G4ClassificationOfNewTrack
    data = getSIMdata(app)
    GetParticleDefinition(track) == data.opticalphoton || return fUrgent
    push!(data.pending, OpticalPhoton(toVec3(GetPosition(track)[]), toVec3(GetMomentumDirection(track)[]),
                                      toVec3(GetPolarization(track)[]), GetTotalEnergy(track),
                                      GetGlobalTime(track), GetParentID(track) == 1))
    push!(data.pendingEvent, length(data.tracks) + 1)
    return fKill
end

function endevent(evt::G4Event, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    gen = app.generator.data
    push!(data.tracks, Track(gen.particle, gen.mass, gen.momentum, gen.direction))
    length(data.tracks) % data.batchSize == 0 && flush!(data, app.detector)
    return
end

function endrun(run::G4Run, app::G4JLApplication)::Nothing
    flush!(getSIMdata(app), app.detector)
    return
end

app = G4JLApplication(; detector = rich,                              ## detector with parameters
                        simdata = RICHSimData(backend, tracer, gputracer),  ## simulation data structure
                        generator = generator,                        ## primary particle generator
                        nthreads = 0,                                 ## # of threads (0 = no MT)
                        physics_type = RICHPhysics,                   ## what physics list to instantiate
                        stackaction_method = stackaction,             ## stacking action
                        endeventaction_method = endevent,             ## end-event action
                        beginrunaction_method = beginrun,             ## begin run action
                        endrunaction_method = endrun                  ## end run action
                      );

configure(app)
initialize(app)

data = getSIMdata(app)
for mode in (:cpu, :gpu)
    data.mode = mode
    beamOn(app, 10)
end

results = Dict{Symbol, Any}()
for mode in (:cpu, :gpu)
    data.mode = mode
    t = @elapsed beamOn(app, 2000)
    results[mode] = (tracks = copy(data.tracks), photons = copy(data.photons), time = t)
    @printf "%-4s: %5.1f s for %d events, %d photons offloaded in %d batches, %.2f s in the propagator, %d hits\n" mode t length(data.tracks) data.nOptPhotons data.nOffloads data.offloadTime length(data.photons)
end

for mode in (:cpu, :gpu), rad in (:aerogel, :gas)
    (; tracks, photons) = results[mode]
    sel = count(p -> p.radiator == rad && p.primary, photons)
    nfired = count(t -> expectedAngle(rad == :aerogel ? rich.aerogelRIndex : nGas, t.mass, t.momentum) > 0, tracks)
    @printf "%-4s %-8s: %6.1f detected photons per track above threshold\n" mode rad sel/max(nfired, 1)
end

fig = Figure(size=(1200, 450))
for (i, (rad, n)) in enumerate(((:aerogel, rich.aerogelRIndex), (:gas, nGas)))
    ax = Axis(fig[1,i], title = rad == :aerogel ? "Aerogel" : "C₄F₁₀",
              xlabel="θc(reconstructed) - θc(expected) [mrad]", ylabel="Photons per track")
    for (mode, color) in ((:cpu, :black), (:gpu, :red))
        (; tracks, photons) = results[mode]
        Δθ = Float64[]
        for p in photons
            (p.radiator == rad && p.primary) || continue
            t = tracks[p.event]
            θ = expectedAngle(n, t.mass, t.momentum)
            isnan(θ) || push!(Δθ, (p.θc - θ)/mrad)
        end
        core = filter(x -> abs(x) < 5, Δθ)
        μ = sum(core)/length(core)
        σ = sqrt(sum(abs2, core .- μ)/length(core))
        h = Hist1D(Δθ; binedges=range(-10, 10, length=101)) * (1/length(tracks))
        stephist!(ax, h, color=color, label=@sprintf("%s (core σ = %.2f mrad)", mode, σ))
    end
    axislegend(ax, position=:rt)
end
display(fig)

data.mode = :collect
beamOn(app, 3000)
sample = copy(data.sample)
sizes = [10^k for k in 2:6 if 10^k <= length(sample)]
backends = [("CPU propagator", nothing, Float64), ("KernelAbstractions CPU()", CPU(), Float64)]
backend isa CPU || push!(backends, ("KernelAbstractions $backend", backend, T))

function timeit(f; repeat=3)
    f()                                            ## compilation and warm-up
    minimum(@elapsed(f()) for _ in 1:repeat)
end

fig = Figure(size=(800, 500))
ax = Axis(fig[1,1], xlabel="Photons per batch", ylabel="Time per photon [ns]", xscale=log10, yscale=log10,
          title="Propagation of optical photons")
for (label, bk, FT) in backends
    times = Float64[]
    for n in sizes
        photons = sample[1:n]
        if isnothing(bk)
            push!(times, timeit(() -> foreach(p -> propagate(tracer, p), photons)))
        else
            gp, gt = [toGPU(FT, p) for p in photons], GPUTracer{FT}(tracer)
            push!(times, timeit(() -> offload(bk, gt, gp)))
        end
    end
    scatterlines!(ax, sizes, times ./ sizes .* 1e9, label=label)
    @printf "%-30s: %s ns/photon\n" label join([@sprintf("%.0f", t/n*1e9) for (t, n) in zip(times, sizes)], " ")
end
axislegend(ax, position=:rt);

display(fig)
