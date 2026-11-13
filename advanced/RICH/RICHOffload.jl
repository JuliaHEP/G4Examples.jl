using Geant4
using Geant4.SystemOfUnits
using Geant4.PhysicalConstants
using FHist, GeometryBasics, Printf
using LinearAlgebra: dot, normalize          ## Geant4 also exports `dot`
using GLMakie

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

mutable struct RICHSimData <: G4JLSimulationData
    ##---Configuration------------------------------------------------------------------------------
    offload::Bool                                  ## offload the optical photons (or track them in Geant4)
    tracer::Union{Nothing, OpticalTracer}          ## the external propagator
    opticalphoton::CxxPtr{G4ParticleDefinition}
    ##---Event data---------------------------------------------------------------------------------
    optphotons::Vector{OpticalPhoton}              ## photons collected by the stacking action
    hits::Vector{Hit}
    ##---Run data-----------------------------------------------------------------------------------
    tracks::Vector{Track}
    photons::Vector{Photon}
    eventHits::Vector{Vector{Hit}}                 ## hits of the first events, for the displays
    nOptPhotons::Int                               ## number of photons offloaded
    offloadTime::Float64                           ## time spent in the propagator
    RICHSimData(; offload=true) = new(offload, nothing, CxxPtr{G4ParticleDefinition}(C_NULL),
                                      OpticalPhoton[], Hit[], Track[], Photon[], [], 0, 0.)
end

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

function beginrun(run::G4Run, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    data.opticalphoton = FindParticle("opticalphoton")
    empty!(data.tracks)
    empty!(data.photons)
    empty!(data.eventHits)
    data.nOptPhotons, data.offloadTime = 0, 0.
    return
end

function beginevent(evt::G4Event, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    empty!(data.hits)
    empty!(data.optphotons)
    return
end

toVec3(v) = Vec3{Float64}(x(v), y(v), z(v))
function stackaction(track::G4Track, app::G4JLApplication)::G4ClassificationOfNewTrack
    data = getSIMdata(app)
    (data.offload && GetParticleDefinition(track) == data.opticalphoton) || return fUrgent
    push!(data.optphotons, OpticalPhoton(toVec3(GetPosition(track)[]), toVec3(GetMomentumDirection(track)[]),
                                         toVec3(GetPolarization(track)[]), GetTotalEnergy(track),
                                         GetGlobalTime(track), GetParentID(track) == 1))
    return fKill
end

radiatorOf(det, zvertex) = zvertex < det.aerogelZ + det.aerogelThickness + 1mm ? :aerogel : :gas

function stepaction(step::G4Step, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    data.offload && return
    track = GetTrack(step)
    GetParticleDefinition(track) == data.opticalphoton || return
    post = GetPostStepPoint(step)
    GetStepStatus(post) == fGeomBoundary || return
    pv = GetPhysicalVolume(post)
    (pv == C_NULL || !startswith(String(GetName(pv)), "Photodetector")) && return
    SetTrackStatus(track, fStopAndKill)
    energy = GetTotalEnergy(track)
    rand() < quantumEfficiency(energy) || return
    p = GetPosition(post)[]
    radiator = radiatorOf(app.detector, z(GetVertexPosition(track)[]))
    push!(data.hits, Hit(Point3{Float64}(x(p), y(p), z(p)), energy, radiator, GetParentID(track) == 1))
    return
end

function offload!(hits::Vector{Hit}, photons::Vector{OpticalPhoton}, tracer::OpticalTracer, det::RICHDetector)
    for photon in photons
        p = propagate(tracer, photon)
        isnothing(p) && continue
        rand() < quantumEfficiency(photon.energy) || continue
        push!(hits, Hit(Point3{Float64}(p...), photon.energy, radiatorOf(det, photon.position[3]), photon.primary))
    end
    return hits
end

function endevent(evt::G4Event, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    if data.offload
        data.nOptPhotons += length(data.optphotons)
        data.offloadTime += @elapsed offload!(data.hits, data.optphotons, data.tracer, app.detector)
    end
    gen = app.generator.data
    track = Track(gen.particle, gen.mass, gen.momentum, gen.direction)
    push!(data.tracks, track)
    ievt = length(data.tracks)
    for hit in data.hits
        push!(data.photons, Photon(ievt, cherenkovAngle(app.detector, track, hit), hit.radiator, hit.primary))
    end
    length(data.eventHits) < 100 && push!(data.eventHits, copy(data.hits))
    return
end

app = G4JLApplication(; detector = rich,                              ## detector with parameters
                        simdata = RICHSimData(),                      ## simulation data structure
                        generator = generator,                        ## primary particle generator
                        nthreads = 0,                                 ## # of threads (0 = no MT)
                        physics_type = RICHPhysics,                   ## what physics list to instantiate
                        stackaction_method = stackaction,             ## stacking action
                        stepaction_method = stepaction,               ## step action method
                        begineventaction_method = beginevent,         ## begin-event action
                        endeventaction_method = endevent,             ## end-event action
                        beginrunaction_method = beginrun              ## begin run action
                      );

configure(app)
initialize(app)

data = getSIMdata(app)
data.tracer = RICHTracer(rich, exportGDML())
data.tracer.model

generator.data.particles, generator.data.pmin, generator.data.pmax = ["pi+"], 10GeV, 10GeV
beamOn(app, 1)
optphotons, hits = copy(data.optphotons), copy(data.hits)
@printf "Photons collected: %d (%d in the aerogel), hits: %d\n" length(optphotons) count(p -> radiatorOf(rich, p.position[3]) == :aerogel, optphotons) length(hits)

fig = Figure(size=(1200, 900))
s = LScene(fig[1,1])
Geant4.draw!(s, GetWorldVolume())
scatter!(s, [Point3{Float64}(p.position) for p in optphotons], markersize=2,
         color=[radiatorOf(rich, p.position[3]) == :aerogel ? :orange : :blue for p in optphotons])
scatter!(s, [h.position for h in hits], markersize=6, color=[h.radiator == :aerogel ? :orange : :blue for h in hits])
display(fig)

generator.data.particles, generator.data.pmin, generator.data.pmax = ["pi+", "kaon+", "proton"], 1GeV, 80GeV
results = Dict{Symbol, Any}()
for mode in (:geant4, :offload)
    data.offload = mode == :offload
    t = @elapsed beamOn(app, 1000)
    results[mode] = (tracks = copy(data.tracks), photons = copy(data.photons), eventHits = copy(data.eventHits), time = t)
    @printf "%-8s: %6.1f s for %d events, %d hits" mode t length(data.tracks) length(data.photons)
    mode == :offload && @printf " (%d photons offloaded, %.1f s in the tracer)" data.nOptPhotons data.offloadTime
    println()
end

for mode in (:geant4, :offload), rad in (:aerogel, :gas)
    (; tracks, photons) = results[mode]
    sel = count(p -> p.radiator == rad && p.primary, photons)
    nfired = count(t -> expectedAngle(rad == :aerogel ? rich.aerogelRIndex : nGas, t.mass, t.momentum) > 0, tracks)
    @printf "%-8s %-8s: %6.1f detected photons per track above threshold\n" mode rad sel/max(nfired, 1)
end

fig = Figure(size=(1200, 450))
for (i, mode) in enumerate((:geant4, :offload))
    ax = Axis(fig[1,i], title = mode == :geant4 ? "Geant4 tracking" : "Offloaded", xlabel="x [cm]", ylabel="y [cm]", aspect=DataAspect())
    xlims!(ax, -rich.detectorHalfX/cm, rich.detectorHalfX/cm)
    ylims!(ax, rich.detectorY ./ cm |> extrema)
    for hits in results[mode].eventHits[1:20]
        h = [hit for hit in hits if hit.position[2] > 0]
        isempty(h) && continue
        scatter!(ax, [p.position[1]/cm for p in h], [p.position[2]/cm for p in h], markersize=5,
                 color=[p.radiator == :aerogel ? :orange : :blue for p in h])
    end
end
display(fig)

fig = Figure(size=(1200, 450))
for (i, (rad, n)) in enumerate(((:aerogel, rich.aerogelRIndex), (:gas, nGas)))
    ax = Axis(fig[1,i], title = rad == :aerogel ? "Aerogel" : "C₄F₁₀",
              xlabel="θc(reconstructed) - θc(expected) [mrad]", ylabel="Photons per track")
    for (mode, color) in ((:geant4, :black), (:offload, :red))
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
