using Geant4
using Geant4.SystemOfUnits:  cm, cm3, mm, pGy, eplus, keV, g, eV
#using GLMakie, Rotations, IGLWrap_jll  # to force loding G4Vis extension

#---Define Detector Parameters struct--------------------------------------------------------------
include(joinpath(@__DIR__, "DetectorB3.jl"))

#---Physics List-----------------------------------------------------------------------------------
struct PhysicsB3a <: G4VUserPhysicsList
    function PhysicsB3a(verbose)
        pl = G4VModularPhysicsList()
        RegisterPhysics(pl, move!(G4DecayPhysics(verbose)))           # Default physics
        RegisterPhysics(pl, move!(G4EmStandardPhysics(verbose)))      # EM physics
        RegisterPhysics(pl, move!(G4RadioactiveDecayPhysics(verbose))) # Radioactive decay
        return pl
    end
end

#---Primary Particle Generator---------------------------------------------------------------------
@with_kw mutable struct GeneratorB3aData <: G4JLGeneratorData
    gun::Union{Nothing, CxxPtr{G4ParticleGun}} = nothing
    ion::Union{Nothing, CxxPtr{G4ParticleDefinition}} = nothing
    Z::Int64 = 9
    A::Int64 = 18
    ionCharge::Float64 = 0eplus
    excitEnergy::Float64 = 0keV
    position::G4ThreeVector = G4ThreeVector(4cm,4cm,4cm)
    direction::G4ThreeVector = G4ThreeVector(1,0,0)
end
function GeneratorB3a(;kwargs...)
    data = GeneratorB3aData(;kwargs...)
    function _init(data::GeneratorB3aData, ::Any)
        gun = data.gun = move!(G4ParticleGun())
        SetParticleMomentumDirection(gun, G4ThreeVector(1,0,0))
        SetParticleEnergy(gun, 1eV)
    end
    function _gen(evt::G4Event, data::GeneratorB3aData)::Nothing
        if isnothing(data.ion)  # late initialize (after physics processes)
            ion = data.ion = GetIon(data.Z, data.A, data.excitEnergy)
            SetParticleDefinition(data.gun, ion)
            SetParticleCharge(data.gun, data.ionCharge)
        end
        position = data.position + G4ThreeVector((rand()-0.5)*1cm, (rand()-0.5)*1cm, (rand()-0.5)*1cm)
        SetParticlePosition(data.gun, position)
        GeneratePrimaryVertex(data.gun, CxxPtr(evt))
    end
    G4JLPrimaryGenerator("GeneratorB3a", data; init_method=_init, generate_method=_gen)
end

#---Event Display----------------------------------------------------------------------------------
#display = G4JLEventDisplay(joinpath(@__DIR__, "VisSettings.jl"))

#---SimulationData---------------------------------------------------------------------------------
mutable struct SimDataB3a <: G4JLSimulationData
    #---Run data
    goodEvents::Int64
    sumDose::Float64
    SimDataB3a() = new(0,0)
end
function add!(x::SimDataB3a, y::SimDataB3a)
    x.goodEvents += y.goodEvents
    x.sumDose += y.sumDose
end

#---Sensitive Detector Crystal---------------------------------------------------------------------
struct CrystalData <: G4JLSDData
    edep::Dict{Int64,Float64} # (CopyNo, Edep)
    CrystalData() = new(Dict{Int64,Float64}())
end
function c_initialize(::G4HCofThisEvent, data::CrystalData)::Nothing
    empty!(data.edep)
    return
end
function c_processHits(step::G4Step, ::G4TouchableHistory, data::CrystalData)::Bool
    edep = step |> GetTotalEnergyDeposit
    edep <  0. && return false
    copy = step |> GetPreStepPoint |> GetTouchable |> GetCopyNumber
    data.edep[copy] = haskey(data.edep, copy) ? data.edep[copy] + edep : edep
    return true
end
#---Create SD instance
crystalSD = G4JLSensitiveDetector("CrystalSD", CrystalData();          # SD name an associated data are mandatory
                                   processhits_method=c_processHits,   # process hist method (also mandatory)
                                   initialize_method=c_initialize)     # intialize method
#---Sensitive Detector Patient---------------------------------------------------------------------
mutable struct PatientData <: G4JLSDData
    dose::Float64
    PatientData() = new(0)
end
function p_initialize(::G4HCofThisEvent, data::PatientData)::Nothing
    data.dose = 0
    return
end
function p_processHits(step::G4Step, ::G4TouchableHistory, data::PatientData)::Bool
    edep = step |> GetTotalEnergyDeposit
    edep <=  0. && return false
    volume  = step |> GetPreStepPoint |> GetTouchable |> GetSolid |> GetCubicVolume
    density = step |> GetPreStepPoint |> GetMaterial |> GetDensity
    data.dose += edep /(density * volume)
    return true
end
#---Create SD instance
patientSD = G4JLSensitiveDetector("PatientSD", PatientData();           # SD name an associated data are mandatory
                                   processhits_method=p_processHits,    # process hist method (also mandatory)
                                   initialize_method=p_initialize)      # intialize method
  
#---User Actions-----------------------------------------------------------------------------------
#---Run Action
function beginrun(run::G4Run, app::G4JLApplication)::Nothing
    data = getSIMdata(app)
    data.goodEvents = 0
    data.sumDose = 0.
    nothing
end
μGy = Geant4.SystemOfUnits.gray/1e6
function endrun(run::G4Run, app::G4JLApplication)::Nothing
    partName = app.generator.data.gun |> GetParticleDefinition |> GetParticleName |> String
    #---end run action is called for each workwer thread and the master one
    if G4Threading!G4GetThreadId() < 0
        data = app.simdata[1]
        #---This is the master thread, so we need to add all the simuation results-----------------
        for d in app.simdata[2:end]
            add!(data, d)
        end
        noEvents = run |> GetNumberOfEvent
        G4JL_println("""
                     --------------------End of Run------------------------------
                      The run was $noEvents $partName Nb of 'good' e+ annihilations: $(data.goodEvents)
                      Total dose in patient : $(data.sumDose/pGy) pGy
                     ------------------------------------------------------------ 
                     """)
    end
end 
#---Event Action
function endevent(evt::G4Event, app::G4JLApplication)
    edep = getSDdata(app, "CrystalSD").edep
    dose = getSDdata(app, "PatientSD").dose
    data = getSIMdata(app)
    if count(>(500keV), values(edep)) == 2
        data.goodEvents += 1
    end
    data.sumDose += dose
    return
end
#---Stacking Action
let G4NeutrinoE, first=true
global function stacking(trk::G4Track, app::G4JLApplication)::G4ClassificationOfNewTrack
    if first  # emulation of C++ static 
        G4NeutrinoE = FindParticle("nu_e")
        first = false
    end
    (trk |> GetParentID) == 0 && return fUrgent               # keep primary particle
    (trk |> GetDefinition) == G4NeutrinoE && return fKill     # kill neutrino
    return fUrgent
end
end

#---Application------------------------------------------------------------------------------------
app = G4JLApplication(; detector = DetectorB3(),                      # detector with parameters
                        simdata = SimDataB3a(),                       # simulation data structure
                        generator = GeneratorB3a(),                   # primary particle generator
                        nthreads = 0,                                 # # of threads (0 = no MT)
                        physics_type = PhysicsB3a,                    # what physics list to instantiate
                        #evtdisplay =  display,                        # set event display 
                        endeventaction_method = endevent,             # end-event action (fill histograms per event)
                        beginrunaction_method = beginrun,             # begin run action
                        endrunaction_method = endrun,                 # end run action
                        stackaction_method = stacking,                # track classification action
                        sdetectors = ["CrystalLV+" => crystalSD,
                                      "PatientLV" => patientSD]       # mapping of LVs to SDs (+ means multiple LVs with same name)
                      );

configure(app)
initialize(app)
beamOn(app, 10000)

#ui`/tracking/verbose 0`

