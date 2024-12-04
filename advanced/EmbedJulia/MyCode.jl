#---Load the needed Julia modules------------------------------------------------------------------
using Geant4
using FHist
using Plots
using Parameters

println("=====> Loading MyCode.jl")

#---Simulation data struct-------------------------------------------------------------------------
@with_kw mutable struct MyData
    edep = 0.0
    edepHist = H1D("Event total Edep distribution", 100, 0., 110.)
end
add!(d::MyData, d2::MyData) = (d.edep += d2.edep; merge!(d.edepHist, d2.edepHist))

const nthreads = ENV["G4NUMTHREADS"] == nothing ? 0 : parse(Int, ENV["G4NUMTHREADS"])
const simdata = [MyData() for i in 1:nthreads+1]   # #workers + 1(master)

function getMyData()
    tid = G4Threading!G4GetThreadId()
    tid < 0 && (tid = -1)   # master thread (-2 for without multi-threading support)
    simdata[tid+2]
end 

#---Actions----------------------------------------------------------------------------------------
function end_of_event_action(event)
    data = getMyData()
    push!(data.edepHist, data.edep)
    return                      # This is mandatory to force to return nothing
end

function begin_of_event_action(event)
    data = getMyData()
    data.edep = 0.0
    return                      # This is mandatory to force to return nothing
end

function stepping_action(step)
    data = getMyData()
    data.edep += GetTotalEnergyDeposit(step)
    return                      # This is mandatory to force to return nothing
end

function begin_of_run_action(run)
    return                      # This is mandatory to force to return nothing
end

function end_of_run_action(run)
    if G4Threading!G4GetThreadId() < 0
        println("=====> End of run")
        data = simdata[1]
        ##---This is the master thread, so we need to add all the simulation results-----------------
        for d in simdata[2:end]
            add!(data, d)
        end
        h = data.edepHist
        img = plot(h.hist, title=h.title)
        savefig(img, "edepHist.png")
        println("=====> edepHist.png saved")
    end    
    return                      # This is mandatory to force to return nothing
end  
