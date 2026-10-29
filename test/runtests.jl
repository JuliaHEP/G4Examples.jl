using Test
using G4Examples
using Literate

#---Generate the script X.jl from X.lit (if it exists) and run it
function run_example(script)
    dir, file = splitdir(script)
    name = splitext(file)[1]
    lit = joinpath(dir, "$name.lit")
    if isfile(lit)
        Literate.script(lit, dir, name = name, keep_comments = false, credit = false)
    end
    @info "Running example $script"
    return run(`julia $script`, devnull, devnull).exitcode
end

@testset "G4Examples" begin
    #---change the working directory
    cd(dirname(dirname(pathof(G4Examples))))
    #---run the examples
    @test run_example("basic/B1/B1.jl")  == 0
    @test run_example("basic/B1/B1vis.jl") == 0
    @test run_example("basic/B2/B2a.jl")  == 0
    @test run_example("basic/B2/B2aVis.jl") == 0
    @test run_example("basic/B3/B3a.jl")  == 0
    @test run_example("basic/B4/B4a.jl")  == 0

    @test run_example("extended/RE03/RE03.jl")  == 0
    @test run_example("extended/GPS/GPS.jl") == 0
    @test run_example("extended/TestEm3/TestEm3.jl")  == 0
    @test run_example("extended/TestEm3/TestEm3-gdml.jl") == 0
    @test run_example("extended/TestEm3/TestEm3Vis.jl") == 0
    @test run_example("extended/field02/field02.jl") == 0
    @test run_example("extended/Solids/Solids.jl") == 0
    @test run_example("extended/OpNovice/OpNovice2.jl") == 0

    @test run_example("advanced/Scintillation/Scintillation.jl")  == 0
    @test run_example("advanced/HBC30/HBC30.jl") == 0
    @test run_example("advanced/AlephTPC/TPCSim.jl") == 0
    @test run_example("advanced/WaterPhantom/WaterPhantom.jl") == 0
    @test run_example("advanced/CustomLib/UserLib.jl") == 0
    @test run_example("advanced/EmbedJulia/JuliaAction.jl") == 0
end
