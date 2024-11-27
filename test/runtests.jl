using Test
using G4Examples

@testset "G4Examples" begin
    #---change the working directory
    cd(dirname(dirname(pathof(G4Examples))))

    @test run(`julia basic/B1/B1.jl`, devnull, devnull).exitcode  == 0
    @test run(`julia basic/B1/B1vis.jl`, devnull, devnull).exitcode == 0
    @test run(`julia basic/B2/B2a.jl`, devnull, devnull).exitcode  == 0
    @test run(`julia basic/B2/B2aVis.jl`, devnull, devnull).exitcode == 0
    @test run(`julia basic/B3/B3a.jl`, devnull, devnull).exitcode  == 0

    @test run(`julia extended/RE03/RE03.jl`, devnull, devnull).exitcode  == 0
    @test run(`julia extended/GPS/GPS.jl`, devnull, devnull).exitcode == 0
    @test run(`julia extended/TestEm3/TestEm3.jl`, devnull, devnull).exitcode  == 0
    @test run(`julia extended/TestEm3/TestEm3Vis.jl`, devnull, devnull).exitcode == 0

    @test run(`julia advanced/Scintillation/Scintillation.jl`, devnull, devnull).exitcode  == 0
    @test run(`julia advanced/HBC30/HBC30.jl`, devnull, devnull).exitcode == 0
    @test run(`julia advanced/AlephTPC/TPCSim.jl`, devnull, devnull).exitcode == 0
    @test run(`julia advanced/WaterPhantom/WaterPhantom.jl`, devnull, devnull).exitcode == 0
end
