using Geant4
using Geant4.SystemOfUnits
using Libdl
using CairoMakie, Rotations, LinearAlgebra, IGLWrap_jll  # to force loading G4Vis extension

prefix = Geant4.Geant4_jll.artifact_dir
dlext = Libdl.dlext
# Compilation of the custom library
Base.run(`c++ -O2 -shared -fPIC -std=c++17 -I$prefix/include/Geant4
         -Wl,-rpath,$prefix/lib -L$prefix/lib
         -lG4geometry -lG4materials -lG4global -lG4clhep
         -o UserLib.$dlext $(@__DIR__)/UserLib.cpp`).exitcode == 0 || error("Compilation failed")

const lib = "./UserLib.$(Libdl.dlext)"
createRoundCube(a,r) = @ccall lib.createRoundCube(a::Float64, r::Float64)::CxxPtr{G4VSolid}
deleteRoundCube(s::CxxPtr{G4VSolid}) = @ccall lib.deleteRoundCube(s::CxxPtr{G4VSolid})::Cvoid
infoRoundCube(s::CxxPtr{G4VSolid}) = (@ccall lib.infoRoundCube(s::CxxPtr{G4VSolid})::Cstring) |> unsafe_string

rcube = createRoundCube(10cm, 1cm)
img = drawDistanceToOut(rcube[], 100000)
display(img)

info = infoRoundCube(rcube)
println(info)

deleteRoundCube(rcube)
