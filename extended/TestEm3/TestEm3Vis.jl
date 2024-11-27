using Geant4
using Geant4.SystemOfUnits
using GLMakie, Rotations, LinearAlgebra, IGLWrap_jll  # to force loading G4Vis extension

include("TestEm3Detector.jl")

det = TestEm3Detector()
TestEm3Construct(det)

world = det.fPhysiWorld
img = draw(world[])
display(img)
