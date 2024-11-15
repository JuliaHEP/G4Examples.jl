using Geant4
using Geant4.SystemOfUnits
using GLMakie, Rotations, LinearAlgebra, IGLWrap_jll  # to force loa ding G4Vis extension

include("DetectorB1.jl")

pv = constructB1Detector(nothing)
img = draw(pv[])
display(img)
