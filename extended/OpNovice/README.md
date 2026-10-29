# Example OpNovice2

Translated to Julia from the **extended/optical/OpNovice2** example of Geant4.

Investigate optical properties and parameters. Details of optical
photon boundary interactions on a surface are recorded. Details
of optical photon generation and transport are recorded.

The C++ UI commands that set the optical properties (`/opnovice2/boxProperty`,
`/opnovice2/surfaceFinish`, ...) are keyword arguments of the `OpDetector` structure
in Julia. The example reproduces four of the C++ macros:

- `electron.mac`: 1 MeV electrons in a scintillating BGO tank. It prints the Cerenkov and
  scintillation photon yields and plots their spectra and the scintillation time profile.
- `fresnel.mac`: polarized optical photons with random direction hit a polished surface.
  It plots the fraction of refracted and reflected photons against the incident angle, and
  compares the Brewster angle and the critical angle with the expected values.
- `coated.mac`: photons hit a surface with a 100 nm coating. The reflectance at normal
  incidence is compared with the thin-film interference formula. In addition to the C++
  macro, a scan in incident angle for both polarizations is plotted against the formula.
- `boundary.mac`: photons hit the surface with each of the 21 configurations of the macro
  (unified, glisur, DAVIS and LUT models with their types and finishes). The fractions of
  the boundary process outcomes are printed and shown as a heatmap.

## Files

- `OpNovice2.lit`: Literate source of the example (source of truth)
- `OpNovice2.jl`, `OpNovice2.ipynb`: script and notebook generated from `OpNovice2.lit`
- `OpDetector.jl`: detector construction with configurable optical properties

## Differences with the C++ example

- The `Run`, `HistoManager` and `TrackInformation` classes are replaced by the
  `OpNovice2SimData` structure, with `FHist` histograms, and a pre-tracking action.
- The run uses one thread.
- The Birks constant of the tank material is not set to 0.126 mm/MeV
  (`G4IonisParamMat::SetBirksConstant` is not wrapped in Geant4.jl). Geant4 uses its
  default for G4_BGO (0.0084 mm/MeV), so there are about 15% more scintillation photons.
  With the same Birks constant, the C++ and Julia versions agree within statistics.
- The group velocity check and the per-step printout of the number of Cerenkov and
  scintillation photons (`G4Cerenkov/G4Scintillation::GetNumPhotons`) are not included.
- `G4Cerenkov` builds its tables only for the materials that exist at the first
  initialization, so the electron study runs first.
- In `boundary.mac` without an angle, the C++ generator chooses a random polarization for
  the first photon and keeps it for the rest of the run. The Julia generator chooses a new
  random polarization for each photon.
