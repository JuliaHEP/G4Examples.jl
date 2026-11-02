# Example XrayTES

Port of the Geant4 advanced example
[xray_TESdetector](https://github.com/Geant4/geant4/tree/master/examples/advanced/xray_TESdetector)
(P. Dondero and R. Stanzani, developed within the ESA AREMBES project). It estimates the particle
background of an X-ray detector in space: a simplified mass model of the X-IFU spectrometer of the
ATHENA mission, from an early configuration.

## Detector

The geometry is defined in `XrayTESDetector.jl`. The parameters of the TES array are keyword
parameters of the `XrayTESDetector` structure. The dimensions of the other elements come from a CAD
model and are fixed.

- **TES array**: 317 pixels with a 0.814 mm pitch in a hexagon of 8.5 mm apothem. Each pixel has
  a 3 μm Bi absorber, a 0.75 μm Si₃N₄ membrane and a piece of the 200 μm Si grid (a
  `G4SubtractionSolid`). The pixels are a `G4PVParameterised` positioned by a Julia function.
  The positions are computed: they are the cells of the square lattice that fit inside a hexagon
  of 8.4 mm apothem. They are the same as in `pixelpos.txt` of the C++ example.
- **Anti-coincidence detector (ACD)**: four 200 μm Si trapezoids 0.75 mm below the array, on a
  Cu plate.
- **Supports**: Cu cage of columns and walls, Si boards.
- **Cryostat**: three thermal shields (Cryoperm, extruded Al shape, thin Al), a ring and the outer
  Al vessel.
- **Filter**: 20 nm Al on a stainless steel mesh (about 2000 trapezoidal bars) in an Al carrier,
  and a Ti aperture cylinder with an Au coating.
- **Spacecraft**: a 1 cm thick Al sphere of 26 cm radius.

The volumes seen by the detector (the array, the wafer, the ACD, the filter and the coating of the
aperture cylinder) form the `InnerRegion` (a `G4Region`), with production cuts of 1 mm.

## Primary particles

Galactic cosmic-ray protons from 10 MeV to 100 GeV with the PAMELA 2009 spectrum at L2. As in
the C++ example, they start on a sphere of 27 cm radius and go towards another random point of
the sphere, which gives an isotropic and uniform flux inside. Only about 1 proton in 1000 hits the
detector, so by default the generator aims at a target sphere of 5 cm around the detector
(`GCRGenerator(target=5cm)`). The flux stays isotropic and uniform inside the target. The
exposure time of the run is computed from the integrated intensity (0.407 protons/(cm² s sr)).

## Physics

QBBC with `G4EmStandardPhysics_option4` and `G4RadioactiveDecayPhysics`, configured with
`G4EmParameters`. In the `InnerRegion`, fluorescence, Auger electrons and PIXE are active and the
multiple scattering is replaced by single scattering (`G4EmStandardSS`). The production cuts are
1 mm everywhere, with a lower limit of 200 eV.

## Output

The stepping action records the tracks that deposit energy in the absorbers, the energy in each
pixel per event and the energy in the ACD. The script shows:

- the generated proton spectrum compared with the GCR spectrum
- the energy and the number of tracks in each pixel
- the particles and creator processes of the background
- the spectra at creation, at the absorber and of the energy deposited (as the `analysis.C` macro
  of the C++ example)
- the background spectrum in counts/(cm² s keV), with and without a 20 keV ACD veto, compared with
  the X-IFU requirement of 5×10⁻³ counts/(cm² s keV) in 2–10 keV

## Differences with the C++ example

- The C++ example uses its own `G4EmStandardPhysics_SpacePhysics`. Its models are very similar to
  those of option 4.
- The C++ example calls `SetAugerCascade(false)` after `SetAuger(true)`. In recent Geant4 versions
  `SetAugerCascade` is the same as `SetAuger`, so this switches the Auger electrons off. Here they
  are switched on.
- The pixel positions are computed instead of read from `pixelpos.txt`.
- The step limits of the C++ example are not used. They have no effect there either, because its
  physics list does not include a step limiter.
- No ROOT file is written. The results are kept in Julia structures and plotted with Makie.
- The geometry can only be built in Julia. Reading the GDML file of the C++ example is not
  supported.

## Running

```sh
julia --project=. advanced/XrayTES/XrayTES.jl
```

The script uses 8 threads. 50000 protons take about 10 s on a laptop and give about 700 tracks in
the absorbers. Increase `nevents` for better statistics, in particular after the ACD veto.
