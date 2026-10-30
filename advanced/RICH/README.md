# Example RICH

A Ring Imaging Cherenkov detector with two radiators (aerogel and C₄F₁₀ gas), similar to the
RICH1 detector of LHCb in its Technical Design Report (see `Schematic-RICH-detector.png`).

Pions, kaons and protons from the origin cross a 5 cm aerogel tile (n = 1.03) and about 85 cm of
C₄F₁₀ (n ≈ 1.0014). Two tilted spherical mirrors (R = 2 m) focus the Cherenkov light onto two
photodetector planes outside the acceptance, where each radiator gives a ring. The example:

- builds the geometry with boolean solids (mirror segments cut from a spherical shell, an aerogel
  tile with a hole for the beam pipe) and a reflecting optical surface on the mirrors
- defines the optical properties as functions of the photon energy: dispersion of C₄F₁₀,
  Rayleigh scattering in aerogel and the quantum efficiency of the photodetectors
- records the detected photons in the stepping action
- reconstructs the Cherenkov angle of each photon from its hit, the track and the mirror geometry
  (finding the reflection point on the mirror and correcting for the refraction at the aerogel exit)
- shows an event display, the rings on the photodetectors, the Cherenkov angle against the
  momentum, the single photon resolution and the reconstructed mass of each track

## Files

- `RICH.lit`: Literate source of the example (source of truth)
- `RICH.jl`, `RICH.ipynb`: script and notebook generated from `RICH.lit`
- `RICHDetector.jl`: detector parameters, materials, optical properties and geometry
- `Schematic-RICH-detector.png`: schematic of the detector
