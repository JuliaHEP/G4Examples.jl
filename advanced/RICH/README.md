# Example RICH

Simulation of a Ring Imaging Cherenkov (RICH) detector with two radiators, silica aerogel and
C₄F₁₀ gas. The layout follows the RICH1 detector of LHCb as described in its Technical Design
Report.

![Schematic of the RICH detector](Schematic-RICH-detector.png)

## Physics

A charged particle that is faster than light in a medium (β > 1/n) emits Cherenkov photons
on a cone with half-angle cos θc = 1/(nβ). If the momentum is measured, the angle
gives the velocity and therefore the mass of the particle. A spherical mirror focuses all photons
with the same direction to the same point, so each cone appears as a ring on the photodetectors.
The radius of the ring is about θc·R/2, where R is the radius of curvature of the mirror.

Each radiator covers a different momentum range:

| Radiator | n | θc max (β = 1) | π threshold | K threshold | p threshold |
|---|---|---|---|---|---|
| Aerogel | 1.03 | 241 mrad | 0.6 GeV/c | 2.0 GeV/c | 3.8 GeV/c |
| C₄F₁₀ (at 3.2 eV) | 1.0014 | 53 mrad | 2.6 GeV/c | 9.3 GeV/c | 17.6 GeV/c |

## Detector

The detector is defined in `RICHDetector.jl`. All dimensions are keyword parameters of the
`RICHDetector` structure.

- **Gas vessel**: a box filled with C₄F₁₀ from z = 1 m to z = 2.1 m. The refractive index
  depends on the wavelength (Sellmeier formula).
- **Beam pipe**: a beryllium tube of 3 cm radius along z.
- **Aerogel**: a 5 cm tile at the entrance of the vessel, with a hole for the beam pipe (a
  `G4SubtractionSolid`). Rayleigh scattering length λ⁴/C, about 5 cm at 400 nm.
- **Mirrors**: two spherical mirrors (top and bottom) with R = 2 m and 90% reflectivity. Each is
  a segment of a spherical shell (a `G4Sphere` with θ and ϕ ranges, its polar axis rotated to the y
  axis) covering one half of the vessel, with a polished `dielectric_metal` optical surface. The
  centres of curvature are at y = ±75 cm, so each mirror is tilted and focuses its light onto the
  photodetector plane of the same half.
- **Photodetectors**: two planes at z ≈ 1 m, between |y| = 45 cm and 135 cm, outside the
  acceptance. Their quantum efficiency peaks at 25% around 400 nm.

## What the example does

1. Shoots one π⁺, K⁺ or proton per event from the origin. The momentum is uniform in log p between
   1 and 80 GeV/c, and the polar angle is between 50 and 250 mrad.
2. Records each optical photon that reaches a photodetector in the stepping action, and keeps
   it with a probability equal to the quantum efficiency.
3. Reconstructs the Cherenkov angle of each detected photon from the hit position and the track.
   The emission point is taken at the middle of the radiator. The reflection point on the mirror
   is found by a bisection in the plane of the emission point, the hit and the centre of
   curvature. The direction of aerogel photons is corrected for the refraction at the exit
   of the tile.
4. Produces the following plots:
   - the optical properties: C₄F₁₀ index, aerogel scattering length, quantum efficiency
   - an event display of a 10 GeV/c pion, with the geometry cut in half
   - the hits of 20 events on the photodetectors, where the aerogel and gas rings are visible
   - the reconstructed θc against the momentum, with the expected curves for π, K and p
   - the single photon resolution for each radiator
   - the mass reconstructed for each track from its average θc and momentum

With 1000 events the example gives about 21 detected photons per track from the aerogel and
27 from the gas. The single photon resolution is about 1 mrad for the aerogel and 2 mrad for the
gas, and π, K and p are well separated in both radiators.

## Running

```sh
julia --project=. advanced/RICH/RICH.jl
```

The script uses GLMakie and opens a window for each figure. It takes about 35 s. The notebook
`RICH.ipynb` uses CairoMakie and shows the figures inline.

To try other configurations, change the detector or generator parameters in `RICH.lit`. For
example:

```julia
rich = RICHDetector(aerogelRIndex = 1.01, mirrorReflectivity = 0.8)
generator = RICHGunGenerator(particles = ["pi+", "kaon+"], pmin = 5GeV, pmax = 30GeV)
```

## Offloading the optical photons

`RICHOffload.lit` is a variant of the example in which Geant4 does not track the optical photons.

1. The **stacking action** (`stackaction_method`) receives each new optical photon, copies its
   position, direction, polarization, energy, time and origin (primary or not) into an
   `OpticalPhoton` (a plain bits structure), and returns `fKill`.
2. At the end of the event, the whole vector of photons is passed to `offload!`, which runs code
   outside Geant4. Here this is the ray tracer in `RICHTracer.jl`. It implements absorption,
   Rayleigh scattering, Fresnel reflection and refraction with polarization, mirror reflectivity,
   detection and absorption at surfaces without a refractive index, using the same algorithms as
   the Geant4 optical processes. The photons that reach a photodetector pass the
   quantum-efficiency test and become hits.
3. The hits are reconstructed as in `RICH.lit`. The same application can also run with the
   photons tracked by Geant4 (`data.offload = false`), so the script runs both modes and compares
   them.

The tracer gets the geometry from [Geom4hep.jl](https://github.com/peremato/Geom4hep), a geometry
package written in Julia. After the initialization of the application, `exportGDML()` writes the
Geant4 geometry to a GDML file with the Geant4 GDML writer, and Geom4hep reads it and builds a flat
geometry model, which provides the navigation: the volume of a point, the distance to the next
boundary and the normal of a surface. The tracer adds:

- the optical properties of the materials, from `richOptics(det)` in `RICHDetector.jl`, the same
  functions used to build the Geant4 material properties tables. They are tabulated on the same
  energy grid and interpolated linearly, as in Geant4.
- the volumes with special surfaces, by name: the two mirrors (with their reflectivity) and the
  photodetectors.

Nothing in the tracer is specific to this detector, so it can be used with other geometries.
Geom4hep is not yet registered: add it with
`julia --project=. -e 'using Pkg; Pkg.add(url="https://github.com/peremato/Geom4hep.git")'`.
For the same reason the offload examples are not run by the tests nor included in the documentation.

When the same Geant4-generated photons are tracked by Geant4 and propagated with the tracer, the
expected numbers of hits agree within 0.1% (50 879 and 50 866 for 1500 events, with the
detection probability of each photon in the tracer estimated from 20 propagations). The two modes
give the same photon yields and the same single-photon resolutions (1.01 mrad for the aerogel and
1.96 mrad for the gas, within 0.01 mrad). The offloaded run takes 0.7 s instead of 1.3 s for 1000
events, of which 0.4 s is spent in the tracer.

The interface between Geant4 and the external code is just the `Vector{OpticalPhoton}`. To use
another propagator (on several threads or on a GPU, or writing the photons to a file), replace
`offload!`.

## Offloading the optical photons to a GPU

`RICHOffloadGPU.lit` runs the same propagator in a
[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) kernel with one thread
per photon (`RICHTracerGPU.jl`), including the Geom4hep navigation: the tracer and the geometry
model are copied to the device with `Adapt.adapt`. It uses `CUDABackend()` on NVIDIA GPUs or
`MetalBackend()` on Apple GPUs, in Float32, and otherwise propagates the photons with Julia threads
in Float64. On data-centre NVIDIA GPUs with fast double precision (A100, H100) you can set
`T = Float64`.

- The tracer is parametric in the floating point type (Apple GPUs do not support Float64).
  Geom4hep uses a larger geometrical tolerance in single precision.
- The random numbers come from a counter-based generator (Philox2x32-10). Its key is a seed for the
  batch and its stream is the index of the photon, so the kernel keeps no generator state in memory,
  and the results of a photon are the same on the CPU and on the GPU.
- The photons of `batchSize` events (100 by default, about 40 000 photons) are offloaded together.
  The end-of-run action offloads the last batch. Larger batches are split into kernel launches of
  at most 200 000 photons: macOS aborts or hangs long-running kernels on the GPU that drives the
  display.

The script compares the CPU and GPU propagators, which give the same photon yields and
resolutions, and measures the time per photon as a function of the batch size. On an Apple M2:

| Photons per batch | 10² | 10³ | 10⁴ | 10⁵ | 10⁶ |
|---|---|---|---|---|---|
| CPU propagator (1 thread) | 1412 ns | 1087 ns | 1112 ns | 1095 ns | 1113 ns |
| CPU propagator (8 threads) | 1268 ns | 369 ns | 249 ns | 218 ns | 208 ns |
| GPU kernel, `MetalBackend()` | 35045 ns | 2028 ns | 414 ns | 280 ns | 273 ns |

The general navigation of Geom4hep is about 3.5 times slower than the hand-written intersections of
the first version of the tracer (about 300 ns per photon on one thread and 60 ns on the GPU). A
photon from the gas takes about 600 ns on one thread, one from the aerogel about 1400 ns: it
scatters and crosses the faces of the tile several times, and each time it leaves the tile it has to
be located again among the volumes of the vessel. Describing the mirrors as `G4Sphere` segments
instead of spherical shells intersected with boxes halved the time of the gas photons. The GPU is
faster than one thread above about 10⁴ photons per batch, and a bit slower than eight threads.

Run it with several threads:

```sh
julia --project=. -t auto advanced/RICH/RICHOffloadGPU.jl
```

## Simplifications

- One particle per event, and no magnetic field.
- The reconstruction takes the track direction and momentum from the generator, as if a tracking
  system measured them perfectly. It also takes the radiator of each photon from the Monte Carlo
  truth.
- The photodetectors have no pixels, so the hit position is exact.
- The aerogel refractive index is constant, so there is no chromatic dispersion for the aerogel.
  Its resolution is therefore better than for a real aerogel.
- The light goes directly from the spherical mirrors to the photodetectors. The final LHCb RICH1
  design adds flat secondary mirrors.
- The run uses one thread.

## Files

- `RICH.lit`: Literate source of the example (source of truth)
- `RICH.jl`, `RICH.ipynb`: script and notebook generated from `RICH.lit`
- `RICHDetector.jl`: detector parameters, materials, optical properties and geometry
- `RICHOffload.lit`, `RICHOffload.jl`, `RICHOffload.ipynb`: variant with the optical photons offloaded
- `RICHOffloadGPU.lit`, `RICHOffloadGPU.jl`, `RICHOffloadGPU.ipynb`: variant with the photons offloaded to a GPU in batches
- `RICHTracer.jl`: stand-alone propagator of optical photons (module `RICHOptics`) using the Geom4hep
  navigation, and `exportGDML`/`RICHTracer` to build it for this detector
- `RICHTracerGPU.jl`: kernel, random number generator and batch transfer to run the propagator on a
  GPU (module `RICHGPU`)
- `Schematic-RICH-detector.png`: schematic of the detector

After you edit `RICH.lit`, regenerate the script and the notebook:

```julia
using Literate
Literate.script("advanced/RICH/RICH.lit", "advanced/RICH", name="RICH", keep_comments=false, documenter=true, credit=false)
Literate.notebook("advanced/RICH/RICH.lit", "advanced/RICH", name="RICH", execute=false, documenter=true)
```
