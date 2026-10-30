# Geant4.jl Examples

[![CI](https://github.com/JuliaHEP/G4Examples.jl/actions/workflows/ci.yaml/badge.svg)](https://github.com/JuliaHEP/G4Examples.jl/actions/workflows/ci.yaml)

Examples for [Geant4.jl](https://github.com/JuliaHEP/Geant4.jl), the Julia bindings to the
[Geant4](https://geant4.web.cern.ch) simulation toolkit. They exercise the functionality of the
bindings, from simple geometries to optical photons, scoring, magnetic fields and event displays.

The examples live in their own repository, outside Geant4.jl, because they depend on extra
packages for visualization and analysis (Makie, Plots, FHist, PYTHIA8, ...). They are also the
source of the example pages in the
[Geant4.jl documentation](https://juliahep.github.io/Geant4.jl/stable/).

## Getting started

Clone the repository and install the dependencies:

```sh
git clone https://github.com/JuliaHEP/G4Examples.jl.git
cd G4Examples.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The first installation downloads the Geant4 libraries and data files, so it takes a while.

Run an example as a script:

```sh
julia --project=. basic/B1/B1.jl
```

You can also open the notebook (`.ipynb`) of an example in Jupyter or VS Code.

The examples with visualization use GLMakie in the scripts and open a window. On a Linux machine
without a display, run them with `xvfb-run`. The notebooks use CairoMakie and show the figures
inline.

Run the test suite, which runs the examples and checks that they finish without errors:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

## How the examples are organised

The examples are grouped in three sections, like the examples of the Geant4 C++ distribution:
`basic`, `extended` and `advanced`. Each example is in its own directory and is written once, in
[Literate.jl](https://github.com/fredrikekre/Literate.jl) format. The other formats are generated
from it.

| File | Content |
|---|---|
| `X.lit` | Source of the example: code with Markdown comments. Edit this file. |
| `X.jl` | Plain Julia script, generated from `X.lit` |
| `X.ipynb` | Jupyter notebook, generated from `X.lit` |
| `XDetector.jl` | Detector construction, in a separate file when it is long or shared |
| `XVisSettings.jl` | Settings of the event display (colours, trajectories, ...) |

After you edit a `.lit` file, regenerate the script and the notebook:

```julia
using Literate
Literate.script("basic/B1/B1.lit", "basic/B1", name="B1", keep_comments=false, documenter=true, credit=false)
Literate.notebook("basic/B1/B1.lit", "basic/B1", name="B1", execute=false, documenter=true)
```

The tests regenerate the script from the `.lit` file before running it.

## Examples

### Basic

| Example | Description |
|---|---|
| [B1](basic/B1) | Simple geometry with a few volumes and standard physics. It uses an interface close to the native C++ one. `B1vis.jl` draws the geometry. |
| [B2](basic/B2) | Fixed-target tracker that collects hits with a sensitive detector. `B2aVis` adds an event display. |
| [B3](basic/B3) | Schematic Positron Emission Tomography (PET) system, with a custom primary generator, a custom physics list and two sensitive detectors. |
| [B4](basic/B4) | Sampling calorimeter with energy deposits accumulated in user actions. `B4aVis.jl` adds an event display. |

### Extended

| Example | Description |
|---|---|
| [Solids](extended/Solids) | Builds and draws all the solid types available in Geant4.jl. |
| [TestEm3](extended/TestEm3) | Sampling calorimeter from *extended/electromagnetic/TestEm3*, based on user actions, with histograms made with FHist and Plots. `TestEm3-gdml.jl` reads the geometry from a GDML file and `TestEm3Vis.jl` adds an event display. |
| [RE03](extended/RE03) | Geant4 built-in scoring, with a scoring mesh defined by a `G4JLScoringMesh`. |
| [GPS](extended/GPS) | Primary particles from the General Particle Source (`G4GeneralParticleSource`). |
| [field02](extended/field02) | Tracking in an electric field and field-dependent electromagnetic processes. |
| [OpNovice](extended/OpNovice) | Port of *extended/optical/OpNovice2*: optical properties and surfaces, Cerenkov and scintillation photons, Fresnel reflection, thin-film coatings and the optical boundary models. |

### Advanced

| Example | Description |
|---|---|
| [WaterPhantom](advanced/WaterPhantom) | Dose in a water phantom with scoring meshes and a medical beam generator. Plots are made after each run. |
| [HBC30](advanced/HBC30) | CERN 30 cm liquid hydrogen bubble chamber in a magnetic field, with an event display of the events that pass a simple trigger. |
| [Scintillation](advanced/Scintillation) | Scintillation detector adapted from [g4-basic-scintillation](https://github.com/settwi/g4-basic-scintillation). Optical photons, optical surfaces and a custom physics list. |
| [AlephTPC](advanced/AlephTPC) | LEP e⁺e⁻ collisions generated with PYTHIA8 and simulated in the ALEPH Time Projection Chamber, with an event display. |
| [RICH](advanced/RICH) | Ring Imaging Cherenkov detector with aerogel and C₄F₁₀ radiators. The Cherenkov angle is reconstructed from the photon hits and used to identify pions, kaons and protons. |
| [CustomLib](advanced/CustomLib) | Builds a user C++ library, here a custom solid, and calls it from Julia, for Geant4 functionality that is not wrapped. |
| [EmbedJulia](advanced/EmbedJulia) | The reverse: a C++ Geant4 application that embeds Julia and calls user actions written in Julia. |

## Versions

The examples use Geant4.jl `0.3` with Geant4 `11.2` or later (`Geant4_jll`). The tests run on
Julia 1.10 and the latest 1.x release.
