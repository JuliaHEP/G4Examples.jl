# Geant4.jl Examples

Examples for [Geant4.jl](https://github.com/JuliaHEP/Geant4.jl), the Julia bindings to the
[Geant4](https://geant4.web.cern.ch) simulation toolkit. They range from simple geometries to
optical photons, scoring, magnetic fields, event displays and offloading to GPUs.

Every page in this documentation is generated from the source of an example in the
[G4Examples.jl](https://github.com/JuliaHEP/G4Examples.jl) repository, and the example is run when
the documentation is built. The pages show the code together with its output and plots, produced
with the current versions of Geant4.jl and Geant4. Each page has links to download the example as a
Jupyter notebook or as a plain Julia script.

## Running the examples

Clone the repository and install the dependencies:

```sh
git clone https://github.com/JuliaHEP/G4Examples.jl.git
cd G4Examples.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The first installation downloads the Geant4 libraries and data files, so it takes a while. Then run
an example as a script, or open its notebook (`.ipynb`) in Jupyter or VS Code:

```sh
julia --project=. basic/B1/B1.jl
```

The scripts with visualization use GLMakie and open a window. On a Linux machine without a
display, run them with `xvfb-run`. The notebooks and these pages use CairoMakie.

## Examples

The examples are grouped like the examples of the Geant4 C++ distribution.

```@contents
Pages = Main.examples_mds
Depth = 1:1
```
