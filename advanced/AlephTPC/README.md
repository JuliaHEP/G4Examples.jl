# ALEPH TPC example

This example shows how easy is to integrate several packages in Julia. The example uses the package `PYTHIA8` to generate LEP collisions of e+ e- together with the `Geant4` package. Only the event display is exercised. 

Figures like this should be produced:
![TPC event](TPCSim-22.png)

The example is written as a `Literate.jl` script. A single source is used to produce a markdown, script or notebook output. There are the commands to generate the output. 

```julia
using Literate
src = "TPCSim_lit.jl"
Literate.markdown(src, pwd(), name = "TPCSim", execute=true)
Literate.notebook(src, pwd(), name = "TPCSim")
Literate.script(src, pwd(), name = "TPCSim")
```



