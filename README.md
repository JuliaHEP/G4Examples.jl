# Geant4.jl Examples Repository

Collection of examples to exercise the available functionality of the Julia bindings to Geant4. The idea was to separate this repository from the Julia package [Geant4.jl](https://github.com/JuliaHEP/Geant4.jl) since the examples carry additional dependencies to packages that provide visualization and data analysis functionality.

The examples are organised in 3 sections following more or less to the organisation of the Geant4 repository: basic, extended and advanced. As soon new functionality will be added to Geant4.jl, new examples will be added in this repository.

All examples are provided as *Julia scripts* and *Jupyter notebooks* at the same time. This is achieved by using the package [`Literate.jl`](https://github.com/fredrikekre/Literate.jl) that uses a single input file to produce the different output formats. The script in `src/literate.jl` does the conversion.

Finally, the examples in this repository as also used to generate the documentation of [`Geant4.jl`](https://juliahep.github.io/Geant4.jl/stable/)

## Basic
### B1/B1
This is most basic example using a more Geant4 native interface. It is creating a geometry with two volumes and standard physcis.

### B2/B2a
Basic example of a calorimeter detector using a sensitive detector to collect 'hits'. 

## B3/B3a
The example simulates schematically a Positron Emitted Tomography system. It uses a custom primary particle generator, a custom physics list and makes use of two sensitive detectors.

## Extended

### Solids
Notebook and script to exercise and visualize all implemented solid types.

### TestEm3
This example comes from *extended/electromagnetic/TestEm3* example. It is fully based on user actions. It comes with additional packages such as FHist.jl and Plots.jl to provide some data analysis and visualization ogf results. 

### RE03
Example using the Geant4 built-in scoring mechanism. A more complete example using scoring is the `WaterPhantom` example in the advanced section.

### GPS
Example using the General Particle Source primary particle generator class `G4GeneralParticleSource`

## Advanced
### WaterPhantom
Example in a notebook format similar to RE03 but with different primary particle generator (MedicalBeam) and using the scoring mechanism. Plots are produced after each run.

### HBC30
Example of a bubble chamber in which we display the particle tracks for an event that passes the trigger. This example introduces a magnetic field.

### Scintillation
Example from the original at https://github.com/settwi/g4-basic-scintillation and adapted to the Geant4.jl interface. Introduces optical photons and a custom physics list. It produces histograms as a result.

### AlephTPC
Example that shows how easy is to integrate several packages in Julia. The example uses the package `PYTHIA8` to generate LEP collisions of e+ e- which are then input as primary particles to a `Geant4` simulation. Only the event display is exercised of the simulation. 
