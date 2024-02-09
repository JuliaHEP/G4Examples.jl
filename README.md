# Geant4.jl Examples Repository

Collection of examples to exercise the available functionality of the Julia bindings to Geant4. The idea was to separate this repository from the Julia package [Geant4.jl](https://github.com/JuliaHEP/Geant4.jl) since the examples carry additional dependencies to packages that provide visualization and data analysis functionality.

The examples are organised in 3 sections following more or less the organisation of the Geant4 repository: basic, extended and advanced. As soon new functionality will be added to Geant4.jl, new examples will be added in this repository.

## Basic
### B1/B1
This is most basic example using a more Geant4 native interface. It is creating a geometry with two volumes and standard physcis.
- **B1.jl** - julia script to run the example. Detector description in a separate file `DetectorB1.jl`
- **B1.ipynb** - the same example in a Jupyter notebook format

### B2/B2a
Basic example of a calorimeter detector using a sensitive detector to collect 'hits'. 
- **B2a.jl** - julia script to run the example. Detector description in a separate file `DetectorB2a.jl`
- **B2aVis.jl** - B2a detector with visualisation
- **B2a.ipynb** - the same example in a Jupyter notebook format

## B3/B3a
The example simulates schematically a Positron Emitted Tomography system. It uses a custom primary particle generator, a custom physics list and makes use of two sensitive detectors.
- **B3a.jl** -  julia script to run the example. Detector description in a separate file `DetectorB3.jl`
- **B3b.jl** - Same with different user actions
- **B3a.ipynb** - the same example in a Jupyter notebook format, with some more detailed explanations 

## Extended

### Solids
- **Solids.ipynb** - Notebook to visualize all implemented solid types

### TestEm3
This example comes from *extended/electromagnetic/TestEm3* example. It is fully based on user actions. It comes with additional packages such as FHist.jl and Plots.jl to provide some data analysis and visualization ogf results. 
- **TestEm3.jl** - julia script to run the example
- **TestEm3-gdml.jl** - the same example but reading the geometry from a GDML file instead of constructing it
- **TestEm3.ipynb** -  the same example in a Jupyter notebook format

### RE03
Example using the Geant4 built-in scoring mechanism
- **RE03.jl** - julia script to setup the scoring and run something simple. The more complete example using scoring if the WaterPhantom example in the advanced section.

### GPS
Example using the General Particle Source primary particle generator class `G4GeneralParticleSource`
- **GPS.jl** - With a simple detector geometry, a GPS is setup, which is then validated with a set of histograms.
- **GPS.ipynb** - same example in form of a note book and more test cases of GPS setups

## Advanced
### WaterPhantom
Example in a notebook format similar to RE03 but with different primary particle generator (MedicalBeam) and using the scoring mechanism. Plots are produced after each run.
- **WaterPhantom.ipynb** - notebook format

### HBC30
Example of a bubble chamber in which we display the particle tracks for an event that passes the trigger. Thi example introduces a magnetic field.
- **HBC30.jl** - script
- **HBC30.ipynb** - in a notebook format

### Scintillation
Example from the original at https://github.com/settwi/g4-basic-scintillation and adapted to the G4JL interface. Introduces optical photons and a custom physics list. It produces histograms as a result.
- **Scintillation.jl** - Main script to run the example. User actions, detector and material description and parameters in separate files.


