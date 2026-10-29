# CLAUDE.md

Examples for [Geant4.jl](https://github.com/JuliaHEP/Geant4.jl), the Julia bindings to Geant4. They live outside Geant4.jl because they pull in extra visualization and analysis packages (Makie, Plots, FHist, PYTHIA8, …). They are also the source of the example pages in the Geant4.jl documentation.

## Layout

- `basic/`, `extended/`, `advanced/`: examples grouped like the Geant4 C++ distribution (B1–B4, TestEm3, RE03, GPS, Solids, field02, WaterPhantom, HBC30, Scintillation, AlephTPC, CustomLib, EmbedJulia). `extended/OpNovice/` is a work-in-progress port of opNovice2 (untracked, only a stub so far).
- `src/G4Examples.jl`: an empty module. It exists only so the repo is a Julia package that `Pkg.test` can run.
- `src/literate.jl`: an older helper that converts `<name>_lit.jl` into a notebook and script with Literate.jl.
- `test/runtests.jl`: runs each example script as a subprocess and checks the exit code. Only the `basic/` examples are enabled. The extended and advanced ones are commented out because they are slow or have heavy dependencies.
- `.github/workflows/ci.yaml`: `julia-runtest` on Julia 1.10 and the latest 1.x, on Ubuntu, under `xvfb-run` (needed for GLMakie).

## Files in an example directory

Each example is written once in Literate.jl format and then generated into other formats:

| File | Role |
|---|---|
| `X.lit` | **Current source of truth.** Literate input with `#md`, `#nb` and `#jl` filter lines. Geant4.jl's docs build reads it (generated `.md` files have `EditURL = "X.lit"`). |
| `X_lit.jl` | Older Literate source, still present in some extended/advanced examples. It may be out of date compared with `X.lit`. |
| `X.jl` | Generated plain script, with comments stripped. This is what the tests run. |
| `X.ipynb` | Generated Jupyter notebook. |
| `X.md` | Generated Documenter page, when present. |
| `XDetector.jl` | Detector construction pulled out so several scripts can `include` it. |
| `XVis.jl`, `XVisSettings.jl` | Event-display variant, plus a NamedTuple of display settings (colors, axes, trajectories) passed to `G4JLEventDisplay(path)`. |

When you change an example, edit the `.lit` file and keep the generated `.jl`/`.ipynb` in sync with it. If you only edit the `.jl`, the next regeneration will overwrite your change. Commits such as "Re-generated basic examples from literate files in Geant4.jl" show that regeneration is done from the Geant4.jl docs tooling.

## Common Geant4.jl idioms in the examples

- `using Geant4, Geant4.SystemOfUnits`. Units are written as multiplications such as `20cm` or `50MeV`.
- Geometry: a construct function `f(det)::CxxPtr{G4VPhysicalVolume}` using `G4NistManager!Instance()`, `FindOrBuildMaterial`, `G4Box`/`G4Tubs`, `G4LogicalVolume` and `G4PVPlacement`. Parameterized detectors are `@with_kw mutable struct ... <: G4JLDetector` (from Parameters.jl) with `Geant4.getConstructor(::MyDet) = f`.
- Application: `G4JLApplication(; detector, generator, physics_type, field, nthreads, simdata, sdetectors, scorers, *_action...)`, then `configure(app)`, `initialize(app)`, `beamOn(app, n)`. Geant4 UI commands use the `` ui`/cmd ...` `` string macro.
- User data is a `mutable struct <: G4JLSimulationData`. With multithreading it needs `Base.merge!` and `Base.empty!`. Hits are collected with `G4JLSensitiveDetector(name, data; processhits_method=...)`. Built-in scoring uses `G4JLScoringMesh`.
- Primaries come from `G4JLGunGenerator`, `G4JLGeneralParticleSource`, or custom generators.
- Visualization: loading `CairoMakie` or `GLMakie` activates Geant4.jl's G4Vis extension. Then `draw(pv[])` draws the geometry, or `G4JLEventDisplay(settingsfile)` is passed as `evtdisplayaction` for event displays.
- Native C++ escape hatches: `advanced/CustomLib` builds a user C++ library and calls it, and `advanced/EmbedJulia` embeds Julia in a C++ Geant4 app.

## Running

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. basic/B1/B1.jl          # run one example
julia --project=. -e 'using Pkg; Pkg.test()'  # run the enabled examples
```

The tests `cd` to the repo root and call `julia <path>` without `--project`. They therefore depend on the active environment that `Pkg.test` sets up. Visualization examples need a display, or `xvfb-run` on headless Linux.

`Manifest.toml` is gitignored. The local one may `dev` Geant4.jl from a sibling checkout (`../Geant4.jl`), so local behavior can differ from the registered Geant4 version (compat `0.3`, Geant4_jll `11.2`).
