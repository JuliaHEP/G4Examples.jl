#---Run one example and generate its documentation page with Literate.jl
#   Usage: julia --project=docs docs/literate.jl <path to X.lit> <output directory>
#   It is called by make.jl in a separate process for each example, because Geant4 can only be
#   initialized once per process.
using Literate

lit, outdir = ARGS
name = splitext(basename(lit))[1]
println("====> Running $name")

# The input is the original .lit file, so the EditURL of the page points to it in the repository.
# The code runs in the output directory, which also has a copy of the files that the example includes.
Literate.markdown(lit, outdir; name = name, execute = true, documenter = true, credit = true)
