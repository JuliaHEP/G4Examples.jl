#---Generate the X.jl script and X.ipynb notebook from every X.lit file in the repository
#   Usage: julia --project=. src/generate.jl [name ...]
#   With no arguments all .lit files are processed, otherwise only the ones with the given names.
using Literate

const rootdir = dirname(@__DIR__)

function litfiles(root)
    files = String[]
    for (dir, _, names) in walkdir(root)
        any(startswith("."), splitpath(relpath(dir, root))) && continue
        append!(files, joinpath(dir, f) for f in names if endswith(f, ".lit"))
    end
    return files
end

function generate(file)
    dir, name = dirname(file), splitext(basename(file))[1]
    @info "Generating $name.jl and $name.ipynb in $(relpath(dir, rootdir))"
    Literate.script(file, dir, name = name, keep_comments = false, credit = false)
    Literate.notebook(file, dir, name = name, execute = false, credit = true)
end

files = litfiles(rootdir)
if !isempty(ARGS)
    files = filter(f -> splitext(basename(f))[1] in ARGS, files)
    isempty(files) && error("No .lit files found for $(join(ARGS, ", "))")
end
foreach(generate, files)
