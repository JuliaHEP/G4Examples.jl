#---Build the documentation of the examples
#   Usage: julia --project=docs docs/make.jl [name ...]
#   Each X.lit file in basic/, extended/ and advanced/ becomes a page with the code and its output,
#   produced by running the example. With names (e.g. `B1 RICH`) only those examples are run, for
#   quick checks; the documentation is then not deployed.
using Documenter

const rootdir     = dirname(@__DIR__)
const examplesdir = joinpath(@__DIR__, "src", "examples")
const sections    = ["basic" => "Basic", "extended" => "Extended", "advanced" => "Advanced"]
const selected    = ARGS

#---Copy the example directories into docs/src/examples--------------------------------------------
#   The pages, the files that the examples include, their images, notebooks and scripts all end up
#   in one flat directory, so that relative links and @ref links between them work.
rm(examplesdir; force = true, recursive = true)
mkpath(examplesdir)
litfiles = Dict(s => String[] for (s, _) in sections)
sources  = Dict{String,String}()                  # file name => original path
for (section, _) in sections, dir in readdir(joinpath(rootdir, section); join = true)
    isdir(dir) || continue
    for f in readdir(dir)
        (startswith(f, ".") || f == "README.md") && continue
        cp(joinpath(dir, f), joinpath(examplesdir, f); force = true)
        sources[f] = joinpath(dir, f)
        endswith(f, ".lit") && push!(litfiles[section], joinpath(dir, f))
    end
end
pagename(lit) = splitext(basename(lit))[1]
isempty(selected) || foreach(s -> filter!(l -> pagename(l) in selected, litfiles[s]), keys(litfiles))

#---Run each example with Literate.jl in its own Julia process--------------------------------------
const julia = `$(Base.julia_cmd()) --project=$(@__DIR__) --threads=auto`
section_mds = map(sections) do (section, title)
    mds = map(litfiles[section]) do lit
        run(`$julia $(joinpath(@__DIR__, "literate.jl")) $lit $examplesdir`)
        "examples/$(pagename(lit)).md"
    end
    title => mds
end
filter!(p -> !isempty(p.second), section_mds)
const examples_mds = vcat(last.(section_mds)...)       # used by the @contents block in index.md

#---Pages with the source of the files that the examples include------------------------------------
#   Any `[file.jl](@ref)` or `[file.cpp](@ref)` link in the selected examples gets a page.
function create_extra(fname)
    name, ext = splitext(fname)
    src = sources[fname]
    open(joinpath(examplesdir, "$name.md"), "w") do io
        println(io, "```@meta\nEditURL = \"$(relpath(src, examplesdir))\"\n```")
        println(io, "# $fname\n")
        println(io, "```", ext == ".cpp" ? "c++" : ext == ".jl" ? "julia" : "")
        write(io, read(src, String))
        println(io, "```")
    end
    return "examples/$name.md"
end
refs = Set{String}()
for (s, _) in sections, lit in litfiles[s], m in eachmatch(r"\[`?([\w-]+\.(?:jl|cpp))`?\]\(@ref\)", read(lit, String))
    push!(refs, m[1])
end
extra_mds = map(create_extra, sort!(collect(refs)))

#---Build the documentation------------------------------------------------------------------------
makedocs(;
    sitename = "G4Examples.jl",
    authors = "Pere Mato",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://juliahep.github.io/G4Examples.jl",
        edit_link = "main",
        size_threshold = 8_000_000,
        size_threshold_warn = 4_000_000,
        search_size_threshold_warn = 2_000_000,
    ),
    pages = [
        "Introduction" => "index.md",
        section_mds...,
        hide.(extra_mds)...,
    ],
    repo = Remotes.GitHub("JuliaHEP", "G4Examples.jl"),
    checkdocs = :none,
    warnonly = !isempty(selected),       # a partial build has dangling links to other examples
)

#---Deploy the documentation-----------------------------------------------------------------------
isempty(selected) && deploydocs(;
    repo = "github.com/JuliaHEP/G4Examples.jl",
    devbranch = "main",
    push_preview = true,
)
