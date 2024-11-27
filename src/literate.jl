import Pkg; Pkg.add("Literate")

using Literate
name = ARGS[1]

file = joinpath(@__DIR__, "..", "$(name)_lit.jl")
isfile(file) || error("File $file does not exist")
destdir =  dirname(@__DIR__)

Literate.notebook(file, destdir, name = name, execute = false, documenter = true, credit = true)
Literate.script(file, destdir, name = name, keep_comments = false, documenter = true, credit = false)

