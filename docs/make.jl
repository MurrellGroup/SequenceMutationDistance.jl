using Documenter, SMD
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/MurrellGroup/SMD.jl.git",
    julia  = "nightly",
    osname = "osx"
makedocs()
