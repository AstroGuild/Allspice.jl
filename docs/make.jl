# Build the manual: julia --project=docs docs/make.jl
using Documenter, Allspice

DocMeta.setdocmeta!(Allspice, :DocTestSetup, :(using Allspice); recursive=true)

makedocs(;
    sitename="Allspice.jl",
    modules=[Allspice],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true"),
    pages=["Home" => "index.md"],
    checkdocs=:exports,
)
