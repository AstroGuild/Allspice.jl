# Build the manual: julia --project=docs docs/make.jl
using Documenter, Allspice

DocMeta.setdocmeta!(Allspice, :DocTestSetup, :(using Allspice); recursive=true)

makedocs(;
    sitename="Allspice.jl",
    modules=[Allspice],
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true",
        size_threshold=600_000, size_threshold_warn=300_000),
    pages=[
        "Home" => "index.md",
        "First session" => "tutorial.md",
        "Tasks" => [
            "howto/contexts.md",
            "howto/derivatives.md",
            "howto/tables.md",
            "howto/events.md",
            "howto/writing.md",
        ],
        "Explanation" => [
            "explanation/cspice.md",
            "explanation/coverage.md",
        ],
        "Reference" => [
            "reference/kernels.md",
            "reference/time.md",
            "reference/states.md",
            "reference/frames.md",
            "reference/geometry.md",
            "reference/events.md",
            "reference/math.md",
            "reference/files.md",
            "reference/ek.md",
        ],
    ],
    checkdocs=:exports,
)
