# Run the doctests of the docstrings only, without a build of the manual:
#     julia --project=docs docs/doctest.jl
using Documenter, Allspice
DocMeta.setdocmeta!(Allspice, :DocTestSetup, :(using Allspice); recursive=true)
doctest(Allspice; manual=false)
