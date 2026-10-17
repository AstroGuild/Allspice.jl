# Allspice.jl

Native Julia port of CSPICE N0067. `CONTRACT.md` records the internal interfaces.

Use the skills whenever you change an exported name, a docstring, the README, the manual or a test. The global skills hold the general rules; the project skills add the rules of this repository:

- `julia-docstrings`, `julia-manual`, `julia-readme`, with `allspice-docs`: docstrings, the Documenter manual, the README.
- `julia-doc-reader-review` (agent `julia-doc-reader`): independent reader review after a documentation change; `allspice-docs` says what the reader gets.
- `julia-tests`, with `allspice-tests`: test design, regression tests, trimming.

Commands:

```sh
julia --startup-file=no --threads=8 test/parallel.jl [group ...]     # tests, 8 processes
julia --project=docs docs/make.jl                                    # manual and doctests
```
