# Allspice.jl

Native Julia port of CSPICE N0067. `CONTRACT.md` records the internal interfaces.

Use the project skills whenever you change an exported name, a docstring, the README, the manual or a test:

- `julia-doc-quality`: README, docstrings, Documenter manual.
- `julia-doc-review`: independent reader review after a documentation change (agent `julia-doc-reader`).
- `julia-test-quality`: test design, regression tests, trimming.

Commands:

```sh
julia --startup-file=no --threads=8 test/parallel.jl [group ...]     # tests, 8 processes
julia --project=docs docs/make.jl                                    # manual and doctests
```
