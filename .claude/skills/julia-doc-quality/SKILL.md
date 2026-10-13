---
name: julia-doc-quality
description: Use when writing or changing the README, a docstring, or a page of the Documenter manual of this Julia package, and when adding or renaming an exported name.
---

# Julia documentation quality

The reader is a competent Julia developer who knows the domain and has never seen this repository. Usefulness to that reader comes before volume: a text is better when the reader understands more and can do more after less reading.

## Where a fact lives

Each fact has one home, and other places link to it.

| Fact | Home |
|---|---|
| What the package is, who needs it, install, one worked example, where to go next | `README.md` |
| A task from start to end, a concept, a numerical assumption, a limit | manual page in `docs/src/` |
| What one call takes, returns, changes and throws | its docstring |
| Internal interfaces and invariants for maintainers | `CONTRACT.md` and private docstrings |

## Prose register

Two registers apply, chosen by what the reader does with the text.

- **Reference text** is scanned or followed step by step: a docstring, the install and example part of the README, the steps of a task page, a table cell. Write short, direct sentences with one fact or one instruction each (20 words at most). The first sentence of a docstring states what the call returns or does.
- **Explanatory text** is followed as an argument: the opening of the README, an explanation page, the paragraph of a task page that says why. Use Simplified Technical English with varied sentence length; join related facts with a conjunction, a colon or a semicolon. Three or more consecutive sentences under about 15 words is a defect there.

Both registers: use the standard Julia and domain terms (do not paraphrase "type stability", "tolerance", "light time"); define a term the reader may not know in one sentence where it first appears; write no sentence that states no fact ("powerful", "seamless", "robust", "leverages"). Measured values go in a table, and every number in the text is one that a command produced: never write an accuracy or a timing from memory.

## Docstring form for an exported name

```julia
"""
    name([ctx,] arg1, arg2; keyword=default) -> result

Return <what>, in <units and frame>, as `name_c` of CSPICE.

<Only what the caller must know: mutation, ordering, defaults, numerical limits,
edge cases. Leave out what the signature already says.>

Throws `SpiceError(:kind)` if <condition>.

# Examples
```jldoctest
julia> name(1.0, 2.0)
3.0
```

See also [`other`](@ref).
"""
```

- One docstring on the generic function covers the methods with and without a context; do not repeat it per method.
- State units on every physical quantity (km, km/s, rad, seconds past J2000 TDB).
- State the differences from CSPICE that a caller can observe: index base, result shape, an error in place of a flag.
- `# Arguments` only when the signature and one sentence cannot carry the meaning.
- An example is wanted for a main entry point and for behaviour that is not obvious; a trivial accessor needs none.
- `jldoctest` only when the example needs no kernel and its output is stable (no random values, no addresses, no timing). An example that needs kernels is a plain `julia` block without invented output; the manual holds the executed version.
- Never write "Contract:" in the docstring of an exported name. Preconditions become "Throws ... if"; postconditions become part of "Return ...".

## Private functions

A private function keeps a docstring only when it states something the code does not: a non-obvious precondition, an invariant, a unit, the CSPICE routine it follows, a known difference. Delete a docstring that only restates the name or the body of a short function. Executable checks (`@assert`) stay.

## Manual

Pages follow reader needs, not `src/` file names: a tutorial, task pages, explanation pages, and a curated reference built with `@docs` blocks grouped by subject. Add a page only when a reader task needs it. Examples run (`@example`, `jldoctest`) on the kernels that `docs/fixture.jl` writes; the repository holds no kernel file. Try each example in a scratch script before it goes into a page, and after the build read the outputs in `docs/build/` once: a result that surprises (a search that finds fewer events than expected) is a fault of the example or of the text beside it.

The reference pages list every exported name once (`checkdocs=:exports` fails the build otherwise); a new export goes into the page of its subject.

## Checks before the work is done

```sh
julia --project=docs docs/make.jl            # build, doctests, cross-references
julia --project=docs docs/doctest.jl         # doctests of docstrings only (fast)
```

Do not pass `warnonly=true`, and do not accept a changed doctest output without reading why it changed. Then hand the changed pages to the `julia-doc-review` skill.
