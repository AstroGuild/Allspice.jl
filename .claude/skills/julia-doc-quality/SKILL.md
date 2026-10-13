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

- One docstring on the generic function covers the methods with and without a context. Write `[ctx,]` in the signature (`name([ctx])` when the context is the only argument) and never show the `KernelSet` method. If methods in two files each have a docstring, keep one, with one signature line per call form.
- An exported type: the signature lines are the constructors a caller uses; list only the fields a caller reads; fold the docstrings of extra constructors into it.
- A family of calls that shares arguments and errors (`spkwNN`, the `gf` searches, calls on one handle type): state the shared part once, on one named member or on the handle type, and link to it from the others.
- State units on every physical quantity (km, km/s, rad, seconds past J2000 TDB).
- State the differences from CSPICE that a caller can observe, in fixed words: "`nothing` takes the place of the `found` flag of CSPICE"; "counts from 1; CSPICE counts from 0"; "This call has no CSPICE counterpart."
- A call that writes a file: say when the bytes reach the disk (at the call or at the close) and whether a call that throws leaves the file unchanged.
- More than three error kinds: write "Throws `SpiceError`:" and a list, one line per kind. A flag or a `nothing` result belongs in the body, not under "Throws".
- `# Arguments` only when the signature and one sentence cannot carry the meaning.
- An example is wanted for a main entry point and for behaviour that is not obvious; a trivial accessor needs none.
- Never write "Contract:" in the docstring of an exported name. Preconditions become "Throws ... if"; postconditions become part of "Return ...".

## Facts are measured, not copied

The largest source of false statements in the old docstrings was "as `other`": a sibling call that was said to throw the same errors and did not. Before a docstring says what a call throws or returns:

- Provoke each error kind with a bad input; the kind often comes from a lower layer and differs between sibling calls.
- Never write "Throws as [`other`](@ref)" without a comparison of the two in the code.
- A kept docstring that states something false is corrected, private or not; report the correction.

## Examples

- `jldoctest` when the example needs no kernel file that it does not write itself and its output is stable. Three ways to need no kernel file: pure mathematics; a fresh `SpiceContext()` filled with `pdpool`, `pcpool`, `lmpool` or `boddef`; a file that the example writes into `mktempdir()` (end that line with `;` so that no path is printed). Never use the default context in a doctest, because all doctests share it.
- The doctest module has only `using Allspice`, so a displayed `SVector` prints with a long header. Show a comparison (`f(x) == [1, 2, 3]`, `≈` for rounded values), a tuple, or `round.(x; digits=6)`. Prefer inputs with exact results. A value that two code paths compute can differ in the last bit: compare with `≈`.
- Write the expected output from a scratch run. Never use `doctest(...; fix=true)`: it turns a failing `== ... true` into a passing `false`.
- An example that needs real kernels is a plain `julia` block without output values. Run it once: against real kernels if the machine has them (use bodies that the named kernel holds), or against `docs/fixture.jl` with `PROBE`, `EARTH` and `SUN` substituted.
- No `$` and no backslash in a docstring that is not `raw`; write units as `km³/s²`. Do not name a variable `first`, `last` or `error`.

## Private functions

A private function keeps a docstring only when it states something the code and the name do not: a non-obvious precondition, an invariant, a unit, an arithmetic rule of f2c, a known difference from CSPICE. "Follows the CSPICE routine X" is a reason only when the name does not already say it; where a file maps names by one rule (`_ekw_<name>` ports `zzek<name>`), state the rule once in the file header. A one-line docstring on a private constant becomes a trailing `#` comment or goes. Executable checks (`@assert`) stay. A name that is public but not exported (`Allspice.det`) gets the public form.

## Manual

Pages follow reader needs, not `src/` file names: a tutorial, task pages, explanation pages, and a curated reference built with `@docs` blocks grouped by subject. Add a page only when a reader task needs it. Examples run (`@example`, `jldoctest`) on the kernels that `docs/fixture.jl` writes; the repository holds no kernel file. Try each example in a scratch script before it goes into a page, and after the build read the outputs in `docs/build/` once: a result that surprises (a search that finds fewer events than expected) is a fault of the example or of the text beside it.

The reference pages list every exported name once (`checkdocs=:exports` fails the build otherwise); a new export goes into the page of its subject.

## Checks before the work is done

```sh
julia --project=docs docs/make.jl            # build, doctests, cross-references
julia --project=docs docs/doctest.jl         # doctests of docstrings only (fast)
```

A change of docstrings only must not change code. Prove it: parse the old and the new file with `Meta.parseall`, replace each `Core.@doc` call by its target, remove the `LineNumberNode`s, and require equality.

Do not pass `warnonly=true`, and do not accept a changed doctest output without reading why it changed. Then hand the changed pages to the `julia-doc-review` skill.
