---
name: julia-test-quality
description: Use when adding, changing, reviewing or trimming tests of this Julia package, and when a fault is fixed and needs a regression test.
---

# Julia test quality

The goal is the smallest set of tests that detects a real change of behaviour. The count of tests, the count of assertions and the coverage percentage are not goals; report them, do not optimise them.

## Before writing tests: the behaviour matrix

Write one row for each scenario before any code, and delete duplicate rows first.

| Behaviour | Regression it detects | Why distinct from other rows | Layer |
|---|---|---|---|

Layers in this repository: differential (result compared with CSPICE through SPICE.jl), analytic (closed-form answer), invariant (round trip, identity), error (kind of `SpiceError`), concurrency, allocation.

## What earns a test

- One representative case for each class of behaviour: each data type, each branch of a method string, each aberration correction family, not each combination of them.
- Boundaries: first and last epoch of coverage, empty input, degenerate geometry, the largest and smallest accepted value.
- Each error kind that a caller can meet, once.
- Each fixed fault, as one focused regression test with a comment that names the fault.
- Differences from CSPICE that are deliberate, pinned by a test that says so.

## What does not

- A Cartesian product of arguments where the factors are independent in the code.
- A random sample larger than what first exposes each branch. Default tier: the smallest sample that reaches every branch; the large sample belongs to the exhaustive tier (`ALLSPICE_TEST_FULL=1`).
- The same comparison repeated through a second entry point that only forwards.
- Behaviour of Base, StaticArrays or SPICE.jl.
- Internal layout that can change without a change of the public result.
- The exact text of an error message, except where the package states that the text equals CSPICE.

## Rules of this repository

- Framework: `Test` with the group files `test/test_<group>.jl`; no test-item migration, no Supposition.
- The oracle is CSPICE. A tolerance is the measured rounding difference with a stated reason, never a value chosen so that the test passes.
- Random input uses the seeded generator of the suite (`ALLSPICE_TEST_SEED`).
- CSPICE keeps state between calls; a test that calls it must reset what it changed (`tsetyr`, `tparch`, loaded kernels) in a `finally` block.
- A reader loop in a concurrency test calls `yield()`.
- Each group has a time budget; a group that grows states its new time. The default suite stays under 120 s in one process.
- User-visible behaviour has a Given/When/Then scenario in `test/features/`; remove or change the scenario with the test.

## Review trigger

More than about 12 new scenarios or 150 new lines for one feature: stop, justify each group of scenarios in one line, and merge what protects the same behaviour.

## Trimming existing tests

For each group: list the behaviours it protects, mark cases that protect the same one, keep the cheapest representative, and move the large sample to the exhaustive tier where it has value there. Report per group: behaviours kept, cases removed and why, assertions and seconds before and after, and line coverage of the source files that the group owns. Coverage may fall only on lines that no remaining behaviour needs; name those lines.

## Commands

```sh
julia --startup-file=no --threads=8 --project=test test/runtests.jl [group ...]
julia --startup-file=no --threads=8 test/parallel.jl [group ...]
julia --startup-file=no --threads=8 --code-coverage=@src --project=test test/runtests.jl
julia --startup-file=no test/coverage_summary.jl
```
