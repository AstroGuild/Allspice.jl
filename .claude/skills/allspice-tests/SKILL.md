---
name: allspice-tests
description: Use when adding, changing, reviewing or trimming tests of Allspice.jl, or when a fault is fixed and needs a regression test; the rules of this repository on top of the global julia-tests skill.
---

# Allspice test rules

The general rules (what earns a test, what does not, the review trigger) are in the global skill `julia-tests`. This file holds only what is particular to Allspice.jl.

## Before writing tests: the behaviour matrix

Write one row for each scenario before any code, and delete duplicate rows first.

| Behaviour | Regression it detects | Why distinct from other rows | Layer |
|---|---|---|---|

Layers in this repository: differential (result compared with CSPICE through SPICE.jl), analytic (closed-form answer), invariant (round trip, identity), error (kind of `SpiceError`), concurrency, allocation.

## Rules of this repository

- Framework: `Test` with the group files `test/test_<group>.jl`; no test-item migration, no Supposition.
- The oracle is CSPICE. A tolerance is the measured rounding difference with a stated reason, never a value chosen so that the test passes.
- Random input uses the seeded generator of the suite (`ALLSPICE_TEST_SEED`).
- CSPICE keeps state between calls; a test that calls it must reset what it changed (`tsetyr`, `tparch`, loaded kernels) in a `finally` block.
- A reader loop in a concurrency test calls `yield()`.
- A count threshold that goes with a random sample (`nhit > 800`) is a fraction of the sample size, so that it holds in both tiers.
- The cases of a fixture share one random stream: add a new case after the existing ones, so that the earlier cases keep their data.
- Each group has a time budget; a group that grows states its new time. The default suite stays under 120 s in one process.

## Trimming existing tests

First measure where the time goes: time each test set and separate the body from the compilation before it. In this repository about 85 % of the time of a group is compilation of the test code and of the package, so a smaller sample removes assertions and not seconds. Trim for regression value per line; to save seconds, reduce the number of distinct compiled types (nesting depth of dual numbers, `Val` orders, closures).

For each group: list the behaviours it protects, mark cases that protect the same one, keep the cheapest representative, and move the large sample to the exhaustive tier where it has value there. Report per group: behaviours kept, cases removed and why, and assertions before and after.

Coverage: compare the set of lines hit before and after, merged over all groups, not the percentage of one group. Julia counts only lines of compiled methods, so the denominator changes between runs and a function that no test calls does not lower the percentage. A line that one group loses and its owner group still runs is not a loss. Find missing tests from the list of exported names that no test calls and from the uncovered lines read as candidate behaviours.

Timing: paired runs without coverage, original and changed alternately, at least twice. One before/after pair on a shared machine is not evidence.

## Commands

```sh
julia --startup-file=no --threads=8 --project=test test/runtests.jl [group ...]
julia --startup-file=no --threads=8 test/parallel.jl [group ...]
julia --startup-file=no --threads=8 --code-coverage=@src --project=test test/runtests.jl
julia --startup-file=no test/coverage_summary.jl
```

Run one coverage measurement at a time in one tree: the summary deletes every `.cov` file under `src/`.
