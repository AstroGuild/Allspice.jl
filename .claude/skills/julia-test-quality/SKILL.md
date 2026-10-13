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
- A random sample larger than what first exposes each branch. The large sample belongs to the exhaustive tier (`nsample(fast, full)`, `ALLSPICE_TEST_FULL=1`). When a smaller sample loses a branch, do not restore the sample: add a fixed case that reaches the branch, with a comment on how rare it is for random input. Do not shrink a sample that is the only oracle for a mode.
- A check that does not depend on the loop variable, repeated for each epoch or each sample (the forms of a call, a result type, `isbits`).
- One `@test` for each sample where one accumulated flag or one maximum gives the same failure; print the first failing samples instead.
- A test of a test fixture (a helper compared with the package through the same oracle).
- The same comparison repeated through a second entry point that only forwards. An entry point forwards if its check lives in a shared routine; keep one case per shared routine and one per entry point, and every entry point whose check takes an argument of its own.
- A syntax or message layer tested again in the group that consumes it. One group owns the layer; the consumer keeps the cases that need its own state.
- Behaviour of Base, StaticArrays or SPICE.jl.
- Internal layout that can change without a change of the public result.
- The exact text of an error message, except where the package states that the text equals CSPICE.

## Rules of this repository

- Framework: `Test` with the group files `test/test_<group>.jl`; no test-item migration, no Supposition.
- The oracle is CSPICE. A tolerance is the measured rounding difference with a stated reason, never a value chosen so that the test passes.
- Random input uses the seeded generator of the suite (`ALLSPICE_TEST_SEED`).
- CSPICE keeps state between calls; a test that calls it must reset what it changed (`tsetyr`, `tparch`, loaded kernels) in a `finally` block.
- A reader loop in a concurrency test calls `yield()`.
- A count threshold that goes with a random sample (`nhit > 800`) is a fraction of the sample size, so that it holds in both tiers.
- The cases of a fixture share one random stream: add a new case after the existing ones, so that the earlier cases keep their data.
- Each group has a time budget; a group that grows states its new time. The default suite stays under 120 s in one process.
- User-visible behaviour has a Given/When/Then scenario in `test/features/`; remove or change the scenario with the test.

## Review trigger

More than about 12 new scenarios or 150 new lines for one feature: stop, justify each group of scenarios in one line, and merge what protects the same behaviour.

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
