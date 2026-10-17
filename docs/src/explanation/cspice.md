# Differences from CSPICE

The package keeps the names, the argument order, the units and the numerical results of CSPICE N0067. It departs from CSPICE where the C interface carries hidden state or output arguments, and this page lists each departure, so that a reader who knows CSPICE or SPICE.jl can port a program with no surprise.

## Call conventions

| CSPICE | This package |
|---|---|
| A name ends in `_c` | The name has no suffix: `spkezr_c` is `spkezr` |
| A call fills or changes an argument | A call that changes its argument has a name that ends in `!`, as Julia code does: `wninsd!` |
| Results come back through output arguments | A call returns its results, several of them as a tuple: `state, lt = spkezr(...)` |
| Vectors and matrices are C arrays | Results are static arrays (`SVector`, `SMatrix`) of StaticArrays.jl; inputs can be any `AbstractVector` or `AbstractMatrix` |
| An error sets a flag that the program must test | An error throws a [`SpiceError`](@ref), whose field `kind` is the short message of CSPICE as a symbol and whose field `msg` is the text; missing SPK data throw a subtype with the bodies and the epoch |
| A "found" flag tells that a lookup failed | A lookup returns `nothing` ([`bodn2c`](@ref) for an unknown name); a call that must return a value throws ([`gdpool`](@ref) for a variable that is not loaded) |
| Indices count from 0 | Indices count from 1: the first interval of a window, the first plate of a DSK, the first row of an EK |
| A file handle is an integer, and several calls work on a hidden "current" file | A handle is an object ([`DAFHandle`](@ref), [`DASHandle`](@ref)) that the caller passes to each call |
| Loaded kernels are one table for the whole process | Loaded kernels belong to a [`SpiceContext`](@ref); see [Contexts and threads](@ref) |
| Output strings need a length argument | A call returns a `String` |
| A cell or a window has a size that the caller declares | A [`Window`](@ref) grows as needed |

The kind of an error is the short message of CSPICE in lower case where the same condition exists there, so `SPICE(INVALIDDEGREE)` becomes `:invaliddegree`.

## No state between calls

CSPICE keeps buffers and saved search results between calls, for speed. This package keeps no such state: a query is a function of the loaded kernels and of its arguments. Three visible consequences follow.

- **Threads.** Queries from many threads need no lock, and a kernel can be loaded while other threads query; see [Contexts and threads](@ref).
- **Settings.** The defaults of [`timdef`](@ref), [`tsetyr`](@ref) and [`tparch`](@ref), the tolerance of [`gfstol`](@ref) and the settings of [`dskstl`](@ref) belong to a context, not to the process.
- **Repeatable answers.** Where the answer of CSPICE depends on an earlier call (for example `ckgp` directly after `ckgpav` on overlapping segments, the dates in a text kernel loaded after certain `str2et` calls, or the four padding bytes of each segment summary of a binary PCK, which CSPICE fills from the last SPK or CK summary it handled), this package returns the answer that CSPICE gives in a new process.

## Numerical results

Each routine follows the arithmetic of its CSPICE counterpart, operation by operation, including the order of the terms of a sum; fused multiply-add is not used. For the same kernels and the same `Float64` input, a result is equal to that of CSPICE bit for bit or differs by rounding; [Coverage and validation](@ref) says how this was checked and gives measured differences.

## What this package adds

| Addition | Calls |
|---|---|
| Derivatives of any order | [`spkder`](@ref), [`pxder`](@ref) |
| Dual-number epochs | all queries; see [Time derivatives](@ref) |
| Tables for inner loops and GPUs | [`ephemeris_table`](@ref), [`rotation_table`](@ref) |
| Aberration correction as a value, parsed once | [`Aberration`](@ref) |
| Scratch files in memory | [`dasops`](@ref), [`ekops`](@ref) |

## What is absent

The toolkit utilities that are programs (`brief`, `mkspk`, `spkmerge` and the others) are not part of the package. [Coverage and validation](@ref) lists the areas of the library that exist.
