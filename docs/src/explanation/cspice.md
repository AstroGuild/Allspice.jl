# Differences from CSPICE

The package keeps the names, the argument order, the units and the numerical results of CSPICE N0067. It departs from CSPICE where the C interface carries hidden state or output arguments, and this page lists each departure, so that a reader who knows CSPICE or SPICE.jl can port a program with no surprise.

## Call conventions

| CSPICE | This package |
|---|---|
| Results come back through output arguments | A call returns its results, several of them as a tuple: `state, lt = spkezr(...)` |
| Vectors and matrices are C arrays | Results are static arrays (`SVector`, `SMatrix`) of StaticArrays.jl; inputs can be any `AbstractVector` or `AbstractMatrix` |
| An error sets a flag that the program must test | An error throws [`SpiceError`](@ref), whose field `kind` is a symbol and whose field `msg` is the text |
| A "found" flag tells that a lookup failed | The call returns `nothing`, as [`bodn2c`](@ref) does for an unknown name, or throws where a value is required |
| Indices count from 0 | Indices count from 1: the first interval of a window, the first plate of a DSK, the first row of an EK |
| A file handle is an integer, and several calls work on a hidden "current" file | A handle is an object ([`DAFHandle`](@ref), [`DASHandle`](@ref)) that the caller passes to each call |
| Loaded kernels are one table for the whole process | Loaded kernels belong to a [`SpiceContext`](@ref); see [Contexts and threads](@ref) |
| Output strings need a length argument | A call returns a `String` |
| A cell or a window has a size that the caller declares | A [`Window`](@ref) grows as needed |

The kind of an error is the short message of CSPICE in lower case where the same condition exists there, so `SPICE(INVALIDDEGREE)` becomes `:invaliddegree`.

## No state between calls

CSPICE keeps buffers and saved search results between calls, for speed, and in a few places an answer depends on the calls that came before. This package keeps no such state: a query is a function of the loaded kernels and of its arguments. Three visible consequences follow.

- **Threads.** Queries from many threads need no lock, and a kernel can be loaded while other threads query.
- **Settings.** The defaults of [`timdef`](@ref) and [`tsetyr`](@ref), the tolerance of [`gfstol`](@ref) and the settings of [`dskstl`](@ref) belong to a context, not to the process.
- **Repeatable answers.** Where the answer of CSPICE depends on an earlier call (for example `ckgp` directly after `ckgpav` on overlapping segments, or the dates in a text kernel loaded after certain `str2et` calls), this package returns the answer that CSPICE gives in a new process.

## Numerical results

Each routine follows the arithmetic of its CSPICE counterpart, operation by operation, including the order of the terms of a sum; fused multiply-add is not used. For the same kernels and the same `Float64` input, most results are equal to those of CSPICE bit for bit, and the others differ by rounding. [Coverage and validation](@ref) gives the measured differences.

Binary kernels that the writers produce are equal byte for byte to the files that the CSPICE writers produce from the same input.

## What this package adds

| Addition | Calls |
|---|---|
| Derivatives of any order | [`spkder`](@ref), [`pxder`](@ref) |
| Dual-number epochs | all queries; see [Time derivatives](@ref) |
| Tables for inner loops and GPUs | [`ephemeris_table`](@ref), [`rotation_table`](@ref) |
| Aberration correction as a value, parsed once | [`Aberration`](@ref) |
| Scratch files in memory | [`dasops`](@ref), [`ekops`](@ref) |

## What is absent

The toolkit utilities that are programs (`brief`, `mkspk`, `spkmerge` and the others) are not part of the package. The geometry finder calls [`gfsubc`](@ref), [`gfsntc`](@ref) and [`gfilum`](@ref) accept ellipsoids and not DSK shapes, as in CSPICE.
