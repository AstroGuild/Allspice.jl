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

Of the 651 routines that the header of CSPICE N0067 declares, 162 have no counterpart of the same name. Each of them is a facility that Julia provides, a routine that NAIF has declared deprecated and replaced, or a part of the C interface that the call conventions above make unnecessary. The table gives the call to use in place of each group, so that a program written for CSPICE or SPICE.jl can be ported without a search.

| CSPICE routines | In this package |
|---|---|
| Error handling: `chkin`, `chkout`, `erract`, `errch`, `errdp`, `errint`, `errprt`, `failed`, `getmsg`, `qcktrc`, `reset`, `return`, `setmsg`, `sigerr`, `trcdep`, `trcnam`, `trcoff` and others | A call throws [`SpiceError`](@ref); catch it with `try` |
| Cells and sets: `insrtd`, `removd`, `elemi`, `union`, `inter`, `diff`, `sdiff`, `scard`, `size`, `valid` and others | `Vector` and `Set` of Julia; windows are a [`Window`](@ref) |
| Arrays: `bsrchd`, `isrchi`, `lstled`, `orderd`, `reordd`, `shelld`, `maxd`, `sumad`, `filld`, `moved`, `mequ` and others | `searchsortedfirst`, `findfirst`, `sortperm`, `sort!`, `maximum`, `sum`, `fill!`, `copy` |
| Strings and parsing: `ucase`, `lcase`, `cmprss`, `lparse`, `nextwd`, `prsdp`, `prsint`, `repmc`, `matchw`, `lx4num` and others | `uppercase`, `lowercase`, `split`, `parse`, `replace`, `occursin`, `Regex` |
| Constants: `pi`, `dpmax`, `dpmin`, `intmax`, `intmin` | `π`, `floatmax`, `-floatmax`, `typemax(Int32)`, `typemin(Int32)` |
| Files and the command line: `exists`, `rdtext`, `prompt`, `getcml`, `putcml`, `ftncls` | `isfile`, `eachline`, `readline`, `ARGS`, `close` |
| Deprecated by NAIF: `bodvar`, `illum`, `spkapp`, `srfxpt`, `subpt`, `subsol`, `dafrda` | [`bodvrd`](@ref) or [`bodvcd`](@ref), [`ilumin`](@ref), [`spkaps`](@ref), [`sincpt`](@ref), [`subpnt`](@ref), [`subslr`](@ref), [`dafgda`](@ref) |
| Loading one file by handle: `spklef`, `spkuef`, `cklpf`, `ckupf`, `pcklof`, `pckuof` | [`furnsh`](@ref) and [`unload`](@ref) |
| Settings of the geometry finder: `gfsstp`, `gfstep`, `gfrepi`, `gfrepu`, `gfrepf`, `gfbail`, `gfinth`, `gfclrh`, `gfseth` | The step argument of each search (`udstep` in [`gfevnt`](@ref), [`gffove`](@ref) and [`gfocce`](@ref)) and the keywords `report` and `bail`; see [`gfdist`](@ref) |
| The current file of DAF and DAS: `dafcs`, `dafgh`, `daswbr` | The handle object passed to each call |

The routines whose names begin with `zz` are private to CSPICE and are left out as well.

The toolkit programs (`brief`, `chronos`, `mkspk`, `spkmerge` and the others) are not part of the package. SPICEApplications.jl runs them from Julia through the NAIF executables. [Coverage and validation](@ref) lists the areas of the library that exist.
