# Differences from CSPICE

The package keeps the names, the argument order, the units and the numerical results of CSPICE N0067. It departs from CSPICE where the C interface carries hidden state, output arguments or arguments for C memory, and in the argument order of a few calls; this page lists each departure.

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
| Some arguments serve only the C memory: counts, sizes, `room`, workspaces | These arguments are absent: [`gdpool`](@ref) takes only the name, and [`gfdist`](@ref) has no `nintvls` |
| `spkobj_c`, `spkcov_c`, `pckfrm_c` and `pckcov_c` take the output cell last | [`spkobj!`](@ref), [`spkcov!`](@ref), [`pckfrm!`](@ref) and [`pckcov!`](@ref) take it first, as SPICE.jl does |
| `spksub_c` takes the handle of the new file last | [`spksub`](@ref) takes it first |
| A cell or a window has a size that the caller declares | A [`Window`](@ref) grows as needed |

The kind of an error is the short message of CSPICE in lower case where the same condition exists there, so `SPICE(INVALIDDEGREE)` becomes `:invaliddegree`.

SPICE.jl already returns its results, names a call that changes its argument with `!`, returns `nothing` for a failed lookup and counts from 1. For a program written with SPICE.jl, the rows that change are the static arrays, which cannot be changed in place where SPICE.jl returns a `Vector` or a `Matrix`; the error type, which is the [`SpiceError`](@ref) of this package with its field `kind` and not `SPICE.SpiceError`; the [`Window`](@ref) in the place of a cell; the context; and the arguments for C memory that SPICE.jl keeps, such as `start` and `room` of `gdpool`, `nintvls` of `gfdist` and `maxn` of `bodvrd`. A call that must return a value also throws where SPICE.jl returns `nothing`, as `gdpool` does for a variable that is not loaded.

## No state between calls

CSPICE keeps buffers and saved search results between calls, for speed. This package keeps no such state: a query is a function of the loaded kernels and of its arguments. Three visible consequences follow.

- **Threads.** Queries from many threads need no lock, and a kernel can be loaded while other threads query; see [Contexts and threads](@ref).
- **Settings.** The defaults of [`timdef`](@ref), [`tsetyr`](@ref) and [`tparch`](@ref), the tolerance of [`gfstol`](@ref) and the settings of [`dskstl`](@ref) belong to a context, not to the process. The tolerance of `gfstol` is the default of the keyword `tol` of each search, and a `tol` given to a search takes its place.
- **Repeatable answers.** Where the answer of CSPICE depends on an earlier call (for example `ckgp` directly after `ckgpav` on overlapping segments, the dates in a text kernel loaded after certain `str2et` calls, or the four padding bytes of each segment summary of a binary PCK, which CSPICE fills from the last SPK or CK summary it handled), this package returns the answer that CSPICE gives in a new process.

## Numerical results

Each routine follows the arithmetic of its CSPICE counterpart, operation by operation, including the order of the terms of a sum; fused multiply-add is not used. For the same kernels and the same `Float64` input, a result is equal to that of CSPICE bit for bit or differs by rounding; [Coverage and validation](@ref) says how this was checked and gives measured differences.

### Exact derivatives and the rates of CSPICE

Order 0 of [`spkder`](@ref) and of [`pxder`](@ref) is the result of `spkpos` and of `pxform`, and every higher order is the derivative of that result as a function of the epoch, with the velocity of the SPK data in the place of the derivative of its position. For an SPK type whose velocity is the derivative of its position, such as type 2 of the planetary ephemerides, the result is therefore the exact derivative. For a type that stores the velocity apart from the position, such as type 9, order 1 of a geometric position is the stored velocity, as in `spkezr`, and the partial of a dual-number epoch is the same velocity.

For a position, let `s` be -1 for reception and +1 for transmission. In the barycentric inertial frame, `spkder` differentiates

```
p(t) = r_targ(t + s ℓ(t)) - r_obs(t)
P(t) = R(t + s |pc(t)| / c) * A(p(t), s * v_obs(t))
```

where `ℓ` is the light time, `A` is the stellar aberration and `R` is the rotation from J2000 to the output frame. For the `LT` corrections `ℓ(t) = |r_targ(t) - r_obs(t)| / c`, the geometric light time, because the single iteration of CSPICE evaluates the target at that light time; for the `CN` corrections, which iterate the light-time equation until it converges, `ℓ` solves `ℓ(t) = |p(t)| / c`, and its derivatives satisfy that equation order by order. `A` is the identity without `+S`, and otherwise `A(p, v) = p + |p| sinφ v̂⊥ + |p| (cosφ - 1) p̂`, with `v⊥` the part of `v` perpendicular to `p` and `sinφ = |v⊥| / c`. `R` is constant for an inertial frame; for another frame, the epoch of `R` is corrected for the light time to the center of the frame, which is the body that the frame is attached to (the Earth for `IAU_EARTH`), and `pc` is `p` with that center in the place of the target.

The velocity of `spkezr` differs from order 1 for two reasons. CSPICE takes the rate of the light time from the converged light-time equation also for the `LT` corrections, which gives a relative difference of about 1e-8 for the Moon seen from the Earth (8.4e-9 measured with `de440s.bsp`). CSPICE also takes the acceleration of the observer for the `+S` corrections from a difference of velocities over ±1 s, which gives a relative difference of about 1e-13.

For a rotation, the table says what `pxder` differentiates and what `sxform` returns in its place.

| Frame | `pxder` differentiates | Derivative block of `sxform` |
|---|---|---|
| Inertial, PCK, fixed offset | the rotation | the same |
| CK | the interpolation of the CK segment, composed with the clock conversion of `sce2c`; the stored angular velocity is not used | built from the stored angular velocity; only segments that have one are used |
| Dynamic, mean equator or mean ecliptic and equinox of date | the precession and obliquity series | the same |
| Dynamic, true equator and equinox of date | the precession, nutation and obliquity series, with all terms | the nutation rates of CSPICE, which omit the drift of the coefficients |
| Dynamic, two-vector | each defining vector with its aberration correction, as in `spkder` | velocities of `spkezr`; a central difference for a velocity vector; zero for a constant vector |
| Dynamic, Euler | the polynomials of the angles | the same |
| Dynamic, frozen at the epoch that `FREEZE_EPOCH` gives in its definition | a constant: zero relative to the base frame | zero |
| Dynamic, treated as inertial by `ROTATION_STATE = 'INERTIAL'` in its definition | the rotating frame | zero |
| Switch, which selects one of several base frames by epoch | the base that `pxform` selects at the epoch | the block of that base |

For data that are polynomials in time (SPK types 2, 3, 8, 9, 12, 13, 14, 18, 20 and the binary PCK), the derivative of any order is the exact derivative of the polynomial. For the other types, and for frames that are defined by formulas, the package evaluates the algorithm on truncated Taylor series in time, which carry each derivative up to the requested order through every operation, and so obtains the derivatives of that algorithm without finite differences. In both cases a derivative has the smoothness of the data and no more: at the boundary between two records of an SPK, where a record is the stretch of time that one set of coefficients or states covers, the geometric position is continuous to rounding and the higher derivatives can jump. An aberration-corrected position also has small steps there, which [Known limits of accuracy](@ref) describes.

## What this package adds

| Addition | Calls |
|---|---|
| Derivatives of any order | [`spkder`](@ref), [`pxder`](@ref) |
| Dual-number epochs | the calls that [Dual-number epochs](@ref) lists |
| Tables for inner loops and GPUs | [`ephemeris_table`](@ref), [`rotation_table`](@ref) |
| Aberration correction as a value, parsed once | [`Aberration`](@ref): `ab = Aberration("LT+S")`, then `spkezr("MOON", et, "J2000", ab, "EARTH")` |
| Scratch files in memory | [`dasops`](@ref), [`ekops`](@ref) |

## What is absent

Of the 651 routines that the header of CSPICE N0067 declares, 162 have no counterpart of the same name. Each of them is a facility that Julia provides, a routine that NAIF has declared deprecated and replaced, or a part of the C interface that the call conventions above make unnecessary. The table gives the call to use in place of each group.

| CSPICE routines | In this package |
|---|---|
| Error handling: `chkin`, `chkout`, `erract`, `errch`, `errdp`, `errint`, `errprt`, `failed`, `getmsg`, `qcktrc`, `reset`, `return`, `setmsg`, `sigerr`, `trcdep`, `trcnam`, `trcoff` and others | A call throws [`SpiceError`](@ref); catch it with `try` |
| Cells and sets: `insrtd`, `removd`, `elemi`, `union`, `inter`, `diff`, `sdiff`, `scard`, `size`, `valid` and others | `Vector` and `Set` of Julia; windows are a [`Window`](@ref) |
| Arrays: `bsrchd`, `isrchi`, `lstled`, `orderd`, `reordd`, `shelld`, `maxd`, `sumad`, `filld`, `moved`, `mequ` and others | `searchsortedfirst`, `findfirst`, `sortperm`, `sort!`, `maximum`, `sum`, `fill!`, `copy` |
| Strings and parsing: `ucase`, `lcase`, `cmprss`, `lparse`, `nextwd`, `prsdp`, `prsint`, `repmc`, `matchw`, `lx4num` and others | `uppercase`, `lowercase`, `split`, `parse`, `replace`, `occursin`, `Regex` |
| Constants: `pi`, `dpmax`, `dpmin`, `intmax`, `intmin` | `π`, `floatmax`, `-floatmax`, `typemax(Int32)`, `typemin(Int32)` |
| Files and the command line: `exists`, `rdtext`, `prompt`, `getcml`, `putcml`, `ftncls` | `isfile`, `eachline`, `readline`, `ARGS`, `close` |
| Deprecated by NAIF: `bodvar`, `illum`, `spkapp`, `srfxpt`, `subpt`, `subsol`, `dafrda` | [`bodvrd`](@ref) or [`bodvcd`](@ref), [`ilumin`](@ref), [`spkaps`](@ref), [`sincpt`](@ref), [`subpnt`](@ref), [`subslr`](@ref), [`dafgda`](@ref) |
| Loading one file by handle: `spklef`, `spkuef`, `cklpf`, `ckupf`, `pcklof`, `pckuof` | [`furnsh`](@ref) and [`unload`](@ref); [`eklef`](@ref) and [`ekuef`](@ref) remain as `furnsh` and `unload` of a file that must be an EK |
| Settings of the geometry finder: `gfsstp`, `gfstep`, `gfrepi`, `gfrepu`, `gfrepf`, `gfbail`, `gfinth`, `gfclrh`, `gfseth` | The step argument of each search (`udstep` in [`gfevnt`](@ref), [`gffove`](@ref) and [`gfocce`](@ref)) and the keywords `report` and `bail`; see [`gfdist`](@ref) |
| The current file of DAF and DAS: `dafcs`, `dafgh`, `daswbr` | The handle object passed to each call |

The routines whose names begin with `zz` are private to CSPICE and are left out as well.

The toolkit programs (`brief`, `chronos`, `mkspk`, `spkmerge` and the others) are not part of the package. [SPICEApplications.jl](https://github.com/jdc-pub/SPICEApplications.jl) runs them from Julia through the NAIF executables. [Coverage and validation](@ref) lists the areas of the library that exist.
