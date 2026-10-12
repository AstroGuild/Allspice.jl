# Internal contract for the prototype

Working note for the first build. Each source file has one owner; the signatures below
are fixed, and everything else inside a file is the owner's choice. `src/types.jl` and
`src/NeoSpice.jl` define the shared types and the include order.

Conventions that hold everywhere:

- Epochs are ET (TDB seconds past J2000). Every evaluator is generic over the epoch
  number type `T<:Real`, so `Float64`, a `ForwardDiff.Dual` and `Taylor{N,Float64}`
  all pass through. `primal(et)::Float64` gives the plain value used for all
  discrete decisions (segment search, record index, branch choice).
- A "jet of order N" is an `NTuple{N+1}` of ordinary time derivatives, orders 0 to N.
- Hot paths allocate nothing and return `StaticArrays` values.
- Errors are `throw(SpiceError(kind, msg))`.
- CSPICE N0067 source is at `/home/alex/mice/src/cspice` for reference; CSPICE is the
  oracle, so match its conventions (units km, km/s, radians in results).

## `taylor.jl`

- `struct Taylor{N,T} <: Real` with `c::NTuple{N+1,T}` holding Taylor coefficients
  (`c[k+1] = f⁽ᵏ⁾/k!`). Arithmetic `+ - *` between Taylor values and with plain reals,
  `/` by a plain real and by a Taylor value, `sin`, `cos`, `sincos`, `sqrt`, `zero`,
  `one`, `convert`, `promote_rule`.
- `taylor_var(t::T, ::Val{N}) -> Taylor{N,T}`: the variable `t` with unit derivative.
- `derivatives(x::Taylor{N,T}) -> NTuple{N+1,T}`: ordinary derivatives (`c[k+1]*k!`).
- `primal(x::Real)`: `Float64(x)` for plain reals, `primal(x.c[1])` for `Taylor`.

## `daf.jl`

- `open_daf(path::AbstractString) -> DAF`: validates the file record, handles both
  byte orders, walks the summary records, fills `Segment` descriptors (including the
  four directory values for SPK types 2, 3, 20 and PCK type 2).

## `spk_eval.jl`

- `find_segment(segs::Vector{Segment}, et::Float64) -> Int`: index of the first
  segment (the vector is in priority order) with `t0 <= et <= t1`, or 0.
- `spk_eval(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) -> NTuple{N+1,Vec3}`:
  position of `seg.target` relative to `seg.center` in frame `seg.frame` and its time
  derivatives to order N. Order 1 is the velocity that CSPICE returns. For types with
  an independent velocity series (3) order k >= 1 is derivative k-1 of the velocity
  series. Types 2 and 3 are required; other types throw `SpiceError(:unsupported, ...)`.
- `spk_eval_position(data, seg, et::T, ::Val{N})`: as above but always derivatives of
  the position series (differs from `spk_eval` only for type 3).
- `pck_eval(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) -> NTuple{N+1,Vec3}`:
  the stored 3-1-3 Euler angles (PHI, DELTA, W) in radians of a binary PCK type 2
  segment, as `pcke02` returns them, and their time derivatives. `PHI = RA + pi/2`,
  `DELTA = pi/2 - DEC`; the matrix from the segment's inertial frame to the body frame
  is `R3(W) * R1(DELTA) * R3(PHI)`.

## `textkernel.jl`

- `kernel_arch(path) -> Symbol`: `:SPK`, `:PCK` (binary), or `:TEXT`.
- `parse_text_kernel(path) -> TextKernel`: full text-kernel grammar (`\begindata`,
  `\begintext`, `=`, `+=`, numeric and string vectors, `D` exponents, `@` dates as TDB
  seconds past J2000, quoted strings with doubled quotes, `//` string continuation).
- `metakernel_files(tk::TextKernel) -> Vector{String}`: paths from `KERNELS_TO_LOAD`
  with `PATH_SYMBOLS`/`PATH_VALUES` substituted; empty if the kernel is not a
  meta-kernel.
- `apply_assignments!(pool::Dict{String,PoolValue}, tk::TextKernel)`.

## `names.jl`

- `normalize_name(s) -> String`: upper case, leading/trailing blanks removed, inner
  blank runs compressed to one.
- `build_names(pool) -> (name2code, code2name)` from `NAIF_BODY_NAME`/`NAIF_BODY_CODE`.
- `body_code(ks::KernelSet, name::AbstractString) -> Int32` (pool first, then the
  built-in table of `zzidmap.c`, then a string that is an integer); throws
  `SpiceError(:unknownbody, ...)`. `body_code(ks, id::Integer) = Int32(id)`.
- `body_name(ks, code) -> Union{String,Nothing}`.

## `frames.jl`

- `build_frames(pool) -> (frames::Dict{Int32,FrameDef}, framename2id)` from the
  `FRAME_<name>`, `FRAME_<id>_NAME/_CLASS/_CLASS_ID/_CENTER` variables.
- `frame_def(ks, id::Integer) -> FrameDef` (pool first, then built-in inertial frames
  of `chgirf.c` and built-in non-inertial frames of `zzfdat.c`); throws
  `SpiceError(:unknownframe, ...)`.
- `frame_id(ks, name::AbstractString) -> Int32`; `frame_id(ks, id::Integer)`.
- `frame_name(ks, id) -> Union{String,Nothing}`.
- `inertial_rotation(id::Integer) -> Mat3{Float64}`: constant matrix that rotates
  vectors from J2000 (frame 1) to built-in inertial frame `id`, as `irfrot` gives it.

## `time.jl`

- `build_lsk(pool) -> Union{Nothing,LeapSeconds}`.
- `deltet(ks, epoch::Real, eptype::AbstractString) -> Float64` (ET − UTC, as `deltet_c`).
- `utc2et(ks, s::AbstractString) -> Float64` and `str2et` for ISO 8601 calendar
  strings (`YYYY-MM-DDThh:mm:ss.fff`, optional time part, optional `UTC`/`TDB` label).
- `et2utc(ks, et, format::AbstractString, prec::Integer) -> String` for `"ISOC"`, `"C"`
  and `"J"`.

## `rotations.jl`

- `build_tk(pool, frames) -> Dict{Int32,TKFrame}` (`TKFRAME_<id>_SPEC` of `MATRIX`,
  `ANGLES`, `QUATERNION`; `_RELATIVE` resolved to a frame ID).
- `build_bodyrot(pool) -> Dict{Int32,BodyRot}`.
- `frame_rotation(ks::KernelSet, from::Int32, to::Int32, et::T) -> Mat3{T}`: matrix
  that rotates vectors from frame `from` to frame `to` at `et`, generic over `T`. With
  `et::Taylor` the result carries all derivative orders. Class 1, class 2 (binary PCK
  segment first if one covers the epoch, text-PCK constants otherwise, as `tisbod`
  does) and class 4 frames; classes 3 and 5 throw `SpiceError(:unsupported, ...)`.
- `rotation_jet(ks, from, to, et::Real, ::Val{N}) -> NTuple{N+1,Mat3}`: the same as
  ordinary derivatives.

## `context.jl`

- `build_kernelset(loaded::Vector{Loaded}) -> KernelSet`; `empty_kernelset()`.
- `SpiceContext()`, `default_context()`, `kernelset(ctx) -> KernelSet` (one atomic read).
- `furnsh([ctx], path)`, `unload([ctx], path)`, `kclear([ctx])`, `ktotal([ctx])`.

## `spk_api.jl`, `frame_api.jl`, `pool_api.jl`

Every public call has the form `f(ctx::SpiceContext, args...)`, plus `f(args...)` on
the default context, plus an inner form on a `KernelSet`.

- `spkgeo(targ::Integer, et, ref, obs::Integer) -> (SVector{6}, lt)`
- `spkezr(targ, et, ref, abcorr, obs) -> (SVector{6}, lt)`; `spkpos -> (SVector{3}, lt)`.
  Only `abcorr = "NONE"` in the prototype; `lt` is the one-way light time `|r|/c`.
- `spkder(targ, et, ref, obs, ::Val{N}) -> NTuple{N+1,SVector{3}}`.
- `pxform(from, to, et) -> SMatrix{3,3}`; `sxform(from, to, et) -> SMatrix{6,6}`;
  `pxder(from, to, et, ::Val{N}) -> NTuple{N+1,SMatrix{3,3}}`.
- `gdpool(name)`, `gcpool(name)`, `gipool(name)`, `bodvrd(body, item)`, `bodn2c`,
  `bodc2n`, `bods2c`, `namfrm`, `frmnam`.

## Test kernels on this machine

- `/home/alex/data/ephemerides/de440s.bsp`
- `/home/alex/.tudat/resource/spice_kernels/`: `naif0012.tls`, `pck00011.tpc`,
  `pck00010.tpc`, `moon_pa_de440_200625.bpc`, `moon_de440_200625.tf`, `moon_080317.tf`
- `/home/alex/malto/data/de405.bsp`
