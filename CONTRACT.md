# Internal contract

This file records, for the maintainers of Allspice.jl, the interfaces between its source
files: the signatures below are what other files rely on, and everything else inside a
file can change without notice. `src/types.jl` and `src/Allspice.jl` define the shared
types and the include order, and the sections below follow that order, so a file uses
only the files of the sections before it.

Conventions that hold everywhere:

- Epochs are ET (ephemeris time: TDB seconds past J2000). Every evaluator is generic
  over the epoch number type `T<:Real`, so `Float64`, a `ForwardDiff.Dual` (a dual
  number, which carries first derivatives) and `Taylor{N,Float64}` all pass through.
  `primal(et)::Float64` gives the plain value used for all discrete decisions (segment
  search, record index, branch choice).
- A "jet of order N" is an `NTuple{N+1}` of ordinary time derivatives, orders 0 to N.
- A "snapshot" is a `KernelSet`: an immutable value that holds everything loaded. A
  query reads the snapshot of its context once and takes no lock, and a change
  (`furnsh`, `unload`, `kclear`, a pool edit) builds a new snapshot under the lock of
  the context and swaps the reference atomically.
- Every public call that needs kernel data has three forms: `f(ks::KernelSet, args...)`,
  `f(ctx::SpiceContext, args...)`, and `f(args...)` on the default context.
- Hot paths allocate nothing and return `StaticArrays` values.
- Errors are `throw(SpiceError(kind, msg))`; `SpiceError` is abstract and this
  constructor returns a `SpiceKindError`. A failure that a caller handles in normal flow
  and that has data of its own gets a subtype with the fields `kind` and `msg` first
  (`SpkInsufficientDataError`, thrown by the helpers of `spk_api.jl`). The `kind` is the
  short message of CSPICE
  in lower case, for example `:invaliddegree` for SPICE(INVALIDDEGREE). A new kind goes
  into `ERROR_KINDS` (`errors.jl`): the constructor refuses a kind that is not there, and
  a test compares the registry with the kinds in the source.
- CSPICE keeps state in saved variables and in process-wide tables; this package keeps
  one, the table of open files of `daf_access.jl`. Other state that CSPICE hides is an
  object that the caller owns (`DAFHandle`, `DASHandle`, `EKResult`, `EKFastLoad`) or an
  atomic field of the `SpiceContext`.
- Where CSPICE computes with 32-bit integers, the port gives the same value: `_cint`
  (the C cast of a double) and `_wrap32` (reduction to the 32-bit range) of `time.jl`,
  with `_dnint` and `_nint`, the rounding functions of f2c (the Fortran to C translator
  that produced most of CSPICE).
- The C source of CSPICE N0067, which NAIF distributes with the toolkit, is the
  reference; CSPICE is the oracle, so match its conventions (units km, km/s, radians in
  results).

## `types.jl`, `dyntypes.jl`, `dsk_types.jl`, `ek/types.jl`

The shared data types. The three small files come before `types.jl` in the include
order because `Loaded`, `KernelSet` and `SpiceContext` have fields of their types.

| Type | File | Content |
|---|---|---|
| `Vec3{T}`, `Mat3{T}` | `types.jl` | `SVector{3,T}` and `SMatrix{3,3,T,9}` |
| `Segment`, `CKSegment` | `types.jl` | plain bits (a value without pointers): one SPK or binary PCK segment, one CK segment |
| `DAF` | `types.jl` | one opened binary kernel: `path`, `kind`, `data`, `segments`, `cksegments` |
| `PoolValue`, `Assignment`, `TextKernel` | `types.jl` | pool value (`Vector{Float64}` or `Vector{String}`) and a parsed text kernel |
| `Loaded` | `types.jl` | one load-list entry: `path`, `source`, `content` (`DAF`, `TextKernel`, `DSKFile` or `EKFile`) |
| `FrameDef`, `TKFrame`, `BodyRot`, `LeapSeconds`, `SclkType1` | `types.jl` | parsed pool data |
| `KernelSet`, `TimeDefaults`, `SpiceContext` | `types.jl` | the snapshot, the defaults of `timdef`, the context |
| `SpiceError`, `SpiceKindError`, `SpkInsufficientDataError`, `ERROR_KINDS` | `errors.jl` | the abstract error type, its two concrete types and the registry of the kinds |
| `DynVector`, `DynFrame`, `DynLevel{K}`, `SwitchFrame`, `PoolWatch` | `dyntypes.jl` | parsed dynamic and switch frames, the recursion link, the watcher state |
| `DLADescr`, `DSKDescr`, `DSK02Header`, `DSKSegment`, `DSKFile`, `DSKTolerances` | `dsk_types.jl` | shape-kernel descriptors, one loaded file, the six tolerances |
| `EKColumnAttr`, `EKSegmentSummary`, `_EKColumn`, `_EKSegment`, `EKFile`, `EKTable` | `ek/types.jl` | table-kernel segments, one loaded file, one catalog entry |

- Every type is immutable after construction, with one exception: the fields `agents`
  and `notify` of `PoolWatch`, which a caller reads and writes only under the lock of
  the context (`version` is atomic and never decreases).
- `KernelSet`, `TimeDefaults`, `DSKTolerances`, `DynFrame` and `SwitchFrame` are declared
  `mutable` with `const` fields only, so that each lives behind one pointer and an
  atomic reference to it is one word. `SpiceContext` has the atomic fields `ks`,
  `timedefs`, `dsktols` and `gftol`, the `lock`, and `watch::PoolWatch`.
- Priority order in `KernelSet.spk`, `pck`, `ck` and `dsk` is the search order of
  CSPICE: later file first, and inside one file later segment first. `withfile(s, file)`
  gives the copy of a `Segment` or `CKSegment` that belongs to file number `file`.

## `taylor.jl`

- `struct Taylor{N,T} <: Real`: the storage is the value `v::T` (coefficient 0) and the
  tuple `d::NTuple{N,T}` of the higher coefficients, because Julia cannot write the
  field type `NTuple{N+1,T}`. The property `x.c::NTuple{N+1,T}` joins them
  (`c[k+1] = f⁽ᵏ⁾/k!`), and `Taylor{N,T}(c)` splits a tuple of `N+1` coefficients.
  Arithmetic `+ - *` between Taylor values and with plain reals, `/` by a plain real
  and by a Taylor value, `inv`, `sin`, `cos`, `sincos`, `sqrt`, `zero`, `one`, `convert`,
  `promote_rule`. Comparisons are not defined: code that branches must use `primal`.
- `taylor_var(t::T, ::Val{N}) -> Taylor{N,T}`: the variable `t` with unit derivative.
- `derivatives(x::Taylor{N,T}) -> NTuple{N+1,T}`: ordinary derivatives (`c[k+1]*k!`).
- `primal(x::Real)`: `Float64(x)` for plain reals, `primal(x.v)` for `Taylor`.

## `daf.jl`

- `open_daf(path::AbstractString) -> DAF`: opens a DAF (double precision array file,
  the container of binary SPK, binary PCK and CK kernels). It validates the file
  record, handles both byte orders, walks the summary records, and fills the
  descriptors in file order: `Segment` for `:SPK` and `:PCK` (including the four
  directory values for SPK types 2, 3, 20 and PCK type 2), `CKSegment` for `:CK`.
  Postcondition: every segment has finite epochs and
  `1 <= baddr <= eaddr <= length(data)`. Throws `SpiceError(:filenotfound)` or `(:baddaf)`.
  `data` is the only memory map of the file that stays after the call.
- `_unmap!(a::Array)` removes the memory map of `a` at once (nothing for an array that
  is not mapped); `a` must not be read afterwards. Windows cannot delete or replace a
  file while a map of it exists, so a map that a call needs only for a while is removed
  with `_unmap!` and not left to the garbage collector: the byte map of `open_daf` and
  `_das_load`, `data` in `_ck_open(f, path, who)` (`ckobj`, `ckcov`) and in
  `kernel_arch`, and the address spaces of a DAS handle at `dascls`
  (`_das_release!(h)`). The maps of a `KernelSet` stay until the collector frees it,
  because a query can still hold the snapshot.
- `_daf_sizeok(size, psize, dsize, offset) -> (ok, n)` mirrors `zzsizeok_`, and
  `_daf_i32(raw, off, little)` reads one 32-bit integer; `daf_access.jl` uses both.

## `spk_eval.jl`

- `find_segment(segs::Vector{Segment}, et::Float64) -> Int`: index of the first
  segment (the vector is in priority order) with `t0 <= et <= t1`, or 0.
- `spk_eval(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) -> NTuple{N+1,Vec3}`:
  position of `seg.target` relative to `seg.center` in frame `seg.frame` and its time
  derivatives to order N. Order 1 is the velocity that CSPICE returns. For types with
  an independent velocity series (3) order k >= 1 is derivative k-1 of the velocity
  series. Types 2 and 3 are in this file; types 1, 5, 8, 9, 12 to 15 and 17 to 21 are in
  `spk_more.jl` and `spk_types/`. Type 10 and unknown types throw
  `SpiceError(:notsupported, ...)`: type 10 needs the leap seconds, which `spk_eval` does
  not have, so only a query through a `KernelSet` evaluates it (see `spk_api.jl`).
  Precondition, not checked: `seg.t0 <= primal(et) <= seg.t1`.
- `pck_eval(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) -> NTuple{N+1,Vec3}`:
  the stored 3-1-3 Euler angles (PHI, DELTA, W) in radians of a binary PCK type 2
  segment, as `pcke02` returns them, and their time derivatives. `PHI = RA + pi/2`,
  `DELTA = pi/2 - DEC`; the matrix from the segment's inertial frame to the body frame
  is `R3(W) * R1(DELTA) * R3(PHI)`.
- `_cheb3(data, a, ncoef, s::T, radius, ::Val{N}) -> NTuple{N+1,Vec3{T}}`: the sum of
  three Chebyshev series and its derivatives, also used by `table.jl`. For a plain
  epoch the value and the first derivative are equal to the last bit to those of
  `chbval_` and `chbint_`, which matters for an angle of many turns.

## `spk_more.jl`, `spk_types/chebyshev.jl`, `conic.jl`, `interp.jl`, `mda.jl`, `numeric.jl`

Evaluation of the SPK types other than 2 and 3. Each type has one function
`state(data, seg, et::T) -> (pos::Vec3{T}, vel::Vec3{T})` that follows the reader and
evaluator pair `spkrNN`/`spkeNN` of CSPICE, so a plain epoch gives the CSPICE state.

| Function | File | SPK types |
|---|---|---|
| `_spk_mda` | `mda.jl` | 1, 21 |
| `_spk5`, `_spk15`, `_spk17` | `conic.jl` | 5, 15, 17 |
| `_spk_equal`, `_spk_unequal`, `_spk18`, `_spk19` | `interp.jl` | 8 and 12, 9 and 13, 18, 19 |
| `_spk14`, `_spk20` | `chebyshev.jl` | 14, 20 |

- `_spk_more(data, seg, et::T, ::Val{N})` is the branch of `spk_eval` for these
  types: order k >= 2 comes from the same `state` function with a Taylor epoch. No
  allocation.
- A type other than 2, 3 and 20 has no directory values in its `Segment`: the
  evaluator reads them from the end of the segment at each call and checks them with
  `_seg_count`, which throws `SpiceError(:baddaf)`. `_lstltd(data, base, n, x)` and
  `_lstled(...)` mirror `lstltd` and `lstled` on a part of the data array.
- `_spk14_layout(data, seg) -> NamedTuple` reads and checks the generic-segment
  directory of a type 14 segment (bases, counts, packet size, Chebyshev degree);
  `_spk14` evaluates from it, and `spksub` reads the records with it.
- `_ad_order(::Type{T}) -> Int` is the number of derivative orders that `T` carries,
  and `_with_primal(x::T, p::Float64) -> T` replaces the plain value of `x` and keeps
  its derivative parts. An iterative solution runs on plain numbers as CSPICE does, and
  Newton steps then give the derivative parts.
- Interpolation kernels in `interp.jl`, shared by the SPK types 8, 9, 12, 13, 18 and
  19, the CK type 5 (`ck.jl`) and the public calls of `math/polynomials.jl`:
  `_neville!(w, n, xa, x, onequal)` (the recurrence of `lgrint`), `_neville_d!(w, dw, n,
  xa, x, onequal)` (of `lgrind`) and `_hermite_table!(w1, w2, n, xa, x, onequal)` (of
  `hrmint`). `xa(i)` is the i-th abscissa, the work arrays hold the values on entry and
  the result in their first place on return, and `onequal(i, j)` is called for two equal
  abscissas: it throws (`:dividebyzero` for the public calls, `:badsegment` for an SPK
  segment) or returns `nothing` for CK type 5, which CSPICE does not check. No
  allocation.
- `_chbigr(data, a, ncof, scale, mid, radius, x)` in `chebyshev.jl` is the integral of a
  Chebyshev expansion for SPK type 20 and for the public `chbigr`; `_eqn_state(et,
  epoch, a, h, k, l, p, q, dlpdt, mldt, nodedt, rapol, decpol) -> (pos, vel)` in
  `conic.jl` is the state of `eqncpv` for SPK type 17 and for the public `eqncpv`. It
  reduces the mean longitude with `_f2c_dmod`, the `d_mod` of f2c, and throws
  `SpiceError(:eccoutofbounds)` for an eccentricity of 0.9 or more, as `kepleq` does.

## `textkernel.jl`

- `kernel_arch(path) -> Symbol`: `:SPK`, `:PCK` or `:CK` for a binary DAF, `:DSK` or
  `:EK` for a DAS file with that ID word, and `:TEXT` for every other file. Throws
  `SpiceError(:filenotfound)` or, for another binary kernel type, `SpiceError(:notsupported)`.
- `parse_text_kernel(path, yearbound=1969, check=false) -> TextKernel`: full text-kernel grammar (`\begindata`,
  `\begintext`, `=`, `+=`, numeric and string vectors, `D` exponents, `@` dates as TDB
  seconds past J2000, quoted strings with doubled quotes). String elements are kept
  raw, so a `//` continuation marker stays in place until `pool_string_joined(values,
  nth, contin="//")` reads it. `parse_text_lines(path, lines, yearbound, check)` does the same on
  lines in memory; `yearbound` and `check` are those of the context for the `@` dates; both throw `SpiceError(:badkernel)` for a malformed assignment.
- `metakernel_files(tk::TextKernel) -> Vector{String}`: paths from `KERNELS_TO_LOAD`
  with `PATH_SYMBOLS`/`PATH_VALUES` substituted; empty if the kernel is not a
  meta-kernel.
- `apply_assignments!(pool::Dict{String,PoolValue}, tk::TextKernel) -> pool`: a stored
  vector is never mutated, so vectors shared with an older snapshot stay valid. Throws
  `SpiceError(:typemismatch)` when `+=` mixes numeric and string values.
- `_nparsd(s) -> Float64`: the double that `nparsd` of CSPICE computes, which is not
  always the correctly rounded one; throws `SpiceError(:badnumber)`.

## `names.jl`

- `normalize_name(s) -> String`: upper case, leading/trailing blanks removed, inner
  blank runs compressed to one.
- `build_names(pool) -> (name2code, code2name)` from `NAIF_BODY_NAME`/`NAIF_BODY_CODE`;
  throws `SpiceError(:badpool)` if the pair is not complete or not consistent.
- `body_code(ks::KernelSet, name::AbstractString) -> Int32` (pool first, then the
  built-in table of `zzidmap.c`, then a string that is an integer); throws
  `SpiceError(:idcodenotfound, ...)`. `body_code(ks, id::Integer) = Int32(id)`.
- `body_name(ks, code) -> Union{String,Nothing}`.

## `frames.jl`

- `build_frames(pool) -> (frames::Dict{Int32,FrameDef}, framename2id)` from the
  `FRAME_<name>`, `FRAME_<id>_NAME/_CLASS/_CLASS_ID/_CENTER` variables.
- `frame_def(ks, id::Integer) -> FrameDef` (built-in inertial frames of `chgirf.c` and
  built-in non-inertial frames of `zzfdat.c` first, then the pool, as `frinfo` does,
  then the rule of `frinfo` that an ID in 13000:13999 is a class 2 frame of the Earth
  with class ID `id - 10000`); throws `SpiceError(:unknownframe, ...)`.
- `frame_id(ks, name::AbstractString) -> Int32`; `frame_id(ks, id::Integer)`.
- `frame_name(ks, id) -> Union{String,Nothing}`.
- `inertial_rotation(id::Integer) -> Mat3{Float64}`: constant matrix that rotates
  vectors from J2000 (frame 1) to built-in inertial frame `id`, as `irfrot` gives it.

## `time.jl`, `timeparse.jl`, `timeout.jl`, `time_api.jl`

Every function with a `ks` argument has a form with a `SpiceContext` and a form on the
default context. The conversions need a leap-second kernel; the kind without one is that of CSPICE for
the call (`:noleapseconds`, `:missingtimeinfo` or `:kernelvarnotfound`).

- `build_lsk(pool) -> Union{Nothing,LeapSeconds}`; `_lsk(ks, kind=:noleapseconds) -> LeapSeconds`.
- `deltet(ks, epoch::Real, eptype::AbstractString) -> Float64` (ET − UTC, as `deltet_c`).
- `unitim(ks, epoch::T, insys, outsys) -> T` between `TAI`, `TDT`, `TT`, `GPS`, `TDB`,
  `ET`, `JDTDT`, `JDTDB`, `JED`; generic over `T`, no allocation.
- `ttrans(ks, from, to, tvec) -> NTuple{7,Float64}` (not exported): the 21
  representations of `ttrans` of SPICELIB. The inner form `_ttrans(l::LeapSeconds,
  pfrom, pto, tv)` takes the codes `_TT_*`, and `_tvec(x...)` builds the time vector.
- `str2et(ks, s, defs=TimeDefaults()) -> Float64`: every string form of `tpartv`.
- `utc2et(ks, s, defs=TimeDefaults()) -> Float64` and `tparse(s, defs) -> Float64`:
  the range checks of `tcheck` run only when `defs.check` is set, as CSPICE makes them
  only after `tparch("YES")`; `str2et` always makes them, as in CSPICE.
- `tpictr(sample) -> String`: the `timout` picture of an example string.
- `timout(ks, et, pictur, defs=TimeDefaults()) -> String`: every mark of `timout`.
- `et2utc(ks, et, format, prec) -> String` for `"C"`, `"D"`, `"J"`, `"ISOC"`, `"ISOD"`.
- `etcal(et) -> String`.
- `timdef([ctx], action, item, value="") -> String`: the defaults `CALENDAR`, `SYSTEM`
  and `ZONE` are an immutable `TimeDefaults` in an atomic field of the `SpiceContext`.
  A reader loads it once with `timedefaults(ctx)` and takes no lock; `SET` builds a new
  object under the lock of the context. `kclear` and `furnsh` do not change it.
- `tsetyr([ctx], year)` sets the field `yearbound` of the `TimeDefaults` of the context:
  the first year of the 100 years that a year of two digits abbreviates (1969 at
  start). `str2et`, `utc2et` and the context form of `tparse` use it; the forms on a
  `KernelSet` take the defaults as an argument, and `tparse(s)` uses 1969. The reader
  of text kernels converts an `@` date with `tparse` and the bound that the context
  has at the load (`furnsh`, `ldpool`, `lmpool`), as `zzrvar_` does, and with the range
  check of `tparch` that the context has at the load.
- `tparch([ctx], type)` sets the field `check` of the `TimeDefaults` of the context:
  `"YES"` (letter case and blanks ignored) turns the range check of `utc2et`, `tparse`
  and the `@` dates on, every other value turns it off, and `""` throws `:emptystring`.
  It is off at start, as in CSPICE; another context is not changed.
- `eklef([ctx], path) -> String` is `furnsh` for a file that is an EK
  (`:wrongarchitecture` otherwise) and returns the path, which `ekuef([ctx], path)`,
  that is `unload`, takes.
- `jul2gr(y, m, d)`, `gr2jul(y, m, d) -> (year, month, day, doy)`.
- Constants: `j1900 j1950 j2000 j2100 jyear tyear spd b1900 b1950 clight halfpi twopi
  rpd dpr`, each a function without arguments.
- Errors: `:noleapseconds`, `:missingtimeinfo`, `:badtime` (string or epoch), `:badtimetype` (system or
  representation), `:badepochtype` (`eptype` of `deltet`), `:badtimedefault` (`timdef`),
  `:badformat` (`et2utc` format, empty picture).

## `rotations.jl`

- `build_tk(pool, frames) -> Dict{Int32,TKFrame}` (`TKFRAME_<id>_SPEC` of `MATRIX`,
  `ANGLES`, `QUATERNION`; `_RELATIVE` resolved to a frame ID), keyed by frame ID, from
  `_tk_parse(pool, frames, classid) -> TKFrame`, which reads the variables as `tkfram_`
  does and throws its errors. A frame that does not parse is left out of the table;
  `_tk_missing(ks, fd)` parses it again at a query and throws the error.
- `build_bodyrot(pool) -> Dict{Int32,BodyRot}`.
- `frame_rotation(ks::KernelSet, from::Int32, to::Int32, et::T) -> Mat3{T}`: matrix
  that rotates vectors from frame `from` to frame `to` at `et`, generic over `T`. With
  `et::Taylor` the result carries all derivative orders. Class 1, class 2 (binary PCK
  segment first if one covers the epoch, text-PCK constants otherwise, as `tisbod`
  does), class 3 (`ck.jl`), class 4, class 5 (`dynframes.jl`) and class 6
  (`switchframes.jl`) frames; a higher class throws `SpiceError(:notsupported, ...)`.
  The two frame chains are followed to their first common frame, as `refchg_` does,
  for at most `MAX_FRAME_CHAIN` (10) links each (`SpiceError(:framechain)` beyond).
- `rotation_jet(ks, from, to, et::Real, ::Val{N}) -> NTuple{N+1,Mat3}`: the same as
  ordinary derivatives.
- `_frame_step(ks, id::Int32, et::T) -> (parent::Int32, rot::Mat3{T})` is one link of a
  chain (requires `id != 1`), and `_frame_parent(ks, id)` gives the parent without any
  evaluation. A class 4 frame without data throws the error of `_tk_missing`, and
  a class 3 frame without a CK segment at `et` throws `SpiceError(:noframeconnect)`.
- `_body_rotation(ks, classid::Int32, et::T) -> Mat3{T}`: J2000 to a class 2 frame, as
  `tisbod_`; throws `SpiceError(:framedatanotfound)`.

## `sclk.jl`

Spacecraft clocks (SCLK) of type 1, as `sc01.f` of CSPICE. An encoded SCLK value counts
"ticks", the units of the least significant clock field, from the start of the first
partition. The data of every clock are read from the pool once, when the snapshot is
built, so the conversions are pure functions of the snapshot.

- `build_sclk(pool) -> Dict{Int32,SclkType1}`: a clock with incomplete data is left
  out, so that a load never fails. `_sclk(ks, sc) -> SclkType1` then throws
  `SpiceError(:kernelvarnotfound)`, or `SpiceError(:notsupported)` for another clock type.
- `_sce2c(ks, clock, et::T) -> T` (`scec01`, continuous ticks) is piecewise linear and
  generic over the number type; it throws `SpiceError(:valueoutofrange)` outside the
  clock range. `_sce2t(ks, clock, et::Float64)` (`scet01`) rounds to whole ticks.
- The public calls `sct2e` and `sce2c` are generic too and do not allocate.

## `ck.jl`, `ck_api.jl`

A CK (orientation kernel) segment gives the C-matrix of an instrument: the matrix that
rotates vectors from the reference frame of the segment to the instrument frame, as a
function of ticks. The evaluators follow the readers `ckrNN` and the evaluators
`ckeNN` of CSPICE for types 1 to 6, but read the data array directly and keep no state.

- `_ck_search(ks, inst::Int32, sclk::T, tol::Float64, needav::Bool) -> (found, cmat,
  av, clkout, ref::Int32)`: the pointing from the first segment in priority order that
  gives one, as `ckbss_`/`cksns_` with `ckpfs_`. Requires `tol >= 0`; no allocation.
- Two evaluation modes, which `dynframes.jl` and `switchframes.jl` share. `_CKExact()`
  gives the rotation in the number type of the epoch (`ckfrot_`). `_CKState()` takes a
  `Float64` epoch and gives a `Mat3{Taylor{1,Float64}}` whose derivative part comes
  from the angular velocity of the segment (`ckfxfm_`), because the state
  transformation of CSPICE is not the time derivative of the rotation for a CK frame.
- `_ck_frame_step(mode, ks, fd::FrameDef, et) -> (ref::Int32, rot)` is one chain link of
  a class 3 frame, and `_ck_chain_to(mode, ks, id, stop, et) -> Mat3` follows a chain
  that may hold class 3 frames, for at most `MAX_FRAME_CHAIN` links, without allocation.
- `_ck_in_chain(ks, id) -> Bool` is true if the chain of `id` toward J2000 holds a frame
  of class 3, 5 or 6. For such a chain `frame_jet(ks, from::Int32, to::Int32,
  et::Float64, Val(1)) -> (R, dR)` of `ck_api.jl`, the jet behind `sxform`, `spkezr` and
  `spkgeo`, is `_ck_state_rotation` (the two blocks of `sxform_c`, as `frmchg_`); for
  every other chain it is `rotation_jet`. `_exact_frame_jet` always gives the exact
  derivatives.
- Public: `ckgp`, `ckgpav` (`nothing` without a pointing; `SpiceError(:valueoutofrange)`
  for `tol < 0`), `ckcov`, `ckcov!` (a `Window`, filled with `wninsd!` interval by
  interval as `ckcov_`), `ckobj`, `ckmeta` (`SpiceError(:unknownckmeta)`), and `ckfrot`,
  `ckfxfm`: one CK link as `_ck_frame_step` evaluates it, `nothing` without CK data,
  without the SCLK kernel (`ckhave_`, `zzsclk_`) or without a segment.
- `ck_records.jl` (after `daf_access.jl`, since it takes a `DAFHandle`): `cknr02`,
  `ckgr02`, `cknr03`, `ckgr03` read the records of one segment through a handle and a
  packed descriptor; `:ckwrongdatatype`, `:cknonexistrec`.
- Deviation from CSPICE, type 6: when `ckr06_` moves to the next mini-segment it rounds
  the clock rate of that mini-segment to a whole number; the stored rate is used here.

## `dynframes.jl`, `switchframes.jl`

Dynamic frames (class 5, the parameterized families of `zzdynfrm.c` and `zzdynrot.c`)
and switch frames (class 6: a frame aligned with the first base frame of a list that
has data at the request time, as `zzswffet_` and `zzswfxfm_`). The definitions are
parsed once when the snapshot is built; the evaluation reads the parsed objects only,
is a pure function of the snapshot and does not allocate.

- `_dyn_levels(base::KernelSet, edited::Bool, bodydefs) -> KernelSet`: the last step of
  `build_kernelset`. It parses the class 5 and class 6 frames (`_dyn_parse`,
  `_sw_parse`) and, if a dynamic frame exists, links a chain of `_DYN_LEVELS` (2) inner
  snapshots with the same data through `inner`; no data are copied.
- Recursion: a dynamic frame of `ks` is evaluated with the states and the frames of
  `ks.inner`. The innermost snapshot has `inner === nothing`, and a dynamic frame that
  is needed there throws `SpiceError(:recursiontoodeep)`, the limit of `zzfrmgt1_`.
- `_dyn_frame_step(mode, ks, fd, et) -> (base::Int32, rot)` and `_sw_frame_step(mode,
  ks, fd, et) -> (parent::Int32, rot)` are the chain links, with the modes of `ck.jl`.
  `_dyn_chain_step` selects between them, and `_dyn_base(ks, id)` gives the base frame
  without any evaluation.
- The parsers do not throw: a definition that is not valid gives an object with
  `errkind !== :none`, and a use of the frame throws that error. A switch frame with no
  base selected throws `SpiceError(:noframeconnect)`.

## `pool_write.jl`

The write interface of the kernel pool (`pdpool`, `pcpool`, `pipool`, `dvpool`,
`lmpool`, `ldpool`, `clpool`), the read calls `expool`, `dtpool`, `gnpool`, `szpool`,
`stpool`, the check `badkpv` (returns `false` or throws the kind of CSPICE), body name definitions (`boddef`) and the watchers (`swpool`, `cvpool`,
`poolversion`). An edit never changes a snapshot: it goes through the context as
`furnsh` does.

- `context.jl` calls three functions under the lock of the context:
  `_pool_x_swap_loaded!(ctx, list, entries)` from `furnsh`,
  `_pool_x_swap_unloaded!(ctx, list)` from `unload`, and `_pool_x_swap_cleared!(ctx)`
  from `kclear`. Each swaps in the new snapshot and notifies the watchers.
- Rule for the load list, which is the rule of CSPICE: an edit stays through a later
  `furnsh` and through the `unload` of a binary kernel; the `unload` of a text kernel
  builds the pool again from the text kernels that stay, so every edit is lost;
  `kclear` removes every edit but keeps the pairs of `boddef`; `clpool` empties the pool
  and leaves the load list as it is.
- Deviation from CSPICE: CSPICE keeps one global agent table, and here each
  `SpiceContext` has its own `PoolWatch`, so `swpool` and `cvpool` are thread safe and
  two contexts do not see the watches of each other. `poolversion(ctx)` is an atomic
  counter of the pool changes that a caller reads without a lock.

## `context.jl`

- `build_kernelset(loaded::Vector{Loaded}, edited=nothing, bodydefs=[]) -> KernelSet`:
  reads no file, since it uses the objects of the load-list entries. `edited` is a pool
  with direct edits that takes the place of the pool of the text kernels.
  `empty_kernelset()`.
- `SpiceContext()`, `default_context()`, `kernelset(ctx) -> KernelSet` (one atomic read).
- `furnsh([ctx], path)`, `unload([ctx], path)`, `kclear([ctx])`, `ktotal([ctx])`. If a
  load fails, the context is unchanged. `furnsh` throws `SpiceError(:blankfilename)` for
  a blank path, `SpiceError(:nosuchfile)` for a missing file and
  `SpiceError(:recursiveloading)` for a meta-kernel that lists a meta-kernel.

## `spk_api.jl`, `frame_api.jl`, `pool_api.jl`

Every public call has the form `f(ctx::SpiceContext, args...)`, plus `f(args...)` on
the default context, plus an inner form on a `KernelSet`.

- `spkgeo(targ::Integer, et, ref, obs::Integer) -> (SVector{6}, lt)`;
  `spkgps(targ::Integer, et, ref, obs::Integer) -> (SVector{3}, lt)` is its position
  part, from the order-0 jet.
- `spkezr(targ, et, ref, abcorr, obs) -> (SVector{6}, lt)`; `spkpos -> (SVector{3}, lt)`.
  `abcorr` is a CSPICE correction string or an `Aberration`; `lt` is the one-way light
  time.
- `spkder(targ, et, ref, obs, ::Val{N}) -> NTuple{N+1,SVector{3}}` (geometric) and
  `spkder(targ, et, ref, abcorr, obs, ::Val{N})`.
- `pxform(from, to, et) -> SMatrix{3,3}`; `sxform(from, to, et) -> SMatrix{6,6}`;
  `pxder(from, to, et, ::Val{N}) -> NTuple{N+1,SMatrix{3,3}}`.
- `gdpool(name)`, `gcpool(name)`, `gipool(name)`, `bodvrd(body, item)`, `bodn2c`,
  `bodc2n`, `bodc2s`, `bods2c`, `namfrm`, `frmnam`. `bodn2c` and `bods2c` return `nothing` for an
  unknown name, and `namfrm` returns 0 for an unknown frame.
- Internal: `spk_jet(ks, targ::Int32, et::T, ref::Int32, obs::Int32, ::Val{N}) ->
  NTuple{N+1,Vec3}` is the geometric position and its derivatives; it does not allocate
  and throws `SpkInsufficientDataError` (`:spkinsuffdata` or `:noloadedfiles`). It is the place where a
  type 10 segment is evaluated, with `_t10_jet` and the leap seconds of `ks`.
  `frame_jet(ks, from, to, et, ::Val{N})` is the rotation jet for any epoch type.
- `check_epoch(et)` throws `SpiceError(:badarg)` for an epoch that is not finite;
  `NameOrID = Union{AbstractString,Integer}`; `CLIGHT` is the value of `clight_c`.

## `aberration.jl`

- `struct Aberration` (plain bits; properties `uselt`, `converged`, `stellar`,
  `transmit`, `geometric`; the first four are the keywords of `Aberration(; ...)`, and
  the field `flags` is private); `Aberration(s::AbstractString)` parses as `zzvalcor` and
  throws `SpiceError(:invalidoption, ...)`. `AbcorrArg` is the union of the two argument forms.
- `abcorr_state(ks, targ::Int32, et, ref::Int32, corr, obs::Int32) -> (SVector{6}, lt)`
  (`spkez`), `abcorr_position(...) -> (SVector{3}, lt)` (`spkezp`),
  `abcorr_jet(..., ::Val{N}) -> NTuple{N+1,Vec3}`.
- Public: `spkssb(targ, et, ref)`, `spkltc(targ, et, ref, abcorr, stobs)`,
  `spkaps(targ, et, ref, abcorr, stobs, accobs)`, `spkacs(targ, et, ref, abcorr, obs)`,
  `spkapo(targ, et, ref, sobs, abcorr)`, `ltime(etobs, obs, dir, targ)`,
  `stelab(pobj, vobs)`, `stlabx(pobj, vobs)`, with the argument order of CSPICE.
- `spkez` and `spkezp` (integer IDs) call `spkezr` and `spkpos`, so their results are
  equal bit for bit.

## `windows.jl`

`windows.jl` comes before `ck_api.jl` in the include order, because `ckcov!` takes a
`Window`.

- `struct Window <: AbstractVector{Float64}`: `data` holds the endpoints of disjoint
  closed intervals in increasing order, so `w[k]`, `length(w)` and iteration give what
  a SPICE double cell gives. `Window(endpoints)` copies and validates with `wnvald!`;
  `_window_of(d)` takes ownership of a valid endpoint vector without a check.
- Public: `card`, `wncard`, `wninsd!`, `wnvald!`, `wnunid`, `wnintd`, `wndifd`, `wncomd`,
  `wnelmd`, `wnincd`, `wnfetd`, `wnsumd`, `wnreld`. The set operations return a new
  window. `wnexpd!`, `wncond!`, `wnextd!`, `wnfild!` and `wnfltd!` change the window in
  place, as their CSPICE routines do, and return it.
- Errors: `:badendpoints` (left > right), `:unmatchendpts` (odd count), `:nointerval`,
  `:invalidendpntspec` (side of `wnextd!`), `:invalidoperation` (operator of `wnreld`).
- Deviation from CSPICE: a window has no size limit, so SPICE(WINDOWEXCESS) cannot
  occur. The coverage calls and every geometry finder search return a `Window`.

## `daf_access.jl`, `daf_write.jl`

Low-level read and write access to DAF files. CSPICE keeps a process-wide table of
open files and one hidden "current array" for each search; here the opened file is a
`DAFHandle` that the caller owns, and the search state lives inside it.

| Group | Functions |
|---|---|
| open and close | `dafopr(path)`, `dafopw(path)`, `dafonw(path, ftype, nd, ni, ifname, resv)`, `dafcls(h)` |
| search and read | `dafbfs`, `dafbbs`, `daffna -> Bool`, `daffpa -> Bool`, `dafgs`, `dafgn`, `dafgda`, `dafgsr` |
| summary | `dafhsf(h) -> (nd, ni)`, `dafps(dc, ic)`, `dafus(sum, nd, ni) -> (dc, ic)`, `dafws`, `dafrs`, `dafrn` |
| add an array | `dafbna(h, sum, name)`, `dafada(h, data)`, `dafena(h)` |
| comments, file record | `dafac`, `dafec -> Vector{String}`, `dafdc`, `dafrfr` |
| architecture | `getfat(path) -> (arch, kind)`, `idw2at(idword) -> (arch, kind)` |

- Thread safety: the only global state is the table of open files, `_OPEN_FILES` under
  `_OPEN_FILES_LOCK`, so two tasks can work on two handles at the same time. One handle
  is not safe for concurrent use, and a caller that shares a handle between tasks must
  hold a lock around every call.
- The table of open files holds, with weak keys, each `DAFHandle` and `DASHandle` that a
  public open call returned (`dafopr`, `dafopw`, `dafonw`, `dasopr`, `dasopw`, `dasonw`,
  `ekopr`, and the calls on top of them), with the device and inode of its file.
  `_open_register!(h, path, writable, who; check=writable)` throws
  `SpiceError(:fileopenconflict)` for a write open of a file that another open handle
  has open, as `zzddhopn` of CSPICE does; a read open is not checked (CSPICE gives
  RWCONFLICT for a file open for write). A close calls `_open_unregister!(h)`. The
  handles of `furnsh` and other internal readers are not in the table.
- A file never shrinks, and the bytes of an array are in place before its summary is,
  so a reader that opens the file between two calls sees a valid DAF. Rule for the
  caller: `furnsh` a file only after `dafcls`, `spkcls` or `ckcls` has returned, and
  unload a kernel before an append, because a loaded map does not see later segments.
- `_daf_add_array!(h, dc, ic, name, data)` adds one complete array; the segment writers
  call it after all input checks. `_daf_check_write(h, who)` throws `:dafnosuchhandle` or
  `:dafinvalidaccess`, and `_daf_check_current(h, who)` throws `:nocurrentarray`.

## `spk_write.jl`, `ck_write.jl`, `pck_write.jl`

The segment writers: `spkopn`, `spkopa`, `spkcls`, `spkpds`, `spkw01` to `spkw21`
(types 1, 2, 3, 5, 8, 9, 12, 13, 15, 17 to 21), `ckopn`, `ckcls`, `ckw01`, `ckw02`,
`ckw03`, `ckw05`, `ckw06`, and `pckopn`, `pckw02`, `pckcls` (binary PCK type 2, whose
records `_seg_cheb_records` lays out as for SPK types 2 and 3). Each writer makes the input checks of the SPICELIB routine in
the same order, so the data area is equal bit for bit to the one that CSPICE writes
from the same input.

- Argument lists follow CSPICE in its order, without the counts and flags that follow
  from the input (`avflag` of `ckw01` and `ckw03` is `avvs !== nothing`; a count of
  `avvs` other than the record count is `:badarg`). A table of states, quaternions or
  packets is a vector of vectors (the form of SPICE.jl, the Julia wrapper of CSPICE), a
  matrix with one row of the table in each column, or a flat vector;
  `_seg_flat(who, what, x, rowlen) -> (flat, n, rowlen)` converts it and throws
  `SpiceError(:badarg)` if the shape does not fit.
- `SegFrame` is a frame name or a frame ID, and `SegContext` (the optional first
  argument `ctx` of every writer, of `spkpds` and of `table_to_spk`; the method without
  it passes `default_context()`) is the `SpiceContext` or `KernelSet` that resolves a
  name; `_seg_refcode(who, ctx, frame)` throws `SpiceError(:invalidrefframe)`.
  `_seg_close(h, what)` closes a file as `spkcls` does: a file that is open for write
  must hold a segment (`:nosegmentsfound`).
- Deviation from SPICE.jl 0.2.5: the record count of `ckw02` and `ckw03` is always the
  length of the input (that version passes one record less to CSPICE), and each element
  of a vector of vectors given to `pckw02` is one record of `3 * (degree + 1)`
  coefficients (that version passes `length(cdata[1]) - 1` as the degree).
- The base frame of `pckw02` is a built-in inertial frame only, found by `irfnum_`
  rules (letter case and blanks do not count) or given by its ID 1 to 21. The padding
  half-word of a binary PCK summary (5 integers) is zero, as a new CSPICE process
  writes it; a CSPICE process that handled an SPK or CK summary before writes there
  the sixth integer of that summary.
- Type 14 has three calls, as in CSPICE: `spk14b([ctx,] handle, segid, body, center,
  frame, first, last, chbdeg) -> SPK14Writer`, `spk14a(writer, [ncsets,] coeffs,
  epochs)` and `spk14e(writer)`. CSPICE keeps the segment in work in a table of files;
  here the `SPK14Writer` holds it, and the file changes only when `spk14e` runs. The data area
  is the one of the generic segment writer `sgwes_` (constant, each record after its
  epoch, epochs, directory, 17 meta data items), equal word for word.

## `das.jl`, `dla.jl`

A DAS file (direct access, segregated) is the container under shape kernels and table
kernels: three address spaces, each numbered from 1, of characters, doubles and 32-bit
integers. The DLA layer (DAS linked array) is a doubly linked list of segments in a
DAS file.

- `DASHandle`: owned by the caller, with the three address spaces as the plain arrays
  `chars`, `dps` and `ints` (`ints[a]` is the integer at DAS address `a`). The arrays
  of a read handle are memory-mapped for a segregated file of native byte order and
  are never modified; a write handle works on owned copies in memory.
- `dascls` of a write handle writes the complete, segregated file and moves it into
  place in one step, so a reader that has the old file mapped keeps a consistent
  image. Windows does not replace a mapped file, so there the move fails
  (`:daswritefail`) while a read handle or a loaded kernel maps the file. `dascls`
  removes the maps of the handle that it closes (`_das_release!`). One handle is not
  safe for concurrent use.
- `dasops()` gives a scratch handle: a write handle with an empty path, which `dascls`
  discards without a file. `ekops()` is the scratch EK on it (CSPICE makes a temporary
  file for both).
- `_das_load(path, who, writable) -> DASHandle` reads a file of either byte order and
  throws `SpiceError(:baddasdirectory)` for a directory that is not consistent.
  `_das_check_open` throws `:dasnosuchhandle`, and `_das_check_write` throws that kind
  or `:dasinvalidaccess`. `_das_check_update(h, who, first, last, n)` throws
  `:invalidaddress` unless both addresses of `dasudc`, `dasudd` or `dasudi` are in
  1:n, also for an empty range, as CSPICE does.
- `_dla_segments(ints, path) -> Vector{DLADescr}` gives all descriptors in file order
  and throws `SpiceError(:baddladescriptor)` if the list has a cycle.

## `dsk.jl`, `dsk_api.jl`, `dsk_write.jl`

Shape kernels (DSK, digital shape kernel) of type 2: a surface of triangular plates
with a grid of volume cells ("voxels") that selects the plates near a ray. The
geometry functions are ports of the CSPICE routines that the intercept needs, in the
order of the floating-point operations of CSPICE; all are pure and none allocates.

- `open_dsk(path) -> DSKFile`: a DAS file with the ID word `DAS/DSK`, of either byte
  order, segregated or not. Throws the kinds of `dasopr`, `:invalidfiletype`,
  `:invalidformat` or `:dasnosuchaddress`. Postcondition: every segment lies inside
  the arrays, so only the pointers inside the data can be out of range later.

| Function | Mirrors | Result | Error kinds |
|---|---|---|---|
| `_dsk_x02(ints, dps, seg, tols, vertex, raydir)` | `dskx02` | `(plid, xpt, found)` | `:rayiszerovector` |
| `_dsk_bux(ks, tols, bodyid, srflst, et, fixfid, vertex, raydir)` | `zzdskbux`, `zzsbfxri` | `(xpt, segidx, plid, found)`; `segidx` indexes `ks.dsk[bodyid]` | `:nodsksegments`, `:dskdatanotfound` |
| `_dsk_bun(ks, tols, bodyid, srflst, et, fixfid, point)` | `zzdskbun`, `zzsbfnrm` | outward unit normal | `:nodsksegments`, `:pointnotinsegment`, `:pointoffsurface` |
| `_dsk_sph(ks, segs, bodyid, srflst)` | `zzdsksph` | `(minrad, maxrad)` | `:dskdatanotfound` |

- All matching segments have the same priority and the nearest intercept wins, so
  `dskxsi` and `dskxv` accept `pri == false` only (`SpiceError(:badpriorityspec)`), as
  CSPICE does. `dskxsi` returns a `DSKSegment` in place of the triple (handle, DLA
  descriptor, DSK descriptor) of CSPICE.
- The tolerances of `dskgtl`/`dskstl` are a `DSKTolerances` in an atomic field of the
  context: `dsktolerances(ctx)` is one atomic read, and a call on a `KernelSet` without
  a context uses `_DSK_DEFAULT_TOLS`.
- `dsk_write.jl` (`dskopn`, `dskmi2`, `dskrb2`, `dskw02`, `dskcls`) has no cross-file
  interface beyond the public API; the spatial index has the layout of CSPICE.

## `surfaces.jl`

Surface names and codes (`srfc2s`, `srfs2c`, `srfcss`, `srfscc`), the computation
method strings, and the surface point calls `latsrf` and `srfnrm`. The surface name
map is read from the snapshot at each query, so there is no hidden state.

- `_dsk_lex_method(method) -> Vector{String}` mirrors `zzlexmet` (`:syntaxerror`), and
  `_dsk_parse_method(ks, bodyid, method) -> (shape, subtyp, pri, srflst)` mirrors
  `zzprsmet`, with `shape` equal to `:ellipsoid` or `:dsk`. The second throws
  `:badmethodsyntax`, `:badpriorityspec`, `:idcodenotfound` or `:toomanysurfaces`.

## `spk_types/sgp4.jl`, `spk_types/type10.jl`, `tle.jl`

SGP4 and SDP4 propagation of two-line element sets, SPK type 10, and the public calls
`getelm`, `evsgp4` and `spkw10`. The propagator gives a state in the TEME frame (true
equator, mean equinox), which `type10.jl` turns to J2000 as `zzteme.f` does.

- CSPICE keeps the initialized propagator in saved variables; here there is none.
  `_sgp4_init(geophs::NTuple{8,Float64}, elems::NTuple{10,Float64}, epoch, opmode) ->
  NamedTuple` returns the constants of one element set as a plain-bits record, and
  `_sgp4_eval(r, t::T) -> (pos, vel)` propagates that record to `t` minutes from the
  element epoch, so a `Taylor` time gives the derivatives of the state.
- Error kinds, those of CSPICE: `_sgp4_init` throws `:suborbital` or `:invalidelements`,
  and `_sgp4_eval` throws `:badmeanmotion`, `:badmeccentricity`, `:badmsemimajor`,
  `:badpeccentricity`, `:badsemilatus` or `:orbitdecay`.
- `_t10_jet(l::LeapSeconds, data, seg, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}` is the
  jet of a type 10 segment in the convention of `spk_eval`; it does not allocate, and
  it throws `SpiceError(:baddaf)` for a segment that is not valid. This type is the
  only one that needs a leap-second kernel, because the element epoch is turned to UTC.
- `_t10_closest(data, lay, x) -> Int` is the index of the reference epoch closest to
  `x`, a tie going to the later one (the index type 4 of the generic segments); the
  evaluator and `spksub` use it.
- `_t10_data(con, refs, packets) -> Vector{Float64}` is the whole data of a type 10
  segment (constants, records of a reference epoch and a 14-double packet, the epochs,
  their directory, the 17 meta items); `spkw10` and `spksub` write it with `dafada`.

## `spk_segments.jl`

`spksfs`, `spkpvn`, `spkuds` and `spksub`: SPK segments named by a file and a
descriptor. `spksub(handle, descr, ident, start, stop, newh)` takes the arguments of
`spksub_c` in its order.

- A loaded file has no handle, so `spksfs(ks, body, et) -> Union{Nothing,Tuple{String,
  Vector{Float64},String}}` returns the path of the file. The segment is
  `find_segment` on `ks.spk[body]`; the descriptor is built with `dafps` and the name is
  read from the name record of the mapped file, in either byte order. `spkpvn` takes
  either a `DAFHandle` or such a path with the kernel set that holds it (otherwise
  `:nosuchhandle`).
- `spkpvn` checks `first <= et <= last` for every type (`:timeoutofbounds`); CSPICE
  checks this for some types only, and the others extrapolate or read outside the
  segment.
- `_spksNN(r::_SegReader, baddr, eaddr, t0, t1) -> Vector{Float64}` port `spksNN` of
  CSPICE line by line, so that the data of the new segment are those of `spksub_c`,
  double for double. `_SegReader` holds the source segment and reads by DAF address,
  with `dafgda` outside the segment. Types 10 and 14 write through their own writers
  (`_spks10` with `_t10_data`, `_spks14` with `spk14b`).
- Defect of CSPICE: `spks19` reads `data[0]` and stops the process when the block of
  100 epochs in which the search for the stop epoch ends begins after that epoch;
  `_spks19` throws `:spkstructureerror` there.

## `kernel_info.jl`

Information calls: `spkobj`, `spkobj!`, `spkcov`, `spkcov!`, `ktotal(kind)`, `kdata`,
`kinfo`, `tkvrsn`.

- `_require_binary_kernel(who, path, kind)` checks with `getfat` and throws
  `:invalidformat` (a transfer file), `:invalidarchtype` or `:invalidfiletype`.
  `_each_summary(f, path)` calls `f(dc, ic, handle)` for every array of a DAF on a
  private read handle that it closes afterwards; `frame_info.jl` uses both helpers.
- Deviation from CSPICE: the handle that `kdata` and `kinfo` return is always 0,
  because a loaded kernel has no DAF handle in this package.

## `math/vectors.jl`, `rotations.jl`, `ellipsoid.jl`, `coordinates.jl`, `conics.jl`, `polynomials.jl`, `hex.jl`

The vector, matrix, rotation, ellipsoid, coordinate, conic, polynomial and
interpolation calls of CSPICE. Each
public function has the CSPICE name, returns its results instead of filling output
arguments, and follows the arithmetic of the CSPICE source, so a `Float64` input gives
the CSPICE result. The functions are generic over the number type, and every discrete
decision is made on `primal`.

- A call that has a "found" flag in CSPICE returns `nothing` when the flag is false. A
  matrix has the mathematical index order, and a Jacobian `J` has
  `J[i, j] = d out_i / d in_j`.
- `SpicePlane{T}` (unit `normal` and `constant >= 0` when made by `nvc2pl`, `nvp2pl` or
  `psv2pl`) and `SpiceEllipse{T}` (orthogonal semi-axes when made by `cgv2el`).
- Kernels without input conversion that the geometry files call: `_m_nearpt(positn, a,
  b, c) -> (npoint, alt)`, `_m_surfpt(positn, u, a, b, c) -> (found, point)`,
  `_m_surfnm`, `_m_cgv2el`, `_m_saelgv`, `_m_latrec`, `_m_reclat`, and vector helpers
  such as `_m_dot`, `_m_cross`, `_m_norm`, `_m_hat`.
- Each 3 by 3 product in the package goes through `_m_mxv`, `_m_mtxv`, `_m_mxm`,
  `_m_mtxm` or `_m_mxmt` and not through `*`: these add the three products of a
  component from left to right, as `mxv_c` and `mxm_c` do, while the `*` of StaticArrays
  lets the compiler fuse a multiplication with an addition, so its result changed with
  the Julia version (Julia 1.10 and 1.13 differed in the last bit).
- `hex.jl` holds `dp2hx` and `hx2dp`, the base 16 strings of SPICELIB; `hx2dp` throws
  `:badnumber` with the `errmsg` of `hx2dp_c` as the message.

## `frame_info.jl`, `lst.jl`

Frame information (`frinfo`, `cidfrm`, `cnmfrm`, `ccifrm`), fixed-offset frame data
(`tkfram`), body orientation (`tipbod`, `tisbod`), the contents of a binary PCK file
(`pckfrm`, `pckfrm!`, `pckcov`, `pckcov!`), the frame lists of a class (`bltfrm`,
`kplfrm`, sorted `Vector{Int}`, `:badframeclass`), body constants (`bodfnd`, `bodvcd`), local
solar time (`et2lst`) and the longitude of the Sun (`lspcn`). These two files have no
cross-file interface beyond the public API: later files call the forms on a
`KernelSet`, for example `frinfo(ks, frcode)`, which returns `(center, class, classid)`
or `nothing`.

## `geometry/`

The geometry calls of CSPICE N0067. `methods.jl` is the shared part, and the other
files hold one family each. The rule for the split of each call is the same: a public
function makes the checks of CSPICE in the order of CSPICE, and a private function
does the computation after the checks, which the geometry finder calls at many epochs.

| File | Public calls | Main functions that other files call |
|---|---|---|
| `methods.jl` | none | `_geo_method`, `_geo_abcorr`, `_geo_body`, `_geo_namfrm`, `_geo_frinfo`, `_geo_fixfid`, `_geo_radii`, `_geo_shape`, `_geo_rayx`, `_geo_spkezp` |
| `constobs.jl` | `spkcpo`, `spkcpt`, `spkcvo`, `spkcvt` | `_GeoFixed`, `_geo_cvo`, `_geo_cvt`, `_geo_pos` |
| `subpoints.jl` | `subpnt`, `subslr` | `_geo_subpnt`, `_geo_same_bodies` |
| `intercept.jl` | `sincpt`, `srfrec` | `_geo_sincpt`, `_geo_sfxcor`, `_geo_sub_in_method` |
| `illumination.jl` | `ilumin`, `illumg`, `illumf`, `phaseq`, `trgsep` | `_geo_illum_states` |
| `occultation.jl` | `occult` | `_GeoOccBody`, `_geo_occ_parse`, `_geo_occ_eval`, `_geo_occed` |
| `limbterm.jl` | `limbpt`, `termpt`, `tangpt` | none |
| `fov.jl` | `getfov`, `getfvn`, `fovray`, `fovtrg` | `_geo_getfov`, `_geo_fovaxi`, `_geo_elvupy`, `_geo_wind2d` |
| `azel.jl` | `azlrec`, `recazl`, `dazldr`, `drdazl`, `azlcpo` | none |

- `_geo_method(ks, bodyid, method) -> _GeoMethod` mirrors `zzprsmet`: the common
  spellings are compared as they are and do not allocate, and every other string goes
  through the lexer and the parser of `surfaces.jl`.
- Error kinds are those of the geometry calls of CSPICE. `_geo_abcorr` (`zzvalcor`)
  throws `:invalidoption` where `Aberration(s)` throws `:badarg`; `_geo_body` throws
  `:idcodenotfound`; `_geo_fixfid` throws `:noframe` (unknown frame) or `:invalidframe`
  (center is not the target); `_geo_radii` (`zzgftreb`) returns three positive values
  and throws `:kernelvarnotfound`, `:typemismatch`, `:invalidcount` or `:badaxislength`.
- Number types: `constobs.jl`, `subpoints.jl`, `intercept.jl` and `illumination.jl` are
  generic over the epoch type for an ellipsoid, while `occultation.jl`, `limbterm.jl`
  and `fov.jl` take a plain epoch because their decisions are discrete. The
  illumination angles need the arcsine, which `Taylor` does not have, so a dual number
  passes through and a `Taylor` epoch does not.

## `gf/`

The geometry finder: the search for the time windows on which a geometric condition
holds. `solver.jl` is the core and knows nothing about geometry; a search is described
by two functions of the epoch, the value of a quantity and the test "is the quantity
decreasing".

- `_GFRun{S,R,B,P}` holds the controls of one search (step, refinement function, `bail`
  and `report` callbacks, tolerance, `threads`), which CSPICE keeps in saved variables.
  `_gf_run(step, tol; refine=gfrefn, bail=nothing, report=nothing, threads=false)`
  builds it and throws `SpiceError(:invalidstep)` or `SpiceError(:invalidtolerance)`.
- `_gf_relx(value, decreasing, relate, refval, adjust, cnfine, run, prefix, suffix) ->
  Window` is the relation driver `zzgfrelx` for `"<"`, `"="`, `">"`, `"LOCMIN"`,
  `"LOCMAX"`, `"ABSMIN"`, `"ABSMAX"`; it throws `:notrecognized`, `:valueoutofrange` or
  `:invalidvalue`. `_gf_scan!(result, cond, w, run, cstep)` is the step-and-bracket
  search (`zzgfsolv`) on each interval of the window `w`; with `run.threads` each
  interval has its own task and the results are inserted in the order of `w`.
- A quantity is a small immutable subtype of `GFQuantity` that holds the snapshot and
  the checked inputs; its constructor is the initializer `zzgf**in` of CSPICE and makes
  the checks in that order. Three methods make it usable by the solver: `_gf_value(q,
  et::Float64) -> Float64`, `_gf_decreasing(q, et::Float64) -> Bool` and `_gf_label(q)`.
  A new quantity is a new type with these methods, and `solver.jl` does not change.

| Type | File | Quantity |
|---|---|---|
| `GFDistance`, `GFSeparation`, `GFRangeRate`, `GFPhase` | `quantities.jl` | distance, angular separation, range rate, phase angle |
| `GFCoordinate{V<:GFVector}` | `quantities.jl` | one coordinate of a position (`GFPositionVector`), of a sub-observer point (`GFSubObserverVector`) or of a surface intercept (`GFInterceptVector`, both in `surface.jl`) |
| `GFIllumination` | `illumination.jl` | phase, incidence or emission angle at a surface point |
| `GFFunction{F,D}` | `quantities.jl` | two functions of the caller (`gfuds`) |

- `_gf_search(q, relate, refval, adjust, cnfine, run) -> Window` is the call of
  `zzgfrelx` that `gfevnt` makes. The method for a `GFCoordinate` (`longitude.jl`)
  confines the search to the window on which the vector exists and solves longitude
  and right ascension with proxy functions, because their values jump at a branch cut.
- The test "is decreasing" is the formula of CSPICE and not the exact derivative of
  the value (range rate by a difference quotient, illumination by `zzilusta`, surface
  points by `zzgfssob` and `zzgfssin`), so that the brackets of the extrema are those of
  CSPICE. Occultation and field-of-view searches have no scalar quantity: each is one
  step-and-bracket search for a binary state.
- Deviations from CSPICE (`api.jl`): the result is a new `Window`, so there is no
  result cell, no workspace and no interval count; the step is an argument, and the
  tolerance, `bail` and `report` are keywords. `gfstol` sets the default tolerance of
  one `SpiceContext` (the atomic field `gftol`, which `kclear` does not change), and a
  search on a `KernelSet` uses `_GF_CNVTOL`, the default of CSPICE.
- Thread safety: a search reads the snapshot once, at its start, so many searches can
  run at the same time and a `furnsh` on another thread does not disturb them.

## `table.jl`, `ext/AllspiceAdaptExt.jl`, `ext/AllspiceForwardDiffExt.jl`

An ephemeris table holds the position of one body relative to another on a fixed
interval of time, as one flat array of Chebyshev coefficients on intervals of equal
length. A query needs no kernel set, so it suits the inner loop of an optimizer.

- `EphemerisTable{T<:AbstractFloat,A<:AbstractVector{T}}` with the fields `coef`, `t0`,
  `dt`, `nint`, `ncoef`, `target`, `observer`, `frame`, `maxerr`. Invariants, which the
  constructor checks (`SpiceError(:badarg)`): `dt > 0`, `nint >= 1`, `ncoef >= 1`,
  `length(coef) == 3 ncoef nint`. The table is an `isbits` value when its array is one.
- `ephemeris_table(targ, t0, t1, ref, abcorr, obs; tol, ncoef, maxintervals)` builds a
  table; `table_jet(tab, et, Val(N)) -> NTuple{N+1,SVector{3}}`, `table_position`,
  `table_state` and `table_span` read it without allocation and without kernel access.
  The range check throws `SpiceError(:spkinsuffdata)` and is a bounds check: under
  `@inbounds` it is removed, and an epoch outside the table is then extrapolated.
- `table_to_spk([ctx,] handle, tab, segid; first, last)` writes the table as one SPK
  segment of type 2 through `spkw02`, with the records of the table as its records.
- `RotationTable{T,A}` (fields `coef`, `t0`, `dt`, `nint`, `ncoef`, `from`, `to`,
  `maxerr`; `length(coef) == 9 ncoef nint`, the nine matrix elements column after
  column) is the counterpart for `pxform`: `rotation_table(from, to, t0, t1; tol,
  ncoef, maxintervals)` builds it, and `table_rotation(tab, et)` and
  `table_rotation_jet(tab, et, Val(N)) -> NTuple{N+1,SMatrix{3,3}}` read it with the
  contract of `table_jet`. `maxerr` is a Frobenius norm; the matrix is not made
  orthogonal again.
- `_tab_fit(f, start, dt, nint, ncoef, Val(C))` fits a function of `3C` components and
  `_tab_refine(...)` doubles the number of intervals; both table types use them.
- `AllspiceAdaptExt` defines `Adapt.adapt_structure` for both table types, which
  moves `coef` to another array type. `AllspiceForwardDiffExt` adds methods for a
  dual-number epoch to `primal`, `spk_jet`, `frame_rotation`, `frame_jet`, `abcorr_jet`,
  `abcorr_position`, `abcorr_state`, `_ad_order` and `_with_primal`: a query takes the
  jet at the plain epoch and assembles the dual result by the chain rule, so the
  partial of a result is the next order of `spkder` or `pxder`.

## `ek/`

The E-kernel (EK) subsystem: a kernel of tables that a query language reads. An EK
file is a DAS file with the ID word `DAS/EK` whose address spaces are divided into
pages; a segment holds rows of one table, and a row has one entry for each column.

| File | Content | Main functions that other files call |
|---|---|---|
| `file.jl` | pages, trees, segment and column descriptors | `open_ek`, `_ek_check`, `_ek_nseg`, `_ek_sdsc`, `_ek_find_column` |
| `columns.jl` | entries of the nine column classes | `_ek_nelt`, `_ek_elt_c`, `_ek_elt_d`, `_ek_elt_i`, `_ek_entry_c`, `_ek_entry_d`, `_ek_entry_i` |
| `query_parse.jl` | scan, parse, normalize, resolve names and times, check | `_ek_encode_query`, `_ek_psel` |
| `query_run.jl` | key selection, joins, union, ORDER BY | `EKResult`, `_ek_search`, `_ek_empty_result` |
| `api.jl` | public read and query calls | `_ek_check_read` |
| `write_tree.jl`, `write_columns.jl` | page allocator, trees, data pages, column entries, indexes | `_ekw_check`, `_ekw_pgan`, `_ekw_trin`, `_ekw_trdl`, the `_ekw_ac*`, `_ekw_ad*` and `_ekw_de*` routines |
| `write.jl` | public write calls, `EKFastLoad` | none |

- The read functions take a "store": any object with the fields `chars`, `dps` and
  `ints` (a `DASHandle` or an `EKFile`). They hold no state, so many tasks can read one
  store at the same time. `open_ek(path) -> EKFile` throws the kinds of `dasopr`,
  `:wrongarchitecture`, `:invalidformat`, `:eknosegments` or `:ekcolnummismatch`.
- `_ek_encode_query(tables, ks, td::TimeDefaults, query) -> (q, failed, errmsg)` is a
  pure function that returns the message of `ekfind_c` for a query that is not valid,
  and throws `SpiceError` only where CSPICE signals a SPICE error. Its docstring lists
  the known differences from CSPICE.
- `_ek_search(ks, q; compat=false) -> EKResult` mirrors `eksrch` and is a pure function
  of the snapshot and the query, so it runs on many threads. Three defects of CSPICE
  N0067 change results (`zzekweed`, `zzekillt`/`zzekille`, `zzekjsqz`; see the header of
  `query_run.jl`); `compat = true` reproduces them for the tests, and the public calls
  use `compat = false`.
- Deviation from CSPICE: `ekfind` returns an `EKResult` that the fetch calls (`ekgc`,
  `ekgd`, `ekgi`, `eknelt`) take as first argument. The object owns its row vectors and
  keeps the data that it points to, so it stays valid after an `unload`.
- The writers follow the C source statement by statement, because the words of a file
  after a sequence of calls must be the words that CSPICE writes. `EKFastLoad` holds
  the state of one fast load, so several loads can be open at once.

## `precompile.jl`

Precompile directives for the common call signatures. The file reads no kernel, since
the directives name argument types only, and it has no cross-file interface.

## Test layout

CSPICE, through SPICE.jl, is the oracle. It has one process-wide kernel state and is
not thread safe, so each group calls `oracle_load` once at its start and never calls
SPICE.jl from a spawned task.

| Item | File | Rule |
|---|---|---|
| Groups | `test/runtests.jl` | `GROUPS` lists 50 names; group `g` is the file `test/test_g.jl`. Run with `julia --startup-file=no --threads=8 --project=test test/runtests.jl [groups]`; an unknown name is an error. |
| Several processes | `test/parallel.jl` | `julia --startup-file=no --threads=8 test/parallel.jl [groups]`; `ALLSPICE_TEST_JOBS` processes (default a quarter of the CPU threads, 8 at most), balanced with the times of the last run in `test/.group_seconds`; exit status 1 if a process failed. |
| Exhaustive tier | `test/helpers.jl` | `ALLSPICE_TEST_FULL=1` sets `FULL`; `nsample(fast, full)` gives `fast` in the default run and `full` in the exhaustive run, and requires `0 <= fast <= full`. |
| Seed | `test/helpers.jl` | `ALLSPICE_TEST_SEED` sets `SEED` (default 20261011), the seed of every random sample. |
| Kernels | `test/kernels.jl` | `fetch_kernels` downloads the missing kernels into `test/data/`, with three attempts for each URL, and checks their SHA-256; `ALLSPICE_KERNEL_DIR`, `ALLSPICE_DE440S`, `ALLSPICE_DE405` override the paths and `ALLSPICE_DOWNLOAD_KERNELS=0` turns the download off. A group whose kernels are missing is skipped with a visible notice, or fails with `ALLSPICE_REQUIRE_KERNELS=1`. |
| Coverage | `test/coverage_summary.jl` | Run the tests with `--code-coverage=@src`, then `julia --startup-file=no test/coverage_summary.jl [--keep] [output file]`; it adds the `.cov` files of all processes, prints one line for each file of `src/` and `ext/` and the total, and deletes the `.cov` files unless `--keep` is given. |

- Time budget: the default run of all groups together has a budget of 120 s, and
  `runtests.jl` prints the time of each group, longest first, at the end of a run. A
  group keeps its large samples, long search spans and large files for the exhaustive
  tier.
- Each group must run alone, because `parallel.jl` gives the groups to the processes in
  any division. A group file therefore includes what it needs with a guard of the form
  `isdefined(@__MODULE__, :NAME) || include(...)`: the fixtures (`ckfixtures.jl`,
  `dskfixtures.jl`, `dynfixtures.jl`, `ekfixtures.jl`, `geomfixtures.jl`), the oracles
  (`gforacle.jl`, `mathoracle.jl`) and `apart.jl`, whose macro `@apart` makes the
  compiler handle each test set on its own.
- The groups `threads`, `gf_threads` and `ek_threads` skip their thread-safety tests
  with one thread and print a notice, so run the suite with `--threads=8`.

## Test kernels

The tests use the NAIF generic kernels `naif0012.tls`, `pck00011.tpc`, `pck00010.tpc`,
`moon_pa_de440_200625.bpc`, `moon_de440_200625.tf`, `moon_080317.tf`, `de440s.bsp` and
the optional `gm_de431.tpc` (the GM checks are skipped if it is absent), and
`de405.bsp` from JPL SSD, which has no binary-file-format string. `SOURCES` in
`test/kernels.jl` gives the URLs (tried in order, three attempts each) and the SHA-256 of
each.
