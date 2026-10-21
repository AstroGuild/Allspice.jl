# Shared data types. Every type here is immutable: a `KernelSet` is a snapshot that
# queries read without a lock, and `furnsh`/`unload` replace it as a whole.

const Vec3{T} = SVector{3,T}
const Mat3{T} = SMatrix{3,3,T,9}

"""
    Segment

Plain-bits descriptor of one DAF segment (SPK or binary PCK).

For an SPK segment `target` moves relative to `center` and the data are in frame
`frame`. For a binary PCK segment `target` is the frame class ID of the body frame,
`frame` is the inertial reference frame and `center` is zero.

`baddr:eaddr` are 1-based DAF double-precision addresses into the data array of file
number `file` of the owning `KernelSet` (`file == 0` until the segment is placed in a
set). `init`, `intlen`, `rsize` and `nrec` are the four directory values that end a
segment of the fixed-interval Chebyshev types (SPK 2, 3, 20 and PCK 2), and are zero
for the other types.
"""
struct Segment
    target::Int32
    center::Int32
    frame::Int32
    dtype::Int32
    t0::Float64
    t1::Float64
    file::Int32
    baddr::Int
    eaddr::Int
    init::Float64
    intlen::Float64
    rsize::Int32
    nrec::Int32
end

withfile(s::Segment, file::Integer) = Segment(s.target, s.center, s.frame, s.dtype, s.t0,
    s.t1, Int32(file), s.baddr, s.eaddr, s.init, s.intlen, s.rsize, s.nrec)

"""
    CKSegment

Plain-bits descriptor of one segment of a CK (orientation) kernel: the orientation of
the structure `inst` relative to the frame `frame`, as data type `dtype`, for the
encoded spacecraft clock times `t0 <= sclk <= t1` (ticks). `avflag` tells whether the
segment holds angular velocity. `baddr:eaddr` and `file` have the meaning that they
have in `Segment`.
"""
struct CKSegment
    inst::Int32
    frame::Int32
    dtype::Int32
    avflag::Bool
    t0::Float64
    t1::Float64
    file::Int32
    baddr::Int
    eaddr::Int
end

withfile(s::CKSegment, file::Integer) = CKSegment(s.inst, s.frame, s.dtype, s.avflag,
    s.t0, s.t1, Int32(file), s.baddr, s.eaddr)

"""
    DAF

One opened binary DAF kernel. `data[a]` is the double at DAF address `a` (the whole
file viewed as doubles, memory-mapped when the byte order is native and converted
into an owned array when it is not). `kind` is `:SPK`, `:PCK` or `:CK`. `segments`
(SPK and PCK) and `cksegments` (CK) are in file order, so a later segment has priority
over an earlier one; the vector that does not fit `kind` is empty.
"""
struct DAF
    path::String
    kind::Symbol
    data::Vector{Float64}
    segments::Vector{Segment}
    cksegments::Vector{CKSegment}
end

DAF(path, kind, data, segments::Vector{Segment}) =
    DAF(path, kind, data, segments, CKSegment[])

"Value of a kernel-pool variable: numeric or character, always a vector."
const PoolValue = Union{Vector{Float64},Vector{String}}

"One assignment of a text kernel: `name = value` or, if `append`, `name += value`."
struct Assignment
    name::String
    append::Bool
    value::PoolValue
end

"One parsed text kernel: its assignments in file order."
struct TextKernel
    path::String
    assignments::Vector{Assignment}
end

"""
    Loaded

One entry of the load list of a `KernelSet`, in load order. `source` is the path of
the meta-kernel that caused the load, or `""` for a kernel loaded directly.
"""
struct Loaded
    path::String
    source::String
    content::Union{DAF,TextKernel,DSKFile,EKFile}
end

"""
    FrameDef

Definition of a reference frame. `class` is the SPICE frame class: 1 inertial,
2 PCK body-fixed, 3 CK, 4 fixed offset (TK), 5 dynamic. `classid` is the ID inside
the class (for class 2 the body or binary PCK frame class ID; for class 4 the ID used
in the `TKFRAME_<classid>_*` variables). `center` is the body at the frame origin.
"""
struct FrameDef
    id::Int32
    name::String
    class::Int32
    classid::Int32
    center::Int32
end

"""
    TKFrame

Fixed-offset frame: `rot` rotates vectors from this frame to frame `relative`
(the meaning of `TKFRAME_*_MATRIX`).
"""
struct TKFrame
    relative::Int32
    rot::Mat3{Float64}
end

"""
    BodyRot

Orientation constants of one body from a text PCK (`BODY<id>_POLE_RA` and so on), as
`tisbod` of CSPICE uses them. Angles are in degrees as stored. `pole_ra`, `pole_dec`
are polynomials in Julian centuries, `pm` is a polynomial in days, the `nut_prec_*`
vectors are the trigonometric amplitudes and `nut_prec_angles` holds the (offset,
rate per century) pairs of the body system, flattened as stored. `ref` is the
inertial reference frame (`BODY<id>_CONSTS_REF_FRAME`, default 1 = J2000) and
`jed_epoch` the reference epoch (`BODY<id>_CONSTS_JED_EPOCH`, default 2451545.0).
`phase_degree` is the polynomial degree of the nutation-precession angles of the body
system (`BODY<sys>_MAX_PHASE_DEGREE`, default 1), stored as given: the evaluation
rejects a value outside 1:3.
"""
struct BodyRot
    pole_ra::Vector{Float64}
    pole_dec::Vector{Float64}
    pm::Vector{Float64}
    nut_prec_ra::Vector{Float64}
    nut_prec_dec::Vector{Float64}
    nut_prec_pm::Vector{Float64}
    nut_prec_angles::Vector{Float64}
    ref::Int32
    jed_epoch::Float64
    phase_degree::Int32
end

"""
    LeapSeconds

Contents of a leap-second kernel: `DELTET/DELTA_T_A`, `K`, `EB`, `M` (two values) and
the `DELTET/DELTA_AT` table as (TAI−UTC seconds, epoch in TDB-like seconds past
J2000 of the UTC day start as the kernel `@date` gives it) pairs, flattened as stored.
"""
struct LeapSeconds
    delta_t_a::Float64
    k::Float64
    eb::Float64
    m::NTuple{2,Float64}
    delta_at::Vector{Float64}
end

"""
    SclkType1

Data of one type 1 spacecraft clock (`SCLK01_*` and `SCLK_PARTITION_*` variables of
clock `id`), as `sc01` of CSPICE uses them; built once by `build_sclk` and never
modified.

- `coef`: the coefficient records, flattened: encoded SCLK (ticks), parallel time
  (seconds past J2000 in the time system `timsys`), rate (seconds per most
  significant count), for each record.
- `pstart`, `pstop`: partition start and stop times, in ticks, as stored.
- `moduli`, `offsets`: modulus and offset of each of the `nfield` clock fields.
- `delim`: code of the output delimiter, 1 to 5 for `.`, `:`, `-`, `,` and blank.
- `timsys`: parallel time system, 1 for TDB and 2 for TDT.
- `tikmsc`: ticks in one count of the most significant field.
- `mxtick`: largest encoded SCLK, the rounded sum of the partition lengths.
"""
struct SclkType1
    id::Int32
    coef::Vector{Float64}
    pstart::Vector{Float64}
    pstop::Vector{Float64}
    moduli::Vector{Float64}
    offsets::Vector{Float64}
    nfield::Int32
    delim::Int32
    timsys::Int32
    tikmsc::Float64
    mxtick::Float64
end

"""
    KernelSet

Immutable snapshot of everything loaded. Built by `build_kernelset` from the load
list; never modified afterwards. Declared `mutable` with `const` fields only so that
it lives behind one pointer: the atomic reference in `SpiceContext` is then one word,
and readers on many threads do not contend.

- `loaded`: the load list, in load order.
- `daf`: data arrays of the DAF files; `Segment.file` indexes this vector.
- `spk`: SPK segments by target body, highest priority first (later file first, and
  inside a file later segment first), which is the search order of CSPICE.
- `pck`: binary PCK segments by frame class ID, in the same priority order.
- `pool`: kernel-pool variables after all assignments, in load order.
- `name2code`, `code2name`: body names defined in the pool (`NAIF_BODY_NAME` and
  `NAIF_BODY_CODE`); keys of `name2code` are normalized with `normalize_name`. The
  built-in table is consulted after these.
- `frames`, `framename2id`: frames defined in the pool; the built-in table is
  consulted before these, as in CSPICE. Keys of `framename2id` are normalized.
- `tk`: fixed-offset frames by frame ID.
- `bodyrot`: text-PCK orientation constants by body ID.
- `lsk`: leap-second data, or `nothing` if no leap-second kernel is loaded.
- `ck`: CK segments by instrument ID, in the priority order of `spk`.
- `sclk`: type 1 spacecraft clocks by clock ID (the spacecraft ID, a negative number).
- `ckmeta`: the pair (`CK_<id>_SCLK`, `CK_<id>_SPK`) for each CK ID that has one of
  these variables in the pool; a missing one of the two holds the default of `ckmeta`.
- `dyn`, `switch`: definitions of the dynamic frames (class 5) and of the switch frames
  (class 6) by frame ID, parsed from the pool (`dynframes.jl`, `switchframes.jl`).
- `inner`: the snapshot with the same data that the evaluation of a dynamic frame of
  this snapshot uses, or `nothing` at the recursion limit of CSPICE (`dynframes.jl`).
- `pooledited`: true if `pool` holds direct edits (`pdpool`, `lmpool`, ...), so that it
  is not the result of the text kernels of `loaded` alone (`pool_write.jl`).
- `bodydefs`: the name and code pairs of `boddef`, oldest first.
- `dskfiles`: the loaded DSK files, in load order; `DSKSegment.file` indexes this vector.
- `dsk`: DSK segments by center body ID, in the search order of CSPICE (later file
  first, and inside a file later segment first).
- `ekfiles`: the loaded EK files, in load order (a file that is loaded twice is here once).
- `ektables`: the tables of the loaded EK files, in load order; `EKTable.segs` indexes
  `ekfiles` (`ek/file.jl`).
"""
mutable struct KernelSet
    const loaded::Vector{Loaded}
    const daf::Vector{Vector{Float64}}
    const spk::Dict{Int32,Vector{Segment}}
    const pck::Dict{Int32,Vector{Segment}}
    const pool::Dict{String,PoolValue}
    const name2code::Dict{String,Int32}
    const code2name::Dict{Int32,String}
    const frames::Dict{Int32,FrameDef}
    const framename2id::Dict{String,Int32}
    const tk::Dict{Int32,TKFrame}
    const bodyrot::Dict{Int32,BodyRot}
    const lsk::Union{Nothing,LeapSeconds}
    const ck::Dict{Int32,Vector{CKSegment}}
    const sclk::Dict{Int32,SclkType1}
    const ckmeta::Dict{Int32,NTuple{2,Int32}}
    const dyn::Dict{Int32,DynFrame}
    const switch::Dict{Int32,SwitchFrame}
    const inner::Union{Nothing,DynLevel{KernelSet}}
    const pooledited::Bool
    const bodydefs::Vector{Tuple{String,Int32}}
    const dskfiles::Vector{DSKFile}
    const dsk::Dict{Int32,Vector{DSKSegment}}
    const ekfiles::Vector{EKFile}
    const ektables::Vector{EKTable}
end

"""
    TimeDefaults

Immutable snapshot of the three defaults of `timdef`: `calendar` is `"GREGORIAN"`,
`"JULIAN"` or `"MIXED"`; `system` is `"UTC"`, `"TDB"`, `"TDT"`, `"TT"` or `""`; `zone`
is `""` or an offset such as `"UTC-7"`. At most one of `system` and `zone` is not
empty. `yearbound` is the first year of the 100 years that a year of two digits
abbreviates (`tsetyr`; 1969 at start, so 69 is 1969 and 68 is 2068). `check` is the
range check of the components of a time string (`tparch`; off at start). Declared
`mutable` with `const` fields only so that it lives behind one pointer, as `KernelSet`
does.
"""
mutable struct TimeDefaults
    const calendar::String
    const system::String
    const zone::String
    const yearbound::Int
    const check::Bool
    function TimeDefaults(calendar, system, zone, yearbound::Integer=1969, check::Bool=false)
        @assert yearbound > 0 "the first year of the interval of abbreviated years is positive"
        return new(calendar, system, zone, yearbound, check)
    end
end

"The defaults of CSPICE at start: Gregorian calendar, UTC, no time zone."
TimeDefaults() = TimeDefaults("GREGORIAN", "UTC", "")

"""
    SpiceContext()

Create an empty kernel context: the loaded kernels and the settings that the queries read.

Every call that reads kernels takes a context as an optional first argument. Without
it the call uses [`default_context`](@ref), one context for the whole process.
Contexts are independent: [`furnsh`](@ref) on one context does not change another.

A context holds:

- the loaded kernels and the kernel pool ([`furnsh`](@ref), [`unload`](@ref),
  [`kclear`](@ref), [`pdpool`](@ref), [`boddef`](@ref));
- the time defaults of [`timdef`](@ref), [`tsetyr`](@ref) and [`tparch`](@ref);
- the pool watchers of [`swpool`](@ref);
- the DSK tolerances of [`dskstl`](@ref) and the search tolerance of [`gfstol`](@ref).

A context can be shared between threads. A query takes no lock. A call that changes
the context (a load, an unload, an edit of the pool or of a setting, a [`cvpool`](@ref)
that clears an update) is atomic: a query
that runs at the same time sees the context before the change or after it, never a part
of it. Two changes on two threads are applied one after the other.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> ktotal(ctx)
0

julia> pdpool(ctx, "MY_VALUES", [1.0, 2.0]);

julia> gdpool(ctx, "MY_VALUES")
2-element Vector{Float64}:
 1.0
 2.0

julia> expool("MY_VALUES")    # the default context is not changed
false
```

See also [`default_context`](@ref), [`furnsh`](@ref).
"""
mutable struct SpiceContext
    @atomic ks::KernelSet
    const lock::ReentrantLock
    @atomic timedefs::TimeDefaults
    const watch::PoolWatch
    @atomic dsktols::DSKTolerances
    @atomic gftol::Float64
    SpiceContext(ks::KernelSet, lock::ReentrantLock) =
        new(ks, lock, TimeDefaults(), PoolWatch(), DSKTolerances(), 1.0e-6)
end

