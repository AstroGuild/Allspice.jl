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

"Copy of `s` that belongs to file number `file` of a `KernelSet`."
withfile(s::Segment, file::Integer) = Segment(s.target, s.center, s.frame, s.dtype, s.t0,
    s.t1, Int32(file), s.baddr, s.eaddr, s.init, s.intlen, s.rsize, s.nrec)

"""
    DAF

One opened binary DAF kernel. `data[a]` is the double at DAF address `a` (the whole
file viewed as doubles, memory-mapped when the byte order is native and converted
into an owned array when it is not). `kind` is `:SPK` or `:PCK`. `segments` are in
file order, so a later segment has priority over an earlier one.
"""
struct DAF
    path::String
    kind::Symbol
    data::Vector{Float64}
    segments::Vector{Segment}
end

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
    content::Union{DAF,TextKernel}
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
end

"""
    TimeDefaults

Immutable snapshot of the three defaults of `timdef`: `calendar` is `"GREGORIAN"`,
`"JULIAN"` or `"MIXED"`; `system` is `"UTC"`, `"TDB"`, `"TDT"`, `"TT"` or `""`; `zone`
is `""` or an offset such as `"UTC-7"`. At most one of `system` and `zone` is not
empty. Declared `mutable` with `const` fields only so that it lives behind one
pointer, as `KernelSet` does.
"""
mutable struct TimeDefaults
    const calendar::String
    const system::String
    const zone::String
end

"The defaults of CSPICE at start: Gregorian calendar, UTC, no time zone."
TimeDefaults() = TimeDefaults("GREGORIAN", "UTC", "")

"""
    SpiceContext()

Holder of the current `KernelSet`. Queries read the reference once and work on that
snapshot; `furnsh`, `unload` and `kclear` build a new snapshot under `lock` and swap
the reference atomically. Readers never take the lock. `timedefs` holds the defaults
of `timdef` in the same way: `timdef` replaces the snapshot under `lock`.
"""
mutable struct SpiceContext
    @atomic ks::KernelSet
    const lock::ReentrantLock
    @atomic timedefs::TimeDefaults
    SpiceContext(ks::KernelSet, lock::ReentrantLock) = new(ks, lock, TimeDefaults())
end

"Error type of the package. `kind` is a short symbol such as `:nocoverage`."
struct SpiceError <: Exception
    kind::Symbol
    msg::String
end
Base.showerror(io::IO, e::SpiceError) = print(io, "SpiceError(", e.kind, "): ", e.msg)
