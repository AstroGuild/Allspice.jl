# Types of the dynamic frames (class 5), the switch frames (class 6) and the kernel-pool
# watchers. They are in a file of their own, included before `types.jl`, because
# `KernelSet` and `SpiceContext` have fields of these types. Every type is immutable
# after construction; the watcher state is the one exception and is guarded by the
# lock of its context.

"""
    DynVector

One of the two defining vectors of a two-vector dynamic frame, parsed from the
`FRAME_<id>_PRI_*` or `FRAME_<id>_SEC_*` variables. Plain bits.

- `kind`: 1 `OBSERVER_TARGET_POSITION`, 2 `OBSERVER_TARGET_VELOCITY`,
  3 `TARGET_NEAR_POINT`, 4 `CONSTANT`; 0 for a vector that is not used.
- `axis`: the frame axis (1, 2, 3) that the vector defines; `negate` is true for the
  negative axis.
- `corr`: flags of the aberration correction (the field of `Aberration`).
- `target`, `observer`: body IDs (for a constant vector only the observer, and only
  with an aberration correction).
- `frame`: frame of the velocity (kind 2), body-fixed frame of the target (kind 3) or
  frame of the constant vector (kind 4); `frclass` and `frcenter` are its class and
  its center.
- `vec`: the constant vector (kind 4) or the radii of the target (kind 3).

Contract: `kind in 0:4`; if `kind != 0` then `axis in 1:3`.
"""
struct DynVector
    kind::Int8
    axis::Int8
    negate::Bool
    corr::UInt8
    target::Int32
    observer::Int32
    frame::Int32
    frclass::Int32
    frcenter::Int32
    vec::SVector{3,Float64}
end

"""
    DynFrame

Definition of one parameterized dynamic frame (frame class 5), parsed once when the
snapshot is built. Declared `mutable` with `const` fields only so that a lookup gives a
reference and copies nothing.

- `id`, `name`, `center`: the frame; `base` is the frame of `FRAME_<id>_RELATIVE`.
- `family`: 1 mean equator and equinox of date, 2 true equator and equinox of date,
  3 mean ecliptic and equinox of date, 4 two-vector, 5 Euler, 6 product.
- `frozen`, `t0`: `FREEZE_EPOCH`; `inertial`: `ROTATION_STATE = 'INERTIAL'`;
  `base_inertial`: the base is one of the built-in inertial frames.
- `v1`, `v2`, `minsep`: two-vector family (`ANGLE_SEP_TOL` in radians).
- `axes`, `epoch`, `scale`, `coef1`, `coef2`, `coef3`: Euler family; `scale` is the
  number of degrees in one unit of the angles.
- `from`, `to`: product family, the frame IDs of `FROM_FRAMES` and `TO_FRAMES`.
- `errkind`, `errmsg`: `:none`, or the `SpiceError` that an evaluation of the frame
  throws because the definition is not valid (CSPICE reports it at that time).

Contract: if `errkind === :none` then `family in 1:6` and the fields of the family are
complete.
"""
mutable struct DynFrame
    const id::Int32
    const name::String
    const center::Int32
    const base::Int32
    const family::Int8
    const frozen::Bool
    const t0::Float64
    const inertial::Bool
    const base_inertial::Bool
    const v1::DynVector
    const v2::DynVector
    const minsep::Float64
    const axes::NTuple{3,Int32}
    const epoch::Float64
    const scale::Float64
    const coef1::Vector{Float64}
    const coef2::Vector{Float64}
    const coef3::Vector{Float64}
    const from::Vector{Int32}
    const to::Vector{Int32}
    const errkind::Symbol
    const errmsg::String
end

"""
    DynLevel{K}

Reference to the snapshot (`K` is `KernelSet`) with which the definitions of dynamic
frames are evaluated, one recursion level below the snapshot that holds it. Mutable, so
that it is passed as one pointer through the typed function pointer of `_dyn_call`.

Contract: `ks` has the data of the snapshot that holds this reference.
"""
mutable struct DynLevel{K}
    const ks::K
end

"""
    SwitchFrame

Definition of one switch frame (frame class 6), parsed once when the snapshot is
built: the base frames in the order of `FRAME_<id>_ALIGNED_WITH` (lowest priority
first), the class of each, and the optional time intervals `starts[i] <= et <=
stops[i]` (ET seconds). `errkind`, `errmsg` are as in `DynFrame`.

Contract: `length(classes) == length(bases)`; if `usetime` then `starts` and `stops`
have that length too and `starts[i] < stops[i]`.
"""
mutable struct SwitchFrame
    const id::Int32
    const bases::Vector{Int32}
    const classes::Vector{Int32}
    const usetime::Bool
    const starts::Vector{Float64}
    const stops::Vector{Float64}
    const errkind::Symbol
    const errmsg::String
end

"""
    PoolWatch

Watcher state of one `SpiceContext`, the per-context form of the agent tables of
`pool.c`. `agents` maps an agent name to the variable names that it watches, `notify`
holds the agents that have an update pending (what `cvpool` reports and clears), and
`version` counts the changes of the kernel pool of the context.

Contract: `agents` and `notify` are read and written only under the lock of the
context; `version` is atomic, never decreases, and can be read without the lock.
"""
mutable struct PoolWatch
    const agents::Dict{String,Vector{String}}
    const notify::Set{String}
    @atomic version::Int
    PoolWatch() = new(Dict{String,Vector{String}}(), Set{String}(), 0)
end
