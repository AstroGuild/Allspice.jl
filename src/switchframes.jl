# Switch frames (frame class 6): a frame that is aligned with the first of a list of
# base frames that has data at the request time, as `zzswffet_` and `zzswfxfm_` of
# CSPICE. The definition is parsed once when the snapshot is built (`_sw_parse`); the
# selection reads the parsed `SwitchFrame` only and does not allocate.

"""
    _sw_times(ks, key) -> Vector{Float64}

Contract: the values of the pool variable `key` as ET seconds: numbers as they are,
strings through `str2et` (which needs a leap-second kernel), as `zzswffet_`.
"""
function _sw_times(ks::KernelSet, key::AbstractString)
    v = ks.pool[key]
    v isa Vector{Float64} && return copy(v)
    return Float64[str2et(ks, s) for s in v]
end

"""
    _sw_parse(ks, fd) -> SwitchFrame

Contract: the definition of the class 6 frame `fd` from the pool of `ks`, with the
checks of `zzswffet_`: the base frame list `FRAME_<id>_ALIGNED_WITH` (names or IDs,
every base with a frame specification) and the optional `FRAME_<id>_START` and
`FRAME_<id>_STOP` (both or none, one value for each base, start before stop). The
function does not throw: a definition that is not valid gives a `SwitchFrame` with
`errkind !== :none`, and a use of the frame throws that error.
"""
function _sw_parse(ks::KernelSet, fd::FrameDef)
    try
        return _sw_parse_checked(ks, fd)
    catch err
        err isa SpiceError || rethrow()
        return SwitchFrame(fd.id, Int32[], Int32[], false, Float64[], Float64[], err.kind,
            err.msg)
    end
end

"Contract: as `_sw_parse`, but throws the `SpiceError` of the first check that fails."
function _sw_parse_checked(ks::KernelSet, fd::FrameDef)
    pool, id = ks.pool, fd.id
    kalign, kstart, kstop = "FRAME_$(id)_ALIGNED_WITH", "FRAME_$(id)_START", "FRAME_$(id)_STOP"
    list = get(pool, kalign, nothing)
    list === nothing && throw(SpiceError(:missingframevar, "kernel variable $kalign, " *
        "needed for the specification of the switch frame $id, is not in the kernel pool"))
    hasstart, hasstop = haskey(pool, kstart), haskey(pool, kstop)
    hasstart == hasstop || throw(SpiceError(:partialframespec, "only one of the kernel " *
        "variables $kstart and $kstop was provided for the switch frame $id; both are " *
        "needed if one is"))
    n = length(list)
    if hasstart && !(length(pool[kstart]) == length(pool[kstop]) == n)
        throw(SpiceError(:countmismatch, "kernel variables of the switch frame $id have " *
            "mismatched sizes: number of base frames = $n, number of start times = " *
            "$(length(pool[kstart])), number of stop times = $(length(pool[kstop]))"))
    end
    bases = Int32[]
    if list isa Vector{Float64}
        append!(bases, round.(Int32, list))
    else
        for name in list
            code = namfrm(ks, name)
            code == 0 && throw(SpiceError(:framenamenotfound, "base frame name $name of " *
                "the switch frame $id could not be translated to a frame ID code"))
            push!(bases, code)
        end
    end
    classes = Int32[]
    for b in bases
        r = frinfo(ks, b)
        r === nothing && throw(SpiceError(:frameinfonotfound,
            "no specification was found for the base frame $b of the switch frame $id"))
        push!(classes, r[2])
    end
    starts, stops = Float64[], Float64[]
    if hasstart
        starts, stops = _sw_times(ks, kstart), _sw_times(ks, kstop)
        k = findfirst(i -> starts[i] >= stops[i], 1:n)
        k === nothing || throw(SpiceError(:badtimebounds, "interval time bounds are not " *
            "strictly increasing at interval index $k for the switch frame $id: " *
            "$(starts[k]) TDB : $(stops[k]) TDB"))
    end
    return SwitchFrame(id, bases, classes, hasstart, starts, stops, :none, "")
end

@noinline _sw_throw_none(fd::FrameDef, et) = throw(SpiceError(:noframeconnect,
    "the switch frame $(fd.id) has no base frame with data at ET $(primal(et))"))
@noinline _sw_throw_missing(id) = throw(SpiceError(:missingframevar,
    "no specification of the switch frame $id is loaded"))

"""
    _sw_ck_step(mode, ks, fd, et) -> (found::Bool, ref::Int32, rot)

The link of the CK frame `fd` as `_ck_frame_step` gives it, but with `found == false`
(and an identity matrix) instead of an error if no CK segment has a pointing at `et`,
which is what lets a switch frame go on to its next base frame. In the state mode only
segments with angular velocity count, as in `ckfxfm_`.
"""
@inline function _sw_ck_step(::_CKExact, ks::KernelSet, fd::FrameDef, et::T) where {T<:Real}
    ticks = _sce2c(ks, _ck_frame_clock(ks, fd, et), et)
    found, cmat, _, _, ref = _ck_search(ks, fd.classid, ticks, 0.0, false)
    found || return false, _J2000_FRAME, one(Mat3{T})
    return true, ref, Mat3{T}(transpose(cmat))
end

@inline function _sw_ck_step(::_CKState, ks::KernelSet, fd::FrameDef, et::Float64)
    ticks = _sce2c(ks, _ck_frame_clock(ks, fd, et), et)
    found, cmat, av, _, ref = _ck_search(ks, fd.classid, ticks, 0.0, true)
    found || return false, _J2000_FRAME, one(Mat3{Taylor{1,Float64}})
    r = Mat3{Float64}(transpose(cmat))
    dr = _m_mxm(_ck_skew(av), r)
    return true, ref, map((a, b) -> Taylor((a, b)), r, dr)
end

"""
    _sw_frame_step(mode, ks, fd, et) -> (parent::Int32, rot)

One link of a frame chain for the class 6 frame `fd`, as `zzswfxfm_`: the base frames
are tried from the highest priority (the last of the list) down; a base whose time
interval does not hold `et` is skipped, and so is a CK base without a pointing at `et`.
For the selected base the result is the link of that base (classes 1 to 4), or the
base itself with an identity matrix (classes 5 and 6). Throws
`SpiceError(:noframeconnect)` if no base is selected, and the error of an invalid
definition. No allocation.
"""
function _sw_frame_step(mode, ks::KernelSet, fd::FrameDef, et::T) where {T<:Real}
    sw = get(ks.switch, fd.id, nothing)
    sw === nothing && _sw_throw_missing(fd.id)
    sw.errkind === :none || _dyn_throw_def(sw.errkind, sw.errmsg)
    S = _ck_numtype(mode, T)
    t = primal(et)
    for i in length(sw.bases):-1:1      # invariant: no base above i was selected
        if sw.usetime
            (@inbounds(sw.starts[i]) <= t <= @inbounds(sw.stops[i])) || continue
        end
        base = @inbounds sw.bases[i]
        class = @inbounds sw.classes[i]
        if class == 3
            found, ref, rot = _sw_ck_step(mode, ks, frame_def(ks, base), et)
            found && return ref, rot::Mat3{S}
        elseif class >= 5
            return base, one(Mat3{S})
        else
            parent, rot = _frame_step(ks, base, _ck_epoch(mode, et))
            return parent, rot::Mat3{S}
        end
    end
    _sw_throw_none(fd, et)
end
