# Dynamic frames (frame class 5): the parameterized families of `zzdynfrm.c` and
# `zzdynrot.c`. The definitions are parsed once when the snapshot is built
# (`_dyn_parse`); the evaluation reads the parsed `DynFrame` only, is a pure function of
# the snapshot and does not allocate.
#
# Two modes, the ones of the CK frames (`ck.jl`):
#
# - `_CKExact()`: the rotation that `pxform` returns, as `zzdynrot_` computes it, in the
#   number type of the epoch. With a `Taylor` epoch the result holds the exact time
#   derivatives of that rotation.
# - `_CKState()`: the state transformation of `zzdynfrm_` at a `Float64` epoch, as a
#   matrix of `Taylor{1,Float64}` values whose derivative part is the derivative block of
#   `sxform`. It follows CSPICE where CSPICE does not differentiate the rotation: the
#   velocities of `spkezr`, the difference quotient of the unit velocity, the nutation
#   rates of `zzwahr_`, a constant vector with zero derivative, `ROTATION_STATE` and
#   `FREEZE_EPOCH`.
#
# Recursion: a dynamic frame of the snapshot `ks` is evaluated with the states and the
# frames of `ks.inner`, which holds the same data one level down. The innermost
# snapshot has `inner === nothing`, and a dynamic frame that is needed there throws
# `SpiceError(:recursiontoodeep)`. A definition can therefore refer to other dynamic
# frames, and those cannot, which is the limit of `zzfrmgt1_`.

"Family codes of `DynFrame.family`."
const _DYN_MEAN_EQ = Int8(1)
const _DYN_TRUE_EQ = Int8(2)
const _DYN_MEAN_ECL = Int8(3)
const _DYN_TWOVEC = Int8(4)
const _DYN_EULER = Int8(5)
const _DYN_PRODUCT = Int8(6)

"Number of snapshots below the one that the user holds: levels 1 and 2 of CSPICE."
const _DYN_LEVELS = 2

"A `DynVector` that is not used."
const _DYN_NOVECTOR = DynVector(0, 1, false, 0x00, 0, 0, 0, 0, 0, zero(SVector{3,Float64}))

# ---------------------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------------------

"""
    _dyn_lookup(pool, id, name, item, required) -> Union{Nothing,Tuple{String,PoolValue}}

Contract: name and value of the frame definition variable `FRAME_<id>_<item>`, or of
`FRAME_<name>_<item>` if the first is absent, with the rules of `zzdynvac_`: a name
longer than 32 characters is not searched for. Throws
`SpiceError(:varnametoolong)` if both names are too long and
`SpiceError(:kernelvarnotfound)` if `required` and no variable exists; returns `nothing`
if not `required` and no variable exists.
"""
function _dyn_lookup(pool::AbstractDict, id::Integer, name::AbstractString,
        item::AbstractString, required::Bool)
    k1 = "FRAME_$(id)_$(item)"
    k2 = "FRAME_$(name)_$(item)"
    if length(k1) <= _TK_MAXNAME
        v = get(pool, k1, nothing)
        v === nothing || return k1, v
    end
    if length(k2) > _TK_MAXNAME
        length(k1) > _TK_MAXNAME && throw(SpiceError(:varnametoolong,
            "kernel variables $k1 and $k2 are longer than $_TK_MAXNAME characters"))
    else
        v = get(pool, k2, nothing)
        v === nothing || return k2, v
    end
    required || return nothing
    throw(SpiceError(:kernelvarnotfound, "at least one of the kernel variables $k1 or " *
        "$k2 is needed to define the frame $name, and neither is in the kernel pool"))
end

@noinline _dyn_badtype(key, name, want) = throw(SpiceError(:badvariabletype,
    "the kernel variable $key of the definition of frame $name must be $want"))
@noinline _dyn_badsize(key, name, maxn, n) = throw(SpiceError(:badvariablesize,
    "the kernel variable $key of the definition of frame $name has $n values; at most " *
    "$maxn are allowed"))

"Contract: the character values (at most `maxn`) of a required definition item, as `zzdynvac_`."
function _dyn_chars(pool, id, name, item, maxn::Int)
    key, v = _dyn_lookup(pool, id, name, item, true)
    v isa Vector{String} || _dyn_badtype(key, name, "character")
    length(v) > maxn && _dyn_badsize(key, name, maxn, length(v))
    return v
end

"Contract: the single character value of a required definition item."
_dyn_char(pool, id, name, item) = _dyn_chars(pool, id, name, item, 1)[1]

"Contract: the single character value of an optional definition item, or `nothing`, as `zzdynoac_`."
function _dyn_optchar(pool, id, name, item)
    r = _dyn_lookup(pool, id, name, item, false)
    r === nothing && return nothing
    key, v = r
    v isa Vector{String} || _dyn_badtype(key, name, "character")
    length(v) > 1 && _dyn_badsize(key, name, 1, length(v))
    return v[1]
end

"Contract: the numeric values (at most `maxn`) of a required definition item, as `zzdynvad_`."
function _dyn_nums(pool, id, name, item, maxn::Int)
    key, v = _dyn_lookup(pool, id, name, item, true)
    v isa Vector{Float64} || _dyn_badtype(key, name, "numeric")
    length(v) > maxn && _dyn_badsize(key, name, maxn, length(v))
    return v
end

"Contract: the single numeric value of an optional definition item, or `nothing`, as `zzdynoad_`."
function _dyn_optnum(pool, id, name, item)
    r = _dyn_lookup(pool, id, name, item, false)
    r === nothing && return nothing
    key, v = r
    v isa Vector{Float64} || _dyn_badtype(key, name, "numeric")
    length(v) > 1 && _dyn_badsize(key, name, 1, length(v))
    return v[1]
end

"""
    _dyn_frameid(ks, id, name, item) -> Int32

Contract: frame ID of a required definition item that holds a frame name or a frame
ID, as `zzdynfid_`. Throws `SpiceError(:notranslation)` for a name that is not a frame
name and not an integer.
"""
function _dyn_frameid(ks::KernelSet, id, name, item)
    key, v = _dyn_lookup(ks.pool, id, name, item, true)
    length(v) > 1 && _dyn_badsize(key, name, 1, length(v))
    v isa Vector{Float64} && return round(Int32, v[1])
    code = namfrm(ks, v[1])
    code != 0 && return Int32(code)
    n = tryparse(Int32, strip(v[1]))
    n === nothing && throw(SpiceError(:notranslation, "the kernel variable $key of the " *
        "definition of frame $name has the value '$(v[1])', which is not a frame name"))
    return n
end

"""
    _dyn_bodyid(ks, id, name, item) -> Int32

Contract: body ID of a required definition item that holds a body name or a body ID,
as `zzdynbid_`. Throws `SpiceError(:notranslation)` for an unknown name.
"""
function _dyn_bodyid(ks::KernelSet, id, name, item)
    key, v = _dyn_lookup(ks.pool, id, name, item, true)
    length(v) > 1 && _dyn_badsize(key, name, 1, length(v))
    v isa Vector{Float64} && return round(Int32, v[1])
    code = bods2c(ks, v[1])
    code === nothing && throw(SpiceError(:notranslation, "the body name '$(v[1])' of the " *
        "kernel variable $key could not be translated to an ID code"))
    return Int32(code)
end

"`s` without blanks and in upper case, as `cmprss_` with zero blanks and `ucase_`."
_dyn_squeeze(s::AbstractString) = uppercase(replace(s, ' ' => ""))

"True if `a` and `b` are equal when blanks and letter case are ignored, as `eqstr_`."
_dyn_eqstr(a::AbstractString, b::AbstractString) = _dyn_squeeze(a) == _dyn_squeeze(b)

"""
    _dyn_degrees(units) -> Float64

Contract: the number of degrees in one angle unit `units`, the factor of `convrt_`.
An angle `x` is converted to radians as `x * f / (180 / pi)`, the two operations of
`convrt_`, so that the result has the bits of CSPICE. Throws `SpiceError(:badunits)`
for a name that `_angle_factor` does not accept.
"""
function _dyn_degrees(units::AbstractString)
    u = uppercase(strip(units))
    u == "RADIANS" && return 180.0 / pi
    u == "DEGREES" && return 1.0
    u == "ARCMINUTES" && return 0.016666666666666666
    u == "ARCSECONDS" && return 2.7777777777777778e-4
    u == "HOURANGLE" && return 15.0
    u == "MINUTEANGLE" && return 0.25
    u == "SECONDANGLE" && return 0.0041666666666666666
    return _angle_factor(units) * (180.0 / pi)
end

"`Aberration(s)`; a value that is not valid throws `SpiceError(:spkinvalidoption)`, as `zzvalcor_`."
function _dyn_parse_abcorr(s::AbstractString)
    try
        return Aberration(s)
    catch e
        e isa SpiceError || rethrow()
        throw(SpiceError(:spkinvalidoption, e.msg))
    end
end

"""
    _dyn_abcorr(s, constant::Bool) -> Aberration

Contract: the correction of a vector definition. For an observer-target vector the
nine values of `spkezr`. For a constant vector also `S` and `XS`, and a value with
both light time and stellar aberration throws `SpiceError(:invalidoption)`, as
`zzdynfrm_` does.
"""
function _dyn_abcorr(s::AbstractString, constant::Bool)
    constant || return _dyn_parse_abcorr(s)
    u = _dyn_squeeze(s)
    u == "S" && return _aberration(_AB_S)
    u == "XS" && return _aberration(_AB_S | _AB_X)
    c = _dyn_parse_abcorr(s)
    c.stellar && throw(SpiceError(:invalidoption, "aberration correction '$s' of a " *
        "constant vector: light time and stellar aberration are mutually exclusive"))
    return c
end

"Class and center of frame `id`; throws `SpiceError(:framedatanotfound)` if the frame is unknown."
function _dyn_frameinfo(ks::KernelSet, id::Int32, name)
    r = frinfo(ks, id)
    r === nothing && throw(SpiceError(:framedatanotfound, "the definition of frame " *
        "$name refers to the frame with ID $id, which has no frame specification"))
    return Int32(r[2]), Int32(r[1])
end

"""
    _dyn_parse_vector(ks, fd, prefix) -> DynVector

Contract: the vector definition `FRAME_<id>_<prefix>*` (`prefix` is `PRI_` or `SEC_`)
of the two-vector frame `fd`, with the checks of `zzdynfrm_` in its order. Throws the
`SpiceError` of the first check that fails.
"""
function _dyn_parse_vector(ks::KernelSet, fd::FrameDef, prefix::AbstractString)
    pool, id, name = ks.pool, fd.id, fd.name
    ax = _dyn_squeeze(_dyn_char(pool, id, name, prefix * "AXIS"))
    negate = startswith(ax, '-')
    ax = replace(ax, '-' => "", '+' => "")
    axis = ax == "X" ? 1 : ax == "Y" ? 2 : ax == "Z" ? 3 : 0
    axis == 0 && throw(SpiceError(:invalidaxis, "the definition of frame $name has the " *
        "axis '$ax'; the valid values are X, -X, Y, -Y, Z, -Z"))
    vdef = _dyn_squeeze(_dyn_char(pool, id, name, prefix * "VECTOR_DEF"))
    zero3 = zero(SVector{3,Float64})
    if vdef == "OBSERVER_TARGET_POSITION"
        targ = _dyn_bodyid(ks, id, name, prefix * "TARGET")
        obs = _dyn_bodyid(ks, id, name, prefix * "OBSERVER")
        corr = _dyn_abcorr(_dyn_char(pool, id, name, prefix * "ABCORR"), false)
        return DynVector(1, axis, negate, corr.flags, targ, obs, 1, 1, 0, zero3)
    elseif vdef == "OBSERVER_TARGET_VELOCITY"
        targ = _dyn_bodyid(ks, id, name, prefix * "TARGET")
        obs = _dyn_bodyid(ks, id, name, prefix * "OBSERVER")
        corr = _dyn_abcorr(_dyn_char(pool, id, name, prefix * "ABCORR"), false)
        frid = _dyn_frameid(ks, id, name, prefix * "FRAME")
        class, center = _dyn_frameinfo(ks, frid, name)
        return DynVector(2, axis, negate, corr.flags, targ, obs, frid, class, center, zero3)
    elseif vdef == "TARGET_NEAR_POINT"
        targ = _dyn_bodyid(ks, id, name, prefix * "TARGET")
        obs = _dyn_bodyid(ks, id, name, prefix * "OBSERVER")
        corr = _dyn_abcorr(_dyn_char(pool, id, name, prefix * "ABCORR"), false)
        fix = cidfrm(ks, targ)
        fix === nothing && throw(SpiceError(:framedatanotfound, "the definition of " *
            "frame $name needs the body-fixed frame of the target body $targ, and no " *
            "frame is associated with this body"))
        radii = bodvcd(ks, targ, "RADII")
        length(radii) == 3 || throw(SpiceError(:invalidcount,
            "only $(length(radii)) axes were found for ID $targ; three axes expected"))
        all(>(0.0), radii) || throw(SpiceError(:badaxislength,
            "an axis of body $targ is negative or zero: BODY$(targ)_RADII = $radii"))
        class, center = _dyn_frameinfo(ks, Int32(fix[1]), name)
        return DynVector(3, axis, negate, corr.flags, targ, obs, fix[1], class, center,
            SVector{3,Float64}(radii))
    elseif vdef == "CONSTANT"
        frid = _dyn_frameid(ks, id, name, prefix * "FRAME")
        cs = _dyn_optchar(pool, id, name, prefix * "ABCORR")
        corr = _dyn_abcorr(cs === nothing ? "NONE" : cs, true)
        class, center, obs = Int32(1), Int32(0), Int32(0)
        if !corr.geometric
            class, center = _dyn_frameinfo(ks, frid, name)
            if corr.stellar || class != 1
                obs = _dyn_bodyid(ks, id, name, prefix * "OBSERVER")
            end
        end
        spec = _dyn_squeeze(_dyn_char(pool, id, name, prefix * "SPEC"))
        if spec == "RECTANGULAR"
            v = _dyn_nums(pool, id, name, prefix * "VECTOR", 3)
            length(v) == 3 || _dyn_badsize("FRAME_$(id)_$(prefix)VECTOR", name, 3, length(v))
            dir = SVector{3,Float64}(v)
        elseif spec == "LATITUDINAL" || spec == "RA/DEC"
            f = _dyn_degrees(_dyn_char(pool, id, name, prefix * "UNITS"))
            lonitem, latitem = spec == "LATITUDINAL" ? ("LONGITUDE", "LATITUDE") : ("RA", "DEC")
            lon = _dyn_nums(pool, id, name, prefix * lonitem, 1)[1] * f / (180.0 / pi)
            lat = _dyn_nums(pool, id, name, prefix * latitem, 1)[1] * f / (180.0 / pi)
            dir = SVector{3,Float64}(cos(lon) * cos(lat), sin(lon) * cos(lat), sin(lat))
        else
            throw(SpiceError(:notsupported, "the definition of frame $name has the " *
                "constant vector specification '$spec', which is not supported"))
        end
        return DynVector(4, axis, negate, corr.flags, 0, obs, frid, class, center, dir)
    end
    throw(SpiceError(:notsupported, "the definition of the two-vector frame $name has " *
        "the vector definition '$vdef', which is not supported"))
end

"""
    _dyn_parse(ks, fd) -> DynFrame

Contract: the definition of the class 5 frame `fd` from the pool of `ks`, with the
checks of `zzdynfrm_` in its order. The function does not throw: a definition that is
not valid gives a `DynFrame` with `errkind !== :none`, and an evaluation of the frame
throws that error, as CSPICE reports it at evaluation time.
"""
function _dyn_parse(ks::KernelSet, fd::FrameDef)
    try
        return _dyn_parse_checked(ks, fd)
    catch err
        err isa SpiceError || rethrow()
        return DynFrame(fd.id, fd.name, fd.center, 0, 0, false, 0.0, false, false,
            _DYN_NOVECTOR, _DYN_NOVECTOR, 0.0, (Int32(0), Int32(0), Int32(0)), 0.0, 1.0,
            Float64[], Float64[], Float64[], Int32[], Int32[], err.kind, err.msg)
    end
end

"Contract: as `_dyn_parse`, but throws the `SpiceError` of the first check that fails."
function _dyn_parse_checked(ks::KernelSet, fd::FrameDef)
    pool, id, name = ks.pool, fd.id, fd.name
    base = _dyn_frameid(ks, id, name, "RELATIVE")
    style = _dyn_char(pool, id, name, "DEF_STYLE")
    _dyn_eqstr(style, "PARAMETERIZED") || throw(SpiceError(:notsupported,
        "dynamic frame style '$style' (in the definition of frame $name) is not supported"))
    family = _dyn_squeeze(_dyn_char(pool, id, name, "FAMILY"))
    fam = family == "MEAN_EQUATOR_AND_EQUINOX_OF_DATE" ? _DYN_MEAN_EQ :
          family == "TRUE_EQUATOR_AND_EQUINOX_OF_DATE" ? _DYN_TRUE_EQ :
          family == "MEAN_ECLIPTIC_AND_EQUINOX_OF_DATE" ? _DYN_MEAN_ECL :
          family == "TWO-VECTOR" ? _DYN_TWOVEC : family == "EULER" ? _DYN_EULER :
          family == "PRODUCT" ? _DYN_PRODUCT : Int8(0)
    ofdate = 1 <= fam <= 3
    freeze = _dyn_optnum(pool, id, name, "FREEZE_EPOCH")
    frozen = freeze !== nothing
    rotsta = _dyn_optchar(pool, id, name, "ROTATION_STATE")
    frozen && rotsta !== nothing && throw(SpiceError(:framedeferror, "the definition " *
        "of frame $name has both FREEZE_EPOCH and ROTATION_STATE; at most one is allowed"))
    ofdate && !frozen && rotsta === nothing && throw(SpiceError(:framedeferror,
        "the definition of frame $name of the family $family has neither FREEZE_EPOCH " *
        "nor ROTATION_STATE; exactly one is required"))
    inertial = false
    if rotsta !== nothing
        inertial = _dyn_eqstr(rotsta, "INERTIAL")
        inertial || _dyn_eqstr(rotsta, "ROTATING") || throw(SpiceError(:notsupported,
            "the definition of frame $name has ROTATION_STATE = '$rotsta'; the valid " *
            "values are ROTATING and INERTIAL"))
    end
    v1 = v2 = _DYN_NOVECTOR
    minsep = 0.0
    axes = (Int32(0), Int32(0), Int32(0))
    epoch, scale = 0.0, 1.0
    coef = (Float64[], Float64[], Float64[])
    from, to = Int32[], Int32[]
    if ofdate
        for (item, model, on) in (("PREC_MODEL", "EARTH_IAU_1976", true),
                ("NUT_MODEL", "EARTH_IAU_1980", fam == _DYN_TRUE_EQ),
                ("OBLIQ_MODEL", "EARTH_IAU_1980", fam == _DYN_MEAN_ECL))
            on || continue
            value = _dyn_char(pool, id, name, item)
            _dyn_eqstr(value, model) || throw(SpiceError(:notsupported, "the definition " *
                "of frame $name has $item = '$value', which is not recognized"))
            fd.center == 399 || throw(SpiceError(:invalidselection, "the definition of " *
                "frame $name has the center $(fd.center) and $item = '$model'; this " *
                "model is not applicable to that body"))
        end
    elseif fam == _DYN_TWOVEC
        v1 = _dyn_parse_vector(ks, fd, "PRI_")
        v2 = _dyn_parse_vector(ks, fd, "SEC_")
        sep = _dyn_optnum(pool, id, name, "ANGLE_SEP_TOL")
        minsep = sep === nothing ? 1e-3 : sep
        v1.axis == v2.axis && throw(SpiceError(:undefinedframe, "the definition of " *
            "frame $name gives the same axis to both vectors"))
    elseif fam == _DYN_EULER
        epoch = _dyn_nums(pool, id, name, "EPOCH", 1)[1]
        key, a = _dyn_lookup(pool, id, name, "AXES", true)
        a isa Vector{Float64} || _dyn_badtype(key, name, "integer")
        length(a) == 3 || _dyn_badsize(key, name, 3, length(a))
        axes = (round(Int32, a[1]), round(Int32, a[2]), round(Int32, a[3]))
        coef = ntuple(i -> copy(_dyn_nums(pool, id, name, "ANGLE_$(i)_COEFFS", 20)), 3)
        scale = _dyn_degrees(_dyn_char(pool, id, name, "UNITS"))
        all(x -> 1 <= x <= 3, axes) || throw(SpiceError(:badaxisnumbers,
            "the Euler axes of frame $name are $axes; each must be 1, 2 or 3"))
    elseif fam == _DYN_PRODUCT
        fnames = _dyn_chars(pool, id, name, "FROM_FRAMES", 10)
        tnames = _dyn_chars(pool, id, name, "TO_FRAMES", 10)
        length(fnames) == length(tnames) || throw(SpiceError(:badframecount,
            "the product frame $name has $(length(fnames)) \"from\" frames and " *
            "$(length(tnames)) \"to\" frames; these counts must match"))
        from = Int32[namfrm(ks, s) for s in fnames]
        to = Int32[namfrm(ks, s) for s in tnames]
        k = findfirst(==(0), vcat(from, to))
        k === nothing || throw(SpiceError(:unknownframe, "the product frame $name " *
            "refers to the frame '$(vcat(fnames, tnames)[k])', which is not known"))
    else
        throw(SpiceError(:notsupported, "dynamic frame family '$family' (in the " *
            "definition of frame $name) is not supported"))
    end
    return DynFrame(id, name, fd.center, base, fam, frozen, frozen ? freeze : 0.0, inertial,
        1 <= base <= N_INERTIAL_FRAMES, v1, v2, minsep, axes, epoch, scale, coef[1],
        coef[2], coef[3], from, to, :none, "")
end

"""
    KernelSet(loaded, daf, spk, pck, pool, name2code, code2name, frames, framename2id, tk,
        bodyrot, lsk, ck, sclk, ckmeta)

Contract: the snapshot with these fields, no dynamic or switch frame definitions, no
inner snapshot, no pool edits and no `boddef` pairs. `_dyn_levels` completes it.
"""
KernelSet(loaded, daf, spk, pck, pool, name2code, code2name, frames, framename2id, tk,
    bodyrot, lsk, ck, sclk, ckmeta) = KernelSet(loaded, daf, spk, pck, pool, name2code,
    code2name, frames, framename2id, tk, bodyrot, lsk, ck, sclk, ckmeta,
    Dict{Int32,DynFrame}(), Dict{Int32,SwitchFrame}(), nothing, false, Tuple{String,Int32}[],
    DSKFile[], Dict{Int32,Vector{DSKSegment}}(), EKFile[], EKTable[])

"""
    _dyn_levels(base, edited, bodydefs) -> KernelSet

Contract: the snapshot with the data of `base`, the parsed definitions of its class 5
and class 6 frames, and the flags `edited` and `bodydefs`. If a dynamic frame exists,
the result has a chain of `_DYN_LEVELS` inner snapshots (`inner.ks`) with the same data,
the innermost with `inner === nothing`; otherwise `inner === nothing`. No data are
copied.
"""
function _dyn_levels(base::KernelSet, edited::Bool, bodydefs::Vector{Tuple{String,Int32}})
    dyn = Dict{Int32,DynFrame}()
    sw = Dict{Int32,SwitchFrame}()
    for fd in values(base.frames)
        haskey(BUILTIN_FRAMES, fd.id) && continue
        fd.class == 5 && (dyn[fd.id] = _dyn_parse(base, fd))
        fd.class == 6 && (sw[fd.id] = _sw_parse(base, fd))
    end
    ks = nothing
    # invariant: `ks` is the snapshot of the level below the one that is built next
    for _ in 0:(isempty(dyn) ? 0 : _DYN_LEVELS)
        ks = KernelSet(base.loaded, base.daf, base.spk, base.pck, base.pool,
            base.name2code, base.code2name, base.frames, base.framename2id, base.tk,
            base.bodyrot, base.lsk, base.ck, base.sclk, base.ckmeta, dyn, sw,
            ks === nothing ? nothing : DynLevel{KernelSet}(ks), edited, bodydefs,
            base.dskfiles, base.dsk, base.ekfiles, base.ektables)
    end
    return ks::KernelSet
end

# ---------------------------------------------------------------------------------------
# Small generic helpers
# ---------------------------------------------------------------------------------------

"Vector of `Taylor{1,Float64}` with the values `r` and the derivatives `v`."
@inline _dyn_pack(r::Vec3{Float64}, v::Vec3{Float64}) = map((a, b) -> Taylor((a, b)), r, v)

"`x` with a zero derivative part (state mode) or `x` itself (plain numbers)."
@inline _dyn_still(x::Taylor{1,Float64}) = Taylor((x.v, 0.0))
@inline _dyn_still(x::Float64) = x

"The constant `x` in the number type `S`."
@inline _dyn_lift(::Type{S}, x::Float64) where {S<:Real} = convert(S, x)
@inline _dyn_lift(::Type{S}, x::Taylor{1,Float64}) where {S<:Real} = convert(S, _dyn_still(x))

"Unit vector along `a`, or the zero vector if `a` is zero, in any number type."
@inline function _dyn_unit(a::Vec3{S}) where {S<:Real}
    n = _norm3(a)
    return primal(n) == 0 ? zero(a) : a / n
end

"Unit vector along `a × b`, as `ucrss_` and `ducrss_`, in any number type."
@inline _dyn_ucross(a::Vec3{S}, b::Vec3{S}) where {S<:Real} = _dyn_unit(cross(a, b))

"`x - m * trunc(x / m)` with the arithmetic of `d_mod` of the f2c library; the subtracted quantity is a constant, so derivatives are not changed."
@inline _dyn_dmod(x::Real, m::Float64) = x - m * trunc(primal(x) / m)

"""
    _dyn_rot(angle, axis) -> Mat3

Contract: the coordinate rotation `[angle]_axis` of `rotate_` for `axis in 1:3`, in the
number type of `angle`.
"""
@inline function _dyn_rot(angle::S, axis::Integer) where {S<:Real}
    s, c = sincos(angle)
    o, z = one(S), zero(S)
    axis == 1 && return Mat3{S}(o, z, z, z, c, -s, z, s, c)
    axis == 2 && return Mat3{S}(c, z, s, z, o, z, -s, z, c)
    return Mat3{S}(c, -s, z, s, c, z, z, z, o)
end

"`[a3]_x3 * [a2]_x2 * [a1]_x1`, as `eul2m_` and (with `Taylor{1}` angles) `eul2xf_`."
@inline _dyn_eul(a3::S, a2::S, a1::S, x3::Integer, x2::Integer, x1::Integer) where {S<:Real} =
    _m_mxm(_dyn_rot(a3, x3), _m_mxm(_dyn_rot(a2, x2), _dyn_rot(a1, x1)))

"""
    _dyn_xf(mode, ks, from, to, et) -> Mat3

Contract: the matrix from frame `from` to frame `to` at `et`: `frame_rotation` in the
exact mode (`zzrefch0_`), `_ck_state_rotation` in the state mode (`zzfrmch0_`).
"""
@inline _dyn_xf(::_CKExact, ks::KernelSet, from::Int32, to::Int32, et::T) where {T<:Real} =
    frame_rotation(ks, from, to, et)::Mat3{T}
@inline _dyn_xf(::_CKState, ks::KernelSet, from::Int32, to::Int32, et::Float64) =
    _ck_state_rotation(ks, from, to, et)::Mat3{Taylor{1,Float64}}

# ---------------------------------------------------------------------------------------
# Earth precession, nutation and obliquity models
# ---------------------------------------------------------------------------------------

"Seconds in a Julian century, `jyear_() * 100`."
const _DYN_CENTURY = 31557600.0 * 100.0
"Radians per arcsecond, `rpd_() / 3600`."
const _DYN_ARCSEC = (pi / 180) / 3600.0

"""
    _dyn_prec(et) -> Mat3

Contract: the 1976 IAU precession matrix from J2000 to the mean equator and equinox
of `et` (TDB seconds), `[-z]_3 [theta]_2 [-zeta]_3` of `zzeprc76_`, in the number type
of `et`.
"""
@inline function _dyn_prec(et::S) where {S<:Real}
    t = et / _DYN_CENTURY
    zeta = t * (t * (t * 0.017998 + 0.30188) + 2306.2181) * _DYN_ARCSEC
    z = t * (t * (t * 0.018203 + 1.09468) + 2306.2181) * _DYN_ARCSEC
    theta = t * (t * (t * -0.041833 - 0.42665) + 2004.3109) * _DYN_ARCSEC
    return _dyn_eul(-z, theta, -zeta, 3, 2, 3)
end

"Contract: the mean obliquity of the ecliptic of date in radians (1980 IAU), as `zzmobliq_`."
@inline function _dyn_mobliq(et::S) where {S<:Real}
    t = et / 31557600.0 / 100.0
    return _DYN_ARCSEC * (t * (t * (t * 0.001813 - 5.9e-4) - 46.815) + 84381.448)
end

"""
    _DYN_WAHR

The 106 terms of the 1980 IAU nutation series as `zzwahr_` stores them: the five
multipliers of the fundamental arguments (l, l', F, D, Omega), then the longitude
coefficient and its rate, then the obliquity coefficient and its rate, in units of
1e-4 arcsecond and 1e-4 arcsecond per 1000 Julian years.
"""
const _DYN_WAHR = (
    (0, 0, 0, 0, 1, -171996, -1742, 92025, 89),
    (0, 0, 0, 0, 2, 2062, 2, -895, 5),
    (-2, 0, 2, 0, 1, 46, 0, -24, 0),
    (2, 0, -2, 0, 0, 11, 0, 0, 0),
    (-2, 0, 2, 0, 2, -3, 0, 1, 0),
    (1, -1, 0, -1, 0, -3, 0, 0, 0),
    (0, -2, 2, -2, 1, -2, 0, 1, 0),
    (2, 0, -2, 0, 1, 1, 0, 0, 0),
    (0, 0, 2, -2, 2, -13187, -16, 5736, -31),
    (0, 1, 0, 0, 0, 1426, -34, 54, -1),
    (0, 1, 2, -2, 2, -517, 12, 224, -6),
    (0, -1, 2, -2, 2, 217, -5, -95, 3),
    (0, 0, 2, -2, 1, 129, 1, -70, 0),
    (2, 0, 0, -2, 0, 48, 0, 1, 0),
    (0, 0, 2, -2, 0, -22, 0, 0, 0),
    (0, 2, 0, 0, 0, 17, -1, 0, 0),
    (0, 1, 0, 0, 1, -15, 0, 9, 0),
    (0, 2, 2, -2, 2, -16, 1, 7, 0),
    (0, -1, 0, 0, 1, -12, 0, 6, 0),
    (-2, 0, 0, 2, 1, -6, 0, 3, 0),
    (0, -1, 2, -2, 1, -5, 0, 3, 0),
    (2, 0, 0, -2, 1, 4, 0, -2, 0),
    (0, 1, 2, -2, 1, 4, 0, -2, 0),
    (1, 0, 0, -1, 0, -4, 0, 0, 0),
    (2, 1, 0, -2, 0, 1, 0, 0, 0),
    (0, 0, -2, 2, 1, 1, 0, 0, 0),
    (0, 1, -2, 2, 0, -1, 0, 0, 0),
    (0, 1, 0, 0, 2, 1, 0, 0, 0),
    (-1, 0, 0, 1, 1, 1, 0, 0, 0),
    (0, 1, 2, -2, 0, -1, 0, 0, 0),
    (0, 0, 2, 0, 2, -2274, -2, 977, -5),
    (1, 0, 0, 0, 0, 712, 1, -7, 0),
    (0, 0, 2, 0, 1, -386, -4, 200, 0),
    (1, 0, 2, 0, 2, -301, 0, 129, -1),
    (1, 0, 0, -2, 0, -158, 0, -1, 0),
    (-1, 0, 2, 0, 2, 123, 0, -53, 0),
    (0, 0, 0, 2, 0, 63, 0, -2, 0),
    (1, 0, 0, 0, 1, 63, 1, -33, 0),
    (-1, 0, 0, 0, 1, -58, -1, 32, 0),
    (-1, 0, 2, 2, 2, -59, 0, 26, 0),
    (1, 0, 2, 0, 1, -51, 0, 27, 0),
    (0, 0, 2, 2, 2, -38, 0, 16, 0),
    (2, 0, 0, 0, 0, 29, 0, -1, 0),
    (1, 0, 2, -2, 2, 29, 0, -12, 0),
    (2, 0, 2, 0, 2, -31, 0, 13, 0),
    (0, 0, 2, 0, 0, 26, 0, -1, 0),
    (-1, 0, 2, 0, 1, 21, 0, -10, 0),
    (-1, 0, 0, 2, 1, 16, 0, -8, 0),
    (1, 0, 0, -2, 1, -13, 0, 7, 0),
    (-1, 0, 2, 2, 1, -10, 0, 5, 0),
    (1, 1, 0, -2, 0, -7, 0, 0, 0),
    (0, 1, 2, 0, 2, 7, 0, -3, 0),
    (0, -1, 2, 0, 2, -7, 0, 3, 0),
    (1, 0, 2, 2, 2, -8, 0, 3, 0),
    (1, 0, 0, 2, 0, 6, 0, 0, 0),
    (2, 0, 2, -2, 2, 6, 0, -3, 0),
    (0, 0, 0, 2, 1, -6, 0, 3, 0),
    (0, 0, 2, 2, 1, -7, 0, 3, 0),
    (1, 0, 2, -2, 1, 6, 0, -3, 0),
    (0, 0, 0, -2, 1, -5, 0, 3, 0),
    (1, -1, 0, 0, 0, 5, 0, 0, 0),
    (2, 0, 2, 0, 1, -5, 0, 3, 0),
    (0, 1, 0, -2, 0, -4, 0, 0, 0),
    (1, 0, -2, 0, 0, 4, 0, 0, 0),
    (0, 0, 0, 1, 0, -4, 0, 0, 0),
    (1, 1, 0, 0, 0, -3, 0, 0, 0),
    (1, 0, 2, 0, 0, 3, 0, 0, 0),
    (1, -1, 2, 0, 2, -3, 0, 1, 0),
    (-1, -1, 2, 2, 2, -3, 0, 1, 0),
    (-2, 0, 0, 0, 1, -2, 0, 1, 0),
    (3, 0, 2, 0, 2, -3, 0, 1, 0),
    (0, -1, 2, 2, 2, -3, 0, 1, 0),
    (1, 1, 2, 0, 2, 2, 0, -1, 0),
    (-1, 0, 2, -2, 1, -2, 0, 1, 0),
    (2, 0, 0, 0, 1, 2, 0, -1, 0),
    (1, 0, 0, 0, 2, -2, 0, 1, 0),
    (3, 0, 0, 0, 0, 2, 0, 0, 0),
    (0, 0, 2, 1, 2, 2, 0, -1, 0),
    (-1, 0, 0, 0, 2, 1, 0, -1, 0),
    (1, 0, 0, -4, 0, -1, 0, 0, 0),
    (-2, 0, 2, 2, 2, 1, 0, -1, 0),
    (-1, 0, 2, 4, 2, -2, 0, 1, 0),
    (2, 0, 0, -4, 0, -1, 0, 0, 0),
    (1, 1, 2, -2, 2, 1, 0, -1, 0),
    (1, 0, 2, 2, 1, -1, 0, 1, 0),
    (-2, 0, 2, 4, 2, -1, 0, 1, 0),
    (-1, 0, 4, 0, 2, 1, 0, 0, 0),
    (1, -1, 0, -2, 0, 1, 0, 0, 0),
    (2, 0, 2, -2, 1, 1, 0, -1, 0),
    (2, 0, 2, 2, 2, -1, 0, 0, 0),
    (1, 0, 0, 2, 1, -1, 0, 0, 0),
    (0, 0, 4, -2, 2, 1, 0, 0, 0),
    (3, 0, 2, -2, 2, 1, 0, 0, 0),
    (1, 0, 2, -2, 0, -1, 0, 0, 0),
    (0, 1, 2, 0, 1, 1, 0, 0, 0),
    (-1, -1, 0, 2, 1, 1, 0, 0, 0),
    (0, 0, -2, 0, 1, -1, 0, 0, 0),
    (0, 0, 2, -1, 2, -1, 0, 0, 0),
    (0, 1, 0, 2, 0, -1, 0, 0, 0),
    (1, 0, -2, -2, 0, -1, 0, 0, 0),
    (0, -1, 2, 0, 1, -1, 0, 0, 0),
    (1, 1, 0, -2, 1, -1, 0, 0, 0),
    (1, 0, -2, 2, 0, -1, 0, 0, 0),
    (2, 0, 0, 2, 0, 1, 0, 0, 0),
    (0, 0, 2, 4, 2, -1, 0, 0, 0),
    (0, 1, 0, 1, 0, 1, 0, 0, 0),
)

"Fundamental arguments of `zzwahr_`: constant (degrees), rate (degrees per day), and the two higher coefficients, for l, l', F, D and Omega."
const _DYN_WAHR_ARGS = (
    (134.96298138888886, 13.064992947243136, 6.5192872572139397e-4, 3.6484365631332527e-7),
    (357.52772333333331, 0.98560028309377146, -1.201414483363923e-5, -6.8408185558748495e-8),
    (93.271910277777778, 13.229350240603848, -2.760338267929901e-4, 6.2707503428852773e-8),
    (297.85036305555559, 12.190749116495549, -1.4348262053484912e-4, 1.0831296046801845e-7),
    (125.04452222222223, -0.052953764841432813, 1.5522608272925558e-4, 4.5605457039165659e-8),
)

"Degrees per radian as `zzwahr_` computes it."
const _DYN_RADIAN = 180.0 / pi
"Divisor that converts the units of the nutation series to radians."
const _DYN_WAHR_FACTOR = _DYN_RADIAN * 3600.0 * 1e4

"""
    _dyn_wahr(et) -> (dpsi, deps)

Contract: nutation in longitude and in obliquity (radians) of the 1980 IAU theory at
`et` (TDB seconds), with the operations of `zzwahr_`, in the number type of `et`. With
a `Taylor` epoch the derivatives are those of the series, including the slow change
of the coefficients that the rates of `zzwahr_` leave out.
"""
function _dyn_wahr(et::S) where {S<:Real}
    dj = et / 86400.0
    dd = dj / 1e4
    dd2 = dd * dd
    t = dj / 365250.0
    ang = map(_DYN_WAHR_ARGS) do c
        _dyn_dmod(c[1] + dj * c[2] + (c[3] + dd * c[4]) * dd2, 360.0) / _DYN_RADIAN
    end
    dpsi = zero(S)
    deps = zero(S)
    for i in 1:106      # invariant: dpsi, deps hold the sums of the terms 1:i-1
        m = @inbounds _DYN_WAHR[i]
        arg = zero(S)
        for j in 1:5
            if m[j] != 0
                arg += m[j] * ang[j]
                arg = _dyn_dmod(arg, _TWOPI)
            end
        end
        cl = zero(S) + m[6]
        m[7] != 0 && (cl += m[7] * t)
        ce = zero(S) + m[8]
        m[9] != 0 && (ce += m[9] * t)
        s, c = sincos(arg)
        dpsi += cl * s / _DYN_WAHR_FACTOR
        deps += ce * c / _DYN_WAHR_FACTOR
    end
    return dpsi, deps
end

"""
    _dyn_wahr_state(et::Float64) -> (dpsi, deps)

Contract: the nutation angles of `_dyn_wahr` as `Taylor{1,Float64}` values whose
derivative parts are the rates that `zzwahr_` returns (radians per second): the rates
of the arguments, without the change of the coefficients.
"""
function _dyn_wahr_state(et::Float64)
    dj = et / 86400.0
    dd = dj / 1e4
    dddj = dd / 1e4
    dd2 = dd * dd
    t = dj / 365250.0
    ang = map(_DYN_WAHR_ARGS) do c
        _dyn_dmod(c[1] + dj * c[2] + (c[3] + dd * c[4]) * dd2, 360.0) / _DYN_RADIAN
    end
    rate = map(_DYN_WAHR_ARGS) do c
        _dyn_dmod(c[2] + dddj * (c[3] * 2.0 + dd * 3.0 * c[4]), 360.0) / _DYN_RADIAN
    end
    dpsi = deps = ddpsi = ddeps = 0.0
    for i in 1:106      # invariant: the four sums hold the terms 1:i-1
        m = @inbounds _DYN_WAHR[i]
        arg = 0.0
        argrt = 0.0
        for j in 1:5
            if m[j] != 0
                arg += m[j] * ang[j]
                argrt += m[j] * rate[j]
                arg = _dyn_dmod(arg, _TWOPI)
            end
        end
        cl = Float64(m[6])
        m[7] != 0 && (cl += m[7] * t)
        ce = Float64(m[8])
        m[9] != 0 && (ce += m[9] * t)
        s, c = sincos(arg)
        dpsi += cl * s / _DYN_WAHR_FACTOR
        deps += ce * c / _DYN_WAHR_FACTOR
        ddpsi += cl * c * argrt / _DYN_WAHR_FACTOR
        ddeps -= ce * s * argrt / _DYN_WAHR_FACTOR
    end
    return Taylor((dpsi, ddpsi / 86400.0)), Taylor((deps, ddeps / 86400.0))
end

"Nutation angles in the number type of the mode: `_dyn_wahr` or `_dyn_wahr_state`."
@inline _dyn_nutation(::_CKExact, et::Real) = _dyn_wahr(et)
@inline _dyn_nutation(::_CKState, et::Float64) = _dyn_wahr_state(et)

"""
    _dyn_ofdate(mode, ks, def, et) -> Mat3

Contract: the matrix from the of-date frame `def` (families 1 to 3) to its base frame
at `et`: the transpose of the precession matrix, of nutation times precession, or of
the obliquity rotation times precession, then the rotation from J2000 to the base.
"""
function _dyn_ofdate(mode, ks::KernelSet, def::DynFrame, et::T) where {T<:Real}
    tv = _ck_epoch(mode, et)
    m = _dyn_prec(tv)
    if def.family == _DYN_TRUE_EQ
        dpsi, deps = _dyn_nutation(mode, et)
        mob = _dyn_mobliq(tv)
        m = _m_mxm(_dyn_eul(-mob - deps, -dpsi, mob, 1, 3, 1), m)
    elseif def.family == _DYN_MEAN_ECL
        m = _m_mxm(_dyn_rot(_dyn_mobliq(tv), 1), m)
    end
    x = typeof(m)(transpose(m))
    def.base == _J2000_FRAME && return x
    return _m_mxm(_dyn_xf(mode, ks, _J2000_FRAME, def.base, et), x)
end

# ---------------------------------------------------------------------------------------
# States and positions of the vector definitions
# ---------------------------------------------------------------------------------------

"""
    _dyn_stelab(p, v) -> Vec3

Contract: the position `p` with the stellar aberration for the observer velocity `v`,
`cos(phi) p + |p| v_perp / c`, the rotation of `stelab_` in closed form (no inverse
sine), for any number type. Requires `p != 0` and `|v| < c`.
"""
@inline function _dyn_stelab(p, v)
    pm = _norm3(p)
    u = p / pm
    vp = v - _dot3(v, u) * u
    c = sqrt(1 - _dot3(vp, vp) / CLIGHT^2)
    return c * p + (pm / CLIGHT) * vp
end

"""
    _dyn_apo(ks, targ, et, ref, corr, obs) -> (r, lt)

Contract: position of `targ` relative to `obs` in frame `ref` with the correction
`corr`, and the light time, as `zzspkzp0_` (`spkpos`). For a `Float64` epoch the values
are those of `spkpos`. For another number type the same algorithm runs on that type
(one light-time iteration for `LT`, three for `CN`, stellar aberration in closed form),
so the derivatives are those of the position that `spkpos` returns. No allocation.
"""
function _dyn_apo(ks::KernelSet, targ::Int32, et::T, ref::Int32, corr,
        obs::Int32) where {T<:Real}
    # The body relative to itself: zero, and no square root of a zero series.
    targ == obs && return zero(Vec3{T}), zero(T)
    if corr.geometric
        r = spk_jet(ks, targ, et, ref, obs, Val(0))[1]::Vec3{T}
        return r, _norm3(r) / CLIGHT
    end
    T === Float64 && return abcorr_position(ks, targ, et, ref, corr, obs)::Tuple{Vec3{T},T}
    inertial = _is_inertial(ks, ref)
    iref = inertial ? ref : _J2000_FRAME
    nos = _aberration(corr.flags & ~_AB_S)
    robs, vobs = _ssb(ks, obs, et, iref)::NTuple{2,Vec3{T}}
    r, lt = _apo(ks, targ, et, iref, nos, robs, vobs)::Tuple{Vec3{T},T}
    if corr.stellar && primal(_dot3(r, r)) > 0
        r = _dyn_stelab(r, corr.transmit ? -vobs : +vobs)
    end
    inertial && return r, lt
    center = frame_def(ks, ref).center
    ltc = center == obs ? zero(lt) : center == targ ? lt :
          _apo(ks, center, et, iref, nos, robs, vobs)[2]::T
    return _m_mxv(_dyn_xf(_CKExact(), ks, _J2000_FRAME, ref, et + _ltsign(corr) * ltc), r), lt
end

"""
    _dyn_spkez(ks, targ, et, ref, corr, obs) -> (r, v, lt)

Contract: position, velocity and light time of `targ` relative to `obs` in frame `ref`
with the correction `corr`, as `zzspkez0_` (`spkezr`), in the number type of `et`. The
velocity is the one of CSPICE (see `spkder` for its difference from the derivative of
the position). No allocation.
"""
function _dyn_spkez(ks::KernelSet, targ::Int32, et::T, ref::Int32, corr,
        obs::Int32) where {T<:Real}
    if corr.geometric
        r, v = spk_jet(ks, targ, et, ref, obs, Val(1))::NTuple{2,Vec3{T}}
        return r, v, _norm3(r) / CLIGHT
    end
    s, lt = abcorr_state(ks, targ, et, ref, corr, obs)::Tuple{SVector{6,T},T}
    return Vec3{T}(s[1], s[2], s[3]), Vec3{T}(s[4], s[5], s[6]), lt
end

"""
    _dyn_nearpt(p, radii) -> Vec3

Contract: the point of the ellipsoid with the semi-axes `radii` (all positive) that is
nearest to the point `p`, in the number type of `p`: the near point of `nearpt`
(`_m_nearpt`), so the plain value is the one of CSPICE and each derivative order comes
from the stationarity condition. Postcondition: the result is on the ellipsoid to
rounding.
"""
@inline _dyn_nearpt(p::Vec3{S}, radii::Vec3{Float64}) where {S<:Real} =
    Vec3{S}(_m_nearpt(p, radii[1], radii[2], radii[3])[1])

"Epoch `et` moved by the light time `lt` in the direction of the correction `corr`, as `zzcorepc_`."
@inline _dyn_corepc(corr, et, lt) = corr.uselt ? et + _ltsign(corr) * lt : et

"""
    _dyn_vector(mode, ks, vd, et) -> Vec3

Contract: the defining vector `vd` of a two-vector frame in J2000 at `et`, before the
sign of the axis is applied.

Exact mode (`zzdynrot_`), in the number type of `et`: the position of `spkpos`; the
unit vector of the velocity of `spkezr`, rotated from its frame at the epoch corrected
for the light time to the center of that frame; the vector from the observer to the
near point on the target ellipsoid, rotated from the body-fixed frame at the epoch
corrected for the light time; the constant vector rotated from its frame, with the
optional stellar aberration.

State mode (`zzdynfrm_`), `Taylor{1,Float64}` values: the state of `spkezr`; the unit
velocity and the central difference of the unit velocity over `±max(1, et * 2^-27)` s;
the state of the near point for the state of `spkezr`; the constant vector with zero
derivative before the rotation. A stellar aberration correction of a constant vector
changes the vector and leaves the derivative as it is, as in CSPICE.
"""
function _dyn_vector(mode::_CKExact, ks::KernelSet, vd::DynVector, et::T) where {T<:Real}
    corr = _aberration(vd.corr)
    if vd.kind == 1
        return _dyn_apo(ks, vd.target, et, _J2000_FRAME, corr, vd.observer)[1]::Vec3{T}
    elseif vd.kind == 2
        _, v, _ = _dyn_spkez(ks, vd.target, et, vd.frame, corr, vd.observer)
        u = _dyn_unit(v)
        vd.frame == _J2000_FRAME && return u::Vec3{T}
        vet = et
        if corr.uselt && vd.frclass != 1
            vet = _dyn_corepc(corr, et, _dyn_apo(ks, vd.frcenter, et, _J2000_FRAME, corr,
                vd.observer)[2])
        end
        return _m_mxv(_dyn_xf(_CKExact(), ks, vd.frame, _J2000_FRAME, vet), u)::Vec3{T}
    elseif vd.kind == 3
        p, lt = _dyn_apo(ks, vd.target, et, vd.frame, corr, vd.observer)
        pobs = -p
        d = _dyn_nearpt(pobs, vd.vec) - pobs
        fet = _dyn_corepc(corr, et, lt)
        return _m_mxv(_dyn_xf(_CKExact(), ks, vd.frame, _J2000_FRAME, fet), d)::Vec3{T}
    end
    # A constant vector.
    fet = et
    if corr.uselt && !corr.stellar && vd.frclass != 1
        fet = _dyn_corepc(corr, et, _dyn_apo(ks, vd.frcenter, et, _J2000_FRAME, corr,
            vd.observer)[2])
    end
    c = Vec3{T}(vd.vec)
    vd.frame == _J2000_FRAME || (c = _m_mxv(_dyn_xf(_CKExact(), ks, vd.frame, _J2000_FRAME, fet), c))
    if corr.stellar
        vobs = spk_jet(ks, vd.observer, et, _J2000_FRAME, Int32(0), Val(1))[2]
        c = _dyn_stelab(c, corr.transmit ? -vobs : +vobs)
    end
    return c::Vec3{T}
end

function _dyn_vector(mode::_CKState, ks::KernelSet, vd::DynVector, et::Float64)
    corr = _aberration(vd.corr)
    if vd.kind == 1
        r, v, _ = _dyn_spkez(ks, vd.target, et, _J2000_FRAME, corr, vd.observer)
        return _dyn_pack(r, v)
    elseif vd.kind == 2
        delta = max(1.0, et * 7.4505805969238281e-9)
        um = _dyn_unit(_dyn_spkez(ks, vd.target, et - delta, vd.frame, corr, vd.observer)[2])
        u0 = _dyn_unit(_dyn_spkez(ks, vd.target, et, vd.frame, corr, vd.observer)[2])
        up = _dyn_unit(_dyn_spkez(ks, vd.target, et + delta, vd.frame, corr, vd.observer)[2])
        s = _dyn_pack(u0, (0.5 / delta) * up + (-0.5 / delta) * um)
        vd.frame == _J2000_FRAME && return s
        vet = et
        if corr.uselt && vd.frclass != 1
            vet = _dyn_corepc(corr, et, _dyn_apo(ks, vd.frcenter, et, _J2000_FRAME, corr,
                vd.observer)[2])
        end
        return _m_mxv(_dyn_xf(_CKState(), ks, vd.frame, _J2000_FRAME, vet), s)
    elseif vd.kind == 3
        r, v, lt = _dyn_spkez(ks, vd.target, et, vd.frame, corr, vd.observer)
        sobs = _dyn_pack(-r, -v)
        d = _dyn_nearpt(sobs, vd.vec) - sobs
        fet = _dyn_corepc(corr, et, lt)
        return _m_mxv(_dyn_xf(_CKState(), ks, vd.frame, _J2000_FRAME, fet), d)
    end
    fet = et
    if corr.uselt && !corr.stellar && vd.frclass != 1
        fet = _dyn_corepc(corr, et, _dyn_apo(ks, vd.frcenter, et, _J2000_FRAME, corr,
            vd.observer)[2])
    end
    c = _dyn_pack(vd.vec, zero(Vec3{Float64}))
    vd.frame == _J2000_FRAME || (c = _m_mxv(_dyn_xf(_CKState(), ks, vd.frame, _J2000_FRAME, fet), c))
    if corr.stellar
        vobs = spk_jet(ks, vd.observer, et, _J2000_FRAME, Int32(0), Val(1))[2]
        p = _stelab(map(x -> x.v, c), corr.transmit ? -vobs : +vobs)
        c = map((a, b) -> Taylor((a, b.d[1])), p, c)
    end
    return c
end

@noinline function _dyn_throw_sep(def::DynFrame, sep, et)
    throw(SpiceError(:degeneratecase, "angular separation of the vectors that define " *
        "the two-vector frame $(def.name) is $sep radians at ET $(primal(et)); the " *
        "minimum allowed difference from 0 or pi is $(def.minsep) radians " *
        "(keyword ANGLE_SEP_TOL)"))
end

"Angular separation of the non-zero vectors `a` and `b` in radians, as `vsep_`."
function _dyn_vsep(a::Vec3{Float64}, b::Vec3{Float64})
    na, nb = norm(a), norm(b)
    (na == 0.0 || nb == 0.0) && return 0.0
    u, v = a / na, b / nb
    d = dot(u, v)
    d > 0 && return 2.0 * asin(0.5 * norm(u - v))
    d < 0 && return pi - 2.0 * asin(0.5 * norm(u + v))
    return pi / 2
end

"""
    _dyn_twovec(mode, ks, def, et) -> Mat3

Contract: the matrix from the two-vector frame `def` to its base frame at `et`. The
columns are the frame axes in J2000 as `zztwovxf_` and `twovec_` build them: the axis
of the primary vector along that vector, the axis of the secondary vector in the
half-plane of the two vectors. Throws `SpiceError(:degeneratecase)` if the sine of the
separation of the vectors is below the sine of `def.minsep`.
"""
function _dyn_twovec(mode, ks::KernelSet, def::DynFrame, et::T) where {T<:Real}
    a = _dyn_vector(mode, ks, def.v1, et)
    b = _dyn_vector(mode, ks, def.v2, et)
    def.v1.negate && (a = -a)
    def.v2.negate && (b = -b)
    sep = _dyn_vsep(map(primal, a), map(primal, b))
    sin(sep) < sin(def.minsep) && _dyn_throw_sep(def, sep, et)
    i1 = Int(def.v1.axis)
    i2 = i1 % 3 + 1
    c1 = _dyn_unit(a)
    if def.v2.axis == i2
        c3 = _dyn_ucross(a, b)
        c2 = _dyn_ucross(c3, a)
    else
        c2 = _dyn_ucross(b, a)
        c3 = _dyn_ucross(a, c2)
    end
    x = i1 == 1 ? hcat(c1, c2, c3) : i1 == 2 ? hcat(c3, c1, c2) : hcat(c2, c3, c1)
    def.base == _J2000_FRAME && return x
    return _m_mxm(_dyn_xf(mode, ks, _J2000_FRAME, def.base, et), x)
end

"""
    _dyn_euler(mode, def, et) -> Mat3

Contract: the matrix from the Euler frame `def` to its base frame at `et`,
`[a1]_x1 [a2]_x2 [a3]_x3` with the angles `a_i` polynomials in `et - def.epoch`.
"""
function _dyn_euler(mode, def::DynFrame, et::T) where {T<:Real}
    delta = _ck_epoch(mode, et) - def.epoch
    a1 = _horner(def.coef1, delta) * def.scale / _DYN_RADIAN
    a2 = _horner(def.coef2, delta) * def.scale / _DYN_RADIAN
    a3 = _horner(def.coef3, delta) * def.scale / _DYN_RADIAN
    return _dyn_eul(a1, a2, a3, def.axes[1], def.axes[2], def.axes[3])
end

"""
    _dyn_product(mode, ks, def, et) -> Mat3

Contract: the matrix from the product frame `def` to its base frame at `et`: the
product, last factor on the left, of the matrices from `def.to[i]` to `def.from[i]`.
"""
function _dyn_product(mode, ks::KernelSet, def::DynFrame, et::T) where {T<:Real}
    x = _dyn_xf(mode, ks, (@inbounds def.to[1]), (@inbounds def.from[1]), et)
    for i in 2:length(def.to)       # invariant: x is the product of the factors 1:i-1
        x = _m_mxm(_dyn_xf(mode, ks, (@inbounds def.to[i]), (@inbounds def.from[i]), et), x)
    end
    return x
end

"""
    _dyn_inert(mode, ks, def, et, x) -> Mat3

Contract: `x` with the derivative part that `ROTATION_STATE = 'INERTIAL'` defines, as
`zzdynfrm_`: zero relative to J2000, so zero if the base is a built-in inertial frame
and otherwise the derivative of the rotation from J2000 to the base, applied to the
fixed orientation. The exact mode has no derivative part to replace and returns `x`.
"""
@inline _dyn_inert(::_CKExact, ks::KernelSet, def::DynFrame, et, x) = x
function _dyn_inert(::_CKState, ks::KernelSet, def::DynFrame, et::Float64,
        x::Mat3{Taylor{1,Float64}})
    def.inertial || return x
    def.base_inertial && return map(_dyn_still, x)
    b = _dyn_xf(_CKState(), ks, def.base, _J2000_FRAME, et)
    out = Mat3{Taylor{1,Float64}}(transpose(b)) * map(_dyn_still, b * x)
    return map((a, o) -> Taylor((a.v, o.d[1])), x, out)
end

"""
    _dyn_eval(mode, ks, def, et) -> Mat3

Contract: the matrix from the dynamic frame `def` to its base frame at `et`, evaluated
with the states and frames of `ks` (the inner snapshot of the one that holds `def`).
Requires `def.errkind === :none`. No allocation.
"""
function _dyn_eval(mode, ks::KernelSet, def::DynFrame, et::T) where {T<:Real}
    S = _ck_numtype(mode, T)
    fam = def.family
    x = fam <= _DYN_MEAN_ECL ? _dyn_ofdate(mode, ks, def, et)::Mat3{S} :
        fam == _DYN_TWOVEC ? _dyn_twovec(mode, ks, def, et)::Mat3{S} :
        fam == _DYN_EULER ? _dyn_euler(mode, def, et)::Mat3{S} :
        _dyn_product(mode, ks, def, et)::Mat3{S}
    return _dyn_inert(mode, ks, def, et, x)::Mat3{S}
end

# ---------------------------------------------------------------------------------------
# The recursion barrier
# ---------------------------------------------------------------------------------------
#
# The evaluation of a dynamic frame calls the state and frame functions of the package,
# and these call the evaluation of a dynamic frame again when a definition names one.
# CSPICE limits this recursion to two levels at run time, but the compiler sees a cycle
# without a limit in which the number types can nest (a series of series), and type
# inference then gives up on some calls of the cycle, also for kernels without any
# dynamic frame. The call from `_dyn_frame_step` to `_dyn_eval` therefore goes through
# a typed C function pointer: the argument and result types are exact, nothing is boxed,
# and the compiler does not follow the call. There is one pointer for each mode and
# number type, made on first use.

"Pointers of `_dyn_call`, by `Tuple{mode type, number type}`. Replaced, never mutated."
const _DYN_POINTERS = Ref(IdDict{Any,Ptr{Cvoid}}())
const _DYN_POINTER_LOCK = ReentrantLock()

"C entry points of `_dyn_call`: `_dyn_eval` with the snapshot of the level `lev`."
_dyn_entry_exact(lev::DynLevel{KernelSet}, def::DynFrame, et::T) where {T<:Real} =
    _dyn_eval(_CKExact(), lev.ks, def, et)
_dyn_entry_state(lev::DynLevel{KernelSet}, def::DynFrame, et::Float64) =
    _dyn_eval(_CKState(), lev.ks, def, et)

"Make the C function pointer of `_dyn_eval` for the mode and number type of `F`."
@generated function _dyn_cfunction(::Type{Tuple{_CKExact,T}}) where {T}
    return :(@cfunction(_dyn_entry_exact, $(Mat3{T}),
        (Ref{DynLevel{KernelSet}}, Ref{DynFrame}, $T)))
end
_dyn_cfunction(::Type{Tuple{_CKState,Float64}}) = @cfunction(_dyn_entry_state,
    Mat3{Taylor{1,Float64}}, (Ref{DynLevel{KernelSet}}, Ref{DynFrame}, Float64))

"""
    _dyn_pointer_slow(F) -> Ptr{Cvoid}

Contract: make the pointer for `F`, publish a new table that has it, and return it. The
table is replaced under a lock, so a reader without the lock sees the old table or the
new one, both complete. The compiler does not follow the call that makes the pointer.
"""
@noinline function _dyn_pointer_slow(@nospecialize(F))
    return @lock _DYN_POINTER_LOCK begin
        p = get(_DYN_POINTERS[], F, C_NULL)::Ptr{Cvoid}
        if p == C_NULL
            p = Base.invokelatest(_dyn_cfunction, F)::Ptr{Cvoid}
            table = copy(_DYN_POINTERS[])
            table[F] = p
            _DYN_POINTERS[] = table
        end
        p
    end
end

"Pointer of `_dyn_call` for `F = Tuple{mode type, number type}`. No allocation after the first use."
@inline function _dyn_pointer(::Type{F}) where {F}
    p = get(_DYN_POINTERS[], F, C_NULL)::Ptr{Cvoid}
    return p == C_NULL ? _dyn_pointer_slow(F) : p
end

"Forget the pointers of an earlier process (they are not valid in a loaded package image)."
function __init__()
    _DYN_POINTERS[] = IdDict{Any,Ptr{Cvoid}}()
    return nothing
end

"""
    _dyn_call(mode, lev, def, et) -> Mat3

Contract: `_dyn_eval(mode, lev.ks, def, et)`, called through the typed function pointer
of the mode and of the number type of `et` (see "The recursion barrier"). A number type
that is not plain bits is evaluated with a direct call. No allocation.
"""
@generated function _dyn_call(::_CKExact, lev::DynLevel{KernelSet}, def::DynFrame,
        et::T) where {T<:Real}
    isbitstype(T) || return :(_dyn_eval(_CKExact(), lev.ks, def, et))
    return :(ccall(_dyn_pointer(Tuple{_CKExact,$T}), $(Mat3{T}),
        (Ref{DynLevel{KernelSet}}, Ref{DynFrame}, $T), lev, def, et))
end
@inline _dyn_call(::_CKState, lev::DynLevel{KernelSet}, def::DynFrame, et::Float64) =
    ccall(_dyn_pointer(Tuple{_CKState,Float64}), Mat3{Taylor{1,Float64}},
        (Ref{DynLevel{KernelSet}}, Ref{DynFrame}, Float64), lev, def, et)

@noinline _dyn_throw_def(kind::Symbol, msg::String) = throw(SpiceError(kind, msg))
@noinline _dyn_throw_missing(id) = throw(SpiceError(:framedatanotfound,
    "no definition of the dynamic frame $id is loaded"))
@noinline _dyn_throw_depth(def::DynFrame) = throw(SpiceError(:recursiontoodeep,
    "the dynamic frame $(def.name) is needed for the evaluation of a dynamic frame " *
    "that is itself needed for the evaluation of a dynamic frame; dynamic frames may " *
    "not be used at recursion level 1"))

"Contract: returns the matrix `x` with each element in the number type `S` (`_dyn_lift`). The type is an argument with a type parameter: a closure that captures a type is not compiled for it in Julia 1.10, and then allocates."
@inline _dyn_liftmat(::Type{S}, x::Mat3) where {S} = Mat3{S}(ntuple(i -> _dyn_lift(S, x[i]), Val(9)))

"""
    _dyn_frame_step(mode, ks, fd, et) -> (base::Int32, rot)

One link of a frame chain for the class 5 frame `fd`: its base frame (`RELATIVE`) and
the matrix from the frame to the base at `et`, in the number type of the mode. A
frozen frame is evaluated at its `FREEZE_EPOCH` and is constant. Throws the error of
an invalid definition, and `SpiceError(:recursiontoodeep)` at the recursion limit.
No allocation. Not inlined, so that the frame functions that hold the call stay small
for the frames of the other classes.
"""
@noinline function _dyn_frame_step(mode, ks::KernelSet, fd::FrameDef, et::T) where {T<:Real}
    def = get(ks.dyn, fd.id, nothing)
    def === nothing && _dyn_throw_missing(fd.id)
    def.errkind === :none || _dyn_throw_def(def.errkind, def.errmsg)
    inner = ks.inner
    inner === nothing && _dyn_throw_depth(def)
    S = _ck_numtype(mode, T)
    if def.frozen
        x = _dyn_call(mode, inner, def, def.t0)
        return def.base, _dyn_liftmat(S, x)
    end
    return def.base, _dyn_call(mode, inner, def, et)::Mat3{S}
end

"""
    _dyn_base(ks, id) -> Int32

Contract: the base frame of the class 5 frame `id`, without any evaluation. Throws as
`_dyn_frame_step` does for a missing or invalid definition.
"""
function _dyn_base(ks::KernelSet, id::Int32)
    def = get(ks.dyn, id, nothing)
    def === nothing && _dyn_throw_missing(id)
    def.errkind === :none || _dyn_throw_def(def.errkind, def.errmsg)
    return def.base
end

"""
    _dyn_chain_step(mode, ks, fd, et) -> (parent::Int32, rot)

One link of a frame chain for a frame of class 5 (`_dyn_frame_step`) or class 6
(`_sw_frame_step`). Throws `SpiceError(:unknownframetype)` for a higher class.
"""
@noinline function _dyn_chain_step(mode, ks::KernelSet, fd::FrameDef, et::T) where {T<:Real}
    fd.class == 5 && return _dyn_frame_step(mode, ks, fd, et)
    fd.class == 6 && return _sw_frame_step(mode, ks, fd, et)
    throw(SpiceError(:unknownframetype,
        "frame $(fd.id) has class $(fd.class), which is not supported"))
end
