# Geometry layer, shared part: the computation method strings (`zzprsmet`), the
# aberration correction with the error kinds of `zzvalcor`, the lookup of bodies, frames
# and radii with the error kinds of the geometry calls of CSPICE, and the target shape
# (`zzsinutl`: `zzsuelin`, `zzsudski`, `zzraysfx`, `zzraynp`, `zzmaxrad`).
#
# Every function is a pure function of the `KernelSet` snapshot and of its arguments.
# The private names of the geometry layer start with `_geo_`, `_Geo` or `_GEO_`.

# --- method strings ---------------------------------------------------------------------

const _GEO_NEARPOINT = 0x01     # sub-point type "NEAR POINT" (legacy ellipsoid forms)
const _GEO_NADIR = 0x02         # sub-point type "NADIR"
const _GEO_INTERCEPT = 0x03     # sub-point type "INTERCEPT"
const _GEO_TANGENT = 0x01       # limb or terminator point definition "TANGENT"
const _GEO_GUIDED = 0x02        # limb or terminator point definition "GUIDED"
const _GEO_UMBRAL = 0x01        # terminator type "UMBRAL"
const _GEO_PENUMBRAL = 0x02     # terminator type "PENUMBRAL"

"The surface list of a method without the `SURFACES` keyword (never modified)."
const _GEO_NOSURF = Int[]

"""
    _GeoMethod

A parsed method string, the outputs of `zzprsmet`: `dsk` (the shape is `DSK`, else
`ELLIPSOID`), `subtyp` (0 or one of `_GEO_NEARPOINT`, `_GEO_NADIR`, `_GEO_INTERCEPT`),
`pntdef` (0, `_GEO_TANGENT`, `_GEO_GUIDED`), `trmtyp` (0, `_GEO_UMBRAL`,
`_GEO_PENUMBRAL`) and the surface ID list `srflst` (empty: all surfaces).
"""
struct _GeoMethod
    dsk::Bool
    subtyp::UInt8
    pntdef::UInt8
    trmtyp::UInt8
    srflst::Vector{Int}
end

@inline _geo_ell(subtyp::UInt8, pntdef::UInt8=0x00, trmtyp::UInt8=0x00) =
    _GeoMethod(false, subtyp, pntdef, trmtyp, _GEO_NOSURF)
@inline _geo_dsk(subtyp::UInt8, pntdef::UInt8=0x00, trmtyp::UInt8=0x00) =
    _GeoMethod(true, subtyp, pntdef, trmtyp, _GEO_NOSURF)

"""
    _geo_method(ks, bodyid, method) -> _GeoMethod

`zzprsmet`: parse the method string `method` for the body `bodyid`.

Contract: the common spellings (upper case, one blank at most between words, no
surface list) are compared as they are and do not allocate; every other string goes
through the lexer and the parser of `surfaces.jl`, which accept and reject what
`zzprsmet` accepts and rejects. Throws `SpiceError` (`:badmethodsyntax`,
`:badpriorityspec`, `:idcodenotfound`, `:toomanysurfaces`, `:syntaxerror`).
"""
@inline function _geo_method(ks::KernelSet, bodyid::Integer, method::AbstractString)
    method == "ELLIPSOID" && return _geo_ell(0x00)
    method == "NEAR POINT/ELLIPSOID" && return _geo_ell(_GEO_NEARPOINT)
    method == "INTERCEPT/ELLIPSOID" && return _geo_ell(_GEO_INTERCEPT)
    method == "ELLIPSOID/NEAR POINT" && return _geo_ell(_GEO_NEARPOINT)
    method == "ELLIPSOID/INTERCEPT" && return _geo_ell(_GEO_INTERCEPT)
    method == "NEAR POINT: ELLIPSOID" && return _geo_ell(_GEO_NEARPOINT)
    method == "INTERCEPT: ELLIPSOID" && return _geo_ell(_GEO_INTERCEPT)
    method == "Near point: ellipsoid" && return _geo_ell(_GEO_NEARPOINT)
    method == "Intercept: ellipsoid" && return _geo_ell(_GEO_INTERCEPT)
    method == "NADIR/ELLIPSOID" && return _geo_ell(_GEO_NADIR)
    method == "TANGENT/ELLIPSOID" && return _geo_ell(0x00, _GEO_TANGENT)
    method == "GUIDED/ELLIPSOID" && return _geo_ell(0x00, _GEO_GUIDED)
    method == "ELLIPSOID/TANGENT" && return _geo_ell(0x00, _GEO_TANGENT)
    method == "ELLIPSOID/GUIDED" && return _geo_ell(0x00, _GEO_GUIDED)
    method == "UMBRAL/TANGENT/ELLIPSOID" && return _geo_ell(0x00, _GEO_TANGENT, _GEO_UMBRAL)
    method == "PENUMBRAL/TANGENT/ELLIPSOID" && return _geo_ell(0x00, _GEO_TANGENT, _GEO_PENUMBRAL)
    method == "UMBRAL/GUIDED/ELLIPSOID" && return _geo_ell(0x00, _GEO_GUIDED, _GEO_UMBRAL)
    method == "PENUMBRAL/GUIDED/ELLIPSOID" && return _geo_ell(0x00, _GEO_GUIDED, _GEO_PENUMBRAL)
    method == "DSK/UNPRIORITIZED" && return _geo_dsk(0x00)
    method == "NADIR/DSK/UNPRIORITIZED" && return _geo_dsk(_GEO_NADIR)
    method == "INTERCEPT/DSK/UNPRIORITIZED" && return _geo_dsk(_GEO_INTERCEPT)
    method == "TANGENT/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_TANGENT)
    method == "GUIDED/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_GUIDED)
    method == "UMBRAL/TANGENT/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_TANGENT, _GEO_UMBRAL)
    method == "PENUMBRAL/TANGENT/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_TANGENT, _GEO_PENUMBRAL)
    method == "UMBRAL/GUIDED/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_GUIDED, _GEO_UMBRAL)
    method == "PENUMBRAL/GUIDED/DSK/UNPRIORITIZED" && return _geo_dsk(0x00, _GEO_GUIDED, _GEO_PENUMBRAL)
    return _geo_method_slow(ks, bodyid, method)
end

"The general case of `_geo_method`."
@noinline function _geo_method_slow(ks::KernelSet, bodyid::Integer, method::AbstractString)
    # The C wrappers of CSPICE reject an empty string before the parser sees it.
    isempty(method) && throw(SpiceError(:emptystring, "the method string is empty"))
    shape, subtyp, _, srflst = _dsk_parse_method(ks, bodyid, method)
    sub = subtyp == "" ? 0x00 : subtyp == "NEAR POINT" ? _GEO_NEARPOINT :
          subtyp == "NADIR" ? _GEO_NADIR : _GEO_INTERCEPT
    pntdef = 0x00
    trmtyp = 0x00
    # The parser has accepted the string, so each of these keywords occurs one time at
    # most, and `zzprsmet` takes a keyword only where it follows a '/'. The legacy
    # ellipsoid forms return before the keywords are read.
    legacy = shape === :ellipsoid && (sub == _GEO_NEARPOINT || (sub == 0x00 &&
             _dsk_eqstr(method, "ELLIPSOID")) || (sub == _GEO_INTERCEPT && (
             _dsk_eqstr(method, "INTERCEPT:ELLIPSOID") || _dsk_eqstr(method, "ELLIPSOID:INTERCEPT") ||
             _dsk_eqstr(method, "INTERCEPT/ELLIPSOID") || _dsk_eqstr(method, "ELLIPSOID/INTERCEPT"))))
    if !legacy
        tok = pushfirst!(_dsk_lex_method(normalize_name(method)), "/")
        for i in 2:length(tok)           # invariant: the keywords of tok[2:i-1] are recorded
            tok[i-1] == "/" || continue
            t = tok[i]
            t == "TANGENT" && (pntdef = _GEO_TANGENT)
            t == "GUIDED" && (pntdef = _GEO_GUIDED)
            t == "UMBRAL" && (trmtyp = _GEO_UMBRAL)
            t == "PENUMBRAL" && (trmtyp = _GEO_PENUMBRAL)
        end
    end
    return _GeoMethod(shape === :dsk, sub, pntdef, trmtyp, isempty(srflst) ? _GEO_NOSURF : srflst)
end

# --- aberration corrections, bodies, frames, radii --------------------------------------

@noinline _geo_bad_abcorr(s) = throw(SpiceError(:invalidoption,
    "the aberration correction \"$s\" is not recognized or not supported; the valid " *
    "values are NONE, LT, LT+S, CN, CN+S, XLT, XLT+S, XCN, XCN+S"))

@noinline function _geo_abcorr_slow(s::AbstractString)
    try
        return _parse_abcorr(s)
    catch e
        e isa SpiceError || rethrow()
        _geo_bad_abcorr(s)
    end
end

"""
    _geo_abcorr(abcorr) -> Aberration

`zzvalcor`: the correction as a plain value. A string that is not one of the nine valid
values (case and blanks ignored) throws `SpiceError(:invalidoption)`, the short message
of CSPICE. No allocation for the common spellings.
"""
@inline function _geo_abcorr(s::AbstractString)
    s == "NONE" && return _aberration(0x00)
    s == "LT" && return _aberration(_AB_LT)
    s == "LT+S" && return _aberration(_AB_LT | _AB_S)
    s == "CN" && return _aberration(_AB_LT | _AB_CN)
    s == "CN+S" && return _aberration(_AB_LT | _AB_CN | _AB_S)
    s == "XLT" && return _aberration(_AB_X | _AB_LT)
    s == "XLT+S" && return _aberration(_AB_X | _AB_LT | _AB_S)
    s == "XCN" && return _aberration(_AB_X | _AB_LT | _AB_CN)
    s == "XCN+S" && return _aberration(_AB_X | _AB_LT | _AB_CN | _AB_S)
    return _geo_abcorr_slow(s)
end
@inline _geo_abcorr(c::Aberration) = c

"`s` of the geometry calls: -1 for reception, +1 for transmission, 0 without light time."
@inline _geo_sign(c::Aberration) = c.uselt ? (c.transmit ? 1.0 : -1.0) : 0.0

@noinline _geo_no_body(what, name) = throw(SpiceError(:idcodenotfound,
    "the $what, '$name', is not a recognized name for an ephemeris object"))

"ID of the body `name` (`zzbods2c`); `SpiceError(:idcodenotfound)` for an unknown name."
@inline _geo_body(::KernelSet, id::Integer, ::AbstractString) = Int32(id)
function _geo_body(ks::KernelSet, name::AbstractString, what::AbstractString)
    code = bods2c(ks, name)
    code === nothing && _geo_no_body(what, name)
    return Int32(code)
end

"ID of the frame `ref`, or 0 if the name is unknown (`namfrm`)."
@inline _geo_namfrm(::KernelSet, ref::Integer) = Int32(ref)
_geo_namfrm(ks::KernelSet, ref::AbstractString) = Int32(namfrm(ks, ref))

"""
    _geo_frinfo(ks, id) -> (center, class, found)

`frinfo` without the class ID. Does not allocate.
"""
@inline function _geo_frinfo(ks::KernelSet, id::Int32)
    id == Int32(0) && return Int32(0), Int32(0), false
    f = get(BUILTIN_FRAMES, id, nothing)
    f === nothing || return f.center, f.class, true
    f = get(ks.frames, id, nothing)
    f === nothing || return f.center, f.class, true
    13000 <= id < 14000 && return Int32(399), Int32(2), true
    return Int32(0), Int32(0), false
end

@noinline _geo_no_frame(ref) = throw(SpiceError(:noframe, "the reference frame $ref is " *
    "not recognized by the frame subsystem; possibly a frame kernel has not been loaded"))
@noinline _geo_bad_center(ref, target, center) = throw(SpiceError(:invalidframe,
    "the reference frame $ref is not centered at the target body $target; the ID code " *
    "of the frame center is $center"))

"""
    _geo_fixfid(ks, fixref, trgcde, target) -> Int32

The ID of the body-fixed frame `fixref` after the checks of the geometry calls: the
frame is known (`SpiceError(:noframe)`) and its center is the target
(`SpiceError(:invalidframe)`).
"""
@inline function _geo_fixfid(ks::KernelSet, fixref::NameOrID, trgcde::Int32, target)
    fixfid = _geo_namfrm(ks, fixref)
    center, _, found = _geo_frinfo(ks, fixfid)
    found || _geo_no_frame(fixref)
    center == trgcde || _geo_bad_center(fixref, target, center)
    return fixfid
end

"The names of the pool variables `BODY<id>_RADII` for the IDs 1 to 999."
const _GEO_RADII_KEYS = String[string("BODY", i, "_RADII") for i in 1:999]

@noinline function _geo_bad_radii(body, v)
    v === nothing && throw(SpiceError(:kernelvarnotfound,
        "the variable BODY$(body)_RADII could not be found in the kernel pool"))
    v isa Vector{Float64} || throw(SpiceError(:typemismatch,
        "the variable BODY$(body)_RADII is not numeric"))
    length(v) == 3 || throw(SpiceError(length(v) > 3 ? :arraytoosmall : :invalidcount,
        "the variable BODY$(body)_RADII has $(length(v)) values; 3 are required"))
    throw(SpiceError(:badaxislength, "a radius of the body $body is not positive: $v"))
end

"""
    _geo_radii(ks, body) -> NTuple{3,Float64}

`zzgftreb`: the three radii of `body` from the kernel pool. Throws `SpiceError`
(`:kernelvarnotfound`, `:typemismatch`, `:invalidcount`, `:badaxislength`).
Postcondition: the three values are positive. No allocation for an ID in 1:999.
"""
@inline function _geo_radii(ks::KernelSet, body::Integer)
    v = get(ks.pool, 1 <= body <= 999 ? (@inbounds _GEO_RADII_KEYS[body]) :
            string("BODY", body, "_RADII"), nothing)
    v isa Vector{Float64} || _geo_bad_radii(body, v)
    length(v) == 3 || _geo_bad_radii(body, v)
    a, b, c = @inbounds v[1], v[2], v[3]
    (a > 0.0 && b > 0.0 && c > 0.0) || _geo_bad_radii(body, v)
    return a, b, c
end

"The vector `a` of length 3 as a `Vec3` (`SpiceError(:badarg)` otherwise)."
@inline _geo_v3(a::Vec3) = a
@inline function _geo_v3(a::AbstractVector)
    length(a) == 3 || throw(SpiceError(:badarg, "a vector must have 3 elements"))
    @inbounds return Vec3(float(a[begin]), float(a[begin+1]), float(a[begin+2]))
end

# --- positions --------------------------------------------------------------------------

"`spkezp_`: position of `targ` relative to `obs` in frame `ref` and the light time."
@inline function _geo_spkezp(ks::KernelSet, targ::Int32, et, ref::Int32, corr::Aberration,
        obs::Int32)
    if corr.geometric || targ == obs
        (r,) = spk_jet(ks, targ, et, ref, obs, Val(0))
        return r, _norm3(r) / CLIGHT
    end
    return abcorr_position(ks, targ, et, ref, corr, obs)
end

"`spkez_`: state of `targ` relative to `obs` in frame `ref` and the light time."
@inline function _geo_spkez(ks::KernelSet, targ::Int32, et, ref::Int32, corr::Aberration,
        obs::Int32)
    if corr.geometric || targ == obs
        r, v = spk_jet(ks, targ, et, ref, obs, Val(1))
        return vcat(r, v), _norm3(r) / CLIGHT
    end
    return abcorr_state(ks, targ, et, ref, corr, obs)
end

# --- the target shape -------------------------------------------------------------------

"""
    _GeoShape

The target shape of a geometry call, the saved values of `zzsinutl`: `dsk` (DSK data,
else the ellipsoid with the radii `a`, `b`, `c`), the target `trgcde`, its body-fixed
frame `fixfid`, the radii of the inner and the outer bounding sphere `minrad` and
`maxrad`, the surface list and the DSK tolerances.
"""
struct _GeoShape
    dsk::Bool
    trgcde::Int32
    fixfid::Int32
    a::Float64
    b::Float64
    c::Float64
    minrad::Float64
    maxrad::Float64
    srflst::Vector{Int}
    tols::DSKTolerances
end

"`zzsuelin`: the ellipsoid of the body `trgcde`."
@inline function _geo_ellshape(ks::KernelSet, trgcde::Int32, fixfid::Int32)
    a, b, c = _geo_radii(ks, trgcde)
    return _GeoShape(false, trgcde, fixfid, a, b, c, min(min(a, b), c), max(max(a, b), c),
        _GEO_NOSURF, _DSK_DEFAULT_TOLS)
end

"`zzdsksph` for the body `trgcde`; `SpiceError(:dskdatanotfound)` without DSK data."
function _geo_dsksph(ks::KernelSet, trgcde::Int32, srflst)
    segs = get(ks.dsk, trgcde, nothing)
    segs === nothing && _dsk_no_segments(trgcde, srflst)
    return _dsk_sph(ks, segs, trgcde, srflst)
end

"`zzsudski`: the DSK surface of the body `trgcde` with the surfaces `srflst`."
function _geo_dskshape(ks::KernelSet, tols::DSKTolerances, trgcde::Int32, fixfid::Int32, srflst)
    minrad, maxrad = _geo_dsksph(ks, trgcde, srflst)
    return _GeoShape(true, trgcde, fixfid, 0.0, 0.0, 0.0, minrad, maxrad, srflst, tols)
end

"The shape of the method `m` (`zzsuelin` or `zzsudski`)."
@inline _geo_shape(ks::KernelSet, tols::DSKTolerances, m::_GeoMethod, trgcde::Int32, fixfid::Int32) =
    m.dsk ? _geo_dskshape(ks, tols, trgcde, fixfid, m.srflst) : _geo_ellshape(ks, trgcde, fixfid)

@noinline _geo_dsk_epoch() = throw(SpiceError(:notsupported,
    "a DSK shape needs a plain (Float64) epoch: the plate model is not differentiated"))

"A plain vector for the DSK calls; `SpiceError(:notsupported)` for another number type."
@inline _geo_plain(v::Vec3{Float64}) = v
@inline _geo_plain(x::Float64) = x
_geo_plain(::Any) = _geo_dsk_epoch()

"""
    _geo_rayx(ks, sh, vertex, raydir, et) -> (found, spoint)

`zzraysfx`: the first point where the ray from `vertex` along `raydir` (body-fixed
frame of the shape, epoch `et` for the DSK segment selection) meets the surface.
Generic over the number type for an ellipsoid; no allocation.
"""
@inline function _geo_rayx(ks::KernelSet, sh::_GeoShape, vertex::SVector{3,T},
        raydir::SVector{3,T}, et) where {T}
    if sh.dsk
        xpt, _, _, found = _dsk_bux(ks, sh.tols, sh.trgcde, sh.srflst, _geo_plain(et),
            sh.fixfid, _geo_plain(vertex), _geo_plain(raydir))
        return found, SVector{3,T}(xpt)
    end
    return _m_surfpt(vertex, raydir, sh.a, sh.b, sh.c)
end

"`zzraynp`: the point of the ray's line nearest to the ellipsoid (or to the outer sphere of DSK data) and the distance."
@inline function _geo_raynp(sh::_GeoShape, vertex::SVector{3}, raydir::SVector{3})
    sh.dsk && return npedln(sh.maxrad, sh.maxrad, sh.maxrad, vertex, raydir)
    return npedln(sh.a, sh.b, sh.c, vertex, raydir)
end

"`zzsbfnrm` with `vhat`: the outward unit normal of the shape at the surface point `spoint`."
@inline function _geo_normal(ks::KernelSet, sh::_GeoShape, spoint::SVector{3,T}, et) where {T}
    sh.dsk || return _m_surfnm(sh.a, sh.b, sh.c, spoint)
    n = _dsk_bun(ks, sh.tols, sh.trgcde, sh.srflst, _geo_plain(et), sh.fixfid,
        _geo_plain(spoint))
    return SVector{3,T}(_dsk_vhat(n))
end
