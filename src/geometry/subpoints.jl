# Sub-observer point and sub-solar point: `subpnt` and `subslr` of CSPICE N0067.
#
# The light time is the one between the observer and the sub-point, not the center of
# the target; the iterations, their counts and their convergence tests are those of the
# C code. The functions are generic over the epoch number type for an ellipsoid; all
# discrete decisions use `primal`.

@noinline _geo_same_bodies(what, a) = throw(SpiceError(:bodiesnotdistinct,
    "in computing the $what, the observing body and the target body are the same: $a"))
@noinline _geo_bad_subtype(method) = throw(SpiceError(:invalidsubtype,
    "the sub-point type is invalid or was not found in the method string <$method>"))
@noinline _geo_no_subpoint(what) = throw(SpiceError(:subpointnotfound,
    "no $what point was found on the surface that the DSK data define"))
@noinline _geo_no_intercept() = throw(SpiceError(:degeneratecase,
    "no intercept of the ray toward the target center was found"))

@inline function _geo_near(m::_GeoMethod, method)
    m.subtyp == 0x00 && _geo_bad_subtype(method)
    near = m.dsk ? m.subtyp == _GEO_NADIR : (m.subtyp == _GEO_NEARPOINT || m.subtyp == _GEO_NADIR)
    near || m.subtyp == _GEO_INTERCEPT || _geo_bad_subtype(method)
    return near
end

"""
    _geo_subpoint(ks, sh, radii, near, pos, et, what) -> (spoint, alt)

One sub-point computation of `subpnt` and `subslr`: the point of the surface below
`pos` (body-fixed) and the distance from `pos` to it. `near` selects the nearest point
of the ellipsoid (for DSK data, the intercept of the ray toward that point); otherwise
the intercept of the ray toward the center.
"""
@inline function _geo_subpoint(ks::KernelSet, sh::_GeoShape, radii::NTuple{3,Float64},
        near::Bool, pos::SVector{3}, et, what)
    if near
        spoint, alt = _m_nearpt(pos, radii[1], radii[2], radii[3])
        sh.dsk || return spoint, alt
        found, spoint = _geo_rayx(ks, sh, pos, spoint - pos, et)
        found || _geo_no_subpoint(what)
        return spoint, _norm3(spoint - pos)
    end
    if sh.dsk
        found, spoint = _geo_rayx(ks, sh, pos, -pos, et)
        found || _geo_no_subpoint(what)
        return spoint, _norm3(spoint - pos)
    end
    found, spoint = _m_surfpt(pos, -pos, radii[1], radii[2], radii[3])
    found || _geo_no_intercept()
    return spoint, _norm3(pos - spoint)
end

"A shape that carries the DSK data of the method, or a placeholder for an ellipsoid."
@inline function _geo_subshape(ks::KernelSet, tols::DSKTolerances, m::_GeoMethod,
        trgcde::Int32, fixfid::Int32)
    m.dsk && return _geo_dskshape(ks, tols, trgcde, fixfid, m.srflst)
    return _GeoShape(false, trgcde, fixfid, 1.0, 1.0, 1.0, 1.0, 1.0, _GEO_NOSURF, tols)
end

function _geo_subpnt(ks::KernelSet, tols::DSKTolerances, method::AbstractString,
        target::NameOrID, et0::Real, fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID)
    corr = _geo_abcorr(abcorr)
    trgcde = _geo_body(ks, target, "target")
    obscde = _geo_body(ks, obsrvr, "observer")
    obscde == trgcde && _geo_same_bodies("sub-observer point", obsrvr)
    fixfid = _geo_fixfid(ks, fixref, trgcde, target)
    m = _geo_method(ks, trgcde, method)
    near = _geo_near(m, method)
    sh = _geo_subshape(ks, tols, m, trgcde, fixfid)
    et = check_epoch(et0)
    s = _geo_sign(corr)
    tpos, lt = _geo_spkezp(ks, trgcde, et, fixfid, corr, obscde)
    obspos = -tpos
    trgepc = et + s * lt
    radii = _geo_radii(ks, trgcde)
    primal(_norm3(obspos)) == 0 && throw(SpiceError(:noseparation,
        "the distance between the observer $obsrvr and the target $target is zero"))
    spoint, alt = _geo_subpoint(ks, sh, radii, near, obspos, trgepc, "sub-observer")
    corr.uselt || return spoint, et, spoint - obspos
    lt = alt / CLIGHT
    trgepc = et + s * lt
    nitr = corr.converged ? 5 : 1
    robs, vobs = _ssb(ks, obscde, et, Int32(1))
    i = 0
    ltdiff = 1.0
    etdiff = 1.0
    prevlt = primal(lt)
    prevet = primal(trgepc)
    # invariant: `spoint` and `lt` are the sub-point and its light time after i iterations
    while i < nitr && ltdiff > abs(primal(lt)) * 1e-17 && etdiff > 0.0
        j2pos = robs - _ssbpos(ks, trgcde, trgepc, Int32(1))
        xform = frame_rotation(ks, Int32(1), fixfid, trgepc)
        obspos = _m_mxv(xform, j2pos)
        if corr.stellar
            # The stellar aberration of the vector from the observer to the sub-point
            # is applied to the observer, with the opposite sign.
            subvec = spoint - obspos
            subvj2 = _m_mtxv(xform, subvec)
            corvj2 = _stelab(subvj2, corr.transmit ? -vobs : vobs)
            obspos = obspos - (_m_mxv(xform, corvj2) - subvec)
        end
        spoint, alt = _geo_subpoint(ks, sh, radii, near, obspos, trgepc, "sub-observer")
        lt = alt / CLIGHT
        trgepc = et + s * lt
        ltdiff = abs(primal(lt) - prevlt)
        etdiff = abs(primal(trgepc) - prevet)
        prevlt = primal(lt)
        prevet = primal(trgepc)
        i += 1
    end
    return spoint, trgepc, spoint - obspos
end

function _geo_subslr(ks::KernelSet, tols::DSKTolerances, method::AbstractString,
        target::NameOrID, et0::Real, fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID)
    corr = _geo_abcorr(abcorr)
    corr.transmit && throw(SpiceError(:notsupported,
        "the aberration correction $(string(corr)) calls for transmission; subslr " *
        "accepts reception corrections only"))
    trgcde = _geo_body(ks, target, "target")
    obscde = _geo_body(ks, obsrvr, "observer")
    trgcde == Int32(10) && throw(SpiceError(:invalidtarget,
        "the target body is the Sun; the sub-solar point is not defined"))
    fixfid = _geo_fixfid(ks, fixref, trgcde, target)
    m = _geo_method(ks, trgcde, method)
    near = _geo_near(m, method)
    sh = _geo_subshape(ks, tols, m, trgcde, fixfid)
    et = check_epoch(et0)
    s = _geo_sign(corr)
    tpos, lt = _geo_spkezp(ks, trgcde, et, fixfid, corr, obscde)
    obspos = -tpos
    trgepc = et + s * lt
    radii = _geo_radii(ks, trgcde)
    spos, _ = _geo_spkezp(ks, Int32(10), trgepc, fixfid, corr, trgcde)
    spoint, _ = _geo_subpoint(ks, sh, radii, near, spos, trgepc, "sub-solar")
    corr.uselt || return spoint, et, spoint - obspos
    lt = _norm3(obspos - spoint) / CLIGHT
    trgepc = et + s * lt
    nitr = corr.converged ? 10 : 1
    robs = _ssbpos(ks, obscde, et, Int32(1))
    i = 0
    ltdiff = 1.0
    etdiff = 1.0
    prevlt = primal(lt)
    prevet = primal(trgepc)
    # invariant: `spoint` and `lt` are the sub-solar point and its light time after i
    # iterations, and `obspos` is the observer at the epoch `trgepc`
    while i < nitr && ltdiff > abs(primal(lt)) * 1e-17 && etdiff > 0.0
        # The Sun as seen from the current sub-solar point at the current target epoch.
        f = _GeoFixed(_geo_frozen(spoint), zero(Vec3{Float64}), 0.0, trgcde, fixfid)
        sunst, _ = _geo_cvo_point(ks, Int32(10), trgepc, fixfid, corr, f, spoint)
        spos = SVector(sunst[1], sunst[2], sunst[3]) + spoint
        spoint, _ = _geo_subpoint(ks, sh, radii, near, spos, trgepc, "sub-solar")
        xform = frame_rotation(ks, Int32(1), fixfid, trgepc)
        obspos = _m_mxv(xform, robs - _ssbpos(ks, trgcde, trgepc, Int32(1)))
        lt = _norm3(obspos - spoint) / CLIGHT
        trgepc = et + s * lt
        ltdiff = abs(primal(lt) - prevlt)
        etdiff = abs(primal(trgepc) - prevet)
        prevlt = primal(lt)
        prevet = primal(trgepc)
        i += 1
        xform = frame_rotation(ks, Int32(1), fixfid, trgepc)
        obspos = _m_mxv(xform, robs - _ssbpos(ks, trgcde, trgepc, Int32(1)))
    end
    if corr.stellar
        f = _GeoFixed(_geo_frozen(spoint), zero(Vec3{Float64}), 0.0, trgcde, fixfid)
        st, _ = _geo_cvt_point(ks, f, et, fixfid, corr, obscde, spoint)
        return spoint, trgepc, SVector(st[1], st[2], st[3])
    end
    return spoint, trgepc, spoint - obspos
end

@inline _geo_frozen(v::SVector{3}) = Vec3{Float64}(primal(v[1]), primal(v[2]), primal(v[3]))

"""
    _geo_cvo_point(ks, targ, et, fixfid, corr, f, point) -> (state, lt)

`spkcpo` for the observer at the surface point `point` of the body `f.ctr`, in the
body-fixed frame `fixfid` with the locus `OBSERVER`. For a plain `point` this is
`_geo_cvo`; for another number type the derivative parts of the point enter the
position at first order (the observer moves with the point).
"""
@inline _geo_cvo_point(ks::KernelSet, targ::Int32, et, fixfid::Int32, corr::Aberration,
    f::_GeoFixed, ::Vec3{Float64}) = _geo_cvo(ks, targ, et, fixfid, _GEO_AT_OBSERVER, corr, f)
@inline function _geo_cvo_point(ks::KernelSet, targ::Int32, et, fixfid::Int32,
        corr::Aberration, f::_GeoFixed, point::SVector{3})
    st, lt = _geo_cvo(ks, targ, et, fixfid, _GEO_AT_OBSERVER, corr, f)
    d = point - f.pos
    return vcat(SVector(st[1], st[2], st[3]) - d, SVector(st[4], st[5], st[6]) - zero(d)), lt
end

"`spkcpt` for the target at the surface point `point` of the body `f.ctr`, with the locus `TARGET`; see `_geo_cvo_point`."
@inline _geo_cvt_point(ks::KernelSet, f::_GeoFixed, et, fixfid::Int32, corr::Aberration,
    obs::Int32, ::Vec3{Float64}) = _geo_cvt(ks, f, et, fixfid, _GEO_AT_TARGET, corr, obs)
@inline function _geo_cvt_point(ks::KernelSet, f::_GeoFixed, et, fixfid::Int32,
        corr::Aberration, obs::Int32, point::SVector{3})
    st, lt = _geo_cvt(ks, f, et, fixfid, _GEO_AT_TARGET, corr, obs)
    d = point - f.pos
    return vcat(SVector(st[1], st[2], st[3]) + d, SVector(st[4], st[5], st[6]) - zero(d)), lt
end

# --- public interface -------------------------------------------------------------------

"""
    subpnt([ctx,] method, target, et, fixref, abcorr, obsrvr) -> (spoint, trgepc, srfvec)

Return the sub-observer point on `target` as seen by `obsrvr` at `et` (seconds past
J2000 TDB), as `subpnt_c` of CSPICE.

- `spoint`: the point, in km, in the body-fixed frame `fixref` of the target.
- `trgepc`: the epoch at the point (seconds past J2000 TDB). With a light time
  correction it is `et` minus (reception) or plus (transmission) the light time between
  the observer and the point. Without one it is `et`.
- `srfvec`: the vector from the observer to the point, in km, in `fixref` at `trgepc`.

`method` selects the definition of the point and the shape of the target:

| `method` | Point |
|---|---|
| `"NEAR POINT/ELLIPSOID"` | point of the reference ellipsoid nearest to the observer |
| `"INTERCEPT/ELLIPSOID"` | where the line from the observer to the target center meets the ellipsoid |
| `"NADIR/DSK/UNPRIORITIZED"` | where the ray from the observer toward the nearest point of the ellipsoid meets the DSK surface |
| `"INTERCEPT/DSK/UNPRIORITIZED"` | where the ray from the observer toward the target center meets the DSK surface |

Case is not significant. `UNPRIORITIZED` means that all loaded DSK segments of the
target are used, whatever their load order. A DSK method can end with
`"/SURFACES = <list>"`, a list of surface names or IDs. The CSPICE forms
`"Near point: ellipsoid"` and `"Intercept: ellipsoid"` are also valid.

`abcorr` is a correction string of [`spkezr`](@ref) or an [`Aberration`](@ref).
Bodies and frames are names or integer IDs. With an ellipsoid method, `et` can be a
`ForwardDiff.Dual`, and the call does not allocate when the bodies and the frame are
integer IDs. A DSK method uses the tolerances that [`dskstl`](@ref) set on the context.

Throws `SpiceError` with the kind
- `:invalidoption` if `abcorr` is not a valid correction;
- `:idcodenotfound` if a body name is not known;
- `:bodiesnotdistinct` if `target` and `obsrvr` are the same body;
- `:noframe` if `fixref` is not known, `:invalidframe` if it is not centered on `target`;
- `:badmethodsyntax` if `method` is malformed, `:badpriorityspec` if a DSK method does
  not have `UNPRIORITIZED`, `:invalidsubtype` if it has no valid point definition;
- `:noseparation` if the observer is at the target center;
- `:subpointnotfound` if the ray misses the DSK surface, `:degeneratecase` if it
  misses the ellipsoid;
- the kinds of [`srfrec`](@ref) if the radii of `target` are missing or not valid,
  `:noloadeddskfiles` if no DSK file is loaded and `:dskdatanotfound` if no DSK data are
  loaded for it (DSK method).

Also throws the errors of [`spkpos`](@ref).

# Examples
```julia
# Needs an SPK, a PCK with the radii and the rotation of Mars, and a leap seconds kernel.
et = str2et("2025 MAR 01 12:00:00 UTC")
spoint, trgepc, srfvec = subpnt("NEAR POINT/ELLIPSOID", "MARS", et, "IAU_MARS", "LT+S", "EARTH")
radius, lon, lat = reclat(spoint)      # km, rad, rad (planetocentric)
distance = vnorm(srfvec)               # km from the Earth to the point
```

See also [`subslr`](@ref), [`sincpt`](@ref), [`gfsubc`](@ref).
"""
function subpnt end

"""
    subslr([ctx,] method, target, et, fixref, abcorr, obsrvr) -> (spoint, trgepc, srfvec)

Return the sub-solar point on `target` as seen by `obsrvr` at `et` (seconds past J2000
TDB), as `subslr_c` of CSPICE.

- `spoint`: the point, in km, in the body-fixed frame `fixref` of the target.
- `trgepc`: the epoch at the point (seconds past J2000 TDB). With a light time
  correction it is `et` minus the light time between the observer and the point.
  Without one it is `et`.
- `srfvec`: the vector from the observer to the point, in km, in `fixref` at `trgepc`.

The methods are those of [`subpnt`](@ref), with the Sun in place of the observer in the
definition of the point. With a correction, the Sun is taken where it appears from the
point at `trgepc`. The other arguments and the notes on `et`, allocation and `ctx` are
those of [`subpnt`](@ref).

Throws `SpiceError(:notsupported)` if `abcorr` is a transmission correction, and
`SpiceError(:invalidtarget)` if `target` is the Sun. The other kinds are those of
[`subpnt`](@ref), except `:bodiesnotdistinct` and `:noseparation`: the observer and the
target can be the same body.

See also [`subpnt`](@ref), [`ilumin`](@ref).
"""
function subslr end

for (f, g) in ((:subpnt, :_geo_subpnt), (:subslr, :_geo_subslr))
    @eval begin
        $f(ks::KernelSet, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
            abcorr::AbcorrArg, obsrvr::NameOrID) =
            $g(ks, _DSK_DEFAULT_TOLS, method, target, et, fixref, abcorr, obsrvr)
        $f(ctx::SpiceContext, method::AbstractString, target::NameOrID, et::Real,
            fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID) =
            $g(kernelset(ctx), dsktolerances(ctx), method, target, et, fixref, abcorr, obsrvr)
        $f(method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
            abcorr::AbcorrArg, obsrvr::NameOrID) =
            $f(default_context(), method, target, et, fixref, abcorr, obsrvr)
    end
end
