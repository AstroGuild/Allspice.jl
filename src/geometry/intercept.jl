# Surface intercept of a ray: `sincpt` of CSPICE N0067 with `zzsfxcor`, and `srfrec`.
#
# The light time is the one between the observer and the intercept point. The steps,
# the iteration counts and the convergence tests are those of `zzsfxcor`. The function
# is generic over the epoch number type for an ellipsoid; all discrete decisions use
# `primal`.

@noinline _geo_sub_in_method(method) = throw(SpiceError(:invalidmethod,
    "the method string <$method> contains a sub-observer point type, which is not valid here"))

"""
    _geo_stl_dir(j2dir, vobs, transmit) -> SVector{3}

The direction that the stellar aberration correction maps to `j2dir`, for an observer
with the velocity `vobs`: the inverse of `stelab` (reception) or of `stlabx`
(transmission), by the fixed-point iteration of `zzsfxcor` (at most five steps, to a
relative error of `1e-17`).
"""
@inline function _geo_stl_dir(j2dir::SVector{3}, vobs::SVector{3}, transmit::Bool)
    vfwd = transmit ? -vobs : vobs      # the velocity of the correction to invert
    stldir = _stelab(j2dir, -vfwd)
    relerr = 1.0
    i = 1
    # invariant: `stldir` is the estimate after i - 1 steps
    while i <= 5 && relerr > 1e-17
        stlerr = j2dir - _stelab(stldir, vfwd)
        stldir = stlerr + stldir
        relerr = primal(_norm3(stlerr)) / primal(_norm3(stldir))
        i += 1
    end
    return stldir
end

"""
    _geo_sfxcor(ks, sh, et, corr, obscde, dfrcde, dvec) -> (found, spoint, trgepc, srfvec)

`zzsfxcor`: the intercept on the shape `sh` of the ray from the body `obscde` along
`dvec` (frame `dfrcde`), with the aberration corrections `corr`.

Contract: `dvec` is not zero. The frame of the ray is evaluated at `et`, or, for a
frame that is not inertial and not centered on the observer, at the epoch corrected for
the light time to its center. `LT` forms make one light-time iteration and `CN` forms
at most ten, with the convergence test of CSPICE. Postcondition: if `found`, `spoint`
is on the surface and `srfvec` is the vector from the observer to it.
"""
function _geo_sfxcor(ks::KernelSet, sh::_GeoShape, et, corr::Aberration, obscde::Int32,
        dfrcde::Int32, dvec::Vec3{Float64})
    trgcde = sh.trgcde
    fixfid = sh.fixfid
    s = _geo_sign(corr)
    loccor = corr.uselt ? Aberration(; lt=true, converged=corr.converged, transmit=corr.transmit) :
             Aberration()
    tpos, lt = _geo_spkezp(ks, trgcde, et, fixfid, loccor, obscde)
    obspos = -tpos
    trgepc = et + s * lt
    dcentr, dclass, _ = _geo_frinfo(ks, dfrcde)
    refepc = et
    if dclass != Int32(1) && corr.uselt && dcentr != obscde
        _, ltcent = _geo_spkezp(ks, dcentr, et, Int32(1), corr, obscde)
        refepc = et + s * ltcent
    end
    j2dir = frame_rotation(ks, dfrcde, Int32(1), refepc) * dvec
    j2geom = j2dir
    j2tmat = frame_rotation(ks, Int32(1), fixfid, trgepc)
    robs, vobs = _ssb(ks, obscde, et, Int32(1))
    corr.stellar && (j2dir = _geo_stl_dir(j2dir, vobs, corr.transmit))
    trgdir = j2tmat * j2dir
    none = (false, zero(obspos), trgepc, zero(obspos))
    range = primal(_norm3(obspos))
    range == 0.0 && throw(SpiceError(:noseparation,
        "the distance between the observer and the target is zero"))
    if range > sh.maxrad * 1.01
        # The ray cannot reach the target if it misses the outer bounding sphere.
        reject = asin(sh.maxrad * 1.01 / range)
        primal(vsep(-obspos, trgdir)) > reject && return none
    end
    found, spoint = _geo_rayx(ks, sh, obspos, trgdir, trgepc)
    (found && !corr.uselt) && return true, spoint, trgepc, spoint - obspos
    xform = j2tmat
    if !found
        # The ray misses at the first epoch: move the epoch to that of the point of the
        # ray nearest to the target and try again.
        nitr = corr.converged ? 3 : 1
        i = 1
        # invariant: no intercept was found in the i - 1 attempts made
        while i <= nitr && !found
            pnear, _ = _geo_raynp(sh, obspos, trgdir)
            lt = _norm3(obspos - pnear) / CLIGHT
            trgepc = et + s * lt
            xform = frame_rotation(ks, Int32(1), fixfid, trgepc)
            obspos = xform * (robs - _ssbpos(ks, trgcde, trgepc, Int32(1)))
            trgdir = xform * j2dir
            found, spoint = _geo_rayx(ks, sh, obspos, trgdir, trgepc)
            i += 1
        end
        found || return none
    end
    nitr = corr.converged ? 10 : 1
    lt = _norm3(obspos - spoint) / CLIGHT
    trgepc = et + s * lt
    prevlt = 0.0
    prevet = primal(trgepc)
    i = 0
    ltdiff = 1.0
    etdiff = 1.0
    # invariant: `spoint` and `lt` are the intercept and its light time after i iterations
    while i < nitr && ltdiff > abs(primal(lt)) * 1e-17 && etdiff > 0.0
        xform = frame_rotation(ks, Int32(1), fixfid, trgepc)
        obspos = xform * (robs - _ssbpos(ks, trgcde, trgepc, Int32(1)))
        trgdir = xform * j2dir
        found, spoint = _geo_rayx(ks, sh, obspos, trgdir, trgepc)
        if !found
            i < nitr - 1 || return none
            pnear, _ = _geo_raynp(sh, obspos, trgdir)
            spoint, _ = nplnpt(obspos, trgdir, pnear)
        end
        lt = _norm3(obspos - spoint) / CLIGHT
        trgepc = et + s * lt
        ltdiff = abs(primal(lt) - prevlt)
        etdiff = abs(primal(trgepc) - prevet)
        prevlt = primal(lt)
        prevet = primal(trgepc)
        i += 1
    end
    udir = _vhat(xform * j2geom)
    return true, spoint, trgepc, _norm3(spoint - obspos) * udir
end

function _geo_sincpt(ks::KernelSet, tols::DSKTolerances, method::AbstractString,
        target::NameOrID, et::Real, fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
        dref::NameOrID, dvec::AbstractVector)
    corr = _geo_abcorr(abcorr)
    trgcde = _geo_body(ks, target, "target")
    obscde = _geo_body(ks, obsrvr, "observer")
    obscde == trgcde && _geo_same_bodies("surface intercept", obsrvr)
    fixfid = _geo_fixfid(ks, fixref, trgcde, target)
    dv = _geo_pos(dvec)
    (dv[1] == 0.0 && dv[2] == 0.0 && dv[3] == 0.0) && throw(SpiceError(:zerovector,
        "the direction vector of the ray is the zero vector"))
    dfrcde = _geo_namfrm(ks, dref)
    _geo_frinfo(ks, dfrcde)[3] || _geo_no_frame(dref)
    m = _geo_method(ks, trgcde, method)
    m.subtyp == 0x00 || _geo_sub_in_method(method)
    sh = _geo_shape(ks, tols, m, trgcde, fixfid)
    found, spoint, trgepc, srfvec = _geo_sfxcor(ks, sh, check_epoch(et), corr, obscde, dfrcde, dv)
    found || return nothing
    return spoint, trgepc, srfvec
end

"""
    sincpt(method, target, et, fixref, abcorr, obsrvr, dref, dvec)
        -> Union{Nothing,Tuple{SVector{3},Real,SVector{3}}}

Return the first point where the ray from `obsrvr` along `dvec` (given in the frame
`dref`) meets the surface of `target`: `(spoint, trgepc, srfvec)`, with `spoint` in the
body-fixed frame `fixref` (km), the epoch `trgepc` at the point and the vector `srfvec`
from the observer to the point, as `sincpt_c`; `nothing` if the ray misses the target.
Forms: on a `KernelSet`, on a `SpiceContext` (which gives the DSK tolerances), and on
the default context.

`method` is `"ELLIPSOID"` or `"DSK/UNPRIORITIZED[/SURFACES = <list>]"`. The light time
is the one to the intercept point, and the ray direction is corrected for stellar
aberration with the `+S` forms.

Contract: `dvec` has 3 elements and is not zero (`SpiceError(:zerovector)`); `fixref`
is centered on `target`. For an ellipsoid the epoch can be a dual number, and the call
does not allocate for integer IDs. Throws `SpiceError`: `:invalidoption`,
`:idcodenotfound`, `:bodiesnotdistinct`, `:noframe`, `:invalidframe`, the kinds of the
method string, `:invalidmethod` (a sub-point type in the method), `:noseparation`, and
the errors of `spkpos`, of the radii lookup and of [`dskxsi`](@ref).
"""
sincpt(ks::KernelSet, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, dref::NameOrID, dvec::AbstractVector) =
    _geo_sincpt(ks, _DSK_DEFAULT_TOLS, method, target, et, fixref, abcorr, obsrvr, dref, dvec)
sincpt(ctx::SpiceContext, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, dref::NameOrID, dvec::AbstractVector) =
    _geo_sincpt(kernelset(ctx), dsktolerances(ctx), method, target, et, fixref, abcorr, obsrvr,
        dref, dvec)
sincpt(method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, dref::NameOrID, dvec::AbstractVector) =
    sincpt(default_context(), method, target, et, fixref, abcorr, obsrvr, dref, dvec)

"""
    srfrec(body, lon, lat) -> SVector{3}

Return the point of the reference ellipsoid of `body` at the planetocentric longitude
`lon` and latitude `lat` (radians), in the body-fixed frame (km), as `srfrec_c`.
Forms: on a `KernelSet`, on a `SpiceContext`, and on the default context.

Contract: the radii of `body` are in the kernel pool and positive (`SpiceError`
`:kernelvarnotfound`, `:invalidcount`, `:badaxislength`). Postcondition: the point is
on the ellipsoid. No allocation for an integer ID in 1:999.
"""
function srfrec(ks::KernelSet, body::NameOrID, lon::Real, lat::Real)
    a, b, c = _geo_radii(ks, _geo_body(ks, body, "body"))
    return edpnt(_m_latrec(1.0, float(lon), float(lat)), a, b, c)
end
srfrec(ctx::SpiceContext, body::NameOrID, lon::Real, lat::Real) =
    srfrec(kernelset(ctx), body, lon, lat)
srfrec(body::NameOrID, lon::Real, lat::Real) = srfrec(default_context(), body, lon, lat)
