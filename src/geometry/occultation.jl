# Occultation state of two bodies: `occult` of CSPICE N0067, with `zzgfocin`,
# `zzgfocst` (`zzgfocu.c`), `zzocced` (two ellipsoids) and `zzasryel` (angular
# separation of a ray and an ellipse).
#
# The decisions are discrete, so the functions take a plain epoch. The ellipsoid and
# the point cases do not allocate.

const _GEO_OCC_POINT = 0x01
const _GEO_OCC_ELLIPSOID = 0x02
const _GEO_OCC_DSK = 0x03

"`dasine_`: the arcsine of `arg` clamped to [-1, 1]; `SpiceError(:inputoutofbounds)` if `|arg| > 1 + tol`."
@inline function _geo_dasine(arg::Float64, tol::Float64)
    abs(arg) - tol > 1.0 && throw(SpiceError(:inputoutofbounds,
        "the argument $arg of the arcsine is outside of [-1, 1] by more than $tol"))
    return asin(max(-1.0, min(1.0, arg)))
end

"""
    _geo_asryel(domin, ellips, vertex, dir) -> (angle, extpt)

`zzasryel`: the smallest (`domin`) or the largest angular separation between the ray
from `vertex` along `dir` and the ellipse `ellips`, and the point of the ellipse where
it occurs. The angle is negative if the ray goes through the plane region that the
ellipse bounds.

Contract: the vertex is not in the plane of the ellipse (`SpiceError(:degeneratecase)`).
The search is the one of CSPICE: a scan of 400 points of the ellipse (320 for a minimum
with the ray outside of the ellipse), then a golden section search of at most 101 steps
to an interval of `1e-9` radians.
"""
function _geo_asryel(domin::Bool, ellips::SpiceEllipse{Float64}, vertex::Vec3{Float64},
        dir::Vec3{Float64})
    center, smajor, sminor = ellips.center, ellips.semimajor, ellips.semiminor
    a = _m_norm(smajor)
    b = _m_norm(sminor)
    (_m_iszero(smajor) || _m_iszero(sminor)) && throw(SpiceError(:invalidaxislength,
        "a semi-axis of the ellipse has the length zero"))
    eplane = psv2pl(center, smajor, sminor)
    _m_iszero(dir) && throw(SpiceError(:zerovector, "the ray direction is the zero vector"))
    vprj = _m_vprjp(vertex, eplane)
    _m_norm(vertex - vprj) == 0.0 && throw(SpiceError(:degeneratecase,
        "the vertex of the ray is in the plane of the ellipse"))
    nxpts, xpt = inrypl(vertex, dir, eplane)
    nxpts == -1 && throw(SpiceError(:degeneratecase, "the ray lies in the plane of the ellipse"))
    npt = 400
    asign = 1.0
    if nxpts != 0
        xoff = xpt - center
        acomp = _m_dot(xoff, smajor) / a
        bcomp = _m_dot(xoff, sminor) / b
        level = acomp * acomp / (a * a) + bcomp * bcomp / (b * b)
        if level <= 1.0
            asign = -1.0
        elseif domin
            npt = 320
        end
    end
    v2 = vertex - center
    extprx = domin ? 2.0 : 0.0
    extidx = 0
    p2 = _M_TWOPI
    delta = p2 / npt
    udir = _m_hat(dir)
    diff = zero(Vec3{Float64})
    for i in 0:npt-1                     # invariant: extidx, extprx are the extreme of the points 0 .. i-1
        theta = i * delta
        diff = -1.0 * v2 + cos(theta) * smajor + sin(theta) * sminor
        proxy = _m_norm(_m_hat(diff) - udir)
        if domin ? proxy < extprx : proxy > extprx
            extidx = i
            extprx = proxy
        end
    end
    gr = (3.0 - sqrt(5.0)) / 2.0
    lower = p2 / npt * (extidx - 1) + p2
    upper = p2 / npt * (extidx + 1) + p2
    btween = p2 / npt * extidx + p2
    btwprx = extprx
    nitr = 0
    # invariant: the extreme is in [lower, upper] and `btween` is the best point found
    while nitr <= 100 && upper - lower > 1e-9
        newpt = btween - lower > upper - btween ? lower + gr * (btween - lower) :
                btween + gr * (upper - btween)
        diff = -1.0 * v2 + cos(newpt) * smajor + sin(newpt) * sminor
        proxy = _m_norm(_m_hat(diff) - udir)
        if newpt > btween
            btween, newpt = newpt, btween
            btwprx, proxy = proxy, btwprx
        end
        if domin ? proxy > btwprx : proxy < btwprx
            lower = newpt
        else
            upper = btween
            btween = newpt
            btwprx = proxy
        end
        nitr += 1
    end
    # As in CSPICE, the result is the last point that the search evaluated.
    return vsep(diff, udir) * asign, center + (diff + v2)
end

"The largest angular radius of the ellipse `(ctr, maj, min)` seen from `view` (an upper bound if the view point is close)."
@inline function _geo_occ_maxang(view::Vec3{Float64}, ctr::Vec3{Float64}, maj::Vec3{Float64},
        mn::Vec3{Float64})
    lplane = psv2pl(ctr, maj, mn)
    vpproj = _m_vprjp(view, lplane)
    vph = _m_norm(view - vpproj)
    ubdist = _m_norm(vpproj - ctr) + _m_norm(maj)
    return atan(ubdist, vph)
end

"""
    _geo_occed(viewpt, centr1, semax1, centr2, semax2) -> Int

`zzocced`: the occultation code of two ellipsoids seen from `viewpt`. The columns of
`semax1` and `semax2` are the semi-axis vectors. The code is 0 (no occultation), ±1
(partial), ±2 (annular) or ±3 (total); it is negative if the second ellipsoid is in
front of the first one.

Contract: positive axis lengths (`SpiceError(:badaxislength)`), orthogonal axes
(`:notarotation`), and a view point outside of both ellipsoids, which do not overlap
(`:notdisjoint`).
"""
function _geo_occed(viewpt::Vec3{Float64}, centr1::Vec3{Float64}, semax1::Mat3{Float64},
        centr2::Vec3{Float64}, semax2::Mat3{Float64})
    ctrs = (centr1, centr2)
    u11, r11 = _m_unorm(semax1[:, 1])
    u12, r12 = _m_unorm(semax1[:, 2])
    u13, r13 = _m_unorm(semax1[:, 3])
    u21, r21 = _m_unorm(semax2[:, 1])
    u22, r22 = _m_unorm(semax2[:, 2])
    u23, r23 = _m_unorm(semax2[:, 3])
    rmat = (hcat(u11, u12, u13), hcat(u21, u22, u23))
    r = ((r11, r12, r13), (r21, r22, r23))
    minrad = (min(min(r11, r12), r13), min(min(r21, r22), r23))
    maxrad = (max(max(r11, r12), r13), max(max(r21, r22), r23))
    (minrad[1] <= 0.0 || minrad[2] <= 0.0) && throw(SpiceError(:badaxislength,
        "the smallest semi-axis lengths of the two ellipsoids are $(minrad[1]) and $(minrad[2])"))
    tpos = (centr1 - viewpt, centr2 - viewpt)
    dist = (_m_norm(tpos[1]), _m_norm(tpos[2]))
    (dist[1] == 0.0 || dist[2] == 0.0) && throw(SpiceError(:notdisjoint,
        "the center of a target coincides with the viewing point"))
    for i in 1:2
        _m_isrot(rmat[i], 1e-14, 1e-12) || throw(SpiceError(:notarotation,
            "the semi-axes of target $i do not form a rotation matrix; the determinant " *
            "is $(_m_det(rmat[i]))"))
    end
    t12pos = tpos[2] - tpos[1]
    ttdist = _m_norm(t12pos)
    ttdist <= minrad[1] + minrad[2] && throw(SpiceError(:notdisjoint,
        "the distance $ttdist between the targets is not larger than the sum of their " *
        "smallest radii"))
    xtpos = (_m_mtxv(rmat[1], tpos[1]), _m_mtxv(rmat[2], tpos[2]))
    for i in 1:2
        d1 = xtpos[i][1] / r[i][1]
        d2 = xtpos[i][2] / r[i][2]
        d3 = xtpos[i][3] / r[i][3]
        d1 * d1 + d2 * d2 + d3 * d3 < 1.0 && throw(SpiceError(:notdisjoint,
            "the viewing point is inside target $i"))
    end
    minang1 = _geo_dasine(minrad[1] / dist[1], 1e-12)
    minang2 = _geo_dasine(minrad[2] / dist[2], 1e-12)
    maxang1 = _geo_occ_bound(r[1], xtpos[1], dist[1], maxrad[1])
    maxang2 = _geo_occ_bound(r[2], xtpos[2], dist[2], maxrad[2])
    trgsep = vsep(tpos[1], tpos[2])
    trgsep > maxang1 + maxang2 && return 0
    bigidx = maxang1 > maxang2 ? 1 : 2
    smlidx = 3 - bigidx
    if minang1 + minang2 > trgsep
        # The inscribed spheres overlap as seen from the view point.
        t2sep = vsep(-tpos[1], t12pos)
        frtidx, s = t2sep < _M_HALFPI - minang1 ? (2, -1) : (1, 1)
        minang1 >= maxang2 + trgsep && return frtidx == 1 ? 3s : 2s
        minang2 >= maxang1 + trgsep && return frtidx == 2 ? 3s : 2s
        (minang1 + trgsep > maxang2 && minang2 + trgsep > maxang1) && return s
    end
    # The limb of the smaller target, in the frame in which the larger one is the unit
    # sphere at the origin of its own axes.
    smlmat = transpose(rmat[smlidx])
    limb = edlimb(r[smlidx][1], r[smlidx][2], r[smlidx][3], _m_mxv(smlmat, viewpt - ctrs[smlidx]))
    smlctr = ctrs[smlidx] + _m_mtxv(smlmat, limb.center)
    smlmaj = _m_mtxv(smlmat, limb.semimajor)
    smlmin = _m_mtxv(smlmat, limb.semiminor)
    rb = r[bigidx]
    sclmat = SMatrix{3,3,Float64,9}(1.0 / rb[1], 0.0, 0.0, 0.0, 1.0 / rb[2], 0.0, 0.0, 0.0, 1.0 / rb[3])
    sclrot = _m_mxmt(sclmat, rmat[bigidx])
    xview = _m_mxv(sclrot, viewpt)
    bigctr = _m_mxv(sclrot, ctrs[bigidx])
    smlctr = _m_mxv(sclrot, smlctr)
    smlmaj, smlmin = _m_saelgv(_m_mxv(sclrot, smlmaj), _m_mxv(sclrot, smlmin))
    xlimb = _m_cgv2el(smlctr, smlmaj, smlmin)
    raydir = bigctr - xview
    minsep, minpt = _geo_asryel(true, xlimb, xview, raydir)
    bigr = _m_norm(raydir)
    uasize = _geo_dasine(1.0 / bigr, 1e-12)
    minsep > uasize && return 0
    sfront = vsep(minpt - bigctr, -raydir) <= _M_HALFPI - uasize
    s = ((sfront && smlidx == 1) || (!sfront && smlidx == 2)) ? 1 : -1
    minsep <= -uasize && return sfront ? 3s : 2s
    xsmlvu = xview - smlctr
    xdist = _m_norm(xsmlvu)
    majlen = _m_norm(smlmaj)
    minlen = _m_norm(smlmin)
    maxang = xdist >= majlen ? _geo_dasine(majlen / xdist, 1e-12) :
             _geo_occ_maxang(xview, smlctr, smlmaj, smlmin)
    tilt = vsep(_m_ucrss(smlmaj, smlmin), xsmlvu)
    tilt > _M_HALFPI && (tilt = Float64(pi) - tilt)
    minang = atan(cos(tilt) * minlen, sin(tilt) * minlen + xdist)
    xasep = vsep(raydir, smlctr - xview)
    xasep + maxang <= uasize && return sfront ? 2s : 3s
    xasep + minang > uasize && return s
    maxsep, _ = _geo_asryel(false, xlimb, xview, raydir)
    abs(maxsep) <= uasize && return sfront ? 2s : 3s
    return s
end

"The upper bound of the angular radius of one ellipsoid of `_geo_occed`."
@inline function _geo_occ_bound(r::NTuple{3,Float64}, xtpos::Vec3{Float64}, dist::Float64,
        maxrad::Float64)
    dist >= maxrad && return _geo_dasine(maxrad / dist, 1e-12)
    xvwtrg = -xtpos
    limb = edlimb(r[1], r[2], r[3], xvwtrg)
    return _geo_occ_maxang(xvwtrg, limb.center, limb.semimajor, limb.semiminor)
end

"""
    _GeoOccBody

One target of `occult`: its ID, the kind of its shape (`_GEO_OCC_POINT`,
`_GEO_OCC_ELLIPSOID`, `_GEO_OCC_DSK`), the parsed method string and the body-fixed
frame (0 for a point).
"""
struct _GeoOccBody
    id::Int32
    kind::UInt8
    m::_GeoMethod
    frame::Int32
end

"`ljust` and `ucase` of a shape string; no allocation for the two common spellings."
@inline function _geo_occ_upper(shape::AbstractString)
    (shape == "ELLIPSOID" || shape == "POINT") && return shape
    return _geo_occ_upper_slow(shape)
end
@noinline _geo_occ_upper_slow(shape::AbstractString) = uppercase(rstrip(lstrip(shape, ' '), ' '))

@noinline _geo_occ_bad_frame(frame, id, why) = throw(SpiceError(:invalidframe,
    "the body-fixed frame <$frame> of the target $id $why"))

"The frame ID of an extended target after the checks of `zzgfocin` (`SpiceError(:invalidframe)`)."
@inline function _geo_occ_frame(ks::KernelSet, frame::NameOrID, id::Int32)
    frame isa AbstractString && all(==(' '), frame) &&
        _geo_occ_bad_frame(frame, id, "is blank; an extended target needs a frame")
    ffrmid = _geo_namfrm(ks, frame)
    center, _, found = _geo_frinfo(ks, ffrmid)
    found || _geo_occ_bad_frame(frame, id, "is not recognized")
    center == id || _geo_occ_bad_frame(frame, id, "is not centered on the target")
    return ffrmid
end

@noinline _geo_occ_bad_shape(shape) = throw(SpiceError(:invalidshape,
    "the shape <$shape> is not supported; the shapes are POINT, ELLIPSOID and DSK/UNPRIORITIZED"))

"""
    _geo_occ_eval(ks, tols, front, fframe, back, bframe, bsrf, obs, corr, et) -> Int

The rest of `zzgfocin` (radii, frames, bounding spheres) and `zzgfocst` for the roles
`front` and `back`. `corr` is the correction without stellar aberration and `bsrf` the
surface list for the bounding sphere of the back target. The result is 0 if `front`
does not occult `back`; otherwise 1, 2 or 3 for a partial, an annular or a total
occultation of two ellipsoids, 3 for an extended `front` and a point `back`, and 2 for
a point `front` and an extended `back` (the codes that `occult` multiplies by the sign
of the roles).
"""
function _geo_occ_eval(ks::KernelSet, tols::DSKTolerances, front::_GeoOccBody, fframe,
        back::_GeoOccBody, bframe, bsrf, obs::Int32, corr::Aberration, et::Float64)
    # Radii, frame and bounding spheres, the front target first.
    frad, fmn, fmx, ffrm = _geo_occ_extent(ks, front, fframe, front.m.srflst)
    brad, bmn, bmx, bfrm = _geo_occ_extent(ks, back, bframe, bsrf)
    frtpos, ltfrnt = _geo_spkezp(ks, front.id, et, Int32(1), corr, obs)
    bckpos, ltback = _geo_spkezp(ks, back.id, et, Int32(1), corr, obs)
    s = _geo_sign(corr)
    if front.kind == _GEO_OCC_ELLIPSOID && back.kind == _GEO_OCC_ELLIPSOID
        mb = frame_rotation(ks, bfrm, Int32(1), et + s * ltback)
        mf = frame_rotation(ks, ffrm, Int32(1), et + s * ltfrnt)
        bsmaxs = hcat(brad[1] * mb[:, 1], brad[2] * mb[:, 2], brad[3] * mb[:, 3])
        fsmaxs = hcat(frad[1] * mf[:, 1], frad[2] * mf[:, 2], frad[3] * mf[:, 3])
        occode = _geo_occed(zero(Vec3{Float64}), bckpos, bsmaxs, frtpos, fsmaxs)
        return occode < 0 ? -occode : 0
    end
    pntocc = back.kind == _GEO_OCC_POINT      # the point is the back target
    bdist = _m_norm(bckpos)
    fdist = _m_norm(frtpos)
    bckfrt = frtpos - bckpos
    mnr = pntocc ? fmn : bmn
    ((pntocc ? fdist : bdist) <= mnr || (pntocc ? bdist : fdist) == 0.0 ||
     _m_norm(bckfrt) <= mnr) && throw(SpiceError(:notdisjoint,
        "the observer or the point target is inside the extended target, or the point " *
        "target coincides with the observer"))
    trgsep = vsep(bckpos, frtpos)
    srad = (pntocc ? fmx : bmx) * 1.01
    tdist = pntocc ? fdist : bdist
    srad < tdist && trgsep > _geo_dasine(srad / tdist, 1e-12) && return 0
    minang = _geo_dasine(mnr * 0.99 / tdist, 1e-12)
    hit = pntocc ? 3 : 2
    if trgsep < minang
        # The point is in the cone of the inner sphere: only the order along the line
        # of sight decides.
        pntocc && return vsep(-frtpos, -bckfrt) > _M_HALFPI - minang ? hit : 0
        return vsep(-bckpos, bckfrt) < _M_HALFPI - minang ? hit : 0
    end
    ext = pntocc ? front : back
    sh = _geo_shape(ks, tols, ext.m, ext.id, pntocc ? ffrm : bfrm)
    found, _, _, srfvec = _geo_sfxcor(ks, sh, et, corr, obs, Int32(1), pntocc ? bckpos : frtpos)
    found || return 0
    return (pntocc ? _m_norm(srfvec) < bdist : _m_norm(srfvec) > fdist) ? hit : 0
end

"Radii, radii of the inner and the outer sphere, and frame of one target of `occult`."
@inline function _geo_occ_extent(ks::KernelSet, b::_GeoOccBody, frame, srflst)
    b.kind == _GEO_OCC_POINT && return (0.0, 0.0, 0.0), 0.0, 0.0, Int32(0)
    if b.kind == _GEO_OCC_ELLIPSOID
        rad = _geo_radii(ks, b.id)
        ffrm = _geo_occ_frame(ks, frame, b.id)
        return rad, min(min(rad[1], rad[2]), rad[3]), max(max(rad[1], rad[2]), rad[3]), ffrm
    end
    ffrm = _geo_occ_frame(ks, frame, b.id)
    mn, mx = _geo_dsksph(ks, b.id, srflst)
    return (0.0, 0.0, 0.0), mn, mx, ffrm
end

"The first step of `occult` for one target: its ID (0 for a point whose name is read later) and its method."
@inline function _geo_occ_parse(ks::KernelSet, targ::NameOrID, shape::AbstractString, what)
    shap = _geo_occ_upper(shape)
    shap == "POINT" && return _GeoOccBody(Int32(0), _GEO_OCC_POINT, _geo_ell(0x00), Int32(0))
    id = _geo_body(ks, targ, what)
    m = _geo_method(ks, id, shap)
    kind = m.dsk ? _GEO_OCC_DSK : (shap == "ELLIPSOID" || _dsk_eqstr(shap, "ELLIPSOID")) ?
           _GEO_OCC_ELLIPSOID : 0x00
    return _GeoOccBody(id, kind, m, Int32(0))
end

function _geo_occult(ks::KernelSet, tols::DSKTolerances, targ1::NameOrID,
        shape1::AbstractString, frame1::NameOrID, targ2::NameOrID, shape2::AbstractString,
        frame2::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, et::Real)
    b1 = _geo_occ_parse(ks, targ1, shape1, "first target")
    b2 = _geo_occ_parse(ks, targ2, shape2, "second target")
    b1 = _GeoOccBody(_geo_body(ks, targ1, "first target"), b1.kind, b1.m, Int32(0))
    b2 = _GeoOccBody(_geo_body(ks, targ2, "second target"), b2.kind, b2.m, Int32(0))
    obs = _geo_body(ks, obsrvr, "observer")
    e = Float64(check_epoch(et))
    # `zzgfocin` makes its body and shape checks before it reads the correction.
    c1 = _geo_occ_checked(ks, tols, b1, shape1, frame1, b2, shape2, frame2, obs, abcorr, e)
    c1 != 0 && return c1
    c2 = _geo_occ_checked(ks, tols, b2, shape2, frame2, b1, shape1, frame1, obs, abcorr, e)
    return -c2
end

"One role assignment of `occult`: the checks of `zzgfocin`, then the state."
@inline function _geo_occ_checked(ks::KernelSet, tols::DSKTolerances, front::_GeoOccBody,
        fshape, fframe, back::_GeoOccBody, bshape, bframe, obs::Int32, abcorr, et::Float64)
    (front.id == back.id || front.id == obs || back.id == obs) && throw(SpiceError(
        :bodiesnotdistinct, "the front target $(front.id), the back target $(back.id) " *
        "and the observer $obs must be three different bodies"))
    front.kind == 0x00 && _geo_occ_bad_shape(fshape)
    back.kind == 0x00 && _geo_occ_bad_shape(bshape)
    return _geo_occ_roles(ks, tols, front, fshape, fframe, back, bshape, bframe, obs, abcorr, et)
end

function _geo_occ_roles(ks::KernelSet, tols::DSKTolerances, front::_GeoOccBody, fshape,
        fframe, back::_GeoOccBody, bshape, bframe, obs::Int32, abcorr, et::Float64)
    bsrf = (back.kind == _GEO_OCC_DSK && !isempty(back.m.srflst)) ?
           _geo_method(ks, front.id, _geo_occ_upper(bshape)).srflst : back.m.srflst
    (front.kind == _GEO_OCC_POINT && back.kind == _GEO_OCC_POINT) && throw(SpiceError(
        :invalidshapecombo, "both targets are points; one target must be extended"))
    ((front.kind == _GEO_OCC_DSK && back.kind != _GEO_OCC_POINT) ||
     (back.kind == _GEO_OCC_DSK && front.kind != _GEO_OCC_POINT)) && throw(SpiceError(
        :invalidshapecombo, "a target with DSK data needs a point as the other target"))
    full = _geo_abcorr(abcorr)
    corr = full.uselt ? Aberration(; lt=true, converged=full.converged, transmit=full.transmit) :
           full
    return _geo_occ_eval(ks, tols, front, fframe, back, bframe, bsrf, obs, corr, et)
end

"""
    occult([ctx,] targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, et) -> Int

Return the occultation state of two targets as seen by `obsrvr` at `et` (seconds past
J2000 TDB), as `occult_c` of CSPICE:

| Code | Meaning |
|---|---|
| -3 | total occultation of the first target by the second |
| -2 | annular occultation of the first target by the second |
| -1 | partial occultation of the first target by the second |
| 0 | no occultation |
| 1 | partial occultation of the second target by the first |
| 2 | annular occultation of the second target by the first |
| 3 | total occultation of the second target by the first |

`shape1` and `shape2` are `"ELLIPSOID"`, `"POINT"` or `"DSK/UNPRIORITIZED"`; the DSK
form can end with `"/SURFACES = <list>"`. At least one target must be extended, and a
DSK target needs a point as the other target. With a point and an extended target, the
only codes are 0, ±2 (the point is in front) and ±3 (the point is behind).

`frame1` and `frame2` are the body-fixed frames of the extended targets, each centered
on its target. The frame of a point target is not used; pass `" "`.

`abcorr` is a correction string of [`spkezr`](@ref) or an [`Aberration`](@ref).
Stellar aberration is not applied: `"LT+S"` acts as `"LT"`, as in CSPICE. `et` cannot
be a `ForwardDiff.Dual`, because the result is discrete. The call does not allocate when
the bodies and the frames are integer IDs and the shapes are ellipsoids or points.
A DSK shape uses the tolerances that [`dskstl`](@ref) set on the context.

Throws `SpiceError` with the kind
- `:idcodenotfound` if a body name is not known;
- `:bodiesnotdistinct` if two of the three bodies are the same;
- `:invalidshape` if a shape is not one of the three values, `:badmethodsyntax` if a
  DSK shape string is malformed, `:badpriorityspec` if it does not have `UNPRIORITIZED`;
- `:invalidshapecombo` if both targets are points, or if a DSK target is paired with
  an extended target;
- `:invalidoption` if `abcorr` is not a valid correction;
- `:invalidframe` if the frame of an extended target is blank, not known or not
  centered on the target;
- `:notdisjoint` if the observer is inside a target, or if the targets overlap;
- the kinds of [`srfrec`](@ref) if the radii of an ellipsoid target are missing or not
  valid, `:dskdatanotfound` if a DSK target has no data loaded.

Also throws the errors of [`spkpos`](@ref).

# Examples
```julia
# Needs an SPK, a PCK with the radii and the rotation of the Moon and the Sun, and a
# leap seconds kernel.
et = str2et("2024 APR 08 18:18:00 UTC")
code = occult("MOON", "ELLIPSOID", "IAU_MOON", "SUN", "ELLIPSOID", "IAU_SUN", "LT", "EARTH", et)
code > 0        # true if the Moon hides a part of the Sun as seen from the Earth center
code == 3       # true if the eclipse is total there
```

See also [`gfoclt`](@ref), [`fovtrg`](@ref).
"""
occult(ks::KernelSet, targ1::NameOrID, shape1::AbstractString, frame1::NameOrID,
    targ2::NameOrID, shape2::AbstractString, frame2::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, et::Real) = _geo_occult(ks, _DSK_DEFAULT_TOLS, targ1, shape1, frame1,
    targ2, shape2, frame2, abcorr, obsrvr, et)
occult(ctx::SpiceContext, targ1::NameOrID, shape1::AbstractString, frame1::NameOrID,
    targ2::NameOrID, shape2::AbstractString, frame2::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, et::Real) = _geo_occult(kernelset(ctx), dsktolerances(ctx), targ1,
    shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, et)
occult(targ1::NameOrID, shape1::AbstractString, frame1::NameOrID, targ2::NameOrID,
    shape2::AbstractString, frame2::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, et::Real) =
    occult(default_context(), targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, et)
