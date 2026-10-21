# Limb points, terminator points and the tangent point of a ray: `limbpt`, `termpt` and
# `tangpt` of CSPICE N0067, with `zztangnt`, `zztanslv`, `zztanutl` (tangent rays on DSK
# surfaces) and `zzedtmpt` (terminator point of an ellipsoid).
#
# The outputs are lists, so the calls allocate the result vectors; the numeric work on
# each point does not allocate. The epoch is a plain number.

const _GEO_LOC_CENTER = 1
const _GEO_LOC_ELLIPSOID = 2
"Size of the result window and of the point buffer of `limbpt` and `termpt` in CSPICE."
const _GEO_MAXWIN = 2000

"""
    _geo_radii3(ks, body, kind) -> NTuple{3,Float64}

`zzbodvcd` for `RADII` with the count test of `limbpt`, `termpt` and `tangpt`: three
values, with no test of their signs. Throws `SpiceError` (`:kernelvarnotfound`,
`:typemismatch`, `:arraytoosmall`, and `kind` for a count that is not 3).
"""
@inline function _geo_radii3(ks::KernelSet, body::Integer, kind::Symbol)
    v = get(ks.pool, 1 <= body <= 999 ? (@inbounds _GEO_RADII_KEYS[body]) :
            string("BODY", body, "_RADII"), nothing)
    v isa Vector{Float64} || _geo_bad_radii3(body, v, kind)
    length(v) == 3 || _geo_bad_radii3(body, v, kind)
    return @inbounds (v[1], v[2], v[3])
end
@noinline function _geo_bad_radii3(body, v, kind::Symbol)
    v === nothing && throw(SpiceError(:kernelvarnotfound,
        "the variable BODY$(body)_RADII could not be found in the kernel pool"))
    v isa Vector{Float64} || throw(SpiceError(:typemismatch,
        "the variable BODY$(body)_RADII is not numeric"))
    throw(SpiceError(length(v) > 3 ? :arraytoosmall : kind,
        "the variable BODY$(body)_RADII has $(length(v)) values; 3 are required"))
end

# --- terminator point of an ellipsoid ---------------------------------------------------

"""
    _geo_edtmpt(umbral, a, b, c, r, axis, plnvec) -> Vec3

`zzedtmpt`: the terminator point of the ellipsoid with the semi-axes `a`, `b`, `c` in
the half plane that the vector `axis` (from the center to a spherical source of the
radius `r`) bounds and that contains `plnvec`. `umbral` selects the umbral terminator,
otherwise the penumbral one.

Contract: positive axes (`SpiceError(:invalidaxislength)`) and radius
(`:invalidradius`), `axis` and `plnvec` not zero (`:zerovector`) and not parallel
(`:degeneratecase`), and a source that does not touch the bounding sphere of the
ellipsoid (`:objectstooclose`). The iteration of CSPICE: at most 20 steps, to an angular
error of `1e-15` (`:noconvergence`).
"""
function _geo_edtmpt(umbral::Bool, a::Float64, b::Float64, c::Float64, r::Float64,
        axis::Vec3{Float64}, plnvec::Vec3{Float64})
    (a <= 0.0 || b <= 0.0 || c <= 0.0) && throw(SpiceError(:invalidaxislength,
        "the target radii must be positive; they are $a, $b, $c"))
    r <= 0.0 && throw(SpiceError(:invalidradius, "the radius of the source must be positive; it is $r"))
    _m_iszero(axis) && throw(SpiceError(:zerovector, "the axis vector is the zero vector"))
    r + max(max(a, b), c) >= _m_norm(axis) && throw(SpiceError(:objectstooclose,
        "the source touches the bounding sphere of the target: the distance between the " *
        "centers is $(_m_norm(axis))"))
    _m_iszero(plnvec) && throw(SpiceError(:zerovector, "the reference vector is the zero vector"))
    # In the scaled space the target is the unit sphere and the source is an ellipsoid.
    ta = 1.0 / a
    tb = 1.0 / b
    tc = 1.0 / c
    xa = ta * r
    xb = tb * r
    xc = tc * r
    tplnvc = Vec3(ta * plnvec[1], tb * plnvec[2], tc * plnvec[3])
    taxis = Vec3(ta * axis[1], tb * axis[2], tc * axis[3])
    utaxis = _m_hat(taxis)
    hplnml = _m_cross(tplnvc, taxis)
    _m_iszero(hplnml) && throw(SpiceError(:degeneratecase,
        "the reference vector and the axis vector are parallel"))
    maxr = max(max(xa, xb), xc)
    d = _m_norm(taxis)
    if umbral
        angle = _geo_dasine((maxr - 1.0) / d, 1e-14)
        targpt = _m_vrotv(utaxis, hplnml, -(_M_HALFPI + angle))
        s = 1.0
    else
        angle = _geo_dasine((maxr + 1.0) / d, 1e-14)
        targpt = _m_vrotv(utaxis, hplnml, angle - _M_HALFPI)
        s = -1.0
    end
    normal = _m_hat(targpt)
    konst = _m_dot(normal, targpt)
    srcpnt = ednmpt(xa, xb, xc, s * normal) + taxis
    h = _m_dot(srcpnt, normal) - konst
    angerr = -h / d
    nitr = 0
    # invariant: `targpt` is on the unit sphere in the half plane; `angerr` is the angle
    # by which the tangent plane at `targpt` misses the source
    while abs(angerr) > 1e-15 && nitr <= 20
        targpt = _m_vrotv(targpt, hplnml, angerr)
        normal = _m_hat(targpt)
        konst = _m_dot(normal, targpt)
        srcpnt = ednmpt(xa, xb, xc, s * normal) + taxis
        h = _m_dot(srcpnt, normal) - konst
        d = _m_norm(_m_perp(srcpnt, hplnml, false) - targpt)
        angerr = -h / d
        nitr += 1
        nitr > 20 && throw(SpiceError(:noconvergence,
            "the terminator point iteration did not converge in 20 steps"))
    end
    return Vec3(a * targpt[1], b * targpt[2], c * targpt[3])
end

# --- tangent rays on a DSK surface ------------------------------------------------------

"""
    _geo_tansta(ks, sh, curve, srcrad, et, vrtx, uaxis, nrml, angle) -> (hit, point)

`zztansta`: the state of the ray at the angle `angle` in the cutting half plane, and
its surface intercept. For the limb (`curve == 0`) the ray leaves the vertex; for the
umbral (`1`) and the penumbral (`2`) terminator it is tangent to the source sphere.
"""
@inline function _geo_tansta(ks::KernelSet, sh::_GeoShape, curve::Int, srcrad::Float64,
        et::Float64, vrtx::Vec3{Float64}, uaxis::Vec3{Float64}, nrml::Vec3{Float64},
        angle::Float64)
    if curve == 0
        return _geo_rayx(ks, sh, vrtx, _m_vrotv(uaxis, nrml, angle), et)
    elseif curve == 1
        raydir = _m_vrotv(uaxis, nrml, angle)
        apex = vrtx + srcrad * _m_vrotv(uaxis, nrml, angle - Float64(pi) / 2)
        return _geo_rayx(ks, sh, apex, raydir, et)
    end
    raydir = _m_vrotv(uaxis, nrml, -angle)
    apex = vrtx + srcrad * _m_vrotv(uaxis, nrml, Float64(pi) * 1.5 - angle)
    return _geo_rayx(ks, sh, apex, raydir, et)
end

@noinline _geo_tan_error(kind::Symbol, msg) = throw(SpiceError(kind, msg))

"""
    _geo_tangnt!(ends, pts, ks, sh, curve, srcrad, et, plnvec, axis, schstp, soltol) -> Int

`zztangnt` with `zztanslv`: the tangent points of rays on the surface `sh` in the half
plane that `axis` bounds and that contains `plnvec`. The points are left in `pts[1:n]`
and `n` is returned; `ends` is scratch space for the transition angles.

Contract: the ray angle goes from `start` to `finish` in steps of `schstp`; every
change of the state (the ray meets the surface or not) is located by bisection to the
tolerance `soltol`, and the intercept that was found last is the tangent point. Throws
`SpiceError`: `:zerovector`, `:degeneratecase`, `:invalidradius`, `:objectstooclose`,
`:invalidtolerance`, `:invalidconststep`, `:noconvergence`, `:arraytoosmall` for more
than $(_GEO_MAXWIN) points, and `:bug` where CSPICE signals `SPICE(BUG)` (a transition
without an intercept in the refinement).
"""
function _geo_tangnt!(ends::Vector{Float64}, pts::Vector{Vec3{Float64}}, ks::KernelSet,
        sh::_GeoShape, curve::Int, srcrad::Float64, et::Float64, plnvec::Vec3{Float64},
        axis::Vec3{Float64}, schstp::Float64, soltol::Float64)
    empty!(ends)
    empty!(pts)
    refvec = curve == 2 ? _m_vrotv(plnvec, axis, Float64(pi)) : plnvec
    # `zztanini`
    _m_iszero(axis) && _geo_tan_error(:zerovector, "the axis vector is the zero vector")
    _m_iszero(refvec) && _geo_tan_error(:zerovector, "the reference vector is the zero vector")
    (curve != 0 && srcrad <= 0.0) && _geo_tan_error(:badsourceradius,
        "the radius of the source must be positive; it is $srcrad")
    nrml = _m_cross(axis, refvec)
    _m_iszero(nrml) && _geo_tan_error(:degeneratecase,
        "the reference vector and the axis vector are parallel")
    nrml = _m_hat(nrml)
    uaxis = _m_hat(axis)
    maxrad = sh.maxrad * 1.001
    maxrad <= 0.0 && _geo_tan_error(:invalidradius, "the radius of the bounding sphere is not positive")
    r = _m_norm(axis)
    pival = Float64(pi)
    if curve == 0
        start = 0.0
        r > maxrad * 1.0000000000010001 && (start = pival - _geo_dasine(maxrad / r, 1e-12))
        finish = pival
    else
        srcrad + maxrad > r && _geo_tan_error(:objectstooclose,
            "the source touches the bounding sphere of the target")
        if curve == 1
            start = pival + _geo_dasine((srcrad - maxrad) / r, 1e-12)
            finish = pival + _geo_dasine(srcrad / r, 1e-12)
        else
            alpha = srcrad / (srcrad + maxrad)
            start = pival - _geo_dasine(srcrad / (alpha * r), 1e-12)
            finish = pival - _geo_dasine(srcrad / r, 1e-12)
        end
    end
    # `zztanslv`
    soltol <= 0.0 && _geo_tan_error(:invalidtolerance, "the tolerance must be positive; it is $soltol")
    start > finish && _geo_tan_error(:boundsoutoforder, "the search bounds are out of order")
    (start - soltol == start || start + soltol == start || finish - soltol == finish ||
     finish + soltol == finish) && _geo_tan_error(:invalidtolerance,
        "the tolerance $soltol is too small for the search bounds")
    schstp <= 0.0 && _geo_tan_error(:invalidconststep, "the search step must be positive; it is $schstp")
    maxmag = max(abs(start), abs(finish))
    maxmag + schstp == maxmag && _geo_tan_error(:invalidconststep,
        "the search step $schstp is too small for the search bounds")
    curx = start
    cursta, xpoint = _geo_tansta(ks, sh, curve, srcrad, et, axis, uaxis, nrml, curx)
    prvpnt = xpoint
    prvset = cursta
    instat = cursta
    openstart = !cursta                # `endflg[0]`
    if cursta
        # The start of the interval is an end point without a tangent point.
        push!(ends, curx)
        push!(pts, xpoint)
    end
    svdx = curx
    savsta = cursta
    # invariant: `ends` holds the transition angles in [start, svdx] and `pts` the
    # intercept found last before each of them
    while svdx < finish
        curx = min(curx + schstp, finish)
        cursta, xpoint = _geo_tansta(ks, sh, curve, srcrad, et, axis, uaxis, nrml, curx)
        cursta && (prvpnt = xpoint; prvset = true)
        while savsta == cursta && svdx < finish
            svdx = curx
            savsta = cursta
            curx = min(curx + schstp, finish)
            cursta, xpoint = _geo_tansta(ks, sh, curve, srcrad, et, axis, uaxis, nrml, curx)
            cursta && (prvpnt = xpoint; prvset = true)
        end
        if savsta != cursta
            state1 = savsta
            x1 = svdx
            x2 = curx
            svdx = curx
            savsta = cursta
            nloop = 0
            # invariant: the state is `state1` at x1 and the other state at x2
            while x2 - x1 > soltol
                nloop += 1
                nloop >= 1000 && _geo_tan_error(:noconvergence,
                    "the bisection of a state change did not converge in 1000 steps")
                t = max(x1, min(x2, x1 * 0.5 + x2 * 0.5))
                if t == x1
                    x2 = t
                elseif t == x2
                    x1 = t
                else
                    st, xpoint = _geo_tansta(ks, sh, curve, srcrad, et, axis, uaxis, nrml, t)
                    st && (prvpnt = xpoint; prvset = true)
                    st == state1 ? (x1 = t) : (x2 = t)
                end
            end
            trnstn = max(x1, min(x2, (x1 + x2) * 0.5))
            instat = !instat
            length(ends) >= _GEO_MAXWIN && _geo_tan_error(:arraytoosmall,
                "more than $_GEO_MAXWIN tangent points were found on one cutting half plane")
            prvset || _geo_tan_error(:bug, "a state change was found without a surface " *
                "intercept in its refinement, as SPICE(BUG) of CSPICE")
            push!(ends, trnstn)
            push!(pts, prvpnt)
            prvset = false
        end
    end
    # `zztangnt`: the ends of the search interval are not tangent points.
    if instat
        push!(ends, finish)
        push!(pts, prvpnt)
    end
    if !isempty(ends) && ends[1] == start && !openstart
        popfirst!(ends)
        popfirst!(pts)
    end
    if !isempty(ends) && ends[end] == finish && instat
        pop!(ends)
        pop!(pts)
    end
    return length(pts)
end

# --- limb points ------------------------------------------------------------------------

@noinline _geo_bad_locus(corloc, values) = throw(SpiceError(:invalidlocus,
    "the aberration correction locus <$corloc> is not recognized; the values are $values"))

"The locus code of `corloc` (`ljucrs`: case, leading blanks and repeated blanks are ignored)."
@inline function _geo_locus(corloc::AbstractString, second::String)
    corloc == "CENTER" && return _GEO_LOC_CENTER
    corloc == second && return _GEO_LOC_ELLIPSOID
    return _geo_locus_slow(corloc, second)
end
@noinline function _geo_locus_slow(corloc::AbstractString, second::String)
    n = normalize_name(corloc)
    n == "CENTER" && return _GEO_LOC_CENTER
    n == second && return _GEO_LOC_ELLIPSOID
    _geo_bad_locus(corloc, "CENTER and $second")
end

"The limb plane normal that points toward the observer, and the limb center."
@inline function _geo_limb_frame(radii::NTuple{3,Float64}, axis::Vec3{Float64})
    limb = edlimb(radii[1], radii[2], radii[3], axis)
    enorml = _m_ucrss(limb.semimajor, limb.semiminor)
    _m_dot(enorml, axis) < 0.0 && (enorml = -1.0 * enorml)
    return limb.center, enorml
end

@noinline _geo_cut_degenerate(i) = throw(SpiceError(:degeneratecase,
    "the cutting half plane $i is parallel to the limb plane"))
@noinline _geo_limb_bug() = throw(SpiceError(:bug,
    "no intercept of the limb ray with the ellipsoid was found"))
@noinline _geo_out_of_room(i, n, room) = throw(SpiceError(:outofroom,
    "out of room in the output arrays: the cutting half plane $i has $n points and the " *
    "room left is $room"))

"The limb point of the ellipsoid in the half plane with the normal `cutnml`."
@inline function _geo_limb_point(radii::NTuple{3,Float64}, center::Vec3{Float64},
        enorml::Vec3{Float64}, cutnml::Vec3{Float64}, i::Int)
    edir = _m_ucrss(cutnml, enorml)
    _m_iszero(edir) && _geo_cut_degenerate(i)
    found, pt = _m_surfpt(center, edir, radii[1], radii[2], radii[3])
    found || _geo_limb_bug()
    return pt, edir
end

function _geo_limbpt(ks::KernelSet, tols::DSKTolerances, method::AbstractString,
        target::NameOrID, et0::Real, fixref::NameOrID, abcorr::AbcorrArg,
        corloc::AbstractString, obsrvr::NameOrID, refvec0::AbstractVector, rolstp0::Real,
        ncuts::Integer, schstp0::Real, soltol0::Real, maxn::Integer)
    corr = _geo_abcorr(abcorr)
    s = _geo_sign(corr)
    trgcde = _geo_body(ks, target, "target")
    obscde = _geo_body(ks, obsrvr, "observer")
    obscde == trgcde && _geo_same_bodies("limb points", obsrvr)
    fixfid = _geo_fixfid(ks, fixref, trgcde, target)
    m = _geo_method(ks, trgcde, method)
    m.pntdef == 0x00 && throw(SpiceError(:invalidlimbtype,
        "the method string <$method> has no limb type; use TANGENT or GUIDED"))
    (m.subtyp != 0x00 || m.trmtyp != 0x00) && throw(SpiceError(:invalidmethod,
        "the method string <$method> has a sub-point type or a terminator type, which " *
        "is not valid here"))
    tangent = m.pntdef == _GEO_TANGENT
    loccde = _geo_locus(corloc, "ELLIPSOID LIMB")
    refvec = _geo_pos(refvec0)
    _m_iszero(refvec) && throw(SpiceError(:zerovector, "the reference vector is the zero vector"))
    maxn < 1 && throw(SpiceError(:invalidsize, "maxn = $maxn; it must be at least 1"))
    (ncuts < 1 || ncuts > maxn) && throw(SpiceError(:invalidcount,
        "ncuts = $ncuts; it must be in 1:$maxn"))
    rolstp, schstp, soltol = Float64(rolstp0), Float64(schstp0), Float64(soltol0)
    if m.dsk
        schstp <= 0.0 && throw(SpiceError(:invalidsearchstep,
            "the search step must be positive; it is $schstp"))
        soltol <= 0.0 && throw(SpiceError(:invalidtolerance,
            "the solution tolerance must be positive; it is $soltol"))
    end
    (ncuts > 1 && rolstp == 0.0) && throw(SpiceError(:invalidrollstep,
        "the roll step is zero and more than one cutting half plane is asked for"))
    sh = m.dsk ? _geo_dskshape(ks, tols, trgcde, fixfid, m.srflst) :
         _GeoShape(false, trgcde, fixfid, 1.0, 1.0, 1.0, 1.0, 1.0, _GEO_NOSURF, tols)
    radii = (!m.dsk || loccde == _GEO_LOC_ELLIPSOID || !tangent) ?
            _geo_radii3(ks, trgcde, :badradiuscount) : (0.0, 0.0, 0.0)
    et = Float64(check_epoch(et0))
    npts = zeros(Int, ncuts)
    points = Vec3{Float64}[]
    epochs = Float64[]
    tangts = Vec3{Float64}[]
    ends = Float64[]
    pntbuf = Vec3{Float64}[]
    pos, lt = _geo_spkezp(ks, trgcde, et, fixfid, corr, obscde)
    _m_iszero(pos) && throw(SpiceError(:noseparation,
        "the distance between the observer and the target is zero"))
    room = Int(maxn)
    if loccde == _GEO_LOC_CENTER
        trgepc = et + s * lt
        axis = -pos
        _m_iszero(_m_cross(axis, refvec)) && throw(SpiceError(:degeneratecase,
            "the reference vector is parallel to the observer-target vector"))
        center = zero(Vec3{Float64})
        enorml = zero(Vec3{Float64})
        (!m.dsk || !tangent) && ((center, enorml) = _geo_limb_frame(radii, axis))
        for i in 1:ncuts                 # invariant: the cuts 1 .. i-1 are in the outputs
            plnvec = _m_vrotv(refvec, axis, (i - 1) * rolstp)
            cutnml = _m_ucrss(axis, plnvec)
            if m.dsk && tangent
                n = _geo_tangnt!(ends, pntbuf, ks, sh, 0, 0.0, trgepc, plnvec, axis, schstp, soltol)
            elseif m.dsk
                edir = _m_ucrss(cutnml, enorml)
                _m_iszero(edir) && _geo_cut_degenerate(i)
                found, pt = _geo_rayx(ks, sh, 1.0 * center + (sh.maxrad * 3.0) * edir, -edir, trgepc)
                empty!(pntbuf)
                found && push!(pntbuf, pt)
                n = length(pntbuf)
            else
                pt, _ = _geo_limb_point(radii, center, enorml, cutnml, i)
                empty!(pntbuf)
                push!(pntbuf, pt)
                n = 1
            end
            npts[i] = n
            n > room && _geo_out_of_room(i, n, room)
            room -= n
            for j in 1:n
                push!(points, pntbuf[j])
                push!(tangts, pntbuf[j] - axis)
                push!(epochs, trgepc)
            end
        end
        return npts, points, epochs, tangts
    end
    tangent || throw(SpiceError(:badlimblocusmix,
        "the limb type GUIDED cannot be used with the locus ELLIPSOID LIMB"))
    robs, vobs = _ssb(ks, obscde, et, Int32(1))
    vstl = corr.transmit ? -vobs : vobs
    stloff = zero(Vec3{Float64})
    for i in 1:ncuts                     # invariant: the cuts 1 .. i-1 are in the outputs
        roll = (i - 1) * rolstp
        local axis, plnvec, epoint, epoch
        if corr.uselt
            numitr = corr.converged ? 5 : 2
            j = 0
            lterr = 1.0
            xform = one(Mat3{Float64})
            cortrg = zero(Vec3{Float64})
            # invariant: `lt` is the light time to the limb point after j iterations
            # (it starts from the value of the previous cut, as in CSPICE)
            while j < numitr && lterr > 1e-18
                epoch = et + s * lt
                ptarg = _ssbpos(ks, trgcde, epoch, Int32(1)) - robs
                cortrg = !corr.stellar ? ptarg : j == 0 ? _stelab(ptarg, vstl) : ptarg + stloff
                xform = frame_rotation(ks, Int32(1), fixfid, epoch)
                axis = -_m_mxv(xform, cortrg)
                plnvec = _m_vrotv(refvec, axis, roll)
                center, enorml = _geo_limb_frame(radii, axis)
                epoint, _ = _geo_limb_point(radii, center, enorml, _m_ucrss(axis, plnvec), i)
                isrfvc = _m_mtxv(xform, epoint) + ptarg
                corr.stellar && (stloff = _stelab(isrfvc, vstl) - isrfvc)
                prvlt = lt
                lt = _m_norm(isrfvc) / CLIGHT
                lterr = abs(lt - prvlt) / max(1.0, abs(epoch))
                j += 1
            end
            axis = -_m_mxv(xform, cortrg)
            epoch = et + s * lt
        else
            axis = -pos
            epoch = et
            plnvec = _m_vrotv(refvec, axis, roll)
            center, enorml = _geo_limb_frame(radii, axis)
            epoint, _ = _geo_limb_point(radii, center, enorml, _m_ucrss(axis, plnvec), i)
        end
        if m.dsk
            n = _geo_tangnt!(ends, pntbuf, ks, sh, 0, 0.0, epoch, plnvec, axis, schstp, soltol)
        else
            empty!(pntbuf)
            push!(pntbuf, epoint)
            n = 1
        end
        npts[i] = n
        n > room && _geo_out_of_room(i, n, room)
        room -= n
        for j in 1:n
            push!(points, pntbuf[j])
            push!(tangts, pntbuf[j] - axis)
            push!(epochs, epoch)
        end
    end
    return npts, points, epochs, tangts
end

# --- terminator points ------------------------------------------------------------------

function _geo_termpt(ks::KernelSet, tols::DSKTolerances, method::AbstractString,
        ilusrc::NameOrID, target::NameOrID, et0::Real, fixref::NameOrID, abcorr::AbcorrArg,
        corloc::AbstractString, obsrvr::NameOrID, refvec0::AbstractVector, rolstp0::Real,
        ncuts::Integer, schstp0::Real, soltol0::Real, maxn::Integer)
    corr = _geo_abcorr(abcorr)
    corr.transmit && throw(SpiceError(:invalidoption,
        "the aberration correction $(string(corr)) calls for transmission; termpt " *
        "accepts reception corrections only"))
    trgcde = _geo_body(ks, target, "target")
    obscde = _geo_body(ks, obsrvr, "observer")
    ilucde = _geo_body(ks, ilusrc, "illumination source")
    obscde == trgcde && _geo_same_bodies("terminator points", obsrvr)
    ilucde == trgcde && throw(SpiceError(:bodiesnotdistinct,
        "the illumination source and the target are the same body: $ilusrc"))
    fixfid = _geo_fixfid(ks, fixref, trgcde, target)
    m = _geo_method(ks, trgcde, method)
    m.pntdef == 0x00 && throw(SpiceError(:invalidtermtype,
        "the method string <$method> has no terminator point definition; use TANGENT or GUIDED"))
    m.trmtyp == 0x00 && throw(SpiceError(:invalidshadow,
        "the method string <$method> has no shadow type; use UMBRAL or PENUMBRAL"))
    m.subtyp == 0x00 || _geo_sub_in_method(method)
    tangent = m.pntdef == _GEO_TANGENT
    umbral = m.trmtyp == _GEO_UMBRAL
    shadow = umbral ? 1 : 2
    loccde = _geo_locus(corloc, "ELLIPSOID TERMINATOR")
    refvec = _geo_pos(refvec0)
    _m_iszero(refvec) && throw(SpiceError(:zerovector, "the reference vector is the zero vector"))
    maxn < 1 && throw(SpiceError(:invalidsize, "maxn = $maxn; it must be at least 1"))
    (ncuts < 0 || ncuts > maxn) && throw(SpiceError(:invalidcount,
        "ncuts = $ncuts; it must be in 0:$maxn"))
    rolstp, schstp, soltol = Float64(rolstp0), Float64(schstp0), Float64(soltol0)
    (ncuts > 1 && rolstp == 0.0) && throw(SpiceError(:invalidrollstep,
        "the roll step is zero and more than one cutting half plane is asked for"))
    sh = m.dsk ? _geo_dskshape(ks, tols, trgcde, fixfid, m.srflst) :
         _GeoShape(false, trgcde, fixfid, 1.0, 1.0, 1.0, 1.0, 1.0, _GEO_NOSURF, tols)
    srad = _geo_radii3(ks, ilucde, :badradiuscount)
    ilurad = max(max(srad[1], srad[2]), srad[3])
    trad = (!m.dsk || loccde == _GEO_LOC_ELLIPSOID || !tangent) ?
           _geo_radii3(ks, trgcde, :badradiuscount) : (0.0, 0.0, 0.0)
    et = Float64(check_epoch(et0))
    npts = zeros(Int, ncuts)
    points = Vec3{Float64}[]
    epochs = Float64[]
    trmvcs = Vec3{Float64}[]
    ends = Float64[]
    pntbuf = Vec3{Float64}[]
    trgpos, lt = _geo_spkezp(ks, trgcde, et, fixfid, corr, obscde)
    _m_iszero(trgpos) && throw(SpiceError(:noseparation,
        "the distance between the observer and the target is zero"))
    room = Int(maxn)
    if loccde == _GEO_LOC_CENTER
        trgepc = et + _geo_sign(corr) * lt
        axis, _ = _geo_spkezp(ks, ilucde, trgepc, fixfid, corr, trgcde)
        _m_iszero(_m_cross(axis, refvec)) && throw(SpiceError(:degeneratecase,
            "the reference vector is parallel to the target-source vector"))
        for i in 1:ncuts                 # invariant: the cuts 1 .. i-1 are in the outputs
            plnvec = _m_vrotv(refvec, axis, (i - 1) * rolstp)
            if m.dsk && tangent
                n = _geo_tangnt!(ends, pntbuf, ks, sh, shadow, ilurad, trgepc, plnvec, axis,
                    schstp, soltol)
            elseif m.dsk
                edir = _m_hat(_geo_edtmpt(umbral, trad[1], trad[2], trad[3], ilurad, axis, plnvec))
                found, pt = _geo_rayx(ks, sh, (sh.maxrad * 3.0) * edir, -edir, trgepc)
                empty!(pntbuf)
                found && push!(pntbuf, pt)
                n = length(pntbuf)
            else
                empty!(pntbuf)
                push!(pntbuf, _geo_edtmpt(umbral, trad[1], trad[2], trad[3], ilurad, axis, plnvec))
                n = 1
            end
            npts[i] = n
            n > room && _geo_out_of_room(i, n, room)
            for j in 1:n
                push!(points, pntbuf[j])
                push!(trmvcs, pntbuf[j] + trgpos)
                push!(epochs, trgepc)
            end
            # As in CSPICE, the room goes down by the count for each point of the cut.
            room -= n * n
        end
        return npts, points, epochs, trmvcs
    end
    tangent || throw(SpiceError(:badtermlocusmix,
        "the point definition GUIDED cannot be used with the locus ELLIPSOID TERMINATOR"))
    robs, vobs = _ssb(ks, obscde, et, Int32(1))
    stloff = zero(Vec3{Float64})
    epoint = zero(Vec3{Float64})
    for i in 1:ncuts                     # invariant: the cuts 1 .. i-1 are in the outputs
        roll = (i - 1) * rolstp
        local axis, epoch
        if corr.uselt
            numitr = corr.converged ? 5 : 1
            j = 0
            lterr = 1.0
            # invariant: `lt` is the light time to the terminator point after j
            # iterations (it starts from the value of the previous cut, as in CSPICE)
            while j < numitr && lterr > 1e-17
                epoch = et - lt
                ptarg = _ssbpos(ks, trgcde, epoch, Int32(1)) - robs
                axis, _ = _geo_spkezp(ks, ilucde, epoch, fixfid, corr, trgcde)
                plnvec = _m_vrotv(refvec, axis, roll)
                epoint = _geo_edtmpt(umbral, trad[1], trad[2], trad[3], ilurad, axis, plnvec)
                trgpos = _m_mxv(frame_rotation(ks, Int32(1), fixfid, epoch), ptarg)
                prvlt = lt
                lt = _m_norm(epoint + trgpos) / CLIGHT
                lterr = abs(lt - prvlt) / max(1.0, abs(epoch))
                j += 1
            end
            epoch = et - lt
            ptarg = _ssbpos(ks, trgcde, epoch, Int32(1)) - robs
            xform = frame_rotation(ks, Int32(1), fixfid, epoch)
            if corr.stellar
                itrmvc = _m_mtxv(xform, epoint) + ptarg
                stloff = _m_mxv(xform, _stelab(itrmvc, vobs) - itrmvc)
            end
            trgpos = _m_mxv(xform, ptarg)
            axis, _ = _geo_spkezp(ks, ilucde, epoch, fixfid, corr, trgcde)
        else
            axis, _ = _geo_spkezp(ks, ilucde, et, fixfid, corr, trgcde)
            epoch = et
        end
        plnvec = _m_vrotv(refvec, axis, roll)
        if m.dsk
            n = _geo_tangnt!(ends, pntbuf, ks, sh, shadow, ilurad, epoch, plnvec, axis, schstp, soltol)
        else
            corr.uselt || (epoint = _geo_edtmpt(umbral, trad[1], trad[2], trad[3], ilurad, axis, plnvec))
            empty!(pntbuf)
            push!(pntbuf, epoint)
            n = 1
        end
        npts[i] = n
        n > room && _geo_out_of_room(i, n, room)
        for j in 1:n
            push!(points, pntbuf[j])
            push!(trmvcs, corr.stellar ? (pntbuf[j] + trgpos) + stloff : pntbuf[j] + trgpos)
            push!(epochs, epoch)
        end
        room -= n * n
    end
    return npts, points, epochs, trmvcs
end

"""
    limbpt([ctx,] method, target, et, fixref, abcorr, corloc, obsrvr, refvec, rolstp,
           ncuts, schstp, soltol, maxn) -> (npts, points, epochs, tangts)

Return limb points of `target` as seen by `obsrvr` at `et` (seconds past J2000 TDB),
as `limbpt_c` of CSPICE.

The points lie in `ncuts` half planes bounded by the line from the target center to
the observer. Half plane `i` contains `refvec` rotated about that line by
`(i - 1) * rolstp` radians. `refvec` has 3 elements, in the body-fixed frame `fixref`.

- `npts`: `npts[i]` is the number of points found in half plane `i` (`ncuts` elements).
- `points`: the limb points, in km in `fixref`, half plane after half plane.
- `epochs`: the epoch at each point, in seconds past J2000 TDB (see `corloc`).
- `tangts`: the vector from the observer to each point, in km in `fixref` at that epoch.

The three lists have `sum(npts)` elements, not `maxn` as the arrays of CSPICE. The
points of half plane `i` start at index `sum(npts[1:i-1]) + 1`. `points` and `tangts`
are of type `Vector{SVector{3,Float64}}`, and `epochs` is a `Vector{Float64}`.

| `method` | Limb points |
|---|---|
| `"TANGENT/ELLIPSOID"` | points of the reference ellipsoid where a ray from the observer is tangent; one per half plane |
| `"TANGENT/DSK/UNPRIORITIZED"` | points of the DSK surface where a ray from the observer is tangent; none, one or several per half plane |
| `"GUIDED/ELLIPSOID"` | as `"TANGENT/ELLIPSOID"` |
| `"GUIDED/DSK/UNPRIORITIZED"` | the outermost DSK surface point on the ray from the center of the ellipsoid limb through each ellipsoid limb point; at most one per half plane |

A DSK method can end with `"/SURFACES = <list>"`.

`corloc` is the place for which the light time and the stellar aberration are
computed: `"CENTER"` (the target center; all points have the same epoch) or
`"ELLIPSOID LIMB"` (each limb point of the reference ellipsoid; `TANGENT` methods only).

`schstp` and `soltol` (radians) are used by `"TANGENT/DSK"` only. The ray is rotated in
the half plane in steps of `schstp`; a tangent point is missed if the ray enters and
leaves the surface within one step. Each tangent ray is then located to within the
angle `soltol`. With an ellipsoid method both values are ignored.

`maxn` is the largest total number of points. `abcorr` is a correction string of
[`spkezr`](@ref) or an [`Aberration`](@ref). `et` cannot be a `ForwardDiff.Dual`. A
DSK method uses the tolerances that [`dskstl`](@ref) set on the context.

Throws `SpiceError` with the kind
- `:invalidoption` if `abcorr` is not a valid correction;
- `:idcodenotfound` if a body name is not known;
- `:bodiesnotdistinct` if `target` and `obsrvr` are the same body;
- `:noframe` if `fixref` is not known, `:invalidframe` if it is not centered on `target`;
- `:badmethodsyntax` if `method` is malformed, `:badpriorityspec` if a DSK method does
  not have `UNPRIORITIZED`, `:invalidlimbtype` if it has neither `TANGENT` nor `GUIDED`,
  `:invalidmethod` if it holds a sub-point or a terminator type;
- `:invalidlocus` if `corloc` is not one of the two values, `:badlimblocusmix` for
  `GUIDED` with `"ELLIPSOID LIMB"`;
- `:zerovector` if `refvec` is zero, `:degeneratecase` if it is parallel to the line
  from the target to the observer;
- `:invalidsize` if `maxn < 1`, `:invalidcount` unless `1 <= ncuts <= maxn`,
  `:invalidrollstep` if `rolstp` is zero and `ncuts > 1`;
- with a DSK method, `:invalidsearchstep` if `schstp` is not positive and
  `:invalidtolerance` if `soltol` is not positive;
- `:invalidconststep` if `schstp` is too small to change the angle;
- `:outofroom` if more than `maxn` points are found, `:arraytoosmall` if one half
  plane has more than 2000;
- `:noseparation` if the observer is at the target center;
- `:kernelvarnotfound` if the radii of `target` are not loaded, `:badradiuscount` if
  there are not 3 of them, `:degeneratecase` if one is not positive;
- `:dskdatanotfound` if no DSK data are loaded for `target`.

Also throws the errors of [`spkpos`](@ref).

# Examples
```julia
# Needs an SPK, a PCK with the radii and the rotation of Mars, and a leap seconds kernel.
et = str2et("2025 MAR 01 12:00:00 UTC")
z = [0.0, 0.0, 1.0]                     # north pole of Mars in IAU_MARS
npts, points, epochs, tangts = limbpt("TANGENT/ELLIPSOID", "MARS", et, "IAU_MARS", "LT+S",
    "CENTER", "EARTH", z, twopi() / 36, 36, 0.0, 0.0, 36)
length(points) == sum(npts) == 36       # one point per half plane on an ellipsoid
lonlat = [reclat(p)[2:3] for p in points]
```

See also [`termpt`](@ref), [`tangpt`](@ref), [`edlimb`](@ref).
"""
limbpt(ks::KernelSet, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbcorrArg, corloc::AbstractString, obsrvr::NameOrID, refvec::AbstractVector,
    rolstp::Real, ncuts::Integer, schstp::Real, soltol::Real, maxn::Integer) =
    _geo_limbpt(ks, _DSK_DEFAULT_TOLS, method, target, et, fixref, abcorr, corloc, obsrvr,
        refvec, rolstp, ncuts, schstp, soltol, maxn)
limbpt(ctx::SpiceContext, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbcorrArg, corloc::AbstractString, obsrvr::NameOrID, refvec::AbstractVector,
    rolstp::Real, ncuts::Integer, schstp::Real, soltol::Real, maxn::Integer) =
    _geo_limbpt(kernelset(ctx), dsktolerances(ctx), method, target, et, fixref, abcorr, corloc,
        obsrvr, refvec, rolstp, ncuts, schstp, soltol, maxn)
limbpt(method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID, abcorr::AbcorrArg,
    corloc::AbstractString, obsrvr::NameOrID, refvec::AbstractVector, rolstp::Real,
    ncuts::Integer, schstp::Real, soltol::Real, maxn::Integer) =
    limbpt(default_context(), method, target, et, fixref, abcorr, corloc, obsrvr, refvec,
        rolstp, ncuts, schstp, soltol, maxn)

"""
    termpt([ctx,] method, ilusrc, target, et, fixref, abcorr, corloc, obsrvr, refvec,
           rolstp, ncuts, schstp, soltol, maxn) -> (npts, points, epochs, trmvcs)

Return terminator points of `target` lit by the source `ilusrc`, as seen by `obsrvr` at
`et` (seconds past J2000 TDB), as `termpt_c` of CSPICE.

The points lie in `ncuts` half planes bounded by the line from the target center to
the source center. Half plane `i` contains `refvec` (body-fixed frame `fixref`)
rotated about that line by `(i - 1) * rolstp` radians. The outputs have the layout of
[`limbpt`](@ref); `trmvcs` holds the vectors from the observer to the points (km).

The source is a sphere with the largest radius of `ilusrc`. `method` joins the shadow,
the point definition and the shape, for example `"UMBRAL/TANGENT/ELLIPSOID"` or
`"PENUMBRAL/TANGENT/DSK/UNPRIORITIZED"`:

- `UMBRAL`: the boundary of the region from which no part of the source is visible.
  `PENUMBRAL`: the boundary of the region from which the whole source is visible.
- `TANGENT`: points where a ray tangent to the source is tangent to the surface.
  `GUIDED`: on a DSK surface, the outermost point on the ray from the target center
  through each terminator point of the reference ellipsoid.
- `ELLIPSOID`, or `DSK/UNPRIORITIZED` with an optional `"/SURFACES = <list>"`.

`corloc` is `"CENTER"` or `"ELLIPSOID TERMINATOR"` (`TANGENT` methods only), with the
meaning it has in [`limbpt`](@ref). `abcorr` is `"NONE"` or a reception correction.
`schstp`, `soltol`, `maxn`, `et` and `ctx` are as in [`limbpt`](@ref).

As in CSPICE, a half plane with `n` points takes `n^2` of the `maxn` places, not `n`.
With an ellipsoid `n` is 1. With a DSK method, `:outofroom` can occur before `maxn`
points are found, so give `maxn` a margin.

Throws the kinds of [`limbpt`](@ref) for the bodies, the frame, `refvec`, `rolstp` and
`maxn`, and `SpiceError` with the kind
- `:invalidoption` if `abcorr` is not valid or is a transmission correction;
- `:bodiesnotdistinct` if `target` is also the observer or the source;
- `:invalidtermtype` if `method` has neither `TANGENT` nor `GUIDED`, `:invalidshadow`
  if it has neither `UMBRAL` nor `PENUMBRAL`, `:invalidmethod` if it holds a
  sub-point type;
- `:invalidlocus` if `corloc` is not one of the two values, `:badtermlocusmix` for
  `GUIDED` with `"ELLIPSOID TERMINATOR"`;
- `:invalidcount` unless `0 <= ncuts <= maxn`;
- `:objectstooclose` if the source touches the bounding sphere of the target, the sphere
  about its center that contains it;
- in a `TANGENT/DSK` search, `:invalidconststep` if `schstp` is not positive and
  `:invalidtolerance` if `soltol` is not positive;
- `:noconvergence` if the search for one terminator point of the ellipsoid takes more
  than 20 iterations;
- `:kernelvarnotfound` if the radii of `target` or `ilusrc` are not loaded,
  `:badradiuscount` if there are not 3 of them, `:invalidaxislength` if one is not
  positive;
- `:dskdatanotfound` if no DSK data are loaded for `target`.

# Examples
```julia
# Needs an SPK, a PCK with the radii of Mars and the Sun and the rotation of Mars, and
# a leap seconds kernel.
et = str2et("2025 MAR 01 12:00:00 UTC")
z = [0.0, 0.0, 1.0]
npts, points, epochs, trmvcs = termpt("UMBRAL/TANGENT/ELLIPSOID", "SUN", "MARS", et,
    "IAU_MARS", "LT+S", "CENTER", "EARTH", z, twopi() / 36, 36, 0.0, 0.0, 36)
lonlat = [reclat(p)[2:3] for p in points]
```

See also [`limbpt`](@ref), [`edterm`](@ref), [`ilumin`](@ref).
"""
termpt(ks::KernelSet, method::AbstractString, ilusrc::NameOrID, target::NameOrID, et::Real,
    fixref::NameOrID, abcorr::AbcorrArg, corloc::AbstractString, obsrvr::NameOrID,
    refvec::AbstractVector, rolstp::Real, ncuts::Integer, schstp::Real, soltol::Real,
    maxn::Integer) = _geo_termpt(ks, _DSK_DEFAULT_TOLS, method, ilusrc, target, et, fixref,
    abcorr, corloc, obsrvr, refvec, rolstp, ncuts, schstp, soltol, maxn)
termpt(ctx::SpiceContext, method::AbstractString, ilusrc::NameOrID, target::NameOrID, et::Real,
    fixref::NameOrID, abcorr::AbcorrArg, corloc::AbstractString, obsrvr::NameOrID,
    refvec::AbstractVector, rolstp::Real, ncuts::Integer, schstp::Real, soltol::Real,
    maxn::Integer) = _geo_termpt(kernelset(ctx), dsktolerances(ctx), method, ilusrc, target,
    et, fixref, abcorr, corloc, obsrvr, refvec, rolstp, ncuts, schstp, soltol, maxn)
termpt(method::AbstractString, ilusrc::NameOrID, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbcorrArg, corloc::AbstractString, obsrvr::NameOrID, refvec::AbstractVector,
    rolstp::Real, ncuts::Integer, schstp::Real, soltol::Real, maxn::Integer) =
    termpt(default_context(), method, ilusrc, target, et, fixref, abcorr, corloc, obsrvr,
        refvec, rolstp, ncuts, schstp, soltol, maxn)

# --- tangent point of a ray -------------------------------------------------------------

"""
    tangpt([ctx,] method, target, et, fixref, abcorr, corloc, obsrvr, dref, dvec)
        -> (tanpt, alt, range, srfpt, trgepc, srfvec)

Return the tangent point of the ray from `obsrvr` along `dvec` on the reference
ellipsoid of `target`, as `tangpt_c` of CSPICE.

- `tanpt`: the point of the ray nearest to the ellipsoid, in km in the body-fixed
  frame `fixref`.
- `alt`: the altitude of `tanpt` above the ellipsoid (km).
- `range`: the distance from the observer to `tanpt` (km).
- `srfpt`: the point of the ellipsoid nearest to `tanpt`, in km in `fixref`.
- `trgepc`: the epoch at the point named by `corloc` (seconds past J2000 TDB). Without
  a correction it is `et`.
- `srfvec`: the vector from the observer to `srfpt`, in km in `fixref` at `trgepc`.

If the ray meets the ellipsoid, `alt` is 0 and `tanpt` equals `srfpt`. If
the ray points away from the target, `tanpt` is the observer position and `range` is 0.

`dvec` is the ray direction (3 elements, not zero) in the frame `dref`. `method` must
be `"ELLIPSOID"`. `corloc` is `"TANGENT POINT"` or `"SURFACE POINT"`: the point for
which the light time and the stellar aberration are computed. `abcorr` is a correction
string of [`spkezr`](@ref) or an [`Aberration`](@ref). `et` is in seconds past J2000
TDB and cannot be a `ForwardDiff.Dual`. The call does not allocate when the bodies and
the frames are integer IDs.

Throws `SpiceError` with the kind
- `:invalidoption` if `abcorr` is not a valid correction;
- `:notsupported` if `method` is not `"ELLIPSOID"` or `corloc` is not one of the two
  values;
- `:zerovector` if `dvec` is zero, `:badarg` if it does not have 3 elements;
- `:idcodenotfound` if a body name is not known;
- `:bodiesnotdistinct` if `target` and `obsrvr` are the same body;
- `:kernelvarnotfound` if the radii of `target` are not loaded, `:invalidcount` if
  there are not 3 of them, `:invalidradii` if one is not positive;
- `:noframe` if `fixref` or `dref` is not known, `:invalidframe` if `fixref` is not
  centered on `target`;
- `:noseparation` if the observer is at the target center, `:invalidgeometry` if it is
  inside or on the ellipsoid.

Also throws the errors of [`spkpos`](@ref).

See also [`limbpt`](@ref), [`sincpt`](@ref), [`npedln`](@ref).
"""
function tangpt(ks::KernelSet, method::AbstractString, target::NameOrID, et0::Real,
        fixref::NameOrID, abcorr::AbcorrArg, corloc::AbstractString, obsrvr::NameOrID,
        dref::NameOrID, dvec0::AbstractVector)
    corr = _geo_abcorr(abcorr)
    s = _geo_sign(corr)
    (method == "ELLIPSOID" || _dsk_eqstr(method, "ELLIPSOID")) || throw(SpiceError(
        :notsupported, "the method <$method> is not supported; use ELLIPSOID"))
    tanloc = _geo_tanloc(corloc)
    dvec = _geo_pos(dvec0)
    _m_iszero(dvec) && throw(SpiceError(:zerovector, "the ray direction is the zero vector"))
    trgcde = _geo_body(ks, target, "target")
    obscde = _geo_body(ks, obsrvr, "observer")
    obscde == trgcde && _geo_same_bodies("tangent point", obsrvr)
    a, b, c = _geo_radii3(ks, trgcde, :invalidcount)
    fixfid = _geo_fixfid(ks, fixref, trgcde, target)
    dfrcde = _geo_namfrm(ks, dref)
    dcentr, dclass, found = _geo_frinfo(ks, dfrcde)
    found || _geo_no_frame(dref)
    et = Float64(check_epoch(et0))
    tpos, lt = _geo_spkezp(ks, trgcde, et, fixfid, corr, obscde)
    _m_iszero(tpos) && throw(SpiceError(:noseparation,
        "the distance between the observer and the target is zero"))
    fixobs = -tpos
    _m_norm(edpnt(fixobs, a, b, c)) >= _m_norm(fixobs) && throw(SpiceError(:invalidgeometry,
        "the observer is inside or on the surface of the reference ellipsoid of the target"))
    trgepc = et + s * lt
    refepc = et
    if !(dclass == Int32(1) || dcentr == obscde || !corr.uselt)
        _, ltcent = _geo_spkezp(ks, dcentr, et, Int32(1), corr, obscde)
        refepc = et + s * ltcent
    end
    j2dir = dfrcde == Int32(1) ? dvec : _m_mxv(frame_rotation(ks, dfrcde, Int32(1), refepc), dvec)
    fixdir = _m_mxv(frame_rotation(ks, Int32(1), fixfid, trgepc), j2dir)
    srfpt, alt = npedln(a, b, c, fixobs, fixdir)
    srfvec = srfpt - fixobs
    tanpt, _ = nplnpt(fixobs, fixdir, srfpt)
    if !corr.uselt
        if _m_dot(tanpt - fixobs, fixdir) <= 0.0
            # The ray points away from the target: the tangent point is the observer.
            srfpt, alt = _m_nearpt(fixobs, a, b, c)
            return fixobs, alt, 0.0, srfpt, et, srfpt - fixobs
        end
        alt == 0.0 && return srfpt, alt, _m_norm(srfvec), srfpt, et, srfvec
        return tanpt, alt, _m_norm(fixobs - tanpt), srfpt, et, srfvec
    end
    nitr = corr.converged ? 10 : 3
    dist = tanloc ? _m_norm(fixobs - tanpt) : _m_norm(srfvec)
    robs, vobs = _ssb(ks, obscde, et, Int32(1))
    vstl = corr.transmit ? -vobs : vobs
    lt = dist / CLIGHT
    trgepc = et + s * lt
    prevlt = 0.0
    prvepc = trgepc
    i = 0
    stlcnv = !corr.stellar
    ltcnv = false
    stlobs = fixobs
    # invariant: `tanpt`, `srfpt` and `lt` are the values after i iterations
    while i < nitr && !(ltcnv && stlcnv)
        prvtan = tanpt
        prvsrf = srfpt
        j2opos = robs - _ssbpos(ks, trgcde, trgepc, Int32(1))
        j2fixm = frame_rotation(ks, Int32(1), fixfid, trgepc)
        stlfix = zero(Vec3{Float64})
        if corr.stellar
            # The stellar aberration offset of the locus point, applied to the observer.
            j2lpos = _m_mtxv(j2fixm, tanloc ? tanpt : srfpt) - j2opos
            stlfix = _m_mxv(j2fixm, _stelab(j2lpos, vstl) - j2lpos)
        end
        fixobs = _m_mxv(j2fixm, j2opos)
        fixdir = _m_mxv(j2fixm, j2dir)
        stlobs = fixobs - stlfix
        srfpt, alt = npedln(a, b, c, stlobs, fixdir)
        tanpt = alt != 0.0 ? nplnpt(stlobs, fixdir, srfpt)[1] : srfpt
        srfvec = srfpt - stlobs
        if _m_dot(tanpt - stlobs, fixdir) <= 0.0
            tanpt = stlobs
            srfpt, alt = _m_nearpt(stlobs, a, b, c)
            srfvec = srfpt - tanpt
            lt = tanloc ? 0.0 : alt / CLIGHT
            trgepc = tanloc ? et : et + s * lt
        else
            dist = tanloc ? _m_norm(fixobs - tanpt) : _m_norm(fixobs - srfpt)
            lt = dist / CLIGHT
            trgepc = et + s * lt
        end
        ltdiff = abs(lt - prevlt)
        epcdif = abs(trgepc - prvepc)
        prevlt = lt
        prvepc = trgepc
        i += 1
        ltcnv = ltdiff < abs(lt) * 1e-16 || epcdif == 0.0
        corr.stellar && (stlcnv = vrel(tanpt, prvtan) < 1e-16 && vrel(srfpt, prvsrf) < 1e-16)
    end
    return tanpt, alt, _m_norm(stlobs - tanpt), srfpt, trgepc, srfvec
end
tangpt(ctx::SpiceContext, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbcorrArg, corloc::AbstractString, obsrvr::NameOrID, dref::NameOrID,
    dvec::AbstractVector) =
    tangpt(kernelset(ctx), method, target, et, fixref, abcorr, corloc, obsrvr, dref, dvec)
tangpt(method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID, abcorr::AbcorrArg,
    corloc::AbstractString, obsrvr::NameOrID, dref::NameOrID, dvec::AbstractVector) =
    tangpt(default_context(), method, target, et, fixref, abcorr, corloc, obsrvr, dref, dvec)

@inline function _geo_tanloc(corloc::AbstractString)
    corloc == "TANGENT POINT" && return true
    corloc == "SURFACE POINT" && return false
    return _geo_tanloc_slow(corloc)
end
@noinline function _geo_tanloc_slow(corloc::AbstractString)
    _dsk_eqstr(corloc, "TANGENTPOINT") && return true
    _dsk_eqstr(corloc, "SURFACEPOINT") && return false
    throw(SpiceError(:notsupported, "the aberration correction locus <$corloc> is not " *
        "recognized; the values are TANGENT POINT and SURFACE POINT"))
end
