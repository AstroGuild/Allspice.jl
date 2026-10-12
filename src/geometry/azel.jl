# Azimuth and elevation: `azlrec`, `recazl`, `dazldr`, `drdazl` and `azlcpo` of CSPICE
# N0067. The coordinate calls are generic over the number type and do not allocate.

"""
    azlrec(range, az, el, azccw, elplsz) -> SVector{3}
    azlrec(range, az, el; azccw=true, elplsz=true)

Return the rectangular coordinates of the point at the distance `range`, the azimuth
`az` and the elevation `el` (radians), as `azlrec_c`. `azccw` tells whether the azimuth
grows counterclockwise about the +Z axis (seen from +Z), and `elplsz` whether the
elevation grows toward +Z.

Contract: no check of the ranges, as in CSPICE. Postcondition: the norm of the result
is `|range|`.
"""
function azlrec(range::Real, az::Real, el::Real, azccw::Bool, elplsz::Bool)
    r, a, e = promote(float(range), float(az), float(el))
    # As in CSPICE, a zero angle keeps its sign.
    return _m_latrec(r, (azccw || primal(a) == 0.0) ? a : -a, (elplsz || primal(e) == 0.0) ? e : -e)
end
azlrec(range::Real, az::Real, el::Real; azccw::Bool=true, elplsz::Bool=true) =
    azlrec(range, az, el, azccw, elplsz)

"""
    recazl(rectan, azccw, elplsz) -> (range, az, el)
    recazl(rectan; azccw=true, elplsz=true)

Return the distance, the azimuth in [0, 2 pi) and the elevation in [-pi/2, pi/2] of the
point `rectan`, as `recazl_c`. The flags are those of [`azlrec`](@ref).

Contract: `rectan` has 3 elements. Postcondition: `azlrec` of the result gives
`rectan` to rounding.
"""
function recazl(rectan::AbstractVector, azccw::Bool, elplsz::Bool)
    range, az, el = recrad(rectan)
    if !azccw && primal(az) > 0.0
        az = 2.0 * Float64(pi) - az
        primal(az) < 0.0 && (az = zero(az))
    end
    (elplsz || primal(el) == 0.0) || (el = -el)
    return range, az, el
end
recazl(rectan::AbstractVector; azccw::Bool=true, elplsz::Bool=true) =
    recazl(rectan, azccw, elplsz)

"""
    dazldr(x, y, z, azccw, elplsz) -> SMatrix{3,3}
    dazldr(x, y, z; azccw=true, elplsz=true)

Return the Jacobian of the transformation from rectangular coordinates to range,
azimuth and elevation at the point `(x, y, z)`, as `dazldr_c`: row `i` holds the
partial derivatives of coordinate `i` of (range, azimuth, elevation).

Contract: the point is not on the Z axis (`SpiceError(:pointonzaxis)`).
"""
function dazldr(x::Real, y::Real, z::Real, azccw::Bool, elplsz::Bool)
    (primal(x) == 0.0 && primal(y) == 0.0) && throw(SpiceError(:pointonzaxis,
        "the Jacobian is not defined for a point on the Z axis"))
    j = dlatdr(x, y, z)
    sa = azccw ? 1.0 : -1.0
    se = elplsz ? 1.0 : -1.0
    return SMatrix{3,3}(j[1, 1], sa * j[2, 1], se * j[3, 1], j[1, 2], sa * j[2, 2],
        se * j[3, 2], j[1, 3], sa * j[2, 3], se * j[3, 3])
end
dazldr(x::Real, y::Real, z::Real; azccw::Bool=true, elplsz::Bool=true) =
    dazldr(x, y, z, azccw, elplsz)

"""
    drdazl(range, az, el, azccw, elplsz) -> SMatrix{3,3}
    drdazl(range, az, el; azccw=true, elplsz=true)

Return the Jacobian of the transformation from range, azimuth and elevation to
rectangular coordinates, as `drdazl_c`: column `j` holds the partial derivatives with
respect to coordinate `j` of (range, azimuth, elevation).

Contract: `range >= 0` (`SpiceError(:valueoutofrange)`).
"""
function drdazl(range::Real, az::Real, el::Real, azccw::Bool, elplsz::Bool)
    primal(range) < 0.0 && throw(SpiceError(:valueoutofrange,
        "the range must not be negative; it is $(primal(range))"))
    sa = azccw ? 1.0 : -1.0
    se = elplsz ? 1.0 : -1.0
    j = drdlat(range, sa * az, se * el)
    return SMatrix{3,3}(j[1, 1], j[2, 1], j[3, 1], sa * j[1, 2], sa * j[2, 2], sa * j[3, 2],
        se * j[1, 3], se * j[2, 3], se * j[3, 3])
end
drdazl(range::Real, az::Real, el::Real; azccw::Bool=true, elplsz::Bool=true) =
    drdazl(range, az, el, azccw, elplsz)

"""
    azlcpo(method, target, et, abcorr, azccw, elplsz, obspos, obsctr, obsref) -> (azlsta, lt)

Return the range, azimuth and elevation of `target` and their rates (km, rad, km/s,
rad/s) as seen by an observer with the constant position `obspos` (km, relative to the
body `obsctr`, in the body-fixed frame `obsref`), and the one-way light time, as
`azlcpo_c`. The angles are in the topocentric frame at the observer: +Z along the
outward normal of the reference ellipsoid of `obsctr` at the point below the observer,
+X toward the north. `method` must be `"ELLIPSOID"`. Forms: on a `KernelSet`, on a
`SpiceContext`, and on the default context. The arguments are in the order of CSPICE.

Contract: `obspos` has 3 elements and `obsref` is centered on `obsctr`. The epoch can
be a dual number. No allocation for integer IDs. Throws `SpiceError`:
`:idcodenotfound`, `:unknownframe`, `:invalidframe`, `:invalidmethod`,
`:invalidoption`, `:pointonzaxis` for a target on the local vertical, and the errors of
[`spkcpo`](@ref) and of the radii lookup.
"""
function azlcpo(ks::KernelSet, method::AbstractString, target::NameOrID, et::Real,
        abcorr::AbcorrArg, azccw::Bool, elplsz::Bool, obspos::AbstractVector, obsctr::NameOrID,
        obsref::NameOrID)
    obscde = _geo_body(ks, obsctr, "observer center")
    fxfcde = _geo_namfrm(ks, obsref)
    center, _, found = _geo_frinfo(ks, fxfcde)
    found || _geo_unknown_frame(obsref)
    center == obscde || _geo_bad_center(obsref, obsctr, center)
    (method == "ELLIPSOID" || _dsk_eqstr(method, "ELLIPSOID")) || throw(SpiceError(
        :invalidmethod, "the method <$method> is not supported; use ELLIPSOID"))
    pos = _geo_pos(obspos)
    xftopo = if pos[1] == 0.0 && pos[2] == 0.0
        # On the polar axis the local +X axis is the +X axis of the body-fixed frame.
        pos[3] < 0.0 ? _m_rotmat(one(Mat3{Float64}), Float64(pi), 1) : one(Mat3{Float64})
    else
        a, b, c = _geo_radii(ks, obscde)
        obsspt, _ = _m_nearpt(pos, a, b, c)
        twovec(_m_surfnm(a, b, c, obsspt), 3, SVector(0.0, 0.0, 1.0), 1)
    end
    corr = _geo_abcorr(abcorr)
    trgcde = _geo_body(ks, target, "target")
    f = _GeoFixed(pos, zero(Vec3{Float64}), 0.0, obscde, fxfcde)
    state, lt = _geo_cvo(ks, trgcde, check_epoch(et), fxfcde, _GEO_AT_OBSERVER, corr, f)
    r, v = _geo_pv(state)
    lr = xftopo * r
    lv = xftopo * v
    range, az, el = recazl(lr, azccw, elplsz)
    rates = dazldr(lr[1], lr[2], lr[3], azccw, elplsz) * lv
    return SVector(range, az, el, rates[1], rates[2], rates[3]), lt
end
azlcpo(ctx::SpiceContext, method::AbstractString, target::NameOrID, et::Real,
    abcorr::AbcorrArg, azccw::Bool, elplsz::Bool, obspos::AbstractVector, obsctr::NameOrID,
    obsref::NameOrID) =
    azlcpo(kernelset(ctx), method, target, et, abcorr, azccw, elplsz, obspos, obsctr, obsref)
azlcpo(method::AbstractString, target::NameOrID, et::Real, abcorr::AbcorrArg, azccw::Bool,
    elplsz::Bool, obspos::AbstractVector, obsctr::NameOrID, obsref::NameOrID) =
    azlcpo(default_context(), method, target, et, abcorr, azccw, elplsz, obspos, obsctr, obsref)
