# Coordinate conversions of CSPICE, their Jacobians, state conversion between coordinate
# systems, and unit conversion.
#
# Angles are in radians. The tuple orders are the argument orders of CSPICE:
# latitudinal (radius, lon, lat), right ascension (range, ra, dec), spherical
# (r, colat, slon), cylindrical (r, clon, z), geodetic and planetographic
# (lon, lat, alt). A Jacobian `J` has `J[i, j] = d out_i / d in_j`.

const _M_TWOPI = 2.0 * Float64(pi)
const _M_HALFPI = Float64(pi) / 2

# --- rectangular, latitudinal, right ascension ----------------------------------------

"""
`reclat` kernel. With `scaled` the longitude comes from the components divided by the
largest one (`reclat_c`); without it from the components themselves (`reclat_`).
"""
@inline function _m_reclat(v::SVector{3}, scaled::Bool)
    vmax = _m_maxabs(v)
    z0 = zero(atan(v[1] / 1.0))
    vmax > 0.0 || return z0, z0, z0
    x1 = v[1] / vmax
    y1 = v[2] / vmax
    z1 = v[3] / vmax
    radius = vmax * sqrt(x1 * x1 + y1 * y1 + z1 * z1)
    lat = atan(z1, sqrt(x1 * x1 + y1 * y1))
    (primal(x1) == 0.0 && primal(y1) == 0.0) && return radius + z0, z0, lat
    return radius + z0, (scaled ? atan(y1, x1) : atan(v[2], v[1])), lat
end

"""
    reclat(rectan) -> (radius, lon, lat)

Return the latitudinal coordinates of a point, as `reclat_c`: `lon` in (-pi, pi], `lat`
in [-pi/2, pi/2]. Postcondition: `lon` is 0 on the z axis, and all three are 0 at the
origin; `latrec(reclat(x)...)` is `x` to rounding.
"""
reclat(rectan::AbstractVector) = _m_reclat(_m_v3(rectan), true)

@inline function _m_latrec(radius, lon, lat)
    slon, clon = sincos(lon)
    slat, clat = sincos(lat)
    return SVector(radius * clon * clat, radius * slon * clat, radius * slat)
end

"""
    latrec(radius, lon, lat) -> SVector{3}

Return the rectangular coordinates of a point from its latitudinal coordinates, as
`latrec_c`.
"""
latrec(radius::Real, lon::Real, lat::Real) = _m_latrec(_m_f(radius), _m_f(lon), _m_f(lat))

"""
    recrad(rectan) -> (range, ra, dec)

Return the range, right ascension in [0, 2 pi) and declination of a point, as
`recrad_c`.
"""
function recrad(rectan::AbstractVector)
    range, ra, dec = _m_reclat(_m_v3(rectan), true)
    return range, (primal(ra) < 0.0 ? ra + _M_TWOPI : ra), dec
end

"""
    radrec(range, ra, dec) -> SVector{3}

Return the rectangular coordinates of a point from its range, right ascension and
declination, as `radrec_c`.
"""
radrec(range::Real, ra::Real, dec::Real) = _m_latrec(_m_f(range), _m_f(ra), _m_f(dec))

# --- spherical ------------------------------------------------------------------------

@inline function _m_recsph(v::SVector{3})
    big = _m_maxabs(v)
    z0 = zero(atan(v[1] / 1.0))
    big > 0.0 || return z0, z0, z0
    x = v[1] / big
    y = v[2] / big
    z = v[3] / big
    r = big * sqrt(x * x + y * y + z * z)
    colat = atan(sqrt(x * x + y * y), z)
    (primal(v[1]) == 0.0 && primal(v[2]) == 0.0) && return r + z0, colat, z0
    return r + z0, colat, atan(v[2], v[1])
end

"""
    recsph(rectan) -> (r, colat, slon)

Return the spherical coordinates of a point, as `recsph_c`: the colatitude in [0, pi]
is the angle from the z axis and `slon` is in (-pi, pi]. Postcondition: `slon` is 0 on
the z axis and all three are 0 at the origin.
"""
recsph(rectan::AbstractVector) = _m_recsph(_m_v3(rectan))

@inline function _m_sphrec(r, colat, slon)
    sl, cl = sincos(slon)
    sc, cc = sincos(colat)
    return SVector(r * cl * sc, r * sl * sc, r * cc)
end

"""
    sphrec(r, colat, slon) -> SVector{3}

Return the rectangular coordinates of a point from its spherical coordinates, as
`sphrec_c`.
"""
sphrec(r::Real, colat::Real, slon::Real) = _m_sphrec(_m_f(r), _m_f(colat), _m_f(slon))

# --- cylindrical ----------------------------------------------------------------------

@inline function _m_reccyl(v::SVector{3})
    big = max(abs(primal(v[1])), abs(primal(v[2])))
    z0 = zero(atan(v[1] / 1.0))
    big == 0.0 && return z0, z0, v[3] + z0
    x = v[1] / big
    y = v[2] / big
    r = big * sqrt(x * x + y * y)
    clon = atan(y, x)
    return r + z0, (primal(clon) < 0.0 ? clon + _M_TWOPI : clon), v[3] + z0
end

"""
    reccyl(rectan) -> (r, clon, z)

Return the cylindrical coordinates of a point, as `reccyl_c`: `clon` in [0, 2 pi).
Postcondition: `clon` is 0 on the z axis.
"""
reccyl(rectan::AbstractVector) = _m_reccyl(_m_v3(rectan))

@inline function _m_cylrec(r, clon, z)
    s, c = sincos(clon)
    return SVector(r * c, r * s, z + zero(r * c))
end

"""
    cylrec(r, clon, z) -> SVector{3}

Return the rectangular coordinates of a point from its cylindrical coordinates, as
`cylrec_c`.
"""
cylrec(r::Real, clon::Real, z::Real) = _m_cylrec(_m_f(r), _m_f(clon), _m_f(z))

"Radius of the point (r, z) with the scaling of `cyllat_c`; zero at the origin."
@inline function _m_hypot(r, z)
    big = max(abs(primal(r)), abs(primal(z)))
    big > 0.0 || return zero(sqrt(one(r) * one(z)))
    x = r / big
    y = z / big
    return big * sqrt(x * x + y * y)
end

"""
    cyllat(r, clon, z) -> (radius, lon, lat)

Convert cylindrical coordinates to latitudinal coordinates, as `cyllat_c`. The
longitude is returned as given.
"""
function cyllat(r::Real, clon::Real, z::Real)
    r, clon, z = promote(_m_f(r), _m_f(clon), _m_f(z))
    rho = _m_hypot(r, z)
    lat = primal(rho) == 0.0 ? zero(atan(z, r)) : atan(z, r)
    return rho + zero(lat), clon + zero(lat), lat
end

"""
    latcyl(radius, lon, lat) -> (r, clon, z)

Convert latitudinal coordinates to cylindrical coordinates, as `latcyl_c`.
"""
function latcyl(radius::Real, lon::Real, lat::Real)
    radius, lon, lat = promote(_m_f(radius), _m_f(lon), _m_f(lat))
    s, c = sincos(lat)
    return radius * c, lon + zero(c), radius * s
end

"""
    cylsph(r, clon, z) -> (radius, colat, slon)

Convert cylindrical coordinates to spherical coordinates, as `cylsph_c`.
"""
function cylsph(r::Real, clon::Real, z::Real)
    r, clon, z = promote(_m_f(r), _m_f(clon), _m_f(z))
    rh = _m_hypot(r, z)
    th = primal(rh) == 0.0 ? zero(atan(r, z)) : atan(r, z)
    return rh + zero(th), th, clon + zero(th)
end

"""
    sphcyl(radius, colat, slon) -> (r, clon, z)

Convert spherical coordinates to cylindrical coordinates, as `sphcyl_c`.
"""
function sphcyl(radius::Real, colat::Real, slon::Real)
    radius, colat, slon = promote(_m_f(radius), _m_f(colat), _m_f(slon))
    s, c = sincos(colat)
    return radius * s, slon + zero(c), radius * c
end

"""
    sphlat(r, colat, slon) -> (radius, lon, lat)

Convert spherical coordinates to latitudinal coordinates, as `sphlat_c`: the latitude
is `pi/2 - colat`.
"""
function sphlat(r::Real, colat::Real, slon::Real)
    r, colat, slon = promote(_m_f(r), _m_f(colat), _m_f(slon))
    return r, slon, _M_HALFPI - colat
end

"""
    latsph(radius, lon, lat) -> (rho, colat, slon)

Convert latitudinal coordinates to spherical coordinates, as `latsph_c`: the colatitude
is `pi/2 - lat`.
"""
function latsph(radius::Real, lon::Real, lat::Real)
    radius, lon, lat = promote(_m_f(radius), _m_f(lon), _m_f(lat))
    return radius, _M_HALFPI - lat, lon
end

# --- geodetic -------------------------------------------------------------------------

@noinline _m_badre(re) = throw(SpiceError(:valueoutofrange,
    "the equatorial radius must be positive; it is $(primal(re))"))
@noinline _m_badflat(f) = throw(SpiceError(:valueoutofrange,
    "the flattening coefficient must be less than 1; it is $(primal(f))"))

@inline function _m_recgeo(v::SVector{3}, re, f)
    primal(re) > 0.0 || _m_badre(re)
    primal(f) < 1.0 || _m_badflat(f)
    c = re - f * re
    base, alt = _m_nearpt(v, re, re, c)
    normal = _m_surfnm(re, re, c, base)
    _, _, lat = _m_reclat(normal, false)
    lon = (primal(v[1]) == 0.0 && primal(v[2]) == 0.0) ? zero(lat) : atan(v[2], v[1]) + zero(lat)
    return lon, lat, alt + zero(lat)
end

"""
    recgeo(rectan, re, f) -> (lon, lat, alt)

Return the geodetic longitude, latitude and altitude of a point relative to the
spheroid with the equatorial radius `re` and the flattening `f` (polar radius
`re (1 - f)`), as `recgeo_c`. The latitude is the one of the outward normal at the
nearest surface point, which `nearpt` gives. Requires `re > 0` and `f < 1`
(`SpiceError(:valueoutofrange)`).
"""
recgeo(rectan::AbstractVector, re::Real, f::Real) = _m_recgeo(_m_v3(rectan), _m_f(re), _m_f(f))

@inline function _m_georec(lon, lat, alt, re, f)
    primal(re) > 0.0 || _m_badre(re)
    primal(f) < 1.0 || _m_badflat(f)
    rp = re - f * re
    sphi, cphi = sincos(lat)
    slmbda, clmbda = sincos(lon)
    big = max(abs(primal(re * cphi)), abs(primal(rp * sphi)))
    x = re * cphi / big
    y = rp * sphi / big
    scale = 1.0 / (big * sqrt(x * x + y * y))
    base = SVector(scale * re * re * clmbda * cphi, scale * re * re * slmbda * cphi,
        scale * rp * rp * sphi)
    normal = _m_surfnm(re, re, rp, base)
    return 1.0 * base + alt * normal
end

"""
    georec(lon, lat, alt, re, f) -> SVector{3}

Return the rectangular coordinates of the point with the geodetic longitude, latitude
and altitude `lon`, `lat`, `alt` relative to the spheroid `re`, `f`, as `georec_c`.
Requires `re > 0` and `f < 1` (`SpiceError(:valueoutofrange)`).
"""
georec(lon::Real, lat::Real, alt::Real, re::Real, f::Real) =
    _m_georec(_m_f(lon), _m_f(lat), _m_f(alt), _m_f(re), _m_f(f))

# --- planetographic -------------------------------------------------------------------

"""
    _m_pgrsense(ks, body) -> Int

+1 if the planetographic longitude of `body` is positive east, -1 if positive west, as
`recpgr_`: the pool variable `BODY<id>_PGR_POSITIVE_LON` (`EAST` or `WEST`) decides if
it is present; otherwise the sign of the prime meridian rate `BODY<id>_PM[2]` does
(retrograde rotation is positive east), with east for the Earth, the Moon and the Sun.
"""
function _m_pgrsense(ks::KernelSet, body::NameOrID)
    id = _m_bodyid(ks, body, :idcodenotfound)
    kv = get(ks.pool, string("BODY", id, "_PGR_POSITIVE_LON"), nothing)
    if kv isa Vector{String} && !isempty(kv)
        s = uppercase(filter(!isspace, kv[1]))
        s = s[1:min(4, lastindex(s))]
        s == "EAST" && return 1
        s == "WEST" && return -1
        throw(SpiceError(:invalidoption,
            "the variable BODY$(id)_PGR_POSITIVE_LON may be EAST or WEST; it is $(kv[1])"))
    end
    id == 399 && return 1
    pm = get(ks.pool, string("BODY", id, "_PM"), nothing)
    (pm isa Vector{Float64} && length(pm) >= 2) || throw(SpiceError(:missingdata,
        "the prime meridian rate BODY$(id)_PM is required but not available"))
    (id == 301 || id == 10) && return 1
    return pm[2] < 0.0 ? 1 : -1
end

"""
    recpgr(body, rectan, re, f) -> (lon, lat, alt)

Return the planetographic longitude in [0, 2 pi], latitude and altitude of a point
relative to the spheroid `re`, `f` of `body` (name or ID), as `recpgr_c`. The sense of
the longitude comes from the kernel pool of the context (see `pgrrec`). Requires a
known body (`SpiceError(:idcodenotfound)`), `re > 0` and `f < 1`
(`SpiceError(:valueoutofrange)`).
"""
function recpgr(ks::KernelSet, body::NameOrID, rectan::AbstractVector, re::Real, f::Real)
    v = _m_v3(rectan)
    id = _m_bodyid(ks, body, :idcodenotfound)
    primal(re) > 0.0 || _m_badre(re)
    primal(f) < 1.0 || _m_badflat(f)
    sense = _m_pgrsense(ks, id)
    lon, lat, alt = _m_recgeo(v, _m_f(re), _m_f(f))
    lon = sense * lon
    primal(lon) < 0.0 && (lon += _M_TWOPI)
    primal(lon) < 0.0 && (lon = zero(lon))
    primal(lon) > _M_TWOPI && (lon = _M_TWOPI + zero(lon))
    return lon, lat, alt
end

"""
    pgrrec(body, lon, lat, alt, re, f) -> SVector{3}

Return the rectangular coordinates of the point with the planetographic coordinates
`lon`, `lat`, `alt` relative to the spheroid `re`, `f` of `body`, as `pgrrec_c`. The
longitude is positive east if the pool variable `BODY<id>_PGR_POSITIVE_LON` says
`EAST`, positive west if it says `WEST`; without the variable it is positive west for
prograde rotation (`BODY<id>_PM` rate positive) and positive east for retrograde
rotation and for the Earth, the Moon and the Sun. Throws
`SpiceError(:idcodenotfound)`, `SpiceError(:valueoutofrange)`,
`SpiceError(:invalidoption)` or `SpiceError(:missingdata)`.
"""
function pgrrec(ks::KernelSet, body::NameOrID, lon::Real, lat::Real, alt::Real, re::Real, f::Real)
    id = _m_bodyid(ks, body, :idcodenotfound)
    primal(re) > 0.0 || _m_badre(re)
    primal(f) < 1.0 || _m_badflat(f)
    sense = _m_pgrsense(ks, id)
    return _m_georec(sense * _m_f(lon), _m_f(lat), _m_f(alt), _m_f(re), _m_f(f))
end

# --- Jacobians ------------------------------------------------------------------------

@noinline _m_onzaxis() = throw(SpiceError(:pointonzaxis,
    "the Jacobian is not defined for a point on the z axis"))

@inline function _m_drdlat(r, lon, lat)
    slon, clon = sincos(lon)
    slat, clat = sincos(lat)
    return SMatrix{3,3}(clon * clat, slon * clat, slat + zero(r),
        -r * slon * clat, r * clon * clat, zero(r * clat),
        -r * clon * slat, -r * slon * slat, r * clat)
end

"""
    drdlat(r, lon, lat) -> SMatrix{3,3}

Return the Jacobian of rectangular coordinates with respect to latitudinal coordinates
(radius, lon, lat), as `drdlat_c`.
"""
drdlat(r::Real, lon::Real, lat::Real) = _m_drdlat(promote(_m_f(r), _m_f(lon), _m_f(lat))...)

"Inverse of a Jacobian with orthogonal columns, with the checks of `invort_`."
@inline function _m_invort(m::SMatrix{3,3})
    unit(i) = begin
        u, len = _m_unorm(_m_col(m, i))
        l0 = primal(len)
        l0 == 0.0 && throw(SpiceError(:zerolengthcolumn, "column $i of the matrix has norm zero"))
        (l0 < 1.0 && l0 * floatmax(Float64) < 1.0) && throw(SpiceError(:columntoosmall,
            "the length $l0 of column $i cannot be inverted"))
        (1.0 / len) * u
    end
    return _m_rows(unit(1), unit(2), unit(3))
end

"""
    dlatdr(x, y, z) -> SMatrix{3,3}

Return the Jacobian of latitudinal coordinates (radius, lon, lat) with respect to
rectangular coordinates, as `dlatdr_c`. Requires a point that is not on the z axis
(`SpiceError(:pointonzaxis)`).
"""
function dlatdr(x::Real, y::Real, z::Real)
    v = SVector(promote(_m_f(x), _m_f(y), _m_f(z))...)
    (primal(v[1]) == 0.0 && primal(v[2]) == 0.0) && _m_onzaxis()
    r, lon, lat = _m_reclat(v, false)
    return _m_invort(_m_drdlat(r, lon, lat))
end

@inline function _m_drdsph(r, colat, slon)
    scolat, ccolat = sincos(colat)
    slong, clong = sincos(slon)
    return SMatrix{3,3}(clong * scolat, slong * scolat, ccolat + zero(r),
        r * clong * ccolat, r * slong * ccolat, -r * scolat,
        -r * slong * scolat, r * clong * scolat, zero(r * scolat))
end

"""
    drdsph(r, colat, slon) -> SMatrix{3,3}

Return the Jacobian of rectangular coordinates with respect to spherical coordinates
(r, colat, slon), as `drdsph_c`.
"""
drdsph(r::Real, colat::Real, slon::Real) = _m_drdsph(promote(_m_f(r), _m_f(colat), _m_f(slon))...)

"""
    dsphdr(x, y, z) -> SMatrix{3,3}

Return the Jacobian of spherical coordinates (r, colat, slon) with respect to
rectangular coordinates, as `dsphdr_c`. Requires a point that is not on the z axis
(`SpiceError(:pointonzaxis)`).
"""
function dsphdr(x::Real, y::Real, z::Real)
    v = SVector(promote(_m_f(x), _m_f(y), _m_f(z))...)
    (primal(v[1]) == 0.0 && primal(v[2]) == 0.0) && _m_onzaxis()
    r, colat, slon = _m_recsph(v)
    return _m_invort(_m_drdsph(r, colat, slon))
end

@inline function _m_drdcyl(r, clon, z)
    s, c = sincos(clon)
    z0 = zero(r * c)
    return SMatrix{3,3}(c + z0, s + z0, z0, -s * r, c * r, z0, z0, z0, one(z0))
end

"""
    drdcyl(r, clon, z) -> SMatrix{3,3}

Return the Jacobian of rectangular coordinates with respect to cylindrical coordinates
(r, clon, z), as `drdcyl_c`.
"""
drdcyl(r::Real, clon::Real, z::Real) = _m_drdcyl(promote(_m_f(r), _m_f(clon), _m_f(z))...)

"""
    dcyldr(x, y, z) -> SMatrix{3,3}

Return the Jacobian of cylindrical coordinates (r, clon, z) with respect to rectangular
coordinates, as `dcyldr_c`. Requires a point that is not on the z axis
(`SpiceError(:pointonzaxis)`).
"""
function dcyldr(x::Real, y::Real, z::Real)
    v = SVector(promote(_m_f(x), _m_f(y), _m_f(z))...)
    (primal(v[1]) == 0.0 && primal(v[2]) == 0.0) && _m_onzaxis()
    r, clon, zz = _m_reccyl(v)
    return _m_invort(_m_drdcyl(r, clon, zz))
end

@inline function _m_drdgeo(lon, lat, alt, re, f)
    primal(f) < 1.0 || _m_badflat(f)
    primal(re) > 0.0 || throw(SpiceError(:badradius,
        "the equatorial radius must be positive; it is $(primal(re))"))
    slat, clat = sincos(lat)
    slon, clon = sincos(lon)
    flat = 1.0 - f
    flat2 = flat * flat
    g = sqrt(clat * clat + flat2 * slat * slat)
    g2 = g * g
    dgdlat = (flat2 - 1.0) * slat * clat / g
    return SMatrix{3,3}(
        -(alt + re / g) * slon * clat,
        (alt + re / g) * clon * clat,
        zero(alt + re / g),
        -re * dgdlat / g2 * clon * clat - (alt + re / g) * clon * slat,
        -re * dgdlat / g2 * slon * clat - (alt + re / g) * slon * slat,
        -flat2 * re * dgdlat / g2 * slat + (alt + flat2 * re / g) * clat,
        clon * clat + zero(alt + re / g),
        slon * clat + zero(alt + re / g),
        slat + zero(alt + re / g))
end

"""
    drdgeo(lon, lat, alt, re, f) -> SMatrix{3,3}

Return the Jacobian of rectangular coordinates with respect to geodetic coordinates
(lon, lat, alt) on the spheroid `re`, `f`, as `drdgeo_c`. Requires `f < 1`
(`SpiceError(:valueoutofrange)`) and `re > 0` (`SpiceError(:badradius)`).
"""
drdgeo(lon::Real, lat::Real, alt::Real, re::Real, f::Real) =
    _m_drdgeo(promote(_m_f(lon), _m_f(lat), _m_f(alt), _m_f(re), _m_f(f))...)

@inline function _m_dgeodr(v::SVector{3}, re, f)
    primal(f) < 1.0 || _m_badflat(f)
    primal(re) > 0.0 || throw(SpiceError(:badradius,
        "the equatorial radius must be positive; it is $(primal(re))"))
    (primal(v[1]) == 0.0 && primal(v[2]) == 0.0) && _m_onzaxis()
    lon, lat, alt = _m_recgeo(v, re, f)
    return _m_invort(_m_drdgeo(promote(lon, lat, alt, re, f)...))
end

"""
    dgeodr(x, y, z, re, f) -> SMatrix{3,3}

Return the Jacobian of geodetic coordinates (lon, lat, alt) with respect to rectangular
coordinates, as `dgeodr_c`. Errors as `drdgeo`, and `SpiceError(:pointonzaxis)` for a
point on the z axis.
"""
dgeodr(x::Real, y::Real, z::Real, re::Real, f::Real) =
    _m_dgeodr(SVector(promote(_m_f(x), _m_f(y), _m_f(z))...), _m_f(re), _m_f(f))

"""
    drdpgr(body, lon, lat, alt, re, f) -> SMatrix{3,3}

Return the Jacobian of rectangular coordinates with respect to the planetographic
coordinates (lon, lat, alt) of `body`, as `drdpgr_c`. Errors as `pgrrec`.
"""
function drdpgr(ks::KernelSet, body::NameOrID, lon::Real, lat::Real, alt::Real, re::Real, f::Real)
    id = _m_bodyid(ks, body, :idcodenotfound)
    primal(re) > 0.0 || _m_badre(re)
    primal(f) < 1.0 || _m_badflat(f)
    sense = _m_pgrsense(ks, id)
    j = _m_drdgeo(promote(sense * _m_f(lon), _m_f(lat), _m_f(alt), _m_f(re), _m_f(f))...)
    return _m_cols(sense * _m_col(j, 1), _m_col(j, 2), _m_col(j, 3))
end

"""
    dpgrdr(body, x, y, z, re, f) -> SMatrix{3,3}

Return the Jacobian of the planetographic coordinates (lon, lat, alt) of `body` with
respect to rectangular coordinates, as `dpgrdr_c`. Errors as `pgrrec`, and
`SpiceError(:pointonzaxis)` for a point on the z axis.
"""
function dpgrdr(ks::KernelSet, body::NameOrID, x::Real, y::Real, z::Real, re::Real, f::Real)
    id = _m_bodyid(ks, body, :idcodenotfound)
    primal(re) > 0.0 || _m_badre(re)
    primal(f) < 1.0 || _m_badflat(f)
    sense = _m_pgrsense(ks, id)
    j = _m_dgeodr(SVector(promote(_m_f(x), _m_f(y), _m_f(z))...), _m_f(re), _m_f(f))
    return _m_rows(sense * _m_row(j, 1), _m_row(j, 2), _m_row(j, 3))
end

# Context and default-context forms of the planetographic calls.
recpgr(ctx::SpiceContext, body::NameOrID, rectan::AbstractVector, re::Real, f::Real) =
    recpgr(kernelset(ctx), body, rectan, re, f)
recpgr(body::NameOrID, rectan::AbstractVector, re::Real, f::Real) =
    recpgr(default_context(), body, rectan, re, f)
for fn in (:pgrrec, :drdpgr, :dpgrdr)
    @eval begin
        $fn(ctx::SpiceContext, body::NameOrID, a::Real, b::Real, c::Real, re::Real, f::Real) =
            $fn(kernelset(ctx), body, a, b, c, re, f)
        $fn(body::NameOrID, a::Real, b::Real, c::Real, re::Real, f::Real) =
            $fn(default_context(), body, a, b, c, re, f)
    end
end

# --- state between coordinate systems -------------------------------------------------

const _M_COSYS = ("RECTANGULAR", "CYLINDRICAL", "LATITUDINAL", "SPHERICAL", "GEODETIC",
    "PLANETOGRAPHIC")

"Index of the coordinate system `name` in `_M_COSYS` (blanks and case ignored), or 0."
function _m_cosys(name::AbstractString)
    s = uppercase(filter(!isspace, name))
    return something(findfirst(==(s), _M_COSYS), 0)
end

"`SpiceError(:numericoverflow)` if a product of `jac` and `v` would overflow, as `xfmsta_`."
@inline function _m_xfm_overflow(jac::SMatrix{3,3}, v::SVector{3})
    toobig = sqrt(floatmax(Float64) / 100.0)
    for i in 1:3, j in 1:3
        sqrt(abs(primal(jac[i, j]))) * sqrt(abs(primal(v[j]))) > toobig &&
            throw(SpiceError(:numericoverflow,
                "a product of a Jacobian element and a velocity component would overflow"))
    end
    return nothing
end

"""
    xfmsta(istate, icosys, ocosys, body) -> SVector{6}

Return the state `istate` (three coordinates and their rates) of the coordinate system
`icosys` in the system `ocosys`, as `xfmsta_c`. The systems are `"RECTANGULAR"`,
`"CYLINDRICAL"`, `"LATITUDINAL"`, `"SPHERICAL"`, `"GEODETIC"` and `"PLANETOGRAPHIC"`
(case and blanks ignored). `body` (name or ID) gives the spheroid for the last two, from
`BODY<id>_RADII` of the context; it is not read otherwise. Throws
`SpiceError(:coordsysnotrec)` for an unknown system, `SpiceError(:idcodenotfound)` for
an unknown body, `SpiceError(:notsupported)` for a body whose first two radii differ,
`SpiceError(:invalidstate)` for a position on the z axis with a velocity off the axis
(the output rates are then not defined), and `SpiceError(:numericoverflow)`.
"""
function xfmsta(ks::KernelSet, istate::AbstractVector, icosys::AbstractString,
        ocosys::AbstractString, body::NameOrID)
    s = _m_v6(istate)
    isys = _m_cosys(icosys)
    osys = _m_cosys(ocosys)
    (isys == 0 || osys == 0) && throw(SpiceError(:coordsysnotrec,
        "a coordinate system of \"$icosys\", \"$ocosys\" is not recognized"))
    isys == osys && return s / 1.0
    T = typeof(atan(s[1] / 1.0))
    re = 1.0
    f = 0.0
    if isys >= 5 || osys >= 5
        id = _m_bodyid(ks, body, :idcodenotfound)
        r1, r2, r3 = _m_radii(ks, id)
        sqrt(abs(r1 - r3)) / sqrt(abs(r1)) >= sqrt(floatmax(Float64)) &&
            throw(SpiceError(:invalidradius, "the radii of body $body give no flattening"))
        r1 == r2 || throw(SpiceError(:notsupported,
            "body $body has two different equatorial radii; a spheroid is required"))
        re = r1
        f = (r1 - r3) / r1
    end
    c1, c2, c3 = s[1], s[2], s[3]
    rates = SVector(s[4], s[5], s[6])
    if isys == 1
        ipos = SVector{3,T}(c1, c2, c3)
        ivel = SVector{3,T}(rates)
    else
        if isys == 2
            ipos = SVector{3,T}(_m_cylrec(c1, c2, c3))
            jac = SMatrix{3,3,T,9}(_m_drdcyl(c1, c2, c3))
        elseif isys == 3
            ipos = SVector{3,T}(_m_latrec(c1, c2, c3))
            jac = SMatrix{3,3,T,9}(_m_drdlat(c1, c2, c3))
        elseif isys == 4
            ipos = SVector{3,T}(_m_sphrec(c1, c2, c3))
            jac = SMatrix{3,3,T,9}(_m_drdsph(c1, c2, c3))
        elseif isys == 5
            ipos = SVector{3,T}(_m_georec(c1, c2, c3, re, f))
            jac = SMatrix{3,3,T,9}(_m_drdgeo(promote(c1, c2, c3, re, f)...))
        else
            ipos = SVector{3,T}(pgrrec(ks, body, c1, c2, c3, re, f))
            jac = SMatrix{3,3,T,9}(drdpgr(ks, body, c1, c2, c3, re, f))
        end
        _m_xfm_overflow(jac, rates)
        ivel = SVector{3,T}(_m_mxv(jac, rates))
    end
    osys == 1 && return _m_state(ipos, ivel)
    if abs(primal(ipos[1])) + abs(primal(ipos[2])) == 0.0
        abs(primal(ivel[1])) + abs(primal(ivel[2])) == 0.0 || throw(SpiceError(:invalidstate,
            "the position is on the z axis and the velocity is not: the output rates are not defined"))
        z0 = zero(T)
        if osys == 2
            return _m_state(SVector{3,T}(_m_reccyl(ipos)), SVector(z0, z0, ivel[3]))
        elseif osys == 3
            return _m_state(SVector{3,T}(_m_reclat(ipos, false)), SVector(ivel[3], z0, z0))
        elseif osys == 4
            return _m_state(SVector{3,T}(_m_recsph(ipos)), SVector(ivel[3], z0, z0))
        elseif osys == 5
            return _m_state(SVector{3,T}(_m_recgeo(ipos, re, f)), SVector(z0, z0, ivel[3]))
        end
        return _m_state(SVector{3,T}(recpgr(ks, body, ipos, re, f)), SVector(z0, z0, ivel[3]))
    end
    if osys == 2
        jac = SMatrix{3,3,T,9}(_m_invort(_m_drdcyl(_m_reccyl(ipos)...)))
        opos = SVector{3,T}(_m_reccyl(ipos))
    elseif osys == 3
        jac = SMatrix{3,3,T,9}(_m_invort(_m_drdlat(_m_reclat(ipos, false)...)))
        opos = SVector{3,T}(_m_reclat(ipos, false))
    elseif osys == 4
        jac = SMatrix{3,3,T,9}(_m_invort(_m_drdsph(_m_recsph(ipos)...)))
        opos = SVector{3,T}(_m_recsph(ipos))
    elseif osys == 5
        jac = SMatrix{3,3,T,9}(_m_dgeodr(ipos, re, f))
        opos = SVector{3,T}(_m_recgeo(ipos, re, f))
    else
        jac = SMatrix{3,3,T,9}(dpgrdr(ks, body, ipos[1], ipos[2], ipos[3], re, f))
        opos = SVector{3,T}(recpgr(ks, body, ipos, re, f))
    end
    _m_xfm_overflow(jac, ivel)
    return _m_state(opos, _m_mxv(jac, ivel))
end
xfmsta(ctx::SpiceContext, istate::AbstractVector, icosys::AbstractString, ocosys::AbstractString,
    body::NameOrID) = xfmsta(kernelset(ctx), istate, icosys, ocosys, body)
xfmsta(istate::AbstractVector, icosys::AbstractString, ocosys::AbstractString, body::NameOrID) =
    xfmsta(default_context(), istate, icosys, ocosys, body)

# --- units ----------------------------------------------------------------------------

"Unit names of `convrt_`, the factor to the base unit of the class, and the class."
const _M_UNITS = (
    ("RADIANS", 180.0 / Float64(pi), :angle), ("DEGREES", 1.0, :angle),
    ("ARCMINUTES", 0.016666666666666666, :angle), ("ARCSECONDS", 2.7777777777777778e-4, :angle),
    ("HOURANGLE", 15.0, :angle), ("MINUTEANGLE", 0.25, :angle),
    ("SECONDANGLE", 0.0041666666666666666, :angle),
    ("METERS", 1.0, :distance), ("KM", 1.0e3, :distance), ("CM", 0.01, :distance),
    ("MM", 0.001, :distance), ("LIGHTSECS", 299792458.0, :distance),
    ("AU", 149597870613.68887, :distance), ("M", 1.0, :distance),
    ("KILOMETERS", 1.0e3, :distance), ("CENTIMETERS", 0.01, :distance),
    ("MILLIMETERS", 0.001, :distance), ("FEET", 0.3048, :distance), ("INCHES", 0.0254, :distance),
    ("STATUTE_MILES", 1609.344, :distance), ("NAUTICAL_MILES", 1852.0, :distance),
    ("YARDS", 0.9144, :distance), ("LIGHTYEARS", 9460730472580800.0, :distance),
    ("PARSECS", 30856775797231604.0, :distance),
    ("SECONDS", 1.0, :time), ("MINUTES", 60.0, :time), ("HOURS", 3600.0, :time),
    ("DAYS", 86400.0, :time), ("JULIAN_YEARS", 31557600.0, :time),
    ("TROPICAL_YEARS", 31556925.976319999, :time), ("YEARS", 31557600.0, :time))

"Index of the unit `name` in `_M_UNITS` (case ignored, trailing blanks ignored), or 0."
function _m_unit(name::AbstractString)
    s = uppercase(rstrip(name, ' '))
    return something(findfirst(u -> u[1] == s, _M_UNITS), 0)
end

"""
    convrt(x, in, out)

Return the measurement `x` in the units `in` converted to the units `out`, as
`convrt_c`. The units are angles (`RADIANS`, `DEGREES`, `ARCMINUTES`, `ARCSECONDS`,
`HOURANGLE`, `MINUTEANGLE`, `SECONDANGLE`), distances (`METERS`, `M`, `KM`,
`KILOMETERS`, `CM`, `CENTIMETERS`, `MM`, `MILLIMETERS`, `FEET`, `INCHES`, `YARDS`,
`STATUTE_MILES`, `NAUTICAL_MILES`, `AU`, `LIGHTSECS`, `LIGHTYEARS`, `PARSECS`) and times
(`SECONDS`, `MINUTES`, `HOURS`, `DAYS`, `JULIAN_YEARS`, `TROPICAL_YEARS`, `YEARS`), in
any case. Throws `SpiceError(:unitsnotrec)` for an unknown unit and
`SpiceError(:incompatibleunits)` for units of two classes.
"""
function convrt(x::Real, in::AbstractString, out::AbstractString)
    i = _m_unit(in)
    j = _m_unit(out)
    (i == 0 || j == 0) && throw(SpiceError(:unitsnotrec,
        "the units \"$in\" or \"$out\" were not recognized"))
    _M_UNITS[i][3] === _M_UNITS[j][3] || throw(SpiceError(:incompatibleunits,
        "cannot convert $in ($(_M_UNITS[i][3])) to $out ($(_M_UNITS[j][3]))"))
    return _m_f(x) * _M_UNITS[i][2] / _M_UNITS[j][2]
end
