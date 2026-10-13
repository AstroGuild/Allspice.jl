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

Return the latitudinal coordinates of the rectangular point `rectan`, as `reclat_c` of
CSPICE.

`radius` has the unit of `rectan`. The angles are in radians: the longitude `lon` in
(-pi, pi], from the x axis toward the y axis, and the latitude `lat` in [-pi/2, pi/2],
from the xy plane. `lon` is 0 on the z axis, and all three are 0 at the origin.

Throws `SpiceError(:badarg)` if `rectan` does not have 3 elements.

# Examples
```jldoctest
julia> reclat([0, 2, 0])
(2.0, 1.5707963267948966, 0.0)

julia> latrec(2, pi / 2, 0) ≈ [0, 2, 0]
true
```

See also [`latrec`](@ref), [`recrad`](@ref), [`recsph`](@ref), [`recgeo`](@ref),
[`dlatdr`](@ref).
"""
reclat(rectan::AbstractVector) = _m_reclat(_m_v3(rectan), true)

@inline function _m_latrec(radius, lon, lat)
    slon, clon = sincos(lon)
    slat, clat = sincos(lat)
    return SVector(radius * clon * clat, radius * slon * clat, radius * slat)
end

"""
    latrec(radius, lon, lat) -> SVector{3}

Return the rectangular coordinates of the point with the latitudinal coordinates
`radius`, `lon`, `lat` (angles in radians), as `latrec_c` of CSPICE.

See also [`reclat`](@ref), which has an example, and [`drdlat`](@ref).
"""
latrec(radius::Real, lon::Real, lat::Real) = _m_latrec(_m_f(radius), _m_f(lon), _m_f(lat))

"""
    recrad(rectan) -> (range, ra, dec)

Return the range, the right ascension and the declination of the rectangular point
`rectan`, as `recrad_c` of CSPICE.

The angles are in radians: `ra` in [0, 2 pi) and `dec` in [-pi/2, pi/2]. They differ
from the result of [`reclat`](@ref) only in the range of the first angle.

Throws `SpiceError(:badarg)` if `rectan` does not have 3 elements.

# Examples
```jldoctest
julia> recrad([0, -2, 0])
(2.0, 4.71238898038469, 0.0)
```

See also [`radrec`](@ref), [`reclat`](@ref).
"""
function recrad(rectan::AbstractVector)
    range, ra, dec = _m_reclat(_m_v3(rectan), true)
    return range, (primal(ra) < 0.0 ? ra + _M_TWOPI : ra), dec
end

"""
    radrec(range, ra, dec) -> SVector{3}

Return the rectangular coordinates of the point with the range `range`, the right
ascension `ra` and the declination `dec` (angles in radians), as `radrec_c` of CSPICE.

See also [`recrad`](@ref).
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

Return the spherical coordinates of the rectangular point `rectan`, as `recsph_c` of
CSPICE.

The angles are in radians: the colatitude `colat` in [0, pi], from the positive z axis,
and the longitude `slon` in (-pi, pi]. `slon` is 0 on the z axis, and all three are 0 at
the origin.

Throws `SpiceError(:badarg)` if `rectan` does not have 3 elements.

# Examples
```jldoctest
julia> recsph([0, 0, -2])
(2.0, 3.141592653589793, 0.0)
```

See also [`sphrec`](@ref), [`reclat`](@ref), [`dsphdr`](@ref).
"""
recsph(rectan::AbstractVector) = _m_recsph(_m_v3(rectan))

@inline function _m_sphrec(r, colat, slon)
    sl, cl = sincos(slon)
    sc, cc = sincos(colat)
    return SVector(r * cl * sc, r * sl * sc, r * cc)
end

"""
    sphrec(r, colat, slon) -> SVector{3}

Return the rectangular coordinates of the point with the spherical coordinates `r`,
`colat`, `slon` (angles in radians), as `sphrec_c` of CSPICE.

See also [`recsph`](@ref), [`drdsph`](@ref).
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

Return the cylindrical coordinates of the rectangular point `rectan`, as `reccyl_c` of
CSPICE.

`r` is the distance from the z axis. The longitude `clon` is in radians in [0, 2 pi),
and it is 0 on the z axis.

Throws `SpiceError(:badarg)` if `rectan` does not have 3 elements.

# Examples
```jldoctest
julia> reccyl([0, -2, 3])
(2.0, 4.71238898038469, 3.0)
```

See also [`cylrec`](@ref), [`dcyldr`](@ref).
"""
reccyl(rectan::AbstractVector) = _m_reccyl(_m_v3(rectan))

@inline function _m_cylrec(r, clon, z)
    s, c = sincos(clon)
    return SVector(r * c, r * s, z - zero(r * c))
end

"""
    cylrec(r, clon, z) -> SVector{3}

Return the rectangular coordinates of the point with the cylindrical coordinates `r`,
`clon`, `z` (`clon` in radians), as `cylrec_c` of CSPICE.

See also [`reccyl`](@ref), [`drdcyl`](@ref).
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

Return the latitudinal coordinates of the point with the cylindrical coordinates `r`,
`clon`, `z`, as `cyllat_c` of CSPICE. Angles are in radians. The longitude is returned
as given.

See also [`latcyl`](@ref), [`reclat`](@ref), [`reccyl`](@ref).
"""
function cyllat(r::Real, clon::Real, z::Real)
    r, clon, z = promote(_m_f(r), _m_f(clon), _m_f(z))
    rho = _m_hypot(r, z)
    lat = primal(rho) == 0.0 ? zero(atan(z, r)) : atan(z, r)
    return rho - zero(lat), clon - zero(lat), lat
end

"""
    latcyl(radius, lon, lat) -> (r, clon, z)

Return the cylindrical coordinates of the point with the latitudinal coordinates
`radius`, `lon`, `lat`, as `latcyl_c` of CSPICE. Angles are in radians. The longitude
is returned as given.

See also [`cyllat`](@ref).
"""
function latcyl(radius::Real, lon::Real, lat::Real)
    radius, lon, lat = promote(_m_f(radius), _m_f(lon), _m_f(lat))
    s, c = sincos(lat)
    return radius * c, lon - zero(c), radius * s
end

"""
    cylsph(r, clon, z) -> (radius, colat, slon)

Return the spherical coordinates of the point with the cylindrical coordinates `r`,
`clon`, `z`, as `cylsph_c` of CSPICE. Angles are in radians. The longitude is returned
as given.

See also [`sphcyl`](@ref), [`recsph`](@ref).
"""
function cylsph(r::Real, clon::Real, z::Real)
    r, clon, z = promote(_m_f(r), _m_f(clon), _m_f(z))
    rh = _m_hypot(r, z)
    th = primal(rh) == 0.0 ? zero(atan(r, z)) : atan(r, z)
    return rh - zero(th), th, clon - zero(th)
end

"""
    sphcyl(radius, colat, slon) -> (r, clon, z)

Return the cylindrical coordinates of the point with the spherical coordinates
`radius`, `colat`, `slon`, as `sphcyl_c` of CSPICE. Angles are in radians. The
longitude is returned as given.

See also [`cylsph`](@ref).
"""
function sphcyl(radius::Real, colat::Real, slon::Real)
    radius, colat, slon = promote(_m_f(radius), _m_f(colat), _m_f(slon))
    s, c = sincos(colat)
    return radius * s, slon - zero(c), radius * c
end

"""
    sphlat(r, colat, slon) -> (radius, lon, lat)

Return the latitudinal coordinates of the point with the spherical coordinates `r`,
`colat`, `slon`, as `sphlat_c` of CSPICE. Angles are in radians. The latitude is
`pi/2 - colat`, and the other two values are returned as given.

See also [`latsph`](@ref).
"""
function sphlat(r::Real, colat::Real, slon::Real)
    r, colat, slon = promote(_m_f(r), _m_f(colat), _m_f(slon))
    return r, slon, _M_HALFPI - colat
end

"""
    latsph(radius, lon, lat) -> (rho, colat, slon)

Return the spherical coordinates of the point with the latitudinal coordinates
`radius`, `lon`, `lat`, as `latsph_c` of CSPICE. Angles are in radians. The colatitude
is `pi/2 - lat`, and the other two values are returned as given.

The longitude is the second latitudinal coordinate and the third spherical one:
```jldoctest
julia> latsph(2.0, 0.5, 0.0)
(2.0, 1.5707963267948966, 0.5)
```

See also [`sphlat`](@ref).
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
    lon = (primal(v[1]) == 0.0 && primal(v[2]) == 0.0) ? zero(lat) : atan(v[2], v[1]) - zero(lat)
    return lon, lat, alt - zero(lat)
end

"""
    recgeo(rectan, re, f) -> (lon, lat, alt)

Return the geodetic coordinates of the rectangular point `rectan` relative to a
spheroid, as `recgeo_c` of CSPICE.

The spheroid has the equatorial radius `re` and the flattening `f`; its polar radius is
`re * (1 - f)`. `lon` and `lat` are in radians. `lat` is the angle between the
equatorial plane and the outward normal at the surface point nearest to `rectan`, which
[`nearpt`](@ref) finds. `alt` is the height above that point, in the unit of `re` and
`rectan`; it is negative inside the spheroid. `lon` is in (-pi, pi], and it is 0 on the
z axis.

Throws `SpiceError(:valueoutofrange)` if `re <= 0` or `f >= 1`, and
`SpiceError(:badarg)` if `rectan` does not have 3 elements.

# Examples
A point 100 km above the equator of the WGS 84 spheroid:
```jldoctest
julia> re, f = 6378.137, 1 / 298.257223563;

julia> round.(recgeo([0, 6478.137, 0], re, f); digits=6)
(1.570796, 0.0, 100.0)

julia> georec(pi / 2, 0, 100, re, f) ≈ [0, 6478.137, 0]
true
```

See also [`georec`](@ref), [`recpgr`](@ref), [`dgeodr`](@ref).
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

Return the rectangular coordinates of the point with the geodetic coordinates `lon`,
`lat` and `alt`, as `georec_c` of CSPICE.

The spheroid has the equatorial radius `re` and the flattening `f`. Angles are in
radians. The result has the unit of `re` and `alt`.

Throws `SpiceError(:valueoutofrange)` if `re <= 0` or `f >= 1`.

See also [`recgeo`](@ref), which has an example, and [`drdgeo`](@ref).
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
    recpgr([ctx,] body, rectan, re, f) -> (lon, lat, alt)

Return the planetographic coordinates of the rectangular point `rectan` relative to the
spheroid of `body`, as `recpgr_c` of CSPICE.

`body` is a name or an ID. `re` is the equatorial radius and `f` the flattening of the
spheroid. `lat` (radians) and `alt` are the geodetic values of [`recgeo`](@ref). `lon`
is in radians in [0, 2 pi], and its sense comes from the kernel pool of the context:

- The variable `BODY<id>_PGR_POSITIVE_LON` (`EAST` or `WEST`) gives the sense, if it is
  present.
- Without it, the longitude is positive west for prograde rotation (the rate
  `BODY<id>_PM[2]` is positive) and positive east for retrograde rotation.
- For the Earth, the Moon and the Sun it is positive east.

Throws:
- `SpiceError(:idcodenotfound)` if `body` is a name without an ID.
- `SpiceError(:valueoutofrange)` if `re <= 0` or `f >= 1`.
- `SpiceError(:invalidoption)` if the sense variable is not `EAST` or `WEST`.
- `SpiceError(:missingdata)` if the sense variable and `BODY<id>_PM` are both absent.
  The Earth needs neither.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pcpool(ctx, "BODY499_PGR_POSITIVE_LON", ["WEST"]);

julia> round.(recpgr(ctx, "MARS", [0, 3496.19, 0], 3396.19, 0.00589); digits=6)
(4.712389, 0.0, 100.0)

julia> pgrrec(ctx, "MARS", 3pi / 2, 0, 100, 3396.19, 0.00589) ≈ [0, 3496.19, 0]
true
```

See also [`pgrrec`](@ref), [`recgeo`](@ref), [`dpgrdr`](@ref).
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
    primal(lon) > _M_TWOPI && (lon = _M_TWOPI - zero(lon))
    return lon, lat, alt
end

"""
    pgrrec([ctx,] body, lon, lat, alt, re, f) -> SVector{3}

Return the rectangular coordinates of the point with the planetographic coordinates
`lon`, `lat` and `alt` of `body`, as `pgrrec_c` of CSPICE.

`re` is the equatorial radius and `f` the flattening of the spheroid. Angles are in
radians. The kernel pool of the context gives the sense of the longitude, by the rules in
[`recpgr`](@ref).

Throws as [`recpgr`](@ref).

See also [`recpgr`](@ref), which has an example, and [`drdpgr`](@ref).
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
    return SMatrix{3,3}(clon * clat, slon * clat, slat - zero(r),
        -r * slon * clat, r * clon * clat, zero(r * clat),
        -r * clon * slat, -r * slon * slat, r * clat)
end

"""
    drdlat(r, lon, lat) -> SMatrix{3,3}

Return the Jacobian of rectangular coordinates with respect to the latitudinal
coordinates `(r, lon, lat)`, as `drdlat_c` of CSPICE.

Element `[i, j]` is the derivative of rectangular coordinate `i` with respect to
latitudinal coordinate `j`. Angles are in radians.

See also [`dlatdr`](@ref), [`latrec`](@ref).
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

Return the Jacobian of the latitudinal coordinates `(radius, lon, lat)` with respect to
rectangular coordinates, at the point `(x, y, z)`, as `dlatdr_c` of CSPICE.

Element `[i, j]` is the derivative of latitudinal coordinate `i` with respect to
rectangular coordinate `j`. Angles are in radians.

Throws `SpiceError(:pointonzaxis)` if the point is on the z axis.

# Examples
The two Jacobians at one point are inverses:
```jldoctest
julia> dlatdr(1, 1, 0) * drdlat(reclat([1, 1, 0])...) ≈ ident()
true
```

See also [`drdlat`](@ref), [`reclat`](@ref).
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
    return SMatrix{3,3}(clong * scolat, slong * scolat, ccolat - zero(r),
        r * clong * ccolat, r * slong * ccolat, -r * scolat,
        -r * slong * scolat, r * clong * scolat, zero(r * scolat))
end

"""
    drdsph(r, colat, slon) -> SMatrix{3,3}

Return the Jacobian of rectangular coordinates with respect to the spherical
coordinates `(r, colat, slon)`, as `drdsph_c` of CSPICE. Angles are in radians.

See also [`dsphdr`](@ref), [`sphrec`](@ref), [`drdlat`](@ref).
"""
drdsph(r::Real, colat::Real, slon::Real) = _m_drdsph(promote(_m_f(r), _m_f(colat), _m_f(slon))...)

"""
    dsphdr(x, y, z) -> SMatrix{3,3}

Return the Jacobian of the spherical coordinates `(r, colat, slon)` with respect to
rectangular coordinates, at the point `(x, y, z)`, as `dsphdr_c` of CSPICE. Angles are
in radians.

Throws `SpiceError(:pointonzaxis)` if the point is on the z axis.

See also [`drdsph`](@ref), [`recsph`](@ref), [`dlatdr`](@ref).
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

Return the Jacobian of rectangular coordinates with respect to the cylindrical
coordinates `(r, clon, z)`, as `drdcyl_c` of CSPICE. `clon` is in radians.

See also [`dcyldr`](@ref), [`cylrec`](@ref), [`drdlat`](@ref).
"""
drdcyl(r::Real, clon::Real, z::Real) = _m_drdcyl(promote(_m_f(r), _m_f(clon), _m_f(z))...)

"""
    dcyldr(x, y, z) -> SMatrix{3,3}

Return the Jacobian of the cylindrical coordinates `(r, clon, z)` with respect to
rectangular coordinates, at the point `(x, y, z)`, as `dcyldr_c` of CSPICE. `clon` is
in radians.

Throws `SpiceError(:pointonzaxis)` if the point is on the z axis.

See also [`drdcyl`](@ref), [`reccyl`](@ref), [`dlatdr`](@ref).
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
        clon * clat - zero(alt + re / g),
        slon * clat - zero(alt + re / g),
        slat - zero(alt + re / g))
end

"""
    drdgeo(lon, lat, alt, re, f) -> SMatrix{3,3}

Return the Jacobian of rectangular coordinates with respect to the geodetic coordinates
`(lon, lat, alt)`, as `drdgeo_c` of CSPICE.

The spheroid has the equatorial radius `re` and the flattening `f`. Angles are in
radians.

Throws `SpiceError(:valueoutofrange)` if `f >= 1`, and `SpiceError(:badradius)` if
`re <= 0`. [`georec`](@ref) throws `SpiceError(:valueoutofrange)` in both cases.

See also [`dgeodr`](@ref), [`georec`](@ref), [`drdlat`](@ref).
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

Return the Jacobian of the geodetic coordinates `(lon, lat, alt)` with respect to
rectangular coordinates, at the point `(x, y, z)`, as `dgeodr_c` of CSPICE.

The spheroid has the equatorial radius `re` and the flattening `f`. Angles are in
radians.

Throws `SpiceError(:valueoutofrange)` if `f >= 1`, `SpiceError(:badradius)` if
`re <= 0`, and `SpiceError(:pointonzaxis)` if the point is on the z axis.

See also [`drdgeo`](@ref), [`recgeo`](@ref), [`dlatdr`](@ref).
"""
dgeodr(x::Real, y::Real, z::Real, re::Real, f::Real) =
    _m_dgeodr(SVector(promote(_m_f(x), _m_f(y), _m_f(z))...), _m_f(re), _m_f(f))

"""
    drdpgr([ctx,] body, lon, lat, alt, re, f) -> SMatrix{3,3}

Return the Jacobian of rectangular coordinates with respect to the planetographic
coordinates `(lon, lat, alt)` of `body`, as `drdpgr_c` of CSPICE. Angles are in
radians. The kernel pool of the context gives the sense of the longitude, by the rules
in [`recpgr`](@ref).

Throws as [`recpgr`](@ref).

See also [`dpgrdr`](@ref), [`pgrrec`](@ref), [`drdgeo`](@ref).
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
    dpgrdr([ctx,] body, x, y, z, re, f) -> SMatrix{3,3}

Return the Jacobian of the planetographic coordinates `(lon, lat, alt)` of `body` with
respect to rectangular coordinates, at the point `(x, y, z)`, as `dpgrdr_c` of CSPICE.
Angles are in radians. The kernel pool of the context gives the sense of the longitude,
by the rules in [`recpgr`](@ref).

Throws as [`recpgr`](@ref), and `SpiceError(:pointonzaxis)` if the point is on the z
axis.

See also [`drdpgr`](@ref), [`recpgr`](@ref), [`dgeodr`](@ref).
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
    xfmsta([ctx,] istate, icosys, ocosys, body) -> SVector{6}

Return the state `istate` converted from the coordinate system `icosys` to the
coordinate system `ocosys`, as `xfmsta_c` of CSPICE.

A state is three coordinates followed by their time derivatives. The names of the
systems ignore case and blanks. Angles are in radians.

| System             | Coordinates          |
|:-------------------|:---------------------|
| `"RECTANGULAR"`    | x, y, z              |
| `"CYLINDRICAL"`    | radius, lon, z       |
| `"LATITUDINAL"`    | radius, lon, lat     |
| `"SPHERICAL"`      | radius, colat, lon   |
| `"GEODETIC"`       | lon, lat, alt        |
| `"PLANETOGRAPHIC"` | lon, lat, alt        |

`body` (name or ID) gives the spheroid of the last two systems, from `BODY<id>_RADII`
in the kernel pool of the context. It is not read for the other systems, so `""` is
enough there.

Throws:
- `SpiceError(:coordsysnotrec)` if a system name is not known.
- `SpiceError(:idcodenotfound)` if `body` is needed and is a name without an ID.
- `SpiceError(:kernelvarnotfound)` if `body` is needed and has no radii in the pool.
- `SpiceError(:notsupported)` if the first two radii of `body` differ.
- `SpiceError(:invalidstate)` if the position is on the z axis and the velocity is
  not along it; the output rates are then not defined.
- `SpiceError(:numericoverflow)` if a velocity component would overflow.
- `SpiceError(:badarg)` if `istate` does not have 6 elements.

# Examples
A point on the x axis that moves along y has a longitude rate:
```jldoctest
julia> s = xfmsta([2, 0, 0, 0, 1, 0], "RECTANGULAR", "LATITUDINAL", "");

julia> s == [2, 0, 0, 0, 0.5, 0]
true
```

See also [`reclat`](@ref), [`recgeo`](@ref), [`recpgr`](@ref), [`dlatdr`](@ref).
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

function _m_unit(name::AbstractString)
    s = uppercase(rstrip(name, ' '))
    return something(findfirst(u -> u[1] == s, _M_UNITS), 0)
end

"""
    convrt(x, in, out)

Return the measurement `x` converted from the unit `in` to the unit `out`, as
`convrt_c` of CSPICE.

| Class    | Units                                                                 |
|:---------|:----------------------------------------------------------------------|
| angle    | `RADIANS`, `DEGREES`, `ARCMINUTES`, `ARCSECONDS`, `HOURANGLE`, `MINUTEANGLE`, `SECONDANGLE` |
| distance | `METERS`, `M`, `KILOMETERS`, `KM`, `CENTIMETERS`, `CM`, `MILLIMETERS`, `MM`, `FEET`, `INCHES`, `YARDS`, `STATUTE_MILES`, `NAUTICAL_MILES`, `AU`, `LIGHTSECS`, `LIGHTYEARS`, `PARSECS` |
| time     | `SECONDS`, `MINUTES`, `HOURS`, `DAYS`, `JULIAN_YEARS`, `TROPICAL_YEARS`, `YEARS` |

The names ignore case and trailing blanks. `YEARS` is `JULIAN_YEARS`.

Throws `SpiceError(:unitsnotrec)` if a unit is not known, and
`SpiceError(:incompatibleunits)` if the two units are of different classes.

# Examples
```jldoctest
julia> convrt(1.5, "HOURS", "MINUTES")
90.0

julia> convrt(1, "km", "METERS")
1000.0
```
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
