# Instrument fields of view: `getfov`, `getfvn`, `fovray` and `fovtrg` of CSPICE N0067,
# with `zzgffvu` (the visibility test), `zzfovaxi` and `zzhullax` (the axis of a polygonal
# field of view), `zzelvupy` (ellipse against pyramid), `zzwind` and `zzwind2d` (winding
# numbers).
#
# The calls build the list of the boundary vectors, so they allocate. The epoch is a
# plain number.

"Largest number of boundary vectors that `fovray` and `fovtrg` accept (`MAXVRT` of CSPICE)."
const _GEO_MAXVRT = 10000

@noinline _geo_fov_error(kind::Symbol, msg) = throw(SpiceError(kind, msg))

@inline _geo_fov_var(ks::KernelSet, instid::Integer, suffix::String) =
    get(ks.pool, string("INS", instid, suffix), nothing)

@inline function _geo_fov_text(ks::KernelSet, instid::Integer, suffix::String)
    v = _geo_fov_var(ks, instid, suffix)
    v isa Vector{String} || return nothing
    isempty(v) && return nothing
    return @inbounds v[1]
end

function _geo_fov_vec3(ks::KernelSet, instid::Integer, suffix::String, missing::Symbol, bad::Symbol)
    v = _geo_fov_var(ks, instid, suffix)
    v === nothing && _geo_fov_error(missing,
        "the variable INS$(instid)$suffix is not in the kernel pool")
    v isa Vector{Float64} || _geo_fov_error(bad,
        "the variable INS$(instid)$suffix must hold 3 numbers")
    length(v) == 3 || _geo_fov_error(bad,
        "the variable INS$(instid)$suffix must hold 3 numbers")
    return @inbounds Vec3{Float64}(v[1], v[2], v[3])
end

"The angle `INS<instid><suffix>` in radians; `missing` is the error kind for a variable that is absent."
function _geo_fov_angle(ks::KernelSet, instid::Integer, suffix::String, missing::Symbol,
        units::Union{Nothing,String})
    v = _geo_fov_var(ks, instid, suffix)
    v isa Vector{Float64} || _geo_fov_error(missing,
        "the variable INS$(instid)$suffix is not in the kernel pool")
    isempty(v) && _geo_fov_error(missing,
        "the variable INS$(instid)$suffix is not in the kernel pool")
    units === nothing && _geo_fov_error(:unitsmissing,
        "the variable INS$(instid)_FOV_ANGLE_UNITS is not in the kernel pool")
    return Float64(convrt(@inbounds(v[1]), units, "RADIANS"))
end

function _geo_getfov(ks::KernelSet, instid::Integer, room::Integer)
    frame = _geo_fov_text(ks, instid, "_FOV_FRAME")
    frame === nothing && _geo_fov_error(:framemissing,
        "the variable INS$(instid)_FOV_FRAME is not in the kernel pool")
    shape0 = _geo_fov_text(ks, instid, "_FOV_SHAPE")
    shape0 === nothing && _geo_fov_error(:shapemissing,
        "the variable INS$(instid)_FOV_SHAPE is not in the kernel pool")
    shape = uppercase(String(strip(shape0, ' ')))
    shape in ("CIRCLE", "ELLIPSE", "POLYGON", "RECTANGLE") || _geo_fov_error(:shapenotsupported,
        "the field of view shape <$shape0> of the instrument $instid is not supported")
    bsight = _geo_fov_vec3(ks, instid, "_BORESIGHT", :boresightmissing, :badboresightspec)
    _m_iszero(bsight) && _geo_fov_error(:zeroboresight,
        "the boresight of the instrument $instid is the zero vector")
    spec = _geo_fov_text(ks, instid, "_FOV_CLASS_SPEC")
    spec === nothing && (spec = "CORNERS")
    bounds = Vec3{Float64}[]
    if _dsk_eqstr(spec, "CORNERS")
        v = _geo_fov_var(ks, instid, "_FOV_BOUNDARY")
        v === nothing && (v = _geo_fov_var(ks, instid, "_FOV_BOUNDARY_CORNERS"))
        v isa Vector{Float64} || _geo_fov_error(:boundarymissing,
            "the variable INS$(instid)_FOV_BOUNDARY is not in the kernel pool")
        n = length(v)
        n > 3 * room && _geo_fov_error(:boundarytoobig,
            "the boundary of the instrument $instid has $n components; the room is for $(3 * room)")
        n % 3 == 0 || _geo_fov_error(:badboundary,
            "the boundary of the instrument $instid has $n components, which is not a multiple of 3")
        n = n ÷ 3
        ((shape == "CIRCLE" && n != 1) || (shape == "ELLIPSE" && n != 2) ||
         (shape == "RECTANGLE" && n != 4) || (shape == "POLYGON" && n < 3)) &&
            _geo_fov_error(:badboundary,
                "the boundary of the instrument $instid has $n vectors, which is not valid for the shape $shape")
        for i in 1:n                     # invariant: bounds holds the vectors 1 .. i-1
            @inbounds push!(bounds, Vec3{Float64}(v[3i-2], v[3i-1], v[3i]))
        end
    elseif _dsk_eqstr(spec, "ANGLES")
        shape == "POLYGON" && _geo_fov_error(:shapenotsupported,
            "the shape POLYGON of the instrument $instid cannot be given by angles")
        refvec = _geo_fov_vec3(ks, instid, "_FOV_REF_VECTOR", :refvectormissing, :badrefvectorspec)
        normal = _m_cross(bsight, refvec)
        _m_norm(normal) == 0.0 && _geo_fov_error(:badrefvectorspec,
            "the reference vector of the instrument $instid is parallel to the boresight")
        units = _geo_fov_text(ks, instid, "_FOV_ANGLE_UNITS")
        refang = _geo_fov_angle(ks, instid, "_FOV_REF_ANGLE", :refanglemissing, units)
        if shape == "CIRCLE"
            room < 1 && _geo_fov_error(:boundarytoobig, "the room is $room; 1 vector is returned")
            push!(bounds, _m_vrotv(bsight, normal, refang))
        elseif shape == "ELLIPSE"
            crsang = _geo_fov_angle(ks, instid, "_FOV_CROSS_ANGLE", :crossanglemissing, units)
            room < 2 && _geo_fov_error(:boundarytoobig, "the room is $room; 2 vectors are returned")
            push!(bounds, _m_vrotv(bsight, normal, refang))
            push!(bounds, _m_vrotv(bsight, _m_cross(bsight, normal), crsang))
        else
            crsang = _geo_fov_angle(ks, instid, "_FOV_CROSS_ANGLE", :crossanglemissing, units)
            room < 4 && _geo_fov_error(:boundarytoobig, "the room is $room; 4 vectors are returned")
            b, bmag = _m_unorm(bsight)
            b1 = _m_hat(_m_perp(refvec, bsight, true))
            b2 = _m_cross(b, b1)
            cosran = cos(refang)
            coscan = cos(crsang)
            (cosran < 1e-15 || coscan < 1e-15) && _geo_fov_error(:badboundary,
                "the half angles of the instrument $instid must be below 90 degrees")
            sinran = sin(refang)
            sincan = sin(crsang)
            n1 = -cosran * b1 + sinran * b
            n2 = -coscan * b2 + sincan * b
            n3 = cosran * b1 + sinran * b
            n4 = coscan * b2 + sincan * b
            for c in (_m_cross(n1, n2), _m_cross(n2, n3), _m_cross(n3, n4), _m_cross(n4, n1))
                push!(bounds, bmag * _m_hat(c))
            end
        end
    else
        _geo_fov_error(:unsupportedspec,
            "the field of view specification <$spec> of the instrument $instid is not supported")
    end
    return shape, String(frame), bsight, bounds
end

"""
    getfov([ctx,] instid, room=10000) -> (shape, frame, bsight, bounds)

Return the field of view of the instrument with the integer ID `instid`, as `getfov_c`
of CSPICE.

- `shape`: `"CIRCLE"`, `"ELLIPSE"`, `"RECTANGLE"` or `"POLYGON"`.
- `frame`: the name of the frame of the vectors.
- `bsight`: the boresight vector, an `SVector{3}`.
- `bounds`: a new `Vector` of boundary vectors (`SVector{3}`): 1 for a circle, 2 for an
  ellipse, 4 for a rectangle, the corners for a polygon.

The vectors are not unit vectors in general; only their directions matter. CSPICE also
returns the number of boundary vectors; here it is `length(bounds)`.

The data come from the pool variables `INS<instid>_FOV_FRAME`, `INS<instid>_FOV_SHAPE`,
`INS<instid>_BORESIGHT` and the other `INS<instid>_FOV_*` variables of an instrument
kernel. `INS<instid>_FOV_CLASS_SPEC` selects the form: `CORNERS` (the default) gives
the boundary vectors, `ANGLES` gives angles from the boresight. `room` is the largest
number of boundary vectors accepted.

Throws `SpiceError` with the kind
- `:framemissing`, `:shapemissing`, `:boresightmissing` or `:boundarymissing` if the
  frame, the shape, the boresight or the boundary variable is not in the kernel pool;
- `:shapenotsupported` if the shape is not one of the four, or is a polygon given by
  angles; `:unsupportedspec` if the specification is not `CORNERS` or `ANGLES`;
- `:badboresightspec` if the boresight does not hold 3 numbers, `:zeroboresight` if it
  is zero;
- `:badboundary` if the boundary does not hold a multiple of 3 numbers, if the number
  of vectors does not fit the shape, or if a half angle of a rectangle is 90 degrees
  or more;
- `:boundarytoobig` if there are more than `room` boundary vectors;
- `:refvectormissing`, `:refanglemissing`, `:crossanglemissing` or `:unitsmissing` if a
  variable of the `ANGLES` form is absent, `:badrefvectorspec` if the reference vector
  does not hold 3 numbers or is parallel to the boresight.

Also throws the errors of [`convrt`](@ref) for the angle unit.

See also [`getfvn`](@ref), [`fovtrg`](@ref), [`fovray`](@ref).
"""
getfov(ks::KernelSet, instid::Integer, room::Integer=_GEO_MAXVRT) = _geo_getfov(ks, instid, room)
getfov(ctx::SpiceContext, instid::Integer, room::Integer=_GEO_MAXVRT) =
    _geo_getfov(kernelset(ctx), instid, room)
getfov(instid::Integer, room::Integer=_GEO_MAXVRT) = getfov(default_context(), instid, room)

"""
    getfvn([ctx,] inst, room=10000) -> (shape, frame, bsight, bounds)

Return the field of view of the instrument with the name `inst`, as `getfvn_c` of
CSPICE. The outputs and the errors are those of [`getfov`](@ref).

`inst` is a name that the loaded kernels map to an ID, or the ID written as a string.

Throws `SpiceError(:idcodenotfound)` if the name is not known.

See also [`getfov`](@ref), [`bodn2c`](@ref).
"""
getfvn(ks::KernelSet, inst::AbstractString, room::Integer=_GEO_MAXVRT) =
    _geo_getfov(ks, _geo_body(ks, inst, "instrument"), room)
getfvn(ctx::SpiceContext, inst::AbstractString, room::Integer=_GEO_MAXVRT) =
    getfvn(kernelset(ctx), inst, room)
getfvn(inst::AbstractString, room::Integer=_GEO_MAXVRT) = getfvn(default_context(), inst, room)

# --- axis of a polygonal field of view --------------------------------------------------

"""
    _geo_hullax(bounds) -> Vec3

`zzhullax`: an axis for the polygonal field of view `bounds`, such that every boundary
vector is less than 90 degrees from it. Throws `SpiceError` (`:invalidcount`,
`:degeneratecase`, `:facenotfound`, `:notsupported`, `:fovtoowide`).
"""
function _geo_hullax(bounds::Vector{Vec3{Float64}})
    n = length(bounds)
    n < 3 && _geo_fov_error(:invalidcount, "a polygonal field of view needs 3 boundary vectors; it has $n")
    halfpi = _M_HALFPI
    found = false
    xidx = 0
    next = 0
    cp = zero(Vec3{Float64})
    # The side of the pyramid that has all other boundary vectors on one side.
    @inline function face(i, nx, cp0)
        c = cp0
        pass1 = true
        ok = true
        m = 1
        # invariant: the vectors 1 .. m-1 are on the inner side of the plane
        while m <= n && ok
            if m != i && m != nx
                sep = vsep(c, bounds[m])
                if pass1
                    if sep > halfpi
                        c = -1.0 * c
                        sep = Float64(pi) - sep
                    end
                    pass1 = false
                end
                ok = sep < halfpi - 1e-12
            end
            ok && (m += 1)
        end
        return ok, c
    end
    i = 1
    while i <= n && !found
        next = i == n ? 1 : i + 1
        cp = _m_cross(bounds[i], bounds[next])
        _m_iszero(cp) && _geo_fov_error(:degeneratecase,
            "the boundary vectors $i and $next of the field of view are parallel")
        ok, cp = face(i, next, cp)
        ok ? (xidx = i; found = true) : (i += 1)
    end
    if !found
        i = 1
        while i <= n && !found
            next = i + 1
            while next <= n && !found
                cp = _m_cross(bounds[i], bounds[next])
                if !_m_iszero(cp)
                    ok, cp = face(i, next, cp)
                    ok && (xidx = i; found = true)
                end
                found || (next += 1)
            end
            found || (i += 1)
        end
    end
    found || _geo_fov_error(:facenotfound,
        "no side of the field of view has all boundary vectors on one side")
    xvec = _m_hat(0.5 * _m_hat(bounds[xidx]) + 0.5 * _m_hat(bounds[next]))
    yvec = _m_hat(cp)
    zvec = _m_ucrss(xvec, yvec)
    minlon = Float64(pi)
    maxlon = 0.0
    for k in 1:n                         # invariant: the extreme longitudes of 1 .. k-1
        (k == xidx || k == next) && continue
        b = bounds[k]
        _, lon, _ = reclat(Vec3{Float64}(_m_dot(xvec, b), _m_dot(yvec, b), _m_dot(zvec, b)))
        lon < minlon && (minlon = lon)
        lon > maxlon && (maxlon = lon)
    end
    minlon < 2e-12 && _geo_fov_error(:notsupported,
        "the field of view is not convex enough: the smallest longitude is $minlon")
    maxlon > Float64(pi) - 2e-12 && _geo_fov_error(:fovtoowide,
        "the field of view is too wide: the largest longitude is $maxlon")
    delta = Float64(pi) - maxlon
    axis = _m_hat(_m_vrotv(yvec, zvec, -delta / 2))
    for k in 1:n
        vsep(bounds[k], axis) > halfpi - 1e-12 && _geo_fov_error(:fovtoowide,
            "the boundary vector $k is 90 degrees or more from the axis of the field of view")
    end
    return axis
end

"`zzfovaxi`: the unit axis of the polygonal field of view `bounds`."
function _geo_fovaxi(bounds::Vector{Vec3{Float64}})
    n = length(bounds)
    n < 3 && _geo_fov_error(:invalidcount, "a polygonal field of view needs 3 boundary vectors; it has $n")
    for i in 1:n
        _m_iszero(_m_cross(bounds[i], bounds[i == n ? 1 : i + 1])) && _geo_fov_error(
            :degeneratecase, "the boundary vectors $i and $(i == n ? 1 : i + 1) are parallel")
    end
    axis = zero(Vec3{Float64})
    for i in 1:n                         # invariant: axis is the sum of the unit vectors 1 .. i-1
        axis = _m_hat(bounds[i]) + axis
    end
    axis = (1.0 / n) * axis
    limit = _M_HALFPI - 1e-12
    for i in 1:n
        if vsep(bounds[i], axis) > limit
            axis = _geo_hullax(bounds)
            break
        end
    end
    return _m_hat(axis)
end

"`frame_`: the unit vector along `x0` and two unit vectors that complete a right-handed frame."
@inline function _geo_frame(x0::Vec3{Float64})
    x = _m_hat(x0)
    a = x[1] * x[1]
    b = x[2] * x[2]
    c = x[3] * x[3]
    a + b + c == 0.0 && return Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 1.0)
    if a <= b && a <= c
        f = sqrt(b + c)
        y = Vec3(0.0, -x[3] / f, x[2] / f)
        return x, y, Vec3(f, -x[1] * y[3], x[1] * y[2])
    elseif b <= a && b <= c
        f = sqrt(a + c)
        y = Vec3(x[3] / f, 0.0, -x[1] / f)
        return x, y, Vec3(x[2] * y[3], f, -x[2] * y[1])
    end
    f = sqrt(a + b)
    y = Vec3(-x[2] / f, x[1] / f, 0.0)
    return x, y, Vec3(-x[3] * y[2], x[3] * y[1], f)
end

# --- winding numbers and the ellipse against the pyramid ----------------------------------

"`zzwind2d`: the winding number of the planar polygon `poly` about `point`."
function _geo_wind2d(poly::Vector{SVector{2,Float64}}, point::SVector{2,Float64})
    n = length(poly)
    n < 3 && _geo_fov_error(:degeneratecase, "a polygon needs 3 vertices; it has $n")
    atotal = 0.0
    rvec = poly[1] - point
    for i in 2:n+1                       # invariant: atotal is the angle swept by the sides 1 .. i-2
        rnext = poly[i <= n ? i : 1] - point
        sep = vsepg(rnext, rvec)
        atotal += (rnext[1] * -rvec[2] + rnext[2] * rvec[1]) >= 0.0 ? sep : -sep
        rvec = rnext
    end
    return round(Int, atotal / _M_TWOPI, RoundNearestTiesAway)
end

"""
    _geo_wind_rays(plane, bounds, point) -> Int

`zzwind` for the polygon whose vertices are the intercepts with `plane` of the rays
from the origin along `bounds`: the winding number about `point`.
`SpiceError(:bug)` if a ray has no single intercept, as in `zzelvupy`.
"""
function _geo_wind_rays(plane::SpicePlane, bounds::Vector{Vec3{Float64}}, point::Vec3{Float64})
    n = length(bounds)
    origin = zero(Vec3{Float64})
    @inline function vertex(i)
        nx, v = inrypl(origin, bounds[i], plane)
        nx == 1 || _geo_fov_error(:bug, "the boundary vector $i does not meet the plane " *
            "normal to the axis of the field of view, as SPICE(BUG) of CSPICE")
        return v
    end
    for i in 1:n
        vertex(i)
    end
    normal, _ = pl2nvc(plane)
    _m_iszero(normal) && _geo_fov_error(:zerovector, "the normal of the plane is the zero vector")
    v1 = vertex(1)
    _m_dot(normal, v1) < 0.0 && (normal = -normal)
    rvec = _m_perp(v1 - point, normal, true)
    atotal = 0.0
    for i in 2:n+1                       # invariant: atotal is the angle swept by the sides 1 .. i-2
        rnext = _m_perp(vertex(i <= n ? i : 1) - point, normal, true)
        sep = vsep(rnext, rvec)
        atotal += _m_dot(rnext, _m_ucrss(normal, rvec)) >= 0.0 ? sep : -sep
        rvec = rnext
    end
    return round(Int, atotal / _M_TWOPI, RoundNearestTiesAway)
end

"""
    _geo_elvupy(ellips, axis, bounds) -> Bool

`zzelvupy` with the vertex at the origin: true if the ellipse `ellips` meets the
pyramid whose edges are the rays from the origin along `bounds` and whose axis is
`axis`. Throws `SpiceError` (`:zerovector`, `:degeneratecase`, `:invalidcount`,
`:invalidfov`, `:bug`).
"""
function _geo_elvupy(ellips::SpiceEllipse{Float64}, axis::Vec3{Float64},
        bounds::Vector{Vec3{Float64}})
    n = length(bounds)
    center, gv1, gv2 = el2cgv(ellips)
    smajor, sminor = _m_saelgv(gv1, gv2)
    (_m_iszero(smajor) || _m_iszero(sminor)) && _geo_fov_error(:zerovector,
        "the limb of the target is a degenerate ellipse")
    scale = 1.0 / max(max(_m_norm(center), _m_norm(smajor)), 0.0)
    center = scale * center
    smajor = scale * smajor
    sminor = scale * sminor
    apex = zero(Vec3{Float64})
    ellscl = _m_cgv2el(center, smajor, sminor)
    a = _m_norm(smajor)
    b = _m_norm(sminor)
    b == 0.0 && _geo_fov_error(:degeneratecase, "the limb of the target is a degenerate ellipse")
    _m_iszero(axis) && _geo_fov_error(:zerovector, "the axis of the field of view is the zero vector")
    (n < 3 || n > _GEO_MAXVRT) && _geo_fov_error(:invalidcount,
        "the field of view has $n boundary vectors; the range is 3:$_GEO_MAXVRT")
    ctrvec = center - apex
    pasize = 0.0
    for i in 1:n                         # invariant: pasize is the largest angle of 1 .. i-1
        asep = vsep(axis, bounds[i])
        asep >= Float64(pi) / 2 && _geo_fov_error(:invalidfov,
            "the boundary vector $i is 90 degrees or more from the axis of the field of view")
        pasize = max(pasize, asep)
        j = i < n ? i + 1 : 1
        if _m_iszero(_m_ucrss(bounds[i], bounds[j]))
            (_m_iszero(bounds[i]) || _m_iszero(bounds[j])) && _geo_fov_error(:zerovector,
                "a boundary vector of the field of view is the zero vector")
            _geo_fov_error(:invalidfov, "the boundary vectors $i and $j are parallel")
        end
    end
    d = _m_norm(center - apex)
    if a < d
        # The cone about the ellipse against the cone about the pyramid.
        easize = asin(a / d)
        vsep(ctrvec, axis) - (easize + pasize) > 0.0 && return false
    end
    if _m_dot(axis, ctrvec) > 0.0
        # The center of the ellipse inside of the pyramid.
        _geo_wind_rays(nvp2pl(axis, ctrvec), bounds, ctrvec) != 0 && return true
    end
    # The first edge of the pyramid through the plane region that the ellipse bounds.
    nxpts, xpt = inrypl(apex, bounds[1], psv2pl(center, smajor, sminor))
    if nxpts == 1
        offset = xpt - center
        x = _m_dot(offset, smajor) / a
        y = _m_dot(offset, sminor) / b
        (x / a)^2 + (y / b)^2 <= 1.0 && return true
    end
    # The ellipse through a side of the pyramid.
    for i in 1:n                         # invariant: the ellipse does not meet the sides 1 .. i-1
        j = i < n ? i + 1 : 1
        nx, xpt1, xpt2 = inelpl(ellscl, psv2pl(apex, bounds[i], bounds[j]))
        nx > 0 || continue
        edge1 = _m_hat(bounds[i])
        edge2 = _m_hat(bounds[j])
        ebsctr = 0.5 * edge1 + 0.5 * edge2
        hafedg = vsep(edge1, edge2) / 2.0
        vxpt1 = _m_hat(xpt1 - apex)
        vxpt2 = _m_hat(xpt2 - apex)
        vbsctr = 0.5 * vxpt1 + 0.5 * vxpt2
        hafsec = vsep(vxpt1, vxpt2) / 2.0
        (vsep(vxpt1, ebsctr) <= hafedg || vsep(vxpt2, ebsctr) <= hafedg ||
         vsep(vbsctr, ebsctr) <= hafedg || vsep(ebsctr, vbsctr) <= hafsec) && return true
    end
    return false
end

# --- visibility in a field of view ------------------------------------------------------

const _GEO_FOV_POINT = 1
const _GEO_FOV_ELLIPSOID = 2
const _GEO_FOV_RAY = 3

"The correction of a ray: `(stellar, transmit)`; a correction with light time is not valid."
function _geo_fov_raycorr(abcorr::AbstractString)
    _dsk_eqstr(abcorr, "NONE") && return false, false
    _dsk_eqstr(abcorr, "S") && return true, false
    _dsk_eqstr(abcorr, "XS") && return true, true
    _geo_fov_error(:invalidoption, "the aberration correction <$abcorr> is not valid " *
        "for a ray: use NONE, S or XS")
end
_geo_fov_raycorr(c::Aberration) = c.uselt ? _geo_fov_error(:invalidoption,
    "a correction with light time is not valid for a ray") : (c.stellar, c.transmit)

function _geo_fov_tshape(tshape::AbstractString)
    s = uppercase(rstrip(lstrip(tshape, ' '), ' '))
    s == "POINT" && return _GEO_FOV_POINT
    s == "ELLIPSOID" && return _GEO_FOV_ELLIPSOID
    s == "RAY" && return _GEO_FOV_RAY
    _geo_fov_error(:invalidshape, "the target shape <$tshape> is not valid: use POINT or ELLIPSOID")
end

"""
    _geo_fovvis(ks, inst, tshape, raydir, target, tframe, abcorr, obsrvr, et) -> Bool

`zzgffvin` followed by `zzgffvst`: the visibility of the target (a point, an ellipsoid
or a ray, by the code `tshape`) in the field of view of `inst` at the epoch `et`. The
checks are made in the order of CSPICE.
"""
function _geo_fovvis(ks::KernelSet, inst::NameOrID, tshape::Int, raydir::Vec3{Float64},
        target::NameOrID, tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, et0::Real)
    obs = _geo_body(ks, obsrvr, "observer")
    uray = tshape == _GEO_FOV_RAY
    etrg = tshape == _GEO_FOV_ELLIPSOID
    targ = Int32(0)
    if !uray
        targ = _geo_body(ks, target, "target")
        targ == obs && _geo_same_bodies("field of view test", obsrvr)
    end
    tfid = Int32(0)
    if uray || etrg
        (tframe isa AbstractString && all(==(' '), tframe)) && _geo_fov_error(:invalidframe,
            "the name of the reference frame is blank")
        tfid = _geo_namfrm(ks, tframe)
        tfid == Int32(0) && _geo_fov_error(:invalidframe,
            "the reference frame $tframe is not recognized")
        center, _, found = _geo_frinfo(ks, tfid)
        found || _geo_fov_error(:frameinfonotfound,
            "no data are available for the reference frame $tframe")
        (etrg && center != targ) && _geo_bad_center(tframe, target, center)
    end
    ustl = false
    xmit = false
    corr = _aberration(0x00)
    if uray
        ustl, xmit = _geo_fov_raycorr(abcorr)
    else
        corr = _geo_abcorr(abcorr)
    end
    trad = etrg ? _geo_radii(ks, targ) : (0.0, 0.0, 0.0)
    rdir = zero(Vec3{Float64})
    if uray
        _m_iszero(raydir) && _geo_fov_error(:zerovector, "the direction of the ray is the zero vector")
        rdir = _m_hat(raydir)
    end
    instid = _geo_body(ks, inst, "instrument")
    shape, ifrname, bsite, bounds = _geo_getfov(ks, instid, _GEO_MAXVRT)
    ifid = _geo_namfrm(ks, ifrname)
    poly = shape == "POLYGON" || shape == "RECTANGLE"
    circle = shape == "CIRCLE"
    faxi = poly ? _geo_fovaxi(bounds) : _m_hat(bsite)
    arad = 0.0
    for b in bounds
        arad = max(arad, vsep(b, faxi))
    end
    arad > _M_HALFPI - 1e-6 && _geo_fov_error(:fovtoowide,
        "the field of view of the instrument $inst has the angular radius $(rad2deg(arad)) " *
        "degrees; the limit is $(rad2deg(_M_HALFPI - 1e-6)) degrees")
    plane = nvc2pl(faxi, 1.0)
    fvct = _m_hat(faxi)
    origin = zero(Vec3{Float64})
    semi1 = zero(Vec3{Float64})
    semi2 = zero(Vec3{Float64})
    xmag1 = 0.0
    xmag2 = 0.0
    if shape == "ELLIPSE" || (circle && etrg)
        b2 = circle ? _m_vrotv(bounds[1], faxi, _M_HALFPI) : bounds[2]
        nx, p = inrypl(origin, bounds[1], plane)
        nx == 1 || _geo_fov_error(:degeneratecase, "the first boundary vector does not meet the plane of the field of view")
        semi1 = p - fvct
        xmag1 = _m_norm(semi1)
        xmag1 == 0.0 && _geo_fov_error(:degeneratecase, "the field of view has a zero semi-axis")
        nx, p = inrypl(origin, b2, plane)
        nx == 1 || _geo_fov_error(:degeneratecase, "the second boundary vector does not meet the plane of the field of view")
        semi2 = p - fvct
        xmag2 = _m_norm(semi2)
        xmag2 == 0.0 && _geo_fov_error(:degeneratecase, "the field of view has a zero semi-axis")
    end
    et = Float64(check_epoch(et0))
    ifid == Int32(0) && _geo_unknown_frame(ifrname)
    if etrg
        pos, lt = _geo_spkezp(ks, targ, et, tfid, corr, obs)
        ettarg = et + _geo_sign(corr) * lt
        insmat = _m_mxm(frame_rotation(ks, Int32(1), ifid, et), frame_rotation(ks, tfid, Int32(1), ettarg))
        if poly
            limb = edlimb(trad[1], trad[2], trad[3], -pos)
            lc, s1, s2 = el2cgv(limb)
            fvlimb = _m_cgv2el(_m_mxv(insmat, pos + lc), _m_mxv(insmat, s1), _m_mxv(insmat, s2))
            return _geo_elvupy(fvlimb, faxi, bounds)
        end
        # The field of view as a small ellipsoid, for the occultation test.
        edct = (xmag1 * xmag1 + 1.0) * fvct
        r3 = xmag1 * sqrt(xmag1^2 + 1.0)
        c3 = r3 * fvct
        c1 = r3 * _m_hat(semi1)
        r2 = xmag2 / xmag1 * r3
        c2 = r2 * _m_ucrss(c3, c1)
        escale = 1e-5 / (_m_norm(edct) + max(r3, r2))
        fsmx = hcat(escale * c1, escale * c2, escale * c3)
        trgsmx = hcat(trad[1] * insmat[:, 1], trad[2] * insmat[:, 2], trad[3] * insmat[:, 3])
        return _geo_occed(origin, escale * edct, fsmx, _m_mxv(insmat, pos), trgsmx) != 0
    end
    if uray
        dir = _m_mxv(frame_rotation(ks, tfid, ifid, et), rdir)
        if ustl
            _, vobs = _ssb(ks, obs, et, Int32(1))
            m1 = frame_rotation(ks, ifid, Int32(1), et)
            dir = _m_mtxv(m1, _stelab(_m_mxv(m1, dir), xmit ? -vobs : vobs))
        end
    else
        v, _ = _geo_spkezp(ks, targ, et, Int32(1), corr, obs)
        dir = _m_mxv(frame_rotation(ks, Int32(1), ifid, et), v)
    end
    sep = vsep(dir, faxi)
    circle && return sep <= vsep(faxi, bounds[1])
    sep > arad && return false
    nx, xpt = inrypl(origin, dir, plane)
    nx == 0 && return false
    nx == 1 || _geo_fov_error(:bug, "the ray is in the plane of the field of view")
    fovpt = xpt - fvct
    if !poly
        c1 = _m_dot(fovpt, semi1) / xmag1
        c2 = _m_dot(fovpt, semi2) / xmag2
        return (c1 / xmag1)^2 + (c2 / xmag2)^2 <= 1.0
    end
    z, x, y = _geo_frame(faxi)
    fpol = Vector{SVector{2,Float64}}(undef, length(bounds))
    for i in eachindex(bounds)           # invariant: fpol[1:i-1] are the vertices in the plane
        nx, p = inrypl(origin, bounds[i], plane)
        nx == 1 || _geo_fov_error(:degeneratecase, "the boundary vector $i does not meet the plane of the field of view")
        q = p - fvct
        fpol[i] = SVector(_m_dot(x, q), _m_dot(y, q))
    end
    return _geo_wind2d(fpol, SVector(_m_dot(x, fovpt), _m_dot(y, fovpt))) != 0
end

"""
    fovtrg([ctx,] inst, target, tshape, tframe, abcorr, obsrvr, et) -> Bool

Return `true` if the body `target` is in the field of view of the instrument `inst` at
`et` (seconds past J2000 TDB), as seen by `obsrvr`, as `fovtrg_c` of CSPICE.

`tshape` is the model of the target:

| `tshape` | Test | `tframe` |
|---|---|---|
| `"ELLIPSOID"` | any part of the reference ellipsoid is in the field of view | body-fixed frame of the target |
| `"POINT"` | the target center is in the field of view | not used; pass `" "` |

`inst` is an instrument name or ID with a field of view in the kernel pool (see
[`getfov`](@ref)). `obsrvr` is the body that carries the instrument. The instrument
frame is evaluated at `et`, so a CK and an SCLK kernel for the spacecraft are needed
if that frame is a CK frame. `abcorr` is a correction string of [`spkezr`](@ref) or an
[`Aberration`](@ref). `et` cannot be a `ForwardDiff.Dual`, because the result is
discrete.

Throws `SpiceError` with the kind
- `:invalidoption` if `tshape` is `"RAY"` (use [`fovray`](@ref)) or if `abcorr` is not a
  valid correction;
- `:idcodenotfound` if a body or instrument name is not known;
- `:invalidshape` if `tshape` is not `"ELLIPSOID"` or `"POINT"`;
- `:bodiesnotdistinct` if `target` and `obsrvr` are the same body;
- `:invalidframe` if `tframe` is blank, not known, or not centered on an ellipsoid
  target; `:frameinfonotfound` if the frame has no definition loaded;
- `:unknownframe` if the instrument frame is not known;
- `:fovtoowide` if the angular radius of the field of view about its axis is not below
  90 degrees minus 1e-6 rad;
- for a rectangle or a polygon: `:degeneratecase` if two neighbor boundary vectors are
  parallel, `:facenotfound` if the boundary vectors are not in one half space,
  `:notsupported` if a boundary vector is within 2e-12 rad of the plane of a side of
  the field of view, `:invalidfov` if a vector is 90 degrees or more from the axis;
- for a circle or an ellipse: `:degeneratecase` if a semi-axis is zero;
- `:notdisjoint` if the observer is inside the target ellipsoid.

Also throws the errors of [`getfov`](@ref) and [`spkpos`](@ref), and the kinds of
[`srfrec`](@ref) if the radii of an ellipsoid target are missing or not valid.

# Examples
```julia
# Needs an SPK, a PCK, a leap seconds kernel, and for the spacecraft an instrument
# kernel, a frame kernel, a CK and an SCLK kernel. The names are those of Cassini.
et = str2et("2004 JUN 11 19:32:00 UTC")
fovtrg("CASSINI_ISS_NAC", "PHOEBE", "ELLIPSOID", "IAU_PHOEBE", "LT+S", "CASSINI", et)
fovtrg("CASSINI_ISS_NAC", "SUN", "POINT", " ", "LT+S", "CASSINI", et)
```

See also [`fovray`](@ref), [`getfov`](@ref), [`gftfov`](@ref), [`occult`](@ref).
"""
function fovtrg(ks::KernelSet, inst::NameOrID, target::NameOrID, tshape::AbstractString,
        tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, et::Real)
    _dsk_eqstr(tshape, "RAY") && _geo_fov_error(:invalidoption,
        "the target shape RAY is not valid in fovtrg; use fovray")
    # The observer is checked before the shape, as in CSPICE.
    _geo_body(ks, obsrvr, "observer")
    return _geo_fovvis(ks, inst, _geo_fov_tshape(tshape), zero(Vec3{Float64}), target, tframe,
        abcorr, obsrvr, et)
end
fovtrg(ctx::SpiceContext, inst::NameOrID, target::NameOrID, tshape::AbstractString,
    tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, et::Real) =
    fovtrg(kernelset(ctx), inst, target, tshape, tframe, abcorr, obsrvr, et)
fovtrg(inst::NameOrID, target::NameOrID, tshape::AbstractString, tframe::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, et::Real) =
    fovtrg(default_context(), inst, target, tshape, tframe, abcorr, obsrvr, et)

"""
    fovray([ctx,] inst, raydir, rframe, abcorr, obsrvr, et) -> Bool

Return `true` if the ray along `raydir` is in the field of view of the instrument `inst`
at `et` (seconds past J2000 TDB), as `fovray_c` of CSPICE.

`raydir` is a direction (3 elements, not zero) in the frame `rframe`, for example the
direction of a star in `"J2000"`. `obsrvr` is the body that carries the instrument; its
velocity gives the stellar aberration. `abcorr` is one of:

| `abcorr` | Correction of the direction |
|---|---|
| `"NONE"` | none |
| `"S"` | stellar aberration, for light that the observer receives |
| `"XS"` | stellar aberration, for light that the observer sends |

An [`Aberration`](@ref) without light time is also valid. `et` cannot be a
`ForwardDiff.Dual`. The notes of [`fovtrg`](@ref) on the instrument apply.

Throws `SpiceError` with the kind
- `:zerovector` if `raydir` is zero, `:badarg` if it does not have 3 elements;
- `:idcodenotfound` if the observer or the instrument name is not known;
- `:invalidframe` if `rframe` is blank or not known, `:frameinfonotfound` if it has no
  definition loaded;
- `:invalidoption` if `abcorr` is another value, such as a correction with light time.

Also throws the errors of [`fovtrg`](@ref) for the field of view.

See also [`fovtrg`](@ref), [`gfrfov`](@ref), [`sincpt`](@ref).
"""
fovray(ks::KernelSet, inst::NameOrID, raydir::AbstractVector, rframe::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, et::Real) =
    _geo_fovvis(ks, inst, _GEO_FOV_RAY, _geo_pos(raydir), 0, rframe, abcorr, obsrvr, et)
fovray(ctx::SpiceContext, inst::NameOrID, raydir::AbstractVector, rframe::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, et::Real) =
    fovray(kernelset(ctx), inst, raydir, rframe, abcorr, obsrvr, et)
fovray(inst::NameOrID, raydir::AbstractVector, rframe::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, et::Real) =
    fovray(default_context(), inst, raydir, rframe, abcorr, obsrvr, et)
