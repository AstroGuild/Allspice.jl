# Ellipsoid, plane and ellipse geometry of CSPICE.
#
# A triaxial ellipsoid is given by its semi-axis lengths `a`, `b`, `c` along x, y, z and
# is centered at the origin. Planes and ellipses are the small immutable structs
# `SpicePlane` and `SpiceEllipse`, which have the invariants of the CSPICE structures of
# the same names. A call that has a "found" flag in CSPICE returns `nothing` when the
# flag is false.

"""
    SpicePlane{T}

Plane of the points `x` with `normal ⋅ x == constant`, as `SpicePlane` of CSPICE.

The fields are `normal::SVector{3,T}` and `constant::T`. A plane from
[`nvc2pl`](@ref), [`nvp2pl`](@ref) or [`psv2pl`](@ref) has a unit `normal` and
`constant >= 0`. `constant` is then the distance of the plane from the origin. The
constructor does not check this.

See also [`pl2nvc`](@ref), [`pl2nvp`](@ref), [`pl2psv`](@ref).
"""
struct SpicePlane{T<:Real}
    normal::SVector{3,T}
    constant::T
end
SpicePlane(normal::SVector{3}, constant::Real) =
    SpicePlane{promote_type(eltype(normal), typeof(constant))}(normal, constant)

"""
    SpiceEllipse{T}

Ellipse of the points `center + cos(θ) * semimajor + sin(θ) * semiminor`, as
`SpiceEllipse` of CSPICE.

The fields `center`, `semimajor` and `semiminor` are `SVector{3,T}` values. An ellipse
from [`cgv2el`](@ref) has orthogonal semi-axes, and `semimajor` is not shorter than
`semiminor`. The constructor does not check this.

See also [`el2cgv`](@ref).
"""
struct SpiceEllipse{T<:Real}
    center::SVector{3,T}
    semimajor::SVector{3,T}
    semiminor::SVector{3,T}
end
function SpiceEllipse(center::SVector{3}, semimajor::SVector{3}, semiminor::SVector{3})
    T = promote_type(eltype(center), eltype(semimajor), eltype(semiminor))
    return SpiceEllipse{T}(center, semimajor, semiminor)
end

@inline _m_abs(x::Real) = primal(x) < 0.0 ? -x : x

@inline _m_max(x, y) = primal(y) > primal(x) ? y : x

@noinline _m_badaxis(a, b, c) = throw(SpiceError(:badaxislength,
    "the ellipsoid semi-axis lengths must be positive; they are $(primal(a)), $(primal(b)), $(primal(c))"))

@inline _m_axes_ok(a, b, c) = primal(a) > 0.0 && primal(b) > 0.0 && primal(c) > 0.0

# --- surface normal, ray intercept, scaling to the surface ------------------------------

@inline function _m_surfnm(a, b, c, point::SVector{3})
    _m_axes_ok(a, b, c) || _m_badaxis(a, b, c)
    m = primal(a) <= primal(b) ? a : b
    m = primal(m) <= primal(c) ? m : c
    a1 = m / a
    b1 = m / b
    c1 = m / c
    return _m_hat(SVector(point[1] * (a1 * a1), point[2] * (b1 * b1), point[3] * (c1 * c1)))
end

"""
    surfnm(a, b, c, point) -> SVector{3}

Return the outward unit normal of an ellipsoid at the surface point `point`, as
`surfnm_c` of CSPICE.

The ellipsoid has its center at the origin and the semi-axes `a`, `b`, `c` along x, y,
z. This holds for all the ellipsoid calls. `point` is not checked to be on the surface.

Throws `SpiceError(:badaxislength)` if a semi-axis is not positive, and
`SpiceError(:badarg)` if `point` does not have 3 elements.

# Examples
```jldoctest
julia> surfnm(1, 2, 3, [0, 0, 3]) == [0, 0, 1]
true
```

See also [`ednmpt`](@ref), [`nearpt`](@ref).
"""
surfnm(a::Real, b::Real, c::Real, point::AbstractVector) =
    _m_surfnm(_m_f(a), _m_f(b), _m_f(c), _m_v3(point))

"`(found, point)` of `surfpt_`; the point is zero when nothing is found."
@inline function _m_surfpt(positn::SVector{3}, u::SVector{3}, a, b, c)
    _m_iszero(u) && throw(SpiceError(:zerovector, "the ray direction is the zero vector"))
    _m_axes_ok(a, b, c) || _m_badaxis(a, b, c)
    x = SVector(u[1] / a, u[2] / b, u[3] / c)
    y = SVector(positn[1] / a, positn[2] / b, positn[3] / c)
    # In the scaled space the ellipsoid is the unit sphere; p is the point of the
    # line that is nearest to the center.
    p = _m_perp(y, x, false)
    yproj = y - p
    ymag = _m_norm(y)
    pmag = _m_norm(p)
    ux = _m_hat(x)
    none = zero(p)
    sign = 1.0
    if primal(ymag) > 1.0
        primal(pmag) > 1.0 && return false, none
        primal(_m_dot(yproj, x)) > 0.0 && return false, none
        primal(pmag) == 1.0 && return true, SVector(p[1] * a, p[2] * b, p[3] * c)
        sign = -1.0
    elseif primal(ymag) == 1.0
        return true, positn + none
    end
    t = 1.0 - pmag * pmag
    scale = primal(t) > 0.0 ? sqrt(t) : zero(t)
    q = 1.0 * p + (sign * scale) * ux
    return true, SVector(q[1] * a, q[2] * b, q[3] * c)
end

"""
    surfpt(positn, u, a, b, c) -> Union{SVector{3},Nothing}

Return the first point where a ray meets an ellipsoid, or `nothing` if the ray misses,
as `surfpt_c` of CSPICE.

The ray starts at `positn` and has the direction `u`. The ellipsoid has the semi-axes
`a`, `b`, `c`. CSPICE returns a `found` flag in place of `nothing`. A start point
inside the ellipsoid gives the point where the ray leaves. A start point on the surface
gives that point.

Throws `SpiceError(:zerovector)` if `u` is zero, `SpiceError(:badaxislength)` if a
semi-axis is not positive, and `SpiceError(:badarg)` if a vector does not have 3
elements.

# Examples
```jldoctest
julia> surfpt([10, 0, 0], [-1, 0, 0], 1, 2, 3) == [1, 0, 0]
true

julia> isnothing(surfpt([10, 0, 0], [1, 0, 0], 1, 2, 3))
true
```

See also [`surfpv`](@ref), [`npedln`](@ref), [`nearpt`](@ref).
"""
function surfpt(positn::AbstractVector, u::AbstractVector, a::Real, b::Real, c::Real)
    found, point = _m_surfpt(_m_v3(positn), _m_v3(u), _m_f(a), _m_f(b), _m_f(c))
    return found ? point : nothing
end

"""
    edpnt(p, a, b, c) -> SVector{3}

Return the point of an ellipsoid on the ray from the center through the point `p`, as
`edpnt_c` of CSPICE. The ellipsoid has the semi-axes `a`, `b`, `c`.

Throws `SpiceError(:invalidradii)` if a semi-axis is not positive, and
`SpiceError(:zerovector)` if `p` is zero. Throws `SpiceError(:pointtoosmall)` if `p` is
so near the center that its scaled square is zero.

See also [`nearpt`](@ref), [`surfpt`](@ref).
"""
function edpnt(p::AbstractVector, a::Real, b::Real, c::Real)
    v = _m_v3(p)
    a, b, c = _m_f(a), _m_f(b), _m_f(c)
    _m_axes_ok(a, b, c) || throw(SpiceError(:invalidradii,
        "the ellipsoid radii must be positive; they are $(primal(a)), $(primal(b)), $(primal(c))"))
    _m_iszero(v) && throw(SpiceError(:zerovector, "the input point is the zero vector"))
    x = v[1] / a
    y = v[2] / b
    z = v[3] / c
    level = x * x + y * y + z * z
    primal(level) > 0.0 || throw(SpiceError(:pointtoosmall,
        "the input point is too close to the origin to be scaled to the ellipsoid"))
    sq = sqrt(level)
    return SVector(v[1] / sq, v[2] / sq, v[3] / sq)
end

"""
    ednmpt(a, b, c, normal) -> SVector{3}

Return the point of an ellipsoid where the outward normal is parallel to `normal`, as
`ednmpt_c` of CSPICE. The ellipsoid has the semi-axes `a`, `b`, `c`.

Throws:
- `SpiceError(:badaxislength)` if a semi-axis is not positive.
- `SpiceError(:axisunderflow)` if the ratio of a semi-axis to the largest one
  underflows.
- `SpiceError(:zerovector)` if `normal` is zero.
- `SpiceError(:degeneratecase)` if the scale factor underflows.

# Examples
```jldoctest
julia> ednmpt(1, 2, 3, [0, 1, 0]) == [0, 2, 0]
true
```

See also [`surfnm`](@ref).
"""
function ednmpt(a::Real, b::Real, c::Real, normal::AbstractVector)
    n = _m_v3(normal)
    a, b, c = _m_f(a), _m_f(b), _m_f(c)
    _m_axes_ok(a, b, c) || _m_badaxis(a, b, c)
    scale = _m_max(_m_max(a, b), c)
    sa = a / scale
    sb = b / scale
    sc = c / scale
    _m_axes_ok(sa, sb, sc) || throw(SpiceError(:axisunderflow,
        "a scaled semi-axis length is not positive"))
    _m_iszero(n) && throw(SpiceError(:zerovector, "the input normal vector is the zero vector"))
    na2 = n[1] * sa * sa
    nb2 = n[2] * sb * sb
    nc2 = n[3] * sc * sc
    arg = na2 * n[1] + nb2 * n[2] + nc2 * n[3]
    primal(arg) > 0.0 || throw(SpiceError(:degeneratecase,
        "the reciprocal of the square of the scale factor is not positive"))
    lambda = (1.0 / sqrt(arg)) * scale
    return SVector(lambda * na2, lambda * nb2, lambda * nc2)
end

# --- nearest point on an ellipsoid ------------------------------------------------------

@noinline _m_nearpt_bug(what::String) = throw(SpiceError(:bug, "nearpt: " * what))

"A coordinate that is negligible against its semi-axis becomes zero, as in `nearpt_`."
@inline _m_nearpt_trim(p::Float64, ax::Float64) =
    (p * 0.5 + ax == ax || p * 0.5 - ax == -ax) ? 0.0 : p

"Lower limit of `lambda` that one coordinate of an interior point gives."
@inline _m_nearpt_low(p::Float64, ax::Float64, top::Float64) =
    p == 0.0 ? -top : ax * (abs(p) * 0.5 - ax)

"One coordinate of the point of the ellipsoid that belongs to `lambda`."
@inline function _m_nearpt_foot(p::Float64, orig::Float64, axsq::Float64, lambda::Float64)
    p == 0.0 && return 0.0
    denom = lambda / axsq + 1.0
    denom <= 0.0 && _m_nearpt_bug("a denominator is not positive")
    return orig / denom
end

"Squares of `point[i] / denom`, with the cut of `nearpt_` for an interior point."
@inline function _m_nearpt_terms(point::NTuple{3,Float64}, axis::NTuple{3,Float64},
        lambda::Float64, inside::Bool)
    return ntuple(Val(3)) do i
        point[i] == 0.0 && return 0.0
        denom = axis[i] + lambda / axis[i]
        trim = abs(point[i]) * 0.5 > denom
        if inside && trim
            return 4.0
        end
        denom == 0.0 && _m_nearpt_bug("a denominator is zero")
        factor = point[i] / denom
        return factor * factor
    end
end

"""
    _m_nearpt(positn, a, b, c) -> (npoint, alt)

`nearpt_` of CSPICE for `Float64` arguments: the same axis ordering, scaling, bracketing
and bisection of the multiplier `lambda`, and the same (at most five) correction steps
that move the solution along the error of the first one.
"""
function _m_nearpt(positn::SVector{3,Float64}, a::Float64, b::Float64, c::Float64)
    (a > 0.0 && b > 0.0 && c > 0.0) || _m_badaxis(a, b, c)
    ax = (a, b, c)
    # `orderd_` for three elements is an insertion sort: ascending and stable.
    i1, i2, i3 = 1, 2, 3
    ax[i1] <= ax[i2] || ((i1, i2) = (i2, i1))
    if !(ax[i2] <= ax[i3])
        i2, i3 = i3, i2
        ax[i1] <= ax[i2] || ((i1, i2) = (i2, i1))
    end
    iorder = (i1, i2, i3)
    scale = 1.0 / ax[i1]
    axis = (scale * ax[i1], scale * ax[i2], scale * ax[i3])
    point = (scale * positn[i1], scale * positn[i2], scale * positn[i3])
    orignl = SVector(point)
    pnorm = _m_norm(orignl)
    toobig = sqrt(floatmax(Float64) / 100.0)
    (axis[2] > toobig || axis[3] > toobig) && throw(SpiceError(:badaxislength,
        "the ratio of the longest to the shortest semi-axis is too large"))
    for i in 1:3
        sqrt(axis[i]) * sqrt(abs(point[i])) > toobig && throw(SpiceError(:inputstoolarge,
            "the product of a semi-axis length and a point component is too large"))
    end
    axisqr = (axis[1] * axis[1], axis[2] * axis[2], axis[3] * axis[3])
    snglpt = 4
    solutn = 1
    solvng = true
    sign = 1.0
    bestpt = orignl
    bestht = 0.0
    olderr = 0.0
    # invariant: bestpt, bestht hold the best solution of the steps before `solutn`
    while solvng
        copy = point
        point = (_m_nearpt_trim(point[1], axis[1]), _m_nearpt_trim(point[2], axis[2]),
            _m_nearpt_trim(point[3], axis[3]))
        inside = false
        d1 = point[1] / axis[1]
        d2 = point[2] / axis[2]
        d3 = point[3] / axis[3]
        q = d1 * d1 + d2 * d2 + d3 * d3 - 1.0
        solutn == 1 && (sign = q >= 0.0 ? 1.0 : -1.0)
        qlower = 0.0
        qupper = 0.0
        lower = 0.0
        upper = 0.0
        lambda = 0.0
        if q == 0.0
            # on the surface: lambda = 0
        elseif q > 0.0
            lower = max(max(max(0.0, (abs(point[1]) * 0.5 - axis[1]) * axis[1]),
                    (abs(point[2]) * 0.5 - axis[2]) * axis[2]), (abs(point[3]) * 0.5 - axis[3]) * axis[3])
            upper = max(max(abs(axis[1] * point[1]), abs(axis[2] * point[2])),
                abs(axis[3] * point[3])) * 2.0
            lambda = upper
        else
            inside = true
            snglpt = point[1] != 0.0 ? 1 : point[2] != 0.0 ? 2 : point[3] != 0.0 ? 3 : 4
            if snglpt <= 3
                lambda = max(max(_m_nearpt_low(point[1], axis[1], axisqr[3]),
                        _m_nearpt_low(point[2], axis[2], axisqr[3])),
                    _m_nearpt_low(point[3], axis[3], axisqr[3]))
                lower = lambda
                upper = max(lower, 0.0)
            else
                # the center: the nearest point is on the shortest axis
                s = copy[1] < 0.0 ? -axis[1] : axis[1]
                point = (s, point[2], point[3])
                copy = (s, 0.0, 0.0)
                q = 0.0
                inside = false
            end
        end
        term = _m_nearpt_terms(point, axis, lambda, inside)
        if !inside
            qlower = q
            qupper = term[1] + term[2] + term[3] - 1.0
        else
            qupper = q
            qlower = term[1] + term[2] + term[3] - 1.0
        end
        qlower = max(0.0, qlower)
        qupper = min(0.0, qupper)
        lambda = upper
        q = qupper
        itr = 0
        # invariant: the root of q(lambda) is in [lower, upper], q(lower) >= 0 >= q(upper)
        while upper - lower > 0.0
            itr += 1
            itr > 2048 && _m_nearpt_bug("the bisection did not converge")
            lower = min(lower, upper)
            qlower = max(0.0, qlower)
            qupper = min(0.0, qupper)
            if q == 0.0
                lower = lambda
                upper = lambda
            else
                if q < 0.0
                    upper = lambda
                    qupper = q
                else
                    lower = lambda
                    qlower = q
                end
                lambda = lower * 0.5 + upper * 0.5
                lambda = max(lower, min(upper, lambda))
            end
            if abs(lambda - lower) <= 1.0e-16 || abs(lambda - upper) <= 1.0e-16
                if abs(qlower) < abs(qupper)
                    upper = lower
                else
                    lower = upper
                end
            end
            if upper - lower > 0.0
                term = _m_nearpt_terms(point, axis, lambda, inside)
                q = term[1] + term[2] + term[3] - 1.0
            end
        end
        lambda = lower
        spoint = SVector(_m_nearpt_foot(point[1], copy[1], axisqr[1], lambda),
            _m_nearpt_foot(point[2], copy[2], axisqr[2], lambda),
            _m_nearpt_foot(point[3], copy[3], axisqr[3], lambda))
        if inside && (snglpt == 2 || snglpt == 3)
            # A point of a principal plane inside the ellipsoid can have its nearest
            # point off that plane.
            extra = false
            epoint = spoint
            if axis[1] != axis[snglpt] && abs(point[2]) <= axis[2] - axisqr[1] / axis[2] &&
               abs(point[3]) <= axis[3] - axisqr[1] / axis[3]
                if axis[1] == axis[2]
                    denom3 = 1.0 - axisqr[1] / axisqr[3]
                    if denom3 > 0.0
                        e2 = 0.0
                        e3 = point[3] / denom3
                        t1 = e2 / axis[2]
                        t2 = e3 / axis[3]
                        temp = 1.0 - t1 * t1 - t2 * t2
                        if temp > 0.0
                            epoint = SVector(axis[1] * sqrt(max(0.0, temp)), e2, e3)
                            extra = true
                        end
                    end
                else
                    denom2 = 1.0 - axisqr[1] / axisqr[2]
                    denom3 = 1.0 - axisqr[1] / axisqr[3]
                    if denom2 > 0.0 && denom3 > 0.0
                        e2 = point[2] / denom2
                        e3 = point[3] / denom3
                        t1 = e2 / axis[2]
                        t2 = e3 / axis[3]
                        temp = 1.0 - t1 * t1 - t2 * t2
                        if temp > 0.0
                            epoint = SVector(axis[1] * sqrt(temp), e2, e3)
                            extra = true
                        end
                    end
                end
            end
            pt = SVector(point)
            if extra && _m_norm(epoint - pt) < _m_norm(spoint - pt)
                spoint = epoint
            end
        end
        if solutn == 1
            point = Tuple(spoint)
            bestpt = spoint
            bestht = _m_norm(bestpt - orignl)
        else
            height = sign * _m_norm(spoint - orignl)
            normal = _m_surfnm(axis[1], axis[2], axis[3], spoint)
            err = SVector(orignl[1] - (spoint[1] + height * normal[1]),
                orignl[2] - (spoint[2] + height * normal[2]),
                orignl[3] - (spoint[3] + height * normal[3]))
            newerr = 0.0
            better = true
            if solutn > 2
                if _m_iszero(err)
                    solvng = false
                else
                    newerr = _m_norm(err)
                end
                better = newerr < olderr
            else
                newerr = _m_norm(err)
            end
            if better
                olderr = newerr
                bestht = height
                bestpt = spoint
                if solutn <= 6
                    errp = _m_perp(err, normal, false)
                    _m_iszero(errp) && (solvng = false)
                    if sign >= 0.0
                        pnorm == 0.0 && _m_nearpt_bug("the norm of the scaled point is zero")
                        solvng && (errp = (_m_norm(spoint) / pnorm) * errp)
                    end
                    np = spoint + errp
                    point = Tuple(np)
                    np == spoint && (solvng = false)
                end
            else
                solvng = false
            end
        end
        solutn += 1
        solvng = solvng && solutn <= 6
    end
    bp = (1.0 / scale) * bestpt
    npoint = SVector(
        iorder[1] == 1 ? bp[1] : iorder[2] == 1 ? bp[2] : bp[3],
        iorder[1] == 2 ? bp[1] : iorder[2] == 2 ? bp[2] : bp[3],
        iorder[1] == 3 ? bp[1] : iorder[2] == 3 ? bp[2] : bp[3])
    return npoint, bestht / scale
end

"""
    _m_nearpt(positn, a, b, c) for another number type

The plain result of the `Float64` method with derivative parts from the stationarity
condition: the near point is `p_i a_i^2 / (a_i^2 + lambda)` where `lambda` is the root
of `sum((p_i a_i / (a_i^2 + lambda))^2) = 1`, refined by Newton steps in the number
type. Where that form does not hold (the center, and the interior points whose near
point is not unique or not smooth) the derivative parts are zero.
"""
function _m_nearpt(positn::SVector{3}, a::Real, b::Real, c::Real)
    p0 = SVector(primal(positn[1]), primal(positn[2]), primal(positn[3]))
    a0, b0, c0 = primal(a), primal(b), primal(c)
    np0, alt0 = _m_nearpt(p0, a0, b0, c0)
    T = typeof(positn[1] * a * b * c / 1.0)
    flat = SVector{3,T}(np0[1], np0[2], np0[3]), convert(T, alt0)
    ax0 = (a0, b0, c0)
    k = abs(p0[1]) >= abs(p0[2]) ? (abs(p0[1]) >= abs(p0[3]) ? 1 : 3) : (abs(p0[2]) >= abs(p0[3]) ? 2 : 3)
    (p0[k] == 0.0 || np0[k] == 0.0) && return flat
    lam0 = ax0[k]^2 * (p0[k] / np0[k] - 1.0)
    size = max(a0, b0, c0)
    for i in 1:3
        den = ax0[i]^2 + lam0
        (den > 0.0 && abs(p0[i] * ax0[i]^2 / den - np0[i]) <= 1.0e-8 * size) || return flat
    end
    axs = (a * one(T), b * one(T), c * one(T))
    lam = convert(T, lam0)
    for _ in 1:_ad_order(T)       # invariant: one more derivative order of lam is exact
        q = -one(T)
        dq = zero(T)
        for i in 1:3
            den = axs[i] * axs[i] + lam
            r = positn[i] * axs[i] / den
            q += r * r
            dq -= 2.0 * r * r / den
        end
        lam -= _with_primal(q, 0.0) / dq
    end
    # `lam` is assigned in the loop above; a closure that captured it would hold it in a box.
    s = let lam = lam
        SVector(ntuple(i -> _with_primal(positn[i] * axs[i] * axs[i] / (axs[i] * axs[i] + lam), np0[i]), Val(3)))
    end
    n = _m_surfnm(axs[1], axs[2], axs[3], s)
    return s, _with_primal(_m_dot(positn - s, n), alt0)
end

"""
    nearpt(positn, a, b, c) -> (npoint::SVector{3}, alt)

Return the point of an ellipsoid that is nearest to the point `positn`, and the
altitude of `positn` above the ellipsoid, as `nearpt_c` of CSPICE.

The ellipsoid has the semi-axes `a`, `b`, `c`. `alt` is negative for a point inside.
The results have the unit of the arguments. The iteration is the one of CSPICE.

The arguments can be `ForwardDiff` or `Taylor` numbers. The value is then the one of
the `Float64` iteration. The derivatives are zero at the center and at an interior
point whose nearest point is not unique.

Throws `SpiceError(:badaxislength)` if a semi-axis is not positive or if the ratio of
the longest to the shortest one is too large. Throws `SpiceError(:inputstoolarge)` if
the product of a semi-axis and a coordinate would overflow.

# Examples
```jldoctest
julia> nearpt([3.5, 0, 0], 1, 2, 3)
([1.0, 0.0, 0.0], 2.5)
```

See also [`dnearp`](@ref), [`recgeo`](@ref), [`npedln`](@ref), [`surfpt`](@ref).
"""
nearpt(positn::AbstractVector, a::Real, b::Real, c::Real) =
    _m_nearpt(_m_v3(positn), _m_f(a), _m_f(b), _m_f(c))

"""
    dnearp(state, a, b, c) -> Union{Tuple{SVector{6},SVector{2}},Nothing}

Return the state of the point of an ellipsoid that is nearest to the position of
`state`, and the altitude with its rate, as `dnearp_c` of CSPICE.

The first result is the position and the velocity of the near point. The second is
`(alt, dalt/dt)`. The result is `nothing` where the velocity of the near point is not
defined, for example at the center of a sphere; CSPICE returns a `found` flag.

Throws as [`nearpt`](@ref), and `SpiceError(:badarg)` if `state` does not have 6
elements.

# Examples
```jldoctest
julia> dnearp([2, 0, 0, 0, 1, 0], 1, 1, 1)
([1.0, 0.0, 0.0, 0.0, 0.5, 0.0], [1.0, 0.0])
```

See also [`nearpt`](@ref), [`surfpv`](@ref).
"""
function dnearp(state::AbstractVector, a::Real, b::Real, c::Real)
    s = _m_v6(state)
    a, b, c = _m_f(a), _m_f(b), _m_f(c)
    pos, vel = _m_pv(s)
    nearp, alt = _m_nearpt(pos, a, b, c)
    g = SVector(1.0 / (a * a), 1.0 / (b * b), 1.0 / (c * c))
    zenith = pos - nearp
    grad = SVector(g[1] * nearp[1], g[2] * nearp[2], g[3] * nearp[3])
    norml, length = _m_unorm(grad)
    l = _m_dot(zenith, norml) / length
    d = SVector(l * g[1] + 1.0, l * g[2] + 1.0, l * g[3] + 1.0)
    (primal(d[1]) == 0.0 || primal(d[2]) == 0.0 || primal(d[3]) == 0.0) && return nothing
    m = SVector(1.0 / d[1], 1.0 / d[2], 1.0 / d[3])
    denom = grad[1] * m[1] * grad[1] + grad[2] * m[2] * grad[2] + grad[3] * m[3] * grad[3]
    primal(denom) == 0.0 && return nothing
    lprime = (grad[1] * m[1] * vel[1] + grad[2] * m[2] * vel[2] + grad[3] * m[3] * vel[3]) / denom
    temp = 1.0 * vel + (-lprime) * grad
    dn = SVector(m[1] * temp[1], m[2] * temp[2], m[3] * temp[3])
    return _m_state(nearp - zero(dn), dn), SVector(alt - zero(l), _m_dot(norml, vel))
end

"""
    surfpv(stvrtx, stdir, a, b, c) -> Union{SVector{6},Nothing}

Return the state of the point where a ray meets an ellipsoid, as `surfpv_c` of CSPICE.

`stvrtx` is the state of the start point of the ray, and `stdir` holds the direction
and its time derivative. The result is the position and the velocity of the intercept
point of [`surfpt`](@ref). The result is `nothing` if the ray misses, if it is tangent
to the surface, or if the velocity would overflow; CSPICE returns a `found` flag.

Throws `SpiceError(:invalidvertex)` if the start point is on the surface and
`SpiceError(:badarg)` if a state does not have 6 elements; other errors as
[`surfpt`](@ref).

# Examples
```jldoctest
julia> surfpv([4, 0, 0, 0, 0, 0], [-1, 0, 0, 0, 0.5, 0], 1, 1, 1) == [1, 0, 0, 0, 1.5, 0]
true
```

See also [`surfpt`](@ref), [`dnearp`](@ref).
"""
function surfpv(stvrtx::AbstractVector, stdir::AbstractVector, a::Real, b::Real, c::Real)
    sv = _m_v6(stvrtx)
    sd = _m_v6(stdir)
    a, b, c = _m_f(a), _m_f(b), _m_f(c)
    v, dv = _m_pv(sv)
    found, x = _m_surfpt(v, _m_pv(sd)[1], a, b, c)
    found || return nothing
    u, du = _m_pv(_m_dvhat(sd, false))
    n = _m_surfnm(a, b, c, x)
    vmx = v - x
    t1 = v[1] / a
    t2 = v[2] / b
    t3 = v[3] / c
    level = t1 * t1 + t2 * t2 + t3 * t3
    primal(level) == 1.0 && throw(SpiceError(:invalidvertex,
        "the ray vertex is on the ellipsoid: the surface velocity is not defined"))
    udn = _m_dot(u, n)
    primal(udn) == 0.0 && return nothing
    second = _m_dot(vmx, n) * du
    dsnum = udn * _m_dot(dv, n) - _m_dot(vmx, n) * _m_dot(du, n)
    third = dsnum * u
    m = max(primal(_m_norm(second)), primal(_m_norm(third)), 1.0)
    m > floatmax(Float64) / 10.0 * (primal(udn) * primal(udn)) && return nothing
    r = 1.0 / udn
    dx = 1.0 * dv + (-r) * second + (-(r * r)) * third
    return _m_state(x - zero(dx), dx)
end

# --- planes -------------------------------------------------------------------------------

@inline _m_plane(normal::SVector{3}, constant) =
    primal(constant) < 0.0 ? SpicePlane(-normal, -constant) : SpicePlane(normal, constant)

"""
    nvc2pl(normal, konst) -> SpicePlane

Return the plane of the points `x` with `normal ⋅ x == konst`, as `nvc2pl_c` of CSPICE.

The result has a unit normal and a constant that is not negative: `normal` is scaled,
and both values change sign if needed.

Throws `SpiceError(:zerovector)` if `normal` is zero.

# Examples
```jldoctest
julia> nvc2pl([0, 0, 2], 4)
SpicePlane{Float64}([0.0, 0.0, 1.0], 2.0)
```

See also [`nvp2pl`](@ref), [`psv2pl`](@ref), [`pl2nvc`](@ref), [`SpicePlane`](@ref).
"""
function nvc2pl(normal::AbstractVector, konst::Real)
    u, mag = _m_unorm(_m_v3(normal))
    primal(mag) == 0.0 && throw(SpiceError(:zerovector, "the normal of a plane must not be zero"))
    return _m_plane(u, _m_f(konst) / mag)
end

"""
    nvp2pl(normal, point) -> SpicePlane

Return the plane through `point` with the normal vector `normal`, as `nvp2pl_c` of
CSPICE. The result has a unit normal and a constant that is not negative.

Throws `SpiceError(:zerovector)` if `normal` is zero.

See also [`nvc2pl`](@ref), [`psv2pl`](@ref), [`pl2nvp`](@ref).
"""
function nvp2pl(normal::AbstractVector, point::AbstractVector)
    n = _m_v3(normal)
    _m_iszero(n) && throw(SpiceError(:zerovector, "the normal of a plane must not be zero"))
    u = _m_hat(n)
    # The zero (which gives the number type of `point` to the normal) is subtracted, not
    # added, because `-0.0 + 0.0` is `0.0`: a zero component keeps its sign, as in CSPICE.
    return _m_plane(u - zero(_m_v3(point)), _m_dot(_m_v3(point), u))
end

"""
    psv2pl(point, span1, span2) -> SpicePlane

Return the plane through `point` that the vectors `span1` and `span2` span, as
`psv2pl_c` of CSPICE. The result has a unit normal and a constant that is not negative.

Throws `SpiceError(:degeneratecase)` if `span1` and `span2` are parallel or one is
zero.

See also [`nvc2pl`](@ref), [`nvp2pl`](@ref), [`pl2psv`](@ref).
"""
function psv2pl(point::AbstractVector, span1::AbstractVector, span2::AbstractVector)
    n = _m_ucrss(_m_v3(span1), _m_v3(span2))
    _m_iszero(n) && throw(SpiceError(:degeneratecase, "the spanning vectors are parallel"))
    p = _m_v3(point)
    return _m_plane(n - zero(p), _m_dot(n, p))      # subtracted: see `nvp2pl`
end

"""
    pl2nvc(plane) -> (normal::SVector{3}, konst)

Return the normal vector and the constant of a plane, as `pl2nvc_c` of CSPICE.

See also [`nvc2pl`](@ref), [`pl2nvp`](@ref), [`pl2psv`](@ref).
"""
pl2nvc(plane::SpicePlane) = plane.normal, plane.constant

"""
    pl2nvp(plane) -> (normal::SVector{3}, point::SVector{3})

Return the normal vector of a plane and the point of the plane nearest to the origin,
as `pl2nvp_c` of CSPICE.

See also [`nvp2pl`](@ref), [`pl2nvc`](@ref), [`pl2psv`](@ref).
"""
pl2nvp(plane::SpicePlane) = plane.normal, plane.constant * plane.normal

"""
Right-handed orthonormal frame `(x, y, z)` whose first axis is along `x`, as `frame_c`:
`y` has a zero in the smallest component of `x`. A zero `x` gives the standard basis.
"""
@inline function _m_frame(x0::SVector{3})
    x = _m_hat(x0)
    a = x[1] * x[1]
    b = x[2] * x[2]
    c = x[3] * x[3]
    z0 = zero(a)
    o = one(a)
    primal(a + b + c) == 0.0 &&
        return SVector(o, z0, z0), SVector(z0, o, z0), SVector(z0, z0, o)
    a0, b0, c0 = primal(a), primal(b), primal(c)
    if a0 <= b0 && a0 <= c0
        f = sqrt(b + c)
        y2 = -x[3] / f
        y3 = x[2] / f
        return x, SVector(z0, y2, y3), SVector(f, -x[1] * y3, x[1] * y2)
    elseif b0 <= a0 && b0 <= c0
        f = sqrt(a + c)
        y3 = -x[1] / f
        y1 = x[3] / f
        return x, SVector(y1, z0, y3), SVector(x[2] * y3, f, -x[2] * y1)
    end
    f = sqrt(a + b)
    y1 = -x[2] / f
    y2 = x[1] / f
    return x, SVector(y1, y2, z0), SVector(-x[3] * y2, x[3] * y1, f)
end

"""
    pl2psv(plane) -> (point::SVector{3}, span1::SVector{3}, span2::SVector{3})

Return the point of a plane nearest to the origin and two orthonormal vectors that
span the plane, as `pl2psv_c` of CSPICE. `span1 × span2` is the normal of the plane.

See also [`psv2pl`](@ref), [`pl2nvc`](@ref), [`pl2nvp`](@ref).
"""
function pl2psv(plane::SpicePlane)
    _, span1, span2 = _m_frame(plane.normal)
    return plane.constant * plane.normal, span1, span2
end

@inline _m_vprjp(v::SVector{3}, plane::SpicePlane) =
    1.0 * v + (plane.constant - _m_dot(v, plane.normal)) * plane.normal

"""
    vprjp(vin, plane) -> SVector{3}

Return the orthogonal projection of the vector `vin` on a plane, as `vprjp_c` of
CSPICE.

Throws `SpiceError(:badarg)` if `vin` does not have 3 elements.

# Examples
```jldoctest
julia> vprjp([1, 2, 3], nvc2pl([0, 0, 1], 0)) == [1, 2, 0]
true
```

See also [`vprjpi`](@ref), [`vproj`](@ref), [`pjelpl`](@ref).
"""
vprjp(vin::AbstractVector, plane::SpicePlane) = _m_vprjp(_m_v3(vin), plane)

"""
    vprjpi(vin, projpl, invpl) -> Union{SVector{3},Nothing}

Return the vector in the plane `invpl` whose orthogonal projection on the plane
`projpl` is `vin`, as `vprjpi_c` of CSPICE.

The result is `nothing` if the planes are orthogonal, or so near to orthogonal that the
result would overflow. CSPICE returns a `found` flag.

Throws `SpiceError(:badarg)` if `vin` does not have 3 elements.

# Examples
```jldoctest
julia> vprjpi([1, 2, 0], nvc2pl([0, 0, 1], 0), nvc2pl([0, 1, 1], 0)) == [1, 2, -2]
true
```

See also [`vprjp`](@ref).
"""
function vprjpi(vin::AbstractVector, projpl::SpicePlane, invpl::SpicePlane)
    v = _m_v3(vin)
    numer = invpl.constant - _m_dot(v, invpl.normal)
    denom = _m_dot(projpl.normal, invpl.normal)
    bound = 10.0 / floatmax(Float64)
    limit = abs(primal(numer)) < 1.0 ? abs(bound) : abs(bound * primal(numer))
    abs(primal(denom)) > limit || return nothing
    return 1.0 * v + (numer / denom) * projpl.normal
end

"""
    inrypl(vertex, dir, plane) -> (nxpts::Int, xpt::SVector{3})

Return the number of points where a ray meets a plane, and the point, as `inrypl_c` of
CSPICE. The ray starts at `vertex` and has the direction `dir`.

| `nxpts` | Case                      | `xpt`             |
|:--------|:--------------------------|:------------------|
| 1       | one point                 | that point        |
| 0       | the ray misses the plane  | the zero vector   |
| -1      | the ray is in the plane   | `vertex`          |

Throws `SpiceError(:zerovector)` if `dir` is zero. Throws `SpiceError(:vectortoobig)`
if the norm of `vertex` or the plane constant is a third of `floatmax(Float64)` or
more.

# Examples
```jldoctest
julia> inrypl([1, 0, 4], [1, 0, -1], nvc2pl([0, 0, 1], 0))
(1, [5.0, 0.0, 0.0])
```

See also [`inelpl`](@ref), [`inedpl`](@ref), [`surfpt`](@ref).
"""
function inrypl(vertex::AbstractVector, dir::AbstractVector, plane::SpicePlane)
    vtx = _m_v3(vertex)
    d = _m_v3(dir)
    toobig = floatmax(Float64) / 3.0
    vnrm = _m_norm(vtx)
    primal(vnrm) >= toobig && throw(SpiceError(:vectortoobig,
        "the vertex of the ray is too far from the origin"))
    normal, constant = plane.normal, plane.constant
    primal(constant) >= toobig && throw(SpiceError(:vectortoobig,
        "the plane is too far from the origin"))
    udir = _m_hat(d)
    _m_iszero(udir) && throw(SpiceError(:zerovector, "the direction of the ray is the zero vector"))
    mscale = _m_max(_m_abs(constant) - zero(vnrm), _m_abs(vnrm) - zero(constant))
    if primal(mscale) != 0.0
        sclvtx = (1.0 / mscale) * vtx
        sclcon = constant / mscale
    else
        sclvtx = vtx - zero(SVector{3,typeof(mscale)})
        sclcon = constant - zero(mscale)
    end
    primal(mscale) > 1.0 && (toobig = toobig / primal(mscale))
    prjdif = sclcon - _m_dot(sclvtx, normal)
    prjdir = _m_dot(udir, normal)
    none = zero(1.0 * sclvtx + prjdif * udir)
    if primal(prjdif) == 0.0
        return (primal(_m_dot(normal, udir)) == 0.0 ? -1 : 1), vtx + none
    end
    pd, pf = primal(prjdir), primal(prjdif)
    samesign = (pd > 0.0 && pf > 0.0) || (pd < 0.0 && pf < 0.0)
    (!samesign || abs(pf) >= abs(pd) * toobig) && return 0, none
    scale = _m_abs(prjdif) / _m_abs(prjdir)
    return 1, mscale * (1.0 * sclvtx + scale * udir)
end

# --- ellipses -----------------------------------------------------------------------------

"""
Roots of `a x^2 + b x + c` as `rquad_c`: `((re1, im1), (re2, im2))`. Requires that `a`
and `b` are not both zero (`SpiceError(:degeneratecase)`).
"""
@inline function _m_rquad(a::T, b::T, c::T) where {T<:Real}
    (primal(a) == 0.0 && primal(b) == 0.0) && throw(SpiceError(:degeneratecase,
        "both the first and the second degree coefficients are zero"))
    scale = max(abs(primal(a)), abs(primal(b)), abs(primal(c)))
    zeroed = (primal(a) != 0.0 && primal(a) / scale == 0.0) ||
             (primal(b) != 0.0 && primal(b) / scale == 0.0) ||
             (primal(c) != 0.0 && primal(c) / scale == 0.0)
    sqr, lin, con = zeroed ? (a, b, c) : (a / scale, b / scale, c / scale)
    z = zero(sqr)
    if primal(sqr) != 0.0
        discrm = lin * lin - 4.0 * sqr * con
        if primal(discrm) >= 0.0
            if primal(lin) < 0.0
                r1 = (-lin + sqrt(discrm)) / (2.0 * sqr)
                return (r1, z), ((con / sqr) / r1, z)
            elseif primal(lin) > 0.0
                r2 = (-lin - sqrt(discrm)) / (2.0 * sqr)
                return ((con / sqr) / r2, z), (r2, z)
            end
            r1 = sqrt(discrm) / (2.0 * sqr)
            return (r1, z), (-r1, z)
        end
        re = -lin / (2.0 * sqr)
        im = sqrt(-discrm) / (2.0 * sqr)
        return (re, im), (re, -im)
    end
    r = -con / lin
    return (r, z), (r, z)
end

"Unit vector of two components with the scaling of `vhatg_c`."
@inline function _m_hat2(x, y)
    scale = max(abs(primal(x)), abs(primal(y)))
    scale == 0.0 && return zero(x) - zero(y), zero(x) - zero(y)
    t1 = x / scale
    t2 = y / scale
    vmag = scale * sqrt(t1 * t1 + t2 * t2)
    return x / vmag, y / vmag
end

"""
Diagonalization of the symmetric matrix `[s11 s12; s12 s22]` as `diags2_c`: the two
diagonal values and the rotation `(c11, c12, c21, c22)` whose columns are the
eigenvectors.
"""
@inline function _m_diags2(s11::T, s12::T, s22::T) where {T<:Real}
    o = one(T)
    z = zero(T)
    primal(s12) == 0.0 && return s11, s22, (o, z, z, o)
    scale = max(abs(primal(s11)), abs(primal(s12)), abs(primal(s22)))
    a = s11 / scale
    b = s12 / scale
    c = s22 / scale
    root1, root2 = _m_rquad(o, -(a + c), a * c - b * b)
    d1 = root1[1]
    d2 = root2[1]
    if abs(primal(d1 - a)) >= abs(primal(d1 - c))
        e1, e2 = _m_hat2(b, _m_max(d1 - a, _m_abs(b)))
        return d2 * scale, d1 * scale, (e2, e1, -e1, e2)
    end
    e1, e2 = _m_hat2(_m_max(d1 - c, _m_abs(b)), b)
    return d1 * scale, d2 * scale, (e1, -e2, e2, e1)
end

@inline function _m_saelgv(vec1::SVector{3}, vec2::SVector{3})
    n1 = _m_norm(vec1)
    n2 = _m_norm(vec2)
    scale = max(abs(primal(n1)), abs(primal(n2)))
    scale == 0.0 && return zero(vec1) - zero(vec2), zero(vec1) - zero(vec2)
    t1 = vec1 / scale - zero(vec2)
    t2 = vec2 / scale - zero(vec1)
    s12 = _m_dot(t1, t2)
    d1, d2, c = _m_diags2(_m_dot(t1, t1), s12, _m_dot(t2, t2))
    # column `major` of the rotation gives the major axis
    if abs(primal(d1)) >= abs(primal(d2))
        return scale * (c[1] * t1 + c[3] * t2), scale * (c[2] * t1 + c[4] * t2)
    end
    return scale * (c[2] * t1 + c[4] * t2), scale * (c[1] * t1 + c[3] * t2)
end

"""
    saelgv(vec1, vec2) -> (smajor::SVector{3}, sminor::SVector{3})

Return the semi-major and the semi-minor axis of the ellipse that two generating
vectors define, as `saelgv_c` of CSPICE.

The ellipse is the set of the points `cos(θ) * vec1 + sin(θ) * vec2`. The results are
orthogonal, and `smajor` is not shorter than `sminor`.

Throws `SpiceError(:badarg)` if a vector does not have 3 elements.

# Examples
```jldoctest
julia> saelgv([0, 1, 0], [2, 0, 0])
([2.0, 0.0, 0.0], [0.0, 1.0, 0.0])
```

See also [`cgv2el`](@ref).
"""
saelgv(vec1::AbstractVector, vec2::AbstractVector) = _m_saelgv(_m_v3(vec1), _m_v3(vec2))

@inline function _m_cgv2el(center::SVector{3}, vec1::SVector{3}, vec2::SVector{3})
    smajor, sminor = _m_saelgv(vec1, vec2)
    return SpiceEllipse(center, smajor, sminor)
end

"""
    cgv2el(center, vec1, vec2) -> SpiceEllipse

Return the ellipse with the center `center` and the generating vectors `vec1` and
`vec2`, as `cgv2el_c` of CSPICE. The semi-axes of the result come from
[`saelgv`](@ref).

Throws `SpiceError(:badarg)` if a vector does not have 3 elements.

See also [`el2cgv`](@ref), [`SpiceEllipse`](@ref).
"""
cgv2el(center::AbstractVector, vec1::AbstractVector, vec2::AbstractVector) =
    _m_cgv2el(_m_v3(center), _m_v3(vec1), _m_v3(vec2))

"""
    el2cgv(ellipse) -> (center, smajor, sminor)

Return the center and the two semi-axes of an ellipse, as `el2cgv_c` of CSPICE.

See also [`cgv2el`](@ref).
"""
el2cgv(ellipse::SpiceEllipse) = ellipse.center, ellipse.semimajor, ellipse.semiminor

"`(found, ellipse)` of `inedpl_c` for positive semi-axes."
@inline function _m_inedpl(a, b, c, plane::SpicePlane)
    point, span1, span2 = pl2psv(plane)
    T = typeof(point[1] * a * b * c)
    none = SpiceEllipse(zero(SVector{3,T}), zero(SVector{3,T}), zero(SVector{3,T}))
    maxrad = max(abs(primal(a)), abs(primal(b)), abs(primal(c)))
    primal(_m_norm(point)) > maxrad && return false, none
    ds = SVector(1.0 / a, 1.0 / b, 1.0 / c)
    inv = SVector(a, b, c)
    # In the distorted space the ellipsoid is the unit sphere and the section a circle.
    dn = _m_ucrss(ds .* span1, ds .* span2)
    _m_iszero(dn) && throw(SpiceError(:degeneratecase, "the spanning vectors are parallel"))
    dplane = _m_plane(dn, _m_dot(dn, ds .* point))
    center, vec1, vec2 = pl2psv(dplane)
    dist = _m_norm(center)
    primal(dist) > 1.0 && return false, none
    t = 1.0 - dist * dist
    rcircl = primal(t) <= 0.0 ? zero(t) : primal(t) >= 1.0 ? one(t) : sqrt(t)
    return true, _m_cgv2el(inv .* center, inv .* (rcircl * vec1), inv .* (rcircl * vec2))
end

"""
    inedpl(a, b, c, plane) -> Union{SpiceEllipse,Nothing}

Return the ellipse where a plane cuts an ellipsoid, or `nothing` if they do not meet,
as `inedpl_c` of CSPICE.

The ellipsoid has the semi-axes `a`, `b`, `c`. CSPICE returns a `found` flag in place
of `nothing`. A tangent plane gives an ellipse of one point, with zero semi-axes.

Throws `SpiceError(:degeneratecase)` if a semi-axis is not positive, as in CSPICE.

# Examples
A plane at the height 0.6 cuts the unit sphere in a circle of radius 0.8:
```jldoctest
julia> e = inedpl(1, 1, 1, nvc2pl([0, 0, 1], 0.6));

julia> e.center == [0, 0, 0.6]
true

julia> vnorm(e.semimajor) ≈ 0.8 && vnorm(e.semiminor) ≈ 0.8
true
```

See also [`edlimb`](@ref), [`inelpl`](@ref), [`inrypl`](@ref).
"""
function inedpl(a::Real, b::Real, c::Real, plane::SpicePlane)
    a, b, c = _m_f(a), _m_f(b), _m_f(c)
    _m_axes_ok(a, b, c) || throw(SpiceError(:degeneratecase,
        "semi-axes: a = $(primal(a)), b = $(primal(b)), c = $(primal(c))"))
    found, ellipse = _m_inedpl(a, b, c, plane)
    return found ? ellipse : nothing
end

"""
    inelpl(ellips, plane) -> (nxpts::Int, xpt1::SVector{3}, xpt2::SVector{3})

Return the number of points where an ellipse meets a plane, and the points, as
`inelpl_c` of CSPICE.

`nxpts` is 0, 1 or 2, or -1 if the ellipse is in the plane. With one point, `xpt1` and
`xpt2` are equal. With no point, and for -1, both are zero.

Throws `SpiceError(:invalidplane)` if the plane has a zero normal or a negative
constant. Throws `SpiceError(:invalidellipse)` if the semi-axes of the ellipse are not
orthogonal.

# Examples
```jldoctest
julia> e = cgv2el([0, 0, 0], [2, 0, 0], [0, 1, 0]);

julia> inelpl(e, nvc2pl([1, 0, 0], 2))
(1, [2.0, 0.0, 0.0], [2.0, 0.0, 0.0])

julia> inelpl(e, nvc2pl([1, 0, 0], 3))
(0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
```

See also [`inedpl`](@ref), [`inrypl`](@ref).
"""
function inelpl(ellips::SpiceEllipse, plane::SpicePlane)
    normal, inpcon = plane.normal, plane.constant
    _m_iszero(normal) && throw(SpiceError(:invalidplane, "the plane has a zero normal vector"))
    primal(inpcon) < 0.0 && throw(SpiceError(:invalidplane,
        "the plane has the negative constant $(primal(inpcon))"))
    center, smajor, sminor = ellips.center, ellips.semimajor, ellips.semiminor
    if !_m_iszero(sminor)
        u1, m1 = _m_unorm(smajor)
        u2, _ = _m_unorm(sminor)
        sep = primal(m1) == 0.0 ? 0.0 : primal(_m_usep(u1, u2))
        abs(sep - Float64(pi) / 2) > 1.0e-9 && throw(SpiceError(:invalidellipse,
            "the ellipse has semi-axes that are not orthogonal; their separation is $sep rad"))
    end
    none = zero(center) - zero(normal)
    if _m_iszero(smajor)
        primal(_m_dot(center, normal)) == primal(inpcon) && return 1, center + none, center + none
        return 0, none, none
    end
    # the plane relative to the center of the ellipse
    point = inpcon * normal - center
    u = _m_hat(normal)
    # The zero is subtracted, not added: `-0.0 + 0.0` is `0.0`, and the sign of a zero
    # component of the normal decides between the angles pi and -pi below.
    trans = _m_plane(u - none, _m_dot(point, u))
    n, constant = trans.normal, trans.constant
    v1 = _m_dot(smajor, n)
    v2 = _m_dot(sminor, n)
    if primal(v1) == 0.0 && primal(v2) == 0.0
        return (primal(constant) == 0.0 ? -1 : 0), none, none
    end
    scale = max(abs(primal(v1)), abs(primal(v2)))
    t1 = v1 / scale
    t2 = v2 / scale
    vn = scale * sqrt(t1 * t1 + t2 * t2)
    primal(vn) < primal(constant) && return 0, none, none
    alpha = _m_acos(constant / vn)
    beta = atan(v2, v1)
    angle1 = beta - alpha
    angle2 = beta + alpha
    nxpts = (_m_iszero(sminor) || primal(angle1) == primal(angle2)) ? 1 : 2
    s1, c1 = sincos(angle1)
    s2, c2 = sincos(angle2)
    return nxpts, 1.0 * center + c1 * smajor + s1 * sminor + none,
    1.0 * center + c2 * smajor + s2 * sminor + none
end

@inline function _m_pjelpl(elin::SpiceEllipse, plane::SpicePlane)
    # `vperp_c` gives zero for a zero second argument
    perp(v) = _m_iszero(plane.normal) ? zero(v) - zero(plane.normal) : _m_perp(v, plane.normal, true)
    return _m_cgv2el(_m_vprjp(elin.center, plane), perp(elin.semimajor), perp(elin.semiminor))
end

"""
    pjelpl(elin, plane) -> SpiceEllipse

Return the orthogonal projection of the ellipse `elin` on a plane, as `pjelpl_c` of
CSPICE.

See also [`vprjp`](@ref), [`cgv2el`](@ref).
"""
pjelpl(elin::SpiceEllipse, plane::SpicePlane) = _m_pjelpl(elin, plane)

"""
    edlimb(a, b, c, viewpt) -> SpiceEllipse

Return the limb of an ellipsoid as seen from the point `viewpt`, as `edlimb_c` of
CSPICE.

The limb is the set of the surface points where the line of sight is tangent to the
ellipsoid. The ellipsoid has the semi-axes `a`, `b`, `c`. A view point on the surface
gives an ellipse of one point.

Throws `SpiceError(:degeneratecase)` if a semi-axis is not positive or if `viewpt` is
inside the ellipsoid.

# Examples
The limb of the unit sphere from a distance of 2:
```jldoctest
julia> limb = edlimb(1, 1, 1, [0, 0, 2]);

julia> limb.center == [0, 0, 0.5]
true

julia> vnorm(limb.semimajor) ≈ sqrt(3) / 2
true
```

See also [`inedpl`](@ref), [`edterm`](@ref).
"""
function edlimb(a::Real, b::Real, c::Real, viewpt::AbstractVector)
    a, b, c = _m_f(a), _m_f(b), _m_f(c)
    vp = _m_v3(viewpt)
    _m_axes_ok(a, b, c) || throw(SpiceError(:degeneratecase,
        "semi-axis lengths: a = $(primal(a)), b = $(primal(b)), c = $(primal(c))"))
    scale = _m_max(_m_max(a, b), c)
    scla = a / scale
    sclb = b / scale
    sclc = c / scale
    scla2 = scla * scla
    sclb2 = sclb * sclb
    sclc2 = sclc * sclc
    (primal(scla2) == 0.0 || primal(sclb2) == 0.0 || primal(sclc2) == 0.0) &&
        throw(SpiceError(:degeneratecase, "a semi-axis is too small"))
    v = (1.0 / scale) * vp
    level = (v[1] * v[1] / scla2) + (v[2] * v[2] / sclb2) + (v[3] * v[3] / sclc2)
    primal(level) < 1.0 && throw(SpiceError(:degeneratecase,
        "the viewing point is inside the ellipsoid"))
    # The limb is in the plane of the points x with (x ⋅ grad at the view point) = 1.
    u, mag = _m_unorm(SVector(v[1] / scla2, v[2] / sclb2, v[3] / sclc2))
    primal(mag) == 0.0 && throw(SpiceError(:zerovector, "the normal of a plane must not be zero"))
    found, limb = _m_inedpl(scla, sclb, sclc, _m_plane(u, 1.0 / mag))
    found || throw(SpiceError(:degeneratecase,
        "the ellipsoid shape and the viewing geometry are too extreme; the limb was not found"))
    return SpiceEllipse(scale * limb.center, scale * limb.semimajor, scale * limb.semiminor)
end

# --- nearest points on lines and ellipses ------------------------------------------------

"""
    nplnpt(linpt, lindir, point) -> (pnear::SVector{3}, dist)

Return the point of a line that is nearest to `point`, and the distance between the
two, as `nplnpt_c` of CSPICE. The line goes through `linpt` and has the direction
`lindir`.

Throws `SpiceError(:zerovector)` if `lindir` is zero.

# Examples
```jldoctest
julia> nplnpt([0, 0, 0], [1, 0, 0], [3, 4, 0])
([3.0, 0.0, 0.0], 4.0)
```

See also [`npelpt`](@ref), [`npedln`](@ref).
"""
function nplnpt(linpt::AbstractVector, lindir::AbstractVector, point::AbstractVector)
    lp = _m_v3(linpt)
    ld = _m_v3(lindir)
    p = _m_v3(point)
    _m_iszero(ld) && throw(SpiceError(:zerovector, "the direction vector must not be zero"))
    pnear = _m_proj(p - lp, ld, true) + lp
    return pnear, _m_norm(pnear - p)
end

@inline function _m_npelpt(point::SVector{3}, ellips::SpiceEllipse)
    center = ellips.center
    minlen = _m_norm(ellips.semiminor)
    majlen = _m_norm(ellips.semimajor)
    min(primal(majlen), primal(minlen)) == 0.0 && throw(SpiceError(:degeneratecase,
        "ellipse semi-axis lengths: $(primal(majlen)) $(primal(minlen))"))
    scale = 1.0 / majlen
    smajor = scale * ellips.semimajor
    sminor = scale * ellips.semiminor
    tmppnt = scale * (point - center)
    # rows of the rotation: the directions of the semi-axes and the ellipse normal
    r1 = _m_hat(smajor)
    r3 = _m_ucrss(smajor, sminor)
    r2 = _m_ucrss(r3, smajor)
    _m_iszero(r2) && throw(SpiceError(:dependentvectors,
        "the semi-axes of the ellipse are linearly dependent"))
    rot = _m_rows(r1, r2, r3)
    t = _m_mxv(rot, tmppnt)
    prjpnt = SVector(t[1], t[2], zero(t[1]))
    pn, _ = _m_nearpt(prjpnt, one(majlen), minlen / majlen, 2.0 * one(majlen))
    pnear = _m_mtxv(rot, majlen * pn) + center
    return pnear, _m_norm(pnear - point)
end

"""
    npelpt(point, ellips) -> (pnear::SVector{3}, dist)

Return the point of an ellipse that is nearest to `point`, and the distance between the
two, as `npelpt_c` of CSPICE.

Throws `SpiceError(:degeneratecase)` if a semi-axis of the ellipse is zero.

See also [`nplnpt`](@ref), [`nearpt`](@ref).
"""
npelpt(point::AbstractVector, ellips::SpiceEllipse) = _m_npelpt(_m_v3(point), ellips)

"""
    npedln(a, b, c, linept, linedr) -> (pnear::SVector{3}, dist)

Return the point of an ellipsoid that is nearest to a line, and the distance between
the two, as `npedln_c` of CSPICE.

The ellipsoid has the semi-axes `a`, `b`, `c`. The line goes through `linept` and has
the direction `linedr`. A line that meets the ellipsoid gives one of the intersection
points and the distance 0.

Throws `SpiceError(:zerovector)` if `linedr` is zero,
`SpiceError(:invalidaxislength)` if a semi-axis is not positive, and
`SpiceError(:degeneratecase)` if the shape is too extreme for the arithmetic.

# Examples
```jldoctest
julia> npedln(1, 1, 1, [2, 0, 0], [0, 0, 1])
([1.0, 0.0, 0.0], 1.0)
```

See also [`nearpt`](@ref), [`surfpt`](@ref), [`nplnpt`](@ref).
"""
function npedln(a::Real, b::Real, c::Real, linept::AbstractVector, linedr::AbstractVector)
    a, b, c = _m_f(a), _m_f(b), _m_f(c)
    lp = _m_v3(linept)
    udir, mag = _m_unorm(_m_v3(linedr))
    primal(mag) == 0.0 && throw(SpiceError(:zerovector, "the line direction is the zero vector"))
    _m_axes_ok(a, b, c) || throw(SpiceError(:invalidaxislength,
        "semi-axis lengths: a = $(primal(a)), b = $(primal(b)), c = $(primal(c))"))
    scale = _m_max(_m_max(a, b), c)
    scla = a / scale
    sclb = b / scale
    sclc = c / scale
    scla2 = scla * scla
    sclb2 = sclb * sclb
    sclc2 = sclc * sclc
    (primal(scla2) == 0.0 || primal(sclb2) == 0.0 || primal(sclc2) == 0.0) &&
        throw(SpiceError(:degeneratecase, "a semi-axis is too small"))
    sclpt = SVector(lp[1] / scale, lp[2] / scale, lp[3] / scale)
    for dir in (udir, -udir)
        found, pt = _m_surfpt(sclpt, dir, scla, sclb, sclc)
        found && return scale * pt, zero(scale * pt[1])
    end
    # The candidate points are where the surface normal is orthogonal to the line.
    normal = SVector(udir[1] / scla2, udir[2] / sclb2, udir[3] / sclc2)
    un, nmag = _m_unorm(normal)
    primal(nmag) == 0.0 && throw(SpiceError(:zerovector, "the normal of a plane must not be zero"))
    candpl = _m_plane(un, zero(nmag) / nmag)
    xfound, cand = _m_inedpl(scla, sclb, sclc, candpl)
    xfound || throw(SpiceError(:degeneratecase, "the candidate ellipse was not found"))
    ud2, dmag = _m_unorm(udir)
    prjpl = _m_plane(ud2, zero(dmag) / dmag)
    prjel = _m_pjelpl(cand, prjpl)
    prjpt = _m_vprjp(sclpt, prjpl)
    prjnpt, dist = _m_npelpt(prjpt, prjel)
    pnear = vprjpi(prjnpt, prjpl, candpl)
    pnear === nothing && throw(SpiceError(:degeneratecase, "the inverse projection was not found"))
    return scale * pnear, dist * scale
end

# --- terminator ---------------------------------------------------------------------------

"""
    _m_edterm(umbral, a, b, c, srcrad, srcpos, npts) -> Vector{SVector{3,Float64}}

Terminator points of `zzedterm_`: for each of `npts` directions about the axis from
the target to the source, the point of the ellipsoid where a plane that is tangent to
the spherical source is tangent to the ellipsoid, found by at most eleven corrections
of the tangent angle.
"""
function _m_edterm(umbral::Bool, a::Float64, b::Float64, c::Float64, srcrad::Float64,
        srcpos::SVector{3,Float64}, npts::Int)
    npts >= 1 || throw(SpiceError(:invalidsize,
        "the number of terminator points must be positive; it is $npts"))
    (a > 0.0 && b > 0.0 && c > 0.0) || throw(SpiceError(:invalidaxislength,
        "semi-axis lengths: a = $a, b = $b, c = $c"))
    srcrad > 0.0 || throw(SpiceError(:invalidradius, "the source radius must be positive; it is $srcrad"))
    d = _m_norm(srcpos)
    rmax = max(a, b, c)
    rmin = min(a, b, c)
    srcrad + rmax >= d && throw(SpiceError(:objectstooclose,
        "the source is too close to the target: distance $d, radii sum $(srcrad + rmax)"))
    z, x, y = _m_frame(-srcpos)
    trans = _m_cols(x, y, z)
    outang = umbral ? asin((srcrad - rmax) / d) : asin((srcrad + rmax) / d)
    inang = umbral ? asin((srcrad - rmin) / d) : asin((srcrad + rmin) / d)
    delta = 2.0 * Float64(pi) / npts
    halfpi = Float64(pi) / 2
    out = Vector{SVector{3,Float64}}(undef, npts)
    for i in 1:npts
        theta = (i - 1) * delta
        st, ct = sincos(theta)
        angle = umbral ? outang : inang
        prvdif = 2.0 * Float64(pi)
        prvang = angle + halfpi
        nitr = 0
        pt = zero(SVector{3,Float64})
        # invariant: |angle - prvang| decreases; at most 11 passes
        while nitr <= 10 && abs(angle - prvang) < prvdif
            nitr += 1
            prvdif = abs(angle - prvang)
            prvang = angle
            sa, ca = sincos(angle)
            local0 = SVector(srcrad * ct * ca, srcrad * st * ca, srcrad * sa)
            vtx = srcpos + _m_mxv(trans, local0)
            dir0 = _m_mxv(trans, local0)
            u = _m_hat(dir0)
            dir = _m_dot(vtx, u) < 0.0 ? -u : u
            v = SVector(a * dir[1], b * dir[2], c * dir[3])
            lambda = 1.0 / _m_norm(v)
            e = SVector(a * v[1], b * v[2], c * v[3])
            pt = lambda * e
            offset = pt - vtx
            u1, m1 = _m_unorm(dir)
            u2, m2 = _m_unorm(offset)
            sep = (m1 == 0.0 || m2 == 0.0) ? 0.0 : _m_usep(u1, u2)
            angerr = sep - halfpi
            s = _m_dot(e, dir) >= 0.0 ? 1.0 : -1.0
            angle = umbral ? angle + s * angerr : angle - s * angerr
        end
        out[i] = pt
    end
    return out
end

"Radii of `body` with the checks of `zzgftreb_`."
function _m_radii(ks::KernelSet, body::Integer)
    v = get(ks.pool, string("BODY", body, "_RADII"), nothing)
    v === nothing && throw(SpiceError(:kernelvarnotfound,
        "the variable BODY$(body)_RADII is not in the kernel pool"))
    v isa Vector{Float64} || throw(SpiceError(:typemismatch,
        "the variable BODY$(body)_RADII is not numeric"))
    length(v) == 3 || throw(SpiceError(length(v) > 3 ? :arraytoosmall : :invalidcount,
        "the variable BODY$(body)_RADII has $(length(v)) values; 3 are required"))
    all(>(0.0), v) || throw(SpiceError(:badaxislength, "a radius of body $body is not positive"))
    return v[1], v[2], v[3]
end

function _m_bodyid(ks::KernelSet, name::NameOrID, kind::Symbol)
    name isa Integer && return Int32(name)
    try
        return body_code(ks, name)
    catch e
        e isa SpiceError && e.kind === :unknownbody &&
            throw(SpiceError(kind, "the body $name has no ID code"))
        rethrow()
    end
end

"""
    edterm([ctx,] trmtyp, source, target, et, fixref, abcorr, obsrvr, npts)
        -> (trgepc, obspos::SVector{3}, trmpts::Vector{SVector{3,Float64}})

Return `npts` points of the terminator of an ellipsoidal body that a spherical body
lights, as `edterm_c` of CSPICE.

`trmtyp` is `"UMBRAL"` or `"PENUMBRAL"`. `source` is the light source and `target` the
lit body. `et` is the epoch at the observer `obsrvr`, in seconds past J2000 TDB.
`fixref` is a body-fixed frame of the target. `abcorr` is the aberration correction.
Bodies and frames are names or IDs.

- `trgepc`: the epoch at the target, `et` corrected for light time if `abcorr` asks
  for it.
- `obspos`: the position of the observer relative to the target center, in `fixref`,
  km.
- `trmpts`: the terminator points relative to the target center, in `fixref`, km.

The call needs SPK data for the three bodies, the orientation of `fixref`, and the
`BODY<id>_RADII` of the source and the target.

Throws:
- `SpiceError(:notranslation)` if a body or frame name is not known.
- `SpiceError(:invalidfixref)` if `fixref` does not have its center at the target.
- `SpiceError(:spkinvalidoption)` if `abcorr` is not known.
- `SpiceError(:notsupported)` for another terminator type.
- `SpiceError(:invalidsize)` if `npts < 1`.
- `SpiceError(:objectstooclose)` if the source touches the target.

# Examples
This call needs kernels:
```julia
trgepc, obspos, trmpts = edterm("UMBRAL", "SUN", "MOON", et, "IAU_MOON", "LT+S", "EARTH", 3)
```

See also [`edlimb`](@ref), [`termpt`](@ref).
"""
function edterm(ks::KernelSet, trmtyp::AbstractString, source::NameOrID, target::NameOrID,
        et::Real, fixref::NameOrID, abcorr::AbstractString, obsrvr::NameOrID, npts::Integer)
    frcode = fixref isa Integer ? Int(fixref) : namfrm(ks, fixref)
    frcode == 0 && throw(SpiceError(:notranslation, "the frame $fixref has no frame ID code"))
    fd = frame_def(ks, frcode)
    trgid = _m_bodyid(ks, target, :notranslation)
    srcid = _m_bodyid(ks, source, :notranslation)
    obsid = _m_bodyid(ks, obsrvr, :notranslation)
    fd.center == trgid || throw(SpiceError(:invalidfixref,
        "the frame $fixref is not centered on the target body $target"))
    trgrad = _m_radii(ks, trgid)
    srcrad = _m_radii(ks, srcid)
    r = max(srcrad[1], srcrad[2], srcrad[3])
    t0 = Float64(et)
    corr = try
        Aberration(abcorr)
    catch e
        e isa SpiceError || rethrow()
        throw(SpiceError(:spkinvalidoption, "the aberration correction $abcorr is not recognized"))
    end
    trgpos, lttarg = spkpos(ks, Int(trgid), t0, frcode, abcorr, Int(obsid))
    trgepc = corr.uselt ? (corr.transmit ? t0 + lttarg : t0 - lttarg) : t0
    srcpos, _ = spkpos(ks, Int(srcid), trgepc, frcode, abcorr, Int(trgid))
    kind = uppercase(strip(trmtyp))
    kind in ("UMBRAL", "PENUMBRAL") || throw(SpiceError(:notsupported,
        "the terminator type $trmtyp is not supported; use UMBRAL or PENUMBRAL"))
    trmpts = _m_edterm(kind == "UMBRAL", trgrad[1], trgrad[2], trgrad[3], r, srcpos, Int(npts))
    return trgepc, -trgpos, trmpts
end
edterm(ctx::SpiceContext, trmtyp::AbstractString, source::NameOrID, target::NameOrID, et::Real,
    fixref::NameOrID, abcorr::AbstractString, obsrvr::NameOrID, npts::Integer) =
    edterm(kernelset(ctx), trmtyp, source, target, et, fixref, abcorr, obsrvr, npts)
edterm(trmtyp::AbstractString, source::NameOrID, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbstractString, obsrvr::NameOrID, npts::Integer) =
    edterm(default_context(), trmtyp, source, target, et, fixref, abcorr, obsrvr, npts)
