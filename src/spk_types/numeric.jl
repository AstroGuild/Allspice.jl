# Numeric helpers of the SPK evaluators: elementary functions that `Taylor` does not
# have yet, the two functions that carry derivatives through an iterative solution,
# and the small vector routines of CSPICE in their own order of operations.

# --- Taylor series of sinh, cosh and atan ---------------------------------------------

"""
    _sinhcosh(x::Taylor) -> (sinh(x), cosh(x))

Coupled recurrences `s_k = (1/k) Σ_{j=1}^{k} j x_j c_{k-j}` and
`c_k = (1/k) Σ_{j=1}^{k} j x_j s_{k-j}`, started from `sinh(x_0)` and `cosh(x_0)`.
"""
@generated function _sinhcosh(x::Taylor{N,T}) where {N,T<:Real}
    ex = _unpack(:a, :x, N)
    push!(ex, :(s0 = sinh(a0)), :(c0 = cosh(a0)))
    for k in 1:N
        ss = _sum([:($j * $(_sym(:a, j)) * $(_sym(:c, k - j))) for j in 1:k])
        cs = _sum([:($j * $(_sym(:a, j)) * $(_sym(:s, k - j))) for j in 1:k])
        push!(ex, :($(_sym(:s, k)) = $ss / $k), :($(_sym(:c, k)) = $cs / $k))
    end
    return Expr(:block, Expr(:meta, :inline), ex..., :(($(_pack(:s, N)), $(_pack(:c, N)))))
end

@inline Base.sinh(x::Taylor) = _sinhcosh(x)[1]
@inline Base.cosh(x::Taylor) = _sinhcosh(x)[2]

"""
    atan(y::Taylor, x::Taylor)

Series of the angle of the point `(x, y)`: `r_0 = atan(y_0, x_0)` and
`r' = (x y' - y x') / (x² + y²)`, by the quotient recurrence on the coefficients of
the derivative. Requires `x_0² + y_0² > 0` (not checked).
"""
@generated function Base.atan(y::Taylor{N,T}, x::Taylor{N,T}) where {N,T<:Real}
    ex = vcat(_unpack(:s, :y, N), _unpack(:c, :x, N))
    push!(ex, :(r0 = atan(s0, c0)))
    for m in 0:N-1
        # q_m: coefficient m of x² + y²; d: coefficient m of x y' - y x';
        # e_m: coefficient m of the derivative of the angle.
        q = _sum([:($(_sym(:c, j)) * $(_sym(:c, m - j)) + $(_sym(:s, j)) * $(_sym(:s, m - j)))
                  for j in 0:m])
        d = _sum([:($(j + 1) * ($(_sym(:c, m - j)) * $(_sym(:s, j + 1)) -
                                $(_sym(:s, m - j)) * $(_sym(:c, j + 1)))) for j in 0:m])
        push!(ex, :($(_sym(:q, m)) = $q))
        if m == 0
            push!(ex, :(e0 = $d / q0))
        else
            sub = _sum([:($(_sym(:q, j)) * $(_sym(:e, m - j))) for j in 1:m])
            push!(ex, :($(_sym(:e, m)) = ($d - $sub) / q0))
        end
        push!(ex, :($(_sym(:r, m + 1)) = $(_sym(:e, m)) / $(m + 1)))
    end
    return Expr(:block, Expr(:meta, :inline), ex..., _pack(:r, N))
end

# --- derivatives through an iterative solution ----------------------------------------
#
# CSPICE solves the two Kepler equations (`prop2b`, `kpsolv`) by bisection, which
# carries no derivative. The evaluators here make the solution with plain `Float64`
# numbers exactly as CSPICE does, and then make Newton steps in the epoch number type
# on the residual with its plain value removed. Such a step does not change the plain
# value of the unknown, and each step makes at least one more derivative order exact;
# `_ad_order` steps are sufficient.

"""
    _ad_order(::Type{T}) -> Int

Number of derivative orders that a value of the number type `T` carries: 0 for a
plain float, `N` plus the order of the coefficient type for `Taylor{N}`, one more
than the value type for a `ForwardDiff.Dual` (package extension). Another `Real`
type gives 8, which is sufficient for a type that carries not more than 8 orders.
"""
@inline _ad_order(::Type{<:AbstractFloat}) = 0
@inline _ad_order(::Type{Taylor{N,T}}) where {N,T} = N + _ad_order(T)
@inline _ad_order(::Type{<:Real}) = 8

"""
    _with_primal(x::T, p::Float64) -> T

`x` with its plain value (`primal(x)`) replaced by `p`; the derivative parts do not
change. Postcondition: `primal(result) == p`.
"""
@inline _with_primal(x::AbstractFloat, p::Float64) = oftype(x, p)
@inline _with_primal(x::Taylor, p::Float64) = _taylor(_with_primal(x.v, p), x.d)
@inline _with_primal(x::Real, p::Float64) = x + (p - primal(x))

# --- vector routines of CSPICE ----------------------------------------------------------

"The value of `pi_()` of CSPICE as a plain double."
const _PI = Float64(pi)

"Dot product in the order of `vdot`."
@inline _dot(a::Vec3, b::Vec3) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

"Cross product in the order of `vcrss`."
@inline _cross(a::Vec3, b::Vec3) = Vec3(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])

@inline _maxabs(v::Vec3) = max(abs(primal(v[1])), abs(primal(v[2])), abs(primal(v[3])))

"Euclidean norm of a plain vector, scaled by its largest component as `vnorm` does."
@inline function _vnorm(v::Vec3{Float64})
    m = _maxabs(v)
    m == 0.0 && return 0.0
    a, b, c = v[1] / m, v[2] / m, v[3] / m
    return m * sqrt(a * a + b * b + c * c)
end

"Unit vector along the plain vector `v`, or the zero vector for a zero `v` (`vhat`)."
@inline function _vhat(v::Vec3{Float64})
    m = _vnorm(v)
    return m > 0.0 ? Vec3{Float64}(v[1] / m, v[2] / m, v[3] / m) : zero(Vec3{Float64})
end

"""
    _vsep(a::Vec3{Float64}, b::Vec3{Float64}) -> Float64

Angle in `[0, pi]` between two plain vectors, by the half-angle formula of `vsep`;
zero if one of them is the zero vector.
"""
function _vsep(a::Vec3{Float64}, b::Vec3{Float64})
    (_vnorm(a) == 0.0 || _vnorm(b) == 0.0) && return 0.0
    u1, u2 = _vhat(a), _vhat(b)
    d = _dot(u1, u2)
    if d > 0.0
        return asin(_vnorm(u1 - u2) * 0.5) * 2.0
    elseif d < 0.0
        return _PI - asin(_vnorm(u1 + u2) * 0.5) * 2.0
    else
        return _PI / 2.0
    end
end

"""
    _vrotv(v::Vec3{T}, axis::Vec3{Float64}, theta::T) -> Vec3{T}

`v` rotated by the angle `theta` about `axis`, as `vrotv`: the part of `v` along the
axis stays, the part normal to it turns. A zero `axis` gives `v`.
"""
@inline function _vrotv(v::Vec3{T}, axis::Vec3{Float64}, theta::T) where {T}
    _vnorm(axis) == 0.0 && return v
    x = _vhat(axis)
    # Projection of `v` on `x`, with the scaling of `vproj`.
    biga = _maxabs(v)
    bigx = _maxabs(x)
    p = if biga == 0.0
        zero(Vec3{T})
    else
        r = x / bigx
        t = v / biga
        (_dot(t, r) * biga / _dot(r, r)) * r
    end
    v1 = v - p
    v2 = _cross(x, v1)
    s, c = sincos(theta)
    return (c * v1 + s * v2) + p
end

"`|a|` with the sign of `b` (`d_sign` of f2c: a zero `b` of either sign gives `+|a|`)."
@inline _dsign(a::Float64, b::Float64) = b >= 0.0 ? abs(a) : -abs(a)
