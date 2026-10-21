# Vector and matrix calls of CSPICE: three-dimensional forms, general-dimension (`g`)
# forms, derivatives of vector functions of a state, and frame construction.
#
# Each public function has the CSPICE name, returns its results instead of filling
# output arguments, and follows the arithmetic of the CSPICE source (the `_c` file where
# CSPICE has native C code, the f2c file otherwise), so that a `Float64` input gives the
# CSPICE result. Fixed-size results are `StaticArrays` values. The functions are generic
# over the number type; every discrete decision is made on `primal`.
#
# Every private helper of the files in `src/math/` has the prefix `_m_`.

# --- input conversion ---------------------------------------------------------------

@inline _m_f(x::Integer) = Float64(x)
@inline _m_f(x::AbstractIrrational) = Float64(x)
@inline _m_f(x::Real) = x

@noinline _m_baddim(what::String, n) =
    throw(SpiceError(:badarg, "$what must have $n elements"))

@inline function _m_v3(a::AbstractVector)
    length(a) == 3 || _m_baddim("a vector", 3)
    @inbounds return SVector(_m_f(a[1]), _m_f(a[2]), _m_f(a[3]))
end
@inline _m_v3(a::SVector{3,<:AbstractFloat}) = a

@inline function _m_v6(a::AbstractVector)
    length(a) == 6 || _m_baddim("a state vector", 6)
    @inbounds return SVector(_m_f(a[1]), _m_f(a[2]), _m_f(a[3]), _m_f(a[4]), _m_f(a[5]),
        _m_f(a[6]))
end

@inline function _m_v4(a::AbstractVector)
    length(a) == 4 || _m_baddim("a quaternion", 4)
    @inbounds return SVector(_m_f(a[1]), _m_f(a[2]), _m_f(a[3]), _m_f(a[4]))
end

@inline function _m_m3(m::AbstractMatrix)
    size(m) == (3, 3) || _m_baddim("a matrix", "3 by 3")
    @inbounds return SMatrix{3,3}(_m_f(m[1, 1]), _m_f(m[2, 1]), _m_f(m[3, 1]), _m_f(m[1, 2]),
        _m_f(m[2, 2]), _m_f(m[3, 2]), _m_f(m[1, 3]), _m_f(m[2, 3]), _m_f(m[3, 3]))
end
@inline _m_m3(m::SMatrix{3,3,<:AbstractFloat,9}) = m

@inline function _m_m6(m::AbstractMatrix)
    size(m) == (6, 6) || _m_baddim("a state transformation", "6 by 6")
    return SMatrix{6,6}(ntuple(k -> @inbounds(_m_f(m[k])), Val(36)))
end
@inline _m_m6(m::SMatrix{6,6,<:AbstractFloat,36}) = m

@inline _m_pv(s::SVector{6}) = SVector(s[1], s[2], s[3]), SVector(s[4], s[5], s[6])

@inline _m_state(p::SVector{3}, v::SVector{3}) = SVector(p[1], p[2], p[3], v[1], v[2], v[3])

@inline _m_maxabs(v::SVector{3}) =
    max(abs(primal(v[1])), abs(primal(v[2])), abs(primal(v[3])))

@inline _m_iszero(v::SVector{3}) =
    primal(v[1]) == 0.0 && primal(v[2]) == 0.0 && primal(v[3]) == 0.0

# --- private kernels (arguments are static vectors and matrices) -----------------------

@inline _m_dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]

@inline _m_cross(a, b) = SVector(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])

"Norm with the scaling by the largest component of `vnorm`."
@inline function _m_norm(v::SVector{3,T}) where {T}
    vmax = _m_maxabs(v)
    vmax == 0.0 && return zero(T)
    a = v[1] / vmax
    b = v[2] / vmax
    c = v[3] / vmax
    return vmax * sqrt(a * a + b * b + c * c)
end

"`(unit vector, norm)` as `unorm`: the zero vector for a zero input."
@inline function _m_unorm(v::SVector{3})
    vmag = _m_norm(v)
    primal(vmag) > 0.0 || return zero(v), vmag
    return SVector(v[1] / vmag, v[2] / vmag, v[3] / vmag), vmag
end

@inline _m_hat(v::SVector{3}) = _m_unorm(v)[1]

"""
Projection of `a` on `b`. With `recip` the scaling is a product with the reciprocal of
the largest component (`vproj_c`); without it the scaling is a division (`vproj_`).
"""
@inline function _m_proj(a::SVector{3}, b::SVector{3}, recip::Bool)
    biga = _m_maxabs(a)
    bigb = _m_maxabs(b)
    (biga == 0.0 || bigb == 0.0) && return zero(a) - zero(b)
    t = recip ? (1.0 / biga) * a : a / biga
    r = recip ? (1.0 / bigb) * b : b / bigb
    scale = _m_dot(t, r) * biga / _m_dot(r, r)
    return scale * r
end

"""
Component of `a` perpendicular to `b`, with the scaling choice of `_m_proj`. A zero `b`
gives `a`, as `vperp_` of SPICELIB and the text of `vperp_c` say.
"""
@inline function _m_perp(a::SVector{3}, b::SVector{3}, recip::Bool)
    biga = _m_maxabs(a)
    bigb = _m_maxabs(b)
    biga == 0.0 && return zero(a) - zero(b)
    bigb == 0.0 && return a - zero(b)
    t = recip ? (1.0 / biga) * a : a / biga
    r = recip ? (1.0 / bigb) * b : b / bigb
    # `vproj_c` and `vproj_` scale their arguments again; t and r have largest
    # component 1, so that scaling changes nothing.
    v = _m_proj(t, r, recip)
    return biga * (t - v)
end

"Unit cross product as `ucrss`: each vector is scaled by its largest component first."
@inline function _m_ucrss(v1::SVector{3}, v2::SVector{3})
    m1 = _m_maxabs(v1)
    m2 = _m_maxabs(v2)
    t1 = m1 != 0.0 ? v1 / m1 : zero(v1)
    t2 = m2 != 0.0 ? v2 / m2 : zero(v2)
    return _m_hat(_m_cross(t1, t2))
end

"Rotation of `v` about `axis` by `theta`, as `vrotv_`."
@inline function _m_vrotv(v::SVector{3}, axis::SVector{3}, theta)
    primal(_m_norm(axis)) == 0.0 && return v - zero(axis) * zero(theta)
    x = _m_hat(axis)
    p = _m_proj(v, x, false)
    v1 = v - p
    v2 = _m_cross(x, v1)
    s, c = sincos(theta)
    return (c * v1 + s * v2) + p
end

"""
    _m_mxv(m, v) -> SVector{3}

Contract: the product of the 3 by 3 matrix `m` and the vector `v`, each component the
sum of its three products from left to right, as `mxv_c`. The product `m * v` of
StaticArrays lets the compiler join a multiplication and an addition into one fused
operation, and compilers differ in where they do so; the result of this function does
not depend on the compiler, which matters where a component is a small difference of
large products.
The forms `_m_mtxv`, `_m_mxm`, `_m_mtxm` and `_m_mxmt` obey the same rule; the package uses
them in place of `*` for each 3 by 3 product.
"""
@inline _m_mxv(m::SMatrix{3,3}, v::SVector{3}) = SVector(
    m[1, 1] * v[1] + m[1, 2] * v[2] + m[1, 3] * v[3],
    m[2, 1] * v[1] + m[2, 2] * v[2] + m[2, 3] * v[3],
    m[3, 1] * v[1] + m[3, 2] * v[2] + m[3, 3] * v[3])

@inline _m_mtxv(m::SMatrix{3,3}, v::SVector{3}) = SVector(
    m[1, 1] * v[1] + m[2, 1] * v[2] + m[3, 1] * v[3],
    m[1, 2] * v[1] + m[2, 2] * v[2] + m[3, 2] * v[3],
    m[1, 3] * v[1] + m[2, 3] * v[2] + m[3, 3] * v[3])

@inline _m_rowcol(k::Int) = (k - 1) % 3 + 1, (k - 1) ÷ 3 + 1

@inline _m_mxm(a::SMatrix{3,3}, b::SMatrix{3,3}) = SMatrix{3,3}(ntuple(Val(9)) do k
    i, j = _m_rowcol(k)
    a[i, 1] * b[1, j] + a[i, 2] * b[2, j] + a[i, 3] * b[3, j]
end)

@inline _m_mtxm(a::SMatrix{3,3}, b::SMatrix{3,3}) = SMatrix{3,3}(ntuple(Val(9)) do k
    i, j = _m_rowcol(k)
    a[1, i] * b[1, j] + a[2, i] * b[2, j] + a[3, i] * b[3, j]
end)

@inline _m_mxmt(a::SMatrix{3,3}, b::SMatrix{3,3}) = SMatrix{3,3}(ntuple(Val(9)) do k
    i, j = _m_rowcol(k)
    a[i, 1] * b[j, 1] + a[i, 2] * b[j, 2] + a[i, 3] * b[j, 3]
end)

@inline _m_xpose(m::SMatrix{3,3}) = SMatrix{3,3}(m[1, 1], m[1, 2], m[1, 3], m[2, 1], m[2, 2],
    m[2, 3], m[3, 1], m[3, 2], m[3, 3])

@inline _m_col(m::SMatrix{3,3}, j::Int) = SVector(m[1, j], m[2, j], m[3, j])
@inline _m_row(m::SMatrix{3,3}, i::Int) = SVector(m[i, 1], m[i, 2], m[i, 3])

@inline _m_cols(a::SVector{3}, b::SVector{3}, c::SVector{3}) =
    SMatrix{3,3}(a[1], a[2], a[3], b[1], b[2], b[3], c[1], c[2], c[3])

@inline _m_rows(a::SVector{3}, b::SVector{3}, c::SVector{3}) =
    SMatrix{3,3}(a[1], b[1], c[1], a[2], b[2], c[2], a[3], b[3], c[3])

"Determinant with the operation order of `det_c`."
@inline _m_det(m::SMatrix{3,3}) =
    (m[1, 1] * (m[2, 2] * m[3, 3] - m[3, 2] * m[2, 3])) -
    (m[1, 2] * (m[2, 1] * m[3, 3] - m[3, 1] * m[2, 3])) +
    (m[1, 3] * (m[2, 1] * m[3, 2] - m[3, 1] * m[2, 2]))

"The index triple `(i1, i2, i3)` of `rotate_` for any integer axis (taken modulo 3)."
@inline function _m_axes(iaxis::Integer)
    t = mod(Int(iaxis), 3)
    return t == 1 ? (1, 2, 3) : t == 2 ? (2, 3, 1) : (3, 1, 2)
end

@inline function _m_rotate(angle::T, iaxis::Integer) where {T<:Real}
    s, c = sincos(angle)
    i1, i2, i3 = _m_axes(iaxis)
    z = zero(s)
    o = one(s)
    return SMatrix{3,3}(ntuple(Val(9)) do k
        i, j = _m_rowcol(k)
        (i == i1 || j == i1) ? (i == j ? o : z) :
        i == j ? c : (i == i2 ? s : -s)
    end)
end

@inline function _m_rotmat(m::SMatrix{3,3}, angle::Real, iaxis::Integer)
    s, c = sincos(angle)
    i1, i2, i3 = _m_axes(iaxis)
    return SMatrix{3,3}(ntuple(Val(9)) do k
        i, j = _m_rowcol(k)
        i == i1 ? m[i1, j] - zero(c) :
        i == i2 ? c * m[i2, j] + s * m[i3, j] : -s * m[i2, j] + c * m[i3, j]
    end)
end

"True if `m` is a rotation within the tolerances, as `isrot_c` (no argument checks)."
@inline function _m_isrot(m::SMatrix{3,3}, ntol::Real, dtol::Real)
    u1, n1 = _m_unorm(_m_col(m, 1))
    u2, n2 = _m_unorm(_m_col(m, 2))
    u3, n3 = _m_unorm(_m_col(m, 3))
    inside(x, tol) = (p = primal(x); 1.0 - tol <= p <= 1.0 + tol)
    normok = inside(n1, ntol) && inside(n2, ntol) && inside(n3, ntol)
    # `isrot_c` takes the determinant of the matrix whose rows are the unit columns.
    return normok && inside(_m_det(_m_rows(u1, u2, u3)), dtol)
end

"""
Unit state of `dvhat`: the unit position and its derivative. With `recip` the
perpendicular part of the velocity is scaled by a product with a reciprocal, as the
C function `dvhat_c` does (through `vperp_c`). Without it the scaling is a division,
as in the f2c routine `dvhat_` (through `vperp_`), which is the routine that `ducrss_`,
`dvsep_`, `zztwovxf_` and `surfpv_` call. The two forms differ by rounding, and the
reciprocal of a velocity below about 5.6e-309 is infinite, so only the division keeps
such a velocity finite. Postcondition: a zero position gives the zero position and the
velocity as it is.
"""
@inline function _m_dvhat(s::SVector{6}, recip::Bool=true)
    p, v = _m_pv(s)
    u, len = _m_unorm(p)
    primal(len) == 0.0 && return _m_state(u, v - zero(u))
    return _m_state(u, (1.0 / len) * _m_perp(v, u, recip))
end

"Cross product of the positions of two states and its derivative, as `dvcrss_c`."
@inline function _m_dvcrss(s1::SVector{6}, s2::SVector{6})
    p1, v1 = _m_pv(s1)
    p2, v2 = _m_pv(s2)
    return _m_state(_m_cross(p1, p2), _m_cross(v1, p2) + _m_cross(p1, v2))
end

"Unit cross product of two states and its derivative, as `ducrss_`."
@inline function _m_ducrss(s1::SVector{6}, s2::SVector{6})
    f1 = _m_maxabs(_m_pv(s1)[1])
    f2 = _m_maxabs(_m_pv(s2)[1])
    a = f1 > 0.0 ? (1.0 / f1) * s1 : s1
    b = f2 > 0.0 ? (1.0 / f2) * s2 : s2
    return _m_dvhat(_m_dvcrss(a, b), false)
end

"""
Quotient `a / b` with the checks of `zzdiv_`: `SpiceError(:dividebyzero)` for a zero
`b`, `SpiceError(:numericoverflow)` if the quotient is above the largest number, and
zero if it is below the smallest one.
"""
@inline function _m_div(a::Real, b::Real)
    a0 = primal(a)
    b0 = primal(b)
    b0 == 0.0 && throw(SpiceError(:dividebyzero, "numerical divide by zero"))
    a0 == 0.0 && return zero(a) * zero(b)
    expnt = trunc(log10(floatmax(Float64)))
    logn = log10(abs(a0))
    logd = log10(abs(b0))
    logn - logd > expnt && throw(SpiceError(:numericoverflow, "numerical overflow in a quotient"))
    logn - logd < -(expnt - 1.0) && return zero(a) * zero(b)
    return a / b
end

# --- three-dimensional vectors -----------------------------------------------------

"""
    vadd(v1, v2) -> SVector{3}

Return the sum `v1 + v2` of two vectors of 3 elements, as `vadd_c` of CSPICE.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

# Examples
```jldoctest
julia> vadd([1, 2, 3], [4, 5, 6]) == [5, 7, 9]
true
```

See also [`vsub`](@ref), [`vlcom`](@ref), [`vaddg`](@ref).
"""
vadd(v1::AbstractVector, v2::AbstractVector) = _m_v3(v1) + _m_v3(v2)

"""
    vsub(v1, v2) -> SVector{3}

Return the difference `v1 - v2` of two vectors of 3 elements, as `vsub_c` of CSPICE.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vadd`](@ref), [`vsubg`](@ref).
"""
vsub(v1::AbstractVector, v2::AbstractVector) = _m_v3(v1) - _m_v3(v2)

"""
    vscl(s, v1) -> SVector{3}

Return the product `s * v1` of a scalar and a vector of 3 elements, as `vscl_c` of
CSPICE.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

See also [`vlcom`](@ref), [`vsclg`](@ref).
"""
vscl(s::Real, v1::AbstractVector) = _m_f(s) * _m_v3(v1)

"""
    vminus(v1) -> SVector{3}

Return `-v1` for a vector of 3 elements, as `vminus_c` of CSPICE.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

See also [`vminug`](@ref).
"""
vminus(v1::AbstractVector) = -_m_v3(v1)

"""
    vdot(v1, v2)

Return the dot product of two vectors of 3 elements, as `vdot_c` of CSPICE.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vcrss`](@ref), [`dvdot`](@ref), [`vdotg`](@ref).
"""
vdot(v1::AbstractVector, v2::AbstractVector) = _m_dot(_m_v3(v1), _m_v3(v2))

"""
    vcrss(v1, v2) -> SVector{3}

Return the cross product `v1 × v2` of two vectors of 3 elements, as `vcrss_c` of CSPICE.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

# Examples
```jldoctest
julia> vcrss([1, 0, 0], [0, 1, 0]) == [0, 0, 1]
true
```

See also [`ucrss`](@ref), [`dvcrss`](@ref), [`vdot`](@ref).
"""
vcrss(v1::AbstractVector, v2::AbstractVector) = _m_cross(_m_v3(v1), _m_v3(v2))

"""
    vnorm(v1)

Return the Euclidean norm of a vector of 3 elements, as `vnorm_c` of CSPICE.

The components are divided by the one largest in magnitude before they are squared, so
the squares do not overflow or underflow. The result is zero only for the zero vector.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

# Examples
```jldoctest
julia> vnorm([3, 4, 0])
5.0

julia> isfinite(vnorm([3e200, 4e200, 0]))
true
```

See also [`vhat`](@ref), [`unorm`](@ref), [`vdist`](@ref), [`vnormg`](@ref).
"""
vnorm(v1::AbstractVector) = _m_norm(_m_v3(v1))

"""
    vhat(v1) -> SVector{3}

Return the unit vector along `v1`, as `vhat_c` of CSPICE. The zero vector gives the
zero vector.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

See also [`unorm`](@ref), [`dvhat`](@ref), [`vhatg`](@ref).
"""
vhat(v1::AbstractVector) = _m_hat(_m_v3(v1))

"""
    unorm(v1) -> (vout::SVector{3}, vmag)

Return the unit vector along `v1` and the norm of `v1`, as `unorm_c` of CSPICE. The
zero vector gives the zero vector and the norm 0.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

# Examples
```jldoctest
julia> unorm([3, 4, 0])
([0.6, 0.8, 0.0], 5.0)
```

See also [`vhat`](@ref), [`vnorm`](@ref), [`unormg`](@ref).
"""
unorm(v1::AbstractVector) = _m_unorm(_m_v3(v1))

"""
    vdist(v1, v2)

Return the distance `vnorm(v1 - v2)` between two points, as `vdist_c` of CSPICE.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vnorm`](@ref), [`vrel`](@ref), [`vdistg`](@ref).
"""
vdist(v1::AbstractVector, v2::AbstractVector) = _m_norm(_m_v3(v1) - _m_v3(v2))

"Angle between two unit vectors `u1` and `u2`, with the two-branch arcsine of `vsep`."
@inline function _m_usep(u1::SVector{3}, u2::SVector{3})
    d = primal(_m_dot(u1, u2))
    if d > 0.0
        return 2.0 * asin(0.5 * _m_norm(u1 - u2))
    elseif d < 0.0
        return Float64(pi) - 2.0 * asin(0.5 * _m_norm(u1 + u2))
    end
    return Float64(pi) / 2 - zero(_m_dot(u1, u2))
end

"""
    vsep(v1, v2)

Return the angle between two vectors of 3 elements, in radians in [0, pi], as `vsep_c`
of CSPICE.

The method keeps the precision near 0 and pi. The result is 0 if a vector is zero.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

# Examples
```jldoctest
julia> vsep([1, 0, 0], [0, 3, 0])
1.5707963267948966
```

See also [`dvsep`](@ref), [`vsepg`](@ref).
"""
function vsep(v1::AbstractVector, v2::AbstractVector)
    u1, m1 = _m_unorm(_m_v3(v1))
    u2, m2 = _m_unorm(_m_v3(v2))
    (primal(m1) == 0.0 || primal(m2) == 0.0) && return zero(asin(_m_dot(u1, u2)))
    return _m_usep(u1, u2)
end

"""
    vrel(v1, v2)

Return the relative difference `|v1 - v2| / max(|v1|, |v2|)` of two vectors of 3
elements, as `vrel_c` of CSPICE. The result is 0 if the vectors are equal, also if both
are zero.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vdist`](@ref), [`vrelg`](@ref).
"""
function vrel(v1::AbstractVector, v2::AbstractVector)
    a = _m_v3(v1)
    b = _m_v3(v2)
    nunorm = _m_norm(a - b)
    primal(nunorm) == 0.0 && return nunorm
    n1 = _m_norm(a)
    n2 = _m_norm(b)
    return nunorm / (primal(n1) >= primal(n2) ? n1 : n2)
end

"""
    vproj(a, b) -> SVector{3}

Return the projection of the vector `a` on the vector `b`, as `vproj_c` of CSPICE. The
result is the zero vector if `a` or `b` is zero.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

# Examples
```jldoctest
julia> vproj([1, 1, 0], [2, 0, 0]) == [1, 0, 0]
true

julia> vperp([1, 1, 0], [2, 0, 0]) == [0, 1, 0]
true
```

See also [`vperp`](@ref), [`vprjp`](@ref), [`vprojg`](@ref).
"""
vproj(a::AbstractVector, b::AbstractVector) = _m_proj(_m_v3(a), _m_v3(b), true)

"""
    vperp(a, b) -> SVector{3}

Return the component of the vector `a` that is perpendicular to the vector `b`, as
`vperp_c` of CSPICE.

For a zero `b` this call returns `a`, as the header of `vperp_c` says; the code of
`vperp_c` in CSPICE version N0067 returns the zero vector.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vproj`](@ref), which has an example.
"""
vperp(a::AbstractVector, b::AbstractVector) = _m_perp(_m_v3(a), _m_v3(b), true)

"""
    vlcom(a, v1, b, v2) -> SVector{3}

Return the linear combination `a*v1 + b*v2` of two vectors of 3 elements, as `vlcom_c`
of CSPICE.

Throws `SpiceError(:badarg)` if a vector does not have 3 elements.

See also [`vlcom3`](@ref), [`vlcomg`](@ref).
"""
vlcom(a::Real, v1::AbstractVector, b::Real, v2::AbstractVector) =
    _m_f(a) * _m_v3(v1) + _m_f(b) * _m_v3(v2)

"""
    vlcom3(a, v1, b, v2, c, v3) -> SVector{3}

Return the linear combination `a*v1 + b*v2 + c*v3` of three vectors of 3 elements, as
`vlcom3_c` of CSPICE.

Throws `SpiceError(:badarg)` if a vector does not have 3 elements.

See also [`vlcom`](@ref).
"""
vlcom3(a::Real, v1::AbstractVector, b::Real, v2::AbstractVector, c::Real, v3::AbstractVector) =
    _m_f(a) * _m_v3(v1) + _m_f(b) * _m_v3(v2) + _m_f(c) * _m_v3(v3)

"""
    vrotv(v, axis, theta) -> SVector{3}

Return the vector `v` rotated about the vector `axis` by `theta` radians, as `vrotv_c`
of CSPICE.

The rotation follows the right-hand rule. `axis` need not have unit length; a zero
`axis` gives `v`.

Throws `SpiceError(:badarg)` if a vector does not have 3 elements.

# Examples
```jldoctest
julia> vrotv([1, 0, 0], [0, 0, 2], pi / 2) ≈ [0, 1, 0]
true
```

See also [`axisar`](@ref), [`rotvec`](@ref).
"""
vrotv(v::AbstractVector, axis::AbstractVector, theta::Real) =
    _m_vrotv(_m_v3(v), _m_v3(axis), _m_f(theta))

"""
    vzero(v) -> Bool

Return `true` if all 3 components of `v` are zero, as `vzero_c` of CSPICE.

Throws `SpiceError(:badarg)` if `v` does not have 3 elements.

See also [`vzerog`](@ref).
"""
vzero(v::AbstractVector) = _m_iszero(_m_v3(v))

"""
    vequ(vin) -> SVector{3}

Return the vector `vin` of 3 elements as an `SVector{3}`, as `vequ_c` of CSPICE.

Throws `SpiceError(:badarg)` if `vin` does not have 3 elements.

See also [`vpack`](@ref), [`vequg`](@ref).
"""
vequ(vin::AbstractVector) = _m_v3(vin)

"""
    vpack(x, y, z) -> SVector{3}

Return the vector with the components `x`, `y`, `z`, as `vpack_c` of CSPICE.

See also [`vupack`](@ref).
"""
vpack(x::Real, y::Real, z::Real) = SVector(_m_f(x), _m_f(y), _m_f(z))

"""
    vupack(v) -> (x, y, z)

Return the three components of the vector `v` as a tuple, as `vupack_c` of CSPICE.

Throws `SpiceError(:badarg)` if `v` does not have 3 elements.

See also [`vpack`](@ref).
"""
vupack(v::AbstractVector) = Tuple(_m_v3(v))

"""
    ucrss(v1, v2) -> SVector{3}

Return the unit vector along the cross product `v1 × v2`, as `ucrss_c` of CSPICE.

Each vector is scaled by its component largest in magnitude before the product, so the
product does not overflow. The result is the zero vector if the vectors are parallel or
one of them is zero.

Throws `SpiceError(:badarg)` if an argument does not have 3 elements.

See also [`vcrss`](@ref), [`ducrss`](@ref).
"""
ucrss(v1::AbstractVector, v2::AbstractVector) = _m_ucrss(_m_v3(v1), _m_v3(v2))

"""
    vtmv(v1, matrix, v2)

Return the bilinear form `v1' * matrix * v2` of two vectors of 3 elements and a 3 by 3
matrix, as `vtmv_c` of CSPICE.

Throws `SpiceError(:badarg)` if a vector does not have 3 elements or the matrix is not
3 by 3.

See also [`vtmvg`](@ref).
"""
function vtmv(v1::AbstractVector, matrix::AbstractMatrix, v2::AbstractVector)
    a = _m_v3(v1)
    m = _m_m3(matrix)
    b = _m_v3(v2)
    val = zero(a[1] * m[1, 1] * b[1])
    for k in 1:3, l in 1:3        # invariant: val is the sum of the terms before (k, l)
        val += a[k] * m[k, l] * b[l]
    end
    return val
end

# --- three-dimensional matrices ------------------------------------------------------

"""
    mxv(m, vin) -> SVector{3}

Return the product `m * vin` of a 3 by 3 matrix and a vector of 3 elements, as `mxv_c`
of CSPICE.

Throws `SpiceError(:badarg)` if `m` is not 3 by 3 or `vin` does not have 3 elements.

# Examples
```jldoctest
julia> m = [0 1 0; -1 0 0; 0 0 1];

julia> mxv(m, [1, 2, 3]) == [2, -1, 3]
true

julia> mtxv(m, [1, 2, 3]) == [-2, 1, 3]
true
```

See also [`mtxv`](@ref), [`mxm`](@ref), [`mxvg`](@ref).
"""
mxv(m::AbstractMatrix, vin::AbstractVector) = _m_mxv(_m_m3(m), _m_v3(vin))

"""
    mtxv(m, vin) -> SVector{3}

Return the product of the transpose of the 3 by 3 matrix `m` and the vector `vin`, as
`mtxv_c` of CSPICE. For a rotation matrix this is the inverse rotation.

Throws `SpiceError(:badarg)` if `m` is not 3 by 3 or `vin` does not have 3 elements.

See also [`mxv`](@ref), which has an example, and [`mtxvg`](@ref).
"""
mtxv(m::AbstractMatrix, vin::AbstractVector) = _m_mtxv(_m_m3(m), _m_v3(vin))

"""
    mxm(m1, m2) -> SMatrix{3,3}

Return the product `m1 * m2` of two 3 by 3 matrices, as `mxm_c` of CSPICE.

Throws `SpiceError(:badarg)` if a matrix is not 3 by 3.

See also [`mtxm`](@ref), [`mxmt`](@ref), [`mxv`](@ref), [`mxmg`](@ref).
"""
mxm(m1::AbstractMatrix, m2::AbstractMatrix) = _m_mxm(_m_m3(m1), _m_m3(m2))

"""
    mtxm(m1, m2) -> SMatrix{3,3}

Return the product of the transpose of `m1` and `m2`, for 3 by 3 matrices, as `mtxm_c`
of CSPICE.

Throws `SpiceError(:badarg)` if a matrix is not 3 by 3.

See also [`mxm`](@ref), [`mxmt`](@ref), [`mtxmg`](@ref).
"""
mtxm(m1::AbstractMatrix, m2::AbstractMatrix) = _m_mtxm(_m_m3(m1), _m_m3(m2))

"""
    mxmt(m1, m2) -> SMatrix{3,3}

Return the product of `m1` and the transpose of `m2`, for 3 by 3 matrices, as `mxmt_c`
of CSPICE.

Throws `SpiceError(:badarg)` if a matrix is not 3 by 3.

See also [`mxm`](@ref), [`mtxm`](@ref), [`mxmtg`](@ref).
"""
mxmt(m1::AbstractMatrix, m2::AbstractMatrix) = _m_mxmt(_m_m3(m1), _m_m3(m2))

"""
    xpose(m1) -> SMatrix{3,3}

Return the transpose of a 3 by 3 matrix, as `xpose_c` of CSPICE.

Throws `SpiceError(:badarg)` if `m1` is not 3 by 3.

See also [`xpose6`](@ref), [`xposeg`](@ref).
"""
xpose(m1::AbstractMatrix) = _m_xpose(_m_m3(m1))

"""
    xpose6(m1) -> SMatrix{6,6}

Return the transpose of a 6 by 6 matrix, as `xpose6_c` of CSPICE.

Throws `SpiceError(:badarg)` if `m1` is not 6 by 6.

See also [`xpose`](@ref), [`xposeg`](@ref).
"""
function xpose6(m1::AbstractMatrix)
    m = _m_m6(m1)
    return SMatrix{6,6}(ntuple(k -> m[(k-1)÷6+1, (k-1)%6+1], Val(36)))
end

"""
    invstm(mat) -> SMatrix{6,6}

Return the inverse of a state transformation matrix, as `invstm_c` of CSPICE.

A state transformation `[R 0; dR R]`, with a rotation `R` and its time derivative `dR`,
as [`sxform`](@ref) returns it, has the inverse
`[R' 0; dR' R']`. The call transposes each of the four 3 by 3 blocks of `mat` in its
place. It does not check that `mat` is a state transformation.

Throws `SpiceError(:badarg)` if `mat` is not 6 by 6.

# Examples
```jldoctest
julia> x = eul2xf([0.1, 0.2, 0.3, 1e-3, 2e-3, 3e-3], 3, 1, 3);

julia> invstm(x) * x ≈ [i == j for i in 1:6, j in 1:6]
true
```

See also [`xpose6`](@ref), [`invert`](@ref).
"""
function invstm(mat::AbstractMatrix)
    m = _m_m6(mat)
    # element k in column order is (i, j); its source is the transposed place in its block
    return SMatrix{6,6}(ntuple(Val(36)) do k
        i, j = (k - 1) % 6, (k - 1) ÷ 6
        bi, bj = 3 * (i ÷ 3), 3 * (j ÷ 3)
        m[bi+j-bj+1, bj+i-bi+1]
    end)
end

"""
    det(m1)

Return the determinant of a 3 by 3 matrix, as `det_c` (expansion along the first row).
Not exported, because `LinearAlgebra` exports the same name: call `Allspice.det`.
"""
det(m1::AbstractMatrix) = _m_det(_m_m3(m1))

"""
    trace(matrix)

Return the trace of a 3 by 3 matrix, as `trace_c` of CSPICE.

Throws `SpiceError(:badarg)` if `matrix` is not 3 by 3.
"""
function trace(matrix::AbstractMatrix)
    m = _m_m3(matrix)
    return m[1, 1] + m[2, 2] + m[3, 3]
end

"""
    invert(m) -> SMatrix{3,3}

Return the inverse of a 3 by 3 matrix, as `invert_c` of CSPICE.

A matrix whose determinant is below `1e-16` in absolute value counts as singular: the
result is then the zero matrix, and no error is thrown.

Throws `SpiceError(:badarg)` if `m` is not 3 by 3.

# Examples
```jldoctest
julia> invert([2 0 0; 0 4 0; 0 0 1]) == [0.5 0 0; 0 0.25 0; 0 0 1]
true

julia> iszero(invert([1 2 3; 2 4 6; 0 0 1]))
true
```

See also [`invort`](@ref), [`xpose`](@ref).
"""
function invert(m1::AbstractMatrix)
    m = _m_m3(m1)
    mdet = _m_det(m)
    abs(primal(mdet)) < 1.0e-16 && return zero(m) * zero(mdet)
    invdet = 1.0 / mdet
    return invdet * SMatrix{3,3}(
        (m[2, 2] * m[3, 3] - m[3, 2] * m[2, 3]),
        -(m[2, 1] * m[3, 3] - m[3, 1] * m[2, 3]),
        (m[2, 1] * m[3, 2] - m[3, 1] * m[2, 2]),
        -(m[1, 2] * m[3, 3] - m[3, 2] * m[1, 3]),
        (m[1, 1] * m[3, 3] - m[3, 1] * m[1, 3]),
        -(m[1, 1] * m[3, 2] - m[3, 1] * m[1, 2]),
        (m[1, 2] * m[2, 3] - m[2, 2] * m[1, 3]),
        -(m[1, 1] * m[2, 3] - m[2, 1] * m[1, 3]),
        (m[1, 1] * m[2, 2] - m[2, 1] * m[1, 2]))
end

"""
    invort(m) -> SMatrix{3,3}

Return the inverse of a 3 by 3 matrix whose columns are orthogonal and have any
non-zero length, as `invort_c` of CSPICE.

The result is the transpose with each row divided by the squared length of the
matching column. The columns are not checked to be orthogonal.

Throws `SpiceError(:zerolengthcolumn)` if a column is zero. Throws
`SpiceError(:columntoosmall)` if the length of a column has no finite reciprocal.
Throws `SpiceError(:badarg)` if `m` is not 3 by 3.

See also [`invert`](@ref).
"""
function invort(m1::AbstractMatrix)
    m = _m_m3(m1)
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
    isrot(m, ntol, dtol) -> Bool

Return `true` if the 3 by 3 matrix `m` is a rotation matrix within two tolerances, as
`isrot_c` of CSPICE.

Each column must have a norm in `1 ± ntol`, and the matrix of the unit columns must
have a determinant in `1 ± dtol`.

Throws `SpiceError(:valueoutofrange)` if `ntol` or `dtol` is negative, and
`SpiceError(:badarg)` if `m` is not 3 by 3.

See also [`sharpr`](@ref).
"""
function isrot(m::AbstractMatrix, ntol::Real, dtol::Real)
    ntol >= 0 || throw(SpiceError(:valueoutofrange, "ntol should be non-negative; it is $ntol"))
    dtol >= 0 || throw(SpiceError(:valueoutofrange, "dtol should be non-negative; it is $dtol"))
    return _m_isrot(_m_m3(m), Float64(ntol), Float64(dtol))
end

"""
    ident() -> SMatrix{3,3,Float64}

Return the 3 by 3 identity matrix, as `ident_c` of CSPICE.
"""
ident() = SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)

"""
    rotate(angle, iaxis) -> SMatrix{3,3}

Return the matrix that rotates a coordinate system by `angle` radians about the axis
`iaxis` (1, 2, 3 for x, y, z), as `rotate_c` of CSPICE.

The product of the result and a vector gives the components of that vector in the
rotated system. The vector itself is not rotated: [`axisar`](@ref) and [`vrotv`](@ref)
do that, and their matrix is the transpose. As in CSPICE, any integer `iaxis` is taken
modulo 3, so 4 is x, 0 is z and -1 is y.

# Examples
```jldoctest
julia> Matrix(round.(rotate(pi / 2, 3); digits=12))
3×3 Matrix{Float64}:
  0.0  1.0  0.0
 -1.0  0.0  0.0
  0.0  0.0  1.0
```

See also [`rotmat`](@ref), [`rotvec`](@ref), [`eul2m`](@ref), [`axisar`](@ref).
"""
rotate(angle::Real, iaxis::Integer) = _m_rotate(_m_f(angle), iaxis)

"""
    rotmat(m1, angle, iaxis) -> SMatrix{3,3}

Return `rotate(angle, iaxis) * m1` for a 3 by 3 matrix `m1`, as `rotmat_c` of CSPICE.
This is the transformation `m1` followed by a rotation of the coordinate system by
`angle` radians about the axis `iaxis`.

Throws `SpiceError(:badarg)` if `m1` is not 3 by 3.

See also [`rotate`](@ref), [`rotvec`](@ref).
"""
rotmat(m1::AbstractMatrix, angle::Real, iaxis::Integer) = _m_rotmat(_m_m3(m1), _m_f(angle), iaxis)

"""
    rotvec(v1, angle, iaxis) -> SVector{3}

Return `rotate(angle, iaxis) * v1`, as `rotvec_c` of CSPICE: the components of the
vector `v1` in a coordinate system rotated by `angle` radians about the axis `iaxis`.

Throws `SpiceError(:badarg)` if `v1` does not have 3 elements.

# Examples
```jldoctest
julia> rotvec([1, 0, 0], pi / 2, 3) ≈ [0, -1, 0]
true
```

See also [`rotate`](@ref), [`rotmat`](@ref), [`vrotv`](@ref).
"""
function rotvec(v1::AbstractVector, angle::Real, iaxis::Integer)
    v = _m_v3(v1)
    s, c = sincos(_m_f(angle))
    i1, i2, i3 = _m_axes(iaxis)
    a = v[i1] - zero(c)
    b = c * v[i2] + s * v[i3]
    d = -s * v[i2] + c * v[i3]
    return i1 == 1 ? SVector(a, b, d) : i1 == 2 ? SVector(d, a, b) : SVector(b, d, a)
end

@noinline function _m_check_twovec(indexa::Integer, indexp::Integer)
    (max(indexp, indexa) > 3 || min(indexp, indexa) < 1) && throw(SpiceError(:badindex,
        "the definition indexes must be in 1:3; they are $indexa and $indexp"))
    indexa == indexp && throw(SpiceError(:undefinedframe,
        "the two definition indexes are both $indexa; they must be different"))
    return nothing
end

"""
    twovec(axdef, indexa, plndef, indexp) -> SMatrix{3,3}

Return the rotation matrix from the input frame to the right-handed frame that two
vectors define, as `twovec_c` of CSPICE.

Axis `indexa` of the new frame is along `axdef`. Axis `indexp` is in the plane of
`axdef` and `plndef`, on the side of `plndef`. The rows of the result are the axes of
the new frame, so the product with a vector gives its components in the new frame.

Throws:
- `SpiceError(:badindex)` if an index is not 1, 2 or 3.
- `SpiceError(:undefinedframe)` if the indexes are equal.
- `SpiceError(:dependentvectors)` if the vectors are parallel or one is zero.
- `SpiceError(:badarg)` if a vector does not have 3 elements.

# Examples
The frame with z along the input x and with x toward the input y:
```jldoctest
julia> Matrix(twovec([2, 0, 0], 3, [1, 1, 0], 1))
3×3 Matrix{Float64}:
 0.0  1.0  0.0
 0.0  0.0  1.0
 1.0  0.0  0.0
```

See also [`twovxf`](@ref).
"""
function twovec(axdef::AbstractVector, indexa::Integer, plndef::AbstractVector, indexp::Integer)
    _m_check_twovec(indexa, indexp)
    a = _m_v3(axdef)
    p = _m_v3(plndef)
    i1 = Int(indexa)
    i2 = i1 % 3 + 1
    r1 = _m_hat(a) - zero(p)
    if indexp == i2
        r3 = _m_ucrss(a, p)
        r2 = _m_ucrss(r3, a)
    else
        r2 = _m_ucrss(p, a)
        r3 = _m_ucrss(a, r2)
    end
    _m_iszero(r2) && throw(SpiceError(:dependentvectors,
        "the vectors axdef and plndef are linearly dependent"))
    return i1 == 1 ? _m_rows(r1, r2, r3) : i1 == 2 ? _m_rows(r3, r1, r2) : _m_rows(r2, r3, r1)
end

"""
    twovxf(axdef, indexa, plndef, indexp) -> SMatrix{6,6}

Return the 6 by 6 state transformation to the frame that two states define, as
`twovxf_c` of CSPICE.

The position parts of the states `axdef` and `plndef` define the frame as in
[`twovec`](@ref); the velocity parts give its time derivative. The upper left and lower
right 3 by 3 blocks of the result are the `twovec` matrix. The lower left block is its
time derivative. The upper right block is zero.

Throws as [`twovec`](@ref), and `SpiceError(:badarg)` if a state does not have 6
elements.
"""
function twovxf(axdef::AbstractVector, indexa::Integer, plndef::AbstractVector, indexp::Integer)
    _m_check_twovec(indexa, indexp)
    a = _m_v6(axdef)
    p = _m_v6(plndef)
    i1 = Int(indexa)
    i2 = i1 % 3 + 1
    s1 = _m_dvhat(a, false) - zero(p)
    if indexp == i2
        s3 = _m_ducrss(a, p)
        s2 = _m_ducrss(s3, a)
    else
        s2 = _m_ducrss(p, a)
        s3 = _m_ducrss(a, s2)
    end
    _m_iszero(_m_pv(s2)[1]) && throw(SpiceError(:dependentvectors,
        "the direction vectors of the states axdef and plndef are linearly dependent"))
    x, y, z = i1 == 1 ? (s1, s2, s3) : i1 == 2 ? (s3, s1, s2) : (s2, s3, s1)
    r = _m_rows(_m_pv(x)[1], _m_pv(y)[1], _m_pv(z)[1])
    d = _m_rows(_m_pv(x)[2], _m_pv(y)[2], _m_pv(z)[2])
    return _m_xform(r, d)
end

@inline function _m_xform(r::SMatrix{3,3}, d::SMatrix{3,3})
    z = zero(r[1, 1] + d[1, 1])
    return SMatrix{6,6}(ntuple(Val(36)) do k
        i = (k - 1) % 6 + 1
        j = (k - 1) ÷ 6 + 1
        j <= 3 ? (i <= 3 ? r[i, j] + z : d[i-3, j] + z) : (i <= 3 ? z : r[i-3, j-3] + z)
    end)
end

@inline function _m_blocks(x::SMatrix{6,6})
    r = SMatrix{3,3}(ntuple(k -> x[(k-1)%3+1, (k-1)÷3+1], Val(9)))
    d = SMatrix{3,3}(ntuple(k -> x[(k-1)%3+4, (k-1)÷3+1], Val(9)))
    return r, d
end

"""
    frame(x) -> (x̂, y, z)

Return a right-handed orthonormal frame whose first axis is along `x`, as `frame_c` of
CSPICE.

`x̂` is the unit vector along `x`, and `y` and `z` complete the frame; the three are
`SVector{3}` values. `y[i] == 0` for the first index `i` of the smallest
`abs(x̂[i])`. The zero vector gives the axes of the input frame. CSPICE overwrites its
argument `x` with `x̂`; this call does not change `x`.

Throws `SpiceError(:badarg)` if `x` does not have 3 elements.

# Examples
```jldoctest
julia> frame([0, 0, 2]) == ([0, 0, 1], [0, -1, 0], [1, 0, 0])
true
```

See also [`twovec`](@ref).
"""
function frame(x::AbstractVector)
    u = _m_hat(_m_v3(x))
    a = u[1] * u[1]
    b = u[2] * u[2]
    c = u[3] * u[3]
    o = zero(a)
    if primal(a + b + c) == 0.0
        return SVector(o + 1, o, o), SVector(o, o + 1, o), SVector(o, o, o + 1)
    end
    # The smallest component (s0) is zero in y; the other two (s1, s2) are swapped, with
    # the sign of the first changed: y[s1] = -x[s2]/f, y[s2] = x[s1]/f, z[s0] = f,
    # z[s1] = -x[s0] y[s2], z[s2] = x[s0] y[s1]. The cyclic order of (s0, s1, s2) makes
    # the frame right-handed.
    pa, pb, pc = primal(a), primal(b), primal(c)
    if pa <= pb && pa <= pc                     # (s0, s1, s2) = (1, 2, 3)
        f = sqrt(b + c)
        y = SVector(o, -u[3] / f, u[2] / f)
        return u, y, SVector(f, -u[1] * y[3], u[1] * y[2])
    elseif pb <= pa && pb <= pc                 # (2, 3, 1)
        f = sqrt(a + c)
        y = SVector(u[3] / f, o, -u[1] / f)
        return u, y, SVector(u[2] * y[3], f, -u[2] * y[1])
    end
    f = sqrt(a + b)                             # (3, 1, 2)
    y = SVector(-u[2] / f, u[1] / f, o)
    return u, y, SVector(-u[3] * y[2], u[3] * y[1], f)
end

# --- derivatives of vector functions of a state -------------------------------------

"""
    dvhat(s1) -> SVector{6}

Return the unit vector along the position of the state `s1`, followed by the time
derivative of that unit vector, as `dvhat_c` of CSPICE.

A state has 6 elements: a position and a velocity. A zero position gives a zero unit
vector and the velocity unchanged.

Throws `SpiceError(:badarg)` if `s1` does not have 6 elements.

# Examples
```jldoctest
julia> dvhat([2, 0, 0, 1, 2, 3]) == [1, 0, 0, 0, 1, 1.5]
true
```

See also [`vhat`](@ref), [`dvnorm`](@ref), [`ducrss`](@ref).
"""
dvhat(s1::AbstractVector) = _m_dvhat(_m_v6(s1))

"""
    dvdot(s1, s2)

Return the time derivative of the dot product of the positions of two states, as
`dvdot_c` of CSPICE.

Throws `SpiceError(:badarg)` if a state does not have 6 elements.

See also [`vdot`](@ref), [`dvcrss`](@ref).
"""
function dvdot(s1::AbstractVector, s2::AbstractVector)
    a = _m_v6(s1)
    b = _m_v6(s2)
    return a[1] * b[4] + a[2] * b[5] + a[3] * b[6] + a[4] * b[1] + a[5] * b[2] + a[6] * b[3]
end

"""
    dvcrss(s1, s2) -> SVector{6}

Return the cross product of the positions of two states, followed by its time
derivative, as `dvcrss_c` of CSPICE.

Throws `SpiceError(:badarg)` if a state does not have 6 elements.

See also [`vcrss`](@ref), [`ducrss`](@ref).
"""
dvcrss(s1::AbstractVector, s2::AbstractVector) = _m_dvcrss(_m_v6(s1), _m_v6(s2))

"""
    ducrss(s1, s2) -> SVector{6}

Return the unit vector along the cross product of the positions of two states, followed
by its time derivative, as `ducrss_c` of CSPICE.

Each state is scaled by its largest position component first. Parallel positions give
a zero unit vector. If a scaled velocity or a product of the scaled states is not
finite, the derivative part is NaN; CSPICE returns numbers without meaning there.

Throws `SpiceError(:badarg)` if a state does not have 6 elements.

See also [`ucrss`](@ref), [`dvcrss`](@ref), [`dvhat`](@ref).
"""
ducrss(s1::AbstractVector, s2::AbstractVector) = _m_ducrss(_m_v6(s1), _m_v6(s2))

"""
    dvnorm(state)

Return the time derivative of the norm of the position of `state`, as `dvnorm_c` of
CSPICE. A zero position gives 0.

Throws `SpiceError(:badarg)` if `state` does not have 6 elements.

See also [`vnorm`](@ref), [`dvhat`](@ref).
"""
function dvnorm(state::AbstractVector)
    p, v = _m_pv(_m_v6(state))
    u, len = _m_unorm(p)
    primal(len) == 0.0 && return zero(_m_dot(v, u))
    return _m_dot(v, u)
end

"""
    dvsep(s1, s2)

Return the time derivative of the angle between the positions of two states, in rad/s
for velocities per second, as `dvsep_c` of CSPICE. Parallel positions, or a zero position,
give 0.

Throws `SpiceError(:numericoverflow)` if the size of the result is above 1e308, and
`SpiceError(:badarg)` if a state does not have 6 elements.

See also [`vsep`](@ref).
"""
function dvsep(s1::AbstractVector, s2::AbstractVector)
    u1, d1 = _m_pv(_m_dvhat(_m_v6(s1), false))
    u2, d2 = _m_pv(_m_dvhat(_m_v6(s2), false))
    pcross = _m_cross(u1, u2)
    numr = _m_dot(u1, d2) + _m_dot(d1, u2)
    _m_iszero(pcross) && return zero(numr)
    return _m_div(-numr, _m_norm(pcross))
end

# --- general-dimension forms ---------------------------------------------------------

@inline function _m_samelen(a::AbstractVector, b::AbstractVector)
    length(a) == length(b) ||
        throw(SpiceError(:badarg, "the vectors have $(length(a)) and $(length(b)) elements"))
    return nothing
end

"""
    vaddg(v1, v2)

Return the sum `v1 + v2` of two vectors of equal length, as `vaddg_c` of CSPICE.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vadd`](@ref), [`vsubg`](@ref), [`vlcomg`](@ref).
"""
vaddg(v1::AbstractVector, v2::AbstractVector) = (_m_samelen(v1, v2); _m_f.(v1) .+ _m_f.(v2))

"""
    vsubg(v1, v2)

Return the difference `v1 - v2` of two vectors of equal length, as `vsubg_c` of
CSPICE.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vsub`](@ref), [`vaddg`](@ref).
"""
vsubg(v1::AbstractVector, v2::AbstractVector) = (_m_samelen(v1, v2); _m_f.(v1) .- _m_f.(v2))

"""
    vsclg(s, v1)

Return the product `s * v1` of a scalar and a vector of any length, as `vsclg_c` of
CSPICE.

See also [`vscl`](@ref).
"""
vsclg(s::Real, v1::AbstractVector) = _m_f(s) .* _m_f.(v1)

"""
    vminug(vin)

Return `-vin` for a vector of any length, as `vminug_c` of CSPICE.

See also [`vminus`](@ref).
"""
vminug(vin::AbstractVector) = .-_m_f.(vin)

"""
    vequg(vin)

Return a copy of the vector `vin` of any length, as `vequg_c` of CSPICE.

See also [`vequ`](@ref).
"""
vequg(vin::AbstractVector) = _m_f.(vin)

"""
    vdotg(v1, v2)

Return the dot product of two vectors of equal length, as `vdotg_c` of CSPICE.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vdot`](@ref).
"""
function vdotg(v1::AbstractVector, v2::AbstractVector)
    _m_samelen(v1, v2)
    dot = zero(_m_f(zero(eltype(v1))) * _m_f(zero(eltype(v2))))
    @inbounds for i in eachindex(v1, v2)       # invariant: dot is the sum of terms before i
        dot += _m_f(v1[i]) * _m_f(v2[i])
    end
    return dot
end

@inline function _m_maxabsg(v::AbstractVector)
    scale = 0.0
    @inbounds for x in v                       # invariant: scale bounds the elements seen
        scale = max(scale, abs(primal(x)))
    end
    return scale
end

"""
    vnormg(v1)

Return the Euclidean norm of a vector of any length, as `vnormg_c` of CSPICE. The
elements are divided by the one largest in magnitude before they are squared, as in
[`vnorm`](@ref).

See also [`unormg`](@ref), [`vdistg`](@ref).
"""
function vnormg(v1::AbstractVector)
    scale = _m_maxabsg(v1)
    norm = zero(_m_f(zero(eltype(v1)))) / 1.0
    scale == 0.0 && return norm
    @inbounds for x in v1                      # invariant: norm is the sum of squares seen
        t = _m_f(x) / scale
        norm += t * t
    end
    return scale * sqrt(norm)
end

"""
    vhatg(v1)

Return the unit vector along a vector of any length, as `vhatg_c` of CSPICE. The zero
vector gives the zero vector.

See also [`vhat`](@ref), [`unormg`](@ref).
"""
vhatg(v1::AbstractVector) = unormg(v1)[1]

"""
    unormg(v1) -> (vout, vmag)

Return the unit vector along a vector of any length and the norm of that vector, as
`unormg_c` of CSPICE. The zero vector gives the zero vector and the norm 0.

See also [`unorm`](@ref), [`vhatg`](@ref).
"""
function unormg(v1::AbstractVector)
    vmag = vnormg(v1)
    primal(vmag) > 0.0 || return _m_f.(v1) .* zero(vmag), vmag
    return _m_f.(v1) ./ vmag, vmag
end

"""
    vdistg(v1, v2)

Return the distance between two points of equal dimension, as `vdistg_c` of CSPICE.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vdist`](@ref).
"""
function vdistg(v1::AbstractVector, v2::AbstractVector)
    _m_samelen(v1, v2)
    scale = 0.0
    @inbounds for i in eachindex(v1, v2)       # invariant: scale bounds the differences seen
        scale = max(scale, abs(primal(_m_f(v1[i]) - _m_f(v2[i]))))
    end
    dist = zero(_m_f(zero(eltype(v1))) - _m_f(zero(eltype(v2)))) / 1.0
    scale == 0.0 && return dist
    @inbounds for i in eachindex(v1, v2)       # invariant: dist is the sum of squares seen
        t = (_m_f(v1[i]) - _m_f(v2[i])) / scale
        dist += t * t
    end
    return scale * sqrt(dist)
end

"""
    vsepg(v1, v2)

Return the angle between two vectors of equal length, in radians in [0, pi], with the
method of [`vsep`](@ref). The result is 0 if a vector is zero.

`vsepg_c` of CSPICE N0067 can lose up to half of the digits near 0 and pi, so the
results of the two differ there.

Throws `SpiceError(:badarg)` if the lengths differ.
"""
function vsepg(v1::AbstractVector, v2::AbstractVector)
    _m_samelen(v1, v2)
    mag1 = vnormg(v1)
    mag2 = vnormg(v2)
    (primal(mag1) == 0.0 || primal(mag2) == 0.0) && return zero(asin(zero(mag1 * mag2)))
    d = primal(vdotg(v1, v2))
    r1 = 1.0 / mag1
    r2 = 1.0 / mag2
    d == 0.0 && return Float64(pi) / 2 - zero(asin(zero(mag1 * mag2)))
    sgn = d > 0.0 ? -1.0 : 1.0
    acc = zero(mag1 * mag2)
    @inbounds for i in eachindex(v1, v2)       # invariant: acc is the sum of squares seen
        t = _m_f(v1[i]) * r1 + sgn * (_m_f(v2[i]) * r2)
        acc += t * t
    end
    half = 2.0 * asin(0.5 * sqrt(acc))
    return d > 0.0 ? half : Float64(pi) - half
end

"""
    vrelg(v1, v2)

Return the relative difference `|v1 - v2| / max(|v1|, |v2|)` of two vectors of equal
length, as `vrelg_c` of CSPICE. The result is 0 for equal vectors.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vrel`](@ref).
"""
function vrelg(v1::AbstractVector, v2::AbstractVector)
    nunorm = vdistg(v1, v2)
    primal(nunorm) == 0.0 && return nunorm
    n1 = vnormg(v1)
    n2 = vnormg(v2)
    return nunorm / (primal(n1) >= primal(n2) ? n1 : n2)
end

"""
    vprojg(a, b)

Return the projection of the vector `a` on the vector `b` of equal length, as
`vprojg_c` of CSPICE. A zero `b` gives the zero vector.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vproj`](@ref).
"""
function vprojg(a::AbstractVector, b::AbstractVector)
    adotb = vdotg(a, b)
    bdotb = vdotg(b, b)
    scale = primal(bdotb) == 0.0 ? zero(adotb / one(bdotb)) : adotb / bdotb
    return scale .* _m_f.(b)
end

"""
    vlcomg(a, v1, b, v2)

Return the linear combination `a*v1 + b*v2` of two vectors of equal length, as
`vlcomg_c` of CSPICE.

Throws `SpiceError(:badarg)` if the lengths differ.

See also [`vlcom`](@ref).
"""
vlcomg(a::Real, v1::AbstractVector, b::Real, v2::AbstractVector) =
    (_m_samelen(v1, v2); _m_f(a) .* _m_f.(v1) .+ _m_f(b) .* _m_f.(v2))

"""
    vzerog(v) -> Bool

Return `true` if every element of the vector `v` of any length is zero, as `vzerog_c`
of CSPICE.

See also [`vzero`](@ref).
"""
vzerog(v::AbstractVector) = all(x -> primal(x) == 0.0, v)

"""
    vtmvg(v1, matrix, v2)

Return the bilinear form `v1' * matrix * v2` for any dimensions that fit, as `vtmvg_c`
of CSPICE.

Throws `SpiceError(:badarg)` if `size(matrix) != (length(v1), length(v2))`.

See also [`vtmv`](@ref).
"""
function vtmvg(v1::AbstractVector, matrix::AbstractMatrix, v2::AbstractVector)
    size(matrix) == (length(v1), length(v2)) ||
        throw(SpiceError(:badarg, "the matrix size does not fit the vectors"))
    val = zero(_m_f(zero(eltype(v1))) * _m_f(zero(eltype(matrix))) * _m_f(zero(eltype(v2))))
    @inbounds for k in 1:length(v1), l in 1:length(v2)   # invariant: sum of terms before (k, l)
        val += _m_f(v1[begin+k-1]) * _m_f(matrix[begin+k-1, begin+l-1]) * _m_f(v2[begin+l-1])
    end
    return val
end

@inline _m_eltype(a, b) = typeof(_m_f(zero(eltype(a))) * _m_f(zero(eltype(b))))

"""
    mxvg(m1, v2)

Return the product `m1 * v2` of a matrix and a vector of any sizes that fit, as
`mxvg_c` of CSPICE.

Throws `SpiceError(:badarg)` if `size(m1, 2) != length(v2)`.

# Examples
```jldoctest
julia> mxvg([1 2 3; 4 5 6], [1, 1, 1])
2-element Vector{Float64}:
  6.0
 15.0
```

See also [`mxv`](@ref), [`mtxvg`](@ref), [`mxmg`](@ref).
"""
function mxvg(m1::AbstractMatrix, v2::AbstractVector)
    size(m1, 2) == length(v2) || throw(SpiceError(:badarg, "the matrix and the vector do not fit"))
    out = zeros(_m_eltype(m1, v2), size(m1, 1))
    @inbounds for i in 1:size(m1, 1), k in 1:size(m1, 2)   # invariant: out[i] sums terms before k
        out[i] += _m_f(m1[begin+i-1, begin+k-1]) * _m_f(v2[begin+k-1])
    end
    return out
end
"Sum of `g(1)`, ..., `g(M)` added to zero in that order: the static forms give the bits of the array forms."
@inline function _m_ordsum(g::G, ::Val{M}) where {G,M}
    t = ntuple(g, Val(M))
    return foldl(+, t; init=zero(t[1]))
end
@inline _m_ordsum(g::G, ::Val{0}) where {G} = 0.0

mxvg(m1::StaticMatrix{N,M}, v2::StaticVector{M}) where {N,M} =
    SVector{N}(ntuple(i -> _m_ordsum(k -> _m_f(m1[i, k]) * _m_f(v2[k]), Val(M)), Val(N)))

"""
    mtxvg(m1, v2)

Return the product of the transpose of `m1` and the vector `v2`, for any sizes that
fit, as `mtxvg_c` of CSPICE.

Throws `SpiceError(:badarg)` if `size(m1, 1) != length(v2)`.

See also [`mtxv`](@ref), [`mxvg`](@ref).
"""
function mtxvg(m1::AbstractMatrix, v2::AbstractVector)
    size(m1, 1) == length(v2) || throw(SpiceError(:badarg, "the matrix and the vector do not fit"))
    out = zeros(_m_eltype(m1, v2), size(m1, 2))
    @inbounds for i in 1:size(m1, 2), k in 1:size(m1, 1)   # invariant: out[i] sums terms before k
        out[i] += _m_f(m1[begin+k-1, begin+i-1]) * _m_f(v2[begin+k-1])
    end
    return out
end
mtxvg(m1::StaticMatrix{N,M}, v2::StaticVector{N}) where {N,M} =
    SVector{M}(ntuple(i -> _m_ordsum(k -> _m_f(m1[k, i]) * _m_f(v2[k]), Val(N)), Val(M)))

"""
    mxmg(m1, m2)

Return the product `m1 * m2` of two matrices of any sizes that fit, as `mxmg_c` of
CSPICE.

Throws `SpiceError(:badarg)` if `size(m1, 2) != size(m2, 1)`.

See also [`mxm`](@ref), [`mtxmg`](@ref), [`mxmtg`](@ref).
"""
function mxmg(m1::AbstractMatrix, m2::AbstractMatrix)
    size(m1, 2) == size(m2, 1) || throw(SpiceError(:badarg, "the matrix sizes do not fit"))
    out = zeros(_m_eltype(m1, m2), size(m1, 1), size(m2, 2))
    @inbounds for i in 1:size(m1, 1), j in 1:size(m2, 2), k in 1:size(m1, 2)
        # invariant: out[i, j] is the sum of the terms before k
        out[i, j] += _m_f(m1[begin+i-1, begin+k-1]) * _m_f(m2[begin+k-1, begin+j-1])
    end
    return out
end
mxmg(m1::StaticMatrix{N,M}, m2::StaticMatrix{M,P}) where {N,M,P} =
    SMatrix{N,P}(ntuple(q -> _m_ordsum(k -> _m_f(m1[(q-1)%N+1, k]) * _m_f(m2[k, (q-1)÷N+1]), Val(M)), Val(N * P)))

"""
    mtxmg(m1, m2)

Return the product of the transpose of `m1` and `m2`, for any sizes that fit, as
`mtxmg_c` of CSPICE.

Throws `SpiceError(:badarg)` if `size(m1, 1) != size(m2, 1)`.

See also [`mtxm`](@ref), [`mxmg`](@ref).
"""
function mtxmg(m1::AbstractMatrix, m2::AbstractMatrix)
    size(m1, 1) == size(m2, 1) || throw(SpiceError(:badarg, "the matrix sizes do not fit"))
    out = zeros(_m_eltype(m1, m2), size(m1, 2), size(m2, 2))
    @inbounds for i in 1:size(m1, 2), j in 1:size(m2, 2), k in 1:size(m1, 1)
        # invariant: out[i, j] is the sum of the terms before k
        out[i, j] += _m_f(m1[begin+k-1, begin+i-1]) * _m_f(m2[begin+k-1, begin+j-1])
    end
    return out
end
mtxmg(m1::StaticMatrix{N,M}, m2::StaticMatrix{N,P}) where {N,M,P} =
    SMatrix{M,P}(ntuple(q -> _m_ordsum(k -> _m_f(m1[k, (q-1)%M+1]) * _m_f(m2[k, (q-1)÷M+1]), Val(N)), Val(M * P)))

"""
    mxmtg(m1, m2)

Return the product of `m1` and the transpose of `m2`, for any sizes that fit, as
`mxmtg_c` of CSPICE.

Throws `SpiceError(:badarg)` if `size(m1, 2) != size(m2, 2)`.

See also [`mxmt`](@ref), [`mxmg`](@ref).
"""
function mxmtg(m1::AbstractMatrix, m2::AbstractMatrix)
    size(m1, 2) == size(m2, 2) || throw(SpiceError(:badarg, "the matrix sizes do not fit"))
    out = zeros(_m_eltype(m1, m2), size(m1, 1), size(m2, 1))
    @inbounds for i in 1:size(m1, 1), j in 1:size(m2, 1), k in 1:size(m1, 2)
        # invariant: out[i, j] is the sum of the terms before k
        out[i, j] += _m_f(m1[begin+i-1, begin+k-1]) * _m_f(m2[begin+j-1, begin+k-1])
    end
    return out
end
mxmtg(m1::StaticMatrix{N,M}, m2::StaticMatrix{P,M}) where {N,M,P} =
    SMatrix{N,P}(ntuple(q -> _m_ordsum(k -> _m_f(m1[(q-1)%N+1, k]) * _m_f(m2[(q-1)÷N+1, k]), Val(M)), Val(N * P)))

"""
    xposeg(matrix)

Return the transpose of a matrix of any size, as `xposeg_c` of CSPICE.

See also [`xpose`](@ref), [`xpose6`](@ref).
"""
xposeg(matrix::AbstractMatrix) = permutedims(_m_f.(matrix))
xposeg(matrix::StaticMatrix) = SMatrix(transpose(_m_f.(matrix)))
