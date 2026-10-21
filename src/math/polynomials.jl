# Polynomial calls of CSPICE: Chebyshev expansions, Lagrange and Hermite interpolation,
# the derivatives of a polynomial, a numerical derivative, the roots of a quadratic and
# the diagonalization of a symmetric 2 by 2 matrix.
#
# The recurrences are those of the SPK and CK evaluators: the Clenshaw steps of
# `spk_eval.jl` (`_clenshaw_step`, `_clenshaw_end`), the integral of `_chbigr` in
# `spk_types/chebyshev.jl`, and the kernels `_neville!`, `_neville_d!` and
# `_hermite_table!` of `spk_types/interp.jl`, and the roots and the diagonalization
# `_m_rquad` and `_m_diags2` of `math/ellipsoid.jl`. Here the data come from the vectors
# of the caller, and a work vector is `similar` to the input, so that a static input
# allocates nothing.

# --- shared input handling --------------------------------------------------------------

"Element type of the result of a call on the numbers `x` and the elements of `vs`."
@inline _m_restype(x, vs::AbstractVector...) =
    typeof(+(_m_f(x), map(v -> _m_f(zero(eltype(v))), vs)...))

"The midpoint and the radius of the interval of a Chebyshev expansion."
@inline function _m_x2s(x2s::AbstractVector)
    length(x2s) == 2 || _m_baddim("x2s", 2)
    @inbounds return _m_f(x2s[begin]), _m_f(x2s[begin+1])
end

@noinline _m_nodegree() = throw(SpiceError(:invaliddegree,
    "the expansion has no coefficient: the degree must be non-negative"))

@noinline _m_badsize(n) = throw(SpiceError(:invalidsize, "array size must be positive; was $n"))

@noinline _m_zerostep() = throw(SpiceError(:invalidstepsize, "the step size is zero"))

@noinline _m_sameabscissa(xvals, i, j) = throw(SpiceError(:dividebyzero,
    "xvals[$i] = xvals[$j] = $(primal(_m_f(xvals[begin+i-1])))"))

"`onequal` of an interpolation on the grid 1, 2, 3, ..., where no two abscissas are equal."
@inline _m_never(_, _) = nothing

# --- Chebyshev expansions ---------------------------------------------------------------

"""
    _m_cheb(cp, s::T, ::Val{N}) -> NTuple{N+1,T}

Value and derivatives with respect to `s`, orders 0 to N, of the Chebyshev series with
the coefficients `cp` at `s`, by the Clenshaw steps of `_cheb3`. Orders 0 and 1 take
the operations of `chbval` and `chbint`. Requires `length(cp) >= 1`.
"""
@inline function _m_cheb(cp::AbstractVector, s::T, ::Val{N}) where {T,N}
    M = Val(N + 1)
    z = ntuple(_ -> zero(T), M)
    s2 = s * 2
    b1, b2 = z, z
    # invariant: b1 and b2 are the Clenshaw terms b_{k+1} and b_{k+2}, with derivatives
    @inbounds for k in length(cp)-1:-1:1
        b1, b2 = _clenshaw_step(_m_f(cp[begin+k]), s2, b1, b2), b1
    end
    return _clenshaw_end(_m_f(@inbounds cp[begin]), s, b1, b2)
end

"Scaled abscissa `s` and radius of `x` for the expansion `cp` on `x2s`, in the result type."
@inline function _m_cheb_s(cp::AbstractVector, x2s::AbstractVector, x::Real)
    mid, rad = _m_x2s(x2s)
    isempty(cp) && _m_nodegree()
    T = _m_restype(x, cp, x2s)
    return convert(T, (_m_f(x) - mid) / rad), rad
end

"""
    chbval(cp, x2s, x) -> p

Return the value at `x` of the Chebyshev expansion with the coefficients `cp`, as
`chbval_c` of CSPICE.

The expansion is `p(x) = sum(cp[k] * T_{k-1}(s))` with `s = (x - x2s[1]) / x2s[2]`:
`x2s` holds the midpoint and the radius of the interval of the expansion. The degree
is `length(cp) - 1`. A value of `x` outside the interval is not an error. As in CSPICE,
`x2s[2]` is not checked: 0 gives NaN, and a negative radius is used as given.

Throws `SpiceError(:invaliddegree)` if `cp` is empty. Throws `SpiceError(:badarg)` if
`x2s` does not have 2 elements.

# Examples
`T_0 + 2 T_1 + 3 T_2` at `s = 0.5`, where `T_2(s) = 2s² - 1`:
```jldoctest
julia> chbval([1, 2, 3], [0, 1], 0.5)
0.5
```

See also [`chbint`](@ref), [`chbder`](@ref), [`chbigr`](@ref).
"""
function chbval(cp::AbstractVector, x2s::AbstractVector, x::Real)
    s, _ = _m_cheb_s(cp, x2s, x)
    return _m_cheb(cp, s, Val(0))[1]
end

"""
    chbint(cp, x2s, x) -> (p, dpdx)

Return the value and the first derivative with respect to `x` of the Chebyshev
expansion of [`chbval`](@ref), as `chbint_c` of CSPICE.

Throws as [`chbval`](@ref).

# Examples
```jldoctest
julia> chbint([1, 2, 3], [0, 1], 0.5)
(0.5, 8.0)
```

See also [`chbder`](@ref).
"""
function chbint(cp::AbstractVector, x2s::AbstractVector, x::Real)
    s, rad = _m_cheb_s(cp, x2s, x)
    p, dp = _m_cheb(cp, s, Val(1))
    return p, dp / rad
end

"""
    chbder(cp, x2s, x, nderiv) -> Vector

Return the value and the derivatives of orders 1 to `nderiv` with respect to `x` of
the Chebyshev expansion of [`chbval`](@ref), as `chbder_c` of CSPICE.

Element `k + 1` of the result is the derivative of order `k`. A negative `nderiv`
gives an empty result, as for [`polyds`](@ref); `chbder_c` signals no error then.

Throws `SpiceError(:invaliddegree)` if `cp` is empty, and `SpiceError(:badarg)` if `x2s`
does not have 2 elements.

# Examples
```jldoctest
julia> chbder([1, 2, 3], [0, 1], 0.5, 3)
4-element Vector{Float64}:
  0.5
  8.0
 12.0
  0.0
```

See also [`chbint`](@ref), [`polyds`](@ref).
"""
function chbder(cp::AbstractVector, x2s::AbstractVector, x::Real, nderiv::Integer)
    s, rad = _m_cheb_s(cp, x2s, x)
    T = typeof(s)
    nd = Int(nderiv)
    nd < 0 && return T[]
    # The recurrence of `_clenshaw_step` for an order that is known only at run time:
    # column m+1 of `b` holds the Clenshaw terms b_k, b_{k+1}, b_{k+2} of derivative m.
    b = zeros(T, 3, nd + 1)
    s2 = s * 2.0
    # invariant: b[1, m+1] and b[2, m+1] are the terms b_{k+1} and b_{k+2} of order m
    @inbounds for k in length(cp):-1:2
        b[3, 1] = b[2, 1]
        b[2, 1] = b[1, 1]
        b[1, 1] = _m_f(cp[begin+k-1]) + (s2 * b[2, 1] - b[3, 1])
        scale = 2.0
        for m in 1:nd
            b[3, m+1] = b[2, m+1]
            b[2, m+1] = b[1, m+1]
            b[1, m+1] = b[2, m] * scale + b[2, m+1] * s2 - b[3, m+1]
            scale += 2.0
        end
    end
    dpdxs = Vector{T}(undef, nd + 1)
    @inbounds dpdxs[1] = _m_f(cp[begin]) + (s * b[1, 1] - b[2, 1])
    scale = 1.0
    @inbounds for m in 1:nd
        dpdxs[m+1] = b[1, m] * scale + b[1, m+1] * s - b[2, m+1]
        scale += 1.0
    end
    # derivative m is divided by rad^m, the power formed by products as in CSPICE
    radm = rad
    @inbounds for m in 1:nd
        dpdxs[m+1] /= radm
        radm = rad * radm
    end
    return dpdxs
end

"""
    chbigr(cp, x2s, x) -> (p, itgrlp)

Return the value at `x` of the Chebyshev expansion of [`chbval`](@ref) and its integral
from the midpoint `x2s[1]` to `x`, as `chbigr_c` of CSPICE.

The integral is in the units of the expansion times the units of `x`.

Throws `SpiceError(:invaliddegree)` if `cp` is empty, `SpiceError(:invalidradius)` if
`x2s[2] <= 0`, and `SpiceError(:badarg)` if `x2s` does not have 2 elements.

# Examples
The integral of `1 + 2s` from 0 to 0.5 on the interval `[-1, 1]`:
```jldoctest
julia> chbigr([1, 2], [0, 1], 0.5)
(2.0, 0.75)
```

See also [`chbval`](@ref).
"""
function chbigr(cp::AbstractVector, x2s::AbstractVector, x::Real)
    mid, rad = _m_x2s(x2s)
    isempty(cp) && _m_nodegree()
    primal(rad) > 0.0 || throw(SpiceError(:invalidradius,
        "the interval radius must be positive; was $(primal(rad))"))
    Base.require_one_based_indexing(cp)
    T = _m_restype(x, cp, x2s)
    return _chbigr(cp, 1, length(cp), 1.0, mid, rad, convert(T, _m_f(x)))
end

# --- Lagrange and Hermite interpolation -------------------------------------------------

"Number of points of an interpolation on `xvals`, `ny` values per point: the checks of CSPICE, then the length of `yvals`."
@inline function _m_npoints(xvals::AbstractVector, yvals::AbstractVector, ny::Int)
    n = length(xvals)
    n >= 1 || _m_badsize(n)
    length(yvals) == ny * n || _m_baddim("yvals", ny * n)
    Base.require_one_based_indexing(xvals, yvals)
    return n
end

"Work vector of the type `T` with the elements of `yvals`."
@inline function _m_work(yvals::AbstractVector, ::Type{T}) where {T}
    w = similar(yvals, T)
    @inbounds for i in eachindex(yvals)
        w[i] = _m_f(yvals[i])
    end
    return w
end

"""
    lgrint(xvals, yvals, x) -> p

Return the value at `x` of the polynomial of degree `n - 1` through the `n` points
`(xvals[i], yvals[i])`, as `lgrint_c` of CSPICE (Lagrange interpolation).

The abscissas need not be ordered. A static `yvals` gives a call without allocation.

Throws `SpiceError(:invalidsize)` if `xvals` is empty, `SpiceError(:dividebyzero)` if
two abscissas are equal, and `SpiceError(:badarg)` if `yvals` does not have the length
of `xvals`.

# Examples
The line through `(1, 3)` and `(2, 5)`:
```jldoctest
julia> lgrint([1, 2], [3, 5], 4)
9.0
```

See also [`lgrind`](@ref), [`lgresp`](@ref), [`hrmint`](@ref).
"""
function lgrint(xvals::AbstractVector, yvals::AbstractVector, x::Real)
    n = _m_npoints(xvals, yvals, 1)
    w = _m_work(yvals, _m_restype(x, xvals, yvals))
    return _neville!(w, n, i -> _m_f(@inbounds xvals[i]), _m_f(x),
        (i, j) -> _m_sameabscissa(xvals, i, j))
end

"""
    lgrind(xvals, yvals, x) -> (p, dp)

Return the value and the derivative at `x` of the interpolating polynomial of
[`lgrint`](@ref), as `lgrind_c` of CSPICE.

Throws as [`lgrint`](@ref).

# Examples
```jldoctest
julia> lgrind([0, 1, 2], [0, 1, 4], 3)    # x²
(9.0, 6.0)
```

See also [`hrmint`](@ref).
"""
function lgrind(xvals::AbstractVector, yvals::AbstractVector, x::Real)
    n = _m_npoints(xvals, yvals, 1)
    T = _m_restype(x, xvals, yvals)
    w = _m_work(yvals, T)
    return _neville_d!(w, similar(w), n, i -> _m_f(@inbounds xvals[i]), _m_f(x),
        (i, j) -> _m_sameabscissa(xvals, i, j))
end

"""
    lgresp(first, step, yvals, x) -> p

Return the value at `x` of the polynomial through the points
`(first + (i - 1) * step, yvals[i])`, as `lgresp_c` of CSPICE (Lagrange interpolation on
equally spaced abscissas).

Throws `SpiceError(:invalidsize)` if `yvals` is empty and
`SpiceError(:invalidstepsize)` if `step == 0`.

# Examples
```jldoctest
julia> lgresp(10, 5, [0, 25, 100], 25)     # (x - 10)²
225.0
```

See also [`lgrint`](@ref), [`hrmesp`](@ref).
"""
function lgresp(first::Real, step::Real, yvals::AbstractVector, x::Real)
    n = length(yvals)
    n >= 1 || _m_badsize(n)
    h = _m_f(step)
    primal(h) == 0.0 && _m_zerostep()
    Base.require_one_based_indexing(yvals)
    newx = (_m_f(x) - _m_f(first)) / h + 1.0
    w = _m_work(yvals, typeof(newx + _m_f(zero(eltype(yvals)))))
    return _neville!(w, n, Float64, newx, _m_never)
end

"""
    hrmint(xvals, yvals, x) -> (f, df)

Return the value and the derivative at `x` of the Hermite polynomial that has, at
`xvals[i]`, the value `yvals[2i-1]` and the derivative `yvals[2i]`, as `hrmint_c` of
CSPICE.

The degree is `2n - 1` for `n` abscissas. The abscissas need not be ordered. A static
`yvals` gives a call without allocation.

Throws `SpiceError(:invalidsize)` if `xvals` is empty, `SpiceError(:dividebyzero)` if
two abscissas are equal, and `SpiceError(:badarg)` if `yvals` does not have twice the
length of `xvals`.

# Examples
`x³`, with the values and the derivatives at 1 and 2:
```jldoctest
julia> hrmint([1, 2], [1, 3, 8, 12], 3)
(27.0, 27.0)
```

See also [`hrmesp`](@ref), [`lgrind`](@ref).
"""
function hrmint(xvals::AbstractVector, yvals::AbstractVector, x::Real)
    n = _m_npoints(xvals, yvals, 2)
    w1 = _m_work(yvals, _m_restype(x, xvals, yvals))
    return _hermite_table!(w1, similar(w1), n, i -> _m_f(@inbounds xvals[i]), _m_f(x),
        (i, j) -> _m_sameabscissa(xvals, i, j))
end

"""
    hrmesp(first, step, yvals, x) -> (f, df)

Return the value and the derivative at `x` of the Hermite polynomial that has, at
`first + (i - 1) * step`, the value `yvals[2i-1]` and the derivative `yvals[2i]`, as
`hrmesp_c` of CSPICE (equally spaced abscissas).

Throws `SpiceError(:invalidsize)` if `yvals` is empty, `SpiceError(:invalidstepsize)`
if `step == 0`, and `SpiceError(:badarg)` if the length of `yvals` is odd.

# Examples
```jldoctest
julia> hrmesp(1, 1, [1, 3, 8, 12], 3)       # x³, as for hrmint
(27.0, 27.0)
```

See also [`hrmint`](@ref), [`lgresp`](@ref).
"""
function hrmesp(first::Real, step::Real, yvals::AbstractVector, x::Real)
    isempty(yvals) && _m_badsize(0)
    iseven(length(yvals)) || throw(SpiceError(:badarg,
        "yvals must hold a value and a derivative for each abscissa: its length is odd"))
    n = length(yvals) ÷ 2
    h = _m_f(step)
    primal(h) == 0.0 && _m_zerostep()
    Base.require_one_based_indexing(yvals)
    newx = (_m_f(x) - _m_f(first)) / h + 1.0
    w1 = _m_work(yvals, typeof(newx + _m_f(zero(eltype(yvals))) * h))
    @inbounds for i in 2:2:2n
        w1[i] = w1[i] * h
    end
    f, df = _hermite_table!(w1, similar(w1), n, Float64, newx, _m_never)
    return f, df / h
end

# --- derivatives of a polynomial --------------------------------------------------------

"""
    polyds(coeffs, nderiv, t) -> Vector

Return the value and the derivatives of orders 1 to `nderiv` at `t` of the polynomial
`coeffs[1] + coeffs[2] t + coeffs[3] t² + ...`, as `polyds_c` of CSPICE.

Element `k + 1` of the result is the derivative of order `k`. An empty `coeffs` is the
zero polynomial. A negative `nderiv` gives an empty result, as CSPICE computes nothing
then.

# Examples
```jldoctest
julia> polyds([1, 2, 3], 3, 2)               # 1 + 2t + 3t²
4-element Vector{Float64}:
 17.0
 14.0
  6.0
  0.0
```

See also [`chbder`](@ref).
"""
function polyds(coeffs::AbstractVector, nderiv::Integer, t::Real)
    tt = _m_f(t)
    T = _m_restype(tt, coeffs)
    nd = Int(nderiv)
    p = zeros(T, max(nd + 1, 0))
    nd < 0 && return p
    # invariant: p[i+1] is derivative i of the polynomial of the coefficients k:end
    @inbounds for k in length(coeffs):-1:1
        scale = Float64(nd)
        for i in nd:-1:1
            p[i+1] = tt * p[i+1] + scale * p[i]
            scale -= 1.0
        end
        p[1] = tt * p[1] + _m_f(coeffs[begin+k-1])
    end
    return p
end

"""
    qderiv(f0, f2, delta) -> dfdt

Return the derivative of a vector function at a point from its values `f0` at
`t - delta` and `f2` at `t + delta`, by the central difference `(f2 - f0) / (2 delta)`,
as `qderiv_c` of CSPICE.

The result has the length of the arguments; static arguments give a static result.

Throws `SpiceError(:dividebyzero)` if `delta == 0`, and `SpiceError(:badarg)` if the
lengths of `f0` and `f2` differ.

# Examples
```jldoctest
julia> qderiv([1, 4], [9, 16], 1) == [4, 6]
true
```
"""
function qderiv(f0::AbstractVector, f2::AbstractVector, delta::Real)
    d = _m_f(delta)
    primal(d) == 0.0 && throw(SpiceError(:dividebyzero,
        "delta abscissa value is zero; a non-zero value is required"))
    return vlcomg(0.5 / d, f2, -0.5 / d, f0)
end

# --- quadratic equations, symmetric 2 by 2 matrices -------------------------------------

"""
    rquad(a, b, c) -> (root1, root2)

Return the roots of the quadratic equation `a x² + b x + c = 0`, as `rquad_c` of
CSPICE. Each root is an `SVector` of its real part and its imaginary part.

`root1` is `(-b + sqrt(b² - 4ac)) / 2a` and `root2` the root with the minus sign. If
`a == 0`, both are the root `-c / b` of the linear equation. The coefficients are
divided by the largest of their magnitudes first, unless that makes one of them zero.

Throws `SpiceError(:degeneratecase)` if `a == b == 0`.

# Examples
```jldoctest
julia> rquad(1, 0, -4)
([2.0, 0.0], [-2.0, 0.0])

julia> rquad(1, 0, 4)
([-0.0, 2.0], [-0.0, -2.0])
```

See also [`diags2`](@ref).
"""
function rquad(a::Real, b::Real, c::Real)
    root1, root2 = _m_rquad(promote(_m_f(a), _m_f(b), _m_f(c))...)
    return SVector(root1), SVector(root2)
end

"""
    diags2(symmat) -> (d, r)

Return the diagonal matrix `d` and the rotation `r` that diagonalize the symmetric 2 by
2 matrix `symmat`, so that `d = r' * symmat * r`, as `diags2_c` of CSPICE.

Only the upper triangle of `symmat` is read. The columns of `r` are unit eigenvectors,
and `d[i, i]` is the eigenvalue of column `i`. Of the rotations that diagonalize
`symmat`, `r` has an angle of the smallest magnitude, so a diagonal `symmat` gives the
identity. If `symmat[1, 1] == symmat[2, 2]`, the angle is ±pi/4 and `d[2, 2]` is the
larger eigenvalue.

Throws `SpiceError(:badarg)` if `symmat` is not 2 by 2.

# Examples
```jldoctest
julia> d, r = diags2([2 1; 1 2]);

julia> round.(d; digits=12) == [1 0; 0 3]
true
```

See also [`rquad`](@ref).
"""
function diags2(symmat::AbstractMatrix)
    size(symmat) == (2, 2) || _m_baddim("a symmetric matrix", "2 by 2")
    @inbounds s11, s12, s22 = promote(_m_f(symmat[begin, begin]), _m_f(symmat[begin, begin+1]),
        _m_f(symmat[begin+1, begin+1]))
    d1, d2, r = _m_diags2(s11, s12, s22)
    z = zero(d1)
    # CSPICE copies the off-diagonal zero of a diagonal input with its sign
    diag = SMatrix{2,2}(d1, z, primal(s12) == 0.0 ? s12 : z, d2)
    return diag, SMatrix{2,2}(r[1], r[3], r[2], r[4])
end
