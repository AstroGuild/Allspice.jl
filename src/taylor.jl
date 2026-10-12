# Truncated Taylor series arithmetic. A `Taylor{N,T}` is the jet of a function of one
# variable (time) to order N; passing one through any evaluator that is generic over
# the epoch type gives all time derivatives to order N in one pass.

"""
    Taylor{N,T} <: Real

Immutable truncated Taylor series of order `N`; `x.c::NTuple{N+1,T}` holds the
coefficients, `x.c[k+1] = f⁽ᵏ⁾/k!` for `k = 0:N`.

Contract: every operation returns the series of the result truncated at order `N`; a
value is isbits when `T` is. Comparisons are not defined: code that branches must use
`primal`.

Julia cannot write the field type `NTuple{N+1,T}` with only the parameters `N` and
`T`, so the storage is the value `v` (coefficient 0) and the tuple `d` of the `N`
higher coefficients; the property `c` joins them and `Taylor{N,T}(c)` splits them.
"""
struct Taylor{N,T<:Real} <: Real
    v::T
    d::NTuple{N,T}
    global _taylor(v::T, d::NTuple{N,T}) where {N,T<:Real} = new{N,T}(v, d)
end

"`Taylor{N,T}(c)`: series with the `N+1` coefficients `c` (requires `length(c) == N+1`)."
@inline function Taylor{N,T}(c::Tuple) where {N,T<:Real}
    length(c) == N + 1 ||
        throw(ArgumentError("Taylor{$N} needs $(N + 1) coefficients, got $(length(c))"))
    return _taylor(convert(T, c[1]), convert(NTuple{N,T}, Base.tail(c)))
end

"`Taylor{N}(c)`: series with coefficients `c` of one element type (requires `length(c) == N+1`)."
@inline Taylor{N}(c::Tuple{T,Vararg{T}}) where {N,T<:Real} = Taylor{N,T}(c)

"`Taylor(c)`: series of order `length(c) - 1` with coefficients `c`."
@inline Taylor(c::Tuple{T,Vararg{T,N}}) where {N,T<:Real} = _taylor(c[1], Base.tail(c))

"`Taylor{N,T}(x)`: the constant `x` as a series (all higher coefficients zero)."
@inline Taylor{N,T}(x::Real) where {N,T<:Real} =
    _taylor(convert(T, x), ntuple(_ -> zero(T), Val(N)))
@inline Taylor{N,T}(x::Taylor{N,T}) where {N,T<:Real} = x
@inline Taylor{N,T}(x::Taylor{N}) where {N,T<:Real} =
    _taylor(convert(T, x.v), convert(NTuple{N,T}, x.d))

"`x.c` is the full coefficient tuple `NTuple{N+1,T}`; `x.v` and `x.d` are the storage."
@inline Base.getproperty(x::Taylor, s::Symbol) =
    s === :c ? (getfield(x, :v), getfield(x, :d)...) : getfield(x, s)
Base.propertynames(::Taylor) = (:c, :v, :d)

"""
    taylor_var(t, Val(N)) -> Taylor{N,typeof(t)}

The independent variable at the value `t`: coefficients `(t, 1, 0, ..., 0)`.
"""
@inline taylor_var(t::T, ::Val{N}) where {T<:Real,N} =
    _taylor(t, ntuple(k -> k == 1 ? one(T) : zero(T), Val(N)))

"""
    derivatives(x::Taylor{N,T}) -> NTuple{N+1,T}

Ordinary derivatives of orders 0 to N: element `k+1` is `x.c[k+1] * k!`.
"""
@generated function derivatives(x::Taylor{N,T}) where {N,T}
    fact = [factorial(big(k)) for k in 1:N]
    all(f -> f <= maxintfloat(Float64), fact) ||
        return :(throw(ArgumentError("order too high for an exact factorial")))
    terms = [:(x.d[$k] * $(Float64(fact[k]))) for k in 1:N]
    return :(tuple(x.v, $(terms...)))
end

"""
    primal(x) -> Float64

Plain value of `x` for discrete decisions: `Float64(x)` for a plain real, the order-0
coefficient (recursively) for a `Taylor`.
"""
@inline primal(x::Real) = Float64(x)
@inline primal(x::Taylor) = primal(x.v)

# --- conversion and promotion -------------------------------------------------------

# Plain number types that mix with `Taylor` values. The list is explicit so that the
# promotion rule cannot become ambiguous with the rule of another `Real` wrapper type.
const TaylorScalar = Union{AbstractFloat,Integer,Rational,AbstractIrrational}

Base.convert(::Type{Taylor{N,T}}, x::Taylor{N,T}) where {N,T<:Real} = x
Base.convert(::Type{Taylor{N,T}}, x::Taylor{N}) where {N,T<:Real} = Taylor{N,T}(x)
Base.convert(::Type{Taylor{N,T}}, x::TaylorScalar) where {N,T<:Real} = Taylor{N,T}(x)

Base.promote_rule(::Type{Taylor{N,T}}, ::Type{S}) where {N,T<:Real,S<:TaylorScalar} =
    Taylor{N,promote_type(T, S)}
Base.promote_rule(::Type{Taylor{N,T}}, ::Type{Taylor{N,S}}) where {N,T<:Real,S<:Real} =
    Taylor{N,promote_type(T, S)}

"`zero` and `one` of a Taylor type are the constants 0 and 1."
Base.zero(::Type{Taylor{N,T}}) where {N,T<:Real} = Taylor{N,T}(zero(T))
Base.one(::Type{Taylor{N,T}}) where {N,T<:Real} = Taylor{N,T}(one(T))
Base.zero(x::Taylor) = zero(typeof(x))
Base.one(x::Taylor) = one(typeof(x))

Base.show(io::IO, x::Taylor{N,T}) where {N,T} = print(io, "Taylor{", N, ",", T, "}", x.c)

# --- linear operations --------------------------------------------------------------

"Sum, difference and negation act coefficient by coefficient."
@inline Base.:+(x::Taylor{N,T}, y::Taylor{N,T}) where {N,T<:Real} =
    _taylor(x.v + y.v, map(+, x.d, y.d))
@inline Base.:-(x::Taylor{N,T}, y::Taylor{N,T}) where {N,T<:Real} =
    _taylor(x.v - y.v, map(-, x.d, y.d))
@inline Base.:-(x::Taylor) = _taylor(-x.v, map(-, x.d))
@inline Base.:+(x::Taylor) = x

# A plain scalar of the coefficient type (or a machine integer) changes only what it
# must; other plain types go through promotion.
const _FastScalar{T} = Union{T,Int}

@inline Base.:+(x::Taylor{N,T}, a::_FastScalar{T}) where {N,T<:AbstractFloat} = _taylor(x.v + a, x.d)
@inline Base.:+(a::_FastScalar{T}, x::Taylor{N,T}) where {N,T<:AbstractFloat} = _taylor(a + x.v, x.d)
@inline Base.:-(x::Taylor{N,T}, a::_FastScalar{T}) where {N,T<:AbstractFloat} = _taylor(x.v - a, x.d)
@inline Base.:-(a::_FastScalar{T}, x::Taylor{N,T}) where {N,T<:AbstractFloat} =
    _taylor(a - x.v, map(-, x.d))
@inline Base.:*(x::Taylor{N,T}, a::_FastScalar{T}) where {N,T<:AbstractFloat} =
    _taylor(x.v * a, map(c -> c * a, x.d))
@inline Base.:*(a::_FastScalar{T}, x::Taylor{N,T}) where {N,T<:AbstractFloat} =
    _taylor(a * x.v, map(c -> a * c, x.d))
@inline Base.:/(x::Taylor{N,T}, a::_FastScalar{T}) where {N,T<:AbstractFloat} =
    _taylor(x.v / a, map(c -> c / a, x.d))

# --- products, quotients and elementary functions -----------------------------------
#
# The bodies below are generated as straight-line code on the coefficients, so that
# they stay allocation-free and type-stable for every N. In the generated code `a0`,
# `a1`, ... are the coefficients of the first argument, `b*` of the second.

_sym(p, k) = Symbol(p, k)

"Expressions that bind `p0 ... pN` to the coefficients of the argument named `arg`."
function _unpack(p, arg, N)
    ex = Any[:($(_sym(p, 0)) = $arg.v)]
    for k in 1:N
        push!(ex, :($(_sym(p, k)) = $arg.d[$k]))
    end
    return ex
end

"Expression that builds the Taylor value with coefficients `p0 ... pN`."
_pack(p, N) = :(_taylor($(_sym(p, 0)), tuple($([_sym(p, k) for k in 1:N]...))))

"Left-to-right sum of the expressions `terms` (at least one)."
_sum(terms) = length(terms) == 1 ? terms[1] : Expr(:call, :+, terms...)

"""
    x * y

Cauchy product truncated at order N: `r_k = Σ_{j=0}^{k} x_j y_{k-j}`.
"""
@generated function Base.:*(x::Taylor{N,T}, y::Taylor{N,T}) where {N,T<:Real}
    ex = vcat(_unpack(:a, :x, N), _unpack(:b, :y, N))
    for k in 0:N
        push!(ex, :($(_sym(:r, k)) = $(_sum([:($(_sym(:a, j)) * $(_sym(:b, k - j))) for j in 0:k]))))
    end
    return Expr(:block, Expr(:meta, :inline), ex..., _pack(:r, N))
end

"""
    x / y

Series quotient: `q_k = (x_k - Σ_{j=1}^{k} y_j q_{k-j}) / y_0`. Requires `y_0 != 0`
(not checked; a zero gives non-finite coefficients, as for plain reals).
"""
@generated function Base.:/(x::Taylor{N,T}, y::Taylor{N,T}) where {N,T<:Real}
    ex = vcat(_unpack(:a, :x, N), _unpack(:b, :y, N))
    push!(ex, :(q0 = a0 / b0))
    for k in 1:N
        s = _sum([:($(_sym(:b, j)) * $(_sym(:q, k - j))) for j in 1:k])
        push!(ex, :($(_sym(:q, k)) = ($(_sym(:a, k)) - $s) / b0))
    end
    return Expr(:block, Expr(:meta, :inline), ex..., _pack(:q, N))
end

"`inv(x)` is `one(x) / x`."
@inline Base.inv(x::Taylor) = one(x) / x

"""
    sincos(x) -> (sin(x), cos(x))

Coupled recurrences `s_k = (1/k) Σ_{j=1}^{k} j x_j c_{k-j}` and
`c_k = -(1/k) Σ_{j=1}^{k} j x_j s_{k-j}`, started from `sincos(x_0)`.
"""
@generated function Base.sincos(x::Taylor{N,T}) where {N,T<:Real}
    ex = _unpack(:a, :x, N)
    push!(ex, :((s0, c0) = sincos(a0)))
    for j in 2:N
        push!(ex, :($(_sym(:ja, j)) = $j * $(_sym(:a, j))))
    end
    ja(j) = j == 1 ? :a1 : _sym(:ja, j)
    for k in 1:N
        ss = _sum([:($(ja(j)) * $(_sym(:c, k - j))) for j in 1:k])
        cs = _sum([:($(ja(j)) * $(_sym(:s, k - j))) for j in 1:k])
        if k == 1
            push!(ex, :(s1 = $ss), :(c1 = -$cs))
        else
            push!(ex, :($(_sym(:s, k)) = $ss / $k), :($(_sym(:c, k)) = -$cs / $k))
        end
    end
    return Expr(:block, Expr(:meta, :inline), ex..., :(($(_pack(:s, N)), $(_pack(:c, N)))))
end

"`sin(x)` and `cos(x)` are the two parts of `sincos(x)`."
@inline Base.sin(x::Taylor) = sincos(x)[1]
@inline Base.cos(x::Taylor) = sincos(x)[2]

"""
    sqrt(x)

Series square root: `r_0 = sqrt(x_0)`, `r_k = (x_k - Σ_{j=1}^{k-1} r_j r_{k-j}) / (2 r_0)`.
Requires `x_0 > 0` for finite higher coefficients (not checked).
"""
@generated function Base.sqrt(x::Taylor{N,T}) where {N,T<:Real}
    ex = _unpack(:a, :x, N)
    push!(ex, :(r0 = sqrt(a0)), :(tr0 = 2 * r0))
    for k in 1:N
        if k == 1
            push!(ex, :(r1 = a1 / tr0))
        else
            s = _sum([:($(_sym(:r, j)) * $(_sym(:r, k - j))) for j in 1:(k - 1)])
            push!(ex, :($(_sym(:r, k)) = ($(_sym(:a, k)) - $s) / tr0))
        end
    end
    return Expr(:block, Expr(:meta, :inline), ex..., _pack(:r, N))
end
