# Evaluation of the SPK data types other than the fixed-interval Chebyshev types 2 and
# 3: types 1, 5, 8, 9, 12, 13, 14, 15, 17, 18, 19, 20 and 21. Type 10 (two-line
# elements) needs the kernel set and is in `spk_types/type10.jl`.
#
# Each type has one function `state(data, seg, et::T) -> (pos::Vec3{T}, vel::Vec3{T})`
# that follows the algorithm of the CSPICE reader and evaluator pair `spkrNN`/`spkeNN`
# and is generic over the epoch number type. With a plain epoch the result is the
# CSPICE state. Higher time derivatives come from the same function with a truncated
# Taylor series as the epoch (`_state_jet`), so that order k >= 2 of `spk_eval` is
# derivative k-1 of the velocity function that CSPICE defines, for every type.
#
# A type other than 2, 3 and 20 has no directory values in its `Segment` descriptor:
# the evaluator reads them from the end of the segment (`data[seg.eaddr - k]`) at each
# call and checks them, because `open_daf` does not validate these types.

@noinline _badsegment(seg::Segment, msg::String) = throw(SpiceError(:baddaf,
    "SPK type $(seg.dtype) segment at $(seg.baddr):$(seg.eaddr): $msg"))

"""
    _seg_count(seg, x::Float64, lo::Int, hi::Int, what) -> Int

Nearest whole number to `x` (as `i_dnnt` of f2c); throws `SpiceError(:baddaf)` that
names `what` unless `x` is finite and the result is in `lo:hi`.
"""
@inline function _seg_count(seg::Segment, x::Float64, lo::Int, hi::Int, what::String)
    (isfinite(x) && lo - 0.5 < x < hi + 0.5) || _badsegment(seg, _count_msg(what, x, lo, hi))
    n = unsafe_trunc(Int, round(x, RoundNearestTiesAway))
    (lo <= n <= hi) || _badsegment(seg, _count_msg(what, x, lo, hi))
    return n
end
@noinline _count_msg(what::String, x::Float64, lo::Int, hi::Int) = "$what = $x is not in $lo:$hi"

"Largest whole number that a segment count may be (the 32-bit address range of a DAF)."
const _MAXCOUNT = Int(typemax(Int32))

"""
    _lstltd(data, base, n, x) -> Int

Number of elements of the ascending range `data[base+1:base+n]` that are less than
`x` (the index of the last such element, as `lstltd` of CSPICE). Requires `n >= 0`
and the range inside `data` (not checked).
"""
@inline function _lstltd(data::Vector{Float64}, base::Int, n::Int, x::Float64)
    lo, hi = 0, n
    # invariant: data[base+i] < x for 1 <= i <= lo, and data[base+i] >= x for i > hi
    @inbounds while lo < hi
        mid = (lo + hi + 1) >>> 1
        if data[base+mid] < x
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

"""
    _lstled(data, base, n, x) -> Int

As `_lstltd`, for the elements that are less than or equal to `x` (`lstled`).
"""
@inline function _lstled(data::Vector{Float64}, base::Int, n::Int, x::Float64)
    lo, hi = 0, n
    # invariant: data[base+i] <= x for 1 <= i <= lo, and data[base+i] > x for i > hi
    @inbounds while lo < hi
        mid = (lo + hi + 1) >>> 1
        if data[base+mid] <= x
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

"Whole part of `x` as an `Int`, with `x` limited to a range that cannot overflow."
@inline _trunc_index(x::Float64) = unsafe_trunc(Int, clamp(x, -1.0e15, 1.0e15))

"Nearest whole number to `x` as an `Int` (ties away from zero, as `i_dnnt`), limited as `_trunc_index`."
@inline _round_index(x::Float64) =
    unsafe_trunc(Int, round(clamp(x, -1.0e15, 1.0e15), RoundNearestTiesAway))

"Three consecutive doubles `data[a:a+2]` as a vector; requires the range inside `data` (not checked)."
@inline _vec3(data::Vector{Float64}, a::Int) = @inbounds Vec3{Float64}(data[a], data[a+1], data[a+2])

include("spk_types/numeric.jl")     # Taylor extras, implicit derivatives, vector helpers
include("spk_types/interp.jl")      # types 8, 9, 12, 13, 18, 19
include("spk_types/mda.jl")         # types 1, 21
include("spk_types/chebyshev.jl")   # types 14, 20
include("spk_types/conic.jl")       # types 5, 15, 17

"""
    _spk_state(data, seg, et::T) -> (pos::Vec3{T}, vel::Vec3{T})

Position (km) and velocity (km/s) of a segment of one of the types of this file at
`et`, by the algorithm of CSPICE; other types throw `SpiceError(:unsupported)`.

Contract: as `spk_eval`. Every discrete decision uses `primal(et)`, so with a
series-valued `et` both results are the series of the CSPICE position and velocity
functions of the selected record or window.
"""
function _spk_state(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    dtype = seg.dtype
    if dtype == 9 || dtype == 13
        return _spk_unequal(data, seg, et)
    elseif dtype == 8 || dtype == 12
        return _spk_equal(data, seg, et)
    elseif dtype == 18
        return _spk18(data, seg, et)
    elseif dtype == 19
        return _spk19(data, seg, et)
    elseif dtype == 1 || dtype == 21
        return _spk_mda(data, seg, et)
    elseif dtype == 14
        return _spk14(data, seg, et)
    elseif dtype == 20
        return _spk20(data, seg, et)
    elseif dtype == 5
        return _spk5(data, seg, et)
    elseif dtype == 15
        return _spk15(data, seg, et)
    elseif dtype == 17
        return _spk17(data, seg, et)
    elseif dtype == 10
        _t10_needs_lsk(seg)     # type 10 needs a kernel set: spk_types/type10.jl
    else
        _unsupported("SPK", dtype)
    end
end

"Plain values of the components of `x`, whose elements are `Taylor{M,T}`."
@inline _values(x::Vec3{Taylor{M,T}}) where {M,T} = Vec3{T}(x[1].v, x[2].v, x[3].v)

"Tuple of the ordinary derivatives of orders 0 to M of `x`, as vectors."
@inline function _vec_derivatives(x::Vec3{Taylor{M,T}}) where {M,T}
    dx, dy, dz = derivatives(x[1]), derivatives(x[2]), derivatives(x[3])
    return ntuple(k -> Vec3{T}(dx[k], dy[k], dz[k]), Val(M + 1))
end

"""
    _spk_more(data, seg, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}

`spk_eval` for the types of this file: order 0 is the CSPICE position, order 1 the
CSPICE velocity, and order k >= 2 is derivative k-1 of the CSPICE velocity function.
No allocation; the epoch type must be a bits type.
"""
@inline _spk_more(data::Vector{Float64}, seg::Segment, et::T, ::Val{0}) where {T<:Real} =
    (_spk_state(data, seg, et)[1],)
@inline _spk_more(data::Vector{Float64}, seg::Segment, et::T, ::Val{1}) where {T<:Real} =
    _spk_state(data, seg, et)
@inline function _spk_more(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) where {T<:Real,N}
    p, v = _spk_state(data, seg, taylor_var(et, Val(N - 1)))
    return (_values(p), _vec_derivatives(v)...)
end

"""
    _spk_more_position(data, seg, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}

`spk_eval_position` for the types of this file: every order is a derivative of the
CSPICE position function. No allocation.
"""
@inline _spk_more_position(data::Vector{Float64}, seg::Segment, et::T, ::Val{0}) where {T<:Real} =
    (_spk_state(data, seg, et)[1],)
@inline function _spk_more_position(data::Vector{Float64}, seg::Segment, et::T,
        ::Val{N}) where {T<:Real,N}
    p, _ = _spk_state(data, seg, taylor_var(et, Val(N)))
    return _vec_derivatives(p)
end
