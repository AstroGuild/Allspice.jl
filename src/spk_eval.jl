# Segment search and evaluation of the fixed-interval Chebyshev data types
# (SPK types 2 and 3, binary PCK type 2), with time derivatives of any order.

"""
    find_segment(segs::Vector{Segment}, et::Float64) -> Int

Index of the first segment of `segs` (which is in priority order) with
`t0 <= et <= t1`, or 0 if no segment covers `et`.
"""
function find_segment(segs::Vector{Segment}, et::Float64)
    @inbounds for i in eachindex(segs)   # invariant: no segment before `i` covers `et`
        s = segs[i]
        (s.t0 <= et <= s.t1) && return i
    end
    return 0
end

"""
    _cheb_record(data, seg, et) -> Int

Address of the first double (MID) of the record of `seg` that covers `primal(et)`, as
`spkr02` selects it: the record index comes from INIT and INTLEN and is clamped to
`1:nrec`, so the segment end epoch falls in the last record. Throws
`SpiceError(:baddaf)` if the record is not inside `data`.
"""
@inline function _cheb_record(data::Vector{Float64}, seg::Segment, et)
    k = unsafe_trunc(Int, (primal(et) - seg.init) / seg.intlen)
    k = clamp(k, 0, Int(seg.nrec) - 1)
    addr = seg.baddr + k * Int(seg.rsize)
    (addr >= 1 && addr + Int(seg.rsize) - 1 <= length(data)) ||
        throw(SpiceError(:baddaf, "segment record is outside the data array"))
    return addr
end

"""
    _clenshaw_step(c, s2::T, b1::NTuple{M,T}, b2::NTuple{M,T}) -> NTuple{M,T}

One step of the Clenshaw recurrence with derivatives: from `b1 = b_{k+1}` and
`b2 = b_{k+2}` (each with derivatives of orders 0 to M-1 with respect to `s`) return
`b_k` for the coefficient `c`, where `s2 = 2s`:
`b_k^(m) = 2m b_{k+1}^(m-1) + 2s b_{k+1}^(m) - b_{k+2}^(m)`, plus `c` for `m = 0`.
Orders 0 and 1 use the operation order of `chbint_` of CSPICE.
"""
@inline _clenshaw_step(c::Real, s2::T, b1::NTuple{M,T}, b2::NTuple{M,T}) where {M,T} =
    ntuple(m -> m == 1 ? c + (s2 * b1[1] - b2[1]) :
        (2 * (m - 1)) * b1[m-1] + b1[m] * s2 - b2[m], Val(M))

"""
    _clenshaw_end(c0, s::T, b1::NTuple{M,T}, b2::NTuple{M,T}) -> NTuple{M,T}

Value and derivatives (orders 0 to M-1 with respect to `s`) of the Chebyshev series
from the last two Clenshaw terms `b1 = b_1` and `b2 = b_2` and the coefficient `c0`:
`f^(m) = m b_1^(m-1) + s b_1^(m) - b_2^(m)`, plus `c0` for `m = 0`.
"""
@inline _clenshaw_end(c0::Real, s::T, b1::NTuple{M,T}, b2::NTuple{M,T}) where {M,T} =
    ntuple(m -> m == 1 ? c0 + (s * b1[1] - b2[1]) :
        (m - 1) * b1[m-1] + s * b1[m] - b2[m], Val(M))

"""
    _cheb3(data, a, ncoef, s::T, radius, ::Val{N}) -> NTuple{N+1,Vec3{T}}

Sum of the three Chebyshev series whose coefficients are `data[a:a+ncoef-1]`,
`data[a+ncoef:a+2ncoef-1]` and `data[a+2ncoef:a+3ncoef-1]` at `s`, with derivatives
to order N; derivative `m` is divided `m` times by `radius` (`d s / d t = 1 / radius`).

Requires `ncoef >= 1` and `a + 3ncoef - 1 <= length(data)` (not checked). One pass
over the coefficients, highest first, with the Clenshaw recurrence. For a plain epoch
the value and the first derivative come from the same floating-point operations as in
`chbval_` and `chbint_` of CSPICE, so they are equal to the last bit; this matters for
an angle of many turns, where a difference of one unit in the last place is large
relative to the reduced angle.
"""
@inline function _cheb3(data::AbstractVector{<:Real}, a::Int, ncoef::Int, s::T, radius::Real,
        ::Val{N}) where {T,N}
    M = Val(N + 1)
    z = ntuple(_ -> zero(T), M)
    s2 = s * 2
    x1, x2 = z, z
    y1, y2 = z, z
    z1, z2 = z, z
    # invariant: (x1, x2), (y1, y2), (z1, z2) are the Clenshaw terms b_{k+1}, b_{k+2}
    # of the three series, with derivatives to order N
    for k in ncoef-1:-1:1
        @inbounds cx, cy, cz = data[a+k], data[a+ncoef+k], data[a+2*ncoef+k]
        x1, x2 = _clenshaw_step(cx, s2, x1, x2), x1
        y1, y2 = _clenshaw_step(cy, s2, y1, y2), y1
        z1, z2 = _clenshaw_step(cz, s2, z1, z2), z1
    end
    @inbounds cx, cy, cz = data[a], data[a+ncoef], data[a+2*ncoef]
    fx = _clenshaw_end(cx, s, x1, x2)
    fy = _clenshaw_end(cy, s, y1, y2)
    fz = _clenshaw_end(cz, s, z1, z2)
    return ntuple(M) do m
        v = Vec3{T}(fx[m], fy[m], fz[m])
        for _ in 2:m      # invariant: `v` has been divided by `radius` once per pass
            v = v / radius
        end
        v
    end
end

"Jet of order N of the first three coefficient sets (`nsets` sets in a record) of the record of `seg` at `et`; `skip` sets are passed over first."
@inline function _cheb_eval(data::Vector{Float64}, seg::Segment, et::T, nsets::Int,
        skip::Int, v::Val{N}) where {T,N}
    addr = _cheb_record(data, seg, et)
    ncoef = (Int(seg.rsize) - 2) ÷ nsets
    @inbounds mid, radius = data[addr], data[addr+1]
    s = (et - mid) / radius
    return _cheb3(data, addr + 2 + skip * ncoef, ncoef, s, radius, v)
end

@inline _spk3(data, seg, et, ::Val{0}) = _cheb_eval(data, seg, et, 6, 0, Val(0))
@inline function _spk3(data, seg, et, ::Val{N}) where {N}
    p = _cheb_eval(data, seg, et, 6, 0, Val(0))
    v = _cheb_eval(data, seg, et, 6, 3, Val(N - 1))
    return (p[1], v...)
end

@noinline _unsupported(what, dtype) =
    throw(SpiceError(:unsupported, "$what data type $dtype is not supported"))

"""
    spk_eval(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}

Position (km) of `seg.target` relative to `seg.center` in frame `seg.frame` at `et`
and its time derivatives to order N; order 1 is the velocity that CSPICE returns.

Contract: `seg` is a descriptor from `open_daf` and `data` is the array of the same
file; `seg.t0 <= primal(et) <= seg.t1` (not checked: outside it the nearest record is
extrapolated). Type 2: all orders are derivatives of the position series. Type 3:
order 0 comes from the position series and order `k >= 1` is derivative `k-1` of the
velocity series. Types 1, 5, 8, 9, 12 to 15 and 17 to 21 (`spk_more.jl`): order
`k >= 2` is derivative `k-1` of the velocity function that CSPICE defines. Type 10
and unknown types throw `SpiceError(:unsupported, ...)`. No allocation.
"""
@inline function spk_eval(data::Vector{Float64}, seg::Segment, et::T, v::Val{N}) where {T<:Real,N}
    if seg.dtype == 2
        return _cheb_eval(data, seg, et, 3, 0, v)
    elseif seg.dtype == 3
        return _spk3(data, seg, et, v)
    else
        return _spk_more(data, seg, et, v)      # the other types: spk_more.jl
    end
end

"""
    spk_eval_position(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}

As `spk_eval`, but every order is a derivative of the position function (the result
differs from `spk_eval` for the types whose velocity is not the derivative of the
position in CSPICE, such as 3, 8, 9 and 14). Same contract as `spk_eval`.
"""
@inline function spk_eval_position(data::Vector{Float64}, seg::Segment, et::T, v::Val{N}) where {T<:Real,N}
    if seg.dtype == 2
        return _cheb_eval(data, seg, et, 3, 0, v)
    elseif seg.dtype == 3
        return _cheb_eval(data, seg, et, 6, 0, v)
    else
        return _spk_more_position(data, seg, et, v)
    end
end

"""
    pck_eval(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}

Euler angles (radians) of a binary PCK type 2 segment at `et` and their time
derivatives to order N, as `pcke02` of CSPICE gives them.

Contract: `seg` is a PCK descriptor from `open_daf` and `data` is the array of the
same file; `seg.t0 <= primal(et) <= seg.t1` (not checked). The three components are
the stored 3-1-3 angles `(PHI, DELTA, W)` with `PHI = RA + pi/2` and
`DELTA = pi/2 - DEC`, so the matrix from the inertial frame `seg.frame` to the body
frame is `R3(W) * R1(DELTA) * R3(PHI)`. Other types throw
`SpiceError(:unsupported, ...)`. No allocation.
"""
@inline function pck_eval(data::Vector{Float64}, seg::Segment, et::T, v::Val{N}) where {T<:Real,N}
    seg.dtype == 2 || _unsupported("binary PCK", seg.dtype)
    return _cheb_eval(data, seg, et, 3, 0, v)
end
