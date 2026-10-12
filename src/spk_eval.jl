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

"Tuple `s .+ c .* t`; `s` and `t` have the same length."
@inline _axpy(s::NTuple{M,T}, c::Float64, t::NTuple{M,T}) where {M,T} = map((a, b) -> a + c * b, s, t)

"""
    _cheb3(data, a, ncoef, tau::T, scale, ::Val{N}) -> NTuple{N+1,Vec3{T}}

Sum of the three Chebyshev series whose coefficients are `data[a:a+ncoef-1]`,
`data[a+ncoef:a+2ncoef-1]` and `data[a+2ncoef:a+3ncoef-1]` at `tau`, with derivatives
to order N; derivative `m` is multiplied by `scale^m` (`scale = d tau / d t`).

Requires `ncoef >= 1` and `a + 3ncoef - 1 <= length(data)` (not checked). One pass
over the coefficients with the recurrence
`T_k^(m) = 2 m T_{k-1}^(m-1) + 2 tau T_{k-1}^(m) - T_{k-2}^(m)`.
"""
@inline function _cheb3(data::Vector{Float64}, a::Int, ncoef::Int, tau::T, scale::Float64,
        ::Val{N}) where {T,N}
    M = Val(N + 1)
    o = one(T)
    z = zero(T)
    tau2 = 2 * tau
    p2 = ntuple(m -> m == 1 ? o : z, M)                      # T_0 and derivatives
    @inbounds cx, cy, cz = data[a], data[a+ncoef], data[a+2*ncoef]
    sx = _axpy(map(zero, p2), cx, p2)
    sy = _axpy(map(zero, p2), cy, p2)
    sz = _axpy(map(zero, p2), cz, p2)
    if ncoef > 1
        p1 = ntuple(m -> m == 1 ? tau : m == 2 ? o : z, M)   # T_1 and derivatives
        @inbounds cx, cy, cz = data[a+1], data[a+ncoef+1], data[a+2*ncoef+1]
        sx = _axpy(sx, cx, p1)
        sy = _axpy(sy, cy, p1)
        sz = _axpy(sz, cz, p1)
        # invariant: p1 = T_{k-1}, p2 = T_{k-2} with derivatives to order N, and the
        # sums hold the terms 0 .. k-1
        for k in 2:ncoef-1
            q1, q2 = p1, p2
            p = ntuple(m -> m == 1 ? tau2 * q1[1] - q2[1] :
                (2 * (m - 1)) * q1[m-1] + tau2 * q1[m] - q2[m], M)
            @inbounds cx, cy, cz = data[a+k], data[a+ncoef+k], data[a+2*ncoef+k]
            sx = _axpy(sx, cx, p)
            sy = _axpy(sy, cy, p)
            sz = _axpy(sz, cz, p)
            p2 = p1
            p1 = p
        end
    end
    fx, fy, fz = sx, sy, sz
    return ntuple(m -> Vec3{T}(fx[m], fy[m], fz[m]) * scale^(m - 1), M)
end

"Jet of order N of the first three coefficient sets (`nsets` sets in a record) of the record of `seg` at `et`; `skip` sets are passed over first."
@inline function _cheb_eval(data::Vector{Float64}, seg::Segment, et::T, nsets::Int,
        skip::Int, v::Val{N}) where {T,N}
    addr = _cheb_record(data, seg, et)
    ncoef = (Int(seg.rsize) - 2) ÷ nsets
    @inbounds mid, radius = data[addr], data[addr+1]
    tau = (et - mid) / radius
    return _cheb3(data, addr + 2 + skip * ncoef, ncoef, tau, 1 / radius, v)
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
velocity series. Other types throw `SpiceError(:unsupported, ...)`. No allocation.
"""
@inline function spk_eval(data::Vector{Float64}, seg::Segment, et::T, v::Val{N}) where {T<:Real,N}
    if seg.dtype == 2
        return _cheb_eval(data, seg, et, 3, 0, v)
    elseif seg.dtype == 3
        return _spk3(data, seg, et, v)
    else
        _unsupported("SPK", seg.dtype)
    end
end

"""
    spk_eval_position(data::Vector{Float64}, seg::Segment, et::T, ::Val{N}) -> NTuple{N+1,Vec3{T}}

As `spk_eval`, but every order is a derivative of the position series (the result
differs from `spk_eval` only for type 3). Same contract as `spk_eval`.
"""
@inline function spk_eval_position(data::Vector{Float64}, seg::Segment, et::T, v::Val{N}) where {T<:Real,N}
    if seg.dtype == 2
        return _cheb_eval(data, seg, et, 3, 0, v)
    elseif seg.dtype == 3
        return _cheb_eval(data, seg, et, 6, 0, v)
    else
        _unsupported("SPK", seg.dtype)
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
