# Interpolation types: Lagrange (8 equal step, 9 unequal step), Hermite (12 equal
# step, 13 unequal step), ESOC/DDID (18) and ESOC/DDID piecewise (19).
#
# The window of states is read directly from the data array and the divided-difference
# tables live in fixed-size stack arrays, so a query does not allocate.

"Largest number of points in one interpolation window (degree 27 Lagrange, `MAXDEG` of CSPICE)."
const _MAXWIN = 28

"Value of `xbase` that selects the abscissas 1, 2, 3, ... of the equal-step routines."
const _GRID = -1

"Abscissa `i` of a window: `data[xbase+i]`, or the whole number `i` if `xbase` is `_GRID`."
@inline _abscissa(data::Vector{Float64}, xbase::Int, i::Int) =
    xbase < 0 ? Float64(i) : @inbounds data[xbase+i]

"""
    _lagrange(seg, data, pbase, n, xbase, x::T) -> SVector{6,T}

Value at `x` of the six Lagrange polynomials through the `n` six-element packets
`data[pbase+6(k-1)+1 : pbase+6k]`, `k = 1:n`, at the abscissas `_abscissa(data, xbase, k)`,
by the Neville recurrence of `lgrint` (of `lgresp` for `xbase == _GRID`, where `x` is
the scaled abscissa `(t - first)/step + 1`).

Requires `1 <= n <= _MAXWIN` and the packets and abscissas inside `data` (not
checked). Throws `SpiceError(:baddaf)` if two abscissas are equal.
"""
function _lagrange(seg::Segment, data::Vector{Float64}, pbase::Int, n::Int, xbase::Int,
        x::T) where {T}
    work = MVector{_MAXWIN,SVector{6,T}}(undef)
    @inbounds for i in 1:n
        a = pbase + 6 * (i - 1)
        work[i] = SVector{6,T}(data[a+1], data[a+2], data[a+3], data[a+4], data[a+5], data[a+6])
    end
    # invariant: after pass j, work[i] is the value at x of the polynomial through the
    # points i .. i+j, for i = 1 .. n-j
    @inbounds for j in 1:n-1
        for i in 1:n-j
            xi = _abscissa(data, xbase, i)
            xij = _abscissa(data, xbase, i + j)
            denom = xi - xij
            denom == 0.0 && _badsegment(seg, "two epochs of one interpolation window are equal")
            c1 = x - xij
            c2 = xi - x
            work[i] = (c1 * work[i] + c2 * work[i+1]) / denom
        end
    end
    return @inbounds work[1]
end

"""
    _hermite(seg, data, pbase, psize, off, n, xbase, x::T, dscale) -> (f::Vec3{T}, df::Vec3{T})

Value `f` and derivative `df` at `x` of the three Hermite polynomials that match, at
the abscissa `_abscissa(data, xbase, k)`, the values `data[a+1:a+3]` and the
derivatives `data[a+4:a+6] * dscale`, with `a = pbase + psize*(k-1) + off`, `k = 1:n`.
The recurrence is the one of `hrmint`; with `xbase == _GRID`, `x` the scaled abscissa
and `dscale` the step it is the one of `hrmesp`, and `df` is then the derivative with
respect to the scaled abscissa.

Requires `1 <= n <= _MAXWIN` and the packets and abscissas inside `data` (not
checked). Throws `SpiceError(:baddaf)` if two abscissas are equal.
"""
function _hermite(seg::Segment, data::Vector{Float64}, pbase::Int, psize::Int, off::Int,
        n::Int, xbase::Int, x::T, dscale::Float64) where {T}
    w1 = MVector{2 * _MAXWIN,Vec3{T}}(undef)      # column 1 of WORK: values
    w2 = MVector{2 * _MAXWIN,Vec3{T}}(undef)      # column 2 of WORK: derivatives
    @inbounds for k in 1:n
        a = pbase + psize * (k - 1) + off
        w1[2k-1] = Vec3{T}(data[a+1], data[a+2], data[a+3])
        w1[2k] = Vec3{T}(data[a+4] * dscale, data[a+5] * dscale, data[a+6] * dscale)
    end
    # First column of the table: the polynomials through two consecutive entries of
    # the doubled abscissa list (x1, x1, x2, x2, ...).
    @inbounds for i in 1:n-1
        xi = _abscissa(data, xbase, i)
        xn = _abscissa(data, xbase, i + 1)
        c1 = xn - x
        c2 = x - xi
        denom = xn - xi
        denom == 0.0 && _badsegment(seg, "two epochs of one interpolation window are equal")
        prev = 2i - 1
        this = prev + 1
        next = this + 1
        w2[prev] = w1[this]
        w2[this] = (w1[next] - w1[prev]) / denom
        temp = w1[this] * (x - xi) + w1[prev]
        w1[this] = (c1 * w1[prev] + c2 * w1[next]) / denom
        w1[prev] = temp
    end
    @inbounds begin
        w2[2n-1] = w1[2n]
        w1[2n-1] = w1[2n] * (x - _abscissa(data, xbase, n)) + w1[2n-1]
    end
    # invariant: after pass j, w1[i] and w2[i] are the value and the derivative at x of
    # the polynomial through the entries i .. i+j of the doubled list
    @inbounds for j in 2:2n-1
        for i in 1:2n-j
            xi = _abscissa(data, xbase, (i + 1) ÷ 2)
            xij = _abscissa(data, xbase, (i + j + 1) ÷ 2)
            c1 = xij - x
            c2 = x - xi
            denom = xij - xi
            denom == 0.0 && _badsegment(seg, "two epochs of one interpolation window are equal")
            w2[i] = (c1 * w2[i] + c2 * w2[i+1] + (w1[i+1] - w1[i])) / denom
            w1[i] = (c1 * w1[i] + c2 * w1[i+1]) / denom
        end
    end
    return @inbounds (w1[1], w2[1])
end

@inline _split6(s::SVector{6,T}) where {T} = (Vec3{T}(s[1], s[2], s[3]), Vec3{T}(s[4], s[5], s[6]))

"""
    _spk_equal(data, seg, et::T) -> (pos, vel)

State of a type 8 (Lagrange) or type 12 (Hermite) segment: states at equal steps, as
`spkr08` with `spke08` or `spke12`. The segment ends with the epoch of the first
state, the step, the window size minus 1 and the number of states.

A window with an odd number of points is centered on the state nearest to the epoch,
one with an even number on the interval that contains it, and the window moves
inward at the two ends of the segment. Type 12 gives the derivative of the position
polynomial as the velocity.
"""
function _spk_equal(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    len >= 10 || _badsegment(seg, "the segment is too short")
    @inbounds start, step = data[e-3], data[e-2]
    (isfinite(start) && isfinite(step) && step != 0.0) ||
        _badsegment(seg, "the first epoch or the step is not valid")
    window = _seg_count(seg, @inbounds(data[e-1]), 0, _MAXWIN - 1, "window size - 1") + 1
    n = _seg_count(seg, @inbounds(data[e]), window, _MAXCOUNT, "state count")
    6 * n + 4 == len || _badsegment(seg, "the state count does not agree with the length")
    degree = window - 1
    q = (primal(et) - start) / step
    k = (isodd(window) ? _round_index(q) : _trunc_index(q)) + 1
    first = min(max(1, k - degree ÷ 2), n - degree)
    pbase = b - 1 + 6 * (first - 1)
    x = (et - (start + (first - 1) * step)) / step + 1.0
    if seg.dtype == 8
        return _split6(_lagrange(seg, data, pbase, window, _GRID, x))
    else
        f, df = _hermite(seg, data, pbase, 6, 0, window, _GRID, x, step)
        return f, df / step
    end
end

"""
    _spk_unequal(data, seg, et::T) -> (pos, vel)

State of a type 9 (Lagrange) or type 13 (Hermite) segment: states at arbitrary
epochs, as `spkr09` with `spke09` or `spke13`. The segment holds the states, the
epochs, an epoch directory, the window size minus 1 and the number of states.
"""
function _spk_unequal(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    len >= 9 || _badsegment(seg, "the segment is too short")
    window = _seg_count(seg, @inbounds(data[e-1]), 0, _MAXWIN - 1, "window size - 1") + 1
    n = _seg_count(seg, @inbounds(data[e]), window, _MAXCOUNT, "state count")
    7 * n + (n - 1) ÷ 100 + 2 == len ||
        _badsegment(seg, "the state count does not agree with the length")
    degree = window - 1
    ebase = b - 1 + 6 * n                       # epoch i is data[ebase+i]
    t = primal(et)
    low = max(_lstltd(data, ebase, n, t), 1)    # last epoch before t, or the first epoch
    center = low
    if isodd(window) && low < n
        @inbounds center = abs(t - data[ebase+low]) < abs(t - data[ebase+low+1]) ? low : low + 1
    end
    first = min(max(center - degree ÷ 2, 1), n - degree)
    pbase = b - 1 + 6 * (first - 1)
    xbase = ebase + first - 1
    if seg.dtype == 9
        return _split6(_lagrange(seg, data, pbase, window, xbase, et))
    else
        return _hermite(seg, data, pbase, 6, 0, window, xbase, et, 1.0)
    end
end

"""
    _packet_window(data, ebase, n, wnd, t) -> (first, count)

Index of the first packet and number of packets of the window of nominal size `wnd`
for the epoch `t` among the `n` epochs `data[ebase+1:ebase+n]`, as `spkr18` and
`spkr19` select it: `wnd/2` packets on each side of the epoch, fewer at an end of
the packet set. Postcondition: `1 <= first` and `first + count - 1 <= n`.
"""
@inline function _packet_window(data::Vector{Float64}, ebase::Int, n::Int, wnd::Int, t::Float64)
    low = max(_lstltd(data, ebase, n, t), 1)
    lsize = min(wnd ÷ 2, low)
    rsize = min(wnd ÷ 2, n - low)
    return low - lsize + 1, lsize + rsize
end

"Packet size of a subtype of types 18 and 19: 12 for Hermite with two derivative sets, else 6."
@inline _packet_size(subtyp::Int) = subtyp == 0 ? 12 : 6

"""
    _packet_state(seg, data, subtyp, pbase, n, xbase, et::T) -> (pos, vel)

State from a window of `n` packets that start after `data[pbase]`, with epochs
`data[xbase+1:xbase+n]`, as `spke18` and `spke19`. Subtype 0: packets of position,
velocity, velocity and acceleration; the position comes from the Hermite polynomial
of the first pair and the velocity from the one of the second pair. Subtype 1:
Lagrange polynomials of the six state components. Subtype 2: Hermite polynomial of
position and velocity, with its derivative as the velocity.
"""
function _packet_state(seg::Segment, data::Vector{Float64}, subtyp::Int, pbase::Int, n::Int,
        xbase::Int, et::T) where {T<:Real}
    if subtyp == 0
        pos, _ = _hermite(seg, data, pbase, 12, 0, n, xbase, et, 1.0)
        vel, _ = _hermite(seg, data, pbase, 12, 6, n, xbase, et, 1.0)
        return pos, vel
    elseif subtyp == 1
        return _split6(_lagrange(seg, data, pbase, n, xbase, et))
    else
        return _hermite(seg, data, pbase, 6, 0, n, xbase, et, 1.0)
    end
end

"""
    _spk18(data, seg, et::T) -> (pos, vel)

State of a type 18 segment, as `spkr18` and `spke18`. The segment holds the packets,
the epochs, an epoch directory, the subtype (0 or 1), the window size and the number
of packets. The window size is even and not above 8 (subtype 0) or 16 (subtype 1).
"""
function _spk18(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    len >= 10 || _badsegment(seg, "the segment is too short")
    subtyp = _seg_count(seg, @inbounds(data[e-2]), 0, 1, "subtype")
    wnd = _seg_count(seg, @inbounds(data[e-1]), 2, subtyp == 0 ? 8 : 16, "window size")
    iseven(wnd) || _badsegment(seg, "the window size is odd")
    n = _seg_count(seg, @inbounds(data[e]), 2, _MAXCOUNT, "packet count")
    psize = _packet_size(subtyp)
    n * (psize + 1) + (n - 1) ÷ 100 + 3 == len ||
        _badsegment(seg, "the packet count does not agree with the length")
    ebase = b - 1 + n * psize
    first, count = _packet_window(data, ebase, n, wnd, primal(et))
    return _packet_state(seg, data, subtyp, b - 1 + (first - 1) * psize, count,
        ebase + first - 1, et)
end

"""
    _spk19(data, seg, et::T) -> (pos, vel)

State of a type 19 segment, as `spkr19` and `spke19`. The segment is a sequence of
mini-segments (each one a packet set with epochs, an epoch directory, a subtype 0, 1
or 2, a window size and a packet count), then the `n+1` bounds of the interpolation
intervals, a directory of the bounds, `n+1` pointers to the mini-segments relative
to the segment start, the selection flag and `n`.

For an epoch that is a bound between two intervals the flag selects the later
interval (1) or the earlier one (0). The window size is even and not above 14
(subtype 0) or 28 (subtypes 1 and 2).
"""
function _spk19(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    len >= 6 || _badsegment(seg, "the segment is too short")
    later = @inbounds(data[e-1]) == 1.0         # the writer of CSPICE stores -1 for "earlier"
    n = _seg_count(seg, @inbounds(data[e]), 1, _MAXCOUNT, "interval count")
    pbas = e - 2 - (n + 1)                      # pointer i is data[pbas+i]
    ivbas = pbas - n ÷ 100 - (n + 1)            # interval bound i is data[ivbas+i]
    # The last pointer is the place after the last mini-segment: the first bound.
    (ivbas >= b && @inbounds(data[pbas+n+1]) == ivbas - b + 2) ||
        _badsegment(seg, "the interval count does not agree with the length")
    t = primal(et)
    mini = later ? min(_lstled(data, ivbas, n + 1, t), n) : _lstltd(data, ivbas, n + 1, t)
    mini = clamp(mini, 1, n)
    @inbounds p1, p2 = data[pbas+mini], data[pbas+mini+1]
    minib = _seg_count(seg, p1, 1, len, "mini-segment pointer") + b - 1
    minie = _seg_count(seg, p2, 1, len, "mini-segment pointer") + b - 2
    (minie - minib >= 3 && minie <= ivbas) ||
        _badsegment(seg, "the pointers of mini-segment $mini are not valid")
    subtyp = _seg_count(seg, @inbounds(data[minie-2]), 0, 2, "subtype")
    wnd = _seg_count(seg, @inbounds(data[minie-1]), 2, subtyp == 0 ? 14 : 28, "window size")
    iseven(wnd) || _badsegment(seg, "the window size is odd")
    npkt = _seg_count(seg, @inbounds(data[minie]), 1, _MAXCOUNT, "packet count")
    psize = _packet_size(subtyp)
    npkt * (psize + 1) + (npkt - 1) ÷ 100 + 3 == minie - minib + 1 ||
        _badsegment(seg, "the packet count of mini-segment $mini does not agree with its length")
    ebase = minib - 1 + npkt * psize
    first, count = _packet_window(data, ebase, npkt, wnd, t)
    return _packet_state(seg, data, subtyp, minib - 1 + (first - 1) * psize, count,
        ebase + first - 1, et)
end
