# Windows: sets of disjoint closed intervals of the real line, as the `wn*` routines of
# CSPICE define them. A window is a plain Julia object; there is no SPICE cell and no
# size limit, so the error SPICE(WINDOWEXCESS) cannot occur.

"""
    Window() ; Window(endpoints)

Window of double-precision intervals: `data` holds the endpoints
`left₁, right₁, left₂, right₂, …` of disjoint closed intervals in increasing order
(`leftᵢ <= rightᵢ < leftᵢ₊₁`). A `Window` is an `AbstractVector{Float64}` of its
endpoints, so `w[k]`, `length(w)` and iteration give what a SPICE double cell gives.

Contract: `Window(endpoints)` copies `endpoints` and validates them with `wnvald!`
(the pairs may be in any order and may overlap); it throws `SpiceError` with kind
`:unmatchendpts` for an odd count and `:badendpoints` for a pair with left > right.
A `Window` is not safe for concurrent mutation; a caller that shares one between
threads must lock it.
"""
struct Window <: AbstractVector{Float64}
    data::Vector{Float64}
    Window() = new(Float64[])
    Window(endpoints::AbstractVector{<:Real}) = wnvald!(new(Float64[x for x in endpoints]))
end

Base.size(w::Window) = size(w.data)
Base.getindex(w::Window, i::Int) = w.data[i]
Base.IndexStyle(::Type{Window}) = IndexLinear()
Base.copy(w::Window) = (c = Window(); append!(c.data, w.data); c)
Base.empty!(w::Window) = (empty!(w.data); w)
Base.:(==)(a::Window, b::Window) = a.data == b.data

"Window that owns the valid endpoint vector `d` (no check, no copy)."
function _window_of(d::Vector{Float64})
    w = Window()
    append!(w.data, d)
    return w
end

"""
    card(x) -> Int

Number of elements of `x`, as `card_c` gives for a SPICE cell: the number of endpoints
of a `Window` (twice the number of intervals) and the length of a vector of IDs.
"""
card(x::AbstractVector) = length(x)

"""
    wncard(window) -> Int

Number of intervals of `window`, as `wncard_c`. Postcondition: `2 * wncard(w) == card(w)`.
"""
wncard(w::Window) = length(w.data) ÷ 2

"""
    wninsd!(window, left, right) -> window

Insert the interval [`left`, `right`] into `window` and merge it with every interval
that it overlaps or touches, as `wninsd_c`. Precondition: `left <= right`, else
`SpiceError(:badendpoints)`. Postcondition: `window` is a valid window that is the
union of its old value and the interval.
"""
function wninsd!(w::Window, left::Real, right::Real)
    l = Float64(left)
    r = Float64(right)
    l > r && throw(SpiceError(:badendpoints,
        "wninsd!: left endpoint $l is greater than right endpoint $r"))
    d = w.data
    card = length(d)
    if card == 0 || l > d[card]
        push!(d, l, r)
        return w
    end
    # First interval whose right endpoint is not less than `l`.
    i = 2
    while i <= card && d[i] < l          # invariant: intervals before i end before l
        i += 2
    end
    if r < d[i-1]
        # The new interval lies strictly between two intervals (or before the first).
        splice!(d, i-1:i-2, (l, r))
        return w
    end
    d[i-1] = min(l, d[i-1])
    d[i] = max(r, d[i])
    j = i + 2
    while j <= card && d[j] <= d[i]      # invariant: intervals i+2 .. j-2 are swallowed
        j += 2
    end
    if j <= card && d[i] >= d[j-1]
        d[i] = d[j]
        j += 2
    end
    while j <= card                      # invariant: d[1:i] is final; j is the next to keep
        i += 2
        d[i-1] = d[j-1]
        d[i] = d[j]
        j += 2
    end
    resize!(d, i)
    return w
end

"""
    wnvald!(window) -> window

Turn the endpoint pairs in `window.data` (any order, overlaps allowed) into a valid
window, as `wnvald_c`. Throws `SpiceError(:unmatchendpts)` for an odd number of
endpoints and `SpiceError(:badendpoints)` for a pair with left > right.
"""
function wnvald!(w::Window)
    pairs = copy(w.data)
    n = length(pairs)
    isodd(n) && throw(SpiceError(:unmatchendpts,
        "wnvald!: the window has an odd number of endpoints ($n)"))
    for i in 1:2:n
        pairs[i] > pairs[i+1] && throw(SpiceError(:badendpoints,
            "wnvald!: left endpoint $(pairs[i]) is greater than right endpoint $(pairs[i+1])"))
    end
    empty!(w.data)
    for i in 1:2:n                       # invariant: w is a valid window of pairs 1 .. i-1
        wninsd!(w, pairs[i], pairs[i+1])
    end
    return w
end

"""
    wnunid(a, b) -> Window

Union of the windows `a` and `b`, as `wnunid_c`.
"""
function wnunid(a::Window, b::Window)
    A, B = a.data, b.data
    acard, bcard = length(A), length(B)
    c = Float64[]
    ap = 1
    bp = 1
    stop = -Inf
    while ap < acard || bp < bcard       # invariant: c is the union of the intervals taken
        usea = bp > bcard || (ap <= acard && A[ap] < B[bp])
        lo, hi = usea ? (A[ap], A[ap+1]) : (B[bp], B[bp+1])
        if lo <= stop && !isempty(c)
            c[end] = max(c[end], hi)
        else
            push!(c, lo, hi)
        end
        stop = c[end]
        usea ? (ap += 2) : (bp += 2)
    end
    return _window_of(c)
end

"""
    wnintd(a, b) -> Window

Intersection of the windows `a` and `b`, as `wnintd_c`. Two intervals that share one
endpoint give the singleton interval at that point.
"""
function wnintd(a::Window, b::Window)
    A, B = a.data, b.data
    acard, bcard = length(A), length(B)
    c = Float64[]
    ap = 1
    bp = 1
    while ap < acard && bp < bcard       # invariant: c holds the common parts found so far
        if A[ap+1] < B[bp+1]
            A[ap+1] >= B[bp] && push!(c, max(B[bp], A[ap]), A[ap+1])
            ap += 2
        else
            B[bp+1] >= A[ap] && push!(c, max(A[ap], B[bp]), B[bp+1])
            bp += 2
        end
    end
    return _window_of(c)
end

"""
    wndifd(a, b) -> Window

Difference of the windows: the points of `a` that are not interior points of `b`, as
`wndifd_c` computes it (a singleton interval of `b` removes nothing, and the result is
a window of closed intervals).
"""
function wndifd(a::Window, b::Window)
    A, B = a.data, b.data
    acard, bcard = length(A), length(B)
    bcard == 0 && return copy(a)
    c = Float64[]
    apb = 1
    bpb = 1
    while apb + 1 <= acard               # invariant: c is the difference for a's intervals before apb
        f = A[apb]
        l = A[apb+1]
        unrslv = bpb + 1 <= bcard
        keep = true
        while unrslv
            bl = B[bpb]
            br = B[bpb+1]
            if l < bl
                unrslv = false
            elseif f > br
                bpb += 2
                unrslv = bpb + 1 <= bcard
            elseif bl <= f && l <= br
                keep = false
                unrslv = false
            elseif bl <= f
                f = br
                bpb += 2
                unrslv = bpb + 1 <= bcard
            elseif f <= bl && l >= br && bl < br
                push!(c, f, bl)
                f = br
                if f == l
                    keep = false
                    unrslv = false
                end
                bpb += 2
                unrslv = unrslv && bpb + 1 <= bcard
            elseif f <= bl && l >= br && bl == br
                bpb += 2
                unrslv = bpb + 1 <= bcard
            else
                l = bl
                unrslv = false
            end
        end
        keep && push!(c, f, l)
        apb += 2
    end
    return _window_of(c)
end

"""
    wncomd(window, left, right) -> Window

Complement of `window` with respect to the interval [`left`, `right`], as `wncomd_c`
(argument order of SPICE.jl). Precondition: `left <= right`, else
`SpiceError(:badendpoints)`.
"""
function wncomd(w::Window, left::Real, right::Real)
    l = Float64(left)
    r = Float64(right)
    l > r && throw(SpiceError(:badendpoints,
        "wncomd: left endpoint $l is greater than right endpoint $r"))
    d = w.data
    card = length(d)
    result = Window()
    if card == 0 || d[1] >= r || d[card] <= l
        return wninsd!(result, l, r)
    end
    i = 2
    while i <= card && d[i] < l          # invariant: intervals before i end before l
        i += 2
    end
    i <= card && d[i-1] > l && wninsd!(result, l, d[i-1])
    while i < card && d[i+1] < r         # invariant: gaps before interval i are in result
        wninsd!(result, d[i], d[i+1])
        i += 2
    end
    i <= card && d[i] < r && wninsd!(result, d[i], r)
    return result
end

"""
    wnelmd(window, point) -> Bool

True if `point` is in one of the intervals of `window`, as `wnelmd_c` (argument order
of SPICE.jl).
"""
function wnelmd(w::Window, point::Real)
    d = w.data
    for i in 1:2:length(d)-1
        d[i] <= point <= d[i+1] && return true
    end
    return false
end

"""
    wnincd(window, left, right) -> Bool

True if the interval [`left`, `right`] is inside one interval of `window`, as
`wnincd_c` (argument order of SPICE.jl).
"""
function wnincd(w::Window, left::Real, right::Real)
    d = w.data
    for i in 1:2:length(d)-1
        (left >= d[i] && right <= d[i+1]) && return true
    end
    return false
end

"""
    wnfetd(window, n) -> (left, right)

Endpoints of interval number `n` (1-based, as in SPICE.jl) of `window`. Precondition:
`1 <= n <= wncard(window)`, else `SpiceError(:nointerval)`.
"""
function wnfetd(w::Window, n::Integer)
    (1 <= n <= wncard(w)) || throw(SpiceError(:nointerval,
        "wnfetd: interval $n is not in a window of $(wncard(w)) interval(s)"))
    return w.data[2n-1], w.data[2n]
end

"""
    wnsumd(window) -> (meas, avg, stddev, idxsml, idxlon)

Summary of `window`, as `wnsumd_c`: total measure, mean and standard deviation of the
interval lengths, and the 1-based index in `window` of the left endpoint of the
shortest and of the longest interval (both 0 for an empty window).
"""
function wnsumd(w::Window)
    d = w.data
    card = length(d)
    card == 0 && return (0.0, 0.0, 0.0, 0, 0)
    sum = 0.0
    sumsqr = 0.0
    idxsml = 1
    idxlon = 1
    mshort = d[2] - d[1]
    mlong = mshort
    for i in 1:2:card                    # invariant: sums and extremes cover intervals before i
        m = d[i+1] - d[i]
        sum += m
        sumsqr += m * m
        if m < mshort
            idxsml = i
            mshort = m
        end
        if m > mlong
            idxlon = i
            mlong = m
        end
    end
    avg = sum * 2.0 / card
    # Rounding can make the variance of equal intervals slightly negative; C gives NaN there.
    stddev = sqrt(max(sumsqr * 2.0 / card - avg * avg, 0.0))
    return (sum, avg, stddev, idxsml, idxlon)
end
