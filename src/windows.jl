# Windows: sets of disjoint closed intervals of the real line, as the `wn*` routines of
# CSPICE define them. A window is a plain Julia object; there is no SPICE cell and no
# size limit, so the error SPICE(WINDOWEXCESS) cannot occur.

"""
    Window()
    Window(endpoints)

Window: a set of disjoint closed intervals of the real line. It takes the place of the
double-precision window cell of CSPICE.

`Window()` is empty. `Window(endpoints)` takes the pairs `left₁, right₁, left₂, right₂, …`.
The pairs can be in any order and can overlap; they are sorted and merged, as
[`wnvald!`](@ref) does. `endpoints` is copied.

A `Window` is an `AbstractVector{Float64}` of its endpoints in increasing order, so
`w[k]`, `length(w)` and iteration work. The field `data` holds the endpoints. A window
has no size limit, so the CSPICE error `SPICE(WINDOWEXCESS)` cannot occur. A `Window`
is not safe for concurrent mutation.

Throws `SpiceError(:unmatchendpts)` if the number of endpoints is odd, and
`SpiceError(:badendpoints)` if a pair has `left > right`.

# Examples
```jldoctest
julia> w = Window([1, 3, 7, 11, 2, 5])
4-element Window:
  1.0
  5.0
  7.0
 11.0

julia> wncard(w)
2
```

See also [`wninsd!`](@ref), [`wnunid`](@ref), [`wnfetd`](@ref), [`wnsumd`](@ref).
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

"Window with a copy of the endpoint vector `d`, which must already be valid (no check)."
function _window_of(d::Vector{Float64})
    w = Window()
    append!(w.data, d)
    return w
end

"""
    card(x) -> Int

Return the number of elements of `x`, as `card_c` of CSPICE gives for a cell.

For a [`Window`](@ref) this is the number of endpoints, which is twice the number of
intervals. For a vector of IDs it is the length.

See also [`wncard`](@ref).
"""
card(x::AbstractVector) = length(x)

"""
    wncard(window) -> Int

Return the number of intervals of `window`, as `wncard_c` of CSPICE.

See also [`card`](@ref), [`wnfetd`](@ref).
"""
wncard(w::Window) = length(w.data) ÷ 2

"""
    wninsd!(window, left, right) -> window

Insert the interval [`left`, `right`] into `window` and return `window`, as `wninsd_c`
of CSPICE.

The new interval is merged with each interval that it overlaps or touches. The window
is the first argument; CSPICE has it last.

Throws `SpiceError(:badendpoints)` if `left > right`.

# Examples
```jldoctest
julia> w = Window([1, 3, 7, 11]);

julia> wninsd!(w, 3, 5) == Window([1, 5, 7, 11])
true
```

See also [`wnunid`](@ref), [`wnvald!`](@ref).
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

Sort and merge the endpoint pairs in `window.data` to make a valid window, and return
`window`, as `wnvald_c` of CSPICE.

The pairs can be in any order and can overlap. `Window(endpoints)` calls this, so it is
necessary only after a direct change of `window.data`.

Throws `SpiceError(:unmatchendpts)` if the number of endpoints is odd, and
`SpiceError(:badendpoints)` if a pair has `left > right`.

See also [`Window`](@ref), [`wninsd!`](@ref).
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

Return the union of the windows `a` and `b`, as `wnunid_c` of CSPICE.

# Examples
```jldoctest
julia> a = Window([1, 3, 7, 11]); b = Window([2, 8]);

julia> wnunid(a, b) == Window([1, 11])
true

julia> wnintd(a, b) == Window([2, 3, 7, 8])
true

julia> wndifd(a, b) == Window([1, 2, 8, 11])
true

julia> wncomd(a, 0, 10) == Window([0, 1, 3, 7])
true
```

See also [`wnintd`](@ref), [`wndifd`](@ref), [`wncomd`](@ref), [`wninsd!`](@ref).
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

Return the intersection of the windows `a` and `b`, as `wnintd_c` of CSPICE.

Two intervals that share only an endpoint give the interval of that one point.

See also [`wnunid`](@ref), which has an example, and [`wndifd`](@ref).
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

Return the difference of the windows `a` and `b`, as `wndifd_c` of CSPICE.

The result holds closed intervals. An interval of `b` that cuts an interval of `a`
leaves the shared endpoints in the result: [1, 5] minus [2, 3] is [1, 2] and [3, 5]. An
interval of `a` that is inside an interval of `b`, endpoints included, is removed
whole. An interval of `b` of one point removes nothing from a longer interval of `a`.

See also [`wnunid`](@ref), which has an example, and [`wncomd`](@ref).
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

Return the complement of `window` with respect to the interval [`left`, `right`], as
`wncomd_c` of CSPICE.

The result holds the closed parts of the interval that `window` does not cover. The
window is the first argument; CSPICE has it after the interval.

Throws `SpiceError(:badendpoints)` if `left > right`.

See also [`wnunid`](@ref), which has an example, and [`wndifd`](@ref).
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

Return `true` if `point` is in an interval of `window`, endpoints included, as
`wnelmd_c` of CSPICE. The window is the first argument; CSPICE has it last.

See also [`wnincd`](@ref).
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

Return `true` if the interval [`left`, `right`] is inside one interval of `window`,
endpoints included, as `wnincd_c` of CSPICE. The window is the first argument; CSPICE
has it last.

See also [`wnelmd`](@ref).
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

Return the endpoints of interval number `n` of `window`, as `wnfetd_c` of CSPICE. `n`
starts at 1; in CSPICE it starts at 0.

Throws `SpiceError(:nointerval)` if `n` is not in `1:wncard(window)`.

# Examples
```jldoctest
julia> wnfetd(Window([1, 3, 7, 11]), 2)
(7.0, 11.0)
```

See also [`wncard`](@ref).
"""
function wnfetd(w::Window, n::Integer)
    (1 <= n <= wncard(w)) || throw(SpiceError(:nointerval,
        "wnfetd: interval $n is not in a window of $(wncard(w)) interval(s)"))
    return w.data[2n-1], w.data[2n]
end

"""
    wnsumd(window) -> (meas, avg, stddev, idxsml, idxlon)

Return a summary of `window`, as `wnsumd_c` of CSPICE.

- `meas`: the sum of the lengths of the intervals.
- `avg`, `stddev`: the mean and the standard deviation of the lengths.
- `idxsml`, `idxlon`: the index in `window` of the left endpoint of the shortest and of
  the longest interval; the first one if several are equal.

The indexes start at 1; in CSPICE they start at 0. An empty window gives 0 for all five
values.

# Examples
```jldoctest
julia> wnsumd(Window([1, 3, 7, 11]))
(6.0, 3.0, 1.0, 1, 3)
```

See also [`wncard`](@ref), [`wnfetd`](@ref).
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
