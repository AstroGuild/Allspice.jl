# Ephemeris tables: the position of one body relative to another on a fixed interval of
# time, as one flat array of Chebyshev coefficients on intervals of equal length. A
# query needs no kernel set, no segment search and no chain of bodies, so it is suited
# to the inner loop of an optimizer and to a GPU kernel (the table is an `isbits` value
# when its array is one, see `ext/AllspiceAdaptExt.jl`).

"""
    EphemerisTable{T,A<:AbstractVector{T}}
    EphemerisTable{T}(tab) -> EphemerisTable{T,Vector{T}}

Position of one body relative to another on a fixed span of time, as piecewise
Chebyshev polynomials that need no kernel to evaluate. [`ephemeris_table`](@ref)
makes one.

A query does no segment search and follows no chain of bodies, so it suits the inner
loop of an optimizer. With Adapt.jl loaded, `adapt` moves the coefficients to another
array type; the table is then an `isbits` value on a GPU.

| Field              | Meaning                                                              |
|:-------------------|:---------------------------------------------------------------------|
| `target`, `observer`, `frame` | NAIF IDs of the bodies and ID of the frame                |
| `t0`, `dt`, `nint` | the table covers `t0 <= et <= t0 + nint * dt`, seconds past J2000 TDB |
| `ncoef`            | Chebyshev coefficients for each component and interval               |
| `coef`             | the coefficients, km                                                 |
| `maxerr`           | largest difference from the source that `ephemeris_table` measured, km |

Interval `i` (from 0) covers `t0 + i dt` to `t0 + (i + 1) dt`. It holds `ncoef`
coefficients for each of x, y and z, in this order, from index `3 ncoef i + 1` of
`coef`: the layout of a record of SPK type 2 without its first two numbers.

`EphemerisTable{T}(tab)` returns a copy with coefficients of type `T` on the CPU. With
`Float32` the position keeps about seven digits (15 km at one astronomical unit);
`maxerr` is not measured again.

Throws `SpiceError(:badarg)` if the constructor gets `dt <= 0`, `nint < 1`,
`ncoef < 1` or `length(coef) != 3 * ncoef * nint`.

# Examples
```jldoctest
julia> tab = EphemerisTable([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], 0.0, 10.0, 1, 2, 301, 399, 1, 0.0);

julia> table_span(tab)
(0.0, 10.0)

julia> table_state(tab, 7.5)    # x = 1 + 2s, y = 3 + 4s, z = 5 + 6s at s = 0.5
6-element StaticArraysCore.SVector{6, Float64} with indices SOneTo(6):
 2.0
 5.0
 8.0
 0.4
 0.8
 1.2
```

See also [`table_position`](@ref), [`table_state`](@ref), [`table_jet`](@ref),
[`table_to_spk`](@ref), [`RotationTable`](@ref).
"""
struct EphemerisTable{T<:AbstractFloat,A<:AbstractVector{T}}
    coef::A
    t0::Float64
    dt::Float64
    nint::Int
    ncoef::Int
    target::Int32
    observer::Int32
    frame::Int32
    maxerr::Float64
    function EphemerisTable(coef::A, t0::Real, dt::Real, nint::Integer, ncoef::Integer,
            target::Integer, observer::Integer, frame::Integer, maxerr::Real) where {T<:AbstractFloat,A<:AbstractVector{T}}
        (dt > 0 && nint >= 1 && ncoef >= 1 && length(coef) == 3 * ncoef * nint) ||
            throw(SpiceError(:badarg, "EphemerisTable: the sizes do not agree"))
        return new{T,A}(coef, t0, dt, nint, ncoef, target, observer, frame, maxerr)
    end
end

EphemerisTable{T}(tab::EphemerisTable) where {T<:AbstractFloat} =
    EphemerisTable(Vector{T}(tab.coef), tab.t0, tab.dt, tab.nint, tab.ncoef, tab.target,
        tab.observer, tab.frame, tab.maxerr)

"""
    table_span(tab) -> (first, last)

Return the first and the last epoch of the table `tab`, in seconds past J2000 TDB.

`tab` is an [`EphemerisTable`](@ref) or a [`RotationTable`](@ref). A query is valid
for `first <= et <= last`.
"""
@inline table_span(tab::EphemerisTable) = (tab.t0, tab.t0 + tab.nint * tab.dt)

Base.show(io::IO, tab::EphemerisTable{T}) where {T} = print(io, "EphemerisTable{", T, "}(",
    tab.target, " relative to ", tab.observer, ", frame ", tab.frame, ", ET ", tab.t0,
    " to ", table_span(tab)[2], ", ", tab.nint, " intervals of ", tab.dt, " s, ",
    tab.ncoef, " coefficients, largest error ", tab.maxerr, " km)")

@noinline _throw_table_range() =
    throw(SpiceError(:nocoverage, "the epoch is outside the ephemeris table"))

@inline _tab_arg(::Type{T}, s::AbstractFloat) where {T} = T(s)
@inline _tab_arg(::Type, s::Real) = s

"""
    table_jet(tab, et, Val(N)) -> NTuple{N+1,SVector{3}}

Return the position from the ephemeris table `tab` at `et`, in km, and its time
derivatives to order `N`.

Element `k + 1` is the derivative of order `k`, in km/s^k; every order is a
derivative of the position polynomial. `et` is in seconds past J2000 TDB and is a
plain number or a dual number. A plain epoch is used as `Float64` also when the
coefficients are `Float32`. The call reads no kernel and does not allocate.

Throws `SpiceError(:nocoverage)` if `et` is outside [`table_span`](@ref)`(tab)`. The
check is a bounds check. Under `@inbounds`, as in a GPU kernel, it is removed, and an
epoch outside the table then gives the extrapolation of the nearest interval.

See also [`table_position`](@ref), [`table_state`](@ref).
"""
Base.@propagate_inbounds function table_jet(tab::EphemerisTable{T}, et::Real, v::Val{N}) where {T,N}
    t = primal(et)
    @boundscheck (tab.t0 <= t <= tab.t0 + tab.nint * tab.dt) || _throw_table_range()
    i = clamp(unsafe_trunc(Int, (t - tab.t0) / tab.dt), 0, tab.nint - 1)
    radius = tab.dt / 2
    s = (et - (tab.t0 + (2i + 1) * radius)) / radius
    return _cheb3(tab.coef, 3 * tab.ncoef * i + 1, tab.ncoef, _tab_arg(T, s),
        _tab_arg(T, radius), v)
end

"""
    table_position(tab, et) -> SVector{3}

Return the position from the ephemeris table `tab` at `et`, in km.

`et` is in seconds past J2000 TDB. The bodies and the frame are `tab.target`,
`tab.observer` and `tab.frame`. The rules for the epoch, for `@inbounds` and for the
number types are those of [`table_jet`](@ref).

Throws `SpiceError(:nocoverage)` if `et` is outside the table.
"""
Base.@propagate_inbounds table_position(tab::EphemerisTable, et::Real) = table_jet(tab, et, Val(0))[1]

"""
    table_state(tab, et) -> SVector{6}

Return the position (km) and the velocity (km/s) from the ephemeris table `tab` at
`et`.

The velocity is the derivative of the position polynomial. For a table with an
aberration correction it is therefore order 1 of [`spkder`](@ref), not the velocity
of [`spkezr`](@ref). The rules for the epoch, for `@inbounds` and for the number
types are those of [`table_jet`](@ref).

Throws `SpiceError(:nocoverage)` if `et` is outside the table.
"""
Base.@propagate_inbounds function table_state(tab::EphemerisTable, et::Real)
    r, v = table_jet(tab, et, Val(1))
    return vcat(r, v)
end

# ---------------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------------

"Segments from `body` to the root of its chain at `t`, or `nothing` if the chain is too long."
function _tab_up(ks::KernelSet, body::Int32, t::Float64)
    list = Segment[]
    while true     # invariant: `list` is the chain from the first body to `body`
        segs = get(ks.spk, body, nothing)
        segs === nothing && return list
        i = find_segment(segs, t)
        i == 0 && return list
        length(list) == MAXCHAIN && return nothing
        push!(list, segs[i])
        body = segs[i].center
    end
end

"""
    _tab_chain(ks, targ, obs, t) -> Union{Nothing,Vector{Segment}}

The segments that `spk_jet` adds and subtracts for `targ` relative to `obs` at `t`,
or `nothing` if the two bodies are not connected at `t`.
"""
function _tab_chain(ks::KernelSet, targ::Int32, obs::Int32, t::Float64)
    a, b = _tab_up(ks, targ, t), _tab_up(ks, obs, t)
    (a === nothing || b === nothing) && return nothing
    root(list, body) = isempty(list) ? body : last(list).center
    root(a, targ) == root(b, obs) || return nothing
    while !isempty(a) && !isempty(b) && last(a) == last(b)   # invariant: the tails removed are equal
        pop!(a)
        pop!(b)
    end
    return append!(a, b)
end

"""
    _tab_native_grid(ks, targ, ref, obs, t0, t1) -> Union{Nothing,Tuple{Float64,Float64,Int,Int}}

`(start, dt, nint, ncoef)` of the grid on which the position of `targ` relative to
`obs` in frame `ref` is one polynomial of degree `ncoef - 1` for each interval, or
`nothing` if there is no such grid. A grid exists when `ref` is inertial, every
segment of the chain is SPK type 2 in an inertial frame, the record boundaries of
all segments fall on the boundaries of the shortest record, and the same segments
serve the whole of `start` to `start + nint * dt`, which contains `t0` to `t1`.
"""
function _tab_native_grid(ks::KernelSet, targ::Int32, ref::Int32, obs::Int32, t0::Float64, t1::Float64)
    _is_inertial(ks, ref) || return nothing
    chain = _tab_chain(ks, targ, obs, (t0 + t1) / 2)
    (chain === nothing || isempty(chain)) && return nothing
    all(s -> s.dtype == Int32(2) && _is_inertial(ks, s.frame), chain) || return nothing
    short = chain[argmin(map(s -> s.intlen, chain))]
    dt, origin = short.intlen, short.init
    whole(x) = abs(x - round(x)) <= 1e-9
    all(s -> whole((s.init - origin) / dt) && whole(s.intlen / dt), chain) || return nothing
    i0 = floor(Int, (t0 - origin) / dt)
    i1 = max(i0 + 1, ceil(Int, (t1 - origin) / dt))
    start, nint = origin + i0 * dt, i1 - i0
    slack = 1e-9 * dt
    all(s -> s.t0 <= start + slack && s.t1 >= start + nint * dt - slack, chain) || return nothing
    for i in 0:nint-1      # invariant: the chain of every interval before `i` is `chain`
        _tab_chain(ks, targ, obs, start + (i + 0.5) * dt) == chain || return nothing
    end
    ncoef = maximum(s -> (Int(s.rsize) - 2) ÷ 3, chain)
    return start, dt, nint, ncoef
end

"""
    _tab_fit(f, start, dt, nint, ncoef) -> (coef::Vector{Float64}, maxerr, scale)

Chebyshev coefficients of `f(et)::SVector{3}` on `nint` intervals of length `dt` from
`start`: interpolation at the `ncoef` zeros of the Chebyshev polynomial of degree
`ncoef`. An epoch near 1e9 s has a resolution of 1e-7 s, so a zero is not an exact
epoch; the function solves the interpolation conditions at the normalized times of
the epochs that `f` received, which keeps a polynomial `f` exact to rounding.
`maxerr` is the largest norm of the difference between `f` and the interpolant at the
`ncoef - 1` points between the zeros (the points where the error of the interpolation
is largest); `scale` is the largest norm of `f` at the zeros.
Requires `dt > 0`, `nint >= 1`, `ncoef >= 2`.
"""
_tab_fit(f::F, start::Float64, dt::Float64, nint::Int, ncoef::Int) where {F} =
    _tab_fit(f, start, dt, nint, ncoef, Val(1))

"The value at `s` of the `C` vector series of one record of `ncoef` coefficients from index `a`, as one vector."
@inline function _tab_value(coef, a::Int, n::Int, s, radius, ::Val{C}) where {C}
    cols = ntuple(c -> _cheb3(coef, a + 3 * n * (c - 1), n, s, radius, Val(0))[1], Val(C))
    return vcat(cols...)
end

"""
    _tab_fit(f, start, dt, nint, ncoef, Val(C)) -> (coef, maxerr, scale)

The fit above for a function `f(et)` of `3C` components (for `C = 3`, the columns of a
matrix one after the other): a record holds `ncoef` coefficients for each component,
in the order of the components, and the two norms are norms of all `3C` components.
"""
function _tab_fit(f::F, start::Float64, dt::Float64, nint::Int, ncoef::Int, vc::Val{C}) where {F,C}
    @assert dt > 0 && nint >= 1 && ncoef >= 2
    n = ncoef
    x = [cospi((j - 0.5) / n) for j in 1:n]
    xc = [cospi(j / n) for j in 1:n-1]
    w = 3 * C * n                       # numbers in one record
    coef = zeros(w * nint)
    vals = Matrix{Float64}(undef, n, 3 * C)
    basis = Matrix{Float64}(undef, n, n)
    radius = dt / 2
    scale = 0.0
    for i in 0:nint-1      # invariant: the coefficients of the intervals before `i` are set
        mid = start + (2i + 1) * radius
        for j in 1:n
            et = mid + radius * x[j]
            s = (et - mid) / radius
            r = f(et)
            vals[j, :] .= r
            scale = max(scale, norm(r))
            basis[j, 1] = 1.0
            basis[j, 2] = s
            for k in 3:n       # invariant: basis[j, m] is the polynomial of degree m - 1 at s
                basis[j, k] = 2 * s * basis[j, k-1] - basis[j, k-2]
            end
        end
        coef[w*i+1:w*(i+1)] .= vec(basis \ vals)
    end
    maxerr = 0.0
    for i in 0:nint-1, sc in xc
        mid = start + (2i + 1) * radius
        et = mid + radius * sc
        fit = _tab_value(coef, w * i + 1, n, (et - mid) / radius, radius, vc)
        maxerr = max(maxerr, norm(fit - f(et)))
    end
    @assert length(coef) == w * nint
    return coef, maxerr, scale
end

"""
    _tab_refine(f, a, b, tol, ncoef, maxintervals, Val(C), who, unit) -> (coef, dt, nint, maxerr)

Fit `f` on `a` to `b` with `_tab_fit` and double the number of intervals until the
measured error is at most `tol`. Without `tol` (`nothing`) the goal is `1e-12` times
the largest norm of `f`, and the loop stops also when the error is small and one more
doubling did not divide it by four: the source then limits the accuracy. Throws
`SpiceError(:noconvergence)` if `maxintervals` intervals do not reach the goal.
"""
function _tab_refine(f::F, a::Float64, b::Float64, tol::Union{Nothing,Real}, ncoef::Int,
        maxintervals::Integer, vc::Val{C}, who::AbstractString, unit::AbstractString) where {F,C}
    nint, previous = 1, Inf
    while true     # invariant: no number of intervals below `nint` was accepted
        dt = (b - a) / nint
        coef, maxerr, scale = _tab_fit(f, a, dt, nint, ncoef, vc)
        goal = tol === nothing ? 1e-12 * scale : Float64(tol)
        floor_reached = tol === nothing && maxerr <= 1e-6 * scale && maxerr > previous / 4
        (maxerr <= goal || floor_reached) && return coef, dt, nint, maxerr
        2 * nint > maxintervals && throw(SpiceError(:noconvergence,
            "$who: $nint intervals give an error of $maxerr$unit, above the tolerance of $goal$unit"))
        nint, previous = 2 * nint, maxerr
    end
end

"""
    ephemeris_table([ctx,] targ, t0, t1, ref, abcorr, obs;
                    tol=nothing, ncoef=16, maxintervals=65536) -> EphemerisTable

Return a table of the position of `targ` relative to `obs` in frame `ref` for
`t0 <= et <= t1`, which [`table_position`](@ref), [`table_state`](@ref) and
[`table_jet`](@ref) evaluate without the kernels. This call has no CSPICE counterpart.

`t0` and `t1` are in seconds past J2000 TDB. `targ`, `ref`, `abcorr` and `obs` are as
for [`spkpos`](@ref). The table is a snapshot: a later change of the kernels does not
change it. The result is an `EphemerisTable{Float64,Vector{Float64}}`, and its field
`maxerr` is the measured error in km.

Two methods give the coefficients:

- **Records of the kernels.** If `abcorr` is `"NONE"`, `ref` is inertial, and every
  segment between the two bodies is SPK type 2 in an inertial frame with record
  boundaries on one common grid, the table takes the record length and the degree of
  the kernels. It then agrees with `spkpos` to rounding, it can start before `t0` and
  end after `t1`, and `tol`, `ncoef` and `maxintervals` are not used.
- **Interpolation.** In all other cases the call interpolates on `t0` to `t1` with
  `ncoef` coefficients for each interval, and it doubles the number of intervals
  until the measured error is at most `tol` (km). Without `tol` the goal is `1e-12`
  times the largest distance in the table. The call then stops also when one more
  doubling does not divide a small error by four: the source limits the accuracy
  (rounding of a rotation angle of many turns, or the small steps of an
  aberration-corrected position at the record boundaries of the kernels).

Throws `SpiceError(:badarg)` if the epochs are not finite with `t0 < t1`, if
`ncoef < 2`, `tol <= 0` or `maxintervals < 1`. Throws `SpiceError(:noconvergence)` if
`maxintervals` intervals do not reach the goal, and the errors of [`spkpos`](@ref) if
the kernels do not cover `t0` to `t1`.

# Examples
```julia
furnsh("naif0012.tls"); furnsh("de440s.bsp")
t0, t1 = str2et("2025-01-01"), str2et("2026-01-01")
tab = ephemeris_table("MARS BARYCENTER", t0, t1, "ECLIPJ2000", "NONE", "SUN")
tab.maxerr                              # km
r = table_position(tab, t0 + 86400.0)   # km, no kernel access
x = table_state(tab, t0 + 86400.0)      # km, km/s
tab32 = EphemerisTable{Float32}(tab)    # smaller, about seven digits
```

See also [`EphemerisTable`](@ref), [`rotation_table`](@ref), [`table_to_spk`](@ref).
"""
function ephemeris_table(ks::KernelSet, targ::NameOrID, t0::Real, t1::Real, ref::NameOrID,
        abcorr::Union{AbstractString,Aberration}, obs::NameOrID;
        tol::Union{Nothing,Real}=nothing, ncoef::Integer=16, maxintervals::Integer=65536)
    a, b = Float64(t0), Float64(t1)
    (isfinite(a) && isfinite(b) && a < b) ||
        throw(SpiceError(:badarg, "ephemeris_table: the epochs must be finite, with t0 < t1"))
    ncoef >= 2 || throw(SpiceError(:badarg, "ephemeris_table: ncoef must be at least 2"))
    (tol === nothing || tol > 0) || throw(SpiceError(:badarg, "ephemeris_table: tol must be positive"))
    maxintervals >= 1 || throw(SpiceError(:badarg, "ephemeris_table: maxintervals must be at least 1"))
    t, f, o = body_code(ks, targ), frame_id(ks, ref), body_code(ks, obs)
    corr = Aberration(abcorr)
    source = corr.geometric ? (et -> spk_jet(ks, t, et, f, o, Val(0))[1]) :
             (et -> abcorr_position(ks, t, et, f, corr, o)[1])
    source(a), source(b)       # data at both ends, or the error of the SPK layer
    grid = corr.geometric ? _tab_native_grid(ks, t, f, o, a, b) : nothing
    if grid !== nothing
        start, dt, nint, n = grid
        coef, maxerr, _ = _tab_fit(source, start, dt, nint, n)
        tab = EphemerisTable(coef, start, dt, nint, n, t, o, f, maxerr)
    else
        coef, dt, nint, maxerr = _tab_refine(source, a, b, tol, Int(ncoef), maxintervals, Val(1),
            "ephemeris_table", " km")
        tab = EphemerisTable(coef, a, dt, nint, ncoef, t, o, f, maxerr)
    end
    lo, hi = table_span(tab)
    @assert lo <= a && b <= hi + 1e-9 * tab.dt "the table does not cover the request"
    return tab
end
ephemeris_table(ctx::SpiceContext, targ::NameOrID, t0::Real, t1::Real, ref::NameOrID,
    abcorr::Union{AbstractString,Aberration}, obs::NameOrID; kw...) =
    ephemeris_table(kernelset(ctx), targ, t0, t1, ref, abcorr, obs; kw...)
ephemeris_table(targ::NameOrID, t0::Real, t1::Real, ref::NameOrID,
    abcorr::Union{AbstractString,Aberration}, obs::NameOrID; kw...) =
    ephemeris_table(default_context(), targ, t0, t1, ref, abcorr, obs; kw...)

# ---------------------------------------------------------------------------------
# Export of a table as an SPK segment
# ---------------------------------------------------------------------------------

"""
    table_to_spk(handle, tab, segid; first=table_span(tab)[1], last=table_span(tab)[2],
                 ctx=default_context()) -> Nothing

Add the ephemeris table `tab` to the SPK file `handle` as one segment of type 2.

The records of the table are the records of the segment. [`spkpos`](@ref) on the new
file therefore gives the position of the table to rounding, and the velocity is its
derivative. `handle` comes from [`spkopn`](@ref) or [`spkopa`](@ref); `segid` is the
segment identifier; `first` and `last` are the bounds of the segment in seconds past
J2000 TDB. The segment is for `tab.target` relative to `tab.observer` in the frame
`tab.frame`. The coefficients are written as `Float64`; a table on a GPU is copied to
the CPU.

A segment has no aberration correction. Write only a table that was made with
`abcorr = "NONE"`; the call cannot check this.

Throws the errors of [`spkw02`](@ref). Among its conditions: `tab.ncoef <= 28`,
`tab.target != 0`, `tab.target != tab.observer`, a known frame, and `first` and `last`
inside the table.

# Examples
```julia
tab = ephemeris_table("MOON", t0, t1, "J2000", "NONE", "EARTH")
h = spkopn("moon.bsp", "moon.bsp", 0)
table_to_spk(h, tab, "Moon from table")
spkcls(h)
```
"""
function table_to_spk(h::DAFHandle, tab::EphemerisTable, segid::AbstractString;
        first::Real=table_span(tab)[1], last::Real=table_span(tab)[2],
        ctx::SegContext=default_context())
    cdata = reshape(Vector{Float64}(Array(tab.coef)), 3 * tab.ncoef, tab.nint)
    spkw02(h, tab.target, tab.observer, Int(tab.frame), first, last, segid, tab.dt, cdata, tab.t0;
        ctx=ctx)
    return nothing
end

# ---------------------------------------------------------------------------------
# Rotation tables
# ---------------------------------------------------------------------------------

"""
    RotationTable{T,A<:AbstractVector{T}}
    RotationTable{T}(tab) -> RotationTable{T,Vector{T}}

Rotation matrix from one frame to another on a fixed span of time, as piecewise
Chebyshev polynomials that need no kernel to evaluate: the counterpart of
[`EphemerisTable`](@ref) for [`pxform`](@ref). [`rotation_table`](@ref) makes one.

| Field              | Meaning                                                              |
|:-------------------|:---------------------------------------------------------------------|
| `from`, `to`       | frame IDs; the matrix rotates a vector from `from` to `to`           |
| `t0`, `dt`, `nint` | the table covers `t0 <= et <= t0 + nint * dt`, seconds past J2000 TDB |
| `ncoef`            | Chebyshev coefficients for each matrix element and interval          |
| `coef`             | the coefficients                                                     |
| `maxerr`           | largest Frobenius norm of the difference from the source that `rotation_table` measured |

Interval `i` (from 0) holds `ncoef` coefficients for each of the nine elements of the
matrix, column after column, from index `9 ncoef i + 1` of `coef`. A matrix of the
table is orthogonal to the accuracy `maxerr`; it is not made orthogonal again.

`RotationTable{T}(tab)` returns a copy with coefficients of type `T` on the CPU;
`maxerr` is not measured again.

Throws `SpiceError(:badarg)` if the constructor gets `dt <= 0`, `nint < 1`,
`ncoef < 1` or `length(coef) != 9 * ncoef * nint`.

See also [`table_rotation`](@ref), [`table_rotation_jet`](@ref), [`table_span`](@ref).
"""
struct RotationTable{T<:AbstractFloat,A<:AbstractVector{T}}
    coef::A
    t0::Float64
    dt::Float64
    nint::Int
    ncoef::Int
    from::Int32
    to::Int32
    maxerr::Float64
    function RotationTable(coef::A, t0::Real, dt::Real, nint::Integer, ncoef::Integer,
            from::Integer, to::Integer, maxerr::Real) where {T<:AbstractFloat,A<:AbstractVector{T}}
        (dt > 0 && nint >= 1 && ncoef >= 1 && length(coef) == 9 * ncoef * nint) ||
            throw(SpiceError(:badarg, "RotationTable: the sizes do not agree"))
        return new{T,A}(coef, t0, dt, nint, ncoef, from, to, maxerr)
    end
end

RotationTable{T}(tab::RotationTable) where {T<:AbstractFloat} =
    RotationTable(Vector{T}(tab.coef), tab.t0, tab.dt, tab.nint, tab.ncoef, tab.from, tab.to,
        tab.maxerr)

@inline table_span(tab::RotationTable) = (tab.t0, tab.t0 + tab.nint * tab.dt)

Base.show(io::IO, tab::RotationTable{T}) where {T} = print(io, "RotationTable{", T, "}(frame ",
    tab.from, " to frame ", tab.to, ", ET ", tab.t0, " to ", table_span(tab)[2], ", ", tab.nint,
    " intervals of ", tab.dt, " s, ", tab.ncoef, " coefficients, largest error ", tab.maxerr, ")")

"""
    table_rotation_jet(tab, et, Val(N)) -> NTuple{N+1,SMatrix{3,3}}

Return the rotation matrix from the rotation table `tab` at `et` and its time
derivatives to order `N`.

Element `k + 1` is the derivative of order `k`, in 1/s^k; every order is a derivative
of the polynomials of the table. Order 1 thus follows order 1 of [`pxder`](@ref). It
is the derivative block of [`sxform`](@ref) for the frames where the two agree, which
excludes CK and dynamic frames. `et` is in seconds past J2000 TDB. The rules for the
epoch, for `@inbounds` and for the number types are those of [`table_jet`](@ref). The
call reads no kernel and does not allocate.

Throws `SpiceError(:nocoverage)` if `et` is outside [`table_span`](@ref)`(tab)`.
"""
Base.@propagate_inbounds function table_rotation_jet(tab::RotationTable{T}, et::Real, v::Val{N}) where {T,N}
    t = primal(et)
    @boundscheck (tab.t0 <= t <= tab.t0 + tab.nint * tab.dt) || _throw_table_range()
    i = clamp(unsafe_trunc(Int, (t - tab.t0) / tab.dt), 0, tab.nint - 1)
    radius = tab.dt / 2
    s = _tab_arg(T, (et - (tab.t0 + (2i + 1) * radius)) / radius)
    r = _tab_arg(T, radius)
    a = 9 * tab.ncoef * i + 1
    c1 = _cheb3(tab.coef, a, tab.ncoef, s, r, v)
    c2 = _cheb3(tab.coef, a + 3 * tab.ncoef, tab.ncoef, s, r, v)
    c3 = _cheb3(tab.coef, a + 6 * tab.ncoef, tab.ncoef, s, r, v)
    return ntuple(k -> hcat(c1[k], c2[k], c3[k]), Val(N + 1))
end

"""
    table_rotation(tab, et) -> SMatrix{3,3}

Return the rotation matrix from frame `tab.from` to frame `tab.to` at `et`, from the
rotation table `tab`.

`et` is in seconds past J2000 TDB. The rules for the epoch, for `@inbounds` and for
the number types are those of [`table_jet`](@ref).

Throws `SpiceError(:nocoverage)` if `et` is outside the table.
"""
Base.@propagate_inbounds table_rotation(tab::RotationTable, et::Real) =
    table_rotation_jet(tab, et, Val(0))[1]

"""
    rotation_table([ctx,] from, to, t0, t1; tol=nothing, ncoef=16, maxintervals=65536)
        -> RotationTable

Return a table of the rotation from frame `from` to frame `to` for `t0 <= et <= t1`,
which [`table_rotation`](@ref) and [`table_rotation_jet`](@ref) evaluate without the
kernels. This call has no CSPICE counterpart.

`t0` and `t1` are in seconds past J2000 TDB; the frames are as for [`pxform`](@ref).
The result is a `RotationTable{Float64,Vector{Float64}}` with
`table_span(result) == (t0, t1)` to rounding, and its field `maxerr` is the measured
error.

The call interpolates the nine elements of the matrix with `ncoef` coefficients for
each interval, and it doubles the number of intervals until the measured error
(Frobenius norm, no unit) is at most `tol`. Without `tol` the goal is `1.7e-12`. The
call then stops also when one more doubling does not divide a small error by four:
the source limits the accuracy (a rotation angle of many turns has an error of
rounding near 1e-11).

Throws `SpiceError(:badarg)` if the epochs are not finite with `t0 < t1`, if
`ncoef < 2`, `tol <= 0` or `maxintervals < 1`. Throws `SpiceError(:noconvergence)` if
`maxintervals` intervals do not reach the goal, and the errors of [`pxform`](@ref) if
the kernels do not cover `t0` to `t1`.

# Examples
```julia
furnsh("naif0012.tls"); furnsh("pck00011.tpc")
t0, t1 = str2et("2025-01-01"), str2et("2025-02-01")
tab = rotation_table("J2000", "IAU_MARS", t0, t1)
R = table_rotation(tab, t0 + 3600.0)
```

See also [`RotationTable`](@ref), [`ephemeris_table`](@ref).
"""
function rotation_table(ks::KernelSet, from::NameOrID, to::NameOrID, t0::Real, t1::Real;
        tol::Union{Nothing,Real}=nothing, ncoef::Integer=16, maxintervals::Integer=65536)
    a, b = Float64(t0), Float64(t1)
    (isfinite(a) && isfinite(b) && a < b) ||
        throw(SpiceError(:badarg, "rotation_table: the epochs must be finite, with t0 < t1"))
    ncoef >= 2 || throw(SpiceError(:badarg, "rotation_table: ncoef must be at least 2"))
    (tol === nothing || tol > 0) || throw(SpiceError(:badarg, "rotation_table: tol must be positive"))
    maxintervals >= 1 || throw(SpiceError(:badarg, "rotation_table: maxintervals must be at least 1"))
    f, t = frame_id(ks, from), frame_id(ks, to)
    source = et -> SVector{9,Float64}(frame_rotation(ks, f, t, et))
    source(a), source(b)       # data at both ends, or the error of the frame layer
    coef, dt, nint, maxerr = _tab_refine(source, a, b, tol, Int(ncoef), maxintervals, Val(3),
        "rotation_table", "")
    tab = RotationTable(coef, a, dt, nint, ncoef, f, t, maxerr)
    @assert table_span(tab)[1] == a && abs(table_span(tab)[2] - b) <= 1e-9 * (b - a)
    return tab
end
rotation_table(ctx::SpiceContext, from::NameOrID, to::NameOrID, t0::Real, t1::Real; kw...) =
    rotation_table(kernelset(ctx), from, to, t0, t1; kw...)
rotation_table(from::NameOrID, to::NameOrID, t0::Real, t1::Real; kw...) =
    rotation_table(default_context(), from, to, t0, t1; kw...)
