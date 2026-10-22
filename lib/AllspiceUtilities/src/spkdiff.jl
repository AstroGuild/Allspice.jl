# spkdiff: compare the trajectories of two bodies from two sets of kernels.

"""
    SpkDiff

The result of [`spkdiff`](@ref): geometric states of the first body from the first set of
kernels and of the second body from the second set, at the same epochs.

- `times`: the epochs, seconds past J2000 TDB.
- `state1`, `state2`: the states, position in km and velocity in km/s.
- `label1`, `label2`: descriptions of the two trajectories for the report.

[`basicstats`](@ref) and [`viewstats`](@ref) summarize the differences; `state1 .-
state2` gives the dump of `spkdiff -t dump`.
"""
struct SpkDiff
    times::Vector{Float64}
    state1::Vector{SVector{6,Float64}}
    state2::Vector{SVector{6,Float64}}
    label1::String
    label2::String
end

"""
    spkdiff(spk1, spk2; body, center, frame="J2000", body2=body, center2=center,
            frame2=frame, kernels=String[], kernels1=String[], kernels2=String[],
            start=nothing, stop=nothing, n=1000, step=nothing) -> SpkDiff

Compare the geometric state of `body` relative to `center` in `frame` from the SPK file
`spk1` with the state of `body2` relative to `center2` in `frame2` from `spk2`, as the
NAIF program `spkdiff`.

Each side has its own kernel context: `kernels`, then `kernels1` (or `kernels2`), then
the SPK file are loaded into it, so that the SPK file has priority. Bodies are names or
NAIF ID codes, frames are names. `start` and `stop` are seconds past J2000 TDB or time
strings (which need a leapseconds kernel in `kernels`); the default is the coverage of
`body` in `spk1` intersected with that of `body2` in `spk2`. The epochs are `n` points at
equal distance over that window or, with `step`, every `step` seconds.

Throws `ArgumentError` for an unknown body or an empty window, and the `SpiceError` of
`Allspice.spkezr` where a state cannot be computed.
"""
function spkdiff(spk1::AbstractString, spk2::AbstractString; body, center,
        frame::AbstractString="J2000", body2=body, center2=center,
        frame2::AbstractString=frame, kernels=String[], kernels1=String[],
        kernels2=String[], start=nothing, stop=nothing, n::Integer=1000, step=nothing)
    ctx1 = _context(kernels, kernels1, spk1)
    ctx2 = _context(kernels, kernels2, spk2)
    b1, c1 = _bodyid(ctx1, body), _bodyid(ctx1, center)
    b2, c2 = _bodyid(ctx2, body2), _bodyid(ctx2, center2)
    window = wnintd(spkcov(spk1, b1), spkcov(spk2, b2))
    if start !== nothing || stop !== nothing
        lo = start === nothing ? -Inf : _et(ctx1, start)
        hi = stop === nothing ? Inf : _et(ctx1, stop)
        lo <= hi || throw(ArgumentError("start is after stop"))
        window = wnintd(window, Window([lo, hi]))
    end
    wncard(window) > 0 || throw(ArgumentError("the two SPK files have no common coverage " *
        "of the bodies in the interval"))
    times = _sample(window, n, step)
    s1 = [spkezr(ctx1, b1, t, frame, "NONE", c1)[1] for t in times]
    s2 = [spkezr(ctx2, b2, t, frame2, "NONE", c2)[1] for t in times]
    label(ctx, b, c, f, spk) = "'$f'-referenced states of '$(_bodyname(ctx, b))' ($b) " *
        "relative to '$(_bodyname(ctx, c))' ($c) from SPK '$spk'"
    return SpkDiff(times, s1, s2, label(ctx1, b1, c1, frame, spk1),
        label(ctx2, b2, c2, frame2, spk2))
end

"""
    basicstats(d::SpkDiff) -> NamedTuple

Return the summary of `spkdiff -t basic`: the maximum and the average over the epochs of
the relative differences `relpos = |r1 - r2| / max(|r1|, |r2|)` and
`relvel = |v1 - v2| / max(|v1|, |v2|)`
and of the absolute differences `pos` (km) and `vel` (km/s), as the fields
`maxrelpos`, `avgrelpos`, `maxrelvel`, `avgrelvel`, `maxpos`, `avgpos`, `maxvel` and
`avgvel`.
"""
function basicstats(d::SpkDiff)
    @assert !isempty(d.times)
    dp = [norm(a[1:3] - b[1:3]) for (a, b) in zip(d.state1, d.state2)]
    dv = [norm(a[4:6] - b[4:6]) for (a, b) in zip(d.state1, d.state2)]
    rp = dp ./ [max(norm(a[1:3]), norm(b[1:3])) for (a, b) in zip(d.state1, d.state2)]
    rv = dv ./ [max(norm(a[4:6]), norm(b[4:6])) for (a, b) in zip(d.state1, d.state2)]
    avg(x) = sum(x) / length(x)
    return (maxrelpos=maximum(rp), avgrelpos=avg(rp), maxrelvel=maximum(rv),
        avgrelvel=avg(rv), maxpos=maximum(dp), avgpos=avg(dp), maxvel=maximum(dv),
        avgvel=avg(dv))
end

"""
    viewstats(d::SpkDiff) -> NamedTuple

Return the summary of `spkdiff -t stats`: the position differences `r1 - r2` in the view
frame of the first state, whose axes are down track (along `v1`), normal to the orbit
plane (along `r1 × v1`) and in the orbit plane (down track × normal).

Each of the fields `mean` (average components), `meanabs` (average absolute
components), `rms`, `atmaxrel` (at the epoch of the largest relative difference) and
`atmaxabs` (at the epoch of the largest absolute difference) is a named tuple
`(downtrack, inplane, normal, dt)` in km, with `dt` the down-track component divided by
`|v1|` in seconds. `epochmaxrel` and `epochmaxabs` are those two epochs.
"""
function viewstats(d::SpkDiff)
    @assert !isempty(d.times)
    comps = map(zip(d.state1, d.state2)) do (a, b)
        r, v = a[1:3], a[4:6]
        x = v / norm(v)
        z = cross(r, v)
        z /= norm(z)
        y = cross(x, z)
        δ = a[1:3] - b[1:3]
        SVector(dot(δ, x), dot(δ, y), dot(δ, z), dot(δ, x) / norm(v))
    end
    tup(c) = (downtrack=c[1], inplane=c[2], normal=c[3], dt=c[4])
    k = length(comps)
    dp = [norm(a[1:3] - b[1:3]) for (a, b) in zip(d.state1, d.state2)]
    rp = dp ./ [max(norm(a[1:3]), norm(b[1:3])) for (a, b) in zip(d.state1, d.state2)]
    irel, iabs = argmax(rp), argmax(dp)
    return (mean=tup(sum(comps) / k), meanabs=tup(sum(c -> abs.(c), comps) / k),
        rms=tup(sqrt.(sum(c -> c .^ 2, comps) / k)), atmaxrel=tup(comps[irel]),
        atmaxabs=tup(comps[iabs]), epochmaxrel=d.times[irel], epochmaxabs=d.times[iabs])
end

function Base.show(io::IO, ::MIME"text/plain", d::SpkDiff)
    s = basicstats(d)
    f(x) = @sprintf("%.13E", x)
    println(io, "#\n# Comparison of $(length(d.times)) ", d.label1, "\n#")
    println(io, "# with $(length(d.times)) ", d.label2, "\n#")
    println(io, "# within the time interval\n#")
    println(io, "#    from '", _cal(first(d.times)), " TDB' (", first(d.times), " TDB seconds)")
    println(io, "#    to   '", _cal(last(d.times)), " TDB' (", last(d.times), " TDB seconds)\n#\n")
    println(io, "Relative differences in state vectors:\n")
    println(io, " "^30, rpad("maximum", 24), "average\n")
    println(io, "  Position:             ", rpad(f(s.maxrelpos), 25), f(s.avgrelpos))
    println(io, "  Velocity:             ", rpad(f(s.maxrelvel), 25), f(s.avgrelvel))
    println(io, "\n\nAbsolute differences in state vectors:\n")
    println(io, " "^30, rpad("maximum", 24), "average\n")
    println(io, "  Position (km):        ", rpad(f(s.maxpos), 25), f(s.avgpos))
    print(io, "  Velocity (km/s):      ", rpad(f(s.maxvel), 25), f(s.avgvel))
end
