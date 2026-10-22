# frmdiff: compare the orientation of two frames from two sets of kernels.

"""
    FrmDiff

The result of [`frmdiff`](@ref): the rotations from the first `from` frame to the first
`to` frame, and from the second `from` frame to the second `to` frame, at the same
epochs.

- `times`: the epochs, seconds past J2000 TDB.
- `rot1`, `rot2`: the rotation matrices, as `Allspice.pxform` gives them.
- `av1`, `av2`: the angular velocities of the `from` frames relative to the `to` frames,
  in rad/s and in the `from` frames, or empty if they were not asked for.
- `avdiff`: `av1 - av2` in the first `from` frame, or rotated into the first `to` frame
  (`avframe=:to`), as `frmdiff -t dumpaa` gives it.
- `angle`: the angle of the rotation `rot1 * rot2'` that takes the second orientation
  to the first, in radians; `axis`: its axis.
"""
struct FrmDiff
    times::Vector{Float64}
    rot1::Vector{SMatrix{3,3,Float64,9}}
    rot2::Vector{SMatrix{3,3,Float64,9}}
    av1::Vector{SVector{3,Float64}}
    av2::Vector{SVector{3,Float64}}
    avdiff::Vector{SVector{3,Float64}}
    angle::Vector{Float64}
    axis::Vector{SVector{3,Float64}}
    label1::String
    label2::String
end

"""
    frmdiff([file1, [file2]]; from, to, from2=from, to2=to, kernels=String[],
            kernels1=String[], kernels2=String[], start=nothing, stop=nothing,
            n=1000, step=nothing, angvel=false, avframe=:from) -> FrmDiff

Compare the rotation from frame `from` to frame `to` computed with the first set of
kernels with the rotation from `from2` to `to2` computed with the second set, as the
NAIF program `frmdiff`.

Each set has its own kernel context: `kernels`, then `kernels1` (or `kernels2`), then
`file1` (or `file2`, which is `file1` if omitted) are loaded into it. Frames are names or
frame IDs. `start` and `stop` are seconds past J2000 TDB or time strings. Without them,
the window is the coverage of the frame `from` in `file1` if `file1` is a binary PCK or
a CK file that holds it. The epochs are `n` points at equal distance over the window or,
with `step`, every `step` seconds.

With `angvel=true` the angular velocities are computed too, and their difference is
expressed in the first `from` frame (`avframe=:from`) or the first `to` frame
(`avframe=:to`).

The NAIF program leaves out the last epoch before `stop` when there are more than two;
the result here keeps it.

Throws `ArgumentError` if there is no window, and the `SpiceError` of `Allspice.pxform`
or `Allspice.sxform` where a rotation cannot be computed.
"""
function frmdiff(files::AbstractString...; from, to, from2=from, to2=to, kernels=String[],
        kernels1=String[], kernels2=String[], start=nothing, stop=nothing,
        n::Integer=1000, step=nothing, angvel::Bool=false, avframe::Symbol=:from)
    length(files) <= 2 || throw(ArgumentError("frmdiff takes at most two kernel files"))
    avframe in (:from, :to) || throw(ArgumentError("avframe is :from or :to, got $avframe"))
    file1 = isempty(files) ? String[] : files[1]
    file2 = length(files) == 2 ? files[2] : file1
    ctx1 = _context(kernels, kernels1, file1)
    ctx2 = _context(kernels, kernels2, file2)
    window = if start !== nothing && stop !== nothing
        lo, hi = _et(ctx1, start), _et(ctx1, stop)
        lo <= hi || throw(ArgumentError("start is after stop"))
        Window([lo, hi])
    else
        w = isempty(files) ? nothing : _framecoverage(ctx1, files[1], from)
        w === nothing && throw(ArgumentError("give start and stop: the coverage of the " *
            "frame cannot be read from the first file"))
        start === nothing ? w : wnintd(w, Window([_et(ctx1, start), Inf]))
    end
    wncard(window) > 0 || throw(ArgumentError("the window is empty"))
    times = _sample(window, n, step)
    nt = length(times)
    rot1 = Vector{SMatrix{3,3,Float64,9}}(undef, nt)
    rot2 = similar(rot1)
    av1 = SVector{3,Float64}[]
    av2 = SVector{3,Float64}[]
    avdiff = SVector{3,Float64}[]
    for (k, t) in enumerate(times)
        if angvel
            r1, w1 = xf2rav(sxform(ctx1, from, to, t))
            r2, w2 = xf2rav(sxform(ctx2, from2, to2, t))
            rot1[k], rot2[k] = r1, r2
            push!(av1, w1)
            push!(av2, w2)
            push!(avdiff, avframe === :to ? r1 * (w1 - w2) : w1 - w2)
        else
            rot1[k] = pxform(ctx1, from, to, t)
            rot2[k] = pxform(ctx2, from2, to2, t)
        end
    end
    axis = Vector{SVector{3,Float64}}(undef, nt)
    angle = Vector{Float64}(undef, nt)
    for k in 1:nt
        axis[k], angle[k] = raxisa(rot1[k] * rot2[k]')
    end
    label(ctx, f, t, file) = "rotations from '$f' to '$t' computed with $(file)"
    return FrmDiff(times, rot1, rot2, av1, av2, avdiff, angle, axis,
        label(ctx1, from, to, isempty(files) ? "the kernels" : files[1]),
        label(ctx2, from2, to2, isempty(files) ? "the kernels" : file2))
end

"Coverage of the frame `from` in `file` if it is a binary PCK or CK file, else `nothing`."
function _framecoverage(ctx, file::AbstractString, from)
    info = frinfo(ctx, _frameid(ctx, from))
    info === nothing && return nothing
    _, class, classid = info
    arch, kind = getfat(file)
    arch == "DAF" || return nothing
    kind == "PCK" && class == 2 && return pckcov(file, classid)
    kind == "CK" && class == 3 && return ckcov(ctx, file, classid, false, "INTERVAL", 0.0, "TDB")
    return nothing
end

"""
    basicstats(d::FrmDiff) -> NamedTuple

Return the summary of `frmdiff -t basic`: the maximum and the average over the epochs of
the rotation angle between the two orientations (`maxangle`, `avgangle`, radians) and,
if angular velocities were computed, of the magnitude of their difference (`maxav`,
`avgav`, rad/s; `NaN` otherwise).
"""
function basicstats(d::FrmDiff)
    @assert !isempty(d.times)
    avg(x) = sum(x) / length(x)
    dav = norm.(d.avdiff)
    return (maxangle=maximum(d.angle), avgangle=avg(d.angle),
        maxav=isempty(dav) ? NaN : maximum(dav), avgav=isempty(dav) ? NaN : avg(dav))
end

function Base.show(io::IO, ::MIME"text/plain", d::FrmDiff)
    s = basicstats(d)
    f(x) = @sprintf("%+.13E", x)
    println(io, "#\n# Comparison of $(length(d.times)) ", d.label1, "\n#")
    println(io, "# with $(length(d.times)) ", d.label2, "\n#")
    println(io, "#    from '", _cal(first(d.times)), "' TDB (", first(d.times), " TDB seconds)")
    println(io, "#    to   '", _cal(last(d.times)), "' TDB (", last(d.times), " TDB seconds)\n#\n")
    println(io, "Absolute difference magnitudes:\n")
    println(io, " "^37, rpad("maximum", 24), "average\n")
    print(io, "  Rotation (rad):             ", rpad(f(s.maxangle), 26), f(s.avgangle))
    if !isempty(d.av1)
        print(io, "\n  Angular velocity (rad/s):   ", rpad(f(s.maxav), 26), f(s.avgav))
    end
end
