# spkmerge: merge and subset SPK files.

"""
    MergeSource(path; bodies=nothing, start=-Inf, stop=Inf, comments=false)

A source SPK file of [`spkmerge`](@ref), as a `SOURCE_SPK_KERNEL` block of the setup
file of the NAIF program.

- `bodies`: the NAIF ID codes of the bodies to take, or `nothing` for all bodies.
- `start`, `stop`: the interval to take, in seconds past J2000 TDB.
- `comments`: `true` copies the comments of the file to the merged file.
"""
struct MergeSource
    path::String
    bodies::Union{Nothing,Vector{Int}}
    start::Float64
    stop::Float64
    comments::Bool
    function MergeSource(path::AbstractString; bodies=nothing, start::Real=-Inf,
            stop::Real=Inf, comments::Bool=false)
        start <= stop || throw(ArgumentError("start $start is after stop $stop"))
        b = bodies === nothing ? nothing : Int[Int(x) for x in bodies]
        return new(String(path), b, Float64(start), Float64(stop), comments)
    end
end

"""
    spkmerge(out, sources; ifname=basename(out), comments=String[]) -> Vector{NamedTuple}

Write the SPK file `out` from the segments of the [`MergeSource`](@ref)s `sources`, as
the NAIF program `spkmerge` writes one `SPK_KERNEL`.

A source listed earlier has priority: for each body, a later source gives data only on
the intervals where no earlier source gave data. Within a source, a segment later in the
file has priority, as when the file is loaded. Each segment of the result is the part of
a source segment on one such interval, written by `Allspice.spksub` with the identifier
of the source segment. Intervals of zero length are not written. The segments of each
source are written in order of priority, the last segment of the file first, as the NAIF
program writes them.

`comments` are the first lines of the comment area of `out`; the comments of each source
with `comments=true` follow. `out` must not exist.

The result has one entry per written segment, in file order, with the fields `body`,
`center`, `frame`, `type`, `start`, `stop`, `ident` and `source` (the path).

Throws `ArgumentError` if `sources` is empty or a source is not an SPK file, and
`SpiceError(:nosegmentsfound)` if no source gives data; `out` is then removed.
"""
function spkmerge(out::AbstractString, sources::AbstractVector{MergeSource};
        ifname::AbstractString=basename(out), comments::AbstractVector{<:AbstractString}=String[])
    isempty(sources) && throw(ArgumentError("spkmerge needs at least one source"))
    for s in sources
        getfat(s.path) == ("DAF", "SPK") || throw(ArgumentError("$(s.path) is not an SPK file"))
    end
    lines = String[comments...]
    for s in sources
        s.comments && append!(lines, AllspiceUtilities.comments(s.path))
    end
    ncomch = isempty(lines) ? 0 : sum(l -> length(l) + 1, lines)
    covered = Dict{Int,Window}()
    written = NamedTuple[]
    newh = spkopn(out, ifname, ncomch)
    ok = false
    try
        isempty(lines) || dafac(newh, lines)
        for s in sources
            h = dafopr(s.path)
            try
                pieces = []
                dafbbs(h)
                while daffpa(h)   # highest priority first
                    descr = dafgs(h)
                    body, center, frame, type, t0, t1 = spkuds(descr)
                    s.bodies === nothing || body in s.bodies || continue
                    lo, hi = max(t0, s.start), min(t1, s.stop)
                    lo <= hi || continue
                    w = Window([lo, hi])
                    cov = get(covered, body, Window())
                    rest = wndifd(w, cov)
                    covered[body] = wnunid(cov, w)
                    push!(pieces, (descr, dafgn(h), rest, (body, center, frame, type)))
                end
                for (descr, ident, rest, (body, center, frame, type)) in pieces
                    for i in 1:wncard(rest)
                        a, b = wnfetd(rest, i)
                        a < b || continue
                        spksub(h, descr, ident, a, b, newh)
                        push!(written, (; body, center, frame, type, start=a, stop=b,
                            ident, source=s.path))
                    end
                end
            finally
                dafcls(h)
            end
        end
        spkcls(newh)
        ok = true
    finally
        if !ok
            dafcls(newh)
            rm(out; force=true)
        end
    end
    @assert all(x -> x.start < x.stop, written)
    return written
end
