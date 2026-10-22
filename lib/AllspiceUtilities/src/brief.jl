# brief: coverage summary of SPK and binary PCK files.

"""
    Coverage

The coverage of one body of an SPK file or one frame of a binary PCK file, as
[`brief`](@ref) returns it.

- `path`: the kernel file.
- `kind`: `:spk` or `:pck`.
- `id`: NAIF ID code of the body (SPK) or frame class ID (PCK).
- `name`: name of the body, `""` for a body without a name and for a PCK frame.
- `window`: the coverage, as a `Window` of seconds past J2000 TDB.
"""
struct Coverage
    path::String
    kind::Symbol
    id::Int
    name::String
    window::Window
end

"""
    brief(paths...; ctx=Allspice.default_context()) -> Vector{Coverage}

Return the coverage of each body of the SPK files and of each frame of the binary PCK
files in `paths`, as the NAIF program `brief`.

The result has one [`Coverage`](@ref) per file and body or frame, in the order of
`paths` and then of the ID codes. Text kernels in `paths` are skipped. Body names come
from `ctx`, so a text kernel loaded there with `NAIF_BODY_NAME` names other bodies.
[`printbrief`](@ref) prints the result in the layout of `brief`.

Throws `ArgumentError` for a binary file that is neither an SPK nor a binary PCK file,
and the `SpiceError` of `Allspice.spkobj` for a file that cannot be read.

# Examples
```julia
julia> cov = brief("de440s.bsp");

julia> cov[1].id, cov[1].name
(1, "MERCURY BARYCENTER")
```
"""
function brief(paths::AbstractString...; ctx::SpiceContext=Allspice.default_context())
    result = Coverage[]
    for p in paths
        arch, kind = getfat(p)
        if arch == "DAF" && kind == "SPK"
            for id in spkobj(p)
                push!(result, Coverage(p, :spk, id, _bodyname(ctx, id), spkcov(p, id)))
            end
        elseif arch == "DAF" && kind == "PCK"
            for id in pckfrm(p)
                push!(result, Coverage(p, :pck, id, "", pckcov(p, id)))
            end
        elseif arch != "KPL"
            throw(ArgumentError("brief reads SPK and binary PCK files; $p is $arch/$kind"))
        end
    end
    return result
end

"""
    printbrief([io::IO=stdout,] cov::AbstractVector{Coverage}) -> Nothing

Print the result of [`brief`](@ref) in the layout of the NAIF program `brief`: for each
file, the bodies or frames with the same coverage together, then the intervals of that
coverage in calendar TDB.
"""
function printbrief(io::IO, cov::AbstractVector{Coverage})
    for path in unique(c.path for c in cov)
        rows = filter(c -> c.path == path, cov)
        println(io, "\nSummary for: ", path, "\n")
        for w in unique(c.window for c in rows)
            group = filter(c -> c.window == w, rows)
            label = group[1].kind === :spk ? "Bodies: " : "Frames: "
            items = [isempty(c.name) ? string(c.id) : "$(c.name) ($(c.id))" for c in group]
            for (k, item) in enumerate(items)
                print(io, k == 1 ? label : (k % 3 == 1 ? " "^length(label) : ""),
                    rpad(item, 24))
                (k % 3 == 0 || k == length(items)) && println(io)
            end
            println(io, " "^8, rpad("Start of Interval (ET)", 36), "End of Interval (ET)")
            println(io, " "^8, rpad("-"^29, 36), "-"^29)
            for i in 1:wncard(w)
                a, b = wnfetd(w, i)
                println(io, " "^8, rpad(_cal(a), 36), _cal(b))
            end
            println(io)
        end
    end
    return nothing
end
printbrief(cov::AbstractVector{Coverage}) = printbrief(stdout, cov)
