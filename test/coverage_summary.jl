# Line coverage of `src/` and `ext/` from the `.cov` files of a run with code coverage.
#
#     julia --startup-file=no --threads=8 --code-coverage=@src --project=test test/runtests.jl [groups]
#     julia --startup-file=no test/coverage_summary.jl [--keep] [output file]
#
# Several `.cov` files of one source file (one for each process) are added. The script
# prints one line for each source file and the total, writes the same lines to the
# output file if one is given, and deletes the `.cov` files unless `--keep` is given.

"Counts of one `.cov` file: `nothing` for a line that is not code, else the number of executions."
function read_counts(path::AbstractString)
    return map(eachline(path)) do line
        field = strip(first(line, min(9, ncodeunits(line))))
        field == "-" || isempty(field) ? nothing : parse(Int, field)
    end
end

function main(args)
    keep = "--keep" in args
    out = filter(a -> a != "--keep", args)
    root = dirname(@__DIR__)
    merged = Dict{String,Vector{Union{Nothing,Int}}}()
    covfiles = String[]
    for dir in ("src", "ext"), (d, _, files) in walkdir(joinpath(root, dir)), f in files
        endswith(f, ".cov") || continue
        path = joinpath(d, f)
        push!(covfiles, path)
        source = relpath(replace(path, r"\.\d+\.cov$" => ""), root)
        counts = read_counts(path)
        old = get(merged, source, nothing)
        merged[source] = old === nothing ? counts :
            [a === nothing ? b : b === nothing ? a : a + b for (a, b) in zip(old, counts)]
    end
    isempty(merged) && error("no .cov file under src/ or ext/; run the tests with --code-coverage=@src")
    lines = String[]
    hit = total = 0
    for source in sort!(collect(keys(merged)))
        c = filter(!isnothing, merged[source])
        h, n = count(>(0), c), length(c)
        hit += h
        total += n
        push!(lines, string(rpad(source, 34), lpad(h, 6), lpad(n, 6), lpad(round(100 * h / max(n, 1); digits=1), 7)))
    end
    push!(lines, string(rpad("all", 34), lpad(hit, 6), lpad(total, 6), lpad(round(100 * hit / total; digits=1), 7)))
    foreach(println, lines)
    isempty(out) || write(out[1], join(lines, "\n") * "\n")
    keep || foreach(rm, covfiles)
    return nothing
end

main(ARGS)
