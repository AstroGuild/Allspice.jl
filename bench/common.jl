# Shared helpers of the benchmark scripts.

using BenchmarkTools, Statistics, Printf, Dates, Random, LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "..", "test", "kernels.jl"))
using .Kernels: Kernels

"Seconds of sampling per benchmark; the environment variable NEOSPICE_BENCH_SECONDS overrides it."
const BENCH_SECONDS = parse(Float64, get(ENV, "NEOSPICE_BENCH_SECONDS", "2"))

const RESULTS_DIR = joinpath(@__DIR__, "results")

"Stop with a visible message if a kernel that the script needs is not on disk."
function require_kernels(paths...)
    gone = Kernels.missing_kernels(paths...)
    isempty(gone) && return nothing
    printstyled("SKIP  missing kernels: ", join(gone, ", "), "\n"; color=:yellow, bold=true)
    exit(0)
end

"""
    anise_epochs() -> Vector{Float64}

The epochs of the ANISE benchmark: 100 epochs evenly spaced from 1900-01-01 12:00 ET to
2099-01-01 12:00 ET, the end excluded, as ET seconds past J2000.
"""
function anise_epochs()
    j2000 = DateTime(2000, 1, 1, 12)
    t0 = Dates.value(DateTime(1900, 1, 1, 12) - j2000) / 1000
    t1 = Dates.value(DateTime(2099, 1, 1, 12) - j2000) / 1000
    step = (t1 - t0) / 100
    epochs = [t0 + step * k for k in 0:99]
    @assert length(epochs) == 100 && epochs[end] < t1 && epochs[1] == t0
    return epochs
end

"One row of a result table. Times are nanoseconds per query."
struct Row
    group::String
    name::String
    min_ns::Float64
    median_ns::Float64
    allocs::Float64
    bytes::Float64
end

"""
    bench_row(group, name, f, args...; n) -> Row

Benchmark `f(args...)`, which makes `n` queries per call, and return the time per query
(minimum and median over the samples) and the allocations per query.
"""
function bench_row(group::AbstractString, name::AbstractString, f, args...; n::Integer)
    @assert n > 0
    f(args...)                                   # compile, and fail early if it throws
    b = @benchmark $f($(args)...) seconds = BENCH_SECONDS
    row = Row(group, name, minimum(b).time / n, median(b).time / n, allocs(b) / n,
        memory(b) / n)
    @printf("  %-34s %-44s min %9.1f ns  median %9.1f ns  allocs %5.2f\n", group, name,
        row.min_ns, row.median_ns, row.allocs)
    return row
end

"Markdown table of `rows`, with the ratio of each minimum to the minimum of row `base`."
function markdown_table(rows::Vector{Row}; base::Union{Nothing,Row}=nothing)
    io = IOBuffer()
    println(io, "| Case | Engine and call | min ns/query | median ns/query | allocs/query | bytes/query |",
        base === nothing ? "" : " ratio to $(base.name) |")
    println(io, "|---|---|---:|---:|---:|---:|", base === nothing ? "" : "---:|")
    for r in rows
        @printf(io, "| %s | %s | %.1f | %.1f | %.2f | %.0f |", r.group, r.name, r.min_ns,
            r.median_ns, r.allocs, r.bytes)
        base === nothing || @printf(io, " %.2f |", r.min_ns / base.min_ns)
        println(io)
    end
    return String(take!(io))
end

"Text that describes the machine and the run, for the head of a result file."
function run_header(title::AbstractString)
    cpu = try
        first(l for l in eachline("/proc/cpuinfo") if startswith(l, "model name"))
    catch
        "model name : unknown"
    end
    return """
    # $title

    - date: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    - Julia $(VERSION), $(Threads.nthreads()) thread(s), $(Sys.CPU_THREADS) logical CPUs
    - CPU $(strip(split(cpu, ":"; limit=2)[2]))
    - BenchmarkTools, $(BENCH_SECONDS) s of samples per row; time per query is the time
      of one loop over the epochs divided by the number of epochs

    """
end

"Write `text` to `bench/results/<name>.md` (the directory is ignored by git) and return the path."
function write_results(name::AbstractString, text::AbstractString)
    mkpath(RESULTS_DIR)
    path = joinpath(RESULTS_DIR, name * ".md")
    write(path, text)
    println("results written to ", path)
    return path
end
