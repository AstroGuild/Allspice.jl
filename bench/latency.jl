# Latency benchmark: time to first result in a fresh Julia process.
#
#     julia --startup-file=no --project=bench bench/latency.jl
#
# The driver starts new Julia processes that run the child script below, one process per
# repetition, and reports the median of three repetitions. One extra process runs first
# and is not counted, so that package precompilation is not part of the numbers. The
# rows are, in the order in which the child measures them:
#
#   using NeoSpice                       load time of the package
#   first furnsh + first spkezr          load de440s.bsp and one state, name form
#   first spkder Val(2)                  first jet of order 2
#   using ForwardDiff                    load time of ForwardDiff and of the extension
#   first spkpos with a Dual epoch       first call through the generic path
#   using SPICE                          for comparison: the CSPICE wrapper
#   first SPICE.furnsh + SPICE.spkezr    for comparison

using Printf, Statistics, Dates

include(joinpath(@__DIR__, "..", "test", "kernels.jl"))
using .Kernels: Kernels

const REPS = parse(Int, get(ENV, "NEOSPICE_BENCH_REPS", "3"))

const CHILD = raw"""
const KERNEL = ARGS[1]
const ET = 1.2345e8
t = time_ns(); using NeoSpice; t_using = (time_ns() - t) / 1e9
t = time_ns()
furnsh(KERNEL)
state, lt = spkezr("EARTH", ET, "J2000", "NONE", "MOON")
t_first = (time_ns() - t) / 1e9
t = time_ns(); jet = spkder("EARTH", ET, "J2000", "MOON", Val(2)); t_der = (time_ns() - t) / 1e9
t = time_ns(); using ForwardDiff; t_fd = (time_ns() - t) / 1e9
t = time_ns()
pos, _ = spkpos(399, ForwardDiff.Dual(ET, 1.0), 1, "NONE", 301)
t_dual = (time_ns() - t) / 1e9
@assert ForwardDiff.value.(pos) == state[1:3]
@assert isapprox(ForwardDiff.partials.(pos, 1), state[4:6]; rtol = 1e-12)
t = time_ns(); import SPICE; t_spice = (time_ns() - t) / 1e9
t = time_ns()
SPICE.furnsh(KERNEL)
cstate, _ = SPICE.spkezr("EARTH", ET, "J2000", "NONE", "MOON")
t_cspice = (time_ns() - t) / 1e9
@assert isapprox(state, cstate; rtol = 1e-12)
println("LATENCY using_neospice ", t_using)
println("LATENCY first_furnsh_spkezr ", t_first)
println("LATENCY first_spkder_2 ", t_der)
println("LATENCY using_forwarddiff ", t_fd)
println("LATENCY first_dual ", t_dual)
println("LATENCY using_spice ", t_spice)
println("LATENCY first_cspice ", t_cspice)
"""

const LABELS = [
    "using_neospice" => "`using NeoSpice`",
    "first_furnsh_spkezr" => "first `furnsh` + first `spkezr` (names)",
    "first_spkder_2" => "first `spkder(..., Val(2))`",
    "using_forwarddiff" => "`using ForwardDiff` (with the extension)",
    "first_dual" => "first `spkpos` with a `ForwardDiff.Dual` epoch",
    "using_spice" => "`using SPICE` (for comparison)",
    "first_cspice" => "first `SPICE.furnsh` + `SPICE.spkezr` (for comparison)",
]

"Run the child once in a new process and return its measurements in seconds."
function run_child(script::AbstractString)
    project = dirname(@__FILE__)
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$project $script $(Kernels.DE440S)`
    out = read(cmd, String)
    result = Dict{String,Float64}()
    for line in split(out, '\n')
        parts = split(line)
        length(parts) == 3 && parts[1] == "LATENCY" || continue
        result[parts[2]] = parse(Float64, parts[3])
    end
    @assert length(result) == length(LABELS) "the child did not report every measurement:\n$out"
    return result
end

function main()
    gone = Kernels.missing_kernels(Kernels.DE440S)
    if !isempty(gone)
        printstyled("SKIP  missing kernels: ", join(gone, ", "), "\n"; color=:yellow, bold=true)
        return
    end
    script, io = mktemp()
    write(io, CHILD); close(io)
    println("warm-up process (not counted) ...")
    run_child(script)
    runs = Dict{String,Float64}[]
    for r in 1:REPS
        push!(runs, run_child(script))
        println("repetition ", r, ": ", join((@sprintf("%s=%.3f", k, runs[end][k]) for (k, _) in LABELS), " "))
    end
    io = IOBuffer()
    println(io, "# Latency: time to first result in a fresh process\n")
    println(io, "- date: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println(io, "- Julia ", VERSION, ", ", REPS, " repetitions, one process each; ",
        "packages precompiled before the first repetition\n")
    println(io, "| Step | median s | min s | max s |")
    println(io, "|---|---:|---:|---:|")
    for (key, label) in LABELS
        v = [r[key] for r in runs]
        @printf(io, "| %s | %.3f | %.3f | %.3f |\n", label, median(v), minimum(v), maximum(v))
    end
    text = String(take!(io))
    println()
    println(text)
    dir = joinpath(@__DIR__, "results")
    mkpath(dir)
    write(joinpath(dir, "latency.md"), text)
    println("results written to ", joinpath(dir, "latency.md"))
    rm(script; force=true)
    return nothing
end

main()
