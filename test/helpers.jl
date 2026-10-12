# Shared helpers of the test groups.
#
# CSPICE (through SPICE.jl) has one process-wide kernel state and is not thread safe.
# Each group therefore calls `oracle_load` once at its start and never calls SPICE.jl
# from a spawned task.

using Test, Random, LinearAlgebra, Printf
using StaticArrays
using NeoSpice
import SPICE

include("kernels.jl")
using .Kernels: Kernels

"Seed of every random sample; change it with the environment variable NEOSPICE_TEST_SEED."
const SEED = parse(Int, get(ENV, "NEOSPICE_TEST_SEED", "20261011"))

"Names of the groups that were skipped, with the reason, for the summary of `runtests.jl`."
const SKIPPED = Tuple{String,String}[]

"Measured values that `runtests.jl` prints again at the end, so thresholds can be tightened."
const MEASURED = String[]

"Print `line` now and keep it for the summary."
function measured(line::AbstractString)
    println(line)
    push!(MEASURED, String(line))
    return nothing
end

"""
    group(body, name, kernels...)

Run `body()` inside `@testset name`, or print a skip notice if one of `kernels` is not
on disk. Postcondition: the group is either run or listed in `SKIPPED`.
"""
function group(body, name::AbstractString, kernels...)
    gone = Kernels.missing_kernels(kernels...)
    if !isempty(gone)
        reason = "missing kernels: " * join(gone, ", ")
        printstyled("SKIP  ", name, "  (", reason, ")\n"; color=:yellow, bold=true)
        push!(SKIPPED, (String(name), reason))
        return nothing
    end
    printstyled("GROUP ", name, "\n"; color=:cyan, bold=true)
    @testset "$name" begin
        body()
    end
    return nothing
end

"Skip notice for a part of a group that cannot run (threads, optional kernel)."
function skip_notice(name::AbstractString, reason::AbstractString)
    printstyled("SKIP  ", name, "  (", reason, ")\n"; color=:yellow, bold=true)
    push!(SKIPPED, (String(name), String(reason)))
    return nothing
end

"Clear CSPICE and load `paths` in order. Postcondition: `SPICE.ktotal(\"ALL\") >= length(paths)`."
function oracle_load(paths)
    SPICE.kclear()
    for p in paths
        SPICE.furnsh(p)
    end
    @assert SPICE.ktotal("ALL") >= length(paths)
    return nothing
end

"New explicit NeoSpice context with `paths` loaded in order."
function neo_context(paths)
    ctx = SpiceContext()
    for p in paths
        furnsh(ctx, p)
    end
    @assert ktotal(ctx) >= length(paths)
    return ctx
end

"""
    quantile_of(v, q) -> Float64

Nearest-rank quantile `q` in [0, 1] of the non-empty vector `v` (which is not modified).
"""
function quantile_of(v::AbstractVector{<:Real}, q::Real)
    @assert !isempty(v) "quantile of an empty sample"
    @assert 0 <= q <= 1
    s = sort(v)
    return Float64(s[clamp(ceil(Int, q * length(s)), 1, length(s))])
end

"Quantiles q50, q75, q99 and the maximum of `v`."
qstats(v) = (q50=quantile_of(v, 0.50), q75=quantile_of(v, 0.75), q99=quantile_of(v, 0.99),
    max=Float64(maximum(v)))

"Print one line of measured quantiles with the label `what` and return the statistics."
function report_quantiles(what::AbstractString, v)
    s = qstats(v)
    measured(@sprintf("  %-58s n=%-7d q50=%.3e q75=%.3e q99=%.3e max=%.3e", what, length(v),
        s.q50, s.q75, s.q99, s.max))
    return s
end

"`n` random epochs, uniform in [`t0`, `t1`] (ET seconds)."
function random_epochs(rng::AbstractRNG, n::Integer, t0::Real, t1::Real)
    @assert t1 > t0 && n > 0
    return [t0 + (t1 - t0) * rand(rng) for _ in 1:n]
end

"""
    record_safe_epochs(rng, n; init, record, span, margin=0.1)

`n` random epochs `init + record * (m + u)` with `u` in [`margin`, 1 - `margin`] and
`|epoch| <= span`. Every Chebyshev record of `de440s.bsp` and of the lunar binary PCK
starts at `init` plus a multiple of `record` (4 days), so `epoch ± h` stays inside one
record of every segment if `h < margin * record`. Finite differences are then not
taken across a record boundary, where the high derivatives of a piecewise polynomial
jump.
"""
function record_safe_epochs(rng::AbstractRNG, n::Integer; init::Float64=-4.734072e9,
        record::Float64=345600.0, span::Float64=6.3e8, margin::Float64=0.1)
    @assert 0 < margin < 0.5
    mlo = ceil(Int, (-span - init) / record)
    mhi = floor(Int, (span - init) / record) - 1
    @assert mhi > mlo
    out = Float64[]
    for _ in 1:n
        m = rand(rng, mlo:mhi)
        u = margin + (1 - 2margin) * rand(rng)
        push!(out, init + record * (m + u))
    end
    @assert all(t -> abs(t) <= span + record, out)
    return out
end

"Common SPK coverage of `de440s.bsp` in ET seconds (all bodies share it), with a margin."
const DE440S_COVERAGE = (-4.734072e9 + 86400.0, 4.735368e9 - 86400.0)
"Coverage of `de405.bsp` in ET seconds, with a margin."
const DE405_COVERAGE = (-1.26248112e10 + 86400.0, 6.3472464e9 - 86400.0)

"Body IDs that `spkobj` of CSPICE reports for the SPK file `path`."
function spk_bodies(path::AbstractString)
    ids = SPICE.spkobj!(SPICE.SpiceIntCell(1000), path)
    return sort!(Int[i for i in ids])
end

"True if `f()` throws a `SpiceError` of NeoSpice with `kind`."
function throws_kind(f, kind::Symbol)
    try
        f()
    catch e
        e isa NeoSpice.SpiceError || rethrow()
        e.kind === kind || @info "SpiceError of another kind" expected = kind got = e.kind msg = e.msg
        return e.kind === kind
    end
    return false
end
