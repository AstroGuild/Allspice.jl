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

"""
True for the exhaustive run (environment variable NEOSPICE_TEST_FULL=1). The default
run has a budget of 120 s for all groups together; a group keeps its large samples,
long search spans and large files for the exhaustive run.
"""
const FULL = get(ENV, "NEOSPICE_TEST_FULL", "0") == "1"

"""
    nsample(fast, full)

Size of a sample: `fast` in the default run, `full` in the exhaustive run
(`NEOSPICE_TEST_FULL=1`). Requires `0 <= fast <= full`.
"""
function nsample(fast::Integer, full::Integer)
    @assert 0 <= fast <= full
    return FULL ? Int(full) : Int(fast)
end

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

# Tolerances of the state comparison, used by the groups states, identities, rotations,
# aberration and management.

"Thresholds of one comparison class. `rel_*` are relative to the vector norm."
struct StateTol
    rel_q75::Float64
    rel_q99::Float64
    max_pos_km::Float64
    max_vel_kms::Float64
end

# The absolute bound `max_pos_km` holds for a component whose vector norm is so small
# that the bound is above the rounding of the norm. One unit in the last place of a
# coordinate of Neptune or Pluto (4.5e9 to 9e9 km) is 9.5e-7 km, so there the bound is
# `ULP_ALLOWANCE` units in the last place of the norm instead: the measured maximum of
# 1.9e-6 km in the first run is two such units, at a relative difference of 1.8e-15.
const ULP_ALLOWANCE = 4

# Inertial frame that is the frame of the segments. The Chebyshev sums use the
# operations of CSPICE, and on the development machine every state was equal to the
# last bit (maximum difference 0 over 109,200 samples); the bound allows one rounding
# for a CSPICE build with different floating-point contraction.
const TOL_J2000 = StateTol(3e-16, 1e-15, 1e-6, 1e-12)
# Constant rotation on top: one matrix product more.
const TOL_INERTIAL = StateTol(3e-16, 1e-15, 1e-6, 1e-12)
# Body-fixed frames. The angle W uses the operations of CSPICE as well, so only the
# matrix products differ in rounding. Measured: relative position 3.3e-16 at most,
# relative velocity 9.0e-15 at most (q99 3.5e-15), velocity 2.9e-11 km/s at most.
const TOL_BODYFIXED = StateTol(3e-16, 2e-14, 1e-6, 1e-9)

"Lowest norm used as the divisor of a relative difference (km, km/s)."
const POS_FLOOR = 1.0
const VEL_FLOOR = 1e-6

"""
    compare_states(ctx, pairs, epochs, ref) -> NamedTuple

Query every (target, observer) pair of `pairs` at every epoch of `epochs` in frame
`ref` with NeoSpice (integer-ID form on `ctx`) and with CSPICE, and return the vectors
of absolute and relative component differences and of light-time differences.
Precondition: CSPICE has the same kernels loaded as `ctx`.
"""
function compare_states(ctx, pairs, epochs, ref::AbstractString)
    @assert !isempty(pairs) && !isempty(epochs)
    abs_pos = Float64[]; abs_vel = Float64[]
    rel_pos = Float64[]; rel_vel = Float64[]
    rel_lt = Float64[]
    excess_pos = Float64[]      # absolute difference above the rounding allowance
    for (targ, obs) in pairs
        st, so = string(targ), string(obs)
        for et in epochs
            ref_state, ref_lt = SPICE.spkezr(st, et, ref, "NONE", so)
            state, lt = spkezr(ctx, targ, et, ref, "NONE", obs)
            nr = max(norm(@view ref_state[1:3]), POS_FLOOR)
            nv = max(norm(@view ref_state[4:6]), VEL_FLOOR)
            for i in 1:3
                d = abs(state[i] - ref_state[i])
                push!(abs_pos, d); push!(rel_pos, d / nr)
                push!(excess_pos, max(d - ULP_ALLOWANCE * eps(nr), 0.0))
                d = abs(state[i+3] - ref_state[i+3])
                push!(abs_vel, d); push!(rel_vel, d / nv)
            end
            push!(rel_lt, abs(lt - ref_lt) / max(ref_lt, POS_FLOOR / 299792.458))
        end
    end
    @assert length(abs_pos) == 3 * length(pairs) * length(epochs)
    return (; abs_pos, abs_vel, rel_pos, rel_vel, rel_lt, excess_pos)
end

"Print the measured quantiles of `d` under `label` and assert them against `tol`."
function assert_states(label::AbstractString, d, tol::StateTol)
    rp = report_quantiles("$label: position, relative", d.rel_pos)
    rv = report_quantiles("$label: velocity, relative", d.rel_vel)
    ap = report_quantiles("$label: position, absolute km", d.abs_pos)
    av = report_quantiles("$label: velocity, absolute km/s", d.abs_vel)
    lt = report_quantiles("$label: light time, relative", d.rel_lt)
    @testset "$label" begin
        @test rp.q75 <= tol.rel_q75
        @test rp.q99 <= tol.rel_q99
        @test maximum(d.excess_pos) <= tol.max_pos_km
        @test rv.q75 <= tol.rel_q75
        @test rv.q99 <= tol.rel_q99
        @test av.max <= tol.max_vel_kms
        # Light time is |r|/c: the same relative error as the position norm, plus one
        # rounding of the square root and one of the division.
        @test lt.q99 <= max(tol.rel_q99, 4eps())
        @test lt.max <= max(10 * tol.rel_q99, 4eps())
    end
    return nothing
end
