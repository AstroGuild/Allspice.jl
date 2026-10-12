# Group 10: windows (sets of disjoint closed intervals) against the SPICE cells of CSPICE.
#
# Feature file: features/windows.feature
#
# Every operation runs on the same random input with NeoSpice (`Window`) and with CSPICE
# (`SPICE.SpiceDoubleCell`), and the endpoint vectors must be equal bit for bit. The
# endpoints are small integers and halves, so that equal endpoints, touching intervals
# and singleton intervals occur often.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of NeoSpice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`NEOSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

"Endpoints of a SPICE double cell as a vector."
cell_values(c) = Float64[c[i] for i in 1:SPICE.card(c)]

"SPICE window (double cell of size `size`) with the valid endpoints `v`."
function cell_window(v::AbstractVector{<:Real}; size::Integer=length(v) + 8)
    c = SPICE.SpiceDoubleCell(size)
    for i in 1:2:length(v)
        SPICE.wninsd!(c, v[i], v[i+1])
    end
    return c
end

"Random pairs `left <= right` (not sorted, overlaps allowed) on a grid of halves in [0, `span`]."
function random_pairs(rng::AbstractRNG, n::Integer; span::Integer=40)
    out = Float64[]
    for _ in 1:n
        a = rand(rng, 0:2span) / 2
        len = rand(rng) < 0.2 ? 0.0 : rand(rng, 0:8) / 2
        push!(out, a, a + len)
    end
    return out
end

"True if `d` holds sorted disjoint intervals (the invariant of a `Window`)."
function valid_window(d::AbstractVector{Float64})
    iseven(length(d)) || return false
    all(d[i] <= d[i+1] for i in 1:2:length(d)) || return false
    return all(d[i] < d[i+1] for i in 2:2:length(d)-1)
end

group("10. windows") do
    rng = MersenneTwister(SEED + 10)
    ntrial = 400

    @apart @testset "wninsd! and wnvald! build the window of CSPICE" begin
        bad = 0
        for _ in 1:ntrial
            pairs = random_pairs(rng, rand(rng, 0:12))
            w = Window()
            c = SPICE.SpiceDoubleCell(length(pairs) + 8)
            for i in 1:2:length(pairs)
                @test wninsd!(w, pairs[i], pairs[i+1]) === w
                SPICE.wninsd!(c, pairs[i], pairs[i+1])
            end
            bad += collect(w) != cell_values(c)
            bad += !valid_window(w.data)
            # `Window(pairs)` and `wnvald!` give the same window from the unsorted pairs.
            bad += Window(pairs) != w
            cv = SPICE.SpiceDoubleCell(length(pairs) + 8)
            foreach(x -> SPICE.appnd(x, cv), pairs)
            SPICE.wnvald!(cv)
            bad += collect(Window(pairs)) != cell_values(cv)
            bad += card(w) != SPICE.card(c)
            bad += wncard(w) != SPICE.wncard(c)
        end
        @test bad == 0
        measured(@sprintf("  %-58s trials=%d mismatches=%d", "windows: wninsd!/wnvald! vs CSPICE", ntrial, bad))
    end

    @apart @testset "set operations are equal to CSPICE" begin
        bad = Dict(k => 0 for k in (:wnunid, :wnintd, :wndifd, :wncomd))
        for _ in 1:ntrial
            a = Window(random_pairs(rng, rand(rng, 0:8)))
            b = Window(random_pairs(rng, rand(rng, 0:8)))
            ca = cell_window(a.data)
            cb = cell_window(b.data)
            u = wnunid(a, b)
            x = wnintd(a, b)
            d = wndifd(a, b)
            bad[:wnunid] += collect(u) != cell_values(SPICE.wnunid(ca, cb))
            bad[:wnintd] += collect(x) != cell_values(SPICE.wnintd(ca, cb))
            bad[:wndifd] += collect(d) != cell_values(SPICE.wndifd(ca, cb))
            @test valid_window(u.data) && valid_window(x.data) && valid_window(d.data)
            # The inputs are not changed.
            @test collect(a) == cell_values(ca) && collect(b) == cell_values(cb)
            left = rand(rng, -4:60) / 2
            right = left + rand(rng, 0:60) / 2
            comp = wncomd(a, left, right)
            bad[:wncomd] += collect(comp) != cell_values(SPICE.wncomd(ca, left, right))
            @test valid_window(comp.data)
        end
        for (k, v) in sort!(collect(bad); by=first)
            @test v == 0
            measured(@sprintf("  %-58s trials=%d mismatches=%d", "windows: $k vs CSPICE", ntrial, v))
        end
    end

    @apart @testset "queries are equal to CSPICE" begin
        bad = 0
        for _ in 1:ntrial
            w = Window(random_pairs(rng, rand(rng, 0:8)))
            c = cell_window(w.data)
            for _ in 1:8
                p = rand(rng, -4:90) / 2
                q = p + rand(rng, 0:6) / 2
                bad += wnelmd(w, p) != SPICE.wnelmd(c, p)
                bad += wnincd(w, p, q) != SPICE.wnincd(c, p, q)
            end
            for n in 1:wncard(w)
                bad += wnfetd(w, n) != SPICE.wnfetd(c, n)
            end
            meas, avg, stddev, isml, ilon = wnsumd(w)
            cmeas, cavg, cstddev, cisml, cilon = SPICE.wnsumd(c)
            bad += (meas, avg, isml, ilon) != (cmeas, cavg, cisml, cilon)
            # CSPICE gives NaN where rounding makes the variance negative; NeoSpice gives 0.
            bad += !(stddev == cstddev || (isnan(cstddev) && stddev == 0.0))
        end
        @test bad == 0
        measured(@sprintf("  %-58s trials=%d mismatches=%d", "windows: wnelmd/wnincd/wnfetd/wnsumd vs CSPICE", ntrial, bad))
    end

    @apart @testset "input checks" begin
        w = Window([1.0, 2.0, 5.0, 6.0])
        @test throws_kind(() -> wninsd!(w, 3.0, 2.0), :badendpoints)
        @test_throws SPICE.SpiceError SPICE.wninsd!(cell_window(w.data), 3.0, 2.0)
        @test collect(w) == [1.0, 2.0, 5.0, 6.0]
        @test throws_kind(() -> wncomd(w, 3.0, 2.0), :badendpoints)
        @test throws_kind(() -> Window([1.0, 2.0, 3.0]), :unmatchendpts)
        @test throws_kind(() -> Window([2.0, 1.0]), :badendpoints)
        @test throws_kind(() -> wnfetd(w, 0), :nointerval)
        @test throws_kind(() -> wnfetd(w, 3), :nointerval)
        @test_throws SPICE.SpiceError SPICE.wnfetd(cell_window(w.data), 3)
        @test wnsumd(Window()) == (0.0, 0.0, 0.0, 0, 0)
        @test card(Window()) == 0 && wncard(Window()) == 0
        # A window is a vector of its endpoints, and `copy` is independent.
        c = copy(w)
        wninsd!(c, 2.0, 5.0)
        @test collect(c) == [1.0, 6.0] && collect(w) == [1.0, 2.0, 5.0, 6.0]
        @test w[3] == 5.0 && length(w) == 4 && card(w) == 4 && card([3, 4, 5]) == 3
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
