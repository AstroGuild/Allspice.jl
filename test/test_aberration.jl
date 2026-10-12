# Group 9: aberration corrections against CSPICE.
#
# Feature file: features/aberration_corrections.feature
#
# Error measure: the Euclidean norm of the difference from the CSPICE vector divided by
# the norm of the CSPICE vector (position and velocity apart), and the relative
# difference of the light time. The table that is printed has the largest value of each
# over all pairs and epochs, per frame and correction, with `NONE` as the baseline.
#
# Why the agreement is near 1e-13 and not 1e-16. A light-time corrected state is the
# difference of two barycentric vectors (target at `et ± lt`, observer at `et`), as in
# `spkltc`. One unit in the last place of a barycentric coordinate (1.5e8 km for the
# Earth and the Moon) is 3e-8 km, which is 8e-14 of the Earth-Moon distance. The two
# libraries sum the Chebyshev series in a different order, so they differ by such units.
#
# Derivatives. `spkder` with a correction returns exact derivatives of the position
# function that CSPICE computes, so it is checked against the central difference of the
# order below (tolerances as in group 3) and against nested `ForwardDiff` of `spkpos`.
# The difference between order 1 and the CSPICE velocity is measured and printed.

import ForwardDiff

const AB_CORRECTIONS = ("LT", "LT+S", "CN", "CN+S", "XLT", "XLT+S", "XCN", "XCN+S")

"Thresholds on the largest relative difference of position, velocity and light time."
# Velocity with `+S`: `spkacs` takes the observer acceleration from the difference of two
# velocities 2 s apart. One unit in the last place of a velocity of 30 km/s moves that
# acceleration by 3e-15 km/s^2, and the derivative of the stellar aberration multiplies
# it by `|r| / c` (15000 s for Neptune): 5e-11 km/s, or 2e-12 of the velocity.
const AB_TOL_INERTIAL = (pos=2e-12, vel=5e-11, lt=2e-12)
# Body-fixed frames: the rounding of the prime-meridian angle far from J2000 (see
# `TOL_BODYFIXED` of group 1) is on top.
const AB_TOL_BODYFIXED = (pos=1e-9, vel=1e-9, lt=2e-12)

const AB_FD_INERTIAL = (h=2000.0, rtol=2e-4)
const AB_FD_FAST = (h=20.0, rtol=2e-5)
"Tolerance of `spkder` against nested ForwardDiff of `spkpos`, orders 0 to 3."
const AB_AD_RTOL = (1e-13, 1e-11, 1e-10, 1e-9)

"Relative difference `|a - b| / |b|`; for `b == 0` the absolute difference."
function ab_reldiff(a, b)
    nb = norm(b)
    return nb == 0 ? norm(a - b) : norm(a - b) / nb
end

"""
    ab_compare(ctx, pairs, epochs, ref, abcorr) -> NamedTuple

Largest relative differences from CSPICE of `spkezr` (position, velocity, light time)
and of `spkpos` (position, light time) over `pairs` and `epochs`.
Precondition: CSPICE has the same kernels loaded as `ctx`.
"""
function ab_compare(ctx, pairs, epochs, ref::AbstractString, abcorr::AbstractString)
    @assert !isempty(pairs) && !isempty(epochs)
    pos = vel = lt = ppos = plt = 0.0
    for (targ, obs) in pairs, et in epochs
        cs, cl = SPICE.spkezr(string(targ), et, ref, abcorr, string(obs))
        s, l = spkezr(ctx, targ, et, ref, abcorr, obs)
        pos = max(pos, ab_reldiff(s[SOneTo(3)], cs[1:3]))
        vel = max(vel, ab_reldiff(s[SVector(4, 5, 6)], cs[4:6]))
        lt = max(lt, abs(l - cl) / cl)
        cp, cpl = SPICE.spkpos(string(targ), et, ref, abcorr, string(obs))
        p, pl = spkpos(ctx, targ, et, ref, abcorr, obs)
        ppos = max(ppos, ab_reldiff(p, cp))
        plt = max(plt, abs(pl - cpl) / cpl)
    end
    return (; pos, vel, lt, ppos, plt)
end

"Nested forward-mode derivatives of `f` at `t`, orders 0 to 3."
function ab_ad_jet(f, t::Float64)
    d1(s) = ForwardDiff.derivative(f, s)
    d2(s) = ForwardDiff.derivative(d1, s)
    d3(s) = ForwardDiff.derivative(d2, s)
    return (f(t), d1(t), d2(t), d3(t))
end

ab_alloc_ezr(ctx, et, ref, ab) = @allocated spkezr(ctx, 301, et, ref, ab, 10)
ab_alloc_pos(ctx, et, ref, ab) = @allocated spkpos(ctx, 301, et, ref, ab, 10)
ab_alloc_der(ctx, et, ref, ab, ::Val{N}) where {N} =
    @allocated spkder(ctx, 301, et, ref, ab, 10, Val(N))
ab_alloc_parse(s) = @allocated Aberration(s)
ab_alloc_ssb(ctx, et) = @allocated spkssb(ctx, 399, et, 1)
ab_alloc_ltc(ctx, et, ab, st) = @allocated spkltc(ctx, 301, et, 1, ab, st)
ab_alloc_aps(ctx, et, ab, st, acc) = @allocated spkaps(ctx, 301, et, 1, ab, st, acc)
ab_alloc_acs(ctx, et, ab) = @allocated spkacs(ctx, 301, et, 1, ab, 399)
ab_alloc_apo(ctx, et, ab, st) = @allocated spkapo(ctx, 301, et, 1, st, ab)
ab_alloc_ltime(ctx, et) = @allocated ltime(ctx, et, 399, "->", 5)
ab_alloc_stelab(p, v) = @allocated stelab(p, v)

"Allocation of the second call of `f(args...)` (the first call compiles)."
function ab_warm(f, args...)
    f(args...)
    return f(args...)
end

group("9. aberration corrections", Kernels.FULL_SET_FILES...) do
    rng = MersenneTwister(SEED + 9)
    set = Kernels.full_set()
    oracle_load(set)
    ctx = neo_context(set)
    # The light time moves the target epoch by up to 5 hours (Neptune), and the stellar
    # aberration needs the observer at ± 1 s: stay one day inside the coverage.
    span = (DE440S_COVERAGE[1] + 86400.0, DE440S_COVERAGE[2] - 86400.0)

    @testset "correction strings" begin
        for (s, uselt, conv, stel, xmit) in (("NONE", false, false, false, false),
                ("LT", true, false, false, false), ("LT+S", true, false, true, false),
                ("CN", true, true, false, false), ("CN+S", true, true, true, false),
                ("XLT", true, false, false, true), ("XLT+S", true, false, true, true),
                ("XCN", true, true, false, true), ("XCN+S", true, true, true, true))
            c = Aberration(s)
            @test (c.uselt, c.converged, c.stellar, c.transmit) == (uselt, conv, stel, xmit)
            @test c.geometric == (s == "NONE")
            @test string(c) == s
            @test isbits(c)
            # Case and blanks are ignored, as `ljucrs` of CSPICE does.
            spaced = " " * join(collect(lowercase(s)), " ") * "  "
            @test Aberration(lowercase(s)) == c
            @test Aberration(spaced) == c
            @test Aberration(SubString(s * "!", 1, length(s))) == c
            @test Aberration(c) === c
            @test Aberration(; lt=uselt, converged=conv, stellar=stel, transmit=xmit) == c
            @test ab_warm(ab_alloc_parse, s) == 0
            @test ab_warm(ab_alloc_parse, spaced) == 0
        end
        @test (@inferred Aberration("cn + s")) == Aberration("CN+S")
        # Values that CSPICE recognizes and rejects, and values that it does not know.
        for bad in ("RL", "RL+S", "XRL", "XRL+S", "S", "XS", "", "   ", "LT+", "LTS", "L",
                "NONE+S", "XNONE", "LT+S+S", "CN+SS", "λT", "LT\t")
            @test throws_kind(() -> Aberration(bad), :badarg)
            @test throws_kind(() -> spkezr(ctx, 301, 0.0, "J2000", bad, 399), :badarg)
            @test throws_kind(() -> spkpos(ctx, 301, 0.0, "J2000", bad, 399), :badarg)
            @test throws_kind(() -> spkder(ctx, 301, 0.0, "J2000", bad, 399, Val(1)), :badarg)
            @test_throws Exception SPICE.spkezr("301", 0.0, "J2000", bad, "399")
        end
        # CSPICE keeps the first five non-blank characters only, so it reads "XLT+SX"
        # as "XLT+S". NeoSpice does not truncate: the value is an error.
        @test SPICE.spkezr("301", 0.0, "J2000", "XLT+SX", "399")[2] > 0
        @test throws_kind(() -> Aberration("XLT+SX"), :badarg)
    end

    pairs = [(301, 399), (399, 301), (10, 399), (399, 10), (5, 301), (4, 0), (3, 199),
        (0, 301), (301, 10), (299, 5), (8, 399), (10, 0)]
    frames = Tuple{String,Any}[("J2000", AB_TOL_INERTIAL), ("ECLIPJ2000", AB_TOL_INERTIAL),
        ("IAU_EARTH", AB_TOL_BODYFIXED), ("IAU_MOON", AB_TOL_BODYFIXED),
        ("MOON_PA_DE440", AB_TOL_BODYFIXED), ("MOON_ME", AB_TOL_BODYFIXED)]

    @testset "states and light time against CSPICE" begin
        epochs = random_epochs(rng, 40, span...)
        measured(@sprintf("  %-14s %-6s %-9s %-9s %-9s %-9s %-9s", "frame", "abcorr",
            "ezr pos", "ezr vel", "ezr lt", "pos pos", "pos lt"))
        for (ref, tol) in frames, ab in ("NONE", AB_CORRECTIONS...)
            d = ab_compare(ctx, pairs, epochs, ref, ab)
            measured(@sprintf("  %-14s %-6s %.2e  %.2e  %.2e  %.2e  %.2e", ref, ab, d.pos,
                d.vel, d.lt, d.ppos, d.plt))
            @test d.pos <= tol.pos
            @test d.vel <= tol.vel
            @test d.lt <= tol.lt
            @test d.ppos <= tol.pos
            @test d.plt <= tol.lt
        end
    end

    @testset "forms of the call" begin
        for et in random_epochs(rng, 5, span...), ab in AB_CORRECTIONS, ref in ("J2000", "IAU_EARTH")
            byid = spkezr(ctx, 301, et, ref, ab, 10)
            @test byid[1] isa SVector{6,Float64} && byid[2] isa Float64
            @test spkezr(ctx, "MOON", et, ref, ab, "SUN") == byid
            @test spkezr(ctx, 301, et, ref, Aberration(ab), 10) == byid
            @test spkezr(ctx, 301, et, ref, " " * lowercase(ab), 10) == byid
            @test spkezr(NeoSpice.kernelset(ctx), 301, et, ref, ab, 10) == byid
            pos, lt = spkpos(ctx, 301, et, ref, ab, 10)
            @test pos isa SVector{3,Float64}
            @test spkpos(ctx, "MOON", et, ref, Aberration(ab), "SUN") == (pos, lt)
            # spkpos follows `spkapo` and spkezr follows `spkltc`: rounding apart.
            @test isapprox(pos, byid[1][SOneTo(3)]; rtol=1e-12)
            @test isapprox(lt, byid[2]; rtol=1e-12)
            # The light time is the norm of the position without stellar aberration,
            # and stellar aberration keeps the norm.
            @test isapprox(lt, norm(pos) / 299792.458; rtol=1e-13)
        end
        # A body relative to itself: zero state and zero light time, as CSPICE.
        for ab in AB_CORRECTIONS, ref in ("J2000", "IAU_MOON")
            state, lt = spkezr(ctx, 399, 1e8, ref, ab, 399)
            @test all(iszero, state) && lt == 0
            @test all(iszero, SPICE.spkezr("399", 1e8, ref, ab, "399")[1])
            @test all(iszero, spkpos(ctx, 399, 1e8, ref, ab, 399)[1])
            @test all(all(iszero, d) for d in spkder(ctx, 399, 1e8, ref, ab, 399, Val(2)))
        end
        # Reception and transmission differ, and the correction has the expected size:
        # the Sun seen from the Earth moves by about 20 arc seconds of aberration.
        geo = spkpos(ctx, 10, 0.0, "J2000", "NONE", 399)[1]
        app = spkpos(ctx, 10, 0.0, "J2000", "LT+S", 399)[1]
        xmt = spkpos(ctx, 10, 0.0, "J2000", "XLT+S", 399)[1]
        @test 15 < acos(clamp(dot(geo, app) / (norm(geo) * norm(app)), -1, 1)) * 206264.8 < 25
        @test norm(app - xmt) > 1e3
    end

    @testset "default context" begin
        kclear()
        foreach(furnsh, set)
        for ab in ("LT+S", "XCN"), ref in ("J2000", "IAU_MOON")
            @test spkezr("MOON", 1e7, ref, ab, "EARTH") == spkezr(ctx, "MOON", 1e7, ref, ab, "EARTH")
            @test spkpos(301, 1e7, ref, Aberration(ab), 399) == spkpos(ctx, 301, 1e7, ref, ab, 399)
            @test spkder(301, 1e7, ref, ab, 399, Val(2)) == spkder(ctx, 301, 1e7, ref, ab, 399, Val(2))
        end
        @test spkssb(399, 1e7, "J2000") == spkssb(ctx, 399, 1e7, "J2000")
        @test spkacs(301, 1e7, "J2000", "CN+S", 399) == spkacs(ctx, 301, 1e7, "J2000", "CN+S", 399)
        @test ltime(1e7, 399, "->", 5) == ltime(ctx, 1e7, 399, "->", 5)
        kclear()
    end

    @testset "low-level calls against CSPICE" begin
        epochs = random_epochs(rng, 15, span...)
        worst = Dict{String,Float64}()
        note(name, x) = (worst[name] = max(get(worst, name, 0.0), x); x)
        tol = 1e-11
        for et in epochs, ref in ("J2000", "ECLIPJ2000"), (targ, obs) in ((301, 399), (5, 301), (399, 10), (10, 0))
            cssb = SPICE.spkssb(obs, et, ref)
            ssb = spkssb(ctx, obs, et, ref)
            @test ssb isa SVector{6,Float64}
            @test note("spkssb", ab_reldiff(ssb, cssb)) <= 1e-13
            @test ssb == spkgeo(ctx, obs, et, ref, 0)[1]
            acc = SVector(1e-6, -2e-6, 5e-7)
            for ab in ("NONE", AB_CORRECTIONS...)
                cs, cl, cd = SPICE.spkltc(targ, et, ref, ab, cssb)
                s, l, d = spkltc(ctx, targ, et, ref, ab, cssb)
                @test note("spkltc state", ab_reldiff(s[SOneTo(3)], cs[1:3]) + ab_reldiff(s[SVector(4, 5, 6)], cs[4:6])) <= tol
                @test note("spkltc lt", abs(l - cl) / cl) <= tol
                @test note("spkltc dlt", abs(d - cd) / max(abs(cd), 1e-9)) <= 1e-8
                @test spkltc(ctx, targ, et, ref, Aberration(ab), ssb) == spkltc(ctx, targ, et, ref, ab, ssb)

                cs, cl, cd = SPICE.spkaps(targ, et, ref, ab, cssb, vcat(acc, zeros(3)))   # SPICE.jl wants 6
                s, l, d = spkaps(ctx, targ, et, ref, ab, cssb, acc)
                @test note("spkaps state", ab_reldiff(s[SOneTo(3)], cs[1:3]) + ab_reldiff(s[SVector(4, 5, 6)], cs[4:6])) <= tol
                @test note("spkaps lt", abs(l - cl) / cl) <= tol
                @test note("spkaps dlt", abs(d - cd) / max(abs(cd), 1e-9)) <= 1e-8

                cs, cl, cd = SPICE.spkacs(targ, et, ref, ab, obs)
                s, l, d = spkacs(ctx, targ, et, ref, ab, obs)
                @test note("spkacs state", ab_reldiff(s[SOneTo(3)], cs[1:3]) + ab_reldiff(s[SVector(4, 5, 6)], cs[4:6])) <= tol
                @test note("spkacs lt", abs(l - cl) / cl) <= tol
                @test note("spkacs dlt", abs(d - cd) / max(abs(cd), 1e-9)) <= 1e-8
                # spkacs in an inertial frame is the state of spkezr. (With NONE, spkacs
                # takes the difference of two barycentric states and spkezr does not.)
                ab == "NONE" || @test (s, l) == spkezr(ctx, targ, et, ref, ab, obs)

                cp, cpl = SPICE.spkapo(targ, et, ref, cssb, ab)
                p, pl = spkapo(ctx, targ, et, ref, cssb, ab)
                @test note("spkapo position", ab_reldiff(p, cp)) <= tol
                @test note("spkapo lt", abs(pl - cpl) / cpl) <= tol
                if ab != "NONE"
                    q, ql = spkpos(ctx, targ, et, ref, ab, obs)
                    @test isapprox(p, q; rtol=1e-12) && isapprox(pl, ql; rtol=1e-12)
                end
            end
            for dir in ("->", "<-")
                cet, cel = SPICE.ltime(et, obs, dir, targ)
                tt, el = ltime(ctx, et, obs, dir, targ)
                @test note("ltime elapsed", abs(el - cel) / cel) <= tol
                @test abs(tt - cet) <= 4 * eps(abs(et)) + tol * cel
            end
            # Stellar aberration of a position, with the true observer velocity and with
            # a velocity near the limit.
            p0 = spkpos(ctx, targ, et, ref, "LT", obs)[1]
            for v in (ssb[SVector(4, 5, 6)], SVector(1e5, -2e5, 3e4), SVector(0.0, 0.0, 0.0))
                capp = SPICE.stelab(collect(p0), collect(v))
                app = stelab(p0, v)
                @test app isa SVector{3,Float64}
                @test note("stelab", ab_reldiff(app, capp)) <= 1e-14
                @test isapprox(norm(app), norm(p0); rtol=1e-14)
                @test note("stlabx", ab_reldiff(stlabx(p0, v), SPICE.stelab(collect(p0), collect(-v)))) <= 1e-14
            end
        end
        for name in sort!(collect(keys(worst)))
            measured(@sprintf("  %-58s max=%.3e", name * ", relative", worst[name]))
        end
        # Errors of the low-level calls.
        st = spkssb(ctx, 399, 0.0, "J2000")
        @test throws_kind(() -> spkltc(ctx, 301, 0.0, "IAU_EARTH", "LT", st), :badarg)
        @test throws_kind(() -> spkaps(ctx, 301, 0.0, "IAU_EARTH", "LT", st, zeros(3)), :badarg)
        @test throws_kind(() -> spkacs(ctx, 301, 0.0, "IAU_EARTH", "LT", 399), :badarg)
        @test throws_kind(() -> spkapo(ctx, 301, 0.0, "IAU_EARTH", st, "LT"), :badarg)
        @test throws_kind(() -> spkltc(ctx, 301, 0.0, "J2000", "LT", zeros(5)), :badarg)
        @test throws_kind(() -> spkaps(ctx, 301, 0.0, "J2000", "LT", st, zeros(2)), :badarg)
        @test throws_kind(() -> spkltc(ctx, 301, 0.0, "J2000", "S", st), :badarg)
        @test throws_kind(() -> ltime(ctx, 0.0, 399, "<>", 301), :badarg)
        @test throws_kind(() -> stelab(SVector(1.0, 0.0, 0.0), SVector(3e5, 0.0, 0.0)), :badarg)
        @test throws_kind(() -> stelab([1.0, 0.0], [0.0, 0.0, 0.0]), :badarg)
        # An epoch whose light-time corrected epoch leaves the coverage.
        @test throws_kind(() -> spkezr(ctx, 8, DE440S_COVERAGE[1] - 86000.0, "J2000", "LT", 399), :nocoverage)
    end

    cases = [(301, 399, "J2000", AB_FD_INERTIAL), (399, 10, "J2000", AB_FD_INERTIAL),
        (10, 301, "ECLIPJ2000", AB_FD_INERTIAL), (5, 399, "J2000", AB_FD_INERTIAL),
        (301, 399, "IAU_EARTH", AB_FD_FAST), (10, 399, "IAU_EARTH", AB_FD_FAST),
        (399, 301, "MOON_PA_DE440", AB_FD_INERTIAL), (10, 301, "IAU_MOON", AB_FD_INERTIAL),
        (4, 399, "IAU_MOON", AB_FD_INERTIAL)]
    depochs = record_safe_epochs(rng, 6)

    @testset "spkder: order 0 and the CSPICE velocity" begin
        measured("  spkder order 1 against the velocity of CSPICE spkezr, largest relative difference")
        for ab in AB_CORRECTIONS
            w0 = 0.0
            w1i = 0.0
            w1b = 0.0
            for (targ, obs, ref, _) in cases, et in depochs
                jet = spkder(ctx, targ, et, ref, ab, obs, Val(1))
                @test jet isa NTuple{2,SVector{3,Float64}}
                state, _ = spkezr(ctx, targ, et, ref, ab, obs)
                cstate, _ = SPICE.spkezr(string(targ), et, ref, ab, string(obs))
                w0 = max(w0, ab_reldiff(jet[1], state[SOneTo(3)]))
                d = ab_reldiff(jet[2], cstate[4:6])
                ref in ("J2000", "ECLIPJ2000") ? (w1i = max(w1i, d)) : (w1b = max(w1b, d))
                # A lower order is a prefix of a higher order.
                full = spkder(ctx, targ, et, ref, ab, obs, Val(3))
                @test all(isapprox(jet[k], full[k]; rtol=1e-12) for k in 1:2)
                @test isapprox(spkder(ctx, targ, et, ref, ab, obs, Val(0))[1], full[1]; rtol=1e-13)
            end
            measured(@sprintf("    %-6s order 0 vs spkezr %.2e; order 1 vs CSPICE velocity: inertial %.2e, body-fixed %.2e",
                ab, w0, w1i, w1b))
            @test w0 <= 1e-12
            # `CN` forms: rounding. `LT` forms: CSPICE takes the light-time rate of the
            # converged equation for a position that is one iteration from the geometric
            # light time, which is a relative difference of order (v/c) * (range rate/c).
            @test w1i <= (occursin("CN", ab) ? 5e-12 : 1e-6)
            @test w1b <= (occursin("CN", ab) ? 5e-11 : 1e-5)
        end
        # NONE through the corrected form is the geometric form.
        @test spkder(ctx, 301, depochs[1], "IAU_EARTH", "NONE", 399, Val(3)) ==
              spkder(ctx, 301, depochs[1], "IAU_EARTH", 399, Val(3))
        @test throws_kind(() -> spkder(ctx, 301, 0.0, "J2000", "LT", 399, Val(-1)), :badarg)
    end

    @testset "spkder: finite differences" begin
        for (targ, obs, ref, fd) in cases, ab in AB_CORRECTIONS
            worst = zeros(3)
            wpos = 0.0
            for et in depochs
                jet(t) = spkder(ctx, targ, t, ref, ab, obs, Val(3))
                j0 = jet(et); jp = jet(et + fd.h); jm = jet(et - fd.h)
                worst .= max.(worst, [ab_reldiff((jp[k] - jm[k]) / (2 * fd.h), j0[k+1]) for k in 1:3])
                # Order 1 against the central difference of the position of spkpos.
                pp = spkpos(ctx, targ, et + fd.h, ref, ab, obs)[1]
                pm = spkpos(ctx, targ, et - fd.h, ref, ab, obs)[1]
                wpos = max(wpos, ab_reldiff((pp - pm) / (2 * fd.h), j0[2]))
            end
            ab in ("LT", "XCN+S") && measured(@sprintf(
                "  FD  %3d rel %3d %-13s %-6s orders 1..3: %s; order 1 vs spkpos %.2e (tol %.0e)",
                targ, obs, ref, ab, join((@sprintf("%.2e", w) for w in worst), " "), wpos, fd.rtol))
            @test all(<=(fd.rtol), worst)
            @test wpos <= fd.rtol
        end
    end

    @testset "spkder: nested ForwardDiff of spkpos (dual-number epoch)" begin
        probe, plt = spkpos(ctx, 301, ForwardDiff.Dual(depochs[1], 1.0), "IAU_EARTH", "CN+S", 399)
        @test eltype(probe) <: ForwardDiff.Dual && plt isa ForwardDiff.Dual
        @test ForwardDiff.value.(probe) == spkpos(ctx, 301, depochs[1], "IAU_EARTH", "CN+S", 399)[1]
        for ab in AB_CORRECTIONS
            worst = zeros(4)
            for (targ, obs, ref, _) in cases, et in depochs[1:3]
                pos(t) = spkpos(ctx, targ, t, ref, ab, obs)[1]
                ad = ab_ad_jet(pos, et)
                jet = spkder(ctx, targ, et, ref, ab, obs, Val(3))
                worst .= max.(worst, [ab_reldiff(jet[k], ad[k]) for k in 1:4])
            end
            measured(@sprintf("  AD  %-6s orders 0..3: %s", ab,
                join((@sprintf("%.2e", w) for w in worst), " ")))
            @test all(worst[k] <= AB_AD_RTOL[k] for k in 1:4)
        end
        # A dual-number epoch through spkezr and through spkder itself.
        et = depochs[2]
        for ab in ("LT+S", "XCN"), ref in ("J2000", "IAU_MOON")
            jet = spkder(ctx, 10, et, ref, ab, 301, Val(3))
            dstate = spkezr(ctx, 10, ForwardDiff.Dual(et, 1.0), ref, ab, 301)[1]
            @test isapprox(ForwardDiff.value.(dstate), spkezr(ctx, 10, et, ref, ab, 301)[1]; rtol=1e-14)
            @test isapprox(ForwardDiff.partials.(dstate, 1)[SOneTo(3)], jet[2]; rtol=1e-6)
            d3 = ForwardDiff.derivative(t -> spkder(ctx, 10, t, ref, ab, 301, Val(2))[3], et)
            @test isapprox(d3, jet[4]; rtol=1e-10)
        end
    end

    @testset "zero allocation and inferred types" begin
        et = 1.2345e8
        refs = (1, 17, Int(namfrm(ctx, "IAU_EARTH")), Int(namfrm(ctx, "MOON_PA_DE440")),
            Int(namfrm(ctx, "MOON_ME")))
        for ref in refs, ab in ("NONE", AB_CORRECTIONS...)
            @test ab_warm(ab_alloc_ezr, ctx, et, ref, ab) == 0
            @test ab_warm(ab_alloc_pos, ctx, et, ref, ab) == 0
            @test ab_warm(ab_alloc_ezr, ctx, et, ref, Aberration(ab)) == 0
            for N in 0:3
                @test ab_warm(ab_alloc_der, ctx, et, ref, ab, Val(N)) == 0
            end
        end
        st = spkssb(ctx, 399, et, 1)
        acc = SVector(1e-6, 2e-6, 3e-6)
        @test ab_warm(ab_alloc_ssb, ctx, et) == 0
        @test ab_warm(ab_alloc_ltime, ctx, et) == 0
        @test ab_warm(ab_alloc_stelab, st[SOneTo(3)], st[SVector(4, 5, 6)]) == 0
        for ab in ("LT", "CN+S", "XLT+S")
            @test ab_warm(ab_alloc_ltc, ctx, et, ab, st) == 0
            @test ab_warm(ab_alloc_aps, ctx, et, ab, st, acc) == 0
            @test ab_warm(ab_alloc_acs, ctx, et, ab) == 0
            @test ab_warm(ab_alloc_apo, ctx, et, ab, st) == 0
        end
        ks = NeoSpice.kernelset(ctx)
        for ref in (1, refs[3]), ab in ("LT", "XCN+S")
            @test (@inferred spkezr(ks, 301, et, ref, ab, 399)) isa Tuple{SVector{6,Float64},Float64}
            @test (@inferred spkpos(ks, 301, et, ref, ab, 399)) isa Tuple{SVector{3,Float64},Float64}
            @test (@inferred spkder(ks, 301, et, ref, ab, 399, Val(3))) isa NTuple{4,SVector{3,Float64}}
            @test (@inferred spkezr(ctx, 301, et, ref, Aberration(ab), 399)) isa Tuple{SVector{6,Float64},Float64}
        end
        @test (@inferred spkssb(ks, 399, et, 1)) isa SVector{6,Float64}
        @test (@inferred spkltc(ks, 301, et, 1, "CN", st)) isa Tuple{SVector{6,Float64},Float64,Float64}
        @test (@inferred spkaps(ks, 301, et, 1, "CN+S", st, acc)) isa Tuple{SVector{6,Float64},Float64,Float64}
        @test (@inferred spkacs(ks, 301, et, 1, "CN+S", 399)) isa Tuple{SVector{6,Float64},Float64,Float64}
        @test (@inferred spkapo(ks, 301, et, 1, st, "CN+S")) isa Tuple{SVector{3,Float64},Float64}
        @test (@inferred ltime(ks, et, 399, "<-", 301)) isa Tuple{Float64,Float64}
        @test (@inferred stelab(st[SOneTo(3)], st[SVector(4, 5, 6)])) isa SVector{3,Float64}
    end
end
