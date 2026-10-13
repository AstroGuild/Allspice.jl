# Geometry finder, quantities from states: `gfdist`, `gfrr`, `gfsep`, `gfpa` and the
# general entry `gfevnt`, against the same calls of CSPICE. For each quantity: every
# relation, several target and observer pairs, every aberration correction that CSPICE
# allows, several steps and adjust values, confinement windows of one and of several
# intervals, searches with no result and searches whose events are on the end points of
# the confinement window. The number of intervals must be equal and each end point must
# agree within 2e-6 s (twice the tolerance); the largest difference is reported for
# each quantity and family of relations.
#
# The default run compares one case in ten of each product of options (`gf_pick`) and
# does not make the table of search times; the exhaustive run (`ALLSPICE_TEST_FULL=1`)
# compares every case and makes the table.
#
# An event can be ill conditioned, and then the two libraries cannot agree to the
# tolerance, because their states differ by rounding. The epoch of an extremum is the
# root of the test "is decreasing", and CSPICE computes that test from numerical
# derivatives whose rounding noise can be large against the slope of the derivative:
#
# - range rate: the acceleration is a difference of two velocities two seconds apart, so
#   the sign of the derivative of the range rate is noise in an interval of up to some
#   hundredths of a second around the extremum (slow bodies have the widest interval);
# - angular separation with stellar aberration: the velocity of the corrected position
#   holds the acceleration of the observer, which CSPICE takes from a difference of
#   velocities too.
#
# Inside that interval the bisection of CSPICE follows its own rounding and the one of
# Allspice follows its own. For these two quantities an end point that differs by more
# than 2e-6 s is accepted only if the quantity has the same value at the two epochs, to
# the rounding of the states (`gf_compare` with `value`).
#
# Behaviour specification: test/features/gf_quantities.feature

isdefined(@__MODULE__, :GFO) || include("gforacle.jl")

const GF_DAY = 86400.0

"Random confinement window of `n` intervals in a span of `span` seconds that starts between 1990 and 2035."
function gf_confinement(rng::AbstractRNG, n::Integer, span::Real)
    t0 = -3.0e8 + 1.4e9 * rand(rng)
    n == 1 && return Window([t0, t0 + span])
    return gf_random_window(rng, n, t0, t0 + span)
end

@gf_light group("gf quantities: distance, range rate, separation, phase angle", Kernels.FULL_SET_FILES...) do
    oracle_load(Kernels.full_set())
    ctx = neo_context(Kernels.full_set())
    ks = Allspice.kernelset(ctx)
    rng = MersenneTwister(SEED + 6200)
    et(s) = SPICE.str2et(s)

    # (target, observer, steps, span of the confinement window)
    pairs = (("MOON", "EARTH", (GF_DAY, 6 * 3600.0), 120 * GF_DAY),
        ("SUN", "EARTH", (5 * GF_DAY, 20 * GF_DAY), 800 * GF_DAY),
        ("MARS BARYCENTER", "EARTH", (10 * GF_DAY, 30 * GF_DAY), 1500 * GF_DAY),
        ("VENUS BARYCENTER", "SUN", (10 * GF_DAY,), 600 * GF_DAY),
        ("EARTH", "MOON", (2 * GF_DAY,), 90 * GF_DAY))

    @testset "gfdist against CSPICE" begin
        i = 0
        for (targ, obs, steps, span) in pairs, abcorr in GF_ABCORRS, (rel, wadj) in GF_VARIANTS
            gf_pick(i += 1, 10) || continue
            cn = gf_confinement(rng, rand(rng, (1, 1, 3)), span)
            step = rand(rng, steps)
            q = Allspice.GFDistance(ks, targ, abcorr, obs)
            refval, adjust = gf_reference(rng, q, cn, rel, wadj)
            w = gfdist(ctx, targ, abcorr, obs, rel, refval, adjust, step, cn)
            ref = GFO.gfdist(targ, abcorr, obs, rel, refval, adjust, step, collect(cn))
            @test gf_compare("gfdist " * gf_family(rel, adjust), w, ref)
        end
        # integer IDs, an `Aberration` value, the form on a snapshot and on the default context
        cn = Window([et("2024 JAN 01"), et("2024 JUL 01")])
        a = gfdist(ctx, "MOON", "LT+S", "EARTH", "=", 3.9e5, 0.0, GF_DAY, cn)
        @test gfdist(ctx, 301, Aberration("LT+S"), 399, "=", 3.9e5, 0.0, GF_DAY, cn) == a
        @test gfdist(ks, "moon", " lt+s ", "399", "=", 3.9e5, 0.0, GF_DAY, cn) == a
        @test !isempty(a)
    end

    @testset "gfrr against CSPICE" begin
        i = 0
        for (targ, obs, steps, span) in pairs, abcorr in GF_ABCORRS, (rel, wadj) in GF_VARIANTS
            gf_pick(i += 1, 10) || continue
            cn = gf_confinement(rng, rand(rng, (1, 1, 3)), span)
            step = rand(rng, steps)
            q = Allspice.GFRangeRate(ks, targ, abcorr, obs)
            refval, adjust = gf_reference(rng, q, cn, rel, wadj)
            w = gfrr(ctx, targ, abcorr, obs, rel, refval, adjust, step, cn)
            ref = GFO.gfrr(targ, abcorr, obs, rel, refval, adjust, step, collect(cn))
            @test gf_compare("gfrr " * gf_family(rel, adjust), w, ref;
                value=t -> Allspice._gf_value(q, t), vtol=64 * eps(30.0))
        end
    end

    @testset "gfsep against CSPICE" begin
        # (target 1, shape 1, frame 1, target 2, shape 2, frame 2, observer, steps, span)
        triples = (("MOON", "POINT", "NULL", "SUN", "POINT", "NULL", "EARTH", (GF_DAY, 3 * GF_DAY), 150 * GF_DAY),
            ("MOON", "SPHERE", "IAU_MOON", "SUN", "SPHERE", "IAU_SUN", "EARTH", (GF_DAY,), 150 * GF_DAY),
            ("VENUS BARYCENTER", "POINT", "NULL", "SUN", "SPHERE", "IAU_SUN", "EARTH", (10 * GF_DAY,), 900 * GF_DAY),
            ("MARS BARYCENTER", "POINT", "NULL", "JUPITER BARYCENTER", "POINT", "NULL", "EARTH", (10 * GF_DAY,), 1500 * GF_DAY),
            ("EARTH", "SPHERE", "IAU_EARTH", "MOON", "SPHERE", "IAU_MOON", "SUN", (2 * GF_DAY,), 120 * GF_DAY))
        i = 0
        for (t1, s1, f1, t2, s2, f2, obs, steps, span) in triples, abcorr in GF_ABCORRS,
            (rel, wadj) in GF_VARIANTS

            gf_pick(i += 1, 10) || continue
            cn = gf_confinement(rng, rand(rng, (1, 1, 3)), span)
            step = rand(rng, steps)
            q = Allspice.GFSeparation(ks, t1, s1, f1, t2, s2, f2, abcorr, obs)
            refval, adjust = gf_reference(rng, q, cn, rel, wadj)
            w = gfsep(ctx, t1, s1, f1, t2, s2, f2, abcorr, obs, rel, refval, adjust, step, cn)
            ref = GFO.gfsep(t1, s1, f1, t2, s2, f2, abcorr, obs, rel, refval, adjust, step, collect(cn))
            @test gf_compare("gfsep " * gf_family(rel, adjust), w, ref;
                value=t -> Allspice._gf_value(q, t), vtol=64 * eps(Float64(pi)))
        end
    end

    @testset "gfpa against CSPICE" begin
        # (target, illuminator, observer, steps, span); reception corrections only
        triples = (("MOON", "SUN", "EARTH", (GF_DAY, 3 * GF_DAY), 150 * GF_DAY),
            ("VENUS BARYCENTER", "SUN", "EARTH", (10 * GF_DAY,), 900 * GF_DAY),
            ("MARS BARYCENTER", "SUN", "EARTH", (10 * GF_DAY, 30 * GF_DAY), 1500 * GF_DAY),
            ("EARTH", "SUN", "MOON", (2 * GF_DAY,), 120 * GF_DAY))
        i = 0
        for (targ, illum, obs, steps, span) in triples, abcorr in ("NONE", "LT", "LT+S", "CN", "CN+S"),
            (rel, wadj) in GF_VARIANTS, trial in 1:2

            gf_pick(i += 1, 11) || continue
            cn = gf_confinement(rng, rand(rng, (1, 1, 3)), span)
            step = rand(rng, steps)
            q = Allspice.GFPhase(ks, targ, illum, abcorr, obs)
            refval, adjust = gf_reference(rng, q, cn, rel, wadj)
            w = gfpa(ctx, targ, illum, abcorr, obs, rel, refval, adjust, step, cn)
            ref = GFO.gfpa(targ, illum, abcorr, obs, rel, refval, adjust, step, collect(cn))
            @test gf_compare("gfpa " * gf_family(rel, adjust), w, ref)
        end
        cn = Window([0.0, 100 * GF_DAY])
        # a transmission correction (one branch of the code: the flag of the parsed correction)
        @test throws_kind(() -> gfpa(ctx, "MOON", "SUN", "XLT+S", "EARTH", "=", 1.0, 0.0, GF_DAY, cn),
            :invalidoption)
        @test GFO.kind(() -> GFO.gfpa("MOON", "SUN", "XLT+S", "EARTH", "=", 1.0, 0.0, GF_DAY,
            collect(cn))) === :invalidoption
    end

    @testset "searches with no result" begin
        cn = Window([et("2020 MAR 01"), et("2020 SEP 01")])
        c = collect(cn)
        for (rel, refval) in (("=", 1.0e6), ("<", 1.0e5), (">", 1.0e6), ("=", 10.0))
            @test isempty(gfdist(ctx, "MOON", "LT", "EARTH", rel, refval, 0.0, GF_DAY, cn))
            @test isempty(GFO.gfdist("MOON", "LT", "EARTH", rel, refval, 0.0, GF_DAY, c))
        end
        @test isempty(gfrr(ctx, "MOON", "NONE", "EARTH", ">", 5.0, 0.0, GF_DAY, cn))
        @test isempty(GFO.gfrr("MOON", "NONE", "EARTH", ">", 5.0, 0.0, GF_DAY, c))
        @test isempty(gfsep(ctx, "MOON", "POINT", "NULL", "SUN", "POINT", "NULL", "NONE", "EARTH",
            ">", 3.2, 0.0, GF_DAY, cn))
        @test isempty(gfpa(ctx, "MOON", "SUN", "LT+S", "EARTH", "<", -0.1, 0.0, GF_DAY, cn))
        # no local extremum of the lunar distance in two days far from the apsides
        peri = gfdist(ctx, "MOON", "NONE", "EARTH", "LOCMIN", 0.0, 0.0, GF_DAY, cn)[1]
        short = Window([peri + 4 * GF_DAY, peri + 6 * GF_DAY])
        for rel in ("LOCMIN", "LOCMAX")
            @test isempty(gfdist(ctx, "MOON", "NONE", "EARTH", rel, 0.0, 0.0, 3600.0, short))
            @test isempty(GFO.gfdist("MOON", "NONE", "EARTH", rel, 0.0, 0.0, 3600.0, collect(short)))
        end
        # the whole window is the result
        @test gfdist(ctx, "MOON", "LT", "EARTH", "<", 1.0e6, 0.0, GF_DAY, cn) == cn
        @test isempty(gfdist(ctx, "MOON", "LT", "EARTH", "<", 1.0e6, 0.0, GF_DAY, Window()))
    end

    @testset "events on the end points of the confinement window" begin
        t0 = et("2022 FEB 01")
        full = Window([t0, t0 + 200 * GF_DAY])
        i = 0
        for abcorr in ("NONE", "CN+S")
            qrr = Allspice.GFRangeRate(ks, "MOON", abcorr, "EARTH")
            r = gfdist(ctx, "MOON", abcorr, "EARTH", "=", 3.85e5, 0.0, GF_DAY, full)
            p = gfdist(ctx, "MOON", abcorr, "EARTH", "LOCMIN", 0.0, 0.0, GF_DAY, full)
            a = gfdist(ctx, "MOON", abcorr, "EARTH", "LOCMAX", 0.0, 0.0, GF_DAY, full)
            @test wncard(r) >= 8 && wncard(p) >= 6 && wncard(a) >= 6
            windows = (Window([r[1], r[7]]), Window([r[3], r[3] + 10 * GF_DAY]),
                Window([r[1] - 5 * GF_DAY, r[5]]), Window([r[3], r[3]]),
                Window([r[1], r[3], r[5], r[7]]), Window([p[1], p[5]]),
                Window(sort([a[1], p[5]])), Window(sort([p[1], a[1], p[3], a[3]])),
                Window([p[3], p[3]]), Window(sort([r[1], p[3], a[3], r[9]])))
            for cn in windows, rel in GF_RELATIONS
                adjust = rel in ("ABSMIN", "ABSMAX") ? 1000.0 : 0.0
                for adj in unique((0.0, adjust))
                    gf_pick(i += 1, 4) || continue
                    w = gfdist(ctx, "MOON", abcorr, "EARTH", rel, 3.85e5, adj, GF_DAY, cn)
                    ref = GFO.gfdist("MOON", abcorr, "EARTH", rel, 3.85e5, adj, GF_DAY, collect(cn))
                    @test gf_compare("gfdist end points, " * gf_family(rel, adj), w, ref)
                    w = gfrr(ctx, "MOON", abcorr, "EARTH", rel, 0.0, adj / 1.0e5, GF_DAY, cn)
                    ref = GFO.gfrr("MOON", abcorr, "EARTH", rel, 0.0, adj / 1.0e5, GF_DAY, collect(cn))
                    @test gf_compare("gfrr end points, " * gf_family(rel, adj), w, ref;
                        value=t -> Allspice._gf_value(qrr, t), vtol=64 * eps(30.0))
                end
            end
        end
    end

    @testset "well-known events" begin
        year = Window([et("2025 JAN 01 00:00 TDB"), et("2026 JAN 01 00:00 TDB")])
        dist(t) = norm(spkpos(ctx, "MOON", t, "J2000", "NONE", "EARTH")[1])
        # Lunar perigees and apogees of 2025.
        peri = gfdist(ctx, "MOON", "NONE", "EARTH", "LOCMIN", 0.0, 0.0, GF_DAY, year)
        apo = gfdist(ctx, "MOON", "NONE", "EARTH", "LOCMAX", 0.0, 0.0, GF_DAY, year)
        @test wncard(peri) in (13, 14) && wncard(apo) in (13, 14)
        @test all(t -> 356300.0 < dist(t) < 370500.0, peri[1:2:end])
        @test all(t -> 404000.0 < dist(t) < 406750.0, apo[1:2:end])
        @test all(t -> dist(t) <= dist(t - 60.0) && dist(t) <= dist(t + 60.0), peri[1:2:end])
        @test all(t -> dist(t) >= dist(t - 60.0) && dist(t) >= dist(t + 60.0), apo[1:2:end])
        near = gfdist(ctx, "MOON", "NONE", "EARTH", "ABSMIN", 0.0, 0.0, GF_DAY, year)
        @test wncard(near) == 1 && dist(near[1]) == minimum(dist, peri[1:2:end])
        @test abs(near[1] - et("2025 NOV 05 22:30 UTC")) < 3600.0      # perigee of 5 November 2025

        # Moon and Sun seen from the Earth: smallest separation at new Moon, largest at
        # full Moon.
        sep(t) = vsep(spkpos(ctx, "MOON", t, "J2000", "LT+S", "EARTH")[1],
            spkpos(ctx, "SUN", t, "J2000", "LT+S", "EARTH")[1])
        newm = gfsep(ctx, "MOON", "POINT", "NULL", "SUN", "POINT", "NULL", "LT+S", "EARTH",
            "LOCMIN", 0.0, 0.0, GF_DAY, year)
        fullm = gfsep(ctx, "MOON", "POINT", "NULL", "SUN", "POINT", "NULL", "LT+S", "EARTH",
            "LOCMAX", 0.0, 0.0, GF_DAY, year)
        @test wncard(newm) in (12, 13) && wncard(fullm) in (12, 13)
        @test all(t -> sep(t) < deg2rad(5.4), newm[1:2:end])
        @test all(t -> sep(t) > deg2rad(174.6), fullm[1:2:end])
        @test abs(newm[1] - et("2025 JAN 29 12:36 UTC")) < 3 * 3600.0  # new Moon of 29 January 2025
        @test abs(fullm[1] - et("2025 JAN 13 22:27 UTC")) < 3 * 3600.0 # full Moon of 13 January 2025
        # The phase angle of the Moon is near 180 degrees at new Moon.
        pa = gfpa(ctx, "MOON", "SUN", "LT+S", "EARTH", "LOCMAX", 0.0, 0.0, GF_DAY, year)
        @test wncard(pa) == wncard(newm)
        @test maximum(abs.(pa[1:2:end] .- newm[1:2:end])) < 600.0

        # The range rate of the Sun relative to the Earth is zero at perihelion and at
        # aphelion, which are the extrema of the distance.
        zero_rr = gfrr(ctx, "SUN", "NONE", "EARTH", "=", 0.0, 0.0, 5 * GF_DAY, year)
        @test wncard(zero_rr) == 2
        @test abs(zero_rr[1] - et("2025 JAN 04 13:28 UTC")) < 2 * 3600.0   # perihelion
        @test abs(zero_rr[3] - et("2025 JUL 03 19:55 UTC")) < 2 * 3600.0   # aphelion
        dmin = gfdist(ctx, "SUN", "NONE", "EARTH", "LOCMIN", 0.0, 0.0, 5 * GF_DAY, year)
        dmax = gfdist(ctx, "SUN", "NONE", "EARTH", "LOCMAX", 0.0, 0.0, 5 * GF_DAY, year)
        @test wncard(dmin) == 1 && wncard(dmax) == 1
        @test abs(dmin[1] - zero_rr[1]) < 1.0 && abs(dmax[1] - zero_rr[3]) < 1.0
        au = 149597870.7
        rsun(t) = norm(spkpos(ctx, "SUN", t, "J2000", "NONE", "EARTH")[1]) / au
        @test abs(rsun(dmin[1]) - 0.98333) < 2.0e-5 && abs(rsun(dmax[1]) - 1.01664) < 2.0e-5
    end

    @testset "gfevnt" begin
        cn = Window([et("2023 JAN 01"), et("2023 JUN 01")])
        nod, noi, nol = Float64[], Int[], Bool[]
        tpar = (["TARGET", "OBSERVER", "ABCORR"], ["MOON", "EARTH", "LT+S"])
        cnames = ["TARGET", "OBSERVER", "ABCORR", "COORDINATE SYSTEM", "COORDINATE", "REFERENCE FRAME",
            "VECTOR DEFINITION", "METHOD", "DVEC", "DREF"]
        # `gfevnt` selects the quantity by its name and gives the relation to the same
        # search for every quantity, so the default run takes each relation once and the
        # quantities in turn (each quantity once at least); the exhaustive run takes the
        # product.
        nev = 0
        for (iv, (rel, wadj)) in enumerate(GF_VARIANTS), iq in 1:5
            FULL || iq == mod1(iv, 5) || continue
            nev += 1
            adjust = wadj ? 2000.0 : 0.0
            if iq == 1
                @test gfevnt(ctx, GF_DAY, gfrefn, "DISTANCE", tpar..., nod, noi, nol, rel, 3.9e5,
                    1.0e-6, adjust, cn) ==
                      gfdist(ctx, "MOON", "LT+S", "EARTH", rel, 3.9e5, adjust, GF_DAY, cn)
            elseif iq == 2
                @test gfevnt(ctx, GF_DAY, gfrefn, " range rate", tpar..., nod, noi, nol, rel, 0.01,
                    1.0e-6, adjust / 1.0e5, cn) ==
                      gfrr(ctx, "MOON", "LT+S", "EARTH", rel, 0.01, adjust / 1.0e5, GF_DAY, cn)
            elseif iq == 3
                @test gfevnt(ctx, GF_DAY, gfrefn, "PHASE ANGLE", ["ILLUM", "TARGET", "OBSERVER", "ABCORR"],
                    ["SUN", "MOON", "EARTH", "LT+S"], nod, noi, nol, rel, 1.0, 1.0e-6, adjust / 1.0e4, cn) ==
                      gfpa(ctx, "MOON", "SUN", "LT+S", "EARTH", rel, 1.0, adjust / 1.0e4, GF_DAY, cn)
            elseif iq == 4
                @test gfevnt(ctx, GF_DAY, gfrefn, "ANGULAR SEPARATION",
                    ["TARGET1", "FRAME1", "SHAPE1", "TARGET2", "FRAME2", "SHAPE2", "OBSERVER", "ABCORR"],
                    ["MOON", "IAU_MOON", "SPHERE", "SUN", "IAU_SUN", "SPHERE", "EARTH", "LT+S"],
                    nod, noi, nol, rel, 1.0, 1.0e-6, adjust / 1.0e4, cn) ==
                      gfsep(ctx, "MOON", "SPHERE", "IAU_MOON", "SUN", "SPHERE", "IAU_SUN", "LT+S", "EARTH",
                    rel, 1.0, adjust / 1.0e4, GF_DAY, cn)
            else
                @test gfevnt(ctx, GF_DAY, gfrefn, "COORDINATE", cnames,
                    ["MOON", "EARTH", "LT+S", "RA/DEC", "DECLINATION", "J2000", "POSITION", " ", " ", " "],
                    [0.0, 0.0, 0.0], noi, nol, rel, 0.1, 1.0e-6, adjust / 1.0e5, cn) ==
                      gfposc(ctx, "MOON", "J2000", "LT+S", "EARTH", "RA/DEC", "DECLINATION", rel, 0.1,
                    adjust / 1.0e5, GF_DAY, cn)
            end
        end
        @test nev == (FULL ? 5 : 1) * length(GF_VARIANTS)
        base = gfdist(ctx, "MOON", "LT+S", "EARTH", "=", 3.9e5, 0.0, GF_DAY, cn)
        # the tolerance is an argument of `gfevnt`
        @test gfevnt(ctx, GF_DAY, gfrefn, "DISTANCE", tpar..., nod, noi, nol, "=", 3.9e5, 1.0e-3,
            0.0, cn) == gfdist(ctx, "MOON", "LT+S", "EARTH", "=", 3.9e5, 0.0, GF_DAY, cn; tol=1.0e-3)
        # a step function and another refinement give the same roots within the tolerance
        varstep = gfevnt(ctx, GFConstStep(GF_DAY), gfrefn, "DISTANCE", tpar..., nod, noi, nol, "=", 3.9e5,
            1.0e-6, 0.0, cn)
        @test varstep == base
        third = gfevnt(ctx, GFConstStep(0.7 * GF_DAY), (t1, t2, s1, s2) -> t1 + (t2 - t1) / 3, "DISTANCE",
            tpar..., nod, noi, nol, "=", 3.9e5, 1.0e-6, 0.0, cn)
        @test length(third) == length(base) && maximum(abs.(collect(third) .- collect(base))) < 2.0e-6
        # the forms on a snapshot and on the default context
        @test gfevnt(ks, GF_DAY, gfrefn, "DISTANCE", tpar..., nod, noi, nol, "=", 3.9e5, 1.0e-6,
            0.0, cn) == base
        # errors
        ev(q, names, vals; op="=") = gfevnt(ctx, GF_DAY, gfrefn, q, names, vals, nod, noi, nol, op,
            1.0, 1.0e-6, 0.0, cn)
        @test throws_kind(() -> ev("SPEED", tpar...), :notrecognized)
        @test throws_kind(() -> ev("DISTANCE", tpar...; op="=="), :notrecognized)
        @test throws_kind(() -> ev("DISTANCE", ["TARGET", "OBSERVER"], ["MOON", "EARTH"]), :missingvalue)
        @test throws_kind(() -> ev("PHASE ANGLE", tpar...), :missingvalue)
        @test throws_kind(() -> ev("DISTANCE", fill("X", 11), fill("Y", 11)), :invalidcount)
        @test throws_kind(() -> gfevnt(ctx, 0.0, gfrefn, "DISTANCE", tpar..., nod, noi, nol, "=",
            1.0, 1.0e-6, 0.0, cn), :invalidstep)
        @test throws_kind(() -> gfevnt(ctx, GFConstStep(-1.0), gfrefn, "DISTANCE", tpar..., nod, noi, nol,
            "=", 1.0, 1.0e-6, 0.0, cn), :invalidstep)
        @test throws_kind(() -> gfevnt(ctx, GF_DAY, gfrefn, "DISTANCE", tpar..., nod, noi, nol, "=",
            1.0, 0.0, 0.0, cn), :invalidtolerance)
        cvals(vd) = ["MOON", "EARTH", "NONE", "LATITUDINAL", "LATITUDE", "IAU_EARTH", vd, " ", " ", " "]
        @test throws_kind(() -> ev("COORDINATE", cnames, cvals("VELOCITY")), :notsupported)
    end

    @testset "errors, with the kinds of CSPICE" begin
        cn = Window([0.0, 30 * GF_DAY])
        c = collect(cn)
        both(kind, neo, ora) = (@test gf_kind(neo) === kind; @test GFO.kind(ora) === kind)
        both(:idcodenotfound, () -> gfdist(ctx, "NOSUCHBODY", "NONE", "EARTH", "=", 1.0, 0.0, GF_DAY, cn),
            () -> GFO.gfdist("NOSUCHBODY", "NONE", "EARTH", "=", 1.0, 0.0, GF_DAY, c))
        both(:idcodenotfound, () -> gfrr(ctx, "MOON", "NONE", "NOSUCHBODY", "=", 1.0, 0.0, GF_DAY, cn),
            () -> GFO.gfrr("MOON", "NONE", "NOSUCHBODY", "=", 1.0, 0.0, GF_DAY, c))
        both(:bodiesnotdistinct, () -> gfdist(ctx, "MOON", "NONE", "301", "=", 1.0, 0.0, GF_DAY, cn),
            () -> GFO.gfdist("MOON", "NONE", "301", "=", 1.0, 0.0, GF_DAY, c))
        both(:bodiesnotdistinct, () -> gfrr(ctx, "EARTH", "NONE", "EARTH", "=", 1.0, 0.0, GF_DAY, cn),
            () -> GFO.gfrr("EARTH", "NONE", "EARTH", "=", 1.0, 0.0, GF_DAY, c))
        both(:bodiesnotdistinct, () -> gfpa(ctx, "MOON", "MOON", "NONE", "EARTH", "=", 1.0, 0.0, GF_DAY, cn),
            () -> GFO.gfpa("MOON", "MOON", "NONE", "EARTH", "=", 1.0, 0.0, GF_DAY, c))
        both(:idcodenotfound, () -> gfpa(ctx, "MOON", "LAMP", "NONE", "EARTH", "=", 1.0, 0.0, GF_DAY, cn),
            () -> GFO.gfpa("MOON", "LAMP", "NONE", "EARTH", "=", 1.0, 0.0, GF_DAY, c))
        sepn(a, sa, b, sb, o) = gfsep(ctx, a, sa, "NULL", b, sb, "NULL", "NONE", o, "=", 1.0, 0.0, GF_DAY, cn)
        sepo(a, sa, b, sb, o) = GFO.gfsep(a, sa, "NULL", b, sb, "NULL", "NONE", o, "=", 1.0, 0.0, GF_DAY, c)
        both(:bodiesnotdistinct, () -> sepn("MOON", "POINT", "MOON", "POINT", "EARTH"),
            () -> sepo("MOON", "POINT", "MOON", "POINT", "EARTH"))
        both(:bodiesnotdistinct, () -> sepn("MOON", "POINT", "SUN", "POINT", "SUN"),
            () -> sepo("MOON", "POINT", "SUN", "POINT", "SUN"))
        both(:idcodenotfound, () -> sepn("MOON", "POINT", "NOSUCHBODY", "POINT", "EARTH"),
            () -> sepo("MOON", "POINT", "NOSUCHBODY", "POINT", "EARTH"))
        both(:notrecognized, () -> sepn("MOON", "CUBE", "SUN", "POINT", "EARTH"),
            () -> sepo("MOON", "CUBE", "SUN", "POINT", "EARTH"))
        # a sphere without radii in the kernel pool
        both(:kernelvarnotfound, () -> sepn("MOON", "POINT", "MARS BARYCENTER", "SPHERE", "EARTH"),
            () -> sepo("MOON", "POINT", "MARS BARYCENTER", "SPHERE", "EARTH"))
        for (kind, (rel, adjust, step)) in ((:invalidstep, ("=", 0.0, 0.0)),
            (:notrecognized, ("=>", 0.0, GF_DAY)), (:valueoutofrange, ("ABSMIN", -1.0, GF_DAY)),
            (:invalidvalue, ("=", 1.0, GF_DAY)))
            both(kind, () -> gfdist(ctx, "MOON", "NONE", "EARTH", rel, 3.9e5, adjust, step, cn),
                () -> GFO.gfdist("MOON", "NONE", "EARTH", rel, 3.9e5, adjust, step, c))
            both(kind, () -> gfpa(ctx, "MOON", "SUN", "NONE", "EARTH", rel, 1.0, adjust, step, cn),
                () -> GFO.gfpa("MOON", "SUN", "NONE", "EARTH", rel, 1.0, adjust, step, c))
        end
        # The step is checked before the bodies, and the relation before the bodies.
        both(:invalidstep, () -> gfdist(ctx, "NOSUCHBODY", "NONE", "EARTH", "=>", 1.0, 0.0, 0.0, cn),
            () -> GFO.gfdist("NOSUCHBODY", "NONE", "EARTH", "=>", 1.0, 0.0, 0.0, c))
        both(:notrecognized, () -> gfdist(ctx, "NOSUCHBODY", "NONE", "EARTH", "=>", 1.0, 0.0, GF_DAY, cn),
            () -> GFO.gfdist("NOSUCHBODY", "NONE", "EARTH", "=>", 1.0, 0.0, GF_DAY, c))
        # An invalid correction and a body without ephemeris data have the kinds that
        # `spkezr` of this package gives them (`:badarg`, `:nopath`); CSPICE signals
        # an error too.
        @test throws_kind(() -> gfdist(ctx, "MOON", "LT+X", "EARTH", "=", 1.0, 0.0, GF_DAY, cn), :badarg)
        @test GFO.kind(() -> GFO.gfdist("MOON", "LT+X", "EARTH", "=", 1.0, 0.0, GF_DAY, c)) !== :none
        @test gf_kind(() -> gfdist(ctx, "PHOBOS", "NONE", "EARTH", "=", 1.0, 0.0, GF_DAY, cn)) in
              (:nocoverage, :nopath)
        @test GFO.kind(() -> GFO.gfdist("PHOBOS", "NONE", "EARTH", "=", 1.0, 0.0, GF_DAY, c)) !== :none
    end

    FULL && @testset "search time against CSPICE (reported, not asserted)" begin
        cn = Window([et("2020 JAN 01"), et("2030 JAN 01")])
        c = collect(cn)
        cases = (
            ("gfdist Moon-Earth, 10 years, step 1 d, '=', NONE",
                () -> gfdist(ctx, "MOON", "NONE", "EARTH", "=", 3.9e5, 0.0, GF_DAY, cn),
                () -> GFO.gfdist("MOON", "NONE", "EARTH", "=", 3.9e5, 0.0, GF_DAY, c)),
            ("gfdist Moon-Earth, 10 years, step 1 d, LOCMIN, LT+S",
                () -> gfdist(ctx, "MOON", "LT+S", "EARTH", "LOCMIN", 0.0, 0.0, GF_DAY, cn),
                () -> GFO.gfdist("MOON", "LT+S", "EARTH", "LOCMIN", 0.0, 0.0, GF_DAY, c)),
            ("gfsep Moon-Sun spheres, 10 years, step 1 d, '<', LT+S",
                () -> gfsep(ctx, "MOON", "SPHERE", "IAU_MOON", "SUN", "SPHERE", "IAU_SUN", "LT+S",
                    "EARTH", "<", 0.5, 0.0, GF_DAY, cn),
                () -> GFO.gfsep("MOON", "SPHERE", "IAU_MOON", "SUN", "SPHERE", "IAU_SUN", "LT+S",
                    "EARTH", "<", 0.5, 0.0, GF_DAY, c)),
            ("gfrr Moon-Earth, 10 years, step 1 d, '=', CN",
                () -> gfrr(ctx, "MOON", "CN", "EARTH", "=", 0.0, 0.0, GF_DAY, cn),
                () -> GFO.gfrr("MOON", "CN", "EARTH", "=", 0.0, 0.0, GF_DAY, c)),
            ("gfpa Moon-Sun-Earth, 10 years, step 1 d, ABSMAX adjust, LT+S",
                () -> gfpa(ctx, "MOON", "SUN", "LT+S", "EARTH", "ABSMAX", 0.0, 0.2, GF_DAY, cn),
                () -> GFO.gfpa("MOON", "SUN", "LT+S", "EARTH", "ABSMAX", 0.0, 0.2, GF_DAY, c)),
        )
        for (name, neo, ora) in cases
            @test gf_compare("timing", neo(), ora())
            tn, to = gf_time(neo), gf_time(ora)
            measured(@sprintf("  %-58s Allspice %8.2f ms  CSPICE %8.2f ms  ratio %.2f", name,
                1e3 * tn, 1e3 * to, to / tn))
        end
    end

    gf_report("gf")
end
