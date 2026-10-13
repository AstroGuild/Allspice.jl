# Geometry finder, illumination angles at a surface point: `gfilum` against `gfilum_c`,
# and the quantity `ILLUMINATION ANGLE` of `gfevnt`.
#
# Three target, source and observer triples (a point of the Moon seen from the Earth
# with the Sun as source, a point of the Earth seen from the Moon with the Sun as
# source, a point of the Earth seen from the Sun with the Moon as source), two random
# surface points for each, the three angles, the seven relations and the two absolute
# extrema with an adjust value, a random correction among `NONE` and the four reception
# corrections, and confinement windows of one and of several intervals.
#
# Tolerance: the convergence tolerance of every search is the default of CSPICE,
# 1e-6 s. The number of intervals must be equal and each end point must agree within
# 2e-6 s, twice the convergence tolerance. A root or an extremum where the angle is
# almost stationary is accepted with a larger difference, up to 1 s, if the angle has
# the same value at the two epochs within 1e-12 rad.
#
# The default run compares one case in four of the product of options (`gf_pick`) and
# does not make the table of search times; the exhaustive run (`ALLSPICE_TEST_FULL=1`)
# compares every case and makes the table.
#
# Behaviour specification: test/features/gf_illum.feature

isdefined(@__MODULE__, :GFO) || include("gforacle.jl")
isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")

"The corrections that `gfilum` accepts."
const GFI_CORRS = ("NONE", "LT", "LT+S", "CN", "CN+S")
"The angles of `gfilum`."
const GFI_ANGLES = ("PHASE", "INCIDENCE", "EMISSION")

@gf_light group("gf illum: gfilum", Kernels.FULL_SET_FILES...) do
    ctx = geo_setup()
    ks = Allspice.kernelset(ctx)
    rng = MersenneTwister(SEED + 6600)
    day = 86400.0
    hour = 3600.0
    start_epoch() = -3.0e8 + 9.0e8 * rand(rng)
    function confinement(span)
        t0 = start_epoch()
        return rand(rng) < 0.7 ? Window([t0, t0 + span]) : gf_random_window(rng, 3, t0, t0 + span)
    end
    surface_point(body) = g_srfrec(GEO_IDS[body], 2pi * rand(rng) - pi, 1.2 * (rand(rng) - 0.5))

    @testset "every angle and relation against CSPICE" begin
        # (target, frame, source, observer, steps, span)
        cases = (("MOON", "IAU_MOON", "SUN", "EARTH", (day, day / 2), 60 * day),
            ("EARTH", "IAU_EARTH", "SUN", "MOON", (2 * hour, 3 * hour), 5 * day),
            ("EARTH", "IAU_EARTH", "MOON", "SUN", (3 * hour,), 4 * day))
        n = 0
        i = 0
        for (targ, fix, src, obs, steps, span) in cases, _ in 1:2, ang in GFI_ANGLES,
            (rel, wadj) in GF_VARIANTS

            gf_pick(i += 1, 4) || continue
            corr = rand(rng, GFI_CORRS)
            step = rand(rng, steps)
            sp = surface_point(targ)
            cn = confinement(span)
            q = Allspice.GFIllumination(ks, "Ellipsoid", ang, targ, src, fix, corr, obs, sp)
            refval, adjust = gf_reference(rng, q, cn, rel, wadj)
            w = gfilum(ctx, "Ellipsoid", ang, targ, src, fix, corr, obs, sp, rel, refval, adjust, step, cn)
            ref = GFO.gfilum("Ellipsoid", ang, targ, src, fix, corr, obs, sp, rel, refval, adjust, step,
                collect(cn))
            ok = gf_compare("gfilum " * lowercase(ang) * " " * gf_family(rel, adjust), w, ref;
                value=t -> Allspice._gf_value(q, t), vtol=1.0e-12)
            ok || @info "gfilum case" targ src obs ang corr sp rel refval adjust step cn
            @test ok
            @test all(t -> wnelmd(cn, t), w)
            n += 1
        end
        @test n == gf_npick(3 * 2 * 3 * 9, 4)
    end

    @testset "the value is the angle of illumg, and a root is a root" begin
        sp = g_srfrec(301, 0.3, 0.2)
        cn = Window([0.0, 60 * day])
        for (i, ang) in enumerate(GFI_ANGLES), corr in ("NONE", "CN+S")
            q = Allspice.GFIllumination(ks, "ELLIPSOID", ang, "MOON", "SUN", "IAU_MOON", corr, "EARTH", sp)
            for t in (1.0e5, 2.0e6, 4.0e6)
                angles = g_illumg("ELLIPSOID", "MOON", "SUN", t, "IAU_MOON", corr, "EARTH", sp)[3:5]
                @test Allspice._gf_value(q, t) ≈ angles[i] atol = 1.0e-11
                # the rate of `zzilusta` against a central difference of the value
                h = 30.0
                fd = (Allspice._gf_value(q, t + h) - Allspice._gf_value(q, t - h)) / (2h)
                @test Allspice._gf_ilurate(q, t) ≈ fd rtol = 1.0e-4 atol = 1.0e-11
            end
            vals = [Allspice._gf_value(q, t) for t in range(cn[1], cn[2]; length=200)]
            refval = (minimum(vals) + maximum(vals)) / 2
            w = gfilum(ctx, "ELLIPSOID", ang, "MOON", "SUN", "IAU_MOON", corr, "EARTH", sp, "=", refval,
                0.0, day, cn)
            @test wncard(w) >= 2
            for k in 1:wncard(w)
                t = w[2k-1]
                @test (Allspice._gf_value(q, t - 5.0) - refval) * (Allspice._gf_value(q, t + 5.0) - refval) < 0
            end
        end
    end

    @testset "gfevnt: ILLUMINATION ANGLE" begin
        sp = g_srfrec(301, -0.4, 0.5)
        cn = Window([0.0, 60 * day])
        names = ["TARGET", "ILLUM", "OBSERVER", "ABCORR", "REFERENCE FRAME", "ANGTYP", "METHOD", "SPOINT"]
        ev(vals, d; op="<", refval=1.2) = gfevnt(ctx, day, gfrefn, "ILLUMINATION ANGLE", names, vals, d,
            Int[], Bool[], op, refval, 1.0e-6, 0.0, cn)
        for ang in GFI_ANGLES
            base = gfilum(ctx, "ELLIPSOID", ang, "MOON", "SUN", "IAU_MOON", "LT+S", "EARTH", sp, "<", 1.2,
                0.0, day, cn)
            @test !isempty(base)
            @test ev(["MOON", "SUN", "EARTH", "LT+S", "IAU_MOON", ang, "ELLIPSOID", " "], sp) == base
            @test ev([" moon", "sun", "earth ", "lt+s", "iau_moon", lowercase(ang), "ellipsoid", " "],
                vcat(sp, 7.0)) == base
        end
        @test gfevnt(ctx, day, gfrefn, "illumination angle", lowercase.(names),
            ["MOON", "SUN", "EARTH", "NONE", "IAU_MOON", "PHASE", "ELLIPSOID", " "], sp, Int[], Bool[],
            "LOCMIN", 0.0, 1.0e-6, 0.0, cn) == gfilum(ctx, "ELLIPSOID", "PHASE", "MOON", "SUN",
            "IAU_MOON", "NONE", "EARTH", sp, "LOCMIN", 0.0, 0.0, day, cn)
        @test throws_kind(() -> ev(["MOON", "SUN", "EARTH", "LT+S", "IAU_MOON", "PHASE", "ELLIPSOID", " "],
                sp[1:2]), :invalidcount)
        @test throws_kind(() -> gfevnt(ctx, day, gfrefn, "ILLUMINATION ANGLE", names[1:7],
                ["MOON", "SUN", "EARTH", "LT+S", "IAU_MOON", "PHASE", "ELLIPSOID"], sp, Int[], Bool[], "<",
                1.2, 1.0e-6, 0.0, cn), :missingvalue)
    end

    @testset "forms, keywords and threads" begin
        sp = g_srfrec(301, 1.0, -0.3)
        cn = Window([0.0, 90 * day])
        many = gf_random_window(rng, 12, 0.0, 200 * day)
        a = ("ELLIPSOID", "INCIDENCE", "MOON", "SUN", "IAU_MOON", "LT+S", "EARTH", sp, "<", 1.4, 0.0, day)
        base = gfilum(ctx, a..., cn)
        @test wncard(base) >= 3
        @test gfilum(ks, a..., cn) == base
        @test gfilum(ks, "ELLIPSOID", "incidence ", 301, 10, namfrm(ctx, "IAU_MOON"), Aberration("LT+S"),
            399, SVector{3}(sp), a[9:end]..., cn) == base
        @test gfilum(ctx, a..., cn; threads=true) == base
        @test gfilum(ctx, a..., many; threads=true) == gfilum(ctx, a..., many)
        @test isempty(gfilum(ctx, a..., Window()))
        loose = gfilum(ctx, a..., cn; tol=1.0e-2)
        @test length(loose) == length(base) && maximum(abs.(collect(loose) .- collect(base))) < 2.0e-2
        gfstol(ctx, 1.0e-2)
        @test gfilum(ctx, a..., cn) == loose
        gfstol(ctx, 1.0e-6)
        log = String[]
        gfilum(ctx, a..., cn; report=GFLog(log), bail=GFStop(false))
        @test log == ["Illumination angle pass 1 of 2", "Illumination angle pass 2 of 2"]
        @test_throws InterruptException gfilum(ctx, a..., cn; report=GFLog(String[]), bail=GFStop(true))
        # many tasks, one context
        serial = [gfilum(ctx, "ELLIPSOID", ang, "MOON", "SUN", "IAU_MOON", corr, "EARTH", sp, "=", 1.0,
            0.0, day, cn) for ang in GFI_ANGLES, corr in GFI_CORRS]
        tasks = [Threads.@spawn gfilum(ctx, "ELLIPSOID", ang, "MOON", "SUN", "IAU_MOON", corr, "EARTH",
            sp, "=", 1.0, 0.0, day, cn) for ang in GFI_ANGLES, corr in GFI_CORRS]
        @test all(fetch(t) == s for (t, s) in zip(tasks, serial))
    end

    @testset "errors, with the kinds of CSPICE" begin
        sp = g_srfrec(301, 0.3, 0.2)
        cn = Window([0.0, 30 * day])
        c = collect(cn)
        both(kind, neo, ora) = (@test gf_kind(neo) === kind; @test GFO.kind(ora) === kind)
        # method, angtyp, target, illmn, fixref, abcorr, obsrvr, spoint, relate, refval,
        # adjust, step
        good = ("ELLIPSOID", "PHASE", "MOON", "SUN", "IAU_MOON", "LT+S", "EARTH", sp, "<", 1.2, 0.0, day)
        with(i, v) = ntuple(k -> k == i ? v : good[k], length(good))
        both(:none, () -> gfilum(ctx, good..., cn), () -> GFO.gfilum(good..., c))
        for (kind, args) in (
            (:invalidmethod, with(1, "DSK/UNPRIORITIZED")),
            (:invalidmethod, with(1, "Near point: ellipsoid")),
            (:notsupported, with(2, "SOLAR INCIDENCE")),
            (:idcodenotfound, with(3, "NOSUCHBODY")),
            (:idcodenotfound, with(4, "NOSUCHBODY")),
            (:idcodenotfound, with(7, "NOSUCHBODY")),
            (:bodiesnotdistinct, with(7, "MOON")),
            (:bodiesnotdistinct, with(4, "MOON")),
            (:unknownframe, with(5, "NOSUCHFRAME")),
            (:invalidframe, with(5, "IAU_EARTH")),
            (:invalidoption, with(6, "LT+X")),
            (:notsupported, with(6, "XLT")),
            (:notsupported, with(6, "XCN+S")),
            (:notrecognized, with(9, "=<")),
            (:valueoutofrange, with(11, -1.0)),
            (:invalidvalue, with(11, 1.0)),
            (:invalidstep, with(12, 0.0)),
            (:invalidstep, with(12, -day)))
            both(kind, () -> gfilum(ctx, args..., cn), () -> GFO.gfilum(args..., c))
        end
        # the source can be the observer; case and blanks of the angle and of the method
        both(:none, () -> gfilum(ctx, with(4, "EARTH")..., cn), () -> GFO.gfilum(with(4, "EARTH")..., c))
        both(:none, () -> gfilum(ctx, " ellipsoid", " emis sion", "moon", "sun", "iau_moon", "cn + s",
                "earth", sp, "locmax", 0.0, 0.0, day, cn),
            () -> GFO.gfilum(" ellipsoid", " emis sion", "moon", "sun", "iau_moon", "cn + s", "earth", sp,
                "locmax", 0.0, 0.0, day, c))
        @test gfilum(ctx, with(4, "EARTH")..., cn) == Window(GFO.gfilum(with(4, "EARTH")..., c))
        # a target without radii
        both(:kernelvarnotfound, () -> gfilum(ctx, "ELLIPSOID", "PHASE", "EARTH BARYCENTER", "SUN",
                "J2000", "NONE", "MOON", sp, "<", 1.2, 0.0, day, cn),
            () -> GFO.gfilum("ELLIPSOID", "PHASE", "EARTH BARYCENTER", "SUN", "J2000", "NONE", "MOON", sp,
                "<", 1.2, 0.0, day, c))
        @test throws_kind(() -> gfilum(ctx, with(8, [1.0, 2.0])..., cn), :badarg)
    end

    FULL && @testset "search time against CSPICE (reported, not asserted)" begin
        t0 = SPICE.str2et("2020 JAN 01")
        sp = g_srfrec(301, 0.3, 0.2)
        spe = g_srfrec(399, -1.2, 0.6)
        cases2 = (
            ("gfilum incidence at a point of the Moon, 1 year, '<', LT+S", 365 * day,
                ("ELLIPSOID", "INCIDENCE", "MOON", "SUN", "IAU_MOON", "LT+S", "EARTH", sp, "<", 1.4, 0.0, day)),
            ("gfilum phase at a point of the Moon, 1 year, LOCMIN, NONE", 365 * day,
                ("ELLIPSOID", "PHASE", "MOON", "SUN", "IAU_MOON", "NONE", "EARTH", sp, "LOCMIN", 0.0, 0.0, day)),
            ("gfilum emission at a point of the Earth, 30 days, '=', CN+S", 30 * day,
                ("ELLIPSOID", "EMISSION", "EARTH", "SUN", "IAU_EARTH", "CN+S", "MOON", spe, "=", 1.0, 0.0,
                    3 * hour)),
        )
        for (name, span, a) in cases2
            cn = Window([t0, t0 + span])
            c = collect(cn)
            @test gf_compare("timing", gfilum(ctx, a..., cn), GFO.gfilum(a..., c))
            tn, to = gf_time(() -> gfilum(ctx, a..., cn)), gf_time(() -> GFO.gfilum(a..., c))
            measured(@sprintf("  %-62s Allspice %8.2f ms  CSPICE %8.2f ms  ratio %.2f", name,
                1e3 * tn, 1e3 * to, to / tn))
        end
    end

    gf_report("gfilum")
end
