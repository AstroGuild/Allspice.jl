# Geometry finder, occultations: `gfoclt` against `gfoclt_c`, and `gfocce`.
#
# Real bodies: the Moon in front of the Sun as seen from the Earth, the Earth in front
# of the Sun as seen from the Moon, the Earth in front of the Moon as seen from the
# Sun, the transits of Mercury and Venus, and Venus behind the Sun; the four occultation
# types, ellipsoids and points, the five corrections that have an effect (`NONE`, `LT`,
# `CN`, `XLT`, `XCN`) and confinement windows of one and of several intervals.
# Synthetic bodies (two triaxial ellipsoids that cross in about an hour, with radii that
# give each occultation type) and a DSK shape of the Moon in front of a point.
#
# Tolerance: the convergence tolerance of every search is the default of CSPICE,
# 1e-6 s. The number of intervals must be equal and each end point must agree within
# 2e-6 s, twice the convergence tolerance: the two libraries make the same bisection of
# a state that is true or false, so the brackets differ only if the state differs at a
# sample, which is possible only within rounding of a contact.
#
# Behaviour specification: test/features/gf_occult.feature

isdefined(@__MODULE__, :GFO) || include("gforacle.jl")
isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")

"The corrections of `gfoclt` that have an effect (stellar aberration is not used)."
const GFOC_CORRS = ("NONE", "LT", "CN", "XLT", "XCN")
const GFOC_TYPES = ("ANY", "PARTIAL", "FULL", "ANNULAR")

group("gf occult: gfoclt, gfocce", Kernels.FULL_SET_FILES...) do
    paths = [geo_kernels(; dsk=true)..., geo_occult_kernels()...]
    oracle_load(paths)
    ctx = neo_context(paths)
    ks = NeoSpice.kernelset(ctx)
    rng = MersenneTwister(SEED + 6700)
    day = 86400.0
    hour = 3600.0
    year = 365.25 * day
    E = "ELLIPSOID"
    P = "POINT"
    found = Dict{String,Int}()
    "Compare one search with CSPICE; `a` holds the arguments before the window."
    function check(label, a, cn; atol=2.0e-6)
        w = gfoclt(ctx, a..., cn)
        ref = GFO.gfoclt(a..., collect(cn))
        ok = gf_compare("gfoclt " * label, w, ref; atol)
        ok || @info "gfoclt case" label a cn
        @test ok
        @test all(t -> wnelmd(cn, t), w)
        found[label] = get(found, label, 0) + wncard(w)
        return w
    end

    @testset "eclipses: the Sun, the Earth and the Moon" begin
        # (label, front, front frame, back, back frame, observer, steps)
        cases = (("Moon in front of the Sun, from the Earth", "MOON", "IAU_MOON", "SUN", "IAU_SUN", "EARTH",
                (2 * hour, 3 * hour)),
            ("Earth in front of the Sun, from the Moon", "EARTH", "IAU_EARTH", "SUN", "IAU_SUN", "MOON",
                (hour, 2 * hour)),
            ("Earth in front of the Moon, from the Sun", "EARTH", "IAU_EARTH", "MOON", "IAU_MOON", "SUN",
                (hour, 2 * hour)))
        for (label, f, ff, b, bf, obs, steps) in cases, typ in GFOC_TYPES, corr in GFOC_CORRS
            t0 = -3.0e8 + 8.0e8 * rand(rng)
            cn = rand(rng) < 0.7 ? Window([t0, t0 + 2 * year]) : gf_random_window(rng, 4, t0, t0 + 3 * year)
            check(label * ", " * typ, (typ, f, E, ff, b, E, bf, corr, obs, rand(rng, steps)), cn)
        end
        # points: the type is ANY
        for corr in GFOC_CORRS
            t0 = -3.0e8 + 4.0e8 * rand(rng)
            cn = Window([t0, t0 + 12 * year])
            check("Moon in front of the Sun as a point", ("ANY", "MOON", E, "IAU_MOON", "SUN", P, " ", corr,
                    "EARTH", 3 * hour), cn)
            check("Moon as a point in front of the Sun", ("ANY", "MOON", P, " ", "SUN", E, "IAU_SUN", corr,
                    "EARTH", 3 * hour), cn)
            check("Earth in front of the Moon as a point", ("any", "earth", "ellipsoid", "iau_earth", "moon",
                    "point", "", corr, "sun", 2 * hour), cn)
        end
        measured(@sprintf("  gfoclt eclipses: intervals found: %s", join(["$k: $v" for (k, v) in
            sort!(collect(found); by=first)], "; ")))
        @test all(v > 0 for (k, v) in found if endswith(k, "ANY") || endswith(k, "PARTIAL") ||
                                               occursin("point", k))
        @test found["Earth in front of the Sun, from the Moon, FULL"] > 0
        @test found["Earth in front of the Moon, from the Sun, FULL"] > 0
    end

    @testset "transits of Mercury and Venus, and Venus behind the Sun" begin
        empty!(found)
        transits = (("MERCURY", "IAU_MERCURY", "2003 MAY 07"), ("MERCURY", "IAU_MERCURY", "2006 NOV 08"),
            ("MERCURY", "IAU_MERCURY", "2016 MAY 09"), ("MERCURY", "IAU_MERCURY", "2019 NOV 11"),
            ("VENUS", "IAU_VENUS", "2004 JUN 08"), ("VENUS", "IAU_VENUS", "2012 JUN 06"))
        for (planet, pf, date) in transits, corr in GFOC_CORRS
            tm = SPICE.str2et(date)
            cn = Window([tm - 4 * day, tm + 4 * day])
            for typ in GFOC_TYPES
                # the partial phases last some minutes
                step = typ == "PARTIAL" ? 60.0 : hour
                check("transit, " * typ, (typ, planet, E, pf, "SUN", E, "IAU_SUN", corr, "EARTH", step), cn)
            end
            check("transit, planet as a point", ("ANY", planet, P, " ", "SUN", E, "IAU_SUN", corr, "EARTH",
                    hour), cn)
        end
        @test found["transit, ANY"] == 30 && found["transit, ANNULAR"] == 30
        measured(@sprintf("  gfoclt transits: intervals found: %s", join(["$k: $v" for (k, v) in
            sort!(collect(found); by=first)], "; ")))
        @test found["transit, PARTIAL"] == 60 && found["transit, FULL"] == 0
        @test found["transit, planet as a point"] == 30
        # Venus behind the Sun: about once in 8 years
        for corr in GFOC_CORRS
            t0 = SPICE.str2et("2000 JAN 01") + 3.0e7 * rand(rng)
            cn = Window([t0, t0 + 10 * year])
            check("Venus behind the Sun, FULL", ("FULL", "SUN", E, "IAU_SUN", "VENUS", E, "IAU_VENUS", corr,
                    "EARTH", 6 * hour), cn)
            check("Venus as a point behind the Sun", ("ANY", "SUN", E, "IAU_SUN", "VENUS", P, " ", corr,
                    "EARTH", 6 * hour), cn)
        end
        @test found["Venus behind the Sun, FULL"] >= 5 && found["Venus as a point behind the Sun"] >= 5
    end

    @testset "synthetic ellipsoids: every type has events" begin
        empty!(found)
        T = GEO_OCC_SPAN - 100.0
        cn = Window([-T, T])
        for radb in ([80.0, 60.0, 50.0], [8.0, 6.0, 5.0], [25.0, 22.0, 12.0])
            both_pdpool(ctx, "BODY-5002_RADII", radb)
            for typ in GFOC_TYPES, corr in GFOC_CORRS, (f, ff, b, bf) in (("GEO_A", "GEO_A_FIXED", "GEO_B",
                    "GEO_B_FIXED"), ("GEO_B", "GEO_B_FIXED", "GEO_A", "GEO_A_FIXED"))
                front = f == "GEO_A" ? "A in front" : "B in front"
                w = check("synthetic, $front, $typ", (typ, f, E, ff, b, E, bf, corr, "GEO_OBS", 20.0), cn)
                # B is behind A as seen from GEO_OBS: no event with B in front
                f == "GEO_B" && @test isempty(w)
            end
            for corr in GFOC_CORRS
                check("synthetic, B as a point", ("ANY", "GEO_A", E, "GEO_A_FIXED", "GEO_B", P, " ", corr,
                        "GEO_OBS", 20.0), cn)
            end
        end
        both_pdpool(ctx, "BODY-5002_RADII", [80.0, 60.0, 50.0])
        for typ in GFOC_TYPES
            @test found["synthetic, A in front, $typ"] > 0
        end
        @test found["synthetic, B as a point"] == 15
        measured(@sprintf("  gfoclt synthetic: intervals found: %s", join(["$k: $v" for (k, v) in
            sort!(collect(found); by=first) if v > 0], "; ")))
    end

    @testset "a DSK shape and a point" begin
        empty!(found)
        T = GEO_OCC_SPAN - 100.0
        cn = Window([-T, T])
        for corr in GFOC_CORRS, m in ("DSK/UNPRIORITIZED", "DSK/UNPRIORITIZED/SURFACES = 2",
                "dsk/unprioritized/surfaces = 1", E)
            check("DSK Moon in front of a point", ("ANY", "MOON", m, "IAU_MOON", "GEO_STAR", P, " ", corr,
                    "GEO_OBS2", 30.0), cn)
            w = check("point in front of the DSK Moon", ("ANY", "GEO_STAR", "point", " ", "MOON", m,
                    "IAU_MOON", corr, "GEO_OBS2", 30.0), cn)
            @test isempty(w)
        end
        @test found["DSK Moon in front of a point"] == 20
    end

    @testset "gfocce: step function and refinement" begin
        t0 = SPICE.str2et("2010 JAN 01")
        cn = Window([t0, t0 + 12 * year])
        a = ("ANY", "MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "EARTH")
        base = gfoclt(ctx, a..., 3 * hour, cn)
        @test wncard(base) >= 4
        @test gfocce(ctx, a..., 1.0e-6, 3 * hour, gfrefn, cn) == base
        @test gfocce(ks, a..., 1.0e-6, et -> 3 * hour, gfrefn, cn) == base
        @test gfocce(ctx, a..., 1.0e-6, 3 * hour, gfrefn, cn; threads=true) == base
        # another refinement gives the same contacts within the tolerance
        quarter(t1, t2, s1, s2) = t1 + 0.25 * (t2 - t1)
        other = gfocce(ctx, a..., 1.0e-6, 3 * hour, quarter, cn)
        @test length(other) == length(base) && maximum(abs.(collect(other) .- collect(base))) <= 2.0e-6
        # a looser tolerance
        loose = gfocce(ctx, a..., 1.0e-1, 3 * hour, gfrefn, cn)
        @test length(loose) == length(base) && 1.0e-6 < maximum(abs.(collect(loose) .- collect(base))) <= 0.2
        @test gfoclt(ctx, a..., 3 * hour, cn; tol=1.0e-1) == loose
        @test throws_kind(() -> gfocce(ctx, a..., 0.0, 3 * hour, gfrefn, cn), :invalidtolerance)
        @test throws_kind(() -> gfocce(ctx, a..., 1.0e-6, et -> 0.0, gfrefn, cn), :invalidstep)
        log = String[]
        gfocce(ctx, a..., 1.0e-6, 3 * hour, gfrefn, cn; report=(e, x...) -> e === :init && push!(log, x[2]))
        @test log == ["Occultation/transit search "]
        @test_throws InterruptException gfocce(ctx, a..., 1.0e-6, 3 * hour, gfrefn, cn; bail=() -> true)
    end

    @testset "forms and threads" begin
        t0 = SPICE.str2et("2010 JAN 01")
        cn = Window([t0, t0 + 3 * year])
        many = gf_random_window(rng, 10, t0, t0 + 6 * year)
        a = ("PARTIAL", "EARTH", E, "IAU_EARTH", "SUN", E, "IAU_SUN", "CN", "MOON", 2 * hour)
        base = gfoclt(ctx, a..., cn)
        @test wncard(base) >= 4
        @test gfoclt(ks, a..., cn) == base
        @test gfoclt(ks, "partial ", 399, E, namfrm(ctx, "IAU_EARTH"), 10, E, namfrm(ctx, "IAU_SUN"),
            Aberration("CN"), 301, 2 * hour, cn) == base
        @test gfoclt(ctx, a..., cn; threads=true) == base
        @test gfoclt(ctx, a..., many; threads=true) == gfoclt(ctx, a..., many)
        @test isempty(gfoclt(ctx, a..., Window()))
        gfstol(ctx, 1.0e-1)
        @test gfoclt(ctx, a..., cn) == gfoclt(ks, a..., cn; tol=1.0e-1) != base
        gfstol(ctx, 1.0e-6)
        serial = [gfoclt(ctx, typ, "MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", corr, "EARTH", 3 * hour, cn)
                  for typ in GFOC_TYPES, corr in GFOC_CORRS]
        tasks = [Threads.@spawn gfoclt(ctx, typ, "MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", corr, "EARTH",
            3 * hour, cn) for typ in GFOC_TYPES, corr in GFOC_CORRS]
        @test all(fetch(t) == s for (t, s) in zip(tasks, serial))
    end

    @testset "errors, with the kinds of CSPICE" begin
        t0 = SPICE.str2et("2019 JUN 20")
        cn = Window([t0, t0 + 30 * day])
        c = collect(cn)
        both(kind, neo, ora) = (@test gf_kind(neo) === kind; @test GFO.kind(ora) === kind)
        # occtyp, front, fshape, fframe, back, bshape, bframe, abcorr, obsrvr, step
        good = ("ANY", "MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "EARTH", 3 * hour)
        with(p...) = (d = Dict(p); ntuple(k -> get(d, k, good[k]), length(good)))
        both(:none, () -> gfoclt(ctx, good..., cn), () -> GFO.gfoclt(good..., c))
        for (kind, args) in (
            (:invalidocctype, with(1 => "TOTAL")),
            (:idcodenotfound, with(2 => "NOSUCHBODY")),
            (:idcodenotfound, with(5 => "NOSUCHBODY")),
            (:idcodenotfound, with(9 => "NOSUCHBODY")),
            (:bodiesnotdistinct, with(5 => "MOON")),
            (:bodiesnotdistinct, with(9 => "MOON")),
            (:bodiesnotdistinct, with(9 => "SUN")),
            (:badmethodsyntax, with(3 => "SPHERE")),
            (:badmethodsyntax, with(6 => "SPHERE")),
            (:invalidshapecombo, with(3 => P, 6 => P)),
            (:badtypeshapecombo, with(1 => "FULL", 6 => P)),
            (:badtypeshapecombo, with(1 => "PARTIAL", 3 => P)),
            (:invalidshapecombo, with(3 => "DSK/UNPRIORITIZED")),
            (:invalidframe, with(4 => "IAU_SUN")),
            (:invalidframe, with(7 => "IAU_MOON")),
            (:invalidframe, with(4 => "NOSUCHFRAME")),
            (:invalidframe, with(4 => " ")),
            (:invalidoption, with(8 => "LT+X")),
            (:invalidstep, with(10 => 0.0)),
            (:invalidstep, with(10 => -1.0)))
            both(kind, () -> gfoclt(ctx, args..., cn), () -> GFO.gfoclt(args..., c))
        end
        # stellar aberration in the correction: accepted and not used
        for corr in ("LT+S", "CN+S", "XLT+S", "XCN+S")
            plain = replace(corr, "+S" => "")
            @test gfoclt(ctx, with(8 => corr)..., cn) == gfoclt(ctx, with(8 => plain)..., cn)
            @test gf_compare("gfoclt correction with +S", gfoclt(ctx, with(8 => corr)..., cn),
                GFO.gfoclt(with(8 => corr)..., c))
        end
        # the frame of a point is not read
        both(:none, () -> gfoclt(ctx, with(6 => P, 7 => "NOSUCHFRAME")..., cn),
            () -> GFO.gfoclt(with(6 => P, 7 => "NOSUCHFRAME")..., c))
        # a body without radii as an ellipsoid
        both(:kernelvarnotfound, () -> gfoclt(ctx, with(2 => "GEO_STAR", 4 => "J2000")..., cn),
            () -> GFO.gfoclt(with(2 => "GEO_STAR", 4 => "J2000")..., c))
    end

    @testset "search time against CSPICE (reported, not asserted)" begin
        t0 = SPICE.str2et("2020 JAN 01")
        cases2 = (
            ("gfoclt Moon in front of the Sun, ANY, 12 years, LT", 12 * year,
                ("ANY", "MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "EARTH", 3 * hour)),
            ("gfoclt Earth in front of the Sun from the Moon, FULL, 3 years, CN", 3 * year,
                ("FULL", "EARTH", E, "IAU_EARTH", "SUN", E, "IAU_SUN", "CN", "MOON", 2 * hour)),
            ("gfoclt Moon as a point in front of the Sun, 12 years, NONE", 12 * year,
                ("ANY", "MOON", P, " ", "SUN", E, "IAU_SUN", "NONE", "EARTH", 3 * hour)),
        )
        for (name, span, a) in cases2
            cn = Window([t0, t0 + span])
            c = collect(cn)
            @test gf_compare("timing", gfoclt(ctx, a..., cn), GFO.gfoclt(a..., c))
            tn, to = gf_time(() -> gfoclt(ctx, a..., cn)), gf_time(() -> GFO.gfoclt(a..., c))
            measured(@sprintf("  %-66s NeoSpice %8.2f ms  CSPICE %8.2f ms  ratio %.2f", name,
                1e3 * tn, 1e3 * to, to / tn))
        end
    end

    n = sum(v for (k, v) in GF_COUNT if startswith(k, "gfoclt"); init=0)
    worst = maximum(v for (k, v) in GF_MAXDIFF if startswith(k, "gfoclt"); init=0.0)
    measured(@sprintf("  gfoclt: %d searches compared with CSPICE, largest end-point difference %.3e s", n, worst))
end
