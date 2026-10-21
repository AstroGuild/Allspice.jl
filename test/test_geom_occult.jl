# Geometry group: occultation state (`occult`) against CSPICE.
#
# The output is a code, so the comparison is exact. A code changes at an instant; where
# the two libraries differ, the test accepts the difference only if CSPICE itself gives
# two codes within 1e-3 s of the epoch (the epoch is on a contact), and it counts these
# cases.
#
# Sources of the cases: scans over real events (solar eclipses of 2017-08-21 and
# 2012-05-20, the Venus transit of 2012-06-06, the lunar eclipse of 2019-01-21 seen from
# the Moon), and synthetic bodies with triaxial shapes, whose ephemerides the tests
# write with the SPK writers of types 8, 9 and 13 (see `geo_occult_kernels`).
#
# The default run has a third to a fifth of the epochs of each scan (21,500 cases and not 72,600);
# the exhaustive run (`ALLSPICE_TEST_FULL=1`) has all. The scans are not random, and
# the default run sees the same codes as the exhaustive run (the test requires it).

isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")

"Compare `occult` over `n + 1` epochs about `et0`; return `(mismatches, on a contact, codes of CSPICE)`."
function occ_scan(ctx, t1, s1, f1, t2, s2, f2, cor, obs, et0, span, n)
    bad = 0
    contact = 0
    codes = Set{Int}()
    for k in 0:n
        et = et0 + span * (k / n - 0.5)
        a = occult(ctx, t1, s1, f1, t2, s2, f2, cor, obs, et)
        b = g_occult(t1, s1, f1, t2, s2, f2, cor, obs, et)
        push!(codes, b)
        a == b && continue
        near = (g_occult(t1, s1, f1, t2, s2, f2, cor, obs, et - 1e-3), g_occult(t1, s1, f1, t2, s2, f2, cor, obs, et + 1e-3))
        if a in near && near[1] != near[2]
            contact += 1
        else
            bad += 1
            bad <= 3 && println("  occult differs: ", (t1, s1, t2, s2, cor, obs, et), " Allspice ", a, ", CSPICE ", b)
        end
    end
    return bad, contact, codes
end

occ_alloc_ee(ks, et, c) = @allocated occult(ks, 301, "ELLIPSOID", 10020, 10, "ELLIPSOID", 10010, c, 399, et)
occ_alloc_pe(ks, et, c) = @allocated occult(ks, 301, "POINT", " ", 10, "ELLIPSOID", 10010, c, 399, et)
occ_alloc_ep(ks, et, c) = @allocated occult(ks, 301, "ELLIPSOID", 10020, 10, "POINT", " ", c, 399, et)

@gf_light group("geom_occult: occult", Kernels.FULL_SET_FILES...) do
    spk, tk = geo_occult_kernels()
    paths = [geo_kernels(; dsk=true)..., spk, tk]
    oracle_load(paths)
    ctx = neo_context(paths)
    rng = MersenneTwister(SEED + 105)
    seen = Set{Int}()
    ncontact = 0

    @testset "real events: eclipses and a transit" begin
        events = (("MOON", "SUN", "EARTH", "2017-08-21T18:26:00", 6 * 3600.0),
            ("MOON", "SUN", "EARTH", "2012-05-20T23:53:00", 6 * 3600.0),
            ("VENUS", "SUN", "EARTH", "2012-06-06T01:30:00", 10 * 3600.0),
            ("EARTH", "SUN", "MOON", "2019-01-21T05:12:00", 8 * 3600.0),
            ("EARTH", "MOON", "SUN", "2019-01-21T05:12:00", 8 * 3600.0))
        nbad = 0
        ncase = 0
        nscan = nsample(50, 150)
        for cor in GEO_CORRS, (a, b, obs, utc, span) in events
            et0 = SPICE.str2et(utc)
            fa = "IAU_" * a
            fb = "IAU_" * b
            for (sa, sb) in (("ELLIPSOID", "ELLIPSOID"), ("POINT", "ELLIPSOID"), ("ELLIPSOID", "POINT")),
                (x, sx, fx, y, sy, fy) in ((a, sa, fa, b, sb, fb), (b, sb, fb, a, sa, fa))

                bad, contact, codes = occ_scan(ctx, x, sx, sx == "POINT" ? " " : fx, y, sy, sy == "POINT" ? " " : fy,
                    cor, obs, et0, span, nscan)
                nbad += bad
                ncontact += contact
                ncase += nscan + 1
                union!(seen, codes)
            end
        end
        measured("  occult, real events: $ncase cases, $nbad differences, codes seen $(sort(collect(seen)))")
        @test nbad == 0
    end

    @testset "synthetic triaxial bodies: all nine codes" begin
        nbad = 0
        ncase = 0
        codes_syn = Set{Int}()
        nscan, npoint = nsample(50, 250), nsample(25, 100)
        # Radii of the back body: larger than the front body (annular), smaller (total),
        # and of the same size (partial only).
        for radb in ([200.0, 150.0, 120.0], [9.0, 7.0, 5.0], [50.0, 40.0, 25.0], [80.0, 60.0, 50.0])
            both_pdpool(ctx, "BODY-5002_RADII", radb)
            for cor in GEO_CORRS
                for (x, fx, y, fy) in (("GEO_A", "GEO_A_FIXED", "GEO_B", "GEO_B_FIXED"), ("GEO_B", "GEO_B_FIXED", "GEO_A", "GEO_A_FIXED"))
                    bad, contact, codes = occ_scan(ctx, x, "ELLIPSOID", fx, y, "ELLIPSOID", fy, cor, "GEO_OBS", 0.0, 2 * GEO_OCC_SPAN * 0.98, nscan)
                    nbad += bad
                    ncontact += contact
                    ncase += nscan + 1
                    union!(codes_syn, codes)
                end
                bad, contact, codes = occ_scan(ctx, "GEO_A", "POINT", " ", "GEO_B", "ELLIPSOID", "GEO_B_FIXED", cor, "GEO_OBS", 0.0, 2 * GEO_OCC_SPAN * 0.98, npoint)
                nbad += bad
                ncontact += contact
                ncase += npoint + 1
                bad, contact, codes = occ_scan(ctx, "GEO_B", "POINT", " ", "GEO_A", "ELLIPSOID", "GEO_A_FIXED", cor, "GEO_OBS", 0.0, 2 * GEO_OCC_SPAN * 0.98, npoint)
                nbad += bad
                ncontact += contact
                ncase += npoint + 1
            end
        end
        both_pdpool(ctx, "BODY-5002_RADII", [80.0, 60.0, 50.0])
        union!(seen, codes_syn)
        measured("  occult, synthetic bodies: $ncase cases, $nbad differences, codes seen $(sort(collect(codes_syn)))")
        @test nbad == 0
        @test codes_syn == Set(-3:3)
        # The forms with IDs.
        @test occult(ctx, -5001, "ellipsoid", "GEO_A_FIXED", -5002, " Ellipsoid", namfrm(ctx, "GEO_B_FIXED"), Aberration("LT"), -5000, 100.0) ==
              occult(ctx, "GEO_A", "ELLIPSOID", "GEO_A_FIXED", "GEO_B", "ELLIPSOID", "GEO_B_FIXED", "LT", "GEO_OBS", 100.0)
    end

    @testset "DSK shape against a point" begin
        nbad = 0
        ncase = 0
        codes_dsk = Set{Int}()
        nscan = nsample(40, 120)
        for cor in GEO_CORRS, m in ("DSK/UNPRIORITIZED", "DSK/UNPRIORITIZED/SURFACES = 2", "dsk/unprioritized/surfaces = 1")
            for (x, sx, fx, y, sy, fy) in (("MOON", m, "IAU_MOON", "GEO_STAR", "POINT", " "), ("GEO_STAR", "point", " ", "MOON", m, "IAU_MOON"))
                bad, contact, codes = occ_scan(ctx, x, sx, fx, y, sy, fy, cor, "GEO_OBS2", 0.0, 2 * GEO_OCC_SPAN * 0.98, nscan)
                nbad += bad
                ncontact += contact
                ncase += nscan + 1
                union!(codes_dsk, codes)
            end
        end
        measured("  occult, DSK shape against a point: $ncase cases, $nbad differences, codes seen $(sort(collect(codes_dsk)))")
        @test nbad == 0
        @test codes_dsk == Set([-3, 0, 3])
    end

    @testset "all codes and contacts" begin
        measured("  occult: epochs on a contact where the code differs: $ncontact")
        @test seen == Set(-3:3)
        @test ncontact <= 5
    end

    @testset "error conditions, kinds of CSPICE" begin
        et = SPICE.str2et("2017-08-21T18:26:00")
        E = "ELLIPSOID"
        for (t1, s1, f1, t2, s2, f2, cor, obs) in (
                ("MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", "FOO", "EARTH"),
                ("MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", "S", "EARTH"),
                ("NOBODY", E, "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "EARTH"),
                ("MOON", E, "IAU_MOON", "NOBODY", E, "IAU_SUN", "LT", "EARTH"),
                ("MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "NOBODY"),
                ("MOON", E, "IAU_MOON", "MOON", E, "IAU_MOON", "LT", "EARTH"),
                ("MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "MOON"),
                ("MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "SUN"),
                ("MOON", "SPHERE", "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "EARTH"),
                ("MOON", E, "IAU_MOON", "SUN", "DISC", "IAU_SUN", "LT", "EARTH"),
                ("MOON", "POINT", " ", "SUN", "POINT", " ", "LT", "EARTH"),
                ("MOON", E, "NOFRAME", "SUN", E, "IAU_SUN", "LT", "EARTH"),
                ("MOON", E, "IAU_MOON", "SUN", E, "NOFRAME", "LT", "EARTH"),
                ("MOON", E, "IAU_EARTH", "SUN", E, "IAU_SUN", "LT", "EARTH"),
                ("MOON", E, "IAU_MOON", "SUN", E, "IAU_MOON", "LT", "EARTH"),
                ("MOON", E, " ", "SUN", E, "IAU_SUN", "LT", "EARTH"),
                ("MOON", "DSK/UNPRIORITIZED", "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "EARTH"),
                ("MOON", "DSK/UNPRIORITIZED", "IAU_MOON", "EARTH", "DSK/UNPRIORITIZED", "IAU_EARTH", "LT", "SUN"),
                ("MOON", "DSK", "IAU_MOON", "SUN", "POINT", " ", "LT", "EARTH"),
                ("MOON", "DSK/UNPRIORITIZED/SURFACES = 9", "IAU_MOON", "SUN", "POINT", " ", "LT", "EARTH"),
                ("EARTH", "DSK/UNPRIORITIZED", "IAU_EARTH", "SUN", "POINT", " ", "LT", "MOON"))
            @test same_kind(() -> g_occult(t1, s1, f1, t2, s2, f2, cor, obs, et),
                () -> occult(ctx, t1, s1, f1, t2, s2, f2, cor, obs, et))
        end
        # Radii that are not valid.
        radii = gdpool(ctx, "BODY301_RADII")
        for bad in ([1737.4, 1737.4], [1737.4, -1.0, 1737.4])
            both_pdpool(ctx, "BODY301_RADII", bad)
            @test same_kind(() -> g_occult("MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "EARTH", et),
                () -> occult(ctx, "MOON", E, "IAU_MOON", "SUN", E, "IAU_SUN", "LT", "EARTH", et))
        end
        both_pdpool(ctx, "BODY301_RADII", radii)
        # An observer inside of a target.
        both_pdpool(ctx, "BODY-5002_RADII", [3000.0, 3000.0, 3000.0])
        @test same_kind(() -> g_occult("GEO_A", E, "GEO_A_FIXED", "GEO_B", E, "GEO_B_FIXED", "NONE", "GEO_OBS", 0.0),
            () -> occult(ctx, "GEO_A", E, "GEO_A_FIXED", "GEO_B", E, "GEO_B_FIXED", "NONE", "GEO_OBS", 0.0))
        both_pdpool(ctx, "BODY-5002_RADII", [80.0, 60.0, 50.0])
    end

    @testset "no allocation and inferred type on ellipsoids and points" begin
        ks = Allspice.kernelset(ctx)
        et0 = SPICE.str2et("2017-08-21T18:26:00")
        for c in GEO_CORRS, dt in (0.0, 3000.0, 86400.0 * 40)
            @test gwarm(occ_alloc_ee, ks, et0 + dt, c) == 0
            @test gwarm(occ_alloc_pe, ks, et0 + dt, c) == 0
            @test gwarm(occ_alloc_ep, ks, et0 + dt, c) == 0
        end
        @test (@inferred occult(ks, 301, "ELLIPSOID", 10020, 10, "ELLIPSOID", 10010, "CN", 399, et0)) isa Int
    end
end
