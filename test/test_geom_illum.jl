# Geometry group: illumination angles, the phase angle and the angular separation of two
# bodies (`ilumin`, `illumg`, `illumf`, `phaseq`, `trgsep`) against CSPICE.
#
# Feature file: features/geometry_illumination.feature
#
# Error measures: the absolute difference of the angles in radians, the norm of the
# difference of the observer-to-point vector over its norm, and the absolute difference
# of the epoch in seconds.
#
# Conditioning. An angle near 0 or pi comes from `vsep`, which loses half of the digits
# there: a relative rounding of 1e-16 in the vectors gives 1e-8 rad at the pole. The
# sample excludes no case; the tolerance on the angles covers it.

isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")
import ForwardDiff

# Measured: 9e-16 rad (angles), 4e-16 (srfvec), 0 (epochs); phaseq 4e-16 rad; trgsep 2.2e-16 rad.
const ILL_TOL = (angle=1e-10, vec=5e-13, epoch=1e-7)

ill_alloc_min(ks, et, c, sp) = @allocated ilumin(ks, "ELLIPSOID", 301, et, 10020, c, 399, sp)
ill_alloc_g(ks, et, c, sp) = @allocated illumg(ks, "ELLIPSOID", 301, 399, et, 10020, c, 10, sp)
ill_alloc_f(ks, et, c, sp) = @allocated illumf(ks, "ELLIPSOID", 301, 10, et, 10020, c, 399, sp)
ill_alloc_ph(ks, et, c) = @allocated phaseq(ks, et, 301, 10, 399, c)
ill_alloc_sep(ks, et, c) = @allocated trgsep(ks, et, 301, "SPHERE", 10020, 10, "POINT", 10010, 399, c)

@gf_light group("geom_illum: ilumin, illumg, illumf, phaseq, trgsep", Kernels.FULL_SET_FILES...) do
    ctx = geo_setup(; dsk=true)
    ks = Allspice.kernelset(ctx)
    rng = MersenneTwister(SEED + 102)
    reception = GEO_CORRS[1:5]

    @testset "ellipsoids, 600 random cases" begin
        st = GeoStats()
        nvis = 0
        nlit = 0
        for k in 1:600
            t = rand(rng, GEO_BODIES)
            o = rand(rng, setdiff(GEO_BODIES, [t]))
            fix = rand(rng, GEO_FIXED[t])
            cor = rand(rng, GEO_CORRS)
            et = geo_epoch(rng)
            sp = g_srfrec(GEO_IDS[t], 2pi * (rand(rng) - 0.5), pi * (rand(rng) - 0.5))
            m = rand(rng, ["ELLIPSOID", "Ellipsoid", " ellipsoid "])
            if t != "SUN"
                a = ilumin(ctx, m, t, et, fix, cor, o, sp)
                b = g_ilumin(m, t, et, fix, cor, o, sp)
                gput!(st, "ilumin trgepc (s)", abs(a[1] - b[1]))
                gput!(st, "ilumin srfvec", grel(a[2], b[2]))
                gput!(st, "ilumin angles (rad)", maximum(abs.(collect(a[3:5]) .- collect(b[3:5]))))
            end
            src = rand(rng, setdiff(GEO_BODIES, [t]))
            a = illumg(ctx, m, t, src, et, fix, cor, o, sp)
            b = g_illumg(m, t, src, et, fix, cor, o, sp)
            gput!(st, "illumg trgepc (s)", abs(a[1] - b[1]))
            gput!(st, "illumg srfvec", grel(a[2], b[2]))
            gput!(st, "illumg angles (rad)", maximum(abs.(collect(a[3:5]) .- collect(b[3:5]))))
            a = illumf(ctx, m, t, src, et, fix, cor, o, sp)
            b = g_illumf(m, t, src, et, fix, cor, o, sp)
            gput!(st, "illumf trgepc (s)", abs(a[1] - b[1]))
            gput!(st, "illumf srfvec", grel(a[2], b[2]))
            gput!(st, "illumf angles (rad)", maximum(abs.(collect(a[3:5]) .- collect(b[3:5]))))
            # A flag can differ only where its angle is at 90 degrees.
            a[6] == b[6] || gput!(st, "illumf visibl differs: |emissn - pi/2|", abs(b[5] - pi / 2))
            a[7] == b[7] || gput!(st, "illumf lit differs: |incdnc - pi/2|", abs(b[4] - pi / 2))
            nvis += b[6]
            nlit += b[7]
            if k % 10 == 0
                c = illumf(ctx, m, GEO_IDS[t], GEO_IDS[src], et, namfrm(ctx, fix), Aberration(cor), GEO_IDS[o], SVector{3}(sp))
                @test c == a
            end
        end
        greport("ellipsoid", st)
        @test 150 < nvis < 450 && 150 < nlit < 450
        for name in ("ilumin", "illumg", "illumf")
            @test gget(st, "$name angles (rad)") < ILL_TOL.angle
            @test gget(st, "$name srfvec") < ILL_TOL.vec
            @test gget(st, "$name trgepc (s)") < ILL_TOL.epoch
        end
        @test !haskey(st.d, "illumf visibl differs: |emissn - pi/2|")
        @test !haskey(st.d, "illumf lit differs: |incdnc - pi/2|")
    end

    @testset "DSK shapes: Moon, 200 random cases" begin
        st = GeoStats()
        nvis = 0
        nlit = 0
        for k in 1:200
            o = rand(rng, ["EARTH", "SUN", "VENUS"])
            cor = rand(rng, GEO_CORRS)
            et = geo_epoch(rng)
            surf = rand(rng, 1:2)
            m = rand(rng, ["DSK/UNPRIORITIZED/SURFACES = $surf", "dsk/unprioritized/surfaces=$surf"])
            sp = c_latsrf(m, "MOON", et, "IAU_MOON", [(2pi * (rand(rng) - 0.5), pi * (rand(rng) - 0.5))])[1]
            a = ilumin(ctx, m, "MOON", et, "IAU_MOON", cor, o, sp)
            b = g_ilumin(m, "MOON", et, "IAU_MOON", cor, o, sp)
            gput!(st, "ilumin trgepc (s)", abs(a[1] - b[1]))
            gput!(st, "ilumin srfvec", grel(a[2], b[2]))
            gput!(st, "ilumin angles (rad)", maximum(abs.(collect(a[3:5]) .- collect(b[3:5]))))
            src = rand(rng, ["SUN", "EARTH"])
            src == o && (src = "SUN"; o = "EARTH")
            a = illumg(ctx, m, "MOON", src, et, "IAU_MOON", cor, o, sp)
            b = g_illumg(m, "MOON", src, et, "IAU_MOON", cor, o, sp)
            gput!(st, "illumg angles (rad)", maximum(abs.(collect(a[3:5]) .- collect(b[3:5]))))
            # `illumf` of CSPICE ignores the surface list for the normal; the DSK file of
            # the test has two surfaces, so the unlisted surface can be the one that is used.
            a = illumf(ctx, m, "MOON", src, et, "IAU_MOON", cor, o, sp)
            b = g_illumf(m, "MOON", src, et, "IAU_MOON", cor, o, sp)
            gput!(st, "illumf srfvec", grel(a[2], b[2]))
            gput!(st, "illumf angles (rad)", maximum(abs.(collect(a[3:5]) .- collect(b[3:5]))))
            (a[6] == b[6] && a[7] == b[7]) || gput!(st, "illumf flags differ", 1.0)
            nvis += b[6]
            nlit += b[7]
        end
        greport("DSK", st)
        @test nvis > 30 && nlit > 30
        for name in ("ilumin", "illumg", "illumf")
            @test gget(st, "$name angles (rad)") < ILL_TOL.angle
        end
        @test gget(st, "ilumin srfvec") < ILL_TOL.vec
        @test gget(st, "illumf srfvec") < ILL_TOL.vec
        @test !haskey(st.d, "illumf flags differ")
    end

    @testset "phaseq, 600 random cases" begin
        worst = 0.0
        for k in 1:600
            t, src, o = shuffle(rng, GEO_BODIES)[1:3]
            cor = rand(rng, reception)
            et = geo_epoch(rng)
            a = phaseq(ctx, et, t, src, o, cor)
            worst = max(worst, abs(a - g_phaseq(et, t, src, o, cor)))
            k % 10 == 0 && @test phaseq(ctx, et, GEO_IDS[t], GEO_IDS[src], GEO_IDS[o], Aberration(cor)) == a
        end
        measured(@sprintf("  %-64s max=%.3e", "phaseq (rad)", worst))
        @test worst < ILL_TOL.angle
    end

    @testset "trgsep, 300 random cases" begin
        worst = 0.0
        frames = Dict("SUN" => "IAU_SUN", "MERCURY" => "IAU_MERCURY", "VENUS" => "IAU_VENUS",
            "EARTH" => "IAU_EARTH", "MOON" => "IAU_MOON")
        for k in 1:300
            t1, t2, o = shuffle(rng, GEO_BODIES)[1:3]
            s1, s2 = rand(rng, ("POINT", "SPHERE")), rand(rng, ("POINT", "SPHERE"))
            cor = rand(rng, GEO_CORRS)
            et = geo_epoch(rng)
            a = trgsep(ctx, et, t1, s1, frames[t1], t2, s2, frames[t2], o, cor)
            worst = max(worst, abs(a - g_trgsep(et, t1, s1, frames[t1], t2, s2, frames[t2], o, cor)))
            k % 10 == 0 && @test trgsep(ctx, et, GEO_IDS[t1], s1, frames[t1], GEO_IDS[t2], s2, frames[t2],
                GEO_IDS[o], Aberration(cor)) == a
        end
        measured(@sprintf("  %-64s max=%.3e", "trgsep (rad)", worst))
        @test worst < ILL_TOL.angle
    end

    @testset "error conditions, kinds of CSPICE" begin
        et = 1.0e8
        sp = g_srfrec(301, 0.3, 0.2)
        for (m, t, src, fix, cor, o, p) in (
                ("ELLIPSOID", "MOON", "SUN", "IAU_MOON", "FOO", "EARTH", sp),
                ("ELLIPSOID", "NOBODY", "SUN", "IAU_MOON", "LT", "EARTH", sp),
                ("ELLIPSOID", "MOON", "NOBODY", "IAU_MOON", "LT", "EARTH", sp),
                ("ELLIPSOID", "MOON", "SUN", "IAU_MOON", "LT", "NOBODY", sp),
                ("ELLIPSOID", "MOON", "SUN", "IAU_MOON", "LT", "MOON", sp),
                ("ELLIPSOID", "MOON", "SUN", "NOFRAME", "LT", "EARTH", sp),
                ("ELLIPSOID", "MOON", "SUN", "IAU_EARTH", "LT", "EARTH", sp),
                ("SPHERE", "MOON", "SUN", "IAU_MOON", "LT", "EARTH", sp),
                ("NEAR POINT/ELLIPSOID", "MOON", "SUN", "IAU_MOON", "LT", "EARTH", sp),
                ("DSK/UNPRIORITIZED", "EARTH", "SUN", "IAU_EARTH", "LT", "MOON", sp),
                ("DSK/UNPRIORITIZED/SURFACES = 9", "MOON", "SUN", "IAU_MOON", "LT", "EARTH", sp),
                ("", "MOON", "SUN", "IAU_MOON", "LT", "EARTH", sp))
            @test same_kind(() -> g_illumg(m, t, src, et, fix, cor, o, p), () -> illumg(ctx, m, t, src, et, fix, cor, o, p))
            @test same_kind(() -> g_illumf(m, t, src, et, fix, cor, o, p), () -> illumf(ctx, m, t, src, et, fix, cor, o, p))
            src == "SUN" && @test same_kind(() -> g_ilumin(m, t, et, fix, cor, o, p), () -> ilumin(ctx, m, t, et, fix, cor, o, p))
        end
        @test neo_kind(() -> ilumin(ctx, "ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "EARTH", [1.0, 2.0])) == :badarg
        # Inputs that CSPICE accepts: a source that is the target, and a method with a
        # limb keyword (the keyword is not used).
        for (m, src) in (("ELLIPSOID", "MOON"), ("TANGENT/ELLIPSOID", "SUN"))
            a = illumf(ctx, m, "MOON", src, et, "IAU_MOON", "LT", "EARTH", sp)
            b = g_illumf(m, "MOON", src, et, "IAU_MOON", "LT", "EARTH", sp)
            @test maximum(abs.(collect(a[3:5]) .- collect(b[3:5]))) < ILL_TOL.angle && a[6:7] == b[6:7]
        end
        for (t, src, o, cor) in (("MOON", "SUN", "EARTH", "FOO"), ("MOON", "SUN", "EARTH", "XLT"),
                ("MOON", "SUN", "EARTH", "XCN+S"), ("NOBODY", "SUN", "EARTH", "LT"),
                ("MOON", "NOBODY", "EARTH", "LT"), ("MOON", "SUN", "NOBODY", "LT"),
                ("MOON", "MOON", "EARTH", "LT"), ("MOON", "SUN", "MOON", "LT"), ("MOON", "SUN", "SUN", "LT"))
            @test same_kind(() -> g_phaseq(et, t, src, o, cor), () -> phaseq(ctx, et, t, src, o, cor))
        end
        # trgsep reads the shape from its first 6 characters, so that "SPHEREX" is a sphere
        # and "POINTX" is not recognized; an observer inside a sphere has no angular radius.
        for (t1, s1, f1, t2, s2, f2, o, cor) in (
                ("MOON", "", "IAU_MOON", "SUN", "POINT", "IAU_SUN", "EARTH", "LT"),
                ("MOON", "POINT", "", "SUN", "POINT", "IAU_SUN", "EARTH", "LT"),
                ("MOON", "POINT", "IAU_MOON", "SUN", "POINT", "IAU_SUN", "EARTH", ""),
                ("NOBODY", "POINT", "IAU_MOON", "SUN", "POINT", "IAU_SUN", "EARTH", "LT"),
                ("MOON", "POINT", "IAU_MOON", "SUN", "POINT", "IAU_SUN", "NOBODY", "LT"),
                ("MOON", "POINT", "IAU_MOON", "MOON", "POINT", "IAU_SUN", "EARTH", "LT"),
                ("MOON", "POINT", "IAU_MOON", "SUN", "POINT", "IAU_SUN", "SUN", "LT"),
                ("MOON", "POINT", "IAU_MOON", "SUN", "POINT", "IAU_SUN", "EARTH", "FOO"),
                ("MOON", "POINTX", "IAU_MOON", "SUN", "POINT", "IAU_SUN", "EARTH", "LT"),
                ("EARTH", "SPHERE", "IAU_EARTH", "SUN", "POINT", "IAU_SUN", "EARTH BARYCENTER", "LT"),
                ("EARTH BARYCENTER", "SPHERE", "IAU_EARTH", "SUN", "POINT", "IAU_SUN", "MOON", "LT"))
            @test same_kind(() -> g_trgsep(et, t1, s1, f1, t2, s2, f2, o, cor),
                () -> trgsep(ctx, et, t1, s1, f1, t2, s2, f2, o, cor))
        end
        @test abs(trgsep(ctx, et, "MOON", "SPHEREX", "IAU_MOON", "SUN", " point junk", "IAU_SUN", "EARTH", "lt + s") -
                  g_trgsep(et, "MOON", "SPHEREX", "IAU_MOON", "SUN", " point junk", "IAU_SUN", "EARTH", "lt + s")) <
              ILL_TOL.angle
    end

    @testset "no allocation and inferred types on ellipsoids" begin
        et = 2.5e8
        sp = srfrec(ks, 301, 0.3, 0.2)
        for c in ("NONE", "LT", "LT+S", "CN", "CN+S")
            @test gwarm(ill_alloc_min, ks, et, c, sp) == 0
            @test gwarm(ill_alloc_g, ks, et, c, sp) == 0
            @test gwarm(ill_alloc_f, ks, et, c, sp) == 0
            @test gwarm(ill_alloc_ph, ks, et, c) == 0
            @test gwarm(ill_alloc_sep, ks, et, c) == 0
        end
        T5 = Tuple{Float64,SVector{3,Float64},Float64,Float64,Float64}
        @test (@inferred ilumin(ks, "ELLIPSOID", 301, et, 10020, "CN+S", 399, sp)) isa T5
        @test (@inferred illumg(ks, "ELLIPSOID", 301, 399, et, 10020, "LT+S", 10, sp)) isa T5
        @test (@inferred illumf(ks, "ELLIPSOID", 301, 10, et, 10020, "LT", 399, sp)) isa
              Tuple{Float64,SVector{3,Float64},Float64,Float64,Float64,Bool,Bool}
        @test (@inferred phaseq(ks, et, 301, 10, 399, "LT+S")) isa Float64
        @test (@inferred trgsep(ks, et, 301, "SPHERE", 10020, 10, "SPHERE", 10010, 399, "CN+S")) isa Float64
        @test (@inferred ilumin(ctx, "ELLIPSOID", "MOON", et, "IAU_MOON", "CN+S", "EARTH", sp)) isa T5
    end

    @testset "epoch as a dual number: derivative against the central difference" begin
        st = GeoStats()
        h = 20.0
        sp = srfrec(ks, 301, 0.3, 0.2)
        spe = srfrec(ks, 399, -1.0, 0.5)
        for et in record_safe_epochs(rng, nsample(6, 24)), cor in ("NONE", "LT", "LT+S", "CN", "CN+S")
            fs = Pair{String,Function}[
                "ilumin angles" => t -> SVector(ilumin(ks, "ELLIPSOID", 301, t, 10020, cor, 399, sp)[3:5]...),
                "ilumin srfvec" => t -> ilumin(ks, "ELLIPSOID", 399, t, 10013, cor, 301, spe)[2],
                "illumg angles" => t -> SVector(illumg(ks, "ELLIPSOID", 301, 399, t, 10020, cor, 10, sp)[3:5]...),
                "phaseq" => t -> SVector(phaseq(ks, t, 301, 10, 399, cor)),
                "trgsep" => t -> SVector(trgsep(ks, t, 301, "SPHERE", 10020, 10, "SPHERE", 10010, 399, cor))]
            for (name, f) in fs
                d = f(ForwardDiff.Dual(et, 1.0))
                @test isapprox(ForwardDiff.value.(d), f(et); rtol=1e-12)
                fd = geo_fd(f, et, h)
                gput!(st, name, norm(ForwardDiff.partials.(d, 1) .- fd) / norm(fd))
            end
        end
        greport("dual epoch against central difference", st)
        for k in keys(st.d)
            @test st.d[k] < 1e-4
        end
        @test neo_kind(() -> ilumin(ks, "DSK/UNPRIORITIZED", 301, ForwardDiff.Dual(1.0e8, 1.0), 10020, "LT", 399, sp)) == :notsupported
    end
end
