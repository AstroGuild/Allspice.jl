# Geometry group: sub-observer point, sub-solar point, surface intercept (`subpnt`,
# `subslr`, `sincpt`, `srfrec`) against CSPICE.
#
# Error measures: the norm of the difference from the CSPICE vector over the norm of
# the CSPICE vector, and the absolute difference of the epoch in seconds.
#
# Why some values are near 1e-10 and not 1e-16. A ray starts at the observer, so one
# unit in the last place of the observer position (3e-8 km at 1 au) moves the ray
# sideways by that amount, which is 1e-11 of the radius of Mercury or of the Moon. The
# intercept of `sincpt`, and the sub-point of a near-point method on DSK data (the ray
# goes from the observer to the near point of the ellipsoid), move with the ray; the
# vector from the observer to the point does not show it, because it is 1 au long. A
# ray that grazes the surface multiplies the shift by `1 / cos(incidence)`. With equal
# inputs the two libraries give equal outputs.

isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")
import ForwardDiff

const SUB_METHODS = ["NEAR POINT/ELLIPSOID", "INTERCEPT/ELLIPSOID", "Near point: ellipsoid",
    "Intercept: ellipsoid", "ELLIPSOID/NEAR POINT", "nadir / ellipsoid", " intercept/ ELLIPSOID "]
const SUB_DSK_METHODS = ["NADIR/DSK/UNPRIORITIZED", "INTERCEPT/DSK/UNPRIORITIZED",
    "NADIR/DSK/UNPRIORITIZED/SURFACES = 1", "INTERCEPT/DSK/UNPRIORITIZED/SURFACES = 2",
    "dsk/unprioritized/nadir/surfaces = 2, 1", "INTERCEPT / DSK / UNPRIORITIZED / SURFACES = 1"]
const INT_DSK_METHODS = ["DSK/UNPRIORITIZED", "DSK/UNPRIORITIZED/SURFACES = 1",
    "DSK/UNPRIORITIZED/SURFACES = 2", "dsk / unprioritized / surfaces = 2, 1"]

# Measured: 5e-16 (sub-points), 3e-15 (vectors), 0 (epochs) on ellipsoids; 7e-11 for
# intercepts of rays and 2e-11 for the near-point methods on DSK data (see the header).
const SUB_TOL = (point=5e-14, vec=5e-13, epoch=1e-7, graze=5e-9, dsk=1e-9)

sub_alloc_pnt(ks, m, et, c) = @allocated subpnt(ks, m, 301, et, 10020, c, 399)
sub_alloc_slr(ks, m, et, c) = @allocated subslr(ks, m, 301, et, 10020, c, 399)
sub_alloc_int(ks, et, c, dv) = @allocated sincpt(ks, "ELLIPSOID", 301, et, 10020, c, 399, 1, dv)
sub_alloc_rec(ks) = @allocated srfrec(ks, 301, 0.1, 0.2)

@gf_light group("geom_subpoints: subpnt, subslr, sincpt, srfrec", Kernels.FULL_SET_FILES...) do
    ctx = geo_setup(; dsk=true)
    ks = Allspice.kernelset(ctx)
    rng = MersenneTwister(SEED + 101)

    @testset "subpnt and subslr on ellipsoids, 600 random cases" begin
        st = GeoStats()
        nslr = 0
        for k in 1:600
            t = rand(rng, GEO_BODIES)
            o = rand(rng, setdiff(GEO_BODIES, [t]))
            fix = rand(rng, GEO_FIXED[t])
            cor = rand(rng, GEO_CORRS)
            m = rand(rng, SUB_METHODS)
            et = geo_epoch(rng)
            a = subpnt(ctx, m, t, et, fix, cor, o)
            b = g_subpnt(m, t, et, fix, cor, o)
            gput!(st, "subpnt spoint", grel(a[1], b[1]))
            gput!(st, "subpnt trgepc (s)", abs(a[2] - b[2]))
            gput!(st, "subpnt srfvec", grel(a[3], b[3]))
            if k % 10 == 0
                # The forms with IDs and a parsed correction give the same values.
                c = subpnt(ctx, m, GEO_IDS[t], et, namfrm(ctx, fix), Aberration(cor), GEO_IDS[o])
                @test c == a
            end
            (t == "SUN" || startswith(cor, "X")) && continue
            nslr += 1
            a = subslr(ctx, m, t, et, fix, cor, o)
            b = g_subslr(m, t, et, fix, cor, o)
            gput!(st, "subslr spoint", grel(a[1], b[1]))
            gput!(st, "subslr trgepc (s)", abs(a[2] - b[2]))
            gput!(st, "subslr srfvec", grel(a[3], b[3]))
        end
        greport("ellipsoid", st)
        @test nslr > 200
        @test gget(st, "subpnt spoint") < SUB_TOL.point
        @test gget(st, "subslr spoint") < SUB_TOL.point
        @test gget(st, "subpnt srfvec") < SUB_TOL.vec
        @test gget(st, "subslr srfvec") < SUB_TOL.vec
        @test gget(st, "subpnt trgepc (s)") < SUB_TOL.epoch
        @test gget(st, "subslr trgepc (s)") < SUB_TOL.epoch
    end

    @testset "sincpt on ellipsoids, 800 random rays" begin
        st = GeoStats()
        nhit = 0
        nmiss = 0
        for k in 1:800
            t = rand(rng, GEO_BODIES)
            o = rand(rng, setdiff(GEO_BODIES, [t]))
            fix = rand(rng, GEO_FIXED[t])
            cor = rand(rng, GEO_CORRS)
            et = geo_epoch(rng)
            dref = rand(rng, ["J2000", fix, "IAU_EARTH", "GEO_SITE", "ECLIPJ2000"])
            r, _ = SPICE.spkpos(t, et, "J2000", cor, o)
            rad = SPICE.bodvrd(t, "RADII", 3)[1]
            ang = asin(rad / norm(r)) * 1.3 * rand(rng)
            dv = SPICE.pxform("J2000", dref, et) * geo_offset_dir(rng, r, ang) * (0.1 + 5rand(rng))
            a = sincpt(ctx, "ELLIPSOID", t, et, fix, cor, o, dref, dv)
            b = g_sincpt("ELLIPSOID", t, et, fix, cor, o, dref, dv)
            if (a === nothing) != (b === nothing)
                # The found flag can differ only for a ray on the limb.
                gput!(st, "found flag differs: distance from the limb (rad)", abs(ang - asin(rad / norm(r))))
                continue
            end
            if a === nothing
                nmiss += 1
                continue
            end
            nhit += 1
            # The intercept of a grazing ray is ill-conditioned: weigh by the incidence.
            cosi = abs(dot(normalize(b[1]), normalize(b[3])))
            gput!(st, cosi > 0.05 ? "sincpt spoint" : "sincpt spoint, grazing rays", grel(a[1], b[1]))
            gput!(st, "sincpt trgepc (s)", abs(a[2] - b[2]))
            gput!(st, cosi > 0.05 ? "sincpt srfvec" : "sincpt srfvec, grazing rays", grel(a[3], b[3]))
        end
        greport("ellipsoid", st)
        @test nhit > 300 && nmiss > 100
        @test gget(st, "sincpt spoint") < SUB_TOL.graze
        @test gget(st, "sincpt srfvec") < SUB_TOL.vec
        @test gget(st, "sincpt trgepc (s)") < SUB_TOL.epoch
        @test !haskey(st.d, "found flag differs: distance from the limb (rad)") ||
              gget(st, "found flag differs: distance from the limb (rad)") < 1e-12
    end

    @testset "DSK shapes: Moon, 250 random cases" begin
        st = GeoStats()
        nhit = 0
        nmiss = 0
        for k in 1:250
            o = rand(rng, ["EARTH", "SUN", "VENUS", "MERCURY"])
            cor = rand(rng, GEO_CORRS)
            et = geo_epoch(rng)
            m = rand(rng, SUB_DSK_METHODS)
            near = occursin("nadir", lowercase(m))
            a = subpnt(ctx, m, "MOON", et, "IAU_MOON", cor, o)
            b = g_subpnt(m, "MOON", et, "IAU_MOON", cor, o)
            gput!(st, near ? "subpnt spoint, nadir" : "subpnt spoint, intercept", grel(a[1], b[1]))
            gput!(st, "subpnt trgepc (s)", abs(a[2] - b[2]))
            gput!(st, "subpnt srfvec", grel(a[3], b[3]))
            if !startswith(cor, "X")
                a = subslr(ctx, m, "MOON", et, "IAU_MOON", cor, o)
                b = g_subslr(m, "MOON", et, "IAU_MOON", cor, o)
                gput!(st, near ? "subslr spoint, nadir" : "subslr spoint, intercept", grel(a[1], b[1]))
                gput!(st, "subslr srfvec", grel(a[3], b[3]))
            end
            m = rand(rng, INT_DSK_METHODS)
            r, _ = SPICE.spkpos("MOON", et, "J2000", cor, o)
            ang = asin(1737.4 / norm(r)) * 1.25 * rand(rng)
            dv = geo_offset_dir(rng, r, ang)
            a = sincpt(ctx, m, "MOON", et, "IAU_MOON", cor, o, "J2000", dv)
            b = g_sincpt(m, "MOON", et, "IAU_MOON", cor, o, "J2000", dv)
            if (a === nothing) != (b === nothing)
                gput!(st, "found flag differs", 1.0)
                continue
            end
            a === nothing && (nmiss += 1; continue)
            nhit += 1
            gput!(st, "sincpt spoint", grel(a[1], b[1]))
            gput!(st, "sincpt trgepc (s)", abs(a[2] - b[2]))
            gput!(st, "sincpt srfvec", grel(a[3], b[3]))
        end
        greport("DSK", st)
        @test nhit > 100 && nmiss > 20
        @test !haskey(st.d, "found flag differs")
        @test gget(st, "subpnt spoint, intercept") < SUB_TOL.point
        @test gget(st, "subpnt spoint, nadir") < SUB_TOL.dsk
        @test gget(st, "subslr spoint, intercept") < SUB_TOL.point
        @test gget(st, "subslr spoint, nadir") < SUB_TOL.dsk
        @test gget(st, "sincpt spoint") < SUB_TOL.dsk
        @test gget(st, "subpnt trgepc (s)") < SUB_TOL.epoch
        @test gget(st, "sincpt trgepc (s)") < SUB_TOL.epoch
    end

    @testset "srfrec" begin
        worst = 0.0
        for k in 1:300
            body = rand(rng, [10, 199, 299, 399, 301, 499, 599])
            lon = 2pi * (rand(rng) - 0.5)
            lat = pi * (rand(rng) - 0.5)
            worst = max(worst, grel(srfrec(ctx, body, lon, lat), g_srfrec(body, lon, lat)))
        end
        measured(@sprintf("  %-64s max=%.3e", "srfrec", worst))
        @test worst < 1e-15
        @test srfrec(ctx, "EARTH", 0.3, 0.2) == srfrec(ctx, 399, 0.3, 0.2)
        @test neo_kind(() -> srfrec(ctx, "NOBODY", 0.0, 0.0)) == :idcodenotfound
        @test same_kind(() -> g_srfrec(3, 0.0, 0.0), () -> srfrec(ctx, 3, 0.0, 0.0); expect=:kernelvarnotfound)
    end

    @testset "error conditions, kinds of CSPICE" begin
        et = 1.0e8
        e3 = [1.0, 0.0, 0.0]
        for (m, t, fix, cor, o) in (
                ("NEAR POINT/ELLIPSOID", "MOON", "IAU_MOON", "FOO", "EARTH"),
                ("NEAR POINT/ELLIPSOID", "MOON", "IAU_MOON", "S", "EARTH"),
                ("NEAR POINT/ELLIPSOID", "NOBODY", "IAU_MOON", "LT", "EARTH"),
                ("NEAR POINT/ELLIPSOID", "MOON", "IAU_MOON", "LT", "NOBODY"),
                ("NEAR POINT/ELLIPSOID", "MOON", "IAU_MOON", "LT", "MOON"),
                ("NEAR POINT/ELLIPSOID", "MOON", "NOFRAME", "LT", "EARTH"),
                ("NEAR POINT/ELLIPSOID", "MOON", "IAU_EARTH", "LT", "EARTH"),
                ("NEAR POINT/ELLIPSOID", "MOON", "J2000", "LT", "EARTH"),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "EARTH"),
                ("TANGENT/ELLIPSOID", "MOON", "IAU_MOON", "LT", "EARTH"),
                ("NEAR POINT/SPHERE", "MOON", "IAU_MOON", "LT", "EARTH"),
                ("NEAR POINT/ELLIPSOID/ELLIPSOID", "MOON", "IAU_MOON", "LT", "EARTH"),
                ("NEAR POINT/DSK/UNPRIORITIZED", "MOON", "IAU_MOON", "LT", "EARTH"),
                ("NADIR/DSK", "MOON", "IAU_MOON", "LT", "EARTH"),
                ("NADIR/DSK/UNPRIORITIZED/SURFACES = NOSURFACE", "MOON", "IAU_MOON", "LT", "EARTH"),
                ("NADIR/DSK/UNPRIORITIZED/SURFACES = 9", "MOON", "IAU_MOON", "LT", "EARTH"),
                ("NADIR/DSK/UNPRIORITIZED", "EARTH", "IAU_EARTH", "LT", "MOON"),
                ("", "MOON", "IAU_MOON", "LT", "EARTH"))
            @test same_kind(() -> g_subpnt(m, t, et, fix, cor, o), () -> subpnt(ctx, m, t, et, fix, cor, o))
            # `subslr` accepts an observer that is the target.
            t == o || @test same_kind(() -> g_subslr(m, t, et, fix, cor, o), () -> subslr(ctx, m, t, et, fix, cor, o))
        end
        @test grel(subslr(ctx, "NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "MOON")[1],
            g_subslr("NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "MOON")[1]) < SUB_TOL.point
        # No ephemeris data for Mars: CSPICE gives SPKINSUFFDATA, Allspice the kind of its SPK layer.
        @test c_kind(() -> g_subpnt("NEAR POINT/ELLIPSOID", "MARS", et, "IAU_MARS", "LT", "EARTH")) == :spkinsuffdata
        @test neo_kind(() -> subpnt(ctx, "NEAR POINT/ELLIPSOID", "MARS", et, "IAU_MARS", "LT", "EARTH")) == :spkinsuffdata
        @test neo_kind(() -> sincpt(ctx, "ELLIPSOID", "MARS", et, "IAU_MARS", "LT", "EARTH", "J2000", e3)) == :spkinsuffdata
        @test same_kind(() -> g_subslr("NEAR POINT/ELLIPSOID", "SUN", et, "IAU_SUN", "LT", "EARTH"),
            () -> subslr(ctx, "NEAR POINT/ELLIPSOID", "SUN", et, "IAU_SUN", "LT", "EARTH"); expect=:invalidtarget)
        @test same_kind(() -> g_subslr("NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "XLT", "EARTH"),
            () -> subslr(ctx, "NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "XLT", "EARTH"); expect=:notsupported)
        for (m, t, fix, cor, o, dref, dv) in (
                ("ELLIPSOID", "MOON", "IAU_MOON", "FOO", "EARTH", "J2000", e3),
                ("ELLIPSOID", "NOBODY", "IAU_MOON", "LT", "EARTH", "J2000", e3),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "NOBODY", "J2000", e3),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "MOON", "J2000", e3),
                ("ELLIPSOID", "MOON", "NOFRAME", "LT", "EARTH", "J2000", e3),
                ("ELLIPSOID", "MOON", "IAU_EARTH", "LT", "EARTH", "J2000", e3),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "EARTH", "NOFRAME", e3),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "EARTH", "J2000", zeros(3)),
                ("NEAR POINT/ELLIPSOID", "MOON", "IAU_MOON", "LT", "EARTH", "J2000", e3),
                ("INTERCEPT/DSK/UNPRIORITIZED", "MOON", "IAU_MOON", "LT", "EARTH", "J2000", e3),
                ("SPHERE", "MOON", "IAU_MOON", "LT", "EARTH", "J2000", e3),
                ("DSK/UNPRIORITIZED", "EARTH", "IAU_EARTH", "LT", "MOON", "J2000", e3),
                ("DSK/UNPRIORITIZED/SURFACES = 9", "MOON", "IAU_MOON", "LT", "EARTH", "J2000", e3),
                ("", "MOON", "IAU_MOON", "LT", "EARTH", "J2000", e3))
            @test same_kind(() -> g_sincpt(m, t, et, fix, cor, o, dref, dv),
                () -> sincpt(ctx, m, t, et, fix, cor, o, dref, dv))
        end
        @test neo_kind(() -> sincpt(ctx, "ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "EARTH", "J2000", [1.0, 0.0])) == :badarg
        # Radii that are not valid.
        radii = gdpool(ctx, "BODY301_RADII")
        for bad in ([1737.4, 1737.4], [1737.4, 1737.4, 1737.4, 1.0], [1737.4, -1.0, 1737.4], [0.0, 1737.4, 1737.4])
            both_pdpool(ctx, "BODY301_RADII", bad)
            @test same_kind(() -> g_subpnt("NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "EARTH"),
                () -> subpnt(ctx, "NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "EARTH"))
            @test same_kind(() -> g_sincpt("ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "EARTH", "J2000", e3),
                () -> sincpt(ctx, "ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "EARTH", "J2000", e3))
        end
        both_dvpool(ctx, "BODY301_RADII")
        @test same_kind(() -> g_subslr("NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "EARTH"),
            () -> subslr(ctx, "NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "EARTH"); expect=:kernelvarnotfound)
        both_pdpool(ctx, "BODY301_RADII", radii)
        ks = Allspice.kernelset(ctx)
        @test grel(subpnt(ctx, "NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "EARTH")[1],
            g_subpnt("NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "EARTH")[1]) < SUB_TOL.point
    end

    @testset "no allocation and inferred types on ellipsoids" begin
        ks = Allspice.kernelset(ctx)
        et = 2.5e8
        dv = SVector{3,Float64}(spkpos(ctx, 301, et, 1, "LT+S", 399)[1])
        for c in ("NONE", "LT", "LT+S", "CN", "CN+S", "XLT+S", "XCN")
            @test gwarm(sub_alloc_pnt, ks, "NEAR POINT/ELLIPSOID", et, c) == 0
            @test gwarm(sub_alloc_pnt, ks, "INTERCEPT/ELLIPSOID", et, c) == 0
            @test gwarm(sub_alloc_int, ks, et, c, dv) == 0
            startswith(c, "X") && continue
            @test gwarm(sub_alloc_slr, ks, "NEAR POINT/ELLIPSOID", et, c) == 0
            @test gwarm(sub_alloc_slr, ks, "INTERCEPT/ELLIPSOID", et, c) == 0
        end
        @test gwarm(sub_alloc_rec, ks) == 0
        T3 = Tuple{SVector{3,Float64},Float64,SVector{3,Float64}}
        @test (@inferred subpnt(ks, "NEAR POINT/ELLIPSOID", 301, et, 10020, "CN+S", 399)) isa T3
        @test (@inferred subslr(ks, "INTERCEPT/ELLIPSOID", 301, et, 10020, "LT+S", 399)) isa T3
        @test (@inferred subpnt(ctx, "NEAR POINT/ELLIPSOID", "MOON", et, "IAU_MOON", "CN+S", "EARTH")) isa T3
        @test (@inferred Nothing sincpt(ks, "ELLIPSOID", 301, et, 10020, "CN+S", 399, 1, dv)) isa T3
        @test (@inferred Nothing sincpt(ks, "ELLIPSOID", 301, et, 10020, "CN+S", 399, 1, -dv)) === nothing
        @test (@inferred srfrec(ks, 301, 0.1, 0.2)) isa SVector{3,Float64}
    end

    @testset "epoch as a dual number: derivative against the central difference" begin
        ks = Allspice.kernelset(ctx)
        st = GeoStats()
        h = 20.0
        for et in record_safe_epochs(rng, nsample(6, 24))
            for cor in ("NONE", "LT", "LT+S", "CN+S", "XCN+S")
                dv = SVector{3,Float64}(spkpos(ctx, 301, et, 1, cor, 399)[1]) + SVector(300.0, -200.0, 100.0)
                fs = Pair{String,Function}[
                    "subpnt near" => t -> vcat(subpnt(ks, "NEAR POINT/ELLIPSOID", 301, t, 10020, cor, 399)[[1, 3]]...),
                    "subpnt intercept" => t -> vcat(subpnt(ks, "INTERCEPT/ELLIPSOID", 399, t, 10013, cor, 301)[[1, 3]]...),
                    "sincpt" => t -> vcat(sincpt(ks, "ELLIPSOID", 301, t, 10020, cor, 399, 1, dv)[[1, 3]]...)]
                startswith(cor, "X") || push!(fs,
                    "subslr near" => t -> vcat(subslr(ks, "NEAR POINT/ELLIPSOID", 301, t, 10020, cor, 399)[[1, 3]]...),
                    "subslr intercept" => t -> vcat(subslr(ks, "INTERCEPT/ELLIPSOID", 399, t, 10013, cor, 10)[[1, 3]]...))
                for (name, f) in fs
                    d = f(ForwardDiff.Dual(et, 1.0))
                    @test isapprox(ForwardDiff.value.(d), f(et); rtol=1e-12)
                    fd = geo_fd(f, et, h)
                    gput!(st, name, norm(ForwardDiff.partials.(d, 1) .- fd) / norm(fd))
                end
                # The epoch at the point follows the light time.
                te = subpnt(ks, "NEAR POINT/ELLIPSOID", 301, ForwardDiff.Dual(et, 1.0), 10020, cor, 399)[2]
                @test abs(ForwardDiff.partials(te, 1) - 1.0) < 1e-5
            end
        end
        greport("dual epoch against central difference", st)
        for k in keys(st.d)
            @test st.d[k] < 1e-4
        end
        # DSK shapes take a plain epoch.
        @test neo_kind(() -> subpnt(ks, "NADIR/DSK/UNPRIORITIZED", 301, ForwardDiff.Dual(1.0e8, 1.0), 10020, "LT", 399)) == :notsupported
    end
end
