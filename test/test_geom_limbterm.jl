# Geometry group: limb points, terminator points and the tangent point of a ray
# (`limbpt`, `termpt`, `tangpt`) against CSPICE.
#
# Error measures: the number of points per cutting half plane must be equal; then the
# norm of the difference from the CSPICE vector over the norm of the CSPICE vector, and
# the absolute difference of the epochs in seconds.
#
# Tangent points on DSK data come from a bisection on the ray angle that stops at the
# angular tolerance `soltol`; the two libraries bisect the same interval, so they agree
# far below `soltol`, but the last step can fall on the other side for inputs that
# differ by one unit in the last place. The tolerance on these points is
# `soltol` times the distance to the observer over the radius.
#
# `tangpt`: the tangent point comes from `npedln`, which finds the near point of an
# ellipse in a plane normal to the ray. One unit in the last place of the ray direction
# changes the result of `npedln_c` itself by 1e-10 of the radius (measured in the test),
# and the two libraries give the same output for the same input.
#
# The default run has 20 random cases on the DSK surface of two parts; the exhaustive
# run (`ALLSPICE_TEST_FULL=1`) has 40. The other samples are the same in the two runs.

isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")

# Measured: 2e-15 (points on ellipsoids), 3e-16 (vectors), 0 (epochs); 2e-9 for tangent
# points on DSK data; tangpt 4e-11 (tangent point and surface point), 1e-15 (range).
const LT_TOL = (point=5e-13, vec=5e-13, epoch=1e-7, dsk=2e-7, tang=2e-8)

"""
    lt_same_tangpt(c, a, rtol) -> Bool

True if every component of the two results of `tangpt` agrees within `rtol` of the size
of its output (the norm of a vector, the absolute value of a number; the epoch within
`rtol` seconds for each second of its absolute value).
"""
function lt_same_tangpt(c, a, rtol::Float64)
    @assert length(c) == length(a) == 6 && rtol > 0
    return all(1:6) do i
        scale = max(norm(a[i]), floatmin())
        all(abs.(collect(c[i]) .- collect(a[i])) .<= rtol * scale)
    end
end

lt_alloc_tang(ks, et, c, loc, dv) = @allocated tangpt(ks, "ELLIPSOID", 301, et, 10020, c, loc, 399, 1, dv)

"Compare two results of `limbpt` or `termpt`; false if the point counts differ."
function lt_compare!(st, name, a, b)
    a[1] == b[1] || return false
    for i in eachindex(a[2])
        gput!(st, "$name points", grel(a[2][i], b[2][:, i]))
        gput!(st, "$name epochs (s)", abs(a[3][i] - b[3][i]))
        gput!(st, "$name vectors", grel(a[4][i], b[4][:, i]))
    end
    return true
end

@gf_light group("geom_limbterm: limbpt, termpt, tangpt", Kernels.FULL_SET_FILES...) do
    paths = [geo_kernels(; dsk=true)..., geo_parts_dsk()]
    oracle_load(paths)
    ctx = neo_context(paths)
    ks = Allspice.kernelset(ctx)
    rng = MersenneTwister(SEED + 106)

    @testset "limbpt and termpt on ellipsoids, 400 random cases" begin
        st = GeoStats()
        ncount = 0
        nterm = 0
        for k in 1:400
            t = rand(rng, GEO_BODIES)
            o = rand(rng, setdiff(GEO_BODIES, [t]))
            fix = rand(rng, GEO_FIXED[t])
            cor = rand(rng, GEO_CORRS)
            et = geo_epoch(rng)
            ref = normalize(randn(rng, 3)) .* (0.5 + 3rand(rng))
            rol = (0.01 + rand(rng)) * rand(rng, (-1, 1))
            nc = rand(rng, 1:8)
            loc = rand(rng, ["CENTER", "ELLIPSOID LIMB", "center", " Ellipsoid  Limb"])
            m = occursin("imb", lowercase(loc)) ? rand(rng, ["TANGENT/ELLIPSOID", "ellipsoid/tangent"]) :
                rand(rng, ["TANGENT/ELLIPSOID", "GUIDED/ELLIPSOID", "Ellipsoid / Guided"])
            a = limbpt(ctx, m, t, et, fix, cor, loc, o, ref, rol, nc, 1e-4, 1e-7, 20)
            b = g_limbpt(m, t, et, fix, cor, loc, o, ref, rol, nc, 1e-4, 1e-7, 20)
            lt_compare!(st, "limbpt", a, b) || (ncount += 1)
            @test all(==(1), a[1]) && length(a[2]) == nc
            if k % 10 == 0
                c = limbpt(ctx, m, GEO_IDS[t], et, namfrm(ctx, fix), Aberration(cor), loc, GEO_IDS[o], SVector{3}(ref), rol, nc, 1e-4, 1e-7, 20)
                @test c == a
            end
            (t == "SUN" || startswith(cor, "X")) && continue
            nterm += 1
            src = rand(rng, ["SUN", "SUN", rand(rng, setdiff(GEO_BODIES, [t]))])
            loc = rand(rng, ["CENTER", "ELLIPSOID TERMINATOR", "ellipsoid terminator"])
            shd = rand(rng, ["UMBRAL", "PENUMBRAL"])
            m = shd * "/" * (loc == "CENTER" ? rand(rng, ["TANGENT", "GUIDED"]) : "TANGENT") * "/ELLIPSOID"
            a = termpt(ctx, m, src, t, et, fix, cor, loc, o, ref, rol, nc, 1e-4, 1e-7, 100)
            b = g_termpt(m, src, t, et, fix, cor, loc, o, ref, rol, nc, 1e-4, 1e-7, 100)
            lt_compare!(st, "termpt", a, b) || (ncount += 1)
        end
        greport("ellipsoid", st)
        @test ncount == 0
        @test nterm > 120
        for name in ("limbpt", "termpt")
            @test gget(st, "$name points") < LT_TOL.point
            @test gget(st, "$name vectors") < LT_TOL.vec
            @test gget(st, "$name epochs (s)") < LT_TOL.epoch
        end
    end

    @testset "limbpt and termpt on DSK data: Moon, 60 random cases" begin
        st = GeoStats()
        ncount = 0
        nmulti = 0
        for k in 1:60
            o = rand(rng, ["EARTH", "SUN", "VENUS"])
            cor = rand(rng, GEO_CORRS)
            et = geo_epoch(rng)
            surf = rand(rng, ["/SURFACES = 1", "/SURFACES = 2", "/SURFACES = 2", ""])
            ref = normalize(randn(rng, 3))
            rol = 0.2 + rand(rng)
            nc = rand(rng, 1:4)
            loc = rand(rng, ["CENTER", "ELLIPSOID LIMB"])
            typ = loc == "CENTER" ? rand(rng, ["TANGENT", "GUIDED"]) : "TANGENT"
            m = typ * "/DSK/UNPRIORITIZED" * surf
            soltol = 1e-9
            a = limbpt(ctx, m, "MOON", et, "IAU_MOON", cor, loc, o, ref, rol, nc, 1e-4, soltol, 400)
            b = g_limbpt(m, "MOON", et, "IAU_MOON", cor, loc, o, ref, rol, nc, 1e-4, soltol, 400)
            lt_compare!(st, typ == "TANGENT" ? "limbpt tangent" : "limbpt guided", a, b) || (ncount += 1)
            nmulti += any(>(1), b[1])
            startswith(cor, "X") && continue
            loc = rand(rng, ["CENTER", "ELLIPSOID TERMINATOR"])
            typ = loc == "CENTER" ? rand(rng, ["TANGENT", "GUIDED"]) : "TANGENT"
            m = rand(rng, ["UMBRAL/", "PENUMBRAL/"]) * typ * "/DSK/UNPRIORITIZED" * surf
            a = termpt(ctx, m, "SUN", "MOON", et, "IAU_MOON", cor, loc, o == "SUN" ? "EARTH" : o, ref, rol, nc, 1e-4, soltol, 400)
            b = g_termpt(m, "SUN", "MOON", et, "IAU_MOON", cor, loc, o == "SUN" ? "EARTH" : o, ref, rol, nc, 1e-4, soltol, 400)
            lt_compare!(st, typ == "TANGENT" ? "termpt tangent" : "termpt guided", a, b) || (ncount += 1)
            nmulti += any(>(1), b[1])
        end
        greport("DSK", st)
        measured("  DSK: cases with more than one point in a half plane: $nmulti; cases with unequal point counts: $ncount")
        @test ncount == 0
        for name in ("limbpt tangent", "termpt tangent")
            @test gget(st, "$name points") < LT_TOL.dsk
            @test gget(st, "$name epochs (s)") < LT_TOL.epoch
        end
        for name in ("limbpt guided", "termpt guided")
            @test gget(st, "$name points") < LT_TOL.point
        end
    end

    @testset "DSK surface of two parts: several points in a half plane" begin
        st = GeoStats()
        ncount = 0
        nmulti = 0
        most = 0
        m3 = "/DSK/UNPRIORITIZED/SURFACES = 3"
        nparts = nsample(20, 40)
        for k in 1:nparts
            o = rand(rng, ["EARTH", "SUN", "VENUS"])
            cor = rand(rng, GEO_CORRS)
            et = geo_epoch(rng)
            # A reference vector near the axis that joins the two parts.
            ref = [0.2 * randn(rng), 0.2 * randn(rng), 1.0]
            loc = rand(rng, ["CENTER", "ELLIPSOID LIMB"])
            a = limbpt(ctx, "TANGENT" * m3, "MOON", et, "IAU_MOON", cor, loc, o, ref, 0.05, 4, 1e-4, 1e-9, 400)
            b = g_limbpt("TANGENT" * m3, "MOON", et, "IAU_MOON", cor, loc, o, ref, 0.05, 4, 1e-4, 1e-9, 400)
            lt_compare!(st, "limbpt", a, b) || (ncount += 1)
            nmulti += count(>(1), b[1])
            most = max(most, maximum(b[1]))
            startswith(cor, "X") && continue
            o == "SUN" && (o = "EARTH")
            loc = rand(rng, ["CENTER", "ELLIPSOID TERMINATOR"])
            shd = rand(rng, ["UMBRAL/", "PENUMBRAL/"])
            a = termpt(ctx, shd * "TANGENT" * m3, "SUN", "MOON", et, "IAU_MOON", cor, loc, o, ref, 0.05, 4, 1e-4, 1e-9, 400)
            b = g_termpt(shd * "TANGENT" * m3, "SUN", "MOON", et, "IAU_MOON", cor, loc, o, ref, 0.05, 4, 1e-4, 1e-9, 400)
            lt_compare!(st, "termpt", a, b) || (ncount += 1)
            nmulti += count(>(1), b[1])
            most = max(most, maximum(b[1]))
        end
        greport("DSK, two parts", st)
        measured("  DSK, two parts: half planes with more than one point: $nmulti (at most $most points); cases with unequal point counts: $ncount")
        @test ncount == 0
        @test nmulti > nparts ÷ 4 && most == 3
        @test gget(st, "limbpt points") < LT_TOL.dsk
        @test gget(st, "termpt points") < LT_TOL.dsk
        # More points than `maxn`: Allspice throws, where `limbpt_c` writes past the arrays.
        et = 1.0e8
        ref = SVector(0.0, 0.0, 1.0)
        n = limbpt(ctx, "TANGENT" * m3, "MOON", et, "IAU_MOON", "NONE", "CENTER", "EARTH", ref, 0.01, 4, 1e-4, 1e-9, 400)[1]
        if sum(n) > 4
            @test neo_kind(() -> limbpt(ctx, "TANGENT" * m3, "MOON", et, "IAU_MOON", "NONE", "CENTER", "EARTH", ref, 0.01, 4, 1e-4, 1e-9, 4)) == :outofroom
        else
            @test neo_kind(() -> limbpt(ctx, "TANGENT" * m3, "MOON", et, "IAU_MOON", "NONE", "CENTER", "EARTH", ref, 0.01, 4, 1e-4, 1e-9, 4)) === nothing
        end
    end

    @testset "tangpt, 800 random rays" begin
        st = GeoStats()
        nbehind = 0
        for k in 1:800
            t = rand(rng, GEO_BODIES)
            o = rand(rng, setdiff(GEO_BODIES, [t]))
            fix = rand(rng, GEO_FIXED[t])
            cor = rand(rng, GEO_CORRS)
            et = geo_epoch(rng)
            loc = rand(rng, ["TANGENT POINT", "SURFACE POINT", "tangent point", "Surface  Point"])
            dref = rand(rng, ["J2000", fix, "IAU_EARTH", "GEO_SITE", "ECLIPJ2000"])
            r, _ = SPICE.spkpos(t, et, "J2000", cor, o)
            rad = SPICE.bodvrd(t, "RADII", 3)[1]
            ang = asin(rad / norm(r)) * (0.2 + 3rand(rng))
            dj = geo_offset_dir(rng, r, ang)
            if k % 20 == 0          # a ray that points away from the target
                dj = -dj
                nbehind += 1
            end
            dv = SPICE.pxform("J2000", dref, et) * dj * (0.1 + 5rand(rng))
            a = tangpt(ctx, "ELLIPSOID", t, et, fix, cor, loc, o, dref, dv)
            b = g_tangpt("ELLIPSOID", t, et, fix, cor, loc, o, dref, dv)
            gput!(st, "tanpt", norm(collect(a[1]) .- b[1]) / max(norm(b[1]), rad))
            gput!(st, "alt over the radius", abs(a[2] - b[2]) / rad)
            gput!(st, "range", abs(a[3] - b[3]) / max(b[3], rad))
            gput!(st, "srfpt", grel(a[4], b[4]))
            gput!(st, "trgepc (s)", abs(a[5] - b[5]))
            gput!(st, "srfvec", grel(a[6], b[6]))
            if k % 10 == 0
                c = tangpt(ctx, "ellipsoid", GEO_IDS[t], et, namfrm(ctx, fix), Aberration(cor), loc, GEO_IDS[o], namfrm(ctx, dref), SVector{3}(dv))
                # The form with a `Vector` and the form with an `SVector` take the ray
                # through code that the compiler can arrange in two ways (it does so
                # with code coverage), so the results can differ by rounding: each
                # output agrees within 1e-12 of its size.
                @test lt_same_tangpt(c, a, 1e-12)
            end
            # The conditioning of `npedln_c`: one unit in the last place of the direction.
            if cor == "NONE" && dref == "J2000"
                p = -SPICE.spkpos(t, et, fix, "NONE", o)[1]
                d = SPICE.pxform("J2000", fix, et) * dv
                rr = SPICE.bodvrd(t, "RADII", 3)
                y = SPICE.npedln(rr[1], rr[2], rr[3], p, d)[1]
                z = SPICE.npedln(rr[1], rr[2], rr[3], p, [nextfloat(d[1]), d[2], prevfloat(d[3])])[1]
                gput!(st, "(npedln_c: change of the output for one unit in the last place)", norm(y - z) / rad)
            end
        end
        greport("tangpt", st)
        @test nbehind == 40
        for k in ("tanpt", "alt over the radius", "srfpt")
            @test gget(st, k) < LT_TOL.tang
        end
        @test gget(st, "range") < LT_TOL.vec
        @test gget(st, "srfvec") < LT_TOL.vec
        @test gget(st, "trgepc (s)") < LT_TOL.epoch
    end

    @testset "error conditions, kinds of CSPICE" begin
        et = 1.0e8
        z = [0.0, 0.0, 1.0]
        T = "TANGENT/ELLIPSOID"
        D = "TANGENT/DSK/UNPRIORITIZED"
        # (method, target, fixref, abcorr, corloc, obsrvr, refvec, rolstp, ncuts, schstp, soltol, maxn)
        for (m, t, fix, cor, loc, o, ref, rol, nc, sch, sol, mx) in (
                (T, "MOON", "IAU_MOON", "FOO", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (T, "NOBODY", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (T, "MOON", "IAU_MOON", "LT", "CENTER", "NOBODY", z, 0.1, 3, 1e-4, 1e-7, 10),
                (T, "MOON", "IAU_MOON", "LT", "CENTER", "MOON", z, 0.1, 3, 1e-4, 1e-7, 10),
                (T, "MOON", "NOFRAME", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (T, "MOON", "IAU_EARTH", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                ("NEAR POINT/ELLIPSOID", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                ("UMBRAL/TANGENT/ELLIPSOID", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                ("TANGENT/SPHERE", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (T, "MOON", "IAU_MOON", "LT", "LIMB", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (T, "MOON", "IAU_MOON", "LT", "ELLIPSOID TERMINATOR", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                ("GUIDED/ELLIPSOID", "MOON", "IAU_MOON", "LT", "ELLIPSOID LIMB", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (T, "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", zeros(3), 0.1, 3, 1e-4, 1e-7, 10),
                (T, "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 0),
                (T, "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 0, 1e-4, 1e-7, 10),
                (T, "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 11, 1e-4, 1e-7, 10),
                (T, "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.0, 3, 1e-4, 1e-7, 10),
                (D, "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 0.0, 1e-7, 10),
                (D, "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 0.0, 10),
                (D, "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-30, 10),
                (D, "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-30, 1e-7, 10),
                (D * "/SURFACES = 9", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (D, "EARTH", "IAU_EARTH", "LT", "CENTER", "MOON", z, 0.1, 3, 1e-4, 1e-7, 10),
                ("", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10))
            @test same_kind(() -> g_limbpt(m, t, 1.0e8, fix, cor, loc, o, ref, rol, nc, sch, sol, mx),
                () -> limbpt(ctx, m, t, 1.0e8, fix, cor, loc, o, ref, rol, nc, sch, sol, mx))
        end
        # A reference vector that is parallel to the observer-target line.
        # (each library gets its own position vector, which differs in the last place)
        r = SPICE.spkpos("MOON", et, "IAU_MOON", "NONE", "EARTH")[1]
        rn = spkpos(ctx, "MOON", et, "IAU_MOON", "NONE", "EARTH")[1]
        @test same_kind(() -> g_limbpt(T, "MOON", et, "IAU_MOON", "NONE", "CENTER", "EARTH", r, 0.1, 3, 1e-4, 1e-7, 10),
            () -> limbpt(ctx, T, "MOON", et, "IAU_MOON", "NONE", "CENTER", "EARTH", rn, 0.1, 3, 1e-4, 1e-7, 10); expect=:degeneratecase)
        U = "UMBRAL/TANGENT/ELLIPSOID"
        UD = "UMBRAL/TANGENT/DSK/UNPRIORITIZED"
        # (method, ilusrc, target, fixref, abcorr, corloc, obsrvr, refvec, rolstp, ncuts, schstp, soltol, maxn)
        for (m, src, t, fix, cor, loc, o, ref, rol, nc, sch, sol, mx) in (
                (U, "SUN", "MOON", "IAU_MOON", "FOO", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "IAU_MOON", "XLT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "NOBODY", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "SUN", "NOBODY", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "NOBODY", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "MOON", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "MOON", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "NOFRAME", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "IAU_EARTH", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (T, "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                ("UMBRAL/ELLIPSOID", "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                ("UMBRAL/TANGENT/NADIR/ELLIPSOID", "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "IAU_MOON", "LT", "TERMINATOR", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "IAU_MOON", "LT", "ELLIPSOID LIMB", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                ("UMBRAL/GUIDED/ELLIPSOID", "SUN", "MOON", "IAU_MOON", "LT", "ELLIPSOID TERMINATOR", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", zeros(3), 0.1, 3, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 0),
                (U, "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, -1, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 11, 1e-4, 1e-7, 10),
                (U, "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.0, 3, 1e-4, 1e-7, 10),
                (UD, "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 0.0, 1e-7, 10),
                (UD, "SUN", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 0.0, 10),
                (UD, "SUN", "EARTH", "IAU_EARTH", "LT", "CENTER", "MOON", z, 0.1, 3, 1e-4, 1e-7, 10))
            @test same_kind(() -> g_termpt(m, src, t, et, fix, cor, loc, o, ref, rol, nc, sch, sol, mx),
                () -> termpt(ctx, m, src, t, et, fix, cor, loc, o, ref, rol, nc, sch, sol, mx))
        end
        e3 = [1.0, 0.0, 0.0]
        for (m, t, fix, cor, loc, o, dref, dv) in (
                ("ELLIPSOID", "MOON", "IAU_MOON", "FOO", "TANGENT POINT", "EARTH", "J2000", e3),
                ("DSK/UNPRIORITIZED", "MOON", "IAU_MOON", "LT", "TANGENT POINT", "EARTH", "J2000", e3),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "CENTER", "EARTH", "J2000", e3),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "TANGENT POINT", "EARTH", "J2000", zeros(3)),
                ("ELLIPSOID", "NOBODY", "IAU_MOON", "LT", "TANGENT POINT", "EARTH", "J2000", e3),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "TANGENT POINT", "NOBODY", "J2000", e3),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "TANGENT POINT", "MOON", "J2000", e3),
                ("ELLIPSOID", "MOON", "NOFRAME", "LT", "TANGENT POINT", "EARTH", "J2000", e3),
                ("ELLIPSOID", "MOON", "IAU_EARTH", "LT", "TANGENT POINT", "EARTH", "J2000", e3),
                ("ELLIPSOID", "MOON", "IAU_MOON", "LT", "TANGENT POINT", "EARTH", "NOFRAME", e3))
            @test same_kind(() -> g_tangpt(m, t, et, fix, cor, loc, o, dref, dv),
                () -> tangpt(ctx, m, t, et, fix, cor, loc, o, dref, dv))
        end
        # An observer inside of the ellipsoid, and radii with a wrong count.
        radii = gdpool(ctx, "BODY301_RADII")
        both_pdpool(ctx, "BODY301_RADII", [5.0e5, 5.0e5, 5.0e5])
        @test same_kind(() -> g_tangpt("ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "TANGENT POINT", "EARTH", "J2000", e3),
            () -> tangpt(ctx, "ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "TANGENT POINT", "EARTH", "J2000", e3); expect=:invalidgeometry)
        both_pdpool(ctx, "BODY301_RADII", [1737.4, 1737.4])
        @test same_kind(() -> g_tangpt("ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "TANGENT POINT", "EARTH", "J2000", e3),
            () -> tangpt(ctx, "ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "TANGENT POINT", "EARTH", "J2000", e3))
        @test same_kind(() -> g_limbpt(T, "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
            () -> limbpt(ctx, T, "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10); expect=:badradiuscount)
        @test same_kind(() -> g_termpt(U, "SUN", "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
            () -> termpt(ctx, U, "SUN", "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10); expect=:badradiuscount)
        both_pdpool(ctx, "BODY301_RADII", radii)
        # A source that touches the sphere about the target.
        sun = gdpool(ctx, "BODY10_RADII")
        both_pdpool(ctx, "BODY10_RADII", [1.6e8, 1.6e8, 1.6e8])
        for m in (U, "PENUMBRAL/TANGENT/ELLIPSOID", UD)
            @test same_kind(() -> g_termpt(m, "SUN", "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10),
                () -> termpt(ctx, m, "SUN", "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 10))
        end
        both_pdpool(ctx, "BODY10_RADII", sun)
    end

    @testset "room in the outputs" begin
        ks = Allspice.kernelset(ctx)
        et = 1.0e8
        z = [0.0, 0.0, 1.0]
        # `limbpt_c` does not count the room that is left (it writes past `maxn` points on
        # DSK data); Allspice throws `:outofroom`. With one point per half plane and
        # `ncuts <= maxn` the two agree.
        a = limbpt(ctx, "TANGENT/ELLIPSOID", "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 10, 1e-4, 1e-7, 10)
        @test sum(a[1]) == 10
        # `termpt_c` takes the number of points of the half plane from the room one time
        # per point, so the room goes down by the square: 3 half planes of 1 point need 3.
        @test c_kind(() -> g_termpt("UMBRAL/TANGENT/ELLIPSOID", "SUN", "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 3)) === nothing
        @test sum(termpt(ctx, "UMBRAL/TANGENT/ELLIPSOID", "SUN", "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 3, 1e-4, 1e-7, 3)[1]) == 3
        # No cutting half plane is a valid request of `termpt`.
        @test c_kind(() -> g_termpt("UMBRAL/TANGENT/ELLIPSOID", "SUN", "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 0, 1e-4, 1e-7, 3)) === nothing
        @test isempty(termpt(ctx, "UMBRAL/TANGENT/ELLIPSOID", "SUN", "MOON", et, "IAU_MOON", "LT", "CENTER", "EARTH", z, 0.1, 0, 1e-4, 1e-7, 3)[2])
    end

    @testset "tangpt: no allocation and inferred type" begin
        et = 2.5e8
        dv = SVector{3,Float64}(spkpos(ctx, 301, et, 1, "LT+S", 399)[1]) + SVector(3000.0, -2000.0, 1000.0)
        for c in GEO_CORRS, loc in ("TANGENT POINT", "SURFACE POINT")
            @test gwarm(lt_alloc_tang, ks, et, c, loc, dv) == 0
            @test gwarm(lt_alloc_tang, ks, et, c, loc, -dv) == 0
        end
        @test (@inferred tangpt(ks, "ELLIPSOID", 301, et, 10020, "CN+S", "TANGENT POINT", 399, 1, dv)) isa
              Tuple{SVector{3,Float64},Float64,Float64,SVector{3,Float64},Float64,SVector{3,Float64}}
        @test (@inferred limbpt(ks, "TANGENT/ELLIPSOID", 301, et, 10020, "LT", "CENTER", 399, dv, 0.1, 3, 1e-4, 1e-7, 10)) isa
              Tuple{Vector{Int},Vector{SVector{3,Float64}},Vector{Float64},Vector{SVector{3,Float64}}}
        @test (@inferred termpt(ks, "UMBRAL/TANGENT/ELLIPSOID", 10, 301, et, 10020, "LT", "CENTER", 399, dv, 0.1, 3, 1e-4, 1e-7, 10)) isa
              Tuple{Vector{Int},Vector{SVector{3,Float64}},Vector{Float64},Vector{SVector{3,Float64}}}
    end
end
