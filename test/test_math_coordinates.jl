# Group 14: coordinate conversions, their Jacobians, state conversion between coordinate
# systems and unit conversion against CSPICE.
#
# Feature file: features/math_coordinates.feature
#
# The random points have magnitudes from 1e-10 to 1e10. The degenerate inputs are the
# origin, the points on the z axis (both poles) and in the x-y plane, points on the
# branch cut of the longitude (negative x axis), and for the geodetic calls points on,
# inside and far outside the spheroid, oblate and prolate.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff
isdefined(@__MODULE__, :CS) || include("mathoracle.jl")

"Points on the axes, in the coordinate planes and on the longitude branch cut."
const SPECIAL_POINTS = (SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -3.0),
    SVector(1.0, 0.0, 0.0), SVector(-2.0, 0.0, 0.0), SVector(0.0, 5.0, 0.0), SVector(0.0, -5.0, 0.0),
    SVector(-1.0, 0.0, 1.0), SVector(-1.0, -0.0, 0.0), SVector(1.0, 1.0, 0.0), SVector(0.0, 1.0, 1.0),
    SVector(-1.0, 1e-300, 0.0), SVector(-1.0, -1e-300, 0.0), SVector(6378.0, 0.0, 0.0),
    SVector(0.0, 0.0, 6356.0))

group("14. math: coordinates", Kernels.PCK11) do
    rng = MersenneTwister(SEED + 14)
    n = 400
    t = Tally()
    oracle_load((Kernels.PCK11,))
    ctx = neo_context((Kernels.PCK11,))
    rpt() = rand(rng) < 0.2 ? rand(rng, SPECIAL_POINTS) : rvec(rng; special=0.0)
    rang() = rand(rng) < 0.15 ? rand(rng, (0.0, pi, -pi, pi / 2, -pi / 2, 2pi)) : 4pi * rand(rng) - 2pi
    rrad() = rand(rng) < 0.05 ? 0.0 : rmag(rng)

    @apart @testset "rectangular to and from the angular systems" begin
        for _ in 1:n
            p = rpt()
            r, a1, a2, z = rrad(), rang(), rang(), rmag(rng) * randn(rng)
            for (name, neo, c) in ((:reclat, reclat, CS.reclat), (:recrad, recrad, CS.recrad),
                (:recsph, recsph, CS.recsph), (:reccyl, reccyl, CS.reccyl))
                versus!(t, name, () -> neo(p), () -> c(p))
            end
            for (name, neo, c) in ((:latrec, latrec, CS.latrec), (:radrec, radrec, CS.radrec),
                (:sphrec, sphrec, CS.sphrec), (:cyllat, cyllat, CS.cyllat),
                (:latcyl, latcyl, CS.latcyl), (:cylsph, cylsph, CS.cylsph),
                (:sphcyl, sphcyl, CS.sphcyl), (:sphlat, sphlat, CS.sphlat), (:latsph, latsph, CS.latsph))
                versus!(t, name, () -> neo(r, a1, a2), () -> c(r, a1, a2))
            end
            versus!(t, :cylrec, () -> cylrec(r, a1, z), () -> CS.cylrec(r, a1, z))
            versus!(t, :cyllat, () -> cyllat(r, a1, z), () -> CS.cyllat(r, a1, z))
            versus!(t, :cylsph, () -> cylsph(r, a1, z), () -> CS.cylsph(r, a1, z))
            # round trips
            tol = 4eps() * max(vnorm(p), floatmin())
            @test vdist(latrec(reclat(p)...), p) <= tol
            @test vdist(radrec(recrad(p)...), p) <= tol
            @test vdist(sphrec(recsph(p)...), p) <= tol
            @test vdist(cylrec(reccyl(p)...), p) <= tol
            @test vdist(latrec(sphlat(recsph(p)...)...), p) <= 2tol
            @test vdist(cylrec(latcyl(reclat(p)...)...), p) <= 2tol
            @test vdist(sphrec(cylsph(reccyl(p)...)...), p) <= 2tol
            @test vdist(latrec(cyllat(reccyl(p)...)...), p) <= 2tol
            @test vdist(cylrec(sphcyl(recsph(p)...)...), p) <= 2tol
            @test vdist(sphrec(latsph(reclat(p)...)...), p) <= 2tol
        end
        # ranges and branch cuts
        @test reclat([0.0, 0.0, 0.0]) == (0.0, 0.0, 0.0)
        @test reclat([0.0, 0.0, 2.0]) == (2.0, 0.0, pi / 2)
        @test reclat([0.0, 0.0, -2.0]) == (2.0, 0.0, -pi / 2)
        @test reclat([-1.0, 0.0, 0.0])[2] == Float64(pi)
        @test reclat([-1.0, -0.0, 0.0])[2] == -Float64(pi) == CS.reclat([-1.0, -0.0, 0.0])[2]
        @test recrad([1.0, -1.0, 0.0])[2] ≈ 7pi / 4
        @test 0.0 <= recrad([-1.0, -1e-300, 0.0])[2] <= 2pi
        @test recsph([0.0, 0.0, -1.0]) == (1.0, Float64(pi), 0.0)
        @test reccyl([0.0, 0.0, 4.0]) == (0.0, 0.0, 4.0)
        @test reccyl([0.0, -1.0, 4.0])[2] ≈ 3pi / 2
    end

    @apart @testset "geodetic and planetographic" begin
        shapes = ((6378.1366, 1 / 298.25642), (3396.19, 0.00589), (1.0, 0.0), (10.0, 0.9), (2.0, -1.5),
            (1e-3, 0.3), (1e8, 0.5))
        bodies = ("EARTH", "MARS", "MOON", "SUN", "JUPITER", "VENUS", "399", "PHOBOS", "NOSUCHBODY", "MERCURY")
        for _ in 1:n
            re, f = rand(rng, shapes)
            rand(rng) < 0.03 && (re = -re)
            rand(rng) < 0.03 && (f = 1.0 + rand(rng))
            # points near the spheroid, far from it, deep inside, and on the axes
            u = rand(rng)
            p = u < 0.15 ? abs(re) * rand(rng, SPECIAL_POINTS) / 6378.0 :
                u < 0.5 ? abs(re) * (0.5 + rand(rng)) * SVector(randn(rng), randn(rng), randn(rng)) :
                u < 0.6 ? abs(re) * 1e-6 * SVector(randn(rng), randn(rng), randn(rng)) : abs(re) * rmag(rng) * 1e-5 * SVector(randn(rng), randn(rng), randn(rng))
            lon, lat = rang(), (rand(rng) < 0.1 ? rand(rng, (pi / 2, -pi / 2, 0.0)) : pi * rand(rng) - pi / 2)
            alt = abs(re) * (rand(rng) < 0.1 ? 0.0 : randn(rng))
            versus!(t, :recgeo, () -> recgeo(p, re, f), () -> CS.recgeo(p, re, f))
            versus!(t, :georec, () -> georec(lon, lat, alt, re, f), () -> CS.georec(lon, lat, alt, re, f))
            versus!(t, :drdgeo, () -> drdgeo(lon, lat, alt, re, f), () -> CS.drdgeo(lon, lat, alt, re, f))
            # deep inside the spheroid the Jacobian has entries of the size of the
            # distance from the center, which carry the rounding error of the radius
            deep = vnorm(p) < 1e-2 * abs(re)
            versus!(t, deep ? :dgeodr_deep : :dgeodr, () -> dgeodr(p..., re, f), () -> CS.dgeodr(p..., re, f))
            b = rand(rng, bodies)
            versus!(t, :recpgr, () -> recpgr(ctx, b, p, re, f), () -> CS.recpgr(b, p, re, f))
            versus!(t, :pgrrec, () -> pgrrec(ctx, b, lon, lat, alt, re, f), () -> CS.pgrrec(b, lon, lat, alt, re, f))
            versus!(t, :drdpgr, () -> drdpgr(ctx, b, lon, lat, alt, re, f), () -> CS.drdpgr(b, lon, lat, alt, re, f))
            versus!(t, deep ? :dpgrdr_deep : :dpgrdr, () -> dpgrdr(ctx, b, p..., re, f), () -> CS.dpgrdr(b, p..., re, f))
            if re > 0 && f < 1 && abs(f) < 0.95
                g = recgeo(p, re, f)
                @test vdist(georec(g..., re, f), p) <= 1e-9 * max(vnorm(p), re)
            end
        end
        @test throws_kind(() -> recgeo([1.0, 2.0, 3.0], -1.0, 0.1), :valueoutofrange)
        @test throws_kind(() -> georec(0.1, 0.2, 3.0, 1.0, 1.0), :valueoutofrange)
        @test throws_kind(() -> drdgeo(0.1, 0.2, 3.0, -1.0, 0.5), :badradius)
        @test throws_kind(() -> dgeodr(0.0, 0.0, 3.0, 1.0, 0.5), :pointonzaxis)
        @test throws_kind(() -> recpgr(ctx, "NOSUCHBODY", [1.0, 2.0, 3.0], 1.0, 0.1), :idcodenotfound)
        # the longitude sense: west for prograde Mars, east for the Earth and the Moon
        p = [100.0, 3000.0, 500.0]
        @test recpgr(ctx, "MARS", p, 3396.19, 0.00589)[1] ≈ 2pi - recgeo(p, 3396.19, 0.00589)[1]
        @test recpgr(ctx, "EARTH", p, 6378.0, 0.003)[1] == recgeo(p, 6378.0, 0.003)[1]
        @test recpgr(ctx, 301, p, 1737.4, 0.0)[1] == recgeo(p, 1737.4, 0.0)[1]
        @test recpgr(ctx, "SUN", p, 695700.0, 0.0)[1] == recgeo(p, 695700.0, 0.0)[1]
        @test recpgr(ctx, "VENUS", p, 6051.8, 0.0)[1] == recgeo(p, 6051.8, 0.0)[1]      # retrograde
        # the pool variable BODY<id>_PGR_POSITIVE_LON overrides the rotation sense
        tk = joinpath(mktempdir(), "pgr.tk")
        write(tk, "\\begindata\nBODY499_PGR_POSITIVE_LON = 'EAST'\nBODY599_PGR_POSITIVE_LON = 'NORTH'\n\\begintext\n")
        SPICE.furnsh(tk)
        ctx2 = neo_context((Kernels.PCK11, tk))
        versus!(t, :recpgr, () -> recpgr(ctx2, "MARS", p, 3396.19, 0.00589), () -> CS.recpgr("MARS", p, 3396.19, 0.00589))
        versus!(t, :recpgr, () -> recpgr(ctx2, "JUPITER", p, 71492.0, 0.06), () -> CS.recpgr("JUPITER", p, 71492.0, 0.06))
        @test recpgr(ctx2, "MARS", p, 3396.19, 0.00589)[1] == recgeo(p, 3396.19, 0.00589)[1]
        @test throws_kind(() -> pgrrec(ctx2, "JUPITER", 0.1, 0.2, 3.0, 71492.0, 0.06), :invalidoption)
        SPICE.unload(tk)
    end

    @apart @testset "Jacobians" begin
        for _ in 1:n
            p = rpt()
            r, a1, a2 = rrad(), rang(), rang()
            versus!(t, :drdlat, () -> drdlat(r, a1, a2), () -> CS.drdlat(r, a1, a2))
            versus!(t, :drdsph, () -> drdsph(r, a1, a2), () -> CS.drdsph(r, a1, a2))
            versus!(t, :drdcyl, () -> drdcyl(r, a1, a2), () -> CS.drdcyl(r, a1, a2))
            versus!(t, :dlatdr, () -> dlatdr(p...), () -> CS.dlatdr(p...))
            versus!(t, :dsphdr, () -> dsphdr(p...), () -> CS.dsphdr(p...))
            versus!(t, :dcyldr, () -> dcyldr(p...), () -> CS.dcyldr(p...))
        end
        # against automatic differentiation of the conversions
        q = SVector(1.3, -0.4, 0.9)
        x = SVector(3100.0, -4200.0, 2500.0)
        re, f = 6378.1366, 1 / 298.25642
        @test drdlat(q...) ≈ ForwardDiff.jacobian(v -> latrec(v...), q) rtol = 1e-14
        @test drdsph(q...) ≈ ForwardDiff.jacobian(v -> sphrec(v...), q) rtol = 1e-14
        @test drdcyl(q...) ≈ ForwardDiff.jacobian(v -> cylrec(v...), q) rtol = 1e-14
        @test drdgeo(q[1], q[2], 300.0, re, f) ≈
              ForwardDiff.jacobian(v -> georec(v[1], v[2], v[3], re, f), SVector(q[1], q[2], 300.0)) rtol = 1e-13
        @test dlatdr(x...) ≈ ForwardDiff.jacobian(v -> SVector(reclat(v)), x) rtol = 1e-13
        @test dsphdr(x...) ≈ ForwardDiff.jacobian(v -> SVector(recsph(v)), x) rtol = 1e-13
        @test dcyldr(x...) ≈ ForwardDiff.jacobian(v -> SVector(reccyl(v)), x) rtol = 1e-13
        # `recgeo` runs the iteration of `nearpt`; its derivative comes from the
        # stationarity condition of the near point
        @test dgeodr(x..., re, f) ≈ ForwardDiff.jacobian(v -> SVector(recgeo(v, re, f)), x) rtol = 1e-10
        @test dpgrdr(ctx, "MARS", x..., 3396.19, 0.00589) ≈
              ForwardDiff.jacobian(v -> SVector(recpgr(ctx, "MARS", v, 3396.19, 0.00589)), x) rtol = 1e-10
        @test drdpgr(ctx, "MARS", q[1], q[2], 300.0, 3396.19, 0.00589) ≈
              ForwardDiff.jacobian(v -> pgrrec(ctx, "MARS", v[1], v[2], v[3], 3396.19, 0.00589),
            SVector(q[1], q[2], 300.0)) rtol = 1e-13
        @test dlatdr(x...) * drdlat(reclat(x)...) ≈ I atol = 1e-12
        @test throws_kind(() -> dlatdr(0.0, 0.0, 1.0), :pointonzaxis)
        @test throws_kind(() -> dsphdr(0.0, 0.0, 1.0), :pointonzaxis)
        @test throws_kind(() -> dcyldr(0.0, 0.0, 1.0), :pointonzaxis)
        # a Taylor series passes through the conversions that need sine and cosine only
        tt = Allspice.taylor_var(0.0, Val(1))
        v = latrec(2.0 + 0.1 * tt, 0.3 + 0.2 * tt, -0.4 + 0.5 * tt)
        @test map(c -> Allspice.derivatives(c)[2], v) ≈ drdlat(2.0, 0.3, -0.4) * [0.1, 0.2, 0.5] rtol = 1e-14
        @test Allspice.derivatives(sphrec(2.0, 0.3 + tt, 0.4)[3])[2] ≈ -2.0 * sin(0.3)
        @test Allspice.derivatives(cylrec(2.0, 0.3 + tt, 0.4)[1])[2] ≈ -2.0 * sin(0.3)
    end

    @apart @testset "state between coordinate systems" begin
        systems = ("RECTANGULAR", "CYLINDRICAL", "LATITUDINAL", "SPHERICAL", "GEODETIC", "PLANETOGRAPHIC")
        bodies = ("EARTH", "MARS", "MOON", "NOSUCHBODY", "JUPITER")
        for _ in 1:2n
            isys, osys = rand(rng, systems), rand(rng, systems)
            rand(rng) < 0.03 && (osys = "POLAR")
            b = rand(rng, bodies)
            vel = rand(rng) < 0.1 ? SVector(0.0, 0.0, randn(rng)) : SVector(randn(rng), randn(rng), randn(rng))
            s = if isys == "RECTANGULAR"
                vcat(rand(rng) < 0.15 ? rand(rng, SPECIAL_POINTS) : 7000.0 * SVector(randn(rng), randn(rng), randn(rng)), vel)
            elseif isys == "CYLINDRICAL"
                SVector(7000.0 * rand(rng), rang(), 7000.0 * randn(rng), randn(rng), 1e-3 * randn(rng), randn(rng))
            elseif isys == "GEODETIC" || isys == "PLANETOGRAPHIC"
                SVector(rang(), pi * rand(rng) - pi / 2, 500.0 * randn(rng), 1e-3 * randn(rng), 1e-3 * randn(rng), randn(rng))
            else
                SVector(7000.0 * rand(rng), rang(), rang() / 2, randn(rng), 1e-3 * randn(rng), 1e-3 * randn(rng))
            end
            # next to the z axis the rate of the longitude is the quotient of two
            # rounding errors, in CSPICE as here
            rect = outcome(() -> xfmsta(ctx, s, isys, "RECTANGULAR", "EARTH"))
            pole = !(rect isa Symbol) && hypot(rect[1], rect[2]) < 1e-6 * hypot(rect[1], rect[2], rect[3])
            versus!(t, pole ? :xfmsta_pole : :xfmsta, () -> xfmsta(ctx, s, isys, osys, b), () -> CS.xfmsta(s, isys, osys, b))
        end
        s = SVector(3100.0, -4200.0, 2500.0, 1.0, 2.0, -3.0)
        for sys in systems[2:end]
            back = xfmsta(ctx, xfmsta(ctx, s, "rectangular", sys, "EARTH"), sys, "RECTANGULAR", "EARTH")
            @test back ≈ s rtol = 1e-11
        end
        # the rates are the time derivatives of the coordinates
        lat = xfmsta(ctx, s, "RECTANGULAR", "LATITUDINAL", "EARTH")
        d = ForwardDiff.derivative(τ -> SVector(reclat(SVector(s[1], s[2], s[3]) + τ * SVector(s[4], s[5], s[6]))), 0.0)
        @test lat[SVector(4, 5, 6)] ≈ d rtol = 1e-13
        @test xfmsta(ctx, s, " rectangular ", "RECTANGULAR", "anything") == s
        @test throws_kind(() -> xfmsta(ctx, s, "RECTANGULAR", "POLAR", "EARTH"), :coordsysnotrec)
        @test throws_kind(() -> xfmsta(ctx, s, "RECTANGULAR", "GEODETIC", "NOSUCHBODY"), :idcodenotfound)
        @test throws_kind(() -> xfmsta(ctx, [0.0, 0, 1, 1, 0, 0], "RECTANGULAR", "SPHERICAL", "EARTH"), :invalidstate)
        @test xfmsta(ctx, [0.0, 0, 2, 0, 0, 3], "RECTANGULAR", "LATITUDINAL", "EARTH") == [2.0, 0.0, pi / 2, 3.0, 0.0, 0.0]
    end

    @apart @testset "unit conversion" begin
        units = ("RADIANS", "DEGREES", "ARCMINUTES", "ARCSECONDS", "HOURANGLE", "MINUTEANGLE", "SECONDANGLE",
            "METERS", "KM", "CM", "MM", "LIGHTSECS", "AU", "M", "KILOMETERS", "CENTIMETERS", "MILLIMETERS",
            "FEET", "INCHES", "STATUTE_MILES", "NAUTICAL_MILES", "YARDS", "LIGHTYEARS", "PARSECS", "SECONDS",
            "MINUTES", "HOURS", "DAYS", "JULIAN_YEARS", "TROPICAL_YEARS", "YEARS", "FURLONGS", "km", "Degrees ")
        for a in units, b in units
            x = rmag(rng) * randn(rng)
            versus!(t, :convrt, () -> convrt(x, a, b), () -> CS.convrt(x, a, b))
        end
        @test convrt(180.0, "DEGREES", "RADIANS") == Float64(pi)
        @test convrt(1.0, "AU", "KM") == CS.convrt(1.0, "AU", "KM") ≈ 149597870.613689
        @test throws_kind(() -> convrt(1.0, "KM", "FURLONGS"), :unitsnotrec)
        @test throws_kind(() -> convrt(1.0, "KM", "SECONDS"), :incompatibleunits)
        @test ForwardDiff.derivative(x -> convrt(x, "HOURS", "SECONDS"), 2.0) == 3600.0
    end

    @apart @testset "differences against CSPICE" begin
        # The trigonometric functions of Julia and of the C library differ by one unit
        # of rounding for some arguments. A geodetic result has the conditioning of the
        # near point: the altitude of a point 1e5 radii away has the rounding error of
        # the distance, and the inverse Jacobians divide by the distance from the axis.
        e = 4eps()
        report!(t, "coordinates"; tol=e, tols=Dict(:recgeo => 1e-14, :recpgr => 1e-14, :georec => 1e-14,
            :pgrrec => 1e-14, :drdgeo => 1e-13, :drdpgr => 1e-13, :dgeodr => 1e-12, :dpgrdr => 1e-12,
            :dgeodr_deep => 1e-1, :dpgrdr_deep => 1e-1, :xfmsta_pole => 1e-3,
            :dlatdr => 1e-14, :dsphdr => 1e-14, :dcyldr => 1e-14, :xfmsta => 1e-11, :convrt => 0.0))
    end

    @apart @testset "inferred types and no allocation for static arguments" begin
        p = SVector(3100.0, -4200.0, 2500.0)
        re, f = 6378.1366, 1 / 298.25642
        T3, V3, M3 = NTuple{3,Float64}, SVector{3,Float64}, SMatrix{3,3,Float64,9}
        calls = Any[(reclat, (p,), T3), (recrad, (p,), T3), (recsph, (p,), T3), (reccyl, (p,), T3),
            (latrec, (1.0, 0.2, 0.3), V3), (radrec, (1.0, 0.2, 0.3), V3), (sphrec, (1.0, 0.2, 0.3), V3),
            (cylrec, (1.0, 0.2, 0.3), V3), (cyllat, (1.0, 0.2, 0.3), T3), (latcyl, (1.0, 0.2, 0.3), T3),
            (cylsph, (1.0, 0.2, 0.3), T3), (sphcyl, (1.0, 0.2, 0.3), T3), (sphlat, (1.0, 0.2, 0.3), T3),
            (latsph, (1.0, 0.2, 0.3), T3), (recgeo, (p, re, f), T3), (georec, (0.1, 0.2, 300.0, re, f), V3),
            (drdlat, (1.0, 0.2, 0.3), M3), (dlatdr, (p[1], p[2], p[3]), M3), (drdsph, (1.0, 0.2, 0.3), M3),
            (dsphdr, (p[1], p[2], p[3]), M3), (drdcyl, (1.0, 0.2, 0.3), M3), (dcyldr, (p[1], p[2], p[3]), M3),
            (drdgeo, (0.1, 0.2, 300.0, re, f), M3), (dgeodr, (p[1], p[2], p[3], re, f), M3),
            (convrt, (1.0, "KM", "AU"), Float64),
            (recpgr, (ctx, 499, p, re, f), T3), (pgrrec, (ctx, 499, 0.1, 0.2, 300.0, re, f), V3),
            (drdpgr, (ctx, 499, 0.1, 0.2, 300.0, re, f), M3), (dpgrdr, (ctx, 499, p[1], p[2], p[3], re, f), M3)]
        for (fn, args, T) in calls
            @test (@inferred fn(args...)) isa T
        end
        # the calls without a kernel pool do not allocate; the pool calls build the
        # name of a pool variable
        for (fn, args, T) in calls[1:24]
            @test alloc_of(fn, args...) == 0
        end
        s = SVector(3100.0, -4200.0, 2500.0, 1.0, 2.0, -3.0)
        @test (@inferred xfmsta(ctx, s, "RECTANGULAR", "LATITUDINAL", "EARTH")) isa SVector{6,Float64}
    end
    SPICE.kclear()
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
