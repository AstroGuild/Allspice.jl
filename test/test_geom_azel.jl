# Geometry group: azimuth and elevation (`azlrec`, `recazl`, `dazldr`, `drdazl`,
# `azlcpo`) against CSPICE.
#
# Error measures: for the coordinate calls, the absolute difference of each output over
# the size of the largest output of the call; for `azlcpo`, the difference of the range
# over the range, of the angles in radians, and of each rate over the size of the rates
# (range rate over the speed, angle rates over the largest angle rate).

isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")
import ForwardDiff

# Measured: 2e-16 (coordinates and Jacobians), 3e-16 (azlcpo range), 2e-15 rad (angles),
# 4e-14 (rates).
const AZ_TOL = (coord=5e-15, range=5e-12, angle=5e-11, rate=5e-10)

az_alloc_cpo(ks, et, c, p, ccw, plsz) = @allocated azlcpo(ks, "ELLIPSOID", 301, et, c, ccw, plsz, p, 399, 10013)
az_alloc_rec(a, b, c) = @allocated azlrec(a, b, c, true, false)
az_alloc_azl(v) = @allocated recazl(v, false, true)
az_alloc_dr(a, b, c) = @allocated dazldr(a, b, c, false, true)
az_alloc_rd(a, b, c) = @allocated drdazl(a, b, c, true, false)

@gf_light group("geom_azel: azlrec, recazl, dazldr, drdazl, azlcpo", Kernels.FULL_SET_FILES...) do
    ctx = geo_setup()
    ks = Allspice.kernelset(ctx)
    rng = MersenneTwister(SEED + 104)

    @testset "coordinates and Jacobians, 1000 random cases" begin
        st = GeoStats()
        for k in 1:1000
            ccw = rand(rng, Bool)
            plsz = rand(rng, Bool)
            range = 10.0^(6 * rand(rng) - 1)
            az = 2pi * rand(rng)
            el = pi * (rand(rng) - 0.5)
            a = azlrec(range, az, el, ccw, plsz)
            b = g_azlrec(range, az, el, ccw, plsz)
            gput!(st, "azlrec", grel(a, b))
            v = range .* normalize(randn(rng, 3))
            a = recazl(v, ccw, plsz)
            b = g_recazl(v, ccw, plsz)
            gput!(st, "recazl", max(abs(a[1] - b[1]) / b[1], abs(a[2] - b[2]), abs(a[3] - b[3])))
            @test 0.0 <= a[2] < 2pi && abs(a[3]) <= pi / 2
            gput!(st, "azlrec of recazl", grel(azlrec(a..., ccw, plsz), v))
            a = dazldr(v..., ccw, plsz)
            b = g_dazldr(v..., ccw, plsz)
            gput!(st, "dazldr", maximum(abs.(a .- b)) / maximum(abs.(b)))
            a = drdazl(range, az, el, ccw, plsz)
            b = g_drdazl(range, az, el, ccw, plsz)
            gput!(st, "drdazl", maximum(abs.(a .- b)) / maximum(abs.(b)))
            if k % 20 == 0
                @test azlrec(range, az, el; azccw=ccw, elplsz=plsz) == azlrec(range, az, el, ccw, plsz)
                @test recazl(v; azccw=ccw, elplsz=plsz) == recazl(v, ccw, plsz)
                @test dazldr(v...; azccw=ccw, elplsz=plsz) == dazldr(v..., ccw, plsz)
                @test drdazl(range, az, el; azccw=ccw, elplsz=plsz) == drdazl(range, az, el, ccw, plsz)
                # The two Jacobians are inverses, and they are the derivatives.
                p = azlrec(range, az, el, ccw, plsz)
                @test isapprox(dazldr(p..., ccw, plsz) * drdazl(range, az, el, ccw, plsz), I; atol=1e-9)
                jf = ForwardDiff.jacobian(x -> azlrec(x[1], x[2], x[3], ccw, plsz), SVector(range, az, el))
                @test isapprox(jf, drdazl(range, az, el, ccw, plsz); rtol=1e-12)
            end
        end
        # Special points: the axes, zero angles with each sense, the origin.
        for v in ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 0.0, -2.0], [-1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0], [0.0, 0.0, 0.0], [1.0, -0.0, 0.0], [1.0, 1e-300, 0.0]), ccw in (true, false),
            plsz in (true, false)

            a = recazl(v, ccw, plsz)
            b = g_recazl(v, ccw, plsz)
            @test a[1] == b[1] && a[2] == b[2] && a[3] == b[3]
            @test collect(azlrec(a..., ccw, plsz)) == g_azlrec(b..., ccw, plsz)
        end
        greport("azimuth and elevation", st)
        for k in keys(st.d)
            @test st.d[k] < AZ_TOL.coord
        end
        @test same_kind(() -> g_dazldr(0.0, 0.0, 1.0, true, true), () -> dazldr(0.0, 0.0, 1.0, true, true); expect=:pointonzaxis)
        @test same_kind(() -> g_drdazl(-1.0, 0.3, 0.2, true, true), () -> drdazl(-1.0, 0.3, 0.2, true, true); expect=:valueoutofrange)
    end

    @testset "azlcpo, 600 random cases" begin
        st = GeoStats()
        for k in 1:600
            ctr = rand(rng, ["EARTH", "MOON", "VENUS"])
            t = rand(rng, setdiff(GEO_BODIES, [ctr]))
            cor = rand(rng, GEO_CORRS)
            et = geo_epoch(rng)
            cref = rand(rng, ctr == "EARTH" ? ["IAU_EARTH", "GEO_SITE"] : GEO_FIXED[ctr])
            rad = SPICE.bodvrd(ctr, "RADII", 3)
            lon = 2pi * rand(rng)
            lat = pi * (rand(rng) - 0.5)
            pos = (1 + 0.1 * rand(rng)) .* [rad[1] * cos(lon) * cos(lat), rad[2] * sin(lon) * cos(lat), rad[3] * sin(lat)]
            k % 50 == 0 && (pos = [0.0, 0.0, rand(rng, (-1, 1)) * rad[3]])     # on the polar axis
            ccw = rand(rng, Bool)
            plsz = rand(rng, Bool)
            a = azlcpo(ctx, "ELLIPSOID", t, et, cor, ccw, plsz, pos, ctr, cref)
            b = g_azlcpo("ELLIPSOID", t, et, cor, ccw, plsz, pos, ctr, cref)
            gput!(st, "range", abs(a[1][1] - b[1][1]) / b[1][1])
            gput!(st, "azimuth and elevation (rad)", max(abs(a[1][2] - b[1][2]), abs(a[1][3] - b[1][3])))
            speed = norm(g_spkcpo(t, et, cref, "OBSERVER", cor, pos, ctr, cref)[1][4:6])
            gput!(st, "range rate over the speed", abs(a[1][4] - b[1][4]) / speed)
            gput!(st, "angle rates over the largest", max(abs(a[1][5] - b[1][5]), abs(a[1][6] - b[1][6])) /
                                                       max(abs(b[1][5]), abs(b[1][6])))
            gput!(st, "light time", abs(a[2] - b[2]) / b[2])
            if k % 10 == 0
                c = azlcpo(ctx, "ellipsoid", GEO_IDS[t], et, Aberration(cor), ccw, plsz, SVector{3}(pos), GEO_IDS[ctr], namfrm(ctx, cref))
                @test c == a
            end
        end
        greport("azlcpo", st)
        @test gget(st, "range") < AZ_TOL.range
        @test gget(st, "azimuth and elevation (rad)") < AZ_TOL.angle
        @test gget(st, "range rate over the speed") < AZ_TOL.rate
        @test gget(st, "angle rates over the largest") < AZ_TOL.rate
        @test gget(st, "light time") < AZ_TOL.range
    end

    @testset "azlcpo: error conditions, kinds of CSPICE" begin
        et = 1.0e8
        p = [6378.0, 0.0, 0.0]
        for (m, t, cor, ctr, cref) in (
                ("ELLIPSOID", "MOON", "FOO", "EARTH", "IAU_EARTH"),
                ("ELLIPSOID", "NOBODY", "LT", "EARTH", "IAU_EARTH"),
                ("ELLIPSOID", "MOON", "LT", "NOBODY", "IAU_EARTH"),
                ("ELLIPSOID", "MOON", "LT", "EARTH", "NOFRAME"),
                ("ELLIPSOID", "MOON", "LT", "EARTH", "IAU_MOON"),
                ("ELLIPSOID", "MOON", "LT", "EARTH", "J2000"),
                ("DSK/UNPRIORITIZED", "MOON", "LT", "EARTH", "IAU_EARTH"),
                ("SPHERE", "MOON", "LT", "EARTH", "IAU_EARTH"))
            @test same_kind(() -> g_azlcpo(m, t, et, cor, true, true, p, ctr, cref),
                () -> azlcpo(ctx, m, t, et, cor, true, true, p, ctr, cref))
        end
        @test neo_kind(() -> azlcpo(ctx, "ELLIPSOID", "MOON", et, "LT", true, true, [1.0], "EARTH", "IAU_EARTH")) == :badarg
    end

    @testset "no allocation and inferred types" begin
        et = 2.5e8
        p = SVector(6378.0 * cos(0.5), 0.0, 6378.0 * sin(0.5))
        for c in GEO_CORRS, ccw in (true, false)
            @test gwarm(az_alloc_cpo, ks, et, c, p, ccw, !ccw) == 0
        end
        @test gwarm(az_alloc_cpo, ks, et, "LT+S", SVector(0.0, 0.0, 6356.0), true, true) == 0
        @test gwarm(az_alloc_rec, 2.0, 0.3, 0.4) == 0
        @test gwarm(az_alloc_azl, p) == 0
        @test gwarm(az_alloc_dr, 1.0, 2.0, 3.0) == 0
        @test gwarm(az_alloc_rd, 2.0, 0.3, 0.4) == 0
        @test (@inferred azlcpo(ks, "ELLIPSOID", 301, et, "CN+S", true, false, p, 399, 10013)) isa Tuple{SVector{6,Float64},Float64}
        @test (@inferred azlrec(2.0, 0.3, 0.4, true, false)) isa SVector{3,Float64}
        @test (@inferred recazl(p, true, false)) isa NTuple{3,Float64}
        @test (@inferred dazldr(1.0, 2.0, 3.0, true, false)) isa SMatrix{3,3,Float64,9}
        @test (@inferred drdazl(2.0, 0.3, 0.4, true, false)) isa SMatrix{3,3,Float64,9}
    end

    @testset "azlcpo: epoch as a dual number" begin
        st = GeoStats()
        p = SVector(6378.0 * cos(0.5), 0.0, 6378.0 * sin(0.5))
        for et in record_safe_epochs(rng, nsample(6, 24)), cor in GEO_CORRS
            f = t -> azlcpo(ks, "ELLIPSOID", 301, t, cor, true, true, p, 399, 10013)[1]
            d = f(ForwardDiff.Dual(et, 1.0))
            v = f(et)
            @test isapprox(ForwardDiff.value.(d), v; rtol=1e-12)
            dp = ForwardDiff.partials.(d, 1)[SOneTo(3)]
            fd = geo_fd(t -> f(t)[SOneTo(3)], et, 20.0)
            w = SVector(1.0 / v[1], 1.0, 1.0)          # range in units of itself, angles in radians
            gput!(st, "d(range, az, el)/dt against central difference", norm(w .* (dp .- fd)) / norm(w .* fd))
            gput!(st, "d(range, az, el)/dt against the returned rates", norm(w .* (dp .- v[SVector(4, 5, 6)])) / norm(w .* v[SVector(4, 5, 6)]))
        end
        greport("azlcpo, dual epoch", st)
        @test gget(st, "d(range, az, el)/dt against central difference") < 1e-4
        @test gget(st, "d(range, az, el)/dt against the returned rates") < 1e-3
    end
end
