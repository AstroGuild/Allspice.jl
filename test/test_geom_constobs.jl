# Geometry group: states relative to an observer or a target with a constant position or
# a constant velocity in a frame (`spkcpo`, `spkcpt`, `spkcvo`, `spkcvt`) against CSPICE.
#
# Feature file: features/geometry_fixed_observer.feature
#
# Error measure: the norm of the difference from the CSPICE vector over the norm of the
# CSPICE vector, for the position and the velocity apart, and the relative difference
# of the light time. The differences are those of the aberration group: a corrected
# state is the difference of two barycentric vectors, and one unit in the last place of
# a barycentric coordinate is 1e-13 of the Earth-Moon distance.

isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")
import ForwardDiff

# Measured: 8e-14 (position), 4e-14 (velocity), 5e-14 (light time).
const CO_TOL = (pos=5e-12, vel=5e-10, lt=5e-12)
const CO_LOCI = ["OBSERVER", "TARGET", "CENTER"]

"A point near the surface of `body` in a body-fixed frame of it, and a slow velocity."
function co_site(rng, body)
    rad = SPICE.bodvrd(body, "RADII", 3)[1]
    lon = 2pi * rand(rng)
    lat = pi * (rand(rng) - 0.5)
    r = (rad + 500 * rand(rng)) .* [cos(lon) * cos(lat), sin(lon) * cos(lat), sin(lat)]
    return r, 0.02 .* randn(rng, 3)
end

co_alloc_cpo(ks, et, c, p) = @allocated spkcpo(ks, 301, et, 10013, "OBSERVER", c, p, 399, 10013)
co_alloc_cpt(ks, et, c, p) = @allocated spkcpt(ks, p, 399, 10013, et, 1, "TARGET", c, 301)
co_alloc_cvo(ks, et, c, s) = @allocated spkcvo(ks, 301, et, 10020, "CENTER", c, s, 1.0e8, 399, 10013)
co_alloc_cvt(ks, et, c, s) = @allocated spkcvt(ks, s, 1.0e8, 399, 1, et, 10013, "OBSERVER", c, 301)

group("geom_constobs: spkcpo, spkcpt, spkcvo, spkcvt", Kernels.FULL_SET_FILES...) do
    ctx = geo_setup()
    ks = NeoSpice.kernelset(ctx)
    rng = MersenneTwister(SEED + 103)

    @testset "600 random cases per call" begin
        st = GeoStats()
        function cmp!(name, a, b)
            gput!(st, "$name position", grel(a[1][1:3], b[1][1:3]))
            gput!(st, "$name velocity", grel(a[1][4:6], b[1][4:6]))
            gput!(st, "$name light time", abs(a[2] - b[2]) / b[2])
        end
        for k in 1:600
            ctr = rand(rng, ["EARTH", "MOON", "VENUS"])
            other = rand(rng, setdiff(GEO_BODIES, [ctr]))
            cor = rand(rng, GEO_CORRS)
            loc = rand(rng, CO_LOCI)
            et = geo_epoch(rng)
            cref = rand(rng, ctr == "EARTH" ? ["IAU_EARTH", "GEO_SITE"] : GEO_FIXED[ctr])
            outref = rand(rng, ["J2000", "ECLIPJ2000", cref, rand(rng, GEO_FIXED[other]), "IAU_EARTH", "GEO_SITE"])
            pos, vel = co_site(rng, ctr)
            # A constant position.
            a = spkcpo(ctx, other, et, outref, loc, cor, pos, ctr, cref)
            b = g_spkcpo(other, et, outref, loc, cor, pos, ctr, cref)
            cmp!("spkcpo", a, b)
            a = spkcpt(ctx, pos, ctr, cref, et, outref, loc, cor, other)
            b = g_spkcpt(pos, ctr, cref, et, outref, loc, cor, other)
            cmp!("spkcpt", a, b)
            if k % 10 == 0
                c = spkcpo(ctx, GEO_IDS[other], et, namfrm(ctx, outref), lowercase(loc), Aberration(cor),
                    SVector{3}(pos), GEO_IDS[ctr], namfrm(ctx, cref))
                @test c == spkcpo(ctx, other, et, outref, loc, cor, pos, ctr, cref)
            end
            # A constant velocity, in a body-fixed frame or in an inertial frame.
            sref = rand(rng, [cref, "J2000", "ECLIPJ2000"])
            sta = vcat(pos, sref == cref ? vel : 100 .* vel)
            epc = et + 2000 * (rand(rng) - 0.5)
            a = spkcvo(ctx, other, et, outref, loc, cor, sta, epc, ctr, sref)
            b = g_spkcvo(other, et, outref, loc, cor, sta, epc, ctr, sref)
            cmp!("spkcvo", a, b)
            a = spkcvt(ctx, sta, epc, ctr, sref, et, outref, loc, cor, other)
            b = g_spkcvt(sta, epc, ctr, sref, et, outref, loc, cor, other)
            cmp!("spkcvt", a, b)
        end
        greport("constant position and velocity", st)
        for name in ("spkcpo", "spkcpt", "spkcvo", "spkcvt")
            @test gget(st, "$name position") < CO_TOL.pos
            @test gget(st, "$name velocity") < CO_TOL.vel
            @test gget(st, "$name light time") < CO_TOL.lt
        end
    end

    @testset "error conditions, kinds of CSPICE" begin
        et = 1.0e8
        p = [6378.0, 0.0, 0.0]
        s = [6378.0, 0.0, 0.0, 0.0, 0.01, 0.0]
        for (t, outref, loc, cor, ctr, cref) in (
                ("MOON", "J2000", "OBSERVER", "FOO", "EARTH", "IAU_EARTH"),
                ("MOON", "J2000", "OBSERVER", "S", "EARTH", "IAU_EARTH"),
                ("NOBODY", "J2000", "OBSERVER", "LT", "EARTH", "IAU_EARTH"),
                ("MOON", "J2000", "OBSERVER", "LT", "NOBODY", "IAU_EARTH"),
                ("MOON", "NOFRAME", "OBSERVER", "LT", "EARTH", "IAU_EARTH"),
                ("MOON", "J2000", "OBSERVER", "LT", "EARTH", "NOFRAME"),
                ("MOON", "J2000", "SOMEWHERE", "LT", "EARTH", "IAU_EARTH"),
                ("MOON", "J2000", "", "LT", "EARTH", "IAU_EARTH"))
            @test same_kind(() -> g_spkcpo(t, et, outref, loc, cor, p, ctr, cref),
                () -> spkcpo(ctx, t, et, outref, loc, cor, p, ctr, cref))
            @test same_kind(() -> g_spkcpt(p, ctr, cref, et, outref, loc, cor, t),
                () -> spkcpt(ctx, p, ctr, cref, et, outref, loc, cor, t))
            @test same_kind(() -> g_spkcvo(t, et, outref, loc, cor, s, et, ctr, cref),
                () -> spkcvo(ctx, t, et, outref, loc, cor, s, et, ctr, cref))
            @test same_kind(() -> g_spkcvt(s, et, ctr, cref, et, outref, loc, cor, t),
                () -> spkcvt(ctx, s, et, ctr, cref, et, outref, loc, cor, t))
        end
        @test neo_kind(() -> spkcpo(ctx, "MOON", et, "J2000", "OBSERVER", "LT", [1.0, 2.0], "EARTH", "IAU_EARTH")) == :badarg
        # Inputs that CSPICE accepts: a target that is the center of the observer, and a
        # frame of the observer that is centered on another body.
        for (t, cref) in (("EARTH", "IAU_EARTH"), ("MOON", "IAU_VENUS"))
            a = spkcpo(ctx, t, et, "J2000", "OBSERVER", "LT+S", p, "EARTH", cref)
            b = g_spkcpo(t, et, "J2000", "OBSERVER", "LT+S", p, "EARTH", cref)
            @test grel(a[1][1:3], b[1][1:3]) < CO_TOL.pos && grel(a[1][4:6], b[1][4:6]) < CO_TOL.vel
        end
        @test neo_kind(() -> spkcvo(ctx, "MOON", et, "J2000", "OBSERVER", "LT", p, et, "EARTH", "IAU_EARTH")) == :badarg
    end

    @testset "no allocation and inferred types" begin
        et = 2.5e8
        p = SVector(6378.0 * cos(0.5), 0.0, 6378.0 * sin(0.5))
        s = SVector(p..., 0.0, 0.01, 0.0)
        for c in GEO_CORRS
            @test gwarm(co_alloc_cpo, ks, et, c, p) == 0
            @test gwarm(co_alloc_cpt, ks, et, c, p) == 0
            @test gwarm(co_alloc_cvo, ks, et, c, s) == 0
            @test gwarm(co_alloc_cvt, ks, et, c, s) == 0
        end
        T = Tuple{SVector{6,Float64},Float64}
        @test (@inferred spkcpo(ks, 301, et, 10013, "OBSERVER", "CN+S", p, 399, 10013)) isa T
        @test (@inferred spkcpt(ks, p, 399, 10013, et, 1, "TARGET", "LT+S", 301)) isa T
        @test (@inferred spkcvo(ks, 301, et, 10020, "CENTER", "XCN", s, 0.0, 399, 10013)) isa T
        @test (@inferred spkcvt(ks, s, 0.0, 399, 1, et, 10013, "OBSERVER", "LT", 301)) isa T
        @test (@inferred spkcpo(ctx, "MOON", et, "IAU_EARTH", "OBSERVER", "CN+S", p, "EARTH", "IAU_EARTH")) isa T
    end

    @testset "epoch as a dual number: derivative against the velocity and the central difference" begin
        st = GeoStats()
        h = 20.0
        p = SVector(6378.0 * cos(0.5), 0.0, 6378.0 * sin(0.5))
        s = SVector(p..., 0.0, 0.01, 0.0)
        for et in record_safe_epochs(rng, 24), cor in GEO_CORRS, loc in CO_LOCI
            fs = Pair{String,Function}[
                "spkcpo" => t -> spkcpo(ks, 301, t, 10013, loc, cor, p, 399, 10013)[1],
                "spkcpt" => t -> spkcpt(ks, p, 399, 10013, t, 1, loc, cor, 301)[1],
                "spkcvo" => t -> spkcvo(ks, 301, t, 10020, loc, cor, s, et - 300.0, 399, 10013)[1],
                "spkcvt" => t -> spkcvt(ks, s, et - 300.0, 399, 1, t, 10013, loc, cor, 301)[1]]
            for (name, f) in fs
                d = f(ForwardDiff.Dual(et, 1.0))
                v = f(et)
                @test isapprox(ForwardDiff.value.(d), v; rtol=1e-12)
                dp = ForwardDiff.partials.(d, 1)[SOneTo(3)]
                fd = geo_fd(t -> f(t)[SOneTo(3)], et, h)
                gput!(st, "$name, d(position)/dt against central difference", norm(dp .- fd) / norm(fd))
                # The velocity of CSPICE is an approximation of this derivative (it neglects
                # terms of the light-time derivative); the difference is reported.
                gput!(st, "$name, d(position)/dt against the returned velocity", norm(dp .- v[SVector(4, 5, 6)]) / norm(v[SVector(4, 5, 6)]))
            end
        end
        greport("dual epoch", st)
        for name in ("spkcpo", "spkcpt", "spkcvo", "spkcvt")
            @test gget(st, "$name, d(position)/dt against central difference") < 1e-4
            @test gget(st, "$name, d(position)/dt against the returned velocity") < 1e-3
        end
    end
end
