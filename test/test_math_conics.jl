# Group 16: two-body calls (Stumpff functions, propagation, conic elements) against
# CSPICE.
#
# Feature file: features/math_conics.feature
#
# The random orbits are ellipses, parabolas and hyperbolas about masses from 1e-3 to
# 1e11 km^3/s^2, with perifocal distances from 1e-2 to 1e8 km. The degenerate inputs
# are the circular and the equatorial orbits (prograde and retrograde), eccentricities
# next to 1, a zero or a parallel position and velocity, and a mass that is not positive.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff
isdefined(@__MODULE__, :CS) || include("mathoracle.jl")

group("16. math: conics") do
    rng = MersenneTwister(SEED + 16)
    n = nsample(100, 400)
    t = Tally()
    r3() = SVector(randn(rng), randn(rng), randn(rng))
    rmu() = rand(rng, (398600.4418, 1.32712440018e11, 4902.800066, 1.0, 1e-3, 3.986e5))
    "Conic elements; with some probability on a boundary of an element range."
    function relts()
        mu = rmu()
        rp = 10.0^(10 * rand(rng) - 2)
        u = rand(rng)
        ecc = u < 0.1 ? 0.0 : u < 0.55 ? rand(rng) : u < 0.65 ? 1.0 : u < 0.7 ? 1.0 - 1e-9 * rand(rng) :
              u < 0.75 ? 1.0 + 1e-9 * rand(rng) : u < 0.95 ? 1.0 + 3 * rand(rng) : 10.0^(2 * rand(rng))
        v = rand(rng)
        inc = v < 0.1 ? 0.0 : v < 0.15 ? Float64(pi) : v < 0.2 ? pi / 2 : pi * rand(rng)
        lnode, argp = 2pi * rand(rng), 2pi * rand(rng)
        rand(rng) < 0.1 && (argp = 0.0)
        rand(rng) < 0.1 && (lnode = 0.0)
        m0 = rand(rng) < 0.1 ? 0.0 : ecc < 1.0 ? 2pi * rand(rng) : 4 * randn(rng)
        t0 = 1e8 * randn(rng)
        w = rand(rng)
        w < 0.02 && (ecc = -ecc - 0.1)
        0.02 <= w < 0.04 && (rp = -rp)
        0.04 <= w < 0.06 && (mu = rand(rng, (0.0, -mu)))
        return SVector(rp, ecc, inc, lnode, argp, m0, t0, mu)
    end
    "Time scale of an orbit: the period of the circular orbit at the perifocal distance."
    scale_of(el) = 2pi * sqrt(abs(el[1])^3 / max(abs(el[8]), 1e-300))
    "State about `mu`: speed from 0.1 to 3 times the circular one; some are degenerate."
    function rst(mu)
        r = 10.0^(10 * rand(rng) - 2)
        pos = r * normalize(r3())
        vc = sqrt(abs(mu) / r)
        u = rand(rng)
        vel = u < 0.1 ? vc * normalize(vcrss(pos, r3())) :                 # circular
              u < 0.15 ? sqrt(2.0) * vc * normalize(vcrss(pos, r3())) :     # parabolic at periapsis
              vc * (0.1 + 2.9 * rand(rng)) * normalize(r3())
        w = rand(rng)
        w < 0.08 && (pos = SVector(pos[1], pos[2], 0.0); vel = SVector(vel[1], vel[2], 0.0))
        0.08 <= w < 0.11 && (pos = SVector(pos[1], 0.0, 0.0); vel = SVector(0.0, -abs(vel[2]) - vc, 0.0))
        0.11 <= w < 0.13 && (vel = 0.3 * vc * normalize(pos))
        0.13 <= w < 0.15 && (vel = SVector(0.0, 0.0, 0.0))
        0.15 <= w < 0.17 && (pos = SVector(0.0, 0.0, 0.0))
        return vcat(pos, vel)
    end

    @apart @testset "Stumpff functions" begin
        for _ in 1:n
            u = rand(rng)
            x = u < 0.3 ? 2.2 * rand(rng) - 1.1 : u < 0.4 ? rand(rng, (0.0, 1.0, -1.0, 1.0 + eps(), -1.0 - eps())) :
                u < 0.9 ? rmag(rng)^0.5 * randn(rng) : -4.0e5 - 2.0e5 * rand(rng)
            versus!(t, :stmp03, () -> stmp03(x), () -> CS.stmp03_(x)[2:end])
        end
        for x in (-30.0, -1.0, -0.3, 0.0, 0.4, 1.0, 55.0)
            c0, c1, c2, c3 = stmp03(x)
            @test c0 ≈ 1 - x * c2 rtol = 1e-14
            @test c1 ≈ 1 - x * c3 rtol = 1e-14
        end
        @test stmp03(0.0) == (1.0, 1.0, 0.5, 1 / 6)
        @test throws_kind(() -> stmp03(-1.0e6), :valueoutofrange)
        d = ForwardDiff.derivative(x -> SVector(stmp03(x)), 0.3)
        h = 1e-6
        @test d ≈ (SVector(stmp03(0.3 + h)) - SVector(stmp03(0.3 - h))) / 2h rtol = 1e-8
    end

    @apart @testset "propagation of a state" begin
        for _ in 1:2n
            mu = rmu()
            rand(rng) < 0.03 && (mu = rand(rng, (0.0, -1.0)))
            s = rst(mu)
            period = 2pi * sqrt(max(vnorm(s[SVector(1, 2, 3)]), 1e-3)^3 / max(abs(mu), 1e-300))
            u = rand(rng)
            dt = u < 0.05 ? 0.0 : u < 0.6 ? period * 2 * randn(rng) : u < 0.9 ? period * 1e-3 * randn(rng) :
                 u < 0.97 ? period * 300 * randn(rng) : rmag(rng)^2 * randn(rng)
            d = versus!(t, :prop2b, () -> prop2b(mu, s, dt), () -> CS.prop2b(mu, s, dt))
            if isfinite(d) && !(outcome(() -> prop2b(mu, s, dt)) isa Symbol)
                # energy and angular momentum do not change
                q = prop2b(mu, s, dt)
                p0, v0 = s[SVector(1, 2, 3)], s[SVector(4, 5, 6)]
                p1, v1 = q[SVector(1, 2, 3)], q[SVector(4, 5, 6)]
                h0 = vcrss(p0, v0)
                if abs(dt) < 10 * period && vnorm(h0) > 1e-3 * vnorm(p0) * vnorm(v0)
                    @test vdist(vcrss(p1, v1), h0) <= 1e-8 * vnorm(h0)
                    e0 = (v0 ⋅ v0) / 2 - mu / vnorm(p0)
                    @test abs((v1 ⋅ v1) / 2 - mu / vnorm(p1) - e0) <= 1e-8 * ((v0 ⋅ v0) / 2 + mu / vnorm(p0))
                end
            end
        end
        # Fixed cases of the bracket search: for this circular orbit the first value of the
        # Kepler function is on the wrong side of the time interval by rounding, so the
        # bound of the bracket is doubled (forward and backward in time). A random state
        # reaches these loops about once in a hundred.
        sc = SVector(1.0e5, 0.0, 0.0, 0.0, sqrt(1.0e-5), 0.0)
        for dt in (7.0, -7.0)
            versus!(t, :prop2b, () -> prop2b(1.0, sc, dt), () -> CS.prop2b(1.0, sc, dt))
        end
        mu = 398600.4418
        s = SVector(7000.0, 100.0, -300.0, 0.3, 7.2, 1.9)
        @test prop2b(mu, s, 0.0) ≈ s rtol = 4eps()
        @test prop2b(mu, prop2b(mu, s, 5000.0), -5000.0) ≈ s rtol = 1e-12
        @test prop2b(mu, prop2b(mu, s, 1234.0), 4321.0) ≈ prop2b(mu, s, 5555.0) rtol = 1e-12
        @test throws_kind(() -> prop2b(0.0, s, 1.0), :nonpositivemass)
        @test throws_kind(() -> prop2b(mu, [0.0, 0, 0, 1, 2, 3], 1.0), :zeroposition)
        @test throws_kind(() -> prop2b(mu, [1.0, 2, 3, 0, 0, 0], 1.0), :zerovelocity)
        @test throws_kind(() -> prop2b(mu, [1.0, 2, 3, 2, 4, 6], 1.0), :nonconicmotion)
        # the derivative with respect to the time is the velocity and the acceleration
        for dt in (-9000.0, 0.0, 10.0, 3.0e4)
            q = prop2b(mu, s, dt)
            p1, v1 = q[SVector(1, 2, 3)], q[SVector(4, 5, 6)]
            d = ForwardDiff.derivative(τ -> prop2b(mu, s, τ), dt)
            @test d[SVector(1, 2, 3)] ≈ v1 rtol = 1e-12
            @test d[SVector(4, 5, 6)] ≈ -mu * p1 / vnorm(p1)^3 rtol = 1e-12
            tt = Allspice.taylor_var(dt, Val(2))
            x = prop2b(mu, s, tt)
            @test map(c -> Allspice.derivatives(c)[1], x) == q
            @test map(c -> Allspice.derivatives(c)[2], x) ≈ d rtol = 1e-12
            @test map(c -> Allspice.derivatives(c)[3], x)[SVector(1, 2, 3)] ≈ -mu * p1 / vnorm(p1)^3 rtol = 1e-11
        end
    end

    @apart @testset "state from elements" begin
        for _ in 1:2n
            el = relts()
            u = rand(rng)
            et = el[7] + (u < 0.1 ? 0.0 : u < 0.7 ? 3 * randn(rng) : 1e-3 * randn(rng)) * scale_of(el)
            # The sine and the cosine of Julia and of the C library differ in the last
            # place for some angles. The energy of an orbit next to a parabola is the
            # small difference of two terms, so that rounding of the periapsis state
            # changes the period by a relative 1e-16 / |1 - ecc|.
            name = abs(el[2] - 1.0) < 1e-6 && el[2] != 1.0 ? :conics_nearparabola : :conics
            versus!(t, name, () -> conics(el, et), () -> CS.conics(el, et))
        end
        el = SVector(7000.0, 0.1, 0.5, 1.0, 2.0, 0.0, 100.0, 398600.4418)
        s = conics(el, 100.0)
        @test vnorm(s[SVector(1, 2, 3)]) ≈ 7000.0 rtol = 1e-14                    # at periapsis
        @test abs(s[SVector(1, 2, 3)] ⋅ s[SVector(4, 5, 6)]) <= 1e-9
        period = 2pi * sqrt((7000.0 / 0.9)^3 / 398600.4418)
        @test conics(el, 100.0 + period) ≈ s rtol = 1e-11
        @test throws_kind(() -> conics(setindex(el, -0.1, 2), 0.0), :badeccentricity)
        @test throws_kind(() -> conics(setindex(el, 0.0, 1), 0.0), :badperiapsevalue)
        @test throws_kind(() -> conics(setindex(el, 0.0, 8), 0.0), :badgm)
        @test throws_kind(() -> conics([1.0, 2.0], 0.0), :badarg)
        d = ForwardDiff.derivative(τ -> conics(el, τ), 2000.0)
        q = conics(el, 2000.0)
        @test d[SVector(1, 2, 3)] ≈ q[SVector(4, 5, 6)] rtol = 1e-12
        @test d[SVector(4, 5, 6)] ≈ -el[8] * q[SVector(1, 2, 3)] / vnorm(q[SVector(1, 2, 3)])^3 rtol = 1e-12
    end

    @apart @testset "elements from a state" begin
        for _ in 1:2n
            mu = rmu()
            rand(rng) < 0.03 && (mu = rand(rng, (0.0, -1.0)))
            s = rst(mu)
            et = 1e8 * randn(rng)
            versus!(t, :oscelt, () -> oscelt(s, et, mu), () -> CS.oscelt(s, et, mu))
            versus!(t, :oscltx, () -> oscltx(s, et, mu), () -> CS.oscltx(s, et, mu))
            el = outcome(() -> oscelt(s, et, mu))
            if !(el isa Symbol)
                @test el[1] > 0.0 && el[2] >= 0.0 && 0.0 <= el[3] <= pi
                @test 0.0 <= el[4] < 2pi + 1e-12 && 0.0 <= el[5] < 2pi + 1e-12
                @test el[7] == et && el[8] == mu
                # the round trip; not for an orbit that is almost rectilinear or parabolic,
                # where the elements do not fix the state to full precision
                p0, v0 = s[SVector(1, 2, 3)], s[SVector(4, 5, 6)]
                if vnorm(vcrss(p0, v0)) > 0.05 * vnorm(p0) * vnorm(v0) && abs(el[2] - 1.0) > 1e-3
                    back = conics(el, et)
                    @test vdist(back[SVector(1, 2, 3)], p0) <= 1e-7 * vnorm(p0)
                    @test vdist(back[SVector(4, 5, 6)], v0) <= 1e-7 * vnorm(v0)
                end
            end
        end
        mu = 398600.4418
        s = SVector(7000.0, 100.0, -300.0, 0.3, 7.2, 1.9)
        el = oscelt(s, 50.0, mu)
        @test conics(el, 50.0) ≈ s rtol = 1e-13
        x = oscltx(s, 50.0, mu)
        @test x[SVector{8}(1:8)] == el
        a = 1 / (2 / vnorm(s[SVector(1, 2, 3)]) - (s[SVector(4, 5, 6)] ⋅ s[SVector(4, 5, 6)]) / mu)
        @test x[10] ≈ a rtol = 1e-13
        @test x[11] ≈ 2pi * sqrt(a^3 / mu) rtol = 1e-13
        @test x[1] ≈ a * (1 - x[2]) rtol = 1e-13
        @test all(iszero, x[SVector{9}(12:20)])
        # a circular equatorial orbit: the node is at the x axis, and the argument of
        # periapsis plus the mean anomaly is the longitude of the body
        c = oscelt([7000.0, 0, 0, 0, sqrt(mu / 7000.0), 0], 0.0, mu)
        @test c[2] < 1e-15 && c[3] == 0.0 && c[4] == 0.0 && abs(rem2pi(c[5] + c[6], RoundNearest)) < 1e-15
        @test c == CS.oscelt([7000.0, 0, 0, 0, sqrt(mu / 7000.0), 0], 0.0, mu)
        @test throws_kind(() -> oscelt(s, 0.0, 0.0), :nonpositivemass)
        @test throws_kind(() -> oscelt([1.0, 2, 3, 2, 4, 6], 0.0, mu), :degeneratecase)
        @test throws_kind(() -> oscelt([0.0, 0, 0, 2, 4, 6], 0.0, mu), :degeneratecase)
        # the derivative of the elements with respect to the state
        jac = ForwardDiff.jacobian(v -> oscelt(v, 50.0, mu), s)
        for k in 1:6
            h = 1e-4 * (k <= 3 ? 1.0 : 1e-3)
            e = SVector{6}(ntuple(i -> i == k ? h : 0.0, 6))
            fd = (oscelt(s + e, 50.0, mu) - oscelt(s - e, 50.0, mu)) / 2h
            @test jac[:, k] ≈ fd rtol = 1e-6 atol = 1e-9
        end
    end

    @apart @testset "differences against CSPICE" begin
        report!(t, "conics"; tol=4eps(), tols=Dict(:prop2b => 1e-13, :conics => 1e-10,
            :conics_nearparabola => 1e-1))
    end

    @apart @testset "inferred types and no allocation for static arguments" begin
        mu = 398600.4418
        s = SVector(7000.0, 100.0, -300.0, 0.3, 7.2, 1.9)
        el = oscelt(s, 50.0, mu)
        calls = Any[(stmp03, (0.3,), NTuple{4,Float64}), (stmp03, (-7.0,), NTuple{4,Float64}),
            (prop2b, (mu, s, 5000.0), SVector{6,Float64}), (conics, (el, 900.0), SVector{6,Float64}),
            (oscelt, (s, 50.0, mu), SVector{8,Float64}), (oscltx, (s, 50.0, mu), SVector{20,Float64})]
        for (fn, args, T) in calls
            @test (@inferred fn(args...)) isa T
            @test alloc_of(fn, args...) == 0
        end
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
