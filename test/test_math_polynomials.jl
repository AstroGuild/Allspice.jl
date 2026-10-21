# Group 20: Chebyshev expansions, Lagrange and Hermite interpolation, derivatives of a
# polynomial, the numerical derivative, the roots of a quadratic and the
# diagonalization of a symmetric 2 by 2 matrix, against CSPICE.
#
# The random inputs are expansions of degree 0 to 20 and point sets of 1 to 12 points,
# with values from 1e-2 to 1e2 and abscissas inside and outside the interval of the
# data; quadratics with a zero leading coefficient, a zero or a negative discriminant
# and coefficients of very different magnitudes; and symmetric matrices with a zero
# off-diagonal element or two equal diagonal elements.
#
# Tolerance: 0 for the calls whose arithmetic is the one of CSPICE and has no library
# function but the square root, which is correctly rounded in both. The measured
# difference of every call of this group is 0.

FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff
isdefined(@__MODULE__, :CS) || include("mathoracle.jl")

group("20. math: polynomials and interpolation") do
    rng = MersenneTwister(SEED + 20)
    n = nsample(150, 1500)
    t = Tally()
    rval() = 10.0^(4 * rand(rng) - 2) * randn(rng)
    "Abscissas: sorted, unsorted, or with a repeated value."
    function rxs(k)
        x = 100 * randn(rng) .+ cumsum(10.0^(2 * rand(rng) - 1) .* (0.1 .+ rand(rng, k)))
        u = rand(rng)
        u < 0.3 && shuffle!(rng, x)
        (u > 0.95 && k >= 2) && (x[rand(rng, 2:k)] = x[1])
        return x
    end
    "An abscissa inside the data, or outside it by up to half its span."
    function rat(xs)
        lo, hi = extrema(xs)
        w = max(hi - lo, 1.0)
        return lo - 0.5w + 2w * rand(rng)
    end

    @apart @testset "Chebyshev expansions" begin
        for _ in 1:n
            cp = [rval() for _ in 1:rand(rng, 1:21)]
            rad = 10.0^(6 * rand(rng) - 3)
            x2s = [1e3 * randn(rng), rad]
            x = x2s[1] + rad * (2.4 * rand(rng) - 1.2)
            versus!(t, :chbval, () -> chbval(cp, x2s, x), () -> CS.chbval(cp, x2s, x))
            versus!(t, :chbint, () -> chbint(cp, x2s, x), () -> CS.chbint(cp, x2s, x))
            k = rand(rng, 0:6)
            versus!(t, :chbder, () -> chbder(cp, x2s, x, k), () -> CS.chbder(cp, x2s, x, k))
            # the checks of chbigr: an empty expansion, a radius that is not positive
            u = rand(rng)
            c = u < 0.05 ? Float64[] : cp
            r = u > 0.95 ? [x2s[1], -rand(rng, (0.0, rad))] : x2s
            versus!(t, :chbigr, () -> chbigr(c, r, x), () -> CS.chbigr(c, r, x))
        end
        # chbint and chbder give the derivative of chbval; chbigr its integral
        cp, x2s = [0.3, -1.2, 0.8, 0.25, -0.1], [2.0, 3.0]
        @test chbint(cp, x2s, 2.7)[2] ≈ ForwardDiff.derivative(x -> chbval(cp, x2s, x), 2.7) rtol = 1e-14
        @test chbder(cp, x2s, 2.7, 2)[3] ≈ ForwardDiff.derivative(x -> chbint(cp, x2s, x)[2], 2.7) rtol = 1e-14
        h = 1e-5
        @test (chbigr(cp, x2s, 3.1 + h)[2] - chbigr(cp, x2s, 3.1 - h)[2]) / 2h ≈ chbval(cp, x2s, 3.1) rtol = 1e-9
        @test chbigr(cp, x2s, 2.0)[2] == 0.0
        @test throws_kind(() -> chbval(Float64[], x2s, 1.0), :invaliddegree)
        @test throws_kind(() -> chbval(cp, [1.0], 1.0), :badarg)
        @test throws_kind(() -> chbder(cp, x2s, 1.0, -1), :badarg)
    end

    @apart @testset "Lagrange and Hermite interpolation" begin
        for _ in 1:n
            k = rand(rng, 1:12)
            xs = rxs(k)
            ys = [rval() for _ in 1:k]
            yd = [rval() for _ in 1:2k]
            x = rat(xs)
            versus!(t, :lgrint, () -> lgrint(xs, ys, x), () -> CS.lgrint(xs, ys, x))
            versus!(t, :lgrind, () -> lgrind(xs, ys, x), () -> CS.lgrind(xs, ys, x))
            versus!(t, :hrmint, () -> hrmint(xs, yd, x), () -> CS.hrmint(xs, yd, x))
            first, step = 100 * randn(rng), rand(rng) < 0.05 ? 0.0 : 10.0^(2 * rand(rng) - 1) * randn(rng)
            x = first + step * (k + 1) * (1.5 * rand(rng) - 0.25)
            versus!(t, :lgresp, () -> lgresp(first, step, ys, x), () -> CS.lgresp(first, step, ys, x))
            versus!(t, :hrmesp, () -> hrmesp(first, step, yd, x), () -> CS.hrmesp(first, step, yd, x))
        end
        # Fixed cases of the checks: no point, and equal abscissas that are neighbours
        # (found in the first column of the Hermite table) or not (found later).
        for (xs, x) in ((Float64[], 1.0), ([1.0, 1.0], 0.5), ([1.0, 2.0, 1.0], 0.5))
            k = length(xs)
            versus!(t, :lgrint, () -> lgrint(xs, ones(k), x), () -> CS.lgrint(xs, ones(k), x))
            versus!(t, :lgrind, () -> lgrind(xs, ones(k), x), () -> CS.lgrind(xs, ones(k), x))
            versus!(t, :hrmint, () -> hrmint(xs, ones(2k), x), () -> CS.hrmint(xs, ones(2k), x))
        end
        versus!(t, :lgresp, () -> lgresp(0.0, 1.0, Float64[], 1.0), () -> CS.lgresp(0.0, 1.0, Float64[], 1.0))
        versus!(t, :hrmesp, () -> hrmesp(0.0, 1.0, Float64[], 1.0), () -> CS.hrmesp(0.0, 1.0, Float64[], 1.0))
        @test throws_kind(() -> lgrint([1.0, 2.0], [1.0], 1.0), :badarg)
        @test throws_kind(() -> hrmint([1.0, 2.0], [1.0, 2.0], 1.0), :badarg)
        @test throws_kind(() -> hrmesp(0.0, 1.0, [1.0, 2.0, 3.0], 1.0), :badarg)
        # the polynomial of degree n - 1 (2n - 1) is reproduced, and lgrind and hrmint
        # give its derivative
        p(x) = 2.0 - x + 0.5x^2 - 0.25x^3
        dp(x) = -1.0 + x - 0.75x^2
        xs = [-1.0, 0.5, 2.0, 3.0]
        @test lgrint(xs, p.(xs), 1.3) ≈ p(1.3) rtol = 1e-14
        @test all(lgrind(xs, p.(xs), 1.3) .≈ (p(1.3), dp(1.3)))
        @test all(hrmint(xs[1:2], [p(-1.0), dp(-1.0), p(0.5), dp(0.5)], 1.3) .≈ (p(1.3), dp(1.3)))
        @test lgresp(-1.0, 1.5, p.([-1.0, 0.5, 2.0, 3.5]), 1.3) ≈ p(1.3) rtol = 1e-14
    end

    @apart @testset "derivatives of a polynomial and numerical derivative" begin
        for _ in 1:n
            c = [rval() for _ in 1:rand(rng, 0:10)]
            k = rand(rng, -1:6)
            x = 3 * randn(rng)
            versus!(t, :polyds, () -> polyds(c, k, x), () -> CS.polyds(c, k, x))
            d = rand(rng, 1:6)
            f0, f2 = [rval() for _ in 1:d], [rval() for _ in 1:d]
            delta = rand(rng) < 0.05 ? 0.0 : 10.0^(4 * rand(rng) - 2)
            versus!(t, :qderiv, () -> qderiv(f0, f2, delta), () -> CS.qderiv(f0, f2, delta))
        end
        @test polyds(Float64[], 2, 1.5) == [0.0, 0.0, 0.0]
        @test isempty(polyds([1.0, 2.0], -1, 1.5))
        @test throws_kind(() -> qderiv([1.0], [1.0, 2.0], 1.0), :badarg)
    end

    @apart @testset "quadratic equations and symmetric 2 by 2 matrices" begin
        resid = 0.0         # largest |rotate' symmat rotate - diag| / |symmat|
        rotdev = 0.0        # largest |det(rotate) - 1|
        for _ in 1:n
            u = rand(rng)
            a = u < 0.1 ? 0.0 : rval()
            b = u > 0.9 ? 0.0 : rval()
            c = rval()
            # a zero or a negative discriminant, and coefficients of very different size
            v = rand(rng)
            v < 0.1 && (c = b * b / (4a))
            v > 0.9 && ((a, c) = (1e-300 * a, 1e300 * c))
            (a == 0.0 && b == 0.0) && (b = 1.0)     # see the note on the degenerate case
            versus!(t, :rquad, () -> rquad(a, b, c), () -> CS.rquad(a, b, c))
            s = [rval() 0.0; 0.0 rval()]
            w = rand(rng)
            s[1, 2] = s[2, 1] = w < 0.15 ? 0.0 : rval()
            w > 0.85 && (s[2, 2] = s[1, 1])
            versus!(t, :diags2, () -> diags2(s), () -> CS.diags2(s))
            d, r = diags2(s)
            resid = max(resid, maximum(abs, r' * s * r - d) / maximum(abs, s))
            rotdev = max(rotdev, abs(r[1, 1] * r[2, 2] - r[1, 2] * r[2, 1] - 1))
        end
        measured(@sprintf("  diags2: residual %.1e, |det(rotate) - 1| %.1e", resid, rotdev))
        # The residual is rounding of the products: up to 4.9e-14 measured in the
        # exhaustive sample, where the two diagonal values differ by many orders.
        @test resid < 1e-13 && rotdev < 1e-15
        # `rquad_c` signals SPICE(DEGENERATECASE) for a == b == 0 (measured); its
        # traceback then has a wrong routine name, so the case is not sent to CSPICE.
        @test throws_kind(() -> rquad(0.0, 0.0, 1.0), :degeneratecase)
        @test diags2([3.0 0.0; 7.0 -1.0]) == ([3.0 0.0; 0.0 -1.0], [1.0 0.0; 0.0 1.0])
        @test throws_kind(() -> diags2(ones(3, 3)), :badarg)
    end

    @apart @testset "differences against CSPICE" begin
        report!(t, "polynomials"; tol=0.0)
    end

    @apart @testset "number types: ForwardDiff and Taylor" begin
        xs, ys = [0.0, 1.0, 2.5, 4.0], [1.0, -2.0, 0.5, 3.0]
        @test ForwardDiff.derivative(x -> lgrint(xs, ys, x), 1.7) ≈ lgrind(xs, ys, 1.7)[2] rtol = 1e-14
        # the derivative with respect to the data: a Lagrange basis polynomial
        g = ForwardDiff.gradient(y -> lgrint(xs, y, 1.7), ys)
        @test g ≈ [lgrint(xs, Float64.(1:4 .== i), 1.7) for i in 1:4] rtol = 1e-14
        yd = [1.0, 0.5, -2.0, 1.0, 0.5, 0.0, 3.0, -1.0]
        @test ForwardDiff.derivative(x -> hrmint(xs, yd, x)[1], 1.7) ≈ hrmint(xs, yd, 1.7)[2] rtol = 1e-13
        tt = Allspice.taylor_var(1.7, Val(2))
        p = lgrind(xs, ys, tt)[1]
        @test Allspice.derivatives(p)[2] ≈ lgrind(xs, ys, 1.7)[2] rtol = 1e-14
        cp, x2s = [0.3, -1.2, 0.8, 0.25], [2.0, 3.0]
        q = chbval(cp, x2s, tt)
        @test collect(Allspice.derivatives(q)) ≈ chbder(cp, x2s, 1.7, 2) rtol = 1e-14
        @test ForwardDiff.derivative(x -> rquad(1.0, -x, 2.0)[1][1], 3.0) ≈ (1 + 3 / sqrt(1.0)) / 2 rtol = 1e-14
    end

    @apart @testset "inferred types and no allocation for static arguments" begin
        xs, ys = SVector(0.0, 1.0, 2.5, 4.0), SVector(1.0, -2.0, 0.5, 3.0)
        yd = SVector(1.0, 0.5, -2.0, 1.0, 0.5, 0.0, 3.0, -1.0)
        cp, x2s = SVector(0.3, -1.2, 0.8, 0.25), SVector(2.0, 3.0)
        s = SMatrix{2,2}(2.0, 1.0, 1.0, 3.0)
        M2 = SMatrix{2,2,Float64,4}
        calls = Any[(lgrint, (xs, ys, 1.7), Float64), (lgrind, (xs, ys, 1.7), NTuple{2,Float64}),
            (lgresp, (0.0, 1.5, ys, 1.7), Float64), (hrmint, (xs, yd, 1.7), NTuple{2,Float64}),
            (hrmesp, (0.0, 1.5, yd, 1.7), NTuple{2,Float64}), (chbval, (cp, x2s, 1.7), Float64),
            (chbint, (cp, x2s, 1.7), NTuple{2,Float64}), (chbigr, (cp, x2s, 1.7), NTuple{2,Float64}),
            (qderiv, (xs, ys, 0.5), SVector{4,Float64}),
            (rquad, (1.0, 2.0, 3.0), NTuple{2,SVector{2,Float64}}), (diags2, (s,), Tuple{M2,M2})]
        for (f, args, T) in calls
            @test (@inferred f(args...)) isa T
            @test alloc_of(f, args...) == 0 skip = STACK_BUFFERS_ESCAPE && f in (lgrint, lgrind)
        end
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
