# Group 12: vector and matrix calls against CSPICE.
#
# Feature file: features/math_vectors.feature
#
# Each function runs on several hundred seeded random inputs with magnitudes from 1e-10
# to 1e10, mixed with the degenerate inputs (zero vector, unit axes, equal, parallel and
# opposite pairs, very small and very large components). The limit of each function is
# 0 where NeoSpice has the operation order of the `_c` function, and a few units of
# rounding where CSPICE runs the f2c routine with another order.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of NeoSpice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`NEOSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0

import ForwardDiff
isdefined(@__MODULE__, :CS) || include("mathoracle.jl")

group("12. math: vectors and matrices") do
    rng = MersenneTwister(SEED + 12)
    n = 400
    t = Tally()
    rs() = rmag(rng) * randn(rng)
    rm3() = rand(rng) < 0.2 ? Matrix(rrot(rng)) : rmag(rng) * randn(rng, 3, 3)
    rax() = rand(rng, -4:7)

    @testset "three-dimensional vectors" begin
        for _ in 1:n
            a, b = rpair(rng)
            c = rvec(rng)
            x, y, z = rs(), rs(), rs()
            th = 20 * rand(rng) - 10
            versus!(t, :vadd, () -> vadd(a, b), () -> CS.vadd(a, b))
            versus!(t, :vsub, () -> vsub(a, b), () -> CS.vsub(a, b))
            versus!(t, :vscl, () -> vscl(x, a), () -> CS.vscl(x, a))
            versus!(t, :vminus, () -> vminus(a), () -> CS.vminus(a))
            versus!(t, :vdot, () -> vdot(a, b), () -> CS.vdot(a, b))
            versus!(t, :vcrss, () -> vcrss(a, b), () -> CS.vcrss(a, b))
            versus!(t, :vnorm, () -> vnorm(a), () -> CS.vnorm(a))
            versus!(t, :vhat, () -> vhat(a), () -> CS.vhat(a))
            versus!(t, :unorm, () -> unorm(a), () -> CS.unorm(a))
            versus!(t, :vdist, () -> vdist(a, b), () -> CS.vdist(a, b))
            versus!(t, :vsep, () -> vsep(a, b), () -> CS.vsep(a, b))
            versus!(t, :vrel, () -> vrel(a, b), () -> CS.vrel(a, b))
            versus!(t, :vproj, () -> vproj(a, b), () -> CS.vproj(a, b))
            # `vperp_c` gives zero for a zero `b` (see the docstring of `vperp`); the
            # f2c routine is the oracle there.
            versus!(t, :vperp, () -> vperp(a, b), () -> iszero(b) ? CS.vperp_(a, b)[2] : CS.vperp(a, b))
            versus!(t, :vlcom, () -> vlcom(x, a, y, b), () -> CS.vlcom(x, a, y, b))
            versus!(t, :vlcom3, () -> vlcom3(x, a, y, b, z, c), () -> CS.vlcom3(x, a, y, b, z, c))
            versus!(t, :vrotv, () -> vrotv(a, b, th), () -> CS.vrotv(a, b, th))
            versus!(t, :vzero, () -> vzero(a), () -> CS.vzero(a) != 0)
            versus!(t, :vequ, () -> vequ(a), () -> CS.vequ(a))
            versus!(t, :vpack, () -> vpack(x, y, z), () -> CS.vpack(x, y, z))
            versus!(t, :vupack, () -> vupack(a), () -> CS.vupack(a))
            versus!(t, :ucrss, () -> ucrss(a, b), () -> CS.ucrss(a, b))
            m = rm3()
            versus!(t, :vtmv, () -> vtmv(a, m, b), () -> CS.vtmv(a, m, b))
        end
        # The stated behaviour at degenerate inputs.
        z3 = SVector(0.0, 0.0, 0.0)
        e1 = SVector(1.0, 0.0, 0.0)
        @test vnorm(z3) == 0.0 && vhat(z3) == z3 && unorm(z3) == (z3, 0.0)
        @test vsep(z3, e1) == 0.0 && vsep(e1, z3) == 0.0
        @test vsep(e1, e1) == 0.0 && vsep(e1, -e1) == Float64(pi) && vsep(e1, SVector(0.0, 2.0, 0.0)) == pi / 2
        @test vproj(e1, z3) == z3 && vproj(z3, e1) == z3
        @test vperp(e1, z3) == e1 && vperp(z3, e1) == z3
        @test CS.vperp(e1, z3) == z3            # the defect of `vperp_c`
        @test ucrss(e1, 5.0 * e1) == z3 && ucrss(z3, e1) == z3
        @test vrotv(e1, z3, 1.0) == e1
        @test vrel(z3, z3) == 0.0 && vrel(e1, e1) == 0.0
        @test vnorm(SVector(1e300, 1e300, 1e300)) == CS.vnorm([1e300, 1e300, 1e300]) < Inf
        @test vnorm(SVector(1e-300, 1e-300, 0.0)) > 0.0
        @test throws_kind(() -> vadd([1.0, 2.0], [1.0, 2.0, 3.0]), :badarg)
        @test vadd([1, 2, 3], [1, 1, 1]) === SVector(2.0, 3.0, 4.0)
    end

    @testset "three-dimensional matrices" begin
        for _ in 1:n
            m1, m2 = rm3(), rm3()
            a = rvec(rng)
            th = 20 * rand(rng) - 10
            versus!(t, :mxv, () -> mxv(m1, a), () -> CS.mxv(m1, a))
            versus!(t, :mtxv, () -> mtxv(m1, a), () -> CS.mtxv(m1, a))
            versus!(t, :mxm, () -> mxm(m1, m2), () -> CS.mxm(m1, m2))
            versus!(t, :mtxm, () -> mtxm(m1, m2), () -> CS.mtxm(m1, m2))
            versus!(t, :mxmt, () -> mxmt(m1, m2), () -> CS.mxmt(m1, m2))
            versus!(t, :xpose, () -> xpose(m1), () -> CS.xpose(m1))
            m6 = rmag(rng) * randn(rng, 6, 6)
            versus!(t, :xpose6, () -> xpose6(m6), () -> CS.xpose6(m6))
            versus!(t, :det, () -> NeoSpice.det(m1), () -> CS.det(m1))
            versus!(t, :trace, () -> trace(m1), () -> CS.trace(m1))
            # singular matrices: a repeated row, and a tiny determinant
            ms = rand(rng) < 0.2 ? vcat(m1[1:2, :], m1[1:1, :]) : rand(rng) < 0.1 ? 1e-6 * m1 / maximum(abs, m1) : m1
            versus!(t, :invert, () -> invert(ms), () -> CS.invert(ms))
            # orthogonal columns of any length, sometimes a zero column
            q = Matrix(rrot(rng)) * Diagonal([rmag(rng), rmag(rng), rand(rng) < 0.1 ? 0.0 : rmag(rng)])
            versus!(t, :invort, () -> invort(q), () -> CS.invort(q))
            # near rotations and tolerances on both sides of the deviation
            r = Matrix(rrot(rng)) + (rand(rng) < 0.5 ? 0.0 : 10.0^(-rand(rng, 1:12))) * randn(rng, 3, 3)
            ntol = rand(rng) < 0.05 ? -1.0 : 10.0^(-rand(rng, 0:14))
            dtol = rand(rng) < 0.05 ? -1.0 : 10.0^(-rand(rng, 0:14))
            versus!(t, :isrot, () -> isrot(r, ntol, dtol), () -> CS.isrot(r, ntol, dtol) != 0)
            ax = rax()
            versus!(t, :rotate, () -> rotate(th, ax), () -> CS.rotate(th, ax))
            versus!(t, :rotmat, () -> rotmat(m1, th, ax), () -> CS.rotmat(m1, th, ax))
            versus!(t, :rotvec, () -> rotvec(a, th, ax), () -> CS.rotvec(a, th, ax))
            p, q2 = rpair(rng)
            ia, ip = rand(rng, 0:4), rand(rng, 0:4)
            versus!(t, :twovec, () -> twovec(p, ia, q2, ip), () -> CS.twovec(p, ia, q2, ip))
            s1, s2 = rstate(rng), rand(rng) < 0.1 ? 2.0 * vcat(p, rvec(rng)) : rstate(rng)
            s1 = rand(rng) < 0.5 ? vcat(p, s1[4:6]) : s1
            versus!(t, :twovxf, () -> twovxf(s1, ia, s2, ip), () -> CS.twovxf(s1, ia, s2, ip))
        end
        versus!(t, :ident, () -> ident(), () -> CS.ident())
        @test invert(zeros(3, 3)) == zeros(3, 3)
        @test throws_kind(() -> invort(zeros(3, 3)), :zerolengthcolumn)
        @test throws_kind(() -> invort(Diagonal([1.0, 1e-320, 1.0])), :columntoosmall)
        @test throws_kind(() -> isrot(ident(), -1.0, 0.0), :valueoutofrange)
        @test isrot(ident(), 0.0, 0.0) && !isrot(2.0 * ident(), 0.5, 0.5)
        @test throws_kind(() -> twovec([1.0, 0, 0], 1, [2.0, 0, 0], 2), :dependentvectors)
        @test throws_kind(() -> twovec([1.0, 0, 0], 1, [0.0, 1, 0], 1), :undefinedframe)
        @test throws_kind(() -> twovec([1.0, 0, 0], 0, [0.0, 1, 0], 1), :badindex)
        @test rotate(0.0, 3) == ident() && rotate(0.3, 4) == rotate(0.3, 1)
    end

    @testset "derivatives of vector functions of a state" begin
        for _ in 1:n
            s1, s2 = rstate(rng), rstate(rng)
            if rand(rng) < 0.1          # parallel positions
                s2 = vcat(-2.0 * s1[1:3], s2[4:6])
            end
            versus!(t, :dvhat, () -> dvhat(s1), () -> CS.dvhat(s1))
            versus!(t, :dvdot, () -> dvdot(s1, s2), () -> CS.dvdot(s1, s2))
            versus!(t, :dvcrss, () -> dvcrss(s1, s2), () -> CS.dvcrss(s1, s2))
            versus!(t, :ducrss, () -> ducrss(s1, s2), () -> CS.ducrss(s1, s2))
            versus!(t, :dvnorm, () -> dvnorm(s1), () -> CS.dvnorm(s1))
            # `dvsep_` divides by the sine of the angle: positions of like magnitude
            # keep the quotient inside the range of the CSPICE overflow test
            versus!(t, :dvsep, () -> dvsep(s1, s2), () -> CS.dvsep(s1, s2))
        end
        @test dvhat([0.0, 0, 0, 1, 2, 3]) == [0.0, 0, 0, 1, 2, 3]
        @test dvnorm([0.0, 0, 0, 1, 2, 3]) == 0.0
        @test dvsep([1.0, 0, 0, 0, 1, 0], [2.0, 0, 0, 0, 0, 1]) == 0.0
    end

    @testset "general-dimension forms" begin
        vsepg_nan = 0
        vsepg_abs = 0.0
        for _ in 1:n
            k = rand(rng, 1:9)
            a = rmag(rng) * randn(rng, k)
            b = rand(rng) < 0.1 ? 2.0 * a : rand(rng) < 0.1 ? -a : rmag(rng) * randn(rng, k)
            rand(rng) < 0.05 && (a = zeros(k))
            x, y = rs(), rs()
            versus!(t, :vaddg, () -> vaddg(a, b), () -> CS.vaddg(a, b))
            versus!(t, :vsubg, () -> vsubg(a, b), () -> CS.vsubg(a, b))
            versus!(t, :vsclg, () -> vsclg(x, a), () -> CS.vsclg(x, a))
            versus!(t, :vminug, () -> vminug(a), () -> CS.vminug(a))
            versus!(t, :vequg, () -> vequg(a), () -> a)
            versus!(t, :vdotg, () -> vdotg(a, b), () -> CS.vdotg(a, b))
            versus!(t, :vnormg, () -> vnormg(a), () -> CS.vnormg(a))
            versus!(t, :vhatg, () -> vhatg(a), () -> CS.vhatg(a))
            versus!(t, :unormg, () -> unormg(a), () -> CS.unormg(a))
            versus!(t, :vdistg, () -> vdistg(a, b), () -> CS.vdistg(a, b))
            # `vsepg_c` has the two arcsine cases exchanged; the f2c routine has not
            versus!(t, :vsepg, () -> vsepg(a, b), () -> CS.vsepg_(a, b))
            sc = CS.vsepg(a, b)
            if isnan(sc)
                vsepg_nan += 1
            else
                vsepg_abs = max(vsepg_abs, abs(vsepg(a, b) - sc))
            end
            versus!(t, :vrelg, () -> vrelg(a, b), () -> CS.vrelg(a, b))
            versus!(t, :vprojg, () -> vprojg(a, b), () -> CS.vprojg(a, b))
            versus!(t, :vlcomg, () -> vlcomg(x, a, y, b), () -> CS.vlcomg(x, a, y, b))
            versus!(t, :vzerog, () -> vzerog(a), () -> CS.vzerog(a))
            l, p = rand(rng, 1:6), rand(rng, 1:6)
            m1 = rmag(rng) * randn(rng, k, l)
            m2 = rmag(rng) * randn(rng, l, p)
            m3 = rmag(rng) * randn(rng, k, p)
            m4 = rmag(rng) * randn(rng, p, l)
            v = rmag(rng) * randn(rng, l)
            versus!(t, :vtmvg, () -> vtmvg(a, m1, v), () -> CS.vtmvg(a, m1, v))
            versus!(t, :mxvg, () -> mxvg(m1, v), () -> CS.mxvg(m1, v))
            versus!(t, :mtxvg, () -> mtxvg(m1, a), () -> CS.mtxvg(m1, a))
            versus!(t, :mxmg, () -> mxmg(m1, m2), () -> CS.mxmg(m1, m2))
            versus!(t, :mtxmg, () -> mtxmg(m1, m3), () -> CS.mtxmg(m1, m3))
            versus!(t, :mxmtg, () -> mxmtg(m1, m4), () -> CS.mxmtg(m1, m4))
            versus!(t, :xposeg, () -> xposeg(m1), () -> CS.xposeg(m1))
        end
        @test throws_kind(() -> vaddg([1.0, 2.0], [1.0]), :badarg)
        @test throws_kind(() -> mxmg(zeros(2, 3), zeros(2, 3)), :badarg)
        # The accuracy defect of `vsepg_c`: for two directions 1e-9 rad apart it loses
        # half of the digits, and NeoSpice does not.
        u = [1.0, 0.0, 0.0, 0.0]
        w = [cos(1e-9), sin(1e-9), 0.0, 0.0]
        @test abs(vsepg(u, w) - 1e-9) < 1e-24
        # Over the random sample: `vsepg_c` is within 1e-7 rad where it gives a number,
        # and it gives NaN for some parallel pairs (the arcsine of a number above 1).
        @test vsepg_abs < 1e-7
        @test isnan(CS.vsepg([-2.790652279880394, -2.349350602390059, 0.1097765487434136,
            0.1663627783663699, 0.6147273839026778], 3.0 .* [-2.790652279880394, -2.349350602390059,
            0.1097765487434136, 0.1663627783663699, 0.6147273839026778]))
        measured(@sprintf("  vsepg against vsepg_c (CSPICE defect): largest difference %.3e rad, NaN from CSPICE in %d of %d",
            vsepg_abs, vsepg_nan, n))
        @test abs(CS.vsepg(u, w) - 1e-9) > 1e-12
        # static arguments give static results
        @test vaddg(SVector(1.0, 2.0), SVector(1.0, 2.0)) isa SVector{2,Float64}
        @test mxmg(SMatrix{2,2}(1.0, 2, 3, 4), SMatrix{2,2}(1.0, 2, 3, 4)) isa SMatrix{2,2,Float64}
    end

    @testset "differences against CSPICE" begin
        # 0 where the arithmetic is the one of the `_c` function; a few units of rounding
        # where a library function (sine, arcsine) or the f2c operation order enters.
        e = 4 * eps()
        report!(t, "vectors"; tol=0.0, tols=Dict(:vsep => e, :vrotv => e, :rotate => e, :rotmat => e,
            :rotvec => e, :twovxf => 1e-13, :ducrss => 1e-13, :dvsep => 1e-12, :vsepg => e))
    end

    @testset "number types: ForwardDiff and Taylor" begin
        a = SVector(1.2, -0.7, 2.1)
        b = SVector(0.3, 1.9, -0.4)
        # gradient of the norm is the unit vector; of the separation, by differences
        @test ForwardDiff.gradient(vnorm, a) ≈ vhat(a) rtol = 1e-14
        @test ForwardDiff.jacobian(vhat, a) ≈ (I - vhat(a) * vhat(a)') / vnorm(a) rtol = 1e-13
        g = ForwardDiff.gradient(x -> vsep(x, b), a)
        h = 1e-6
        fd = [(vsep(a + h * SVector(Float64(i == 1), Float64(i == 2), Float64(i == 3)), b) -
               vsep(a - h * SVector(Float64(i == 1), Float64(i == 2), Float64(i == 3)), b)) / 2h for i in 1:3]
        @test g ≈ fd rtol = 1e-8
        @test ForwardDiff.derivative(th -> vrotv(a, b, th), 0.4) ≈ vcrss(vhat(b), vrotv(a, b, 0.4)) rtol = 1e-13
        @test ForwardDiff.derivative(th -> rotate(th, 3), 0.4) ≈
              [-sin(0.4) cos(0.4) 0; -cos(0.4) -sin(0.4) 0; 0 0 0] rtol = 1e-14
        # the unit state of `dvhat` is the derivative of `vhat` along the velocity
        s = SVector(1.2, -0.7, 2.1, 0.3, 0.5, -0.9)
        dv = ForwardDiff.derivative(τ -> vhat(SVector(s[1], s[2], s[3]) + τ * SVector(s[4], s[5], s[6])), 0.0)
        @test dvhat(s)[SVector(4, 5, 6)] ≈ dv rtol = 1e-14
        s2 = SVector(-0.4, 0.8, 0.1, 0.2, -0.3, 0.6)
        pos(x, τ) = SVector(x[1], x[2], x[3]) + τ * SVector(x[4], x[5], x[6])
        @test dvsep(s, s2) ≈ ForwardDiff.derivative(τ -> vsep(pos(s, τ), pos(s2, τ)), 0.0) rtol = 1e-13
        @test dvnorm(s) ≈ ForwardDiff.derivative(τ -> vnorm(pos(s, τ)), 0.0) rtol = 1e-14
        @test dvdot(s, s2) ≈ ForwardDiff.derivative(τ -> vdot(pos(s, τ), pos(s2, τ)), 0.0) rtol = 1e-14
        @test ducrss(s, s2)[SVector(4, 5, 6)] ≈
              ForwardDiff.derivative(τ -> ucrss(pos(s, τ), pos(s2, τ)), 0.0) rtol = 1e-13
        # a Taylor time series passes through the smooth calls: r(t) = a + b t
        tt = NeoSpice.taylor_var(0.0, Val(2))
        r = SVector(a[1] + b[1] * tt, a[2] + b[2] * tt, a[3] + b[3] * tt)
        nrm = vnorm(r)
        @test nrm isa NeoSpice.Taylor{2,Float64}
        @test NeoSpice.derivatives(nrm)[1] == vnorm(a)
        @test NeoSpice.derivatives(nrm)[2] ≈ vdot(a, b) / vnorm(a) rtol = 1e-14
        @test NeoSpice.derivatives(vhat(r)[1])[2] ≈ dvhat(vcat(a, b))[4] rtol = 1e-14
        @test NeoSpice.derivatives(rotate(tt, 1)[2, 3])[2] == 1.0
        @test NeoSpice.derivatives(mxv(rotate(0.3 * tt, 3), r)[1])[1] == a[1]
        @test NeoSpice.derivatives(vcrss(r, b)[1])[1] == vcrss(a, b)[1]
    end

    @testset "inferred types and no allocation for static arguments" begin
        a = SVector(1.2, -0.7, 2.1)
        b = SVector(0.3, 1.9, -0.4)
        m = rrot(rng)
        s = SVector(1.2, -0.7, 2.1, 0.3, 0.5, -0.9)
        s2 = SVector(-0.4, 0.8, 0.1, 0.2, -0.3, 0.6)
        m6 = SMatrix{6,6}(randn(rng, 6, 6))
        V3, M3, V6, M6 = SVector{3,Float64}, SMatrix{3,3,Float64,9}, SVector{6,Float64}, SMatrix{6,6,Float64,36}
        calls = Any[(vadd, (a, b), V3), (vsub, (a, b), V3), (vscl, (2.0, a), V3), (vminus, (a,), V3),
            (vdot, (a, b), Float64), (vcrss, (a, b), V3), (vnorm, (a,), Float64), (vhat, (a,), V3),
            (unorm, (a,), Tuple{V3,Float64}), (vdist, (a, b), Float64), (vsep, (a, b), Float64),
            (vrel, (a, b), Float64), (vproj, (a, b), V3), (vperp, (a, b), V3),
            (vlcom, (1.0, a, 2.0, b), V3), (vlcom3, (1.0, a, 2.0, b, 3.0, a), V3),
            (vrotv, (a, b, 0.3), V3), (vzero, (a,), Bool), (vequ, (a,), V3),
            (vpack, (1.0, 2.0, 3.0), V3), (vupack, (a,), NTuple{3,Float64}), (ucrss, (a, b), V3),
            (vtmv, (a, m, b), Float64), (mxv, (m, a), V3), (mtxv, (m, a), V3), (mxm, (m, m), M3),
            (mtxm, (m, m), M3), (mxmt, (m, m), M3), (xpose, (m,), M3), (xpose6, (m6,), M6),
            (NeoSpice.det, (m,), Float64), (trace, (m,), Float64), (invert, (m,), M3),
            (invort, (m,), M3), (isrot, (m, 1e-10, 1e-10), Bool), (ident, (), M3),
            (rotate, (0.3, 2), M3), (rotmat, (m, 0.3, 2), M3), (rotvec, (a, 0.3, 2), V3),
            (twovec, (a, 1, b, 2), M3), (twovxf, (s, 1, s2, 2), M6), (dvhat, (s,), V6),
            (dvdot, (s, s2), Float64), (dvcrss, (s, s2), V6), (ducrss, (s, s2), V6),
            (dvnorm, (s,), Float64), (dvsep, (s, s2), Float64),
            (vaddg, (a, b), V3), (vsubg, (a, b), V3), (vsclg, (2.0, a), V3), (vminug, (a,), V3),
            (vequg, (a,), V3), (vdotg, (a, b), Float64), (vnormg, (a,), Float64), (vhatg, (a,), V3),
            (unormg, (a,), Tuple{V3,Float64}), (vdistg, (a, b), Float64), (vsepg, (a, b), Float64),
            (vrelg, (a, b), Float64), (vprojg, (a, b), V3), (vlcomg, (1.0, a, 2.0, b), V3),
            (vzerog, (a,), Bool), (vtmvg, (a, m, b), Float64), (mxvg, (m, a), V3), (mtxvg, (m, a), V3),
            (mxmg, (m, m), M3), (mtxmg, (m, m), M3), (mxmtg, (m, m), M3), (xposeg, (m,), M3)]
        for (f, args, T) in calls
            @test (@inferred f(args...)) isa T
            @test alloc_of(f, args...) == 0
        end
        # the scalar g forms do not allocate for ordinary vectors either
        va, vb = collect(a), collect(b)
        for f in (vdotg, vdistg, vsepg, vrelg)
            @test alloc_of(f, va, vb) == 0
        end
        @test alloc_of(vnormg, va) == 0 && alloc_of(vzerog, va) == 0
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
