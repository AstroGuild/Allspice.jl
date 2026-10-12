# Group 13: rotation calls (Euler angles, axis and angle, quaternions, state
# transformations) against CSPICE.
#
# Feature file: features/math_rotations.feature
#
# The random rotations come from Euler angles over the whole circle. The degenerate
# inputs are the identity, the half turns about each axis, Euler sequences with the
# middle angle at an end of its range, matrices that are not rotations, and axis numbers
# outside 1:3.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of NeoSpice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`NEOSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0

import ForwardDiff
isdefined(@__MODULE__, :CS) || include("mathoracle.jl")

"Rotations that sit on the branch boundaries of `m2q`, `m2eul` and `raxisa`."
function special_rotations()
    out = Any[Matrix(1.0I, 3, 3)]
    for d in ([1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [-1.0, -1.0, 1.0])
        push!(out, Matrix(Diagonal(d)))                    # half turns about x, y, z
    end
    for ax in 1:3, ang in (pi / 2, -pi / 2, pi, 1e-9, pi - 1e-9, 2pi / 3)
        push!(out, Matrix(CS.rotate(ang, ax)))
    end
    for (a, b, c) in ((3, 1, 3), (1, 2, 3), (3, 2, 1), (2, 3, 2), (1, 3, 2)), mid in (0.0, pi, pi / 2, -pi / 2)
        push!(out, Matrix(CS.eul2m(0.7, mid, -1.1, a, b, c)))
    end
    push!(out, [0.0 1.0 0.0; 0.0 0.0 1.0; 1.0 0.0 0.0])     # third of a turn about (1,1,1)
    push!(out, [0.0 1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 -1.0])    # half turn about (1,1,0)
    return out
end

group("13. math: rotations") do
    rng = MersenneTwister(SEED + 13)
    n = 400
    t = Tally()
    rang() = 2pi * rand(rng) - pi
    specials = special_rotations()
    "A rotation: random, special, or (rarely) not a rotation at all."
    function rmat()
        u = rand(rng)
        u < 0.2 && return rand(rng, specials)
        u < 0.25 && return Matrix(rrot(rng)) .* (1.0 + 0.3 * randn(rng))    # scaled
        u < 0.3 && return randn(rng, 3, 3)
        u < 0.35 && return Matrix(rrot(rng)) + 1e-9 * randn(rng, 3, 3)      # slightly off
        return Matrix(rrot(rng))
    end

    @testset "Euler angles" begin
        for _ in 1:n
            a3, a2, a1 = rang(), rang(), rang()
            rand(rng) < 0.15 && (a2 = rand(rng, (0.0, pi, -pi, pi / 2, -pi / 2)))
            # mostly valid axes; sometimes one outside 1:3
            ax = rand(rng) < 0.1 ? rand(rng, 0:4, 3) : rand(rng, 1:3, 3)
            versus!(t, :eul2m, () -> eul2m(a3, a2, a1, ax...), () -> CS.eul2m(a3, a2, a1, ax...))
            r = rmat()
            versus!(t, :m2eul, () -> m2eul(r, ax...), () -> CS.m2eul(r, ax...))
            # a matrix built with the same axes, so that the degenerate branch is met
            if all(in(1:3), ax) && ax[2] != ax[1] && ax[2] != ax[3]
                r2 = CS.eul2m(a3, a2, a1, ax...)
                versus!(t, :m2eul, () -> m2eul(r2, ax...), () -> CS.m2eul(r2, ax...))
                back = eul2m(m2eul(r2, ax...)..., ax...)
                @test maximum(abs, back - r2) < 1e-14
            end
            e = SVector(a3, a2, a1, randn(rng), randn(rng), randn(rng))
            # `eul2xf_` reads a table with an axis number outside 1:3 before its
            # check, so CSPICE is called with valid axes only
            axv = rand(rng, 1:3, 3)
            versus!(t, :eul2xf, () -> eul2xf(e, axv...), () -> CS.eul2xf(e, axv...))
            x = Matrix(CS.eul2xf(e, axv...))
            versus!(t, :xf2eul, () -> xf2eul(x, ax...), () -> CS.xf2eul(x, ax...))
            if axv[2] != axv[1] && axv[2] != axv[3]
                versus!(t, :xf2eul, () -> xf2eul(x, axv...), () -> CS.xf2eul(x, axv...))
                e2, unique = xf2eul(x, axv...)
                @test !unique || maximum(abs, eul2xf(e2, axv...) - x) < 1e-9 * max(1.0, maximum(abs, x))
            end
            xr = Matrix(CS.rav2xf(rmat(), randn(rng, 3)))
            versus!(t, :xf2eul, () -> xf2eul(xr, ax...), () -> CS.xf2eul(xr, ax...))
        end
        @test throws_kind(() -> eul2m(0.1, 0.2, 0.3, 3, 1, 4), :badaxisnumbers)
        @test throws_kind(() -> m2eul(ident(), 3, 3, 1), :badaxisnumbers)
        @test throws_kind(() -> m2eul(2.0 * ident(), 3, 1, 3), :notarotation)
        @test throws_kind(() -> eul2xf(zeros(6), 0, 1, 3), :badaxisnumbers)
        # the degenerate case sets the first angle to zero and keeps the rotation
        @test m2eul(eul2m(0.7, 0.0, -1.1, 3, 1, 3), 3, 1, 3)[1] == 0.0
        @test m2eul(ident(), 3, 1, 3) == (0.0, 0.0, 0.0)
        @test m2eul(ident(), 1, 2, 3) == (0.0, 0.0, 0.0)
        @test xf2eul(eul2xf([0.7, 0.0, -1.1, 0.1, 0.2, 0.3], 3, 1, 3), 3, 1, 3)[2] == false
        @test xf2eul(eul2xf([0.7, 0.4, -1.1, 0.1, 0.2, 0.3], 3, 1, 3), 3, 1, 3)[2] == true
    end

    @testset "axis and angle, quaternions" begin
        for _ in 1:n
            ax = rvec(rng)
            th = 20 * rand(rng) - 10
            versus!(t, :axisar, () -> axisar(ax, th), () -> CS.axisar(ax, th))
            r = rmat()
            versus!(t, :raxisa, () -> raxisa(r), () -> CS.raxisa(r))
            versus!(t, :m2q, () -> m2q(r), () -> CS.m2q(r))
            versus!(t, :sharpr, () -> sharpr(r), () -> CS.sharpr_(r))
            q1 = rand(rng) < 0.05 ? zeros(4) : rmag(rng) * randn(rng, 4)
            q2 = rand(rng) < 0.3 ? Vector(CS.m2q(rrot(rng))) : randn(rng, 4)
            dq = randn(rng, 4)
            versus!(t, :q2m, () -> q2m(q1), () -> CS.q2m(q1))
            versus!(t, :q2m, () -> q2m(q2), () -> CS.q2m(q2))
            versus!(t, :qxq, () -> qxq(q1, q2), () -> CS.qxq(q1, q2))
            versus!(t, :qdq2av, () -> qdq2av(q1, dq), () -> CS.qdq2av(q1, dq))
            av = rvec(rng)
            versus!(t, :rav2xf, () -> rav2xf(r, av), () -> CS.rav2xf(r, av))
            x = rmag(rng) * randn(rng, 6, 6)
            versus!(t, :xf2rav, () -> xf2rav(x), () -> CS.xf2rav(x))
            # round trips on a true rotation
            rr = Matrix(rrot(rng))
            @test maximum(abs, q2m(m2q(rr)) - rr) < 1e-14
            axis, angle = raxisa(rr)
            @test maximum(abs, axisar(axis, angle) - rr) < 1e-14
            @test 0.0 <= angle <= pi
            w = randn(rng, 3)
            r2, w2 = xf2rav(rav2xf(rr, w))
            @test r2 == rr && maximum(abs, w2 - w) < 1e-14
            @test abs(vnormg(m2q(rr)) - 1.0) < 4eps() && m2q(rr)[1] >= 0.0
            qa, qb = m2q(rr), m2q(rrot(rng))
            @test maximum(abs, q2m(qxq(qa, qb)) - q2m(qa) * q2m(qb)) < 1e-14
        end
        for r in specials
            versus!(t, :m2q, () -> m2q(r), () -> CS.m2q(r))
            versus!(t, :raxisa, () -> raxisa(r), () -> CS.raxisa(r))
            for ax in ((3, 1, 3), (1, 2, 3), (3, 2, 1), (2, 1, 2), (2, 3, 1), (1, 3, 2))
                versus!(t, :m2eul, () -> m2eul(r, ax...), () -> CS.m2eul(r, ax...))
            end
        end
        @test m2q(ident()) == [1.0, 0.0, 0.0, 0.0]
        # a half turn has the scalar part 0, which is not "> 0": CSPICE then negates
        @test m2q(Diagonal([1.0, -1.0, -1.0])) == [0.0, -1.0, 0.0, 0.0] == CS.m2q(Diagonal([1.0, -1.0, -1.0]))
        @test m2q(Diagonal([-1.0, 1.0, -1.0])) == [0.0, 0.0, -1.0, 0.0]
        @test m2q(Diagonal([-1.0, -1.0, 1.0])) == [0.0, 0.0, 0.0, -1.0]
        @test raxisa(ident()) == ([0.0, 0.0, 1.0], 0.0)
        @test raxisa(Diagonal([1.0, -1.0, -1.0])) == ([-1.0, 0.0, 0.0], Float64(pi))
        @test axisar([0.0, 0.0, 0.0], 1.0) == ident()
        @test throws_kind(() -> m2q(zeros(3, 3)), :notarotation)
        @test throws_kind(() -> raxisa(3.0 * ident()), :notarotation)
        @test throws_kind(() -> q2m([1.0, 0.0, 0.0]), :badarg)
        @test q2m(zeros(4)) == ident()
    end

    @testset "differences against CSPICE" begin
        # The sine, cosine, arcsine and arctangent of Julia and of the C library differ
        # by one unit of rounding for some arguments, and a matrix product carries that
        # on. An Euler angle near a degenerate case has the error of the matrix divided
        # by the sine of the middle angle.
        report!(t, "rotations"; tol=0.0, tols=Dict(:eul2m => 1e-15, :m2eul => 4e-15, :eul2xf => 4e-15,
            :xf2eul => 1e-13, :axisar => 1e-15, :raxisa => 1e-15, :qdq2av => 1e-15))
    end

    @testset "number types: ForwardDiff and Taylor" begin
        e0 = SVector(0.7, 0.4, -1.1)
        de = SVector(0.1, -0.2, 0.3)
        for ax in ((3, 1, 3), (1, 2, 3), (3, 2, 1), (2, 1, 2))
            # the lower left block of `eul2xf` is the time derivative of `eul2m`
            dr = ForwardDiff.derivative(τ -> eul2m((e0 + τ * de)..., ax...), 0.0)
            x = eul2xf(vcat(e0, de), ax...)
            @test x[SVector(4, 5, 6), SVector(1, 2, 3)] ≈ dr rtol = 1e-13
            # the same derivative from a Taylor series of the angles
            tt = NeoSpice.taylor_var(0.0, Val(1))
            rt = eul2m((e0[1] + de[1] * tt), (e0[2] + de[2] * tt), (e0[3] + de[3] * tt), ax...)
            @test map(v -> NeoSpice.derivatives(v)[2], rt) ≈ dr rtol = 1e-13
            # the Jacobian of `m2eul` inverts the one of `eul2m`
            j1 = ForwardDiff.jacobian(v -> vec(eul2m(v..., ax...)), e0)
            j2 = ForwardDiff.jacobian(m -> SVector(m2eul(reshape(m, 3, 3), ax...)), vec(Matrix(eul2m(e0..., ax...))))
            @test j2 * j1 ≈ I atol = 1e-12
        end
        # angular velocity: the quaternion rate of a rotation about a fixed axis
        axis = vhat([0.3, -0.5, 0.8])
        rate = 0.37
        qf(τ) = m2q(xpose(axisar(axis, 0.6 + rate * τ)))      # frame rotation: the transpose
        q0 = qf(0.0)
        dq = ForwardDiff.derivative(qf, 0.0)
        @test qdq2av(q0, dq) ≈ rate * axis rtol = 1e-13
        rf(τ) = xpose(axisar(axis, 0.6 + rate * τ))
        x = NeoSpice._m_xform(rf(0.0), ForwardDiff.derivative(rf, 0.0))
        @test xf2rav(x)[2] ≈ rate * axis rtol = 1e-13
        @test rav2xf(rf(0.0), rate * axis) ≈ x rtol = 1e-13
        # Taylor series through the quaternion calls
        tt = NeoSpice.taylor_var(0.0, Val(2))
        qt = m2q(eul2m(0.3 + 0.1 * tt, 0.2 * tt, -0.4 + 0.5 * tt, 3, 1, 3))
        @test qt isa SVector{4,NeoSpice.Taylor{2,Float64}}
        @test map(v -> NeoSpice.derivatives(v)[1], qt) == m2q(eul2m(0.3, 0.0, -0.4, 3, 1, 3))
        back = q2m(qt)
        ref = eul2m(0.3 + 0.1 * tt, 0.2 * tt, -0.4 + 0.5 * tt, 3, 1, 3)
        @test maximum(abs(NeoSpice.derivatives(back[k])[j] - NeoSpice.derivatives(ref[k])[j]) for k in 1:9, j in 1:3) < 1e-14
        @test NeoSpice.derivatives(qxq(qt, qt)[1])[1] ≈ 2 * m2q(eul2m(0.3, 0.0, -0.4, 3, 1, 3))[1]^2 - 1
    end

    @testset "inferred types and no allocation for static arguments" begin
        r = rrot(rng)
        q = m2q(r)
        v = SVector(0.3, -0.5, 0.8)
        e = SVector(0.7, 0.4, -1.1, 0.1, -0.2, 0.3)
        x = eul2xf(e, 3, 1, 3)
        V3, V4, V6 = SVector{3,Float64}, SVector{4,Float64}, SVector{6,Float64}
        M3, M6 = SMatrix{3,3,Float64,9}, SMatrix{6,6,Float64,36}
        calls = Any[(eul2m, (0.1, 0.2, 0.3, 3, 1, 3), M3), (m2eul, (r, 3, 1, 3), NTuple{3,Float64}),
            (m2eul, (r, 1, 2, 3), NTuple{3,Float64}), (eul2xf, (e, 3, 1, 3), M6),
            (xf2eul, (x, 3, 1, 3), Tuple{V6,Bool}), (axisar, (v, 0.3), M3),
            (raxisa, (r,), Tuple{V3,Float64}), (m2q, (r,), V4), (q2m, (q,), M3), (qxq, (q, q), V4),
            (qdq2av, (q, q), V3), (rav2xf, (r, v), M6), (xf2rav, (x,), Tuple{M3,V3}), (sharpr, (r,), M3)]
        for (f, args, T) in calls
            @test (@inferred f(args...)) isa T
            @test alloc_of(f, args...) == 0
        end
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
