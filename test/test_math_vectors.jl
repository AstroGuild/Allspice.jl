# Group 12: vector and matrix calls, and the base 16 strings of doubles, against CSPICE.
#
# Feature file: features/math_vectors.feature
#
# Each function runs on several hundred seeded random inputs with magnitudes from 1e-10
# to 1e10, mixed with the degenerate inputs (zero vector, unit axes, equal, parallel and
# opposite pairs, very small and very large components). The limit of each function is
# 0 where Allspice has the operation order of the `_c` function, and a few units of
# rounding where CSPICE runs the f2c routine with another order.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff
isdefined(@__MODULE__, :CS) || include("mathoracle.jl")

"""
    dvsep_overflows(s1, s2) -> Bool

True if `dvsep_c` must not be called for the states `s1` and `s2`. The f2c routine
`dvsep_` divides a numerator (the sum of two dot products of the unit positions and
their derivatives) by the norm of the cross product of the unit positions. If the
numerator is infinite, `zzdiv_` signals `SPICE(NUMERICOVERFLOW)` and formats the
numerator for the message with `errdp_` (zzdiv.c, line 406). `dpstr_` cannot format an
infinite number: it converts the logarithm of the number to an integer and uses the
result as a subscript of the table `ipower` (dpstr.c, lines 368 to 397), and the
subscript check of f2c then stops the whole process. Contract: the function makes the
same CSPICE calls as `dvsep_` before the quotient (`dvhat_c`, `vcrss_c`, `vdot_c`), which
signal no error, and it does not call Allspice.
"""
function dvsep_overflows(s1, s2)
    u1, u2 = CS.dvhat(s1), CS.dvhat(s2)
    iszero(CS.vcrss(u1[1:3], u2[1:3])) && return false    # `dvsep_` returns 0 before the quotient
    return isinf(CS.vdot(u1[1:3], u2[4:6]) + CS.vdot(u1[4:6], u2[1:3]))
end

"""
    ducrss_overflows(s1, s2) -> Bool

True if the derivative part of `ducrss_c` is not defined for the states `s1` and `s2`:
`ducrss_` divides each state by its largest position component and takes the cross
product and its derivative with `dvcrss_`, and that intermediate result has an element
that is not finite (a velocity above about 1e8 with the position 1e-300 becomes
infinite, and infinity times zero is NaN). CSPICE then returns numbers without meaning:
the macro `max` of f2c (f2c.h, line 392) drops a NaN operand, so `vperp_` (vperp.c, lines
178 and 182) can take a NaN velocity for a zero velocity. Allspice returns NaN.
Contract: the function repeats the scaling of `ducrss_` and calls `dvcrss_c` only.
"""
function ducrss_overflows(s1, s2)
    function scaled(s)
        f = maximum(abs, s[1:3])
        return f > 0.0 ? (1.0 / f) .* collect(Float64, s) : collect(Float64, s)
    end
    return !all(isfinite, CS.dvcrss(scaled(s1), scaled(s2)))
end

"""
    twovxf_overflows(axdef, indexa, plndef, indexp) -> Bool

True if one of the two calls of `ducrss_` in `zztwovxf_` has an intermediate result that
is not finite (see `ducrss_overflows`). Axis numbers that are not valid give `false`:
`twovxf_c` signals an error for them before any arithmetic.
"""
function twovxf_overflows(axdef, indexa, plndef, indexp)
    (indexa in 1:3 && indexp in 1:3 && indexa != indexp) || return false
    if indexp == indexa % 3 + 1
        return ducrss_overflows(axdef, plndef) || ducrss_overflows(CS.ducrss(axdef, plndef), axdef)
    end
    return ducrss_overflows(plndef, axdef) || ducrss_overflows(axdef, CS.ducrss(plndef, axdef))
end

"""
    oracle_hx2dp(s) -> Union{Float64,Symbol,String}

`hx2dp_c` of CSPICE: the number, the short message of a signalled error as a lower-case
symbol, or the text that CSPICE writes into `errmsg` with its error flag.
"""
function oracle_hx2dp(s::AbstractString)
    num = Ref{Cdouble}(NaN)
    err = Ref{Cint}(0)
    msg = zeros(UInt8, 256)
    ccall((:hx2dp_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Ref{Cdouble}, Ref{Cint},
        Ptr{UInt8}), s, length(msg), num, err, msg)
    try
        CS.check()
    catch e
        e isa CS.CError || rethrow()
        return e.kind
    end
    return err[] != 0 ? GC.@preserve(msg, unsafe_string(pointer(msg))) : num[]
end

"`hx2dp(s)` in the form of `oracle_hx2dp`: the number, an error kind, or the message of `:badnumber`."
function neo_hx2dp(s::AbstractString)
    try
        return hx2dp(s)
    catch e
        e isa SpiceError || rethrow()
        return e.kind === :badnumber ? e.msg : e.kind
    end
end

"""
A random string in the neighbourhood of the base 16 form: an optional sign, up to 33
digits (with now and then a character that is not a digit), the caret or not, and an
exponent of up to 4 digits, with blanks around. The exponent is kept short because
CSPICE scales by one factor of 16 per unit of the exponent.
"""
function random_hex(rng::AbstractRNG)
    pick(chars) = chars[rand(rng, 1:length(chars))]
    digit() = rand(rng) < 0.03 ? pick("G -.") : pick("0123456789ABCDEFabcdef")
    io = IOBuffer()
    rand(rng) < 0.2 && print(io, " "^rand(rng, 1:2))
    rand(rng) < 0.4 && print(io, pick("+-"))
    for _ in 1:(rand(rng) < 0.1 ? rand(rng, 30:33) : rand(rng, 0:14))
        print(io, digit())
    end
    if rand(rng) < 0.95
        print(io, '^')
        rand(rng) < 0.1 && print(io, ' ')
        rand(rng) < 0.5 && print(io, pick("+-"))
        for _ in 1:rand(rng, 0:4)
            print(io, digit())
        end
    end
    rand(rng) < 0.2 && print(io, " "^rand(rng, 1:2))
    return String(take!(io))
end

group("12. math: vectors and matrices") do
    rng = MersenneTwister(SEED + 12)
    n = nsample(100, 400)
    t = Tally()
    rs() = rmag(rng) * randn(rng)
    rm3() = rand(rng) < 0.2 ? Matrix(rrot(rng)) : rmag(rng) * randn(rng, 3, 3)
    rax() = rand(rng, -4:7)

    @apart @testset "three-dimensional vectors" begin
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
            # `frame_c` overwrites its input with the unit vector; SPICE.jl passes a copy.
            versus!(t, :frame, () -> frame(a), () -> SPICE.frame(collect(Float64, a)))
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
        @test frame(z3) == (e1, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0)) == SPICE.frame(zeros(3))
        @test throws_kind(() -> frame([1.0, 2.0]), :badarg)
        @test vrel(z3, z3) == 0.0 && vrel(e1, e1) == 0.0
        @test vnorm(SVector(1e300, 1e300, 1e300)) == CS.vnorm([1e300, 1e300, 1e300]) < Inf
        @test vnorm(SVector(1e-300, 1e-300, 0.0)) > 0.0
        @test throws_kind(() -> vadd([1.0, 2.0], [1.0, 2.0, 3.0]), :badarg)
        @test vadd([1, 2, 3], [1, 1, 1]) === SVector(2.0, 3.0, 4.0)
    end

    @apart @testset "three-dimensional matrices" begin
        nskip = 0
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
            versus!(t, :invstm, () -> invstm(m6), () -> CS.invstm(m6))
            versus!(t, :det, () -> Allspice.det(m1), () -> CS.det(m1))
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
            # Not compared where an intermediate result of CSPICE is not finite.
            if twovxf_overflows(s1, ia, s2, ip)
                nskip += 1
            else
                versus!(t, :twovxf, () -> twovxf(s1, ia, s2, ip), () -> CS.twovxf(s1, ia, s2, ip))
            end
        end
        measured("  twovxf: $nskip of $n inputs not compared (intermediate result of CSPICE not finite)")
        # Fixed case of the operation order: `zztwovxf_` calls the f2c routines `dvhat_`
        # and `ducrss_`, which scale by a division.
        sa = [-2.9, 8.2, -11.9, 2.1, -12.4, -6.1]
        sb = [-11.6, -0.3, 8.3, 7.8, 12.3, -15.7]
        @test twovxf(sa, 1, sb, 2) == CS.twovxf(sa, 1, sb, 2)
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
        # the inverse of a state transformation
        x = eul2xf([0.4, -1.1, 2.0, 1e-3, -2e-4, 5e-5], 3, 1, 3)
        @test invstm(x) * x ≈ I atol = 1e-15
    end

    @apart @testset "derivatives of vector functions of a state" begin
        nguard = 0
        nskip = 0
        for _ in 1:n
            s1, s2 = rstate(rng), rstate(rng)
            if rand(rng) < 0.1          # parallel positions
                s2 = vcat(-2.0 * s1[1:3], s2[4:6])
            end
            versus!(t, :dvhat, () -> dvhat(s1), () -> CS.dvhat(s1))
            versus!(t, :dvdot, () -> dvdot(s1, s2), () -> CS.dvdot(s1, s2))
            versus!(t, :dvcrss, () -> dvcrss(s1, s2), () -> CS.dvcrss(s1, s2))
            # Where an intermediate result of CSPICE is not finite, the derivative part
            # has no meaning and the unit vector only is compared.
            last = ducrss_overflows(s1, s2) ? 3 : 6
            nskip += last == 3
            versus!(t, :ducrss, () -> collect(ducrss(s1, s2))[1:last], () -> collect(CS.ducrss(s1, s2))[1:last])
            versus!(t, :dvnorm, () -> dvnorm(s1), () -> CS.dvnorm(s1))
            # A velocity of 1e300 with a small position gives an infinite numerator in
            # `dvsep_`. CSPICE then stops the process (see `dvsep_overflows`), so the
            # oracle is not called and the expected result is the error that `zzdiv_`
            # signals for this input.
            guard = dvsep_overflows(s1, s2)
            nguard += guard
            versus!(t, :dvsep, () -> dvsep(s1, s2), () -> guard ? throw(CS.CError(:numericoverflow)) : CS.dvsep(s1, s2))
        end
        measured("  dvsep: $nguard of $n inputs not given to dvsep_c (infinite numerator, CSPICE stops the process)")
        measured("  ducrss: derivative of $nskip of $n inputs not compared (intermediate result of CSPICE not finite)")
        # Fixed cases of the operation order: `ducrss_` and `dvsep_` call the f2c routine
        # `dvhat_`, which scales by a division, and not `dvhat_c`, which multiplies by a
        # reciprocal. The reciprocal of a velocity below 5.6e-309 is infinite (first
        # case); the two orders differ by rounding (second and third case), and a
        # position of 1e-300 makes that difference as large as the result (fourth case).
        sden1 = [-6.354774078047717e-7, 8.943378680402032e-6, 1.6011695882103261e-6, 0.0, 0.0, 0.0]
        sden2 = [1e300, 1e300, 1e300, -2.572684309046176e-10, -4.428686661044401e-10, -2.7522373404770423e-10]
        sa = [-2.9, 8.2, -11.9, 2.1, -12.4, -6.1]
        sb = [-11.6, -0.3, 8.3, 7.8, 12.3, -15.7]
        sc = [10.3, 1.3, -6.5, 8.4, -21.3, 10.2]
        sd = [13.7, 4.7, -6.6, -7.4, 11.8, -4.1]
        sbig = [2.353009784756466e-6, 8.590299545730118e-6, -4.739435602996215e-6, 3.9803597213088137e8,
            2.7987192714291748e7, -3.951825994570644e8]
        stiny = [1e-300, 0.0, 0.0, 1.0, 0.0, 0.0]
        @test all(isfinite, ducrss(sden1, sden2))
        @test all(p -> !ducrss_overflows(p...) && ducrss(p...) == CS.ducrss(p...), ((sden1, sden2), (sa, sb), (sbig, stiny)))
        @test dvsep(sc, sd) == CS.dvsep(sc, sd)
        # Fixed case of the guard: the scaled velocity of the first state is infinite.
        sover = [1e-300, 0.0, 0.0, 6.340909641677687e7, -6.894680105462036e7, 7.179757811367899e8]
        sunit = [0.0, 0.0, -1.0, 1e-300, 0.0, 0.0]
        @test ducrss_overflows(sover, sunit) && twovxf_overflows(sover, 2, sunit, 1) &&
              ducrss(sover, sunit)[1:3] == CS.ducrss(sover, sunit)[1:3] &&
              all(isnan, ducrss(sover, sunit)[4:6]) && all(iszero, CS.ducrss(sover, sunit)[4:6])
        # Fixed case of the guard: the derivative of the unit position is infinite.
        sinf = [1e-300, 0.0, 0.0, 0.0, 1e300, 0.0]
        sfin = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
        @test dvsep_overflows(sinf, sfin) && dvsep_overflows(sfin, sinf) &&
              !dvsep_overflows([1.0, 0, 0, 0, 1, 0], [2.0, 0, 0, 0, 0, 1])
        @test all(p -> throws_kind(() -> dvsep(p...), :numericoverflow), ((sinf, sfin), (sfin, sinf)))
        @test dvhat([0.0, 0, 0, 1, 2, 3]) == [0.0, 0, 0, 1, 2, 3]
        @test dvnorm([0.0, 0, 0, 1, 2, 3]) == 0.0
        @test dvsep([1.0, 0, 0, 0, 1, 0], [2.0, 0, 0, 0, 0, 1]) == 0.0
    end

    @apart @testset "general-dimension forms" begin
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
        # half of the digits, and Allspice does not.
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

    @apart @testset "differences against CSPICE" begin
        # 0 where the arithmetic is the one of the CSPICE function; a few units of
        # rounding where a library function (sine, arcsine) enters.
        e = 4 * eps()
        report!(t, "vectors"; tol=0.0, tols=Dict(:vsep => e, :vrotv => e, :rotate => e, :rotmat => e,
            :rotvec => e, :vsepg => e))
    end

    @apart @testset "number types: ForwardDiff and Taylor" begin
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
        tt = Allspice.taylor_var(0.0, Val(2))
        r = SVector(a[1] + b[1] * tt, a[2] + b[2] * tt, a[3] + b[3] * tt)
        nrm = vnorm(r)
        @test nrm isa Allspice.Taylor{2,Float64}
        @test Allspice.derivatives(nrm)[1] == vnorm(a)
        @test Allspice.derivatives(nrm)[2] ≈ vdot(a, b) / vnorm(a) rtol = 1e-14
        @test Allspice.derivatives(vhat(r)[1])[2] ≈ dvhat(vcat(a, b))[4] rtol = 1e-14
        @test Allspice.derivatives(rotate(tt, 1)[2, 3])[2] == 1.0
        @test Allspice.derivatives(mxv(rotate(0.3 * tt, 3), r)[1])[1] == a[1]
        @test Allspice.derivatives(vcrss(r, b)[1])[1] == vcrss(a, b)[1]
    end

    @apart @testset "base 16 strings of doubles: dp2hx and hx2dp" begin
        specials = [0.0, -0.0, 1.0, -1.0, 672.0, -11.0, 0.1, 16.0^-3, nextfloat(0.0),
            -nextfloat(0.0), floatmin(), prevfloat(floatmin()), floatmax(), -floatmax(), NaN]
        nbad = 0
        for x in vcat(specials, [rmag(rng) * randn(rng) for _ in 1:n])
            h = dp2hx(x)
            nbad += h != SPICE.dp2hx(x)
            # Every finite double has an exact string, so the round trip is the identity;
            # the string of -0.0 is the one of 0.0.
            isnan(x) || (nbad += hx2dp(h) !== (iszero(x) ? 0.0 : x))
        end
        @test nbad == 0
        @test throws_kind(() -> dp2hx(Inf), :valueoutofrange)   # dp2hx_c does not return
        # The error messages of CSPICE, each limit of the exponent and of the mantissa, and
        # a random sample around the form.
        strings = ["2A^3", "-B^1", "", " ", "   1^1  ", "1^ 2", "1^", "1", "^1", "-^0", "- 1^1",
            "1 ^1", "+1^1", "1^-", "1^G", "G^1", "1\t^1", "0123456789ABCDEF0123456789ABCDE^0",
            "0123456789ABCDEF0123456789ABCDEF^0", "1^7FFFFFFF", "1^80000000", "1^-80000001",
            "-1^7FFFFFFF", "FFFFFFFFFFFFF8^100", "FFFFFFFFFFFFF8^101", "-FFFFFFFFFFFFF8^100",
            "-FFFFFFFFFFFFF8^101", "4^-10C", "8^-10E", "abc^-1f", "-8^-80000000"]
        nbad = 0
        for str in vcat(strings, [random_hex(rng) for _ in 1:4n])
            got, want = neo_hx2dp(str), oracle_hx2dp(str)
            same = got isa Float64 && want isa Float64 ? got === want : isequal(got, want)
            same || @info "hx2dp" str got want
            nbad += !same
        end
        @test nbad == 0
        # CSPICE scales a zero 2^31 times here; the result is zero by arithmetic.
        @test hx2dp("0^7FFFFFFF") === 0.0 && hx2dp("1^-7FFFFFFF") === 0.0
    end

    @apart @testset "inferred types and no allocation for static arguments" begin
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
            (mtxm, (m, m), M3), (mxmt, (m, m), M3), (xpose, (m,), M3), (xpose6, (m6,), M6), (invstm, (m6,), M6),
            (Allspice.det, (m,), Float64), (trace, (m,), Float64), (invert, (m,), M3),
            (invort, (m,), M3), (isrot, (m, 1e-10, 1e-10), Bool), (ident, (), M3),
            (rotate, (0.3, 2), M3), (rotmat, (m, 0.3, 2), M3), (rotvec, (a, 0.3, 2), V3),
            (twovec, (a, 1, b, 2), M3), (twovxf, (s, 1, s2, 2), M6), (frame, (a,), NTuple{3,V3}),
            (dvhat, (s,), V6),
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
