# Group 3: time derivatives of states and rotations.
#
# Feature file: features/derivative_orders.feature
#
# Three independent checks of `spkder` and `pxder`:
#   (a) orders 0 and 1 are the `spkezr` state (the `pxform` matrix and the derivative
#       block of `sxform`);
#   (b) order k+1 is the central finite difference of order k;
#   (c) every order is the nested `ForwardDiff.derivative` of the position function
#       (of `pxform`), evaluated through the generic path with a `Dual` epoch. The
#       public calls answer a `Dual` epoch from the analytical derivatives (the
#       ForwardDiff extension), so the generic evaluators are reached with `invoke`
#       (`generic_pos`, `generic_rot`); the public calls are checked as well.
#
# Tolerance of (b). The central difference (f(t+h) - f(t-h)) / 2h has the truncation
# error h^2/6 * |f'''| and the rounding error eps * |f| / h. With tau the time scale of
# the motion (|f^(k)| ~ |f| / tau^k) the relative errors are (h/tau)^2 / 6 and
# eps * tau / h:
#
#   case                         tau (s)   h (s)   truncation   rounding   tolerance
#   inertial frame, MOON_PA      >= 3.7e5  2000    <= 5e-6      <= 4e-14   2e-4
#   IAU_EARTH (tau = 1 / omega)  1.4e4     20      4e-7         1.4e-13    2e-5
#
# The tolerance has a margin of 40 above the truncation error, because the high
# derivatives of a Chebyshev record carry the fit error of the record amplified by
# (degree^2 / record length)^k. All three points t-h, t, t+h lie in one record
# (`record_safe_epochs`), so no difference is taken across a record boundary.
#
# Tolerance of (c). Both paths differentiate the same polynomial, so they differ only in
# rounding, which the same amplification raises with the order. Starting values per
# order 0..4, relative to the norm of the derivative: see `AD_RTOL`.

import ForwardDiff

const FD_CASES_INERTIAL = (h=2000.0, rtol=2e-4)
const FD_CASES_FAST = (h=20.0, rtol=2e-5)
const AD_RTOL = (2e-15, 1e-13, 1e-11, 1e-9, 1e-7)
const MAXORDER = 4

"Relative difference `|a - b| / |b|`; for `b == 0` the absolute difference."
function reldiff(a, b)
    nb = norm(b)
    return nb == 0 ? norm(a - b) : norm(a - b) / nb
end

"""
    fd_errors(jet, et, h, N) -> Vector{Float64}

For `jet(t) -> NTuple{N+1}` return, for k = 0..N-1, the relative difference between
order k+1 at `et` and the central difference of order k with step `h`.
"""
function fd_errors(jet, et::Float64, h::Float64, N::Int)
    @assert h > 0 && N >= 1
    j0 = jet(et); jp = jet(et + h); jm = jet(et - h)
    @assert length(j0) == N + 1
    return [reldiff((jp[k+1] - jm[k+1]) / (2h), j0[k+2]) for k in 0:N-1]
end

"Nested forward-mode derivatives of `f` at `t`, orders 0 to 4, as a tuple."
function ad_jet(f, t::Float64)
    d1(s) = ForwardDiff.derivative(f, s)
    d2(s) = ForwardDiff.derivative(d1, s)
    d3(s) = ForwardDiff.derivative(d2, s)
    d4(s) = ForwardDiff.derivative(d3, s)
    return (f(t), d1(t), d2(t), d3(t), d4(t))
end

"Position through the generic evaluators (no analytical shortcut for a `Dual` epoch)."
generic_pos(ctx, targ, t, ref, obs) = invoke(NeoSpice.spk_jet,
    Tuple{NeoSpice.KernelSet,Int32,Real,Int32,Int32,Val{0}}, NeoSpice.kernelset(ctx),
    Int32(targ), t, NeoSpice.frame_id(NeoSpice.kernelset(ctx), ref), Int32(obs), Val(0))[1]

"Rotation through the generic evaluators (no analytical shortcut for a `Dual` epoch)."
function generic_rot(ctx, from, to, t)
    ks = NeoSpice.kernelset(ctx)
    return invoke(NeoSpice.frame_rotation, Tuple{NeoSpice.KernelSet,Int32,Int32,Real}, ks,
        NeoSpice.frame_id(ks, from), NeoSpice.frame_id(ks, to), t)
end

group("3. derivatives", Kernels.FULL_SET_FILES...) do
    rng = MersenneTwister(SEED + 3)
    set = Kernels.full_set()
    oracle_load(set)
    ctx = neo_context(set)
    epochs = record_safe_epochs(rng, 12)

    state_cases = [
        (301, 399, "J2000", FD_CASES_INERTIAL), (399, 10, "J2000", FD_CASES_INERTIAL),
        (10, 301, "J2000", FD_CASES_INERTIAL), (5, 399, "ECLIPJ2000", FD_CASES_INERTIAL),
        (199, 4, "J2000", FD_CASES_INERTIAL),
        (301, 399, "IAU_EARTH", FD_CASES_FAST), (10, 399, "IAU_EARTH", FD_CASES_FAST),
        (399, 301, "MOON_PA", FD_CASES_INERTIAL), (10, 301, "MOON_PA", FD_CASES_INERTIAL),
        (399, 301, "MOON_ME", FD_CASES_INERTIAL),
    ]

    @testset "spkder: result types and consistency between orders" begin
        et = epochs[1]
        for N in 0:MAXORDER
            jet = spkder(ctx, 301, et, "J2000", 399, Val(N))
            @test jet isa NTuple{N + 1,SVector{3,Float64}}
        end
        # A lower order is a prefix of a higher order.
        for (targ, obs, ref, _) in state_cases, et in epochs[1:4]
            full = spkder(ctx, targ, et, ref, obs, Val(MAXORDER))
            for N in 0:MAXORDER-1
                part = spkder(ctx, targ, et, ref, obs, Val(N))
                @test all(isapprox(part[k], full[k]; rtol=1e-13) for k in 1:N+1)
            end
            # Name form equals the integer-ID form.
            @test spkder(ctx, string(targ), et, ref, string(obs), Val(2)) ==
                  spkder(ctx, targ, et, SPICE.namfrm(ref), obs, Val(2))
        end
    end

    @testset "spkder: orders 0 and 1 are the spkezr state" begin
        worst = 0.0
        for (targ, obs, ref, _) in state_cases, et in epochs
            state, _ = spkezr(ctx, targ, et, ref, "NONE", obs)
            cstate, _ = SPICE.spkezr(string(targ), et, ref, "NONE", string(obs))
            r, v = spkder(ctx, targ, et, ref, obs, Val(1))
            @test r == state[SOneTo(3)]
            @test v == state[SVector(4, 5, 6)]
            (r0,) = spkder(ctx, targ, et, ref, obs, Val(0))
            @test isapprox(r0, r; rtol=1e-14)
            worst = max(worst, reldiff(r, cstate[1:3]), reldiff(v, cstate[4:6]))
        end
        measured(@sprintf("  %-58s max=%.3e", "spkder orders 0, 1 against CSPICE, relative", worst))
        @test worst <= 1e-9
    end

    @testset "spkder: order k+1 against the central difference of order k" begin
        for (targ, obs, ref, fd) in state_cases
            worst = zeros(MAXORDER)
            for et in epochs
                e = fd_errors(t -> spkder(ctx, targ, t, ref, obs, Val(MAXORDER)), et, fd.h,
                    MAXORDER)
                worst .= max.(worst, e)
            end
            measured(@sprintf("  FD  %3d rel %3d %-11s h=%-6g orders 1..4: %s (tol %.0e)",
                targ, obs, ref, fd.h, join((@sprintf("%.2e", w) for w in worst), " "),
                fd.rtol))
            @test all(<=(fd.rtol), worst)
        end
    end

    @testset "spkder: nested ForwardDiff through the generic path" begin
        # The epoch that reaches NeoSpice is a Dual (nested up to four deep).
        probe = spkpos(ctx, 301, ForwardDiff.Dual(epochs[1], 1.0), "J2000", "NONE", 399)[1]
        @test eltype(probe) <: ForwardDiff.Dual
        @test ForwardDiff.value.(probe) == spkpos(ctx, 301, epochs[1], "J2000", "NONE", 399)[1]
        for (targ, obs, ref, _) in state_cases
            paths = ("generic" => t -> generic_pos(ctx, targ, t, ref, obs),
                "public" => t -> spkpos(ctx, targ, t, ref, "NONE", obs)[1])
            for (label, pos) in paths
                worst = zeros(MAXORDER + 1)
                for et in epochs[1:6]
                    ad = ad_jet(pos, et)
                    jet = spkder(ctx, targ, et, ref, obs, Val(MAXORDER))
                    worst .= max.(worst, [reldiff(jet[k], ad[k]) for k in 1:MAXORDER+1])
                end
                measured(@sprintf("  AD %-7s %3d rel %3d %-11s orders 0..4: %s", label, targ,
                    obs, ref, join((@sprintf("%.2e", w) for w in worst), " ")))
                @test all(worst[k] <= AD_RTOL[k] for k in 1:MAXORDER+1)
            end
        end
        # The velocity part of spkezr with a Dual epoch carries the acceleration.
        et = epochs[2]
        acc = ForwardDiff.derivative(t -> spkezr(ctx, 301, t, "J2000", "NONE", 399)[1][SVector(4, 5, 6)], et)
        @test isapprox(acc, spkder(ctx, 301, et, "J2000", 399, Val(2))[3]; rtol=1e-10)
    end

    rot_cases = [
        ("J2000", "IAU_EARTH", FD_CASES_FAST), ("IAU_EARTH", "J2000", FD_CASES_FAST),
        ("J2000", "IAU_MOON", FD_CASES_INERTIAL), ("J2000", "MOON_PA", FD_CASES_INERTIAL),
        ("MOON_ME", "ECLIPJ2000", FD_CASES_INERTIAL), ("IAU_EARTH", "MOON_PA", FD_CASES_FAST),
        ("IAU_MARS", "IAU_JUPITER", FD_CASES_FAST),
    ]

    @testset "pxder: orders 0 and 1 are pxform and the sxform block" begin
        worst = 0.0
        for (from, to, _) in rot_cases, et in epochs
            for N in 0:MAXORDER
                @test pxder(ctx, from, to, et, Val(N)) isa NTuple{N + 1,SMatrix{3,3,Float64,9}}
            end
            R, dR = pxder(ctx, from, to, et, Val(1))
            X = sxform(ctx, from, to, et)
            @test isapprox(R, pxform(ctx, from, to, et); atol=4eps())
            @test dR == X[SVector(4, 5, 6), SOneTo(3)]
            Xc = SPICE.sxform(from, to, et)
            worst = max(worst, reldiff(dR, Xc[4:6, 1:3]))
        end
        measured(@sprintf("  %-58s max=%.3e", "pxder order 1 against CSPICE sxform, relative", worst))
        @test worst <= 1e-7
        # Between inertial frames every derivative is zero.
        jet = pxder(ctx, "J2000", "GALACTIC", epochs[1], Val(3))
        @test all(all(iszero, jet[k]) for k in 2:4)
    end

    @testset "pxder: order k+1 against the central difference of order k" begin
        for (from, to, fd) in rot_cases
            worst = zeros(MAXORDER)
            for et in epochs
                e = fd_errors(t -> pxder(ctx, from, to, t, Val(MAXORDER)), et, fd.h, MAXORDER)
                worst .= max.(worst, e)
            end
            measured(@sprintf("  FD  %-11s to %-11s h=%-6g orders 1..4: %s (tol %.0e)", from,
                to, fd.h, join((@sprintf("%.2e", w) for w in worst), " "), fd.rtol))
            @test all(<=(fd.rtol), worst)
        end
    end

    @testset "pxder: nested ForwardDiff through the generic path" begin
        for (from, to, _) in rot_cases
            paths = ("generic" => t -> generic_rot(ctx, from, to, t),
                "public" => t -> pxform(ctx, from, to, t))
            for (label, rot) in paths
                worst = zeros(MAXORDER + 1)
                for et in epochs[1:6]
                    ad = ad_jet(rot, et)
                    jet = pxder(ctx, from, to, et, Val(MAXORDER))
                    worst .= max.(worst, [reldiff(jet[k], ad[k]) for k in 1:MAXORDER+1])
                end
                measured(@sprintf("  AD %-7s %-11s to %-11s orders 0..4: %s", label, from, to,
                    join((@sprintf("%.2e", w) for w in worst), " ")))
                @test all(worst[k] <= AD_RTOL[k] for k in 1:MAXORDER+1)
            end
        end
        # pxder with a Dual epoch: the partial of order k is order k+1.
        et = epochs[3]
        jd = pxder(ctx, "J2000", "IAU_EARTH", ForwardDiff.Dual(et, 1.0), Val(2))
        jp = pxder(ctx, "J2000", "IAU_EARTH", et, Val(3))
        @test all(ForwardDiff.value.(jd[k]) == jp[k] for k in 1:3)
        @test all(ForwardDiff.partials.(jd[k], 1) == jp[k+1] for k in 1:3)
    end
end
