# Group 3: time derivatives of states and rotations.
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
#
# Size of (c). A `Dual` nested k deep through the generic evaluators is a new compilation
# of the SPK and frame code for each k, and the cost grows fast with k (about 2 s, 3 s,
# 6 s and 16 s for the position). The default run therefore nests the generic path to
# `AD_GENERIC_ORDER` = 1 and the public calls to `AD_PUBLIC_ORDER` = 2; the exhaustive run
# (`ALLSPICE_TEST_FULL=1`) nests both to order 4. Depth 2 is the first depth at which the
# extension assembles a dual number from dual numbers; a deeper nesting runs the same
# methods on another type. Check (b) covers orders 1 to 4 of `spkder` and `pxder` in both
# runs.
#
# Size of (a). The per-epoch assertions of (a) use `EPOCHS_EXACT` epochs (2 in the default
# run, all 12 in the exhaustive run): the relation between the orders does not depend on
# the epoch, and the checks (b) and (c) use the larger sets.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff

const FD_CASES_INERTIAL = (h=2000.0, rtol=2e-4)
const FD_CASES_FAST = (h=20.0, rtol=2e-5)
const AD_RTOL = (2e-15, 1e-13, 1e-11, 1e-9, 1e-7)
const MAXORDER = 4
"Nesting depth of ForwardDiff through the generic evaluators, and through the public calls."
const AD_GENERIC_ORDER = nsample(1, MAXORDER)
const AD_PUBLIC_ORDER = nsample(2, MAXORDER)
"Number of epochs of the per-epoch assertions of check (a)."
const EPOCHS_EXACT = nsample(2, 12)

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

"The function `t -> d^N f / dt^N` by nested forward mode; `f` itself for `N == 0`."
ad_nth(f, ::Val{0}) = f
function ad_nth(f, ::Val{N}) where {N}
    g = ad_nth(f, Val(N - 1))
    return s -> ForwardDiff.derivative(g, s)
end

"""
    ad_jet(f, t, Val(N)) -> NTuple{N+1}

Nested forward-mode derivatives of `f` at `t`, orders 0 to `N`. Requires `0 <= N`.
"""
function ad_jet(f, t::Float64, ::Val{N}) where {N}
    @assert N isa Int && N >= 0
    jet = ntuple(k -> ad_nth(f, Val(k - 1))(t), Val(N + 1))
    @assert length(jet) == N + 1
    return jet
end

"Position through the generic evaluators (no analytical shortcut for a `Dual` epoch)."
generic_pos(ctx, targ, t, ref, obs) = invoke(Allspice.spk_jet,
    Tuple{Allspice.KernelSet,Int32,Real,Int32,Int32,Val{0}}, Allspice.kernelset(ctx),
    Int32(targ), t, Allspice.frame_id(Allspice.kernelset(ctx), ref), Int32(obs), Val(0))[1]

"Rotation through the generic evaluators (no analytical shortcut for a `Dual` epoch)."
function generic_rot(ctx, from, to, t)
    ks = Allspice.kernelset(ctx)
    return invoke(Allspice.frame_rotation, Tuple{Allspice.KernelSet,Int32,Int32,Real}, ks,
        Allspice.frame_id(ks, from), Allspice.frame_id(ks, to), t)
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

    @apart @testset "spkder: result types and consistency between orders" begin
        et = epochs[1]
        for N in 0:MAXORDER
            jet = spkder(ctx, 301, et, "J2000", 399, Val(N))
            @test jet isa NTuple{N + 1,SVector{3,Float64}}
        end
        # A lower order is a prefix of a higher order, in each class of frame.
        for (targ, obs, ref, _) in state_cases, et in epochs[1:nsample(1, 4)]
            full = spkder(ctx, targ, et, ref, obs, Val(MAXORDER))
            for N in 0:MAXORDER-1
                part = spkder(ctx, targ, et, ref, obs, Val(N))
                @test all(isapprox(part[k], full[k]; rtol=1e-13) for k in 1:N+1)
            end
        end
        # Name form equals the integer-ID form (the name form only resolves the names).
        @test spkder(ctx, "301", et, "IAU_EARTH", "399", Val(2)) ==
              spkder(ctx, 301, et, SPICE.namfrm("IAU_EARTH"), 399, Val(2))
    end

    @apart @testset "spkder: orders 0 and 1 are the spkezr state" begin
        worst = 0.0
        for (targ, obs, ref, _) in state_cases, et in epochs[1:EPOCHS_EXACT]
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

    @apart @testset "spkder: order k+1 against the central difference of order k" begin
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

    @apart @testset "spkder: nested ForwardDiff through the generic path" begin
        # The epoch that reaches Allspice is a Dual (nested up to `AD_PUBLIC_ORDER` deep).
        probe = spkpos(ctx, 301, ForwardDiff.Dual(epochs[1], 1.0), "J2000", "NONE", 399)[1]
        @test eltype(probe) <: ForwardDiff.Dual
        @test ForwardDiff.value.(probe) == spkpos(ctx, 301, epochs[1], "J2000", "NONE", 399)[1]
        for (targ, obs, ref, _) in state_cases
            generic = t -> generic_pos(ctx, targ, t, ref, obs)
            public = t -> spkpos(ctx, targ, t, ref, "NONE", obs)[1]
            for (label, order) in (("generic", AD_GENERIC_ORDER), ("public", AD_PUBLIC_ORDER))
                worst = zeros(order + 1)
                for et in epochs[1:6]
                    ad = label == "generic" ? collect(ad_jet(generic, et, Val(AD_GENERIC_ORDER))) :
                         collect(ad_jet(public, et, Val(AD_PUBLIC_ORDER)))
                    jet = spkder(ctx, targ, et, ref, obs, Val(MAXORDER))
                    worst .= max.(worst, [reldiff(jet[k], ad[k]) for k in 1:order+1])
                end
                measured(@sprintf("  AD %-7s %3d rel %3d %-11s orders 0..%d: %s", label, targ,
                    obs, ref, order, join((@sprintf("%.2e", w) for w in worst), " ")))
                @test all(worst[k] <= AD_RTOL[k] for k in 1:order+1)
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

    @apart @testset "pxder: orders 0 and 1 are pxform and the sxform block" begin
        worst = 0.0
        for N in 0:MAXORDER
            @test pxder(ctx, "J2000", "IAU_EARTH", epochs[1], Val(N)) isa NTuple{N + 1,SMatrix{3,3,Float64,9}}
        end
        for (from, to, _) in rot_cases, et in epochs[1:EPOCHS_EXACT]
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

    @apart @testset "pxder: order k+1 against the central difference of order k" begin
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

    @apart @testset "pxder: nested ForwardDiff through the generic path" begin
        for (from, to, _) in rot_cases
            generic = t -> generic_rot(ctx, from, to, t)
            public = t -> pxform(ctx, from, to, t)
            for (label, order) in (("generic", AD_GENERIC_ORDER), ("public", AD_PUBLIC_ORDER))
                worst = zeros(order + 1)
                for et in epochs[1:6]
                    ad = label == "generic" ? collect(ad_jet(generic, et, Val(AD_GENERIC_ORDER))) :
                         collect(ad_jet(public, et, Val(AD_PUBLIC_ORDER)))
                    jet = pxder(ctx, from, to, et, Val(MAXORDER))
                    worst .= max.(worst, [reldiff(jet[k], ad[k]) for k in 1:order+1])
                end
                measured(@sprintf("  AD %-7s %-11s to %-11s orders 0..%d: %s", label, from, to,
                    order, join((@sprintf("%.2e", w) for w in worst), " ")))
                @test all(worst[k] <= AD_RTOL[k] for k in 1:order+1)
            end
        end
        # pxder with a Dual epoch: the partial of order k is order k+1.
        et = epochs[3]
        jd = pxder(ctx, "J2000", "IAU_EARTH", ForwardDiff.Dual(et, 1.0), Val(2))
        jp = pxder(ctx, "J2000", "IAU_EARTH", et, Val(3))
        @test all(ForwardDiff.value.(jd[k]) == jp[k] for k in 1:3)
        @test all(ForwardDiff.partials.(jd[k], 1) == jp[k+1] for k in 1:3)
    end

    @apart @testset "calls without derivatives refuse a Dual epoch" begin
        # (fault fixed: limbpt, termpt, tangpt, occult, et2lst and fovtrg threw a
        # MethodError for a Dual epoch; they take a plain epoch, and a Dual now throws
        # SpiceError(:badarg))
        t = ForwardDiff.Dual(epochs[1], 1.0)
        z = [0.0, 0.0, 1.0]
        calls = (
            () -> limbpt(ctx, "TANGENT/ELLIPSOID", 399, t, "IAU_EARTH", "NONE", "CENTER", 301,
                z, 1.0, 1, 0.0, 0.0, 1),
            () -> termpt(ctx, "UMBRAL/TANGENT/ELLIPSOID", 10, 399, t, "IAU_EARTH", "NONE",
                "CENTER", 301, z, 1.0, 1, 0.0, 0.0, 1),
            () -> tangpt(ctx, "ELLIPSOID", 399, t, "IAU_EARTH", "NONE", "TANGENT POINT", 301,
                "J2000", z),
            () -> occult(ctx, 301, "ELLIPSOID", "IAU_MOON", 10, "ELLIPSOID", "IAU_SUN", "NONE",
                399, t),
            () -> et2lst(ctx, t, 399, 0.0, "PLANETOCENTRIC"),
            () -> fovtrg(ctx, -999, 301, "POINT", " ", "NONE", 399, t))
        @test all(f -> throws_kind(f, :badarg), calls)
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
