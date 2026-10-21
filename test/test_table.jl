# Group: ephemeris and rotation tables (`ephemeris_table`, `table_position`, `table_state`,
# `table_jet`, `table_to_spk`, `rotation_table`, `table_rotation`, `table_rotation_jet`).
#
# The oracle is CSPICE at random epochs of the table. Thresholds:
#   native grid (geometric, inertial, SPK type 2): the table is the polynomial of the
#     kernels, so the position agrees to a few units of rounding (measured 6e-16
#     relative) and the velocity to rounding times degree^2 (measured 8e-13);
#   interpolation: the position agrees within `maxerr` of the table plus the
#     difference between Allspice and CSPICE.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff
import Adapt

struct ToFloat32 end
Adapt.adapt_storage(::ToFloat32, x::Array{Float64}) = Float32.(x)

table_unchecked(tab, et) = @inbounds table_position(tab, et)
table_alloc(tab, et) = @allocated table_state(tab, et)

group("ephemeris table", Kernels.FULL_SET_FILES...) do
    rng = MersenneTwister(SEED + 41)
    set = Kernels.full_set()
    oracle_load(set)
    ctx = neo_context(set)
    t0, t1 = 7.0e8 + 1234.5, 7.0e8 + 200 * 86400.0
    rel(a, b) = norm(a - b) / norm(b)

    @apart @testset "native grid" begin
        # The Sun has records of 16 days in de440s.bsp; the Earth and the Moon, 4 days.
        for (targ, obs, ref, days) in [(301, 399, "J2000", 4), (399, 10, "ECLIPJ2000", 4),
                (4, 399, "J2000", 4), (10, 0, "J2000", 16), (399, 301, "ECLIPJ2000", 4)]
            tab = ephemeris_table(ctx, targ, t0, t1, ref, "NONE", obs)
            lo, hi = table_span(tab)
            @test lo <= t0 && t1 <= hi
            @test tab.dt == days * 86400.0 && 8 <= tab.ncoef <= 13 && length(tab.coef) == 3 * tab.ncoef * tab.nint
            @test (tab.target, tab.observer) == (targ, obs)
            ep = ev = 0.0
            for et in random_epochs(rng, 2000, lo, hi)
                ref_state, _ = SPICE.spkezr(string(targ), et, ref, "NONE", string(obs))
                x = table_state(tab, et)
                ep = max(ep, rel(x[1:3], ref_state[1:3]))
                ev = max(ev, rel(x[4:6], ref_state[4:6]))
            end
            measured(@sprintf("  table native %3d rel %3d %-10s nint=%d  position %.2e  velocity %.2e (relative)  maxerr %.2e km",
                targ, obs, ref, tab.nint, ep, ev, tab.maxerr))
            @test ep <= 5e-15 && ev <= 1e-11
            @test tab.maxerr <= 5e-15 * norm(table_position(tab, lo))
        end
        # Names and the default context give the same table.
        a = ephemeris_table(ctx, "MOON", t0, t1, "J2000", "NONE", "EARTH")
        @test a.coef == ephemeris_table(ctx, 301, t0, t1, "J2000", Aberration("NONE"), 399).coef
    end

    @apart @testset "interpolation" begin
        tab = ephemeris_table(ctx, 399, t0, t1, "MOON_PA", "NONE", 301; tol=1e-6)
        @test table_span(tab) == (t0, t0 + tab.nint * tab.dt) && table_span(tab)[2] ≈ t1
        @test tab.maxerr <= 1e-6 && tab.ncoef == 16
        e = maximum(et -> norm(table_position(tab, et) - SPICE.spkpos("399", et, "MOON_PA", "NONE", "301")[1]),
            random_epochs(rng, 2000, t0, t1))
        measured(@sprintf("  table MOON_PA tol=1e-6: nint=%d maxerr %.2e km, against CSPICE %.2e km", tab.nint, tab.maxerr, e))
        @test e <= 3e-6
        # Without a tolerance the function stops at the noise of the source.
        for (targ, obs, ref, ab, limit) in [(301, 399, "J2000", "LT+S", 1e-4), (4, 399, "J2000", "CN", 1e-3),
                (10, 399, "IAU_EARTH", "NONE", 2e-2), (301, 399, "IAU_EARTH", "XLT", 1e-4)]
            tab = ephemeris_table(ctx, targ, t0, t1, ref, ab, obs)
            e = maximum(et -> norm(table_position(tab, et) - SPICE.spkpos(string(targ), et, ref, ab, string(obs))[1]),
                random_epochs(rng, 2000, t0, t1))
            measured(@sprintf("  table %3d rel %3d %-9s %-4s nint=%d maxerr %.2e km, against CSPICE %.2e km",
                targ, obs, ref, ab, tab.nint, tab.maxerr, e))
            @test tab.maxerr <= limit && e <= 3 * limit
        end
        tab2 = ephemeris_table(ctx, 301, t0, t1, "IAU_EARTH", "NONE", 399; tol=1e-3, ncoef=8)
        @test tab2.ncoef == 8 && tab2.maxerr <= 1e-3
    end

    @apart @testset "errors" begin
        @test throws_kind(() -> ephemeris_table(ctx, 301, t0, t1, "J2000", "LT+S", 399; tol=1e-9, maxintervals=64), :noconvergence)
        @test throws_kind(() -> ephemeris_table(ctx, 301, t1, t0, "J2000", "NONE", 399), :badarg)
        @test throws_kind(() -> ephemeris_table(ctx, 301, t0, t0, "J2000", "NONE", 399), :badarg)
        @test throws_kind(() -> ephemeris_table(ctx, 301, t0, Inf, "J2000", "NONE", 399), :badarg)
        @test throws_kind(() -> ephemeris_table(ctx, 301, t0, t1, "J2000", "NONE", 399; ncoef=1), :badarg)
        @test throws_kind(() -> ephemeris_table(ctx, 301, t0, t1, "J2000", "NONE", 399; tol=0.0), :badarg)
        @test throws_kind(() -> ephemeris_table(ctx, 301, t0, t1, "J2000", "NONE", 399; maxintervals=0), :badarg)
        @test throws_kind(() -> ephemeris_table(ctx, 301, t0, t1, "J2000", "BAD", 399), :invalidoption)
        @test_throws SpiceError ephemeris_table(ctx, 301, t0, 1e11, "J2000", "NONE", 399)
        @test_throws SpiceError ephemeris_table(ctx, 499, t0, t1, "J2000", "NONE", 399)
        @test throws_kind(() -> Allspice.EphemerisTable(zeros(5), 0.0, 1.0, 1, 2, 1, 2, 1, 0.0), :badarg)
        tab = ephemeris_table(ctx, 301, t0, t1, "J2000", "NONE", 399)
        lo, hi = table_span(tab)
        @test throws_kind(() -> table_position(tab, lo - 1.0), :spkinsuffdata)
        @test throws_kind(() -> table_state(tab, hi + 1.0), :spkinsuffdata)
        @test throws_kind(() -> table_jet(tab, NaN, Val(2)), :spkinsuffdata)
        @test table_position(tab, lo) isa SVector{3,Float64} && table_position(tab, hi) isa SVector{3,Float64}
        if Base.JLOptions().check_bounds != 1
            @test all(isfinite, table_unchecked(tab, hi + 1.0))   # extrapolation, by contract
        end
        @test occursin("EphemerisTable{Float64}", sprint(show, tab))
    end

    @apart @testset "derivatives, allocation, types" begin
        tab = ephemeris_table(ctx, 301, t0, t1, "J2000", "NONE", 399)
        et = t0 + 4.56e6
        jet = table_jet(tab, et, Val(3))
        @test jet[1] == table_position(tab, et) && vcat(jet[1], jet[2]) == table_state(tab, et)
        ref_jet = spkder(ctx, 301, et, "J2000", 399, Val(3))
        @test all(k -> rel(jet[k], ref_jet[k]) <= (1e-14, 1e-11, 1e-9, 1e-7)[k], 1:4)
        d1 = ForwardDiff.derivative(t -> table_position(tab, t), et)
        d2 = ForwardDiff.derivative(s -> ForwardDiff.derivative(t -> table_position(tab, t), s), et)
        @test rel(d1, jet[2]) <= 1e-13 && rel(d2, jet[3]) <= 1e-11
        @test (@inferred table_state(tab, et)) isa SVector{6,Float64}
        table_alloc(tab, et)
        @test table_alloc(tab, et) == 0

        single = Allspice.EphemerisTable{Float32}(tab)
        @test single.coef isa Vector{Float32} && table_span(single) == table_span(tab)
        @test table_position(single, et) isa SVector{3,Float32}
        @test rel(table_position(single, et), jet[1]) <= 1e-5

        moved = Adapt.adapt(ToFloat32(), tab)
        @test moved.coef isa Vector{Float32} && moved.coef == single.coef
        @test (moved.t0, moved.dt, moved.nint, moved.ncoef, moved.maxerr) ==
              (tab.t0, tab.dt, tab.nint, tab.ncoef, tab.maxerr)
        small = Allspice.EphemerisTable(SVector{6,Float64}(1, 2, 3, 4, 5, 6), 0.0, 2.0, 1, 2, 1, 2, 1, 0.0)
        @test isbits(small) && table_position(small, 1.5) == SVector(2.0, 5.0, 8.0)
    end

    @apart @testset "export as an SPK segment" begin
        dir = mktempdir()
        tab = ephemeris_table(ctx, 301, t0, t1, "ECLIPJ2000", "NONE", 399)
        lo, hi = table_span(tab)
        p = joinpath(dir, "table.bsp")
        h = spkopn(p, "TABLE", 0)
        table_to_spk(ctx, h, tab, "moon from the table")
        fit = ephemeris_table(ctx, 399, t0, t1, "MOON_PA", "NONE", 301; tol=1e-6)
        @test throws_kind(() -> table_to_spk(ctx, h, fit, "x"; first=t0 - 10.0), :coveragegap)
        # Regression: table_to_spk wrote a frame ID that no context knows into the segment;
        # spkw02_c signals SPICE(INVALIDREFFRAME) for a frame that it cannot resolve.
        unknown = Allspice.EphemerisTable(fit.coef, fit.t0, fit.dt, fit.nint, fit.ncoef,
            fit.target, fit.observer, -999, fit.maxerr)
        hc = SPICE.spkopn(joinpath(dir, "table_c.bsp"), "TABLE", 0)
        ccall((:spkw02_c, SPICE.libcspice), Cvoid, (Cint, Cint, Cint, Cstring, Cdouble, Cdouble,
            Cstring, Cdouble, Cint, Cint, Ptr{Cdouble}, Cdouble), hc, -9, 399, "NO_SUCH_FRAME",
            0.0, 2.0, "x", 2.0, 1, 1, collect(1.0:6.0), 0.0)
        @test cspice_kind() === :invalidrefframe
        SPICE.dafcls(hc)
        @test throws_kind(() -> table_to_spk(ctx, h, unknown, "x"), :invalidrefframe)
        # The frame of `fit`, MOON_PA, is defined by a frame kernel of `ctx` only.
        @test throws_kind(() -> table_to_spk(SpiceContext(), h, fit, "x"), :invalidrefframe)
        spkcls(h)
        # A second file: the Earth from the Moon in the same file would close a loop of bodies.
        p2 = joinpath(dir, "table2.bsp")
        h2 = spkopn(p2, "TABLE", 0)
        table_to_spk(ctx, h2, Allspice.EphemerisTable{Float32}(fit), "earth, single precision")
        spkcls(h2)
        own = neo_context([p])
        SPICE.kclear(); SPICE.furnsh(Kernels.LSK); SPICE.furnsh(p)
        ep = ev = ec = 0.0
        for et in random_epochs(rng, 500, lo, hi)
            x = table_state(tab, et)
            y = spkezr(own, 301, et, "ECLIPJ2000", "NONE", 399)[1]
            c = SPICE.spkezr("301", et, "ECLIPJ2000", "NONE", "399")[1]
            ep = max(ep, rel(y[1:3], x[1:3])); ev = max(ev, rel(y[4:6], x[4:6])); ec = max(ec, rel(c, y))
        end
        measured(@sprintf("  table as SPK type 2: position %.1e, velocity %.1e against the table; %.1e against CSPICE on the file", ep, ev, ec))
        @test ep <= 5e-15 && ev <= 1e-11 && ec <= 1e-13
        # The segment keeps the bodies, the frame and the span of the table.
        seg = only(Allspice.kernelset(own).spk[Int32(301)])
        @test (seg.center, seg.dtype, seg.t0, seg.t1) == (399, 2, lo, hi)
        own2 = neo_context([p2])
        seg2 = only(Allspice.kernelset(own2).spk[Int32(399)])
        @test (seg2.center, seg2.frame, seg2.dtype) == (301, fit.frame, 2)
        d = Allspice.kernelset(own2).daf[seg2.file]
        @test d[seg2.baddr:seg2.baddr+3*fit.ncoef+1][3:end] == Float64.(Float32.(fit.coef[1:3*fit.ncoef]))
        oracle_load(set)
    end

    @apart @testset "rotation table" begin
        oracle_load(set)
        days = FULL ? 30 : 6
        r0, r1 = t0, t0 + days * 86400.0
        frob(a, b) = norm(a - b)
        for (from, to, limit) in [("J2000", "IAU_EARTH", 1e-9), ("IAU_MOON", "J2000", 1e-10),
                ("MOON_PA", "IAU_EARTH", 1e-9), ("ECLIPJ2000", "J2000", 1e-14)]
            tab = rotation_table(ctx, from, to, r0, r1)
            @test table_span(tab)[1] == r0 && table_span(tab)[2] ≈ r1
            @test length(tab.coef) == 9 * tab.ncoef * tab.nint && tab.maxerr <= limit
            e = eo = 0.0
            for et in random_epochs(rng, 500, r0, r1)
                m = table_rotation(tab, et)
                e = max(e, frob(m, SPICE.pxform(from, to, et)))
                eo = max(eo, frob(m * m', I))
            end
            measured(@sprintf("  rotation table %-10s to %-9s nint=%d maxerr %.1e, against CSPICE %.1e, orthogonality %.1e",
                from, to, tab.nint, tab.maxerr, e, eo))
            @test e <= 3 * limit && eo <= 6 * limit
        end
        tab = rotation_table(ctx, "J2000", "IAU_EARTH", r0, r1; tol=1e-8)
        @test tab.maxerr <= 1e-8
        et = r0 + 0.37 * (r1 - r0)
        jet = table_rotation_jet(tab, et, Val(2))
        refj = pxder(ctx, "J2000", "IAU_EARTH", et, Val(2))
        @test jet[1] == table_rotation(tab, et)
        @test frob(jet[2], refj[2]) <= 1e-9 && frob(jet[3], refj[3]) <= 1e-11
        d1 = ForwardDiff.derivative(t -> table_rotation(tab, t), et)
        @test frob(d1, jet[2]) <= 1e-18
        @test (@inferred table_rotation(tab, et)) isa SMatrix{3,3,Float64,9}
        rot_alloc(tb, t) = @allocated table_rotation_jet(tb, t, Val(1))
        rot_alloc(tab, et)
        @test rot_alloc(tab, et) == 0
        single = Allspice.RotationTable{Float32}(tab)
        @test table_rotation(single, et) isa SMatrix{3,3,Float32,9}
        @test frob(table_rotation(single, et), jet[1]) <= 1e-5
        moved = Adapt.adapt(ToFloat32(), tab)
        @test moved.coef == single.coef && (moved.from, moved.to, moved.nint) == (tab.from, tab.to, tab.nint)
        @test occursin("RotationTable{Float64}", sprint(show, tab))
        @test rotation_table(ctx, 1, "IAU_EARTH", r0, r1; tol=1e-8).coef == tab.coef
        @test throws_kind(() -> table_rotation(tab, r1 + 1.0), :spkinsuffdata)
        @test throws_kind(() -> rotation_table(ctx, "J2000", "IAU_EARTH", r1, r0), :badarg)
        @test throws_kind(() -> rotation_table(ctx, "J2000", "IAU_EARTH", r0, r1; ncoef=1), :badarg)
        @test throws_kind(() -> rotation_table(ctx, "J2000", "IAU_EARTH", r0, r1; tol=-1.0), :badarg)
        @test throws_kind(() -> rotation_table(ctx, "J2000", "IAU_EARTH", r0, r1; maxintervals=0), :badarg)
        @test throws_kind(() -> rotation_table(ctx, "J2000", "IAU_EARTH", r0, r1; tol=1e-15, maxintervals=8), :noconvergence)
        @test throws_kind(() -> Allspice.RotationTable(zeros(5), 0.0, 1.0, 1, 2, 1, 2, 0.0), :badarg)
        @test_throws SpiceError rotation_table(ctx, "J2000", "NO_SUCH_FRAME", r0, r1)
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
