# Group 9: spacecraft clocks (SCLK type 1) against CSPICE.
#
# The SCLK kernel is written by `test/ckfixtures.jl`: clock -77 (two fields, two
# partitions, four coefficient records, parallel time TDB) and clock -78 (three fields
# with a first-field offset, three partitions, parallel time TDT).
#
# Error measures: seconds for ET, ticks for encoded SCLK; strings must be equal.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

isdefined(@__MODULE__, :CKFix) || include("ckfixtures.jl")

# ET of the clocks is below 4e9 s (one unit in the last place 5e-7 s) at the end of the
# last partition and below 2e5 s (3e-11 s) in the part that the CK tests use. The two
# engines do the same operations in the same order, so the measured difference is 0;
# the thresholds leave room for a `sin` that differs in the last place (TDT clock).
const SCLK_TOL_ET = 1e-7          # seconds, over the whole clock range
const SCLK_TOL_TICKS = 1e-6       # ticks, continuous conversion over the whole clock range

"Largest encoded SCLK of the clock with partitions `pstart`, `pstop`."
sclk_mxtick(pstart, pstop) = sum(round.(pstop .- pstart))

"True if both engines throw, or both return and `same(a, b)` holds."
function agree(same, neo, oracle)
    a = try
        neo()
    catch e
        e isa Allspice.SpiceError || rethrow()
        nothing
    end
    b = try
        oracle()
    catch e
        e isa SPICE.SpiceError || rethrow()
        nothing
    end
    (a === nothing || b === nothing) && return a === nothing && b === nothing
    return same(a, b)
end

alloc_sct2e(ks, sc, x) = @allocated sct2e(ks, sc, x)
alloc_sce2c(ks, sc, x) = @allocated sce2c(ks, sc, x)
alloc_sce2t(ks, sc, x) = @allocated sce2t(ks, sc, x)

group("9. spacecraft clocks against CSPICE", Kernels.LSK) do
    rng = MersenneTwister(SEED + 9)
    set = [Kernels.LSK, CKFix.sclk_kernel()]
    oracle_load(set)
    ctx = neo_context(set)
    ks = Allspice.kernelset(ctx)
    clocks = ((CKFix.SC_TDB, CKFix.TDB_PSTART, CKFix.TDB_PSTOP, CKFix.TDB_COEF),
        (CKFix.SC_TDT, CKFix.TDT_PSTART, CKFix.TDT_PSTOP, CKFix.TDT_COEF))

    @apart @testset "clock data" begin
        for (sc, pstart, pstop, _) in clocks
            # CSPICE has no C wrapper of `sctype_`; the Fortran-style function is called.
            @test sctype(ctx, sc) == ccall((:sctype_, SPICE.libcspice), Cint, (Ref{Cint},), sc) == 1
            ps, pe = scpart(ctx, sc)
            cs, ce = SPICE.scpart(sc)
            @test ps == cs == pstart
            @test pe == ce == pstop
        end
        @test sort(collect(keys(ks.sclk))) == [CKFix.SC_TDT, CKFix.SC_TDB]
    end

    for (sc, pstart, pstop, coef) in clocks
        mx = sclk_mxtick(pstart, pstop)
        # Random ticks over the whole range and over the part with the coefficient
        # records, and the special values: range ends, partition totals, record ticks.
        ticks = Float64[mx * rand(rng) for _ in 1:1500]
        append!(ticks, [coef[end, 1] * 1.2 * rand(rng) for _ in 1:1500])
        totals = cumsum(round.(pstop .- pstart))
        for s in (0.0, mx, totals..., coef[:, 1]...), d in (-1.0, -0.25, 0.0, 0.25, 1.0)
            0 <= s + d <= mx && push!(ticks, s + d)
        end
        det = Float64[]; dc = Float64[]; dt = Float64[]; rt = Float64[]
        strings_ok = true
        for t in ticks
            et = sct2e(ctx, sc, t)
            etc = SPICE.sct2e(sc, t)
            push!(det, abs(et - etc))
            push!(dc, abs(sce2c(ctx, sc, etc) - SPICE.sce2c(sc, etc)))
            push!(dt, abs(sce2t(ctx, sc, etc) - SPICE.sce2t(sc, etc)))
            push!(rt, abs(sce2c(ctx, sc, et) - t))
            s = scdecd(ctx, sc, t)
            strings_ok &= s == SPICE.scdecd(sc, t)
            strings_ok &= sce2s(ctx, sc, etc) == SPICE.sce2s(sc, etc)
            strings_ok &= scencd(ctx, sc, s) == SPICE.scencd(sc, s)
            strings_ok &= scs2e(ctx, sc, s) == SPICE.scs2e(sc, s)
            fields = last(split(s, '/'))
            strings_ok &= sctiks(ctx, sc, fields) == SPICE.sctiks(sc, fields)
            # `scfmt` of the tick count inside the partition
            k = sctiks(ctx, sc, fields)
            strings_ok &= scfmt(ctx, sc, k) == SPICE.scfmt(sc, k) == fields
        end
        a = report_quantiles("clock $sc: sct2e, seconds", det)
        b = report_quantiles("clock $sc: sce2c, ticks", dc)
        c = report_quantiles("clock $sc: sce2t, ticks", dt)
        r = report_quantiles("clock $sc: sce2c(sct2e(t)) - t, ticks", rt)
        @test a.max <= SCLK_TOL_ET
        @test b.max <= SCLK_TOL_TICKS
        @test c.max == 0.0
        # The round trip is exact to rounding, except in the 50 s jump of the parallel
        # time of clock -77 at its partition boundary, where two tick values share an ET.
        @test r.q50 <= SCLK_TOL_TICKS
        @test strings_ok
    end

    @apart @testset "clock strings as a user writes them" begin
        cases = ((CKFix.SC_TDB, ["1/100.5", "100.5", "1/0000000100.005", "2/2000.000",
                "2/ 2000 . 17", "  1/  77 12", "9999.255", "2/9999:255", "1/100", "100,3",
                "2/2000-1", "1/.7", "1/7.", "2000.001", "1.5/100.5"]),
            (CKFix.SC_TDT, ["1/1:0:0", "1:59:49", "1/000000854:20:01", "2/3:00:00",
                "3/1:0:0", "2/5", "1/7 30", "1/7::3", "1/7:5:", "12:30"]))
        for (sc, list) in cases, s in list
            @test agree(==, () -> scencd(ctx, sc, s), () -> SPICE.scencd(sc, s))
            @test agree(==, () -> scs2e(ctx, sc, s), () -> SPICE.scs2e(sc, s))
            f = String(last(split(s, '/')))
            @test agree(==, () -> sctiks(ctx, sc, f), () -> SPICE.sctiks(sc, f))
        end
        # At least one string of each clock is valid, so the comparison is not vacuous.
        @test scencd(ctx, CKFix.SC_TDB, "2/2000.001") == SPICE.scencd(CKFix.SC_TDB, "2/2000.001")
        @test scencd(ctx, CKFix.SC_TDT, "1/1:0:0") == 0.0
    end

    @apart @testset "errors" begin
        sc = CKFix.SC_TDB
        mx = sclk_mxtick(CKFix.TDB_PSTART, CKFix.TDB_PSTOP)
        @test throws_kind(() -> sct2e(ctx, sc, -1.0), :valueoutofrange)
        @test throws_kind(() -> sct2e(ctx, sc, mx + 1.0), :valueoutofrange)
        @test_throws SPICE.SpiceError SPICE.sct2e(sc, -1.0)
        @test_throws SPICE.SpiceError SPICE.sct2e(sc, mx + 1.0)
        @test throws_kind(() -> sce2c(ctx, sc, -10.0), :valueoutofrange)
        @test_throws SPICE.SpiceError SPICE.sce2c(sc, -10.0)
        @test throws_kind(() -> sce2c(ctx, sc, 1e12), :valueoutofrange)
        @test throws_kind(() -> sce2t(ctx, sc, 1e12), :valueoutofrange)
        # `sce2t` extends the first record to earlier epochs; the result is then negative.
        @test throws_kind(() -> sce2t(ctx, sc, -10.0), :valueoutofrange)
        @test throws_kind(() -> scdecd(ctx, sc, -1.0), :valueoutofrange)
        @test throws_kind(() -> scdecd(ctx, sc, mx + 1.0), :valueoutofrange)
        @test throws_kind(() -> scfmt(ctx, sc, -1.0), :valueoutofrange)
        @test throws_kind(() -> sctiks(ctx, sc, "  "), :invalidsclkstring)
        @test throws_kind(() -> sctiks(ctx, sc, "1.2.3"), :invalidsclkstring)
        @test throws_kind(() -> sctiks(ctx, sc, "12.x"), :invalidsclkstring)
        @test throws_kind(() -> sctiks(ctx, CKFix.SC_TDT, "0:0:0"), :invalidsclkstring)
        @test throws_kind(() -> scencd(ctx, sc, "/100.5"), :badpartnumber)
        @test throws_kind(() -> scencd(ctx, sc, "x/100.5"), :badpartnumber)
        @test throws_kind(() -> scencd(ctx, sc, "3/100.5"), :badpartnumber)
        @test throws_kind(() -> scencd(ctx, sc, "1/"), :invalidsclkstring)
        # A character of several bytes in a field is an error of the string, not of the code.
        @test all(s -> try scencd(ctx, sc, s); false catch e; e isa Allspice.SpiceError end,
            ("1/é:2", "1/100:4é", "é"))
        @test throws_kind(() -> scencd(ctx, sc, "1/10001.0"), :notinpart)
        @test throws_kind(() -> scencd(ctx, CKFix.SC_TDT, "2/1:0:0"), :notinpart)
        @test throws_kind(() -> scencd(ctx, CKFix.SC_TDT, "900000001:0:0"), :nopartition)
        @test throws_kind(() -> sct2e(ctx, -5, 0.0), :kernelvarnotfound)
        @test throws_kind(() -> sctype(ctx, -5), :kernelvarnotfound)
        @test throws_kind(() -> scpart(ctx, -5), :kernelvarnotfound)
        # The TDT clock needs the leap-second kernel; the TDB clock does not.
        bare = neo_context([CKFix.sclk_kernel()])
        @test throws_kind(() -> sct2e(bare, CKFix.SC_TDT, 1000.0), :nolsk)
        @test sct2e(bare, sc, 1000.0) == sct2e(ctx, sc, 1000.0)
    end

    @apart @testset "derivatives pass through the conversions" begin
        for (sc, _, _, coef) in clocks
            tikmsc = sc == CKFix.SC_TDB ? 256.0 : CKFix.TDT_TIKMSC
            for k in 1:size(coef, 1)-1
                t = coef[k, 1] + (coef[k+1, 1] - coef[k, 1]) * (0.2 + 0.6 * rand(rng))
                x = sct2e(ctx, sc, Allspice.taylor_var(t, Val(2)))
                e0, e1, e2 = Allspice.derivatives(x)
                @test e0 == sct2e(ctx, sc, t)
                h = 64.0
                fd = (sct2e(ctx, sc, t + h) - sct2e(ctx, sc, t - h)) / 2h
                @test isapprox(e1, fd; rtol=1e-9)
                # seconds per tick: the record rate (per count) over the ticks in a count
                @test isapprox(e1, coef[k, 3] / tikmsc; rtol=sc == CKFix.SC_TDB ? 1e-15 : 1e-7)
                sc == CKFix.SC_TDB && @test e2 == 0.0
                y = sce2c(ctx, sc, Allspice.taylor_var(e0, Val(2)))
                c0, c1, _ = Allspice.derivatives(y)
                @test isapprox(c0, t; atol=1e-6)
                @test isapprox(c1 * e1, 1.0; rtol=1e-12)
            end
        end
    end

    @apart @testset "no allocation, inferred types" begin
        for (sc, _, _, coef) in clocks
            t = coef[2, 1] * 0.7
            et = sct2e(ks, sc, t)
            tv = Allspice.taylor_var(t, Val(2))
            ev = Allspice.taylor_var(et, Val(2))
            @test CKFix.warm(alloc_sct2e, ks, sc, t) == 0
            @test CKFix.warm(alloc_sce2c, ks, sc, et) == 0
            @test CKFix.warm(alloc_sce2t, ks, sc, et) == 0
            @test CKFix.warm(alloc_sct2e, ks, sc, tv) == 0
            @test CKFix.warm(alloc_sce2c, ks, sc, ev) == 0
            @test (@inferred sct2e(ks, sc, t)) isa Float64
            @test (@inferred sce2c(ks, sc, et)) isa Float64
            @test (@inferred sce2t(ks, sc, et)) isa Float64
            @test (@inferred sct2e(ks, sc, tv)) isa Allspice.Taylor{2,Float64}
            @test (@inferred sce2c(ks, sc, ev)) isa Allspice.Taylor{2,Float64}
            @test (@inferred sct2e(ctx, sc, 1000)) isa Float64
        end
    end

    @apart @testset "forms on the default context" begin
        kclear()
        foreach(furnsh, set)
        sc = CKFix.SC_TDT
        @test sct2e(sc, 5000.0) == sct2e(ctx, sc, 5000.0)
        @test sce2c(sc, 100.0) == sce2c(ctx, sc, 100.0)
        @test sce2t(sc, 100.0) == sce2t(ctx, sc, 100.0)
        @test sce2s(sc, 100.0) == sce2s(ctx, sc, 100.0)
        @test scs2e(sc, "1/7:30:00") == scs2e(ctx, sc, "1/7:30:00")
        @test scencd(sc, "1/7:30:00") == scencd(ctx, sc, "1/7:30:00")
        @test scdecd(sc, 5000.0) == scdecd(ctx, sc, 5000.0)
        @test sctiks(sc, "7:30:00") == sctiks(ctx, sc, "7:30:00")
        @test scfmt(sc, 5000.0) == scfmt(ctx, sc, 5000.0)
        @test scpart(sc) == scpart(ctx, sc)
        @test sctype(sc) == 1
        kclear()
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
