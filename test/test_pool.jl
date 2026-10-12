# Group 5: kernel pool, names and time conversions against CSPICE.
#
# Conventions for an absent item, as the prototype source documents them: `gdpool`,
# `gcpool`, `gipool` and `bodvrd` throw `SpiceError(:notfound)` (SPICE.jl returns
# `nothing`, or throws for `bodvrd`); `bodn2c`, `bods2c` and `bodc2n` return `nothing`;
# `namfrm` returns 0 and `frmnam` returns "", as CSPICE does.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of NeoSpice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`NEOSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

"Pool variable names of CSPICE that match `template`."
function oracle_names(template::AbstractString)
    out = String[]
    start = 1
    while true
        part = SPICE.gnpool(template, start, 500)
        part === nothing && break
        append!(out, part)
        length(part) < 500 && break
        start += 500
    end
    return unique!(out)
end

"Number of time conversions of NeoSpice that threw, and the first messages."
const NEO_TIME_ERRORS = Ref(0)

"""
    neo_time(f, what) -> Float64

Value of `f()`, or `NaN` if it throws a `SpiceError`. The error is counted and the first
three are printed, so that one defect does not hide the other measurements of the group.
"""
function neo_time(f, what)
    try
        return Float64(f())
    catch e
        e isa NeoSpice.SpiceError || rethrow()
        NEO_TIME_ERRORS[] += 1
        NEO_TIME_ERRORS[] <= 3 && println("  NeoSpice threw for ", repr(what), ": ", sprint(showerror, e))
        return NaN
    end
end

"UTC calendar string `YYYY-MM-DDThh:mm:ss.fff` of the given fields (no validation of the day)."
iso(y, mo, d, h, mi, s) = @sprintf("%04d-%02d-%02dT%02d:%02d:%06.3f", y, mo, d, h, mi, s)

group("5. kernel pool and names against CSPICE", Kernels.TEXT_SET_FILES...,
    Kernels.MOON_TF_OLD) do
    rng = MersenneTwister(SEED + 5)
    set = [Kernels.text_set()..., Kernels.MOON_TF_OLD]
    have_gm = isfile(Kernels.GM_TPC)
    have_gm ? push!(set, Kernels.GM_TPC) :
    skip_notice("5. GM values", "missing kernel: " * Kernels.GM_TPC)
    oracle_load(set)
    ctx = neo_context(set)

    @apart @testset "every pool variable of CSPICE" begin
        names = oracle_names("*")
        println("  ", length(names), " pool variables in CSPICE")
        @test length(names) > 500
        nnum = 0; nchr = 0; bad = String[]
        for name in names
            n, type = SPICE.dtpool(name)
            if type === :N
                ref = SPICE.gdpool(name; room=n)
                got = try gdpool(ctx, name) catch e; e end
                (got isa Vector{Float64} && got == ref) || push!(bad, name)
                nnum += 1
            else
                ref = SPICE.gcpool(name; room=n, lenout=1024)
                got = try gcpool(ctx, name) catch e; e end
                (got isa Vector{String} && got == ref) || push!(bad, name)
                nchr += 1
            end
        end
        println("  numeric: ", nnum, ", character: ", nchr, ", different: ", length(bad))
        isempty(bad) || println("  first different variables: ", first(bad, 10))
        @test isempty(bad)
        @test nnum > 0 && nchr > 0
    end

    @apart @testset "gdpool, gcpool, gipool samples" begin
        @test gdpool(ctx, "BODY399_RADII") == [6378.1366, 6378.1366, 6356.7519]
        @test gdpool(ctx, "BODY399_RADII") isa Vector{Float64}
        for name in ("BODY301_PM", "BODY499_POLE_RA", "BODY5_NUT_PREC_ANGLES",
                "DELTET/DELTA_AT", "DELTET/DELTA_T_A", "TKFRAME_31009_ANGLES",
                "BODY399_N_GEOMAG_CTR_DIPOLE_LAT")
            @test gdpool(ctx, name) == SPICE.gdpool(name; room=1000)
        end
        for name in ("FRAME_31008_NAME", "TKFRAME_31009_SPEC", "TKFRAME_31009_UNITS",
                "TKFRAME_31009_RELATIVE", "TKFRAME_31000_RELATIVE", "FRAME_31006_NAME")
            @test gcpool(ctx, name) == SPICE.gcpool(name)
        end
        for name in ("FRAME_MOON_PA_DE440", "TKFRAME_31009_AXES", "FRAME_31008_CLASS",
                "BODY399_RADII")
            @test gipool(ctx, name) == SPICE.gipool(name)
        end
        # The result is a copy: changing it does not change the pool.
        v = gdpool(ctx, "BODY399_RADII"); v[1] = -1.0
        @test gdpool(ctx, "BODY399_RADII")[1] == 6378.1366
        # Absent variable and wrong type.
        @test SPICE.gdpool("NO_SUCH_VARIABLE") === nothing
        @test throws_kind(() -> gdpool(ctx, "NO_SUCH_VARIABLE"), :notfound)
        @test throws_kind(() -> gcpool(ctx, "NO_SUCH_VARIABLE"), :notfound)
        @test throws_kind(() -> gipool(ctx, "NO_SUCH_VARIABLE"), :notfound)
        @test throws_kind(() -> gcpool(ctx, "BODY399_RADII"), :badtype)
        @test throws_kind(() -> gdpool(ctx, "TKFRAME_31009_SPEC"), :badtype)
    end

    @apart @testset "bodvrd" begin
        for body in ("EARTH", "MOON", "SUN", "MARS", "JUPITER", "399", "earth")
            @test bodvrd(ctx, body, "RADII") == SPICE.bodvrd(body, "RADII")
        end
        @test bodvrd(ctx, 399, "RADII") == SPICE.bodvrd("399", "RADII")
        for (body, item) in (("EARTH", "POLE_RA"), ("MOON", "NUT_PREC_PM"), ("MARS", "PM"),
                ("JUPITER", "NUT_PREC_RA"), ("SATURN", "POLE_DEC"))
            @test bodvrd(ctx, body, item) == SPICE.bodvrd(body, item, 100)
        end
        if have_gm
            for body in ("EARTH", "MOON", "SUN", "EARTH BARYCENTER", "JUPITER BARYCENTER", "3")
                @test bodvrd(ctx, body, "GM") == SPICE.bodvrd(body, "GM")
            end
            @test bodvrd(ctx, "EARTH", "GM") isa Vector{Float64}
        end
        @test_throws SPICE.SpiceError SPICE.bodvrd("EARTH", "NO_SUCH_ITEM")
        # A blank name is an error in CSPICE; NeoSpice reports "not found".
        @test bodn2c(ctx, "") === nothing
        @test namfrm(ctx, "") == 0
        @test throws_kind(() -> bodvrd(ctx, "EARTH", "NO_SUCH_ITEM"), :notfound)
        @test throws_kind(() -> bodvrd(ctx, "NO_SUCH_BODY", "RADII"), :unknownbody)
    end

    @apart @testset "bodn2c, bods2c, bodc2n" begin
        names = ["EARTH", "earth", "  Earth  ", "MOON", "SUN", "SSB", "SOLAR SYSTEM BARYCENTER",
            "EMB", "EARTH BARYCENTER", "Earth-Moon   Barycenter", "MARS", "PHOBOS", "IO",
            "JUPITER BARYCENTER", "PLUTO", "CHARON", "CERES", "CASSINI", "VOYAGER 1", "MRO",
            "JUICE", "DSS-14", "TEMPEL 1", "NO_SUCH_BODY", "399", "-82"]
        for name in names
            @test bodn2c(ctx, name) == SPICE.bodn2c(name)
            @test bods2c(ctx, name) == SPICE.bods2c(name)
        end
        codes = [0, 1, 3, 10, 199, 299, 301, 399, 401, 499, 501, 599, 699, 799, 899, 999, 901,
            2000001, -82, -31, -32, -74, 399014, 1000005, 123456789, -99999]
        for code in codes
            @test bodc2n(ctx, code) == SPICE.bodc2n(code)
        end
        # Round trip for a random sample of the built-in codes.
        for code in rand(rng, 1:999, 200)
            name = SPICE.bodc2n(code)
            @test bodc2n(ctx, code) == name
            name === nothing || @test bodn2c(ctx, name) == SPICE.bodn2c(name)
        end
    end

    @apart @testset "namfrm, frmnam" begin
        names = ["J2000", "j2000", "ECLIPJ2000", "B1950", "FK4", "GALACTIC", "MARSIAU",
            "ECLIPB1950", "DE-118", "IAU_EARTH", "IAU_MOON", "IAU_MARS", "IAU_JUPITER",
            "IAU_SUN", "IAU_PHOBOS", "ITRF93", "EARTH_FIXED", "MOON_PA", "MOON_ME",
            "MOON_PA_DE440", "MOON_ME_DE440_ME421", "MOON_PA_DE421", "MOON_ME_DE421",
            "NO_SUCH_FRAME"]
        for name in names
            @test namfrm(ctx, name) == SPICE.namfrm(name)
        end
        for id in [1:21..., 10001:10020..., 10081, 13000, 31000, 31001, 31006, 31007, 31008,
                31009, 0, -5, 123456]
            @test frmnam(ctx, id) == SPICE.frmnam(id)
        end
    end

    @apart @testset "utc2et, et2utc, deltet" begin
        # Random calendar dates from 1972 to 2035 with millisecond resolution.
        worst = 0.0; same_str2et = true
        for _ in 1:300
            s = iso(rand(rng, 1972:2035), rand(rng, 1:12), rand(rng, 1:28), rand(rng, 0:23),
                rand(rng, 0:59), rand(rng, 0:59999) / 1000)
            ref = SPICE.utc2et(s)
            got = neo_time(() -> utc2et(ctx, s), s)
            worst = max(worst, abs(got - ref))
            isnan(got) || (same_str2et &= str2et(ctx, s) == got)
        end
        @test same_str2et
        measured(@sprintf("  %-58s max=%.3e", "utc2et, random dates, absolute s", worst))
        # One unit in the last place of 1e9 s is 1.2e-7 s.
        @test worst <= 5e-7

        # Around every leap second of the kernel, and inside the leap second itself.
        at = SPICE.gdpool("DELTET/DELTA_AT"; room=200)
        nleap = 0; wleap = 0.0
        for i in 4:2:length(at)                 # the first entry is the initial offset
            # The entry is the first instant with the new offset: 1 Jan or 1 Jul 00:00.
            day = SPICE.et2utc(at[i] + 100.0, "ISOC", 0)
            y = parse(Int, day[1:4]); mo = parse(Int, day[6:7])
            (py, pmo, pd) = mo == 1 ? (y - 1, 12, 31) : (y, 6, 30)
            for s in (iso(py, pmo, pd, 23, 59, 58.25), iso(py, pmo, pd, 23, 59, 59.75),
                    iso(py, pmo, pd, 23, 59, 60.0), iso(py, pmo, pd, 23, 59, 60.5),
                    iso(y, mo, 1, 0, 0, 0.0), iso(y, mo, 1, 0, 0, 1.5))
                wleap = max(wleap, abs(neo_time(() -> utc2et(ctx, s), s) - SPICE.utc2et(s)))
            end
            nleap += 1
        end
        println("  leap seconds checked: ", nleap)
        measured(@sprintf("  %-58s max=%.3e", "utc2et, around leap seconds, absolute s", wleap))
        @test nleap >= 26
        @test wleap <= 5e-7
        # The leap second is one second long: 23:59:60.5 lies between its neighbours.
        u(s) = neo_time(() -> utc2et(ctx, s), s)
        @test u("2016-12-31T23:59:60.500") - u("2016-12-31T23:59:59.500") ≈ 1.0 atol = 1e-6
        @test u("2017-01-01T00:00:00.500") - u("2016-12-31T23:59:60.500") ≈ 1.0 atol = 1e-6
        println("  utc2et calls of NeoSpice that threw: ", NEO_TIME_ERRORS[])
        @test NEO_TIME_ERRORS[] == 0

        # Forms of the input string.
        for s in ("2000-01-01T12:00:00", "2000-01-01", "2024-02-29T23:59:59.999",
                "1999-12-31T23:59:59.123456", "2020-06-15T06:07:08 UTC", "2000-01-01 12:00:00 TDB",
                "2030-03-04 05:06:07.5 TDB")
            ref = try SPICE.str2et(s) catch; nothing end
            ref === nothing && continue
            @test isapprox(neo_time(() -> str2et(ctx, s), s), ref; atol=5e-7)
        end

        # et2utc: equal strings in the three formats.
        ets = random_epochs(rng, 200, SPICE.utc2et("1972-01-02T00:00:00"),
            SPICE.utc2et("2035-01-01T00:00:00"))
        for i in 4:2:length(at)
            t = SPICE.utc2et(SPICE.et2utc(at[i] + 100.0, "ISOC", 0)[1:10] * "T00:00:00")
            append!(ets, (t - 2.5, t - 1.5, t - 0.5, t, t + 0.5))
        end
        nbad = 0
        for et in ets, (fmt, prec) in (("ISOC", 3), ("ISOC", 0), ("ISOC", 6), ("C", 3), ("J", 6))
            ref = SPICE.et2utc(et, fmt, prec)
            got = et2utc(ctx, et, fmt, prec)
            if got != ref
                nbad += 1
                nbad <= 5 && println("  et2utc differs: et=", et, " ", fmt, " ", prec, ": ",
                    repr(got), " against ", repr(ref))
            end
        end
        println("  et2utc strings compared: ", 5 * length(ets), ", different: ", nbad)
        @test nbad == 0
        # Round trip.
        for et in ets[1:50]
            @test abs(neo_time(() -> utc2et(ctx, et2utc(ctx, et, "ISOC", 6)), "round trip") - et) <= 1e-6
        end

        # deltet in both directions.
        wd = 0.0
        for et in ets
            wd = max(wd, abs(deltet(ctx, et, "ET") - SPICE.deltet(et, "ET")),
                abs(deltet(ctx, et, "UTC") - SPICE.deltet(et, "UTC")))
        end
        measured(@sprintf("  %-58s max=%.3e", "deltet, absolute s", wd))
        @test wd <= 1e-9
    end

    @apart @testset "no leap-second kernel" begin
        empty = SpiceContext()
        @test_throws NeoSpice.SpiceError utc2et(empty, "2000-01-01T12:00:00")
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
