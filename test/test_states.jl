# Group 1: states against CSPICE.
#
# Feature file: features/load_then_query.feature
#
# Error measure, in the manner of the ANISE validation: for each query the absolute
# difference per component, |Allspice - CSPICE|, and the same difference divided by the
# Euclidean norm of the CSPICE vector (position norm for position components, velocity
# norm for velocity components). The assertions are on the quantiles q75 and q99 of the
# relative difference and on the maximum absolute difference.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

# `StateTol`, the tolerances `TOL_*` and the floors are in helpers.jl, because other
# groups use them.

# `compare_states` and `assert_states` are in helpers.jl as well.

"All ordered pairs (a, b), a != b, of `bodies`."
ordered_pairs(bodies) = [(a, b) for a in bodies for b in bodies if a != b]

group("1. states against CSPICE", Kernels.FULL_SET_FILES...) do
    rng = MersenneTwister(SEED)
    set = Kernels.full_set()
    oracle_load(set)
    ctx = neo_context(set)

    bodies = spk_bodies(Kernels.DE440S)
    @test length(bodies) >= 14
    pairs = ordered_pairs(bodies)

    @apart @testset "de440s, all ordered pairs, J2000" begin
        epochs = random_epochs(rng, 200, DE440S_COVERAGE...)
        println("  bodies: ", bodies, "; ", length(pairs), " ordered pairs x ",
            length(epochs), " epochs")
        d = compare_states(ctx, pairs, epochs, "J2000")
        assert_states("de440s J2000", d, TOL_J2000)
    end

    # A smaller set in the other frames: a mix of short and long chains.
    small = [(301, 399), (399, 301), (10, 399), (399, 10), (301, 10), (5, 399),
        (4, 301), (199, 299), (3, 0), (0, 301)]
    small = [p for p in small if (p[1] == 0 || p[1] in bodies) && (p[2] == 0 || p[2] in bodies)]
    epochs50 = random_epochs(rng, 50, DE440S_COVERAGE...)
    @apart @testset "de440s, frame $ref" for (ref, tol) in (("ECLIPJ2000", TOL_INERTIAL),
            ("IAU_EARTH", TOL_BODYFIXED), ("IAU_MOON", TOL_BODYFIXED),
            ("MOON_PA", TOL_BODYFIXED), ("MOON_ME", TOL_BODYFIXED))
        d = compare_states(ctx, small, epochs50, ref)
        assert_states("de440s $ref", d, tol)
    end

    @apart @testset "forms of the call" begin
        epochs = random_epochs(rng, 20, DE440S_COVERAGE...)
        named = [("MOON", 301, "EARTH", 399), ("EARTH", 399, "SUN", 10),
            ("SUN", 10, "MOON", 301), ("EARTH BARYCENTER", 3, "SOLAR SYSTEM BARYCENTER", 0),
            ("JUPITER BARYCENTER", 5, "VENUS", 299), ("MERCURY", 199, "EMB", 3)]
        for (tn, ti, on, oi) in named, et in epochs
            byid = spkezr(ctx, ti, et, "J2000", "NONE", oi)
            # Name form, integer-ID form, integer-string form, integer frame ID and a
            # name with other case and blanks: one and the same answer, bit for bit.
            @test spkezr(ctx, tn, et, "J2000", "NONE", on) == byid
            @test spkezr(ctx, string(ti), et, "J2000", "NONE", string(oi)) == byid
            @test spkezr(ctx, ti, et, 1, "NONE", oi) == byid
            @test spkezr(ctx, "  " * lowercase(tn) * " ", et, "j2000", "NONE", on) == byid
            # spkpos is the position part with the same light time.
            pos, lt = spkpos(ctx, tn, et, "J2000", "NONE", on)
            @test pos == byid[1][SOneTo(3)]
            @test lt == byid[2]
            @test spkpos(ctx, ti, et, 1, "NONE", oi) == (pos, lt)
            # spkgeo is the same geometric state.
            @test spkgeo(ctx, ti, et, "J2000", oi) == byid
            # CSPICE answers for the name form and for spkpos and spkgeo.
            cstate, clt = SPICE.spkezr(tn, et, "J2000", "NONE", on)
            @test isapprox(byid[1], cstate; rtol=1e-13)
            @test isapprox(lt, clt; rtol=1e-13)
            cpos, cplt = SPICE.spkpos(tn, et, "J2000", "NONE", on)
            @test isapprox(pos, cpos; rtol=1e-13)
            cgeo, cglt = SPICE.spkgeo(ti, et, "J2000", oi)
            @test isapprox(spkgeo(ctx, ti, et, "J2000", oi)[1], cgeo; rtol=1e-13)
            @test isapprox(lt, cglt; rtol=1e-13)
            # Light time is the one-way time |r|/c.
            @test isapprox(lt, norm(pos) / 299792.458; rtol=4eps())
        end
        # The result types are the ones the contract gives.
        state, lt = spkezr(ctx, 301, 0.0, "J2000", "NONE", 399)
        @test state isa SVector{6,Float64}
        @test lt isa Float64
        @test spkpos(ctx, 301, 0.0, "J2000", "NONE", 399)[1] isa SVector{3,Float64}
        # A body relative to itself is the zero state.
        @test all(iszero, spkezr(ctx, 399, 0.0, "J2000", "NONE", 399)[1])
    end

    @apart @testset "default context" begin
        kclear()
        foreach(furnsh, set)
        for et in random_epochs(rng, 10, DE440S_COVERAGE...)
            @test spkezr("MOON", et, "J2000", "NONE", "EARTH") ==
                  spkezr(ctx, "MOON", et, "J2000", "NONE", "EARTH")
            @test spkpos(301, et, "IAU_EARTH", "NONE", 399) ==
                  spkpos(ctx, 301, et, "IAU_EARTH", "NONE", 399)
        end
        kclear()
    end
end

group("1b. states against CSPICE, de405", Kernels.DE405, Kernels.LSK) do
    rng = MersenneTwister(SEED + 1)
    set = [Kernels.LSK, Kernels.DE405]
    oracle_load(set)
    ctx = neo_context(set)
    # Body 1000000 of this file has the degenerate coverage [0, 0]; leave it out.
    bodies = [b for b in spk_bodies(Kernels.DE405) if b < 1000]
    @test length(bodies) >= 14
    epochs = random_epochs(rng, 20, DE405_COVERAGE...)
    println("  bodies: ", bodies, "; ", length(epochs), " epochs")
    d = compare_states(ctx, ordered_pairs(bodies), epochs, "J2000")
    assert_states("de405 J2000", d, TOL_J2000)
    d = compare_states(ctx, [(301, 399), (399, 10), (499, 301)], epochs, "ECLIPJ2000")
    assert_states("de405 ECLIPJ2000", d, TOL_INERTIAL)
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
