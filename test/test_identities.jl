# Group 4: algebraic identities at random epochs. No oracle is needed.
#
# Tolerances are relative to the largest vector of the identity: the two sides sum the
# same segment values in a different order, so a few units in the last place of the
# largest term is the expected difference.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

const IDENT_RTOL = 1e-14
const IDENT_ROT_ATOL = 1e-14

group("4. algebraic identities", Kernels.FULL_SET_FILES...) do
    rng = MersenneTwister(SEED + 4)
    ctx = neo_context(Kernels.full_set())
    bodies = [0, 1, 2, 3, 4, 5, 10, 199, 299, 301, 399]
    epochs = random_epochs(rng, 25, DE440S_COVERAGE...)
    state(a, b, et, ref) = spkezr(ctx, a, et, ref, "NONE", b)[1]

    @apart @testset "state(A,B) = -state(B,A), frame $ref" for ref in ("J2000", "ECLIPJ2000",
            "IAU_EARTH", "MOON_PA")
        worst = 0.0
        for et in epochs, _ in 1:8
            a, b = rand(rng, bodies), rand(rng, bodies)
            ab, ba = state(a, b, et, ref), state(b, a, et, ref)
            scale = max(norm(ab[SOneTo(3)]), POS_FLOOR)
            vscale = max(norm(ab[SVector(4, 5, 6)]), VEL_FLOOR)
            worst = max(worst, norm(ab[SOneTo(3)] + ba[SOneTo(3)]) / scale,
                norm(ab[SVector(4, 5, 6)] + ba[SVector(4, 5, 6)]) / vscale)
        end
        measured(@sprintf("  %-58s max=%.3e", "antisymmetry, $ref, relative", worst))
        @test worst <= IDENT_RTOL
    end

    @apart @testset "state(A,C) = state(A,B) + state(B,C), frame $ref" for ref in ("J2000",
            "ECLIPJ2000", "IAU_EARTH", "MOON_PA")
        worst = 0.0
        for et in epochs, _ in 1:8
            a, b, c = rand(rng, bodies), rand(rng, bodies), rand(rng, bodies)
            ac, ab, bc = state(a, c, et, ref), state(a, b, et, ref), state(b, c, et, ref)
            scale = max(norm(ac[SOneTo(3)]), norm(ab[SOneTo(3)]), norm(bc[SOneTo(3)]), POS_FLOOR)
            vscale = max(norm(ac[SVector(4, 5, 6)]), norm(ab[SVector(4, 5, 6)]),
                norm(bc[SVector(4, 5, 6)]), VEL_FLOOR)
            res = ac - (ab + bc)
            worst = max(worst, norm(res[SOneTo(3)]) / scale, norm(res[SVector(4, 5, 6)]) / vscale)
        end
        measured(@sprintf("  %-58s max=%.3e", "triangle identity, $ref, relative", worst))
        @test worst <= IDENT_RTOL
    end

    @apart @testset "the same identities for every derivative order" begin
        worst = 0.0
        for et in epochs[1:10], _ in 1:4
            a, b, c = rand(rng, bodies), rand(rng, bodies), rand(rng, bodies)
            ac = spkder(ctx, a, et, "J2000", c, Val(3))
            ab = spkder(ctx, a, et, "J2000", b, Val(3))
            bc = spkder(ctx, b, et, "J2000", c, Val(3))
            for k in 1:4
                scale = max(norm(ac[k]), norm(ab[k]), norm(bc[k]))
                scale == 0 && continue
                worst = max(worst, norm(ac[k] - (ab[k] + bc[k])) / scale)
            end
        end
        measured(@sprintf("  %-58s max=%.3e", "triangle identity, orders 0..3, relative", worst))
        # Orders 2 and 3 amplify the rounding of the Chebyshev coefficients.
        @test worst <= 1e-9
    end

    frames = ["J2000", "ECLIPJ2000", "B1950", "GALACTIC", "IAU_EARTH", "IAU_MOON", "IAU_MARS",
        "IAU_JUPITER", "MOON_PA", "MOON_ME"]
    @apart @testset "rotation identities" begin
        wchain = 0.0; wtrans = 0.0; wself = 0.0
        for et in epochs, _ in 1:12
            a, b, c = rand(rng, frames), rand(rng, frames), rand(rng, frames)
            # A vector goes a -> b -> c: R(b->c) * R(a->b) = R(a->c). In the notation of
            # the task, with pxform(x, y) the map from x to y read right to left.
            wchain = max(wchain, maximum(abs.(pxform(ctx, b, c, et) * pxform(ctx, a, b, et) -
                                               pxform(ctx, a, c, et))))
            wtrans = max(wtrans, maximum(abs.(pxform(ctx, a, b, et)' - pxform(ctx, b, a, et))))
            wself = max(wself, maximum(abs.(pxform(ctx, a, a, et) - I)))
        end
        measured(@sprintf("  %-58s max=%.3e", "R(b,c) R(a,b) - R(a,c)", wchain))
        measured(@sprintf("  %-58s max=%.3e", "R(a,b)' - R(b,a)", wtrans))
        @test wchain <= IDENT_ROT_ATOL
        @test wtrans <= IDENT_ROT_ATOL
        @test wself <= IDENT_ROT_ATOL
    end

    @apart @testset "a state in a rotated frame is sxform times the J2000 state" begin
        worst = 0.0
        for et in epochs, ref in ("ECLIPJ2000", "IAU_EARTH", "MOON_PA")
            s = state(301, 399, et, "J2000")
            sr = state(301, 399, et, ref)
            res = sxform(ctx, "J2000", ref, et) * s - sr
            worst = max(worst, norm(res[SOneTo(3)]) / norm(s[SOneTo(3)]),
                norm(res[SVector(4, 5, 6)]) / norm(sr[SVector(4, 5, 6)]))
        end
        measured(@sprintf("  %-58s max=%.3e", "sxform * state(J2000) - state(ref), relative", worst))
        @test worst <= 1e-13
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
