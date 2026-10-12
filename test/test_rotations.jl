# Group 2: rotations against CSPICE.
#
# Error measure: absolute difference per matrix element. Elements of a rotation matrix
# are at most 1 in magnitude, so the absolute difference is also the relative one. For
# the derivative block of `sxform` the difference is divided by `RATE_SCALE`, an upper
# bound of the angular rates involved, so both blocks share one threshold.

const ROT_FRAMES = ["J2000", "ECLIPJ2000", "B1950", "GALACTIC", "IAU_EARTH", "IAU_MOON",
    "IAU_MARS", "IAU_JUPITER", "MOON_PA", "MOON_ME"]
const ROT_INERTIAL = ("J2000", "ECLIPJ2000", "B1950", "GALACTIC")

# Inertial to inertial: constant matrices, a few roundings of difference at most.
const ROT_TOL_INERTIAL = (q99=2e-15, max=5e-15)
# With a body-fixed frame: see the note on TOL_BODYFIXED in test_states.jl; the angle W of
# Jupiter reaches 8e5 rad at the ends of the sampled range (unit in the last place 1e-10).
const ROT_TOL_BODYFIXED = (q99=1e-9, max=5e-9)
# Largest angular rate among the frames, rad/s (Jupiter 1.76e-4, pairs can add).
const RATE_SCALE = 4e-4
const ORTHO_TOL = 1e-14

group("2. rotations against CSPICE", Kernels.TEXT_SET_FILES..., Kernels.MOON_BPC) do
    rng = MersenneTwister(SEED + 2)
    set = [Kernels.text_set()..., Kernels.MOON_BPC]
    oracle_load(set)
    ctx = neo_context(set)
    epochs = random_epochs(rng, 40, DE440S_COVERAGE...)

    din = Float64[]; dbf = Float64[]          # pxform differences
    sin_ = Float64[]; sbf = Float64[]         # sxform derivative-block differences, scaled
    ortho = Float64[]
    structure_ok = true
    for from in ROT_FRAMES, to in ROT_FRAMES
        from == to && continue
        inertial = from in ROT_INERTIAL && to in ROT_INERTIAL
        for et in epochs
            R = pxform(ctx, from, to, et)
            X = sxform(ctx, from, to, et)
            Rc = SPICE.pxform(from, to, et)
            Xc = SPICE.sxform(from, to, et)
            append!(inertial ? din : dbf, abs.(R .- Rc))
            append!(inertial ? din : dbf, abs.(X[SOneTo(3), SOneTo(3)] .- Xc[1:3, 1:3]))
            append!(inertial ? sin_ : sbf,
                abs.(X[SVector(4, 5, 6), SOneTo(3)] .- Xc[4:6, 1:3]) ./ RATE_SCALE)
            push!(ortho, maximum(abs.(R * R' - I)))
            push!(ortho, abs(det(R) - 1))
            # sxform = [R 0; dR R], with R the pxform matrix to within rounding (the
            # two come from the order 1 and the order 0 evaluation; measured 1.1e-16).
            structure_ok &= maximum(abs.(X[SOneTo(3), SOneTo(3)] - R)) <= 4eps()
            structure_ok &= X[SVector(4, 5, 6), SVector(4, 5, 6)] == X[SOneTo(3), SOneTo(3)]
            structure_ok &= all(iszero, X[SOneTo(3), SVector(4, 5, 6)])
            inertial && (structure_ok &= all(iszero, X[SVector(4, 5, 6), SOneTo(3)]))
        end
    end
    a = report_quantiles("pxform, inertial pairs, absolute", din)
    b = report_quantiles("pxform, body-fixed pairs, absolute", dbf)
    c = report_quantiles("sxform dR block / $(RATE_SCALE), body-fixed pairs", sbf)
    o = report_quantiles("orthogonality |R R' - I|, |det R - 1|", ortho)
    @test a.q99 <= ROT_TOL_INERTIAL.q99
    @test a.max <= ROT_TOL_INERTIAL.max
    @test b.q99 <= ROT_TOL_BODYFIXED.q99
    @test b.max <= ROT_TOL_BODYFIXED.max
    @test c.q99 <= ROT_TOL_BODYFIXED.q99
    @test c.max <= ROT_TOL_BODYFIXED.max
    @test maximum(sin_) == 0.0
    @test o.max <= ORTHO_TOL
    @test structure_ok

    @testset "forms and types" begin
        et = epochs[1]
        R = pxform(ctx, "J2000", "IAU_EARTH", et)
        @test R isa SMatrix{3,3,Float64,9}
        @test sxform(ctx, "J2000", "IAU_EARTH", et) isa SMatrix{6,6,Float64,36}
        @test pxform(ctx, 1, SPICE.namfrm("IAU_EARTH"), et) == R
        @test pxform(ctx, "j2000", " iau_earth ", et) == R
        @test pxform(ctx, "IAU_MARS", "IAU_MARS", et) == SMatrix{3,3,Float64}(I)
        # The alias frames give the DE440 frames exactly as CSPICE resolves them.
        @test isapprox(pxform(ctx, "MOON_PA", "MOON_PA_DE440", et), I; atol=1e-15)
        @test isapprox(pxform(ctx, "MOON_ME", "MOON_ME_DE440_ME421", et), I; atol=1e-15)
        @test isapprox(pxform(ctx, "MOON_PA", "MOON_ME", et),
            SPICE.pxform("MOON_PA", "MOON_ME", et); atol=1e-15)
    end

    # The same frames near J2000, where the angle W is small: this isolates the model
    # (pole, nutation-precession terms, Chebyshev angles) from the rounding of W.
    @testset "near J2000" begin
        near = random_epochs(rng, 40, -86400.0 * 30, 86400.0 * 30)
        dn = Float64[]
        for from in ROT_FRAMES, to in ROT_FRAMES, et in near
            from == to && continue
            append!(dn, abs.(pxform(ctx, from, to, et) .- SPICE.pxform(from, to, et)))
        end
        n = report_quantiles("pxform, all pairs within 30 d of J2000", dn)
        @test n.q99 <= 1e-12
        @test n.max <= 1e-11
    end
end
