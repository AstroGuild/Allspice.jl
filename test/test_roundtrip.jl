# Group: files written by Allspice and read back by Allspice.
#
# The writers and the evaluators of the same data type were built apart, so this group
# closes the loop: a type 13 SPK segment and a type 3 CK segment (the two that a
# trajectory export needs) are written with Allspice, loaded with Allspice and with
# CSPICE, and both engines must return the same states and the same pointing.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

group("roundtrip", Kernels.LSK, Kernels.DE440S) do
    rng = MersenneTwister(SEED + 40)
    dir = mktempdir()
    # --- SPK type 13: a spacecraft near the Earth, states on a smooth curve -------------
    n = 60
    epochs = collect(range(1.0e8, 1.0e8 + 5900.0, length=n))
    a, w = 7000.0, 1.1e-3
    states = [[a * cos(w * (t - epochs[1])), a * sin(w * (t - epochs[1])), 30.0 * sin(2w * (t - epochs[1])),
        -a * w * sin(w * (t - epochs[1])), a * w * cos(w * (t - epochs[1])), 60.0 * w * cos(2w * (t - epochs[1]))]
        for t in epochs]
    bsp = joinpath(dir, "roundtrip.bsp")
    h = spkopn(bsp, "roundtrip", 0)
    spkw13(h, -901, 399, "J2000", epochs[1], epochs[end], "type 13 round trip", 7, states, epochs)
    spkcls(h)
    @test getfat(bsp) == ("DAF", "SPK")
    @test spkobj(bsp) == [-901]

    # --- CK type 3 with a type 1 clock --------------------------------------------------
    tsc = joinpath(dir, "roundtrip.tsc")
    write(tsc, """
        \\begindata
        SCLK_KERNEL_ID            = ( @2026-10-11 )
        SCLK_DATA_TYPE_901        = ( 1 )
        SCLK01_TIME_SYSTEM_901    = ( 1 )
        SCLK01_N_FIELDS_901       = ( 2 )
        SCLK01_MODULI_901         = ( 4294967296 1000 )
        SCLK01_OFFSETS_901        = ( 0 0 )
        SCLK01_OUTPUT_DELIM_901   = ( 1 )
        SCLK_PARTITION_START_901  = ( 0.0 )
        SCLK_PARTITION_END_901    = ( 4.294967296E+12 )
        SCLK01_COEFFICIENTS_901   = ( 0.0  0.0  1.0 )
        \\begintext
        """)
    tf = joinpath(dir, "roundtrip.tf")
    write(tf, """
        \\begindata
        FRAME_RT_BUS              = -901000
        FRAME_-901000_NAME        = 'RT_BUS'
        FRAME_-901000_CLASS       = 3
        FRAME_-901000_CLASS_ID    = -901000
        FRAME_-901000_CENTER      = -901
        CK_-901000_SCLK           = -901
        CK_-901000_SPK            = -901
        \\begintext
        """)
    ticks = [1000.0 * t for t in epochs]          # 1000 ticks per second, ET = seconds
    quats = map(epochs) do t
        ang = 2.0e-3 * (t - epochs[1])
        axis = normalize(SVector(0.3, -0.5, 0.8))
        [cos(ang / 2), (sin(ang / 2) * axis)...]
    end
    avvs = [collect(2.0e-3 * normalize(SVector(0.3, -0.5, 0.8))) for _ in epochs]
    bc = joinpath(dir, "roundtrip.bc")
    h = ckopn(bc, "roundtrip", 0)
    ckw03(h, ticks[1], ticks[end], -901000, "J2000", "type 3 round trip", ticks, quats,
        avvs, [ticks[1]])
    ckcls(h)
    @test getfat(bc) == ("DAF", "CK")

    set = (Kernels.LSK, Kernels.DE440S, bsp, tsc, tf, bc)
    oracle_load(set)
    ctx = neo_context(set)
    ets = vcat(epochs[1], epochs[end], epochs[7], random_epochs(rng, nsample(30, 200), epochs[1], epochs[end]))

    @apart @testset "type 13 SPK written and read by Allspice" begin
        worst = 0.0
        for et in ets, (obs, ref) in ((399, "J2000"), (301, "J2000"), (10, "ECLIPJ2000"))
            s, lt = spkezr(ctx, -901, et, ref, "NONE", obs)
            sc, ltc = SPICE.spkezr("-901", et, ref, "NONE", string(obs))
            worst = max(worst, norm(s - sc) / norm(sc))
        end
        measured(@sprintf("  %-58s max=%.3e", "round trip type 13, state against CSPICE, relative", worst))
        @test worst <= 1e-14
        # The states at the knots are the states that were written.
        for k in (1, 7, n)
            @test isapprox(spkezr(ctx, -901, epochs[k], "J2000", "NONE", 399)[1], states[k]; rtol=1e-14)
        end
        @test spkcov(bsp, -901)[1] == epochs[1]
    end

    @apart @testset "type 3 CK written and read by Allspice" begin
        worst = 0.0
        allfound = true
        for et in ets
            R = pxform(ctx, "J2000", "RT_BUS", et)
            worst = max(worst, maximum(abs.(R - SPICE.pxform("J2000", "RT_BUS", et))))
            allfound &= ckgpav(ctx, -901000, 1000.0 * et, 0.0, "J2000") !== nothing
        end
        @test allfound
        measured(@sprintf("  %-58s max=%.3e", "round trip type 3 CK, pxform against CSPICE, absolute", worst))
        @test worst <= 1e-14
        s = spkezr(ctx, 301, ets[5], "RT_BUS", "LT+S", -901)[1]
        sc = SPICE.spkezr("301", ets[5], "RT_BUS", "LT+S", "-901")[1]
        @test isapprox(s, sc; rtol=1e-10)
    end
    # CSPICE keeps a loaded binary kernel open, and a context keeps it mapped until the
    # garbage collector frees the snapshot; Windows deletes neither kind of file.
    SPICE.kclear()
    kclear(ctx)
    remove_temp(dir)
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
