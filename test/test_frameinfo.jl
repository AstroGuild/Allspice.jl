# Group 19: frame information, the lists of frames, body constants, body-fixed rotations,
# binary PCK coverage, local solar time and L_s against CSPICE.
#
# Conventions for an absent item: `frinfo`, `cidfrm`, `cnmfrm` and `ccifrm` return
# `nothing` where CSPICE returns `found = false`; `tkfram` throws where `tkfram_c` signals
# an error. `tkfram_` has no wrapper in SPICE.jl and is called through `ccall`.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

isdefined(@__MODULE__, :CKFix) || include("ckfixtures.jl")
isdefined(@__MODULE__, :DynFix) || include("dynfixtures.jl")

"`tkfram_` of CSPICE: `(rot, relative frame ID)`, the error kind, or `nothing` if not found."
function oracle_tkfram(id::Integer)
    rot = zeros(9)
    frame = Ref{Cint}(0)
    found = Ref{Cint}(0)
    ccall((:tkfram_, SPICE.libcspice), Cint, (Ref{Cint}, Ptr{Float64}, Ref{Cint}, Ref{Cint}),
        Cint(id), rot, frame, found)
    kind = cspice_kind()
    kind === :none || return kind
    found[] == 0 && return nothing
    return SMatrix{3,3,Float64,9}(rot), Int(frame[])
end

"`tkfram(ctx, id)`, or the kind of the `SpiceError` that it throws."
neo_tkfram(ctx, id) = try tkfram(ctx, id) catch e; e isa Allspice.SpiceError ? e.kind : rethrow() end

group("19. frame information, body constants, local solar time", Kernels.FULL_SET_FILES...) do
    rng = MersenneTwister(SEED + 19)
    set = [Kernels.full_set()..., CKFix.frame_kernel(), DynFix.frame_kernel()]
    extra = joinpath(mktempdir(), "objects.tf")
    write(extra, """
    KPL/FK

    \\begindata
       OBJECT_MARS_FRAME      = 'IAU_MARS'
       OBJECT_599_FRAME       = 'J2000'
       OBJECT_EARTH_FRAME     = 13000
       OBJECT_-77_FRAME       = 'TESTSC_T3'
       OBJECT_301_FRAME       = 'MOON_ME'
       OBJECT_299_FRAME       = 31000
       NAIF_BODY_NAME        += ( 'TEST_PROBE' )
       NAIF_BODY_CODE        += ( -77 )
    \\begintext
    """)
    # The second pass adds the `OBJECT_*_FRAME` variables. Only `cidfrm` and `cnmfrm` read
    # them, so the other calls are compared in the first pass only.
    for with_objects in (false, true)
        with_objects && push!(set, extra)
        oracle_load(set)
        ctx = neo_context(set)
        label = with_objects ? " (with OBJECT_*_FRAME variables)" : ""

        frames = vcat(collect(-5:25), collect(10001:10150), collect(13000:13010),
            [1850000 + i for i in 0:25], [31000, 31001, 31002, 31003, 31004, 31005, 31006,
                31007, 31008, -77001, -77003, -77900, -77099, -99001, 1400499, 99999999])
        with_objects || @testset "frinfo" begin
            nfound = 0
            for id in frames
                want = SPICE.frinfo(id)
                got = frinfo(ctx, id)
                got == want || @info "frinfo" id got want
                @test got == want
                want === nothing || (nfound += 1)
            end
            println("  frinfo: ", length(frames), " frame IDs, ", nfound, " found")
            @test nfound > 150
            @test (@inferred Union{Nothing,NTuple{3,Int}} frinfo(ctx, 13000)) == SPICE.frinfo(13000)
        end

        @testset "cidfrm and cnmfrm$label" begin
            bodies = vcat(collect(-3:12), collect(199:100:999), collect(301:1:305),
                collect(401:402), collect(501:516), collect(601:618), [10, 9511010, 2000001,
                    2000004, 2000433, 1000005, -77, -78, -31, 3, 4, 399, 301, 123456789])
            nfound = 0
            for b in bodies
                want = SPICE.cidfrm(b)
                got = cidfrm(ctx, b)
                got == want || @info "cidfrm" b got want
                @test got == want
                want === nothing || (nfound += 1)
            end
            println("  cidfrm: ", length(bodies), " bodies, ", nfound, " found")
            @test nfound > 50
            names = ["EARTH", "earth", "  Moon ", "MARS", "JUPITER", "SUN", "PHOBOS", "TEST_PROBE", "VENUS",
                "SOLAR SYSTEM BARYCENTER", "MARS BARYCENTER", "NOT A BODY", "TITAN", "CERES",
                "PLUTO", "IO", "Earth  Barycenter", "399", "-77", "499"]
            for n in names
                want = SPICE.cnmfrm(n)
                got = cnmfrm(ctx, n)
                got == want || @info "cnmfrm" n got want
                @test got == want
            end
        end

        with_objects || @testset "ccifrm" begin
            n = 0
            for id in frames
                info = SPICE.frinfo(id)
                info === nothing && continue
                want = SPICE.ccifrm(info[2], info[3])
                got = ccifrm(ctx, info[2], info[3])
                got == want || @info "ccifrm" id info got want
                @test got == want
                n += 1
            end
            for (c, k) in ((1, 0), (1, 22), (2, 398), (2, 31006), (3, -77001), (3, 5), (4, -77900),
                (5, 1850001), (6, 1850021), (7, 1), (0, 0), (2, 10), (4, 31001))
                want = SPICE.ccifrm(c, k)
                got = ccifrm(ctx, c, k)
                got == want || @info "ccifrm" c k got want
                @test got == want
            end
            println("  ccifrm: ", n, " (class, class ID) pairs of loaded frames")
        end

        with_objects || @testset "tkfram" begin
            worst = 0.0
            nfound = 0
            for id in (-77900, 31000, 31001, 31002, 31003, 31006, 31007, 1850001, 10013, 5, -5)
                want = oracle_tkfram(id)
                got = neo_tkfram(ctx, id)
                want isa Tuple || (@test got === want; continue)
                nfound += 1
                @test got[2] == want[2]
                worst = max(worst, maximum(abs.(got[1] .- want[1])))
            end
            @test nfound >= 3
            @test worst <= 4 * eps()
            measured(@sprintf("  %-58s n=%-7d max=%.3e", "tkfram: element of the matrix",
                nfound, worst))
        end
    end

    @testset "tkfram and pxform: a TK frame definition that is not valid" begin
        # Regression: tkfram returned `nothing`, and pxform threw :framedatanotfound, for
        # an incomplete or inconsistent TK frame definition; tkfram_c and pxform_c signal
        # the error of tkfram_ for each. One case for each check of tkfram_.
        head = ["FRAME_TKA = -1000", "FRAME_-1000_NAME = 'TKA'", "FRAME_-1000_CLASS = 4",
            "FRAME_-1000_CLASS_ID = -1000", "FRAME_-1000_CENTER = 399"]
        rel(name) = "TKFRAME_-1000_RELATIVE = '$name'"
        spec(s) = "TKFRAME_-1000_SPEC = '$s'"
        angles(units, axes) = [rel("J2000"); spec("ANGLES"); "TKFRAME_-1000_ANGLES = (1, 2, 3)";
            "TKFRAME_-1000_AXES = $axes"; "TKFRAME_-1000_UNITS = '$units'"]
        unit = "TKFRAME_-1000_MATRIX = (1, 0, 0, 0, 1, 0, 0, 0, 1)"
        cases = (
            (-1000, head),                                                  # variablenotfound
            (-1000, [head; rel("NO_SUCH_FRAME"); spec("MATRIX"); unit]),    # badframespec
            (-1000, [head; rel("TKA"); spec("MATRIX"); unit]),              # badframespec2
            (-1000, [head; rel("J2000"); "TKFRAME_TKA_RELATIVE = 'J2000'"; spec("MATRIX"); unit]),
            (-1000, [head; rel("J2000"); spec("EULER")]),                   # unknownframespec
            (-1000, [head; rel("J2000"); spec("MATRIX"); "TKFRAME_-1000_MATRIX = (1, 0, 0)"]),
            (-1000, [head; angles("FURLONGS", "(3, 1, 3)")]),               # unitsnotrec
            (-1000, [head; angles("DEGREES", "(3, 4, 3)")]),                # badaxisnumbers
            (-2001, ["FRAME_TKB = -1001", "FRAME_-1001_NAME = 'TKB'", "FRAME_-1001_CLASS = 4",
                "FRAME_-1001_CLASS_ID = -2001", "FRAME_-1001_CENTER = 399",
                "TKFRAME_-2001_RELATIVE = 'J2000'", "TKFRAME_-2001_SPEC = 'QUATERNION'",
                "TKFRAME_-2001_Q = (1, 0, 0, 0)"]),                         # incompleteframe
        )
        r = zeros(9)
        try
            for (id, lines) in cases
                SPICE.kclear()
                SPICE.lmpool(lines)
                tctx = SpiceContext()
                lmpool(tctx, lines)
                want = oracle_tkfram(id)
                @test want isa Symbol && neo_tkfram(tctx, id) === want
                name = id == -1000 ? "TKA" : "TKB"
                ccall((:pxform_c, SPICE.libcspice), Cvoid, (Cstring, Cstring, Cdouble,
                    Ptr{Float64}), name, "J2000", 0.0, r)
                @test throws_kind(() -> pxform(tctx, name, "J2000", 0.0), cspice_kind())
            end
            @test neo_tkfram(SpiceContext(), 0) === oracle_tkfram(0) === :zeroframeid
        finally
            oracle_load(set)
        end
    end

    set = [Kernels.full_set()..., CKFix.frame_kernel(), DynFix.frame_kernel()]
    oracle_load(set)
    ctx = neo_context(set)
    ks = Allspice.kernelset(ctx)

    @apart @testset "tipbod and tisbod" begin
        dp = Float64[]; ds = Float64[]; dd = Float64[]
        structure_ok = true
        # 31008 is the frame class of the lunar binary PCK; the others use text PCK constants.
        for body in (399, 301, 499, 599, 10, 199, 299, 699, 799, 899, 401, 502, 606, 31008),
            ref in ("J2000", "ECLIPJ2000", "B1950", "FK4", "GALACTIC", 1, 17)
            for et in random_epochs(rng, 6, -6.0e8, 6.0e8)
                refname = ref isa Integer ? frmnam(ctx, ref) : ref
                a = tipbod(ctx, ref, body, et)
                b = SPICE.tipbod(refname, body, et)
                push!(dp, maximum(abs.(a .- b)))
                s = tisbod(ctx, ref, body, et)
                t = SPICE.tisbod(refname, body, et)
                push!(ds, maximum(abs.(s[1:3, 1:3] .- t[1:3, 1:3])))
                scale = maximum(abs.(t[4:6, 1:3]))
                push!(dd, maximum(abs.(s[4:6, 1:3] .- t[4:6, 1:3])) / scale)
                structure_ok &= s[1:3, 4:6] == zeros(3, 3) && s[4:6, 4:6] == s[1:3, 1:3]
            end
        end
        sp = report_quantiles("tipbod: element of the matrix", dp)
        ss = report_quantiles("tisbod: element of the rotation block", ds)
        sd = report_quantiles("tisbod: derivative block (relative)", dd)
        @test sp.max <= 1e-13
        @test ss.max <= 1e-13
        @test sd.max <= 1e-12
        @test structure_ok          # [R 0; dR R]
        @test (@inferred tipbod(ks, 1, 399, 1.0e8)) isa SMatrix{3,3,Float64,9}
        @test (@inferred tisbod(ks, 1, 399, 1.0e8)) isa SMatrix{6,6,Float64,36}
        @test throws_kind(() -> tipbod(ctx, "IAU_EARTH", 399, 0.0), :irfnotrec)
        @test throws_kind(() -> tipbod(ctx, "J2000", 123456, 0.0), :framedatanotfound)
    end

    @apart @testset "pckfrm and pckcov" begin
        want = sort!(Int[i for i in SPICE.pckfrm!(SPICE.SpiceIntCell(100), Kernels.MOON_BPC)])
        @test pckfrm(Kernels.MOON_BPC) == want
        @test pckfrm!([7, 3], Kernels.MOON_BPC) == sort!(unique!(vcat([7, 3], want)))
        # Regression: `spkobj!` did not sort its input, so that an input that was not
        # sorted came back with duplicates; `spkobj_c` inserts into a set.
        @test spkobj!([399, 3, 399], Kernels.DE440S) == sort!(unique!(vcat([399, 3], spkobj(Kernels.DE440S))))
        @test !isempty(want)
        for id in vcat(want, [399, 0])
            cell = SPICE.pckcov!(SPICE.SpiceDoubleCell(2000), Kernels.MOON_BPC, id)
            w = Float64[x for x in cell]
            got = pckcov(Kernels.MOON_BPC, id)
            @test got isa Window
            @test collect(got) == w
        end
        pre = Window([-1.0e12, -0.9e12])
        got = pckcov!(pre, Kernels.MOON_BPC, want[1])
        cell = SPICE.SpiceDoubleCell(2000)
        SPICE.wninsd!(cell, -1.0e12, -0.9e12)
        SPICE.pckcov!(cell, Kernels.MOON_BPC, want[1])
        @test collect(got) == Float64[x for x in cell]
        @test throws_kind(() -> pckfrm(Kernels.DE440S), :invalidfiletype)
        @test throws_kind(() -> pckcov(Kernels.LSK, 301), :invalidarchtype)
    end

    @apart @testset "bodfnd and bodvcd" begin
        n = 0
        for body in (399, 301, 499, 10, 599, 1, 3, 9, -77, 2000004, 606),
            item in ("RADII", "PM", "POLE_RA", "POLE_DEC", "GM", "NUT_PREC_RA", "LONG_AXIS",
                "radii", "NOTHING", "NUT_PREC_ANGLES", "MAX_PHASE_DEGREE")
            want = SPICE.bodfnd(body, item)
            @test bodfnd(ctx, body, item) == want
            want || continue
            @test bodvcd(ctx, body, item) == SPICE.bodvcd(body, item)
            n += 1
        end
        println("  bodvcd: ", n, " items equal to CSPICE")
        @test n > 30
        @test throws_kind(() -> bodvcd(ctx, 399, "NOTHING"), :kernelvarnotfound)
        @test (@inferred bodfnd(ctx, 399, "RADII")) === true
    end

    @apart @testset "et2lst" begin
        n = 0
        dsec = Int[]
        # Bodies with a segment in de440s.bsp; Venus rotates in the retrograde direction.
        for body in (399, 301, 299, 199, 10),
            typ in ("PLANETOCENTRIC", "PLANETOGRAPHIC", "planetographic", " Planetocentric ")
            for et in random_epochs(rng, 12, -6.0e8, 6.0e8)
                lon = (2 * rand(rng) - 1) * 2pi
                want = SPICE.et2lst(et, body, lon, strip(typ))
                got = et2lst(ctx, et, body, lon, typ)
                n += 1
                if got == want
                    push!(dsec, 0)
                    continue
                end
                # The second is truncated: one unit of rounding in the angle can move it.
                a = 3600 * got[1] + 60 * got[2] + got[3]
                b = 3600 * want[1] + 60 * want[2] + want[3]
                d = min(mod(a - b, 86400), mod(b - a, 86400))
                push!(dsec, d)
                d <= 1 || @info "et2lst" body typ et lon got want
            end
        end
        nexact = count(==(0), dsec)
        measured(@sprintf("  %-58s n=%-7d equal=%d max=%d s", "et2lst: all five outputs", n,
            nexact, maximum(dsec)))
        @test maximum(dsec) <= 1
        @test nexact >= n - 3
        @test et2lst(ctx, 0.0, 10, 1.0, "PLANETOCENTRIC") == SPICE.et2lst(0.0, 10, 1.0, "PLANETOCENTRIC")
        @test throws_kind(() -> et2lst(ctx, 0.0, 399, 0.0, "GEODETIC"), :unknownsystem)
        # `BODY<id>_PGR_POSITIVE_LON` sets the sense of a planetographic longitude.
        try
            for sense in ("EAST", "WEST")
                SPICE.pcpool("BODY299_PGR_POSITIVE_LON", [sense])
                pcpool(ctx, "BODY299_PGR_POSITIVE_LON", [sense])
                @test et2lst(ctx, 1.0e8, 299, 0.7, "PLANETOGRAPHIC") ==
                      SPICE.et2lst(1.0e8, 299, 0.7, "PLANETOGRAPHIC")
            end
            pcpool(ctx, "BODY299_PGR_POSITIVE_LON", ["NORTH"])
            @test throws_kind(() -> et2lst(ctx, 1.0e8, 299, 0.7, "PLANETOGRAPHIC"), :invalidoption)
        finally
            SPICE.dvpool("BODY299_PGR_POSITIVE_LON")
            dvpool(ctx, "BODY299_PGR_POSITIVE_LON")
        end
        @test throws_kind(() -> et2lst(ctx, 0.0, 1234567, 0.0, "PLANETOCENTRIC"), :cantfindframe)
        r = @inferred et2lst(ctx, 1.0e8, 399, 0.3, "PLANETOCENTRIC")
        @test r isa Tuple{Int,Int,Int,String,String}
    end

    @apart @testset "lspcn" begin
        d = Float64[]
        # Bodies with a segment in de440s.bsp and a pole in the PCK.
        for body in ("EARTH", 399, "MOON", 301, "VENUS", 199, "mercury"),
            ab in ("NONE", "LT", "LT+S", "CN", "CN+S")
            for et in random_epochs(rng, 10, -6.0e8, 6.0e8)
                want = SPICE.lspcn(string(body isa Integer ? bodc2n(ctx, body) : body), et, ab)
                got = lspcn(ctx, body, et, ab)
                x = abs(got - want)
                push!(d, min(x, 2pi - x))
            end
        end
        s = report_quantiles("lspcn: angle (rad)", d)
        @test s.max <= 2e-12
        @test 0 <= lspcn(ctx, "EARTH", 0.0, "NONE") < 2pi
        @test (@inferred lspcn(ks, 399, 1.0e8, "LT+S")) isa Float64
    end

    @apart @testset "bltfrm and kplfrm" begin
        for cls in (-1, 1, 2, 3, 4, 5, 6)
            @test bltfrm(cls) == Int[i for i in SPICE.bltfrm(cls)]
            @test kplfrm(ctx, cls) == Int[i for i in SPICE.kplfrm(cls, 1000)]
        end
        @test length(kplfrm(ctx, -1)) > 20
        for cls in (0, 7, -2)
            @test throws_kind(() -> bltfrm(cls), :badframeclass)
            @test throws_kind(() -> kplfrm(ctx, cls), :badframeclass)
            @test_throws SPICE.SpiceError SPICE.kplfrm(cls, 10)
        end
        # The rules of `kplfrm_` on a pool of its own. CSPICE reads the variables, not the
        # frame definitions: KP_A has its class under its name; KP_B has no FRAME_KP_B and
        # KP_C no class, so both are left out; FRAME_KP_D is rounded to the ID of OTHER,
        # whose name does not match and which is listed; KP_E has a name variable in
        # lower case, which the template matches and the lookup of the ID does not.
        lines = ["FRAME_KP_A = 1700001", "FRAME_1700001_NAME = 'KP_A'", "FRAME_KP_A_CLASS = 3",
            "FRAME_1700002_NAME = 'KP_B'", "FRAME_1700002_CLASS = 4",
            "FRAME_KP_C = 1700003", "FRAME_1700003_NAME = 'KP_C'",
            "FRAME_1700005_NAME = 'KP_D'", "FRAME_KP_D = 1700004.4",
            "FRAME_1700004_NAME = 'OTHER'", "FRAME_1700004_CLASS = 4",
            "frame_1700006_name = 'KP_E'", "FRAME_KP_E = 1700006", "FRAME_1700006_CLASS = 4"]
        try
            SPICE.kclear()
            SPICE.lmpool(lines)
            kc = SpiceContext()
            lmpool(kc, lines)
            for cls in (-1, 3, 4)
                @test kplfrm(kc, cls) == Int[i for i in SPICE.kplfrm(cls, 100)]
            end
            @test kplfrm(kc, -1) == [1700001, 1700004]
        finally
            SPICE.kclear()
        end
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
