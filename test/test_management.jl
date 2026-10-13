# Group 6: kernel management.
#
# Feature files: features/unload.feature, features/load_order_precedence.feature,
# features/coverage_error.feature

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

"Write a meta-kernel that loads `files` (all in one directory `dir`) and return its path."
function write_metakernel(tmp::AbstractString, dir::AbstractString, files)
    @assert all(f -> dirname(f) == dir, files)
    @assert length(dir) <= 80 "PATH_VALUES entry longer than the CSPICE limit"
    path = joinpath(tmp, "allspice_test.tm")
    list = join(("'\$K/" * basename(f) * "'" for f in files), ",\n                        ")
    write(path, """
    KPL/MK

    Meta-kernel written by test_management.jl.

    \\begindata

       PATH_VALUES     = ( '$dir' )
       PATH_SYMBOLS    = ( 'K' )
       KERNELS_TO_LOAD = ( $list )

    \\begintext
    """)
    return path
end

group("6. kernel management", Kernels.FULL_SET_FILES..., Kernels.DE405, Kernels.PCK10) do
    rng = MersenneTwister(SEED + 6)
    # Epochs inside both ephemerides.
    epochs = random_epochs(rng, 20, DE440S_COVERAGE...)
    query(ctx, et) = spkezr(ctx, 301, et, "J2000", "NONE", 10)
    cquery(et) = SPICE.spkezr("301", et, "J2000", "NONE", "10")

    @apart @testset "furnsh then unload restores the earlier answers" begin
        ctx = neo_context([Kernels.LSK, Kernels.DE440S])
        @test ktotal(ctx) == 2
        before = [query(ctx, et) for et in epochs]
        furnsh(ctx, Kernels.DE405)
        @test ktotal(ctx) == 3
        during = [query(ctx, et) for et in epochs]
        @test all(during[i] != before[i] for i in eachindex(epochs))
        unload(ctx, Kernels.DE405)
        @test ktotal(ctx) == 2
        @test [query(ctx, et) for et in epochs] == before
        # The same for a text kernel: the pool goes back to the earlier values.
        furnsh(ctx, Kernels.PCK11)
        r11 = bodvrd(ctx, "NEPTUNE", "PM")
        R11 = pxform(ctx, "J2000", "IAU_NEPTUNE", epochs[1])
        M11 = pxform(ctx, "J2000", "IAU_MARS", epochs[1])
        furnsh(ctx, Kernels.PCK10)
        @test bodvrd(ctx, "NEPTUNE", "PM") != r11
        @test pxform(ctx, "J2000", "IAU_NEPTUNE", epochs[1]) != R11
        # pck00010 after pck00011 leaves the 26 Mars amplitudes of pck00011 with the one
        # angle pair of pck00010; CSPICE signals SPICE(INSUFFICIENTANGLES) for that pool.
        @test throws_kind(() -> pxform(ctx, "J2000", "IAU_MARS", epochs[1]), :insufficientangles)
        unload(ctx, Kernels.PCK10)
        @test bodvrd(ctx, "NEPTUNE", "PM") == r11
        @test pxform(ctx, "J2000", "IAU_NEPTUNE", epochs[1]) == R11
        @test pxform(ctx, "J2000", "IAU_MARS", epochs[1]) == M11
        # Unloading a kernel that is not loaded changes nothing, as in CSPICE.
        unload(ctx, "/no/such/kernel.bsp")
        @test ktotal(ctx) == 3
        @test [query(ctx, et) for et in epochs] == before
    end

    @apart @testset "kclear" begin
        ctx = neo_context(Kernels.full_set())
        @test ktotal(ctx) == length(Kernels.full_set())
        kclear(ctx)
        @test ktotal(ctx) == 0
        @test_throws Allspice.SpiceError query(ctx, epochs[1])
        @test_throws Allspice.SpiceError pxform(ctx, "J2000", "MOON_PA", epochs[1])
        @test throws_kind(() -> gdpool(ctx, "BODY399_RADII"), :kernelvarnotfound)
        # Built-in names and frames need no kernel.
        @test bodn2c(ctx, "EARTH") == 399
        @test namfrm(ctx, "ECLIPJ2000") == 17
        @test pxform(ctx, "J2000", "ECLIPJ2000", 0.0) isa SMatrix{3,3,Float64,9}
        # The context can be used again.
        furnsh(ctx, Kernels.DE440S)
        @test ktotal(ctx) == 1
        @test query(ctx, epochs[1])[1] isa SVector{6,Float64}
    end

    @apart @testset "loading order decides precedence" begin
        for order in ([Kernels.DE405, Kernels.DE440S], [Kernels.DE440S, Kernels.DE405])
            oracle_load(order)
            ctx = neo_context(order)
            label = join(basename.(order), " then ")
            d = compare_states(ctx, [(301, 10), (399, 10), (301, 399), (5, 0), (4, 399)],
                epochs, "J2000")
            assert_states("precedence, $label", d, TOL_J2000)
        end
        first405 = neo_context([Kernels.DE405, Kernels.DE440S])     # de440s wins
        first440 = neo_context([Kernels.DE440S, Kernels.DE405])     # de405 wins
        only440 = neo_context([Kernels.DE440S])
        only405 = neo_context([Kernels.DE405])
        for et in epochs
            @test query(first405, et) == query(only440, et)
            @test query(first440, et) == query(only405, et)
            @test query(first405, et) != query(first440, et)
        end
        # Outside the coverage of the later file the earlier file answers.
        late = DE405_COVERAGE[2] - 86400.0 * 100          # year 2201: only in de405
        @test late > DE440S_COVERAGE[2]
        @test query(first405, late) == query(only405, late)
        oracle_load([Kernels.DE405, Kernels.DE440S])
        @test isapprox(query(first405, late)[1], cquery(late)[1]; rtol=1e-13)
        # Text kernels: the later assignment of a variable wins.
        a = neo_context([Kernels.PCK10, Kernels.PCK11])
        b = neo_context([Kernels.PCK11, Kernels.PCK10])
        oracle_load([Kernels.PCK10, Kernels.PCK11])
        for body in ("NEPTUNE", "MERCURY", "MARS", "JUPITER")
            @test bodvrd(a, body, "PM") == SPICE.bodvrd(body, "PM")
            @test isapprox(pxform(a, "J2000", "IAU_" * body, 1e8),
                SPICE.pxform("J2000", "IAU_" * body, 1e8); atol=1e-10)
        end
        oracle_load([Kernels.PCK11, Kernels.PCK10])
        for body in ("NEPTUNE", "MERCURY", "JUPITER")
            @test bodvrd(b, body, "PM") == SPICE.bodvrd(body, "PM")
            @test isapprox(pxform(b, "J2000", "IAU_" * body, 1e8),
                SPICE.pxform("J2000", "IAU_" * body, 1e8); atol=1e-10)
        end
        @test bodvrd(a, "NEPTUNE", "PM") != bodvrd(b, "NEPTUNE", "PM")
        @test maximum(abs.(pxform(a, "J2000", "IAU_NEPTUNE", 1e8) - pxform(b, "J2000", "IAU_NEPTUNE", 1e8))) > 1e-3
        # The mixed Mars constants are an error in CSPICE and in Allspice.
        @test_throws SPICE.SpiceError SPICE.pxform("J2000", "IAU_MARS", 1e8)
        @test throws_kind(() -> pxform(b, "J2000", "IAU_MARS", 1e8), :insufficientangles)
    end

    @apart @testset "meta-kernel in a temporary directory" begin
        mktempdir() do tmp
            files = [Kernels.LSK, Kernels.PCK11, Kernels.MOON_TF, Kernels.MOON_BPC]
            mk = write_metakernel(tmp, Kernels.DIR, files)
            oracle_load([mk, Kernels.DE440S])
            ctx = neo_context([mk, Kernels.DE440S])
            direct = neo_context([files..., Kernels.DE440S])
            # CSPICE counts the meta-kernel and each file that it lists.
            @test ktotal(ctx) == SPICE.ktotal("ALL")
            @test ktotal(ctx) == length(files) + 2
            # The meta-kernel variables do not stay in the pool.
            @test SPICE.gcpool("KERNELS_TO_LOAD") === nothing
            @test throws_kind(() -> gcpool(ctx, "KERNELS_TO_LOAD"), :kernelvarnotfound)
            @test throws_kind(() -> gcpool(ctx, "PATH_VALUES"), :kernelvarnotfound)
            for et in epochs[1:5]
                @test spkezr(ctx, "MOON", et, "MOON_PA_DE440", "NONE", "EARTH") ==
                      spkezr(direct, "MOON", et, "MOON_PA_DE440", "NONE", "EARTH")
                cs, _ = SPICE.spkezr("MOON", et, "MOON_PA_DE440", "NONE", "EARTH")
                @test isapprox(spkezr(ctx, "MOON", et, "MOON_PA_DE440", "NONE", "EARTH")[1],
                    cs; rtol=1e-9)
            end
            @test utc2et(ctx, "2020-01-01T00:00:00") ≈ SPICE.utc2et("2020-01-01T00:00:00") atol = 1e-6
            @test bodvrd(ctx, "EARTH", "RADII") == SPICE.bodvrd("EARTH", "RADII")
            # Unloading the meta-kernel unloads everything that it loaded.
            unload(ctx, mk)
            SPICE.unload(mk)
            @test ktotal(ctx) == SPICE.ktotal("ALL") == 1
            @test throws_kind(() -> gdpool(ctx, "BODY399_RADII"), :kernelvarnotfound)
            @test spkezr(ctx, 301, epochs[1], "J2000", "NONE", 399)[1] isa SVector{6,Float64}
            @test_throws Allspice.SpiceError pxform(ctx, "J2000", "MOON_PA_DE440", epochs[1])
            # A missing file gives an error and leaves the context unchanged.
            @test throws_kind(() -> furnsh(ctx, joinpath(tmp, "absent.bsp")), :nosuchfile)
            @test ktotal(ctx) == 1
            # A meta-kernel that lists a meta-kernel is an error, as in CSPICE.
            outer = joinpath(tmp, "outer.tm")
            write(outer, "KPL/MK\n\\begindata\nKERNELS_TO_LOAD = ( '$mk' )\n\\begintext\n")
            @test_throws SPICE.SpiceError SPICE.furnsh(outer)
            @test throws_kind(() -> furnsh(ctx, outer), :recursiveloading)
            @test ktotal(ctx) == 1
        end
    end

    @apart @testset "errors" begin
        ctx = neo_context(Kernels.full_set())
        lo, hi = -4.734072e9, 4.735368e9
        # Out of coverage: on both sides, and after the end in every form of the call.
        @test throws_kind(() -> spkezr(ctx, 301, lo - 1e6, "J2000", "NONE", 399), :spkinsuffdata)
        et = hi + 1e6
        @test throws_kind(() -> spkezr(ctx, 301, et, "J2000", "NONE", 399), :spkinsuffdata)
        @test throws_kind(() -> spkezr(ctx, "MOON", et, "J2000", "NONE", "EARTH"), :spkinsuffdata)
        @test throws_kind(() -> spkpos(ctx, 10, et, "J2000", "NONE", 0), :spkinsuffdata)
        @test throws_kind(() -> spkgeo(ctx, 399, et, "J2000", 10), :spkinsuffdata)
        @test throws_kind(() -> spkder(ctx, 301, et, "J2000", 399, Val(2)), :spkinsuffdata)
        # The ends of the coverage are inside it, as in CSPICE.
        oracle_load(Kernels.full_set())
        for et in (lo, hi)
            cs, _ = SPICE.spkezr("301", et, "J2000", "NONE", "399")
            @test isapprox(spkezr(ctx, 301, et, "J2000", "NONE", 399)[1], cs; rtol=1e-12)
        end
        @test_throws SPICE.SpiceError SPICE.spkezr("301", hi + 1e6, "J2000", "NONE", "399")
        # Unknown body, by name; a body with no ephemeris data, by ID.
        @test throws_kind(() -> spkezr(ctx, "NO_SUCH_BODY", 0.0, "J2000", "NONE", "EARTH"), :idcodenotfound)
        @test throws_kind(() -> spkezr(ctx, "EARTH", 0.0, "J2000", "NONE", "NO_SUCH_BODY"), :idcodenotfound)
        @test throws_kind(() -> spkder(ctx, "NO_SUCH_BODY", 0.0, "J2000", "EARTH", Val(1)), :idcodenotfound)
        @test_throws Allspice.SpiceError spkezr(ctx, 599, 0.0, "J2000", "NONE", 399)
        @test_throws Allspice.SpiceError spkezr(ctx, "PHOBOS", 0.0, "J2000", "NONE", "EARTH")
        # Unknown frame.
        @test throws_kind(() -> spkezr(ctx, 301, 0.0, "NO_SUCH_FRAME", "NONE", 399), :unknownframe)
        @test throws_kind(() -> pxform(ctx, "J2000", "NO_SUCH_FRAME", 0.0), :unknownframe)
        @test throws_kind(() -> sxform(ctx, "NO_SUCH_FRAME", "J2000", 0.0), :unknownframe)
        @test throws_kind(() -> pxder(ctx, "J2000", "NO_SUCH_FRAME", 0.0, Val(1)), :unknownframe)
        @test throws_kind(() -> pxform(ctx, 1, 987654, 0.0), :unknownframe)       # by frame ID
        # A correction string that is not valid (group 9 has the full list).
        @test throws_kind(() -> spkezr(ctx, 301, 0.0, "J2000", "RL+S", 399), :invalidoption)
        @test spkezr(ctx, 301, 0.0, "J2000", " none ", 399) == spkezr(ctx, 301, 0.0, "J2000", "NONE", 399)
        # A frame whose data are not loaded (binary PCK of DE421 is not on disk).
        if isfile(Kernels.MOON_TF_OLD)
            old = neo_context([Kernels.text_set()..., Kernels.MOON_TF_OLD, Kernels.MOON_BPC])
            @test_throws Allspice.SpiceError pxform(old, "J2000", "MOON_PA_DE421", 0.0)
        end
        # The error is printable and names the kind.
        err = try spkezr(ctx, 301, 1e11, "J2000", "NONE", 399) catch e; e end
        @test err isa Allspice.SpiceError
        @test occursin("spkinsuffdata", sprint(showerror, err))
        @test !isempty(err.msg)
        # Missing SPK data carry the query and the body without data. `loaded` tells an
        # epoch outside the data (the kind `:nocoverage` of earlier versions) from a body
        # that is in no file (`:nopath`).
        caught(f) = try f(); nothing catch e; e end
        data(e) = (e.kind, e.target, e.observer, e.et, e.body, e.loaded)
        err = caught(() -> spkezr(ctx, 301, 1e11, "J2000", "NONE", 399))
        @test err isa SpkInsufficientDataError
        @test data(err) == (:spkinsuffdata, 301, 399, 1e11, 399, true)
        err = caught(() -> spkpos(ctx, "JUPITER", 0.0, "J2000", "NONE", "EARTH"))
        @test err isa SpkInsufficientDataError
        @test data(err) == (:spkinsuffdata, 599, 399, 0.0, 599, false)
        err = caught(() -> spkgeo(SpiceContext(), 301, 0.0, "J2000", 399))
        @test err isa SpkInsufficientDataError
        @test data(err) == (:noloadedfiles, 301, 399, 0.0, nothing, false)
        # Every other failure has the kind and the message only.
        err = caught(() -> spkezr(ctx, 301, 0.0, "NO_SUCH_FRAME", "NONE", 399))
        @test err isa SpiceKindError && err isa SpiceError
        @test SpiceError(:badarg, "text") isa SpiceKindError
        @test sprint(showerror, SpiceError(:badarg, "text")) == "SpiceError(badarg): text"
    end

    @apart @testset "the kind of an error is the short message of CSPICE" begin
        # Regression: these calls had kinds of their own (`:nofile`, `:badarg`,
        # `:unknownbody`, `:notfound`, `:badtype`, `:baddaf`, `:unsupported`), so that the
        # same failure had two kinds in sibling calls. One case for each kind that changed.
        function ckind(f)                       # short message of the direct CSPICE call `f`
            f()
            ccall((:failed_c, SPICE.libcspice), Cint, ()) == 0 && return :none
            buf = zeros(UInt8, 64)
            ccall((:getmsg_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Ptr{UInt8}), "SHORT", 64, buf)
            ccall((:reset_c, SPICE.libcspice), Cvoid, ())
            return Symbol(lowercase(replace(unsafe_string(pointer(buf)), "SPICE(" => "", ")" => "")))
        end
        nkind(f) = try f(); :none catch e; e isa Allspice.SpiceError ? e.kind : rethrow() end
        oracle_load(Kernels.full_set())
        ctx = neo_context(Kernels.full_set())
        SPICE.pcpool("BODY399_MGMT_TEXT", ["a"]); pcpool(ctx, "BODY399_MGMT_TEXT", ["a"])
        none = joinpath(mktempdir(), "none.bc")
        st, lt, n = zeros(6), Ref(0.0), Ref{Cint}(0)
        ids = SPICE.SpiceIntCell(10); cover = SPICE.SpiceDoubleCell(10)
        ezr(targ, corr) = ccall((:spkezr_c, SPICE.libcspice), Cvoid, (Cstring, Cdouble, Cstring, Cstring,
            Cstring, Ptr{Cdouble}, Ptr{Cdouble}), targ, 0.0, "J2000", corr, "EARTH", st, lt)
        pos(targ, corr) = ccall((:spkpos_c, SPICE.libcspice), Cvoid, (Cstring, Cdouble, Cstring, Cstring,
            Cstring, Ptr{Cdouble}, Ptr{Cdouble}), targ, 0.0, "J2000", corr, "EARTH", st, lt)
        value(body, item) = ccall((:bodvrd_c, SPICE.libcspice), Cvoid, (Cstring, Cstring, Cint, Ptr{Cint},
            Ptr{Cdouble}), body, item, 6, n, st)
        cases = (
            (:nosuchfile, () -> furnsh(ctx, none), () -> ccall((:furnsh_c, SPICE.libcspice), Cvoid, (Cstring,), none)),
            (:blankfilename, () -> furnsh(ctx, " "), () -> ccall((:furnsh_c, SPICE.libcspice), Cvoid, (Cstring,), " ")),
            (:fileopenfailed, () -> ldpool(ctx, none), () -> ccall((:ldpool_c, SPICE.libcspice), Cvoid, (Cstring,), none)),
            (:filenotfound, () -> eklef(ctx, none),
                () -> ccall((:eklef_c, SPICE.libcspice), Cvoid, (Cstring, Ptr{Cint}), none, n)),
            (:filenotfound, () -> ckobj(none),
                () -> ccall((:ckobj_c, SPICE.libcspice), Cvoid, (Cstring, Ref{SPICE.Cell{SPICE.SpiceInt}}), none, ids.cell)),
            (:invalidarchtype, () -> ckobj(Kernels.LSK), () -> ccall((:ckobj_c, SPICE.libcspice), Cvoid,
                (Cstring, Ref{SPICE.Cell{SPICE.SpiceInt}}), Kernels.LSK, ids.cell)),
            (:idcodenotfound, () -> spkezr(ctx, "NO_SUCH_BODY", 0.0, "J2000", "NONE", "EARTH"),
                () -> ezr("NO_SUCH_BODY", "NONE")),
            (:invalidoption, () -> spkezr(ctx, "MOON", 0.0, "J2000", "RL+S", "EARTH"),
                () -> ezr("MOON", "RL+S")),
            (:spkinvalidoption, () -> spkpos(ctx, "MOON", 0.0, "J2000", "RL+S", "EARTH"),
                () -> pos("MOON", "RL+S")),
            (:spkinvalidoption, () -> spkapo(ctx, 301, 0.0, "J2000", zeros(6), "RL+S"),
                () -> ccall((:spkapo_c, SPICE.libcspice), Cvoid, (Cint, Cdouble, Cstring, Ptr{Cdouble}, Cstring,
                    Ptr{Cdouble}, Ptr{Cdouble}), 301, 0.0, "J2000", zeros(6), "RL+S", st, lt)),
            (:notranslation, () -> bodvrd(ctx, "NO_SUCH_BODY", "RADII"), () -> value("NO_SUCH_BODY", "RADII")),
            (:kernelvarnotfound, () -> bodvrd(ctx, "EARTH", "NO_SUCH_ITEM"), () -> value("EARTH", "NO_SUCH_ITEM")),
            (:typemismatch, () -> bodvrd(ctx, "EARTH", "MGMT_TEXT"), () -> value("EARTH", "MGMT_TEXT")),
        )
        for (kind, neo, c) in cases
            @test nkind(neo) === ckind(c) === kind
        end
        SPICE.dvpool("BODY399_MGMT_TEXT")
        # Regression: missing SPK data and missing leap seconds had the kinds `:nocoverage`,
        # `:nopath` and `:nolsk` of this package. With SPK files loaded CSPICE has one
        # message for an epoch outside the data and for a body in no file; with no SPK
        # file it has another; without leap seconds the message depends on the call.
        geo(targ, et) = ccall((:spkgeo_c, SPICE.libcspice), Cvoid, (Cint, Cdouble, Cstring, Cint,
            Ptr{Cdouble}, Ptr{Cdouble}), targ, et, "J2000", 399, st, lt)
        @test nkind(() -> spkgeo(ctx, 301, 1.0e11, "J2000", 399)) === ckind(() -> geo(301, 1.0e11)) === :spkinsuffdata
        @test nkind(() -> spkgeo(ctx, 599, 0.0, "J2000", 399)) === ckind(() -> geo(599, 0.0)) === :spkinsuffdata
        SPICE.kclear()
        bare = SpiceContext()
        buf = zeros(UInt8, 40)
        bare_cases = (
            (:noloadedfiles, () -> spkgeo(bare, 301, 0.0, "J2000", 399), () -> geo(301, 0.0)),
            (:noleapseconds, () -> str2et(bare, "2000 JAN 1"),
                () -> ccall((:str2et_c, SPICE.libcspice), Cvoid, (Cstring, Ptr{Cdouble}), "2000 JAN 1", lt)),
            (:missingtimeinfo, () -> et2utc(bare, 0.0, "C", 3), () -> ccall((:et2utc_c, SPICE.libcspice), Cvoid,
                (Cdouble, Cstring, Cint, Cint, Ptr{UInt8}), 0.0, "C", 3, 40, buf)),
            (:missingtimeinfo, () -> timout(bare, 0.0, "YYYY"), () -> ccall((:timout_c, SPICE.libcspice), Cvoid,
                (Cdouble, Cstring, Cint, Ptr{UInt8}), 0.0, "YYYY", 40, buf)),
            (:missingtimeinfo, () -> unitim(bare, 0.0, "TDB", "TAI"), () -> ccall((:unitim_c, SPICE.libcspice),
                Cdouble, (Cdouble, Cstring, Cstring), 0.0, "TDB", "TAI")),
            (:kernelvarnotfound, () -> deltet(bare, 0.0, "ET"), () -> ccall((:deltet_c, SPICE.libcspice), Cvoid,
                (Cdouble, Cstring, Ptr{Cdouble}), 0.0, "ET", lt)),
        )
        for (kind, neo, c) in bare_cases
            @test nkind(neo) === ckind(c) === kind
        end
    end

    @apart @testset "the registry of kinds holds the kinds that the source throws" begin
        # The kinds are read from the syntax trees of `src/`: each lower-case symbol in
        # the first argument of `SpiceError(...)` or of the constructor of a subtype (a
        # literal or a conditional of literals),
        # in an argument of a call that a function hands on to `SpiceError` (directly or
        # through another such function), and on the right of `kind = ...`.
        # Not seen: a kind that reaches `SpiceError` from a table or a field. The
        # constructor refuses a kind that is not registered, so the suite fails where
        # such a kind is thrown.
        name(f) = f isa Symbol ? f : f isa Expr && f.head === :. && f.args[end] isa QuoteNode ?
                                     f.args[end].value : nothing
        islower(s) = all(c -> islowercase(c) || isdigit(c), string(s))
        function symbols!(out, x)
            x isa QuoteNode && x.value isa Symbol && islower(x.value) && push!(out, x.value)
            x isa Expr && foreach(a -> symbols!(out, a), x.head === :. ? x.args[1:1] : x.args)
            return out
        end
        mentions(x, s) = x === s || (x isa Expr && any(a -> mentions(a, s), x.args))
        positional(x) = filter(a -> !(a isa Expr && a.head === :parameters), x.args[2:end])
        argname(a) = a isa Symbol ? a : a isa Expr && a.head in (:(::), :kw) && !isempty(a.args) ?
                                        argname(a.args[1]) : nothing
        # the positions of the arguments of a call of `x.args[1]` that hold a kind
        errtypes = (:SpiceError, :SpiceKindError, :SpkInsufficientDataError)
        kindargs(x, helpers) = name(x.args[1]) in errtypes ? (1,) : get(helpers, name(x.args[1]), ())
        hands_on(x, arg, helpers) = x isa Expr && ((x.head === :call && any(p -> p <= length(positional(x)) &&
            mentions(positional(x)[p], arg), kindargs(x, helpers))) || any(a -> hands_on(a, arg, helpers), x.args))
        function helpers!(helpers, x)             # true if a helper or a position was added
            x isa Expr || return false
            added = false
            if x.head in (:function, :(=)) && x.args[1] isa Expr && length(x.args) == 2
                sig = x.args[1]
                while sig isa Expr && sig.head in (:where, :(::))
                    sig = sig.args[1]
                end
                if sig isa Expr && sig.head === :call && !(name(sig.args[1]) in errtypes)
                    for (pos, a) in enumerate(positional(sig))
                        an = argname(a)
                        (an === nothing || pos in get(helpers, name(sig.args[1]), ())) && continue
                        if hands_on(x.args[2], an, helpers)
                            helpers[name(sig.args[1])] = (get(helpers, name(sig.args[1]), ())..., pos)
                            added = true
                        end
                    end
                end
            end
            for a in x.args
                added |= helpers!(helpers, a)
            end
            return added
        end
        function kinds!(out, other, x, helpers)
            x isa Expr || return out
            if x.head === :call
                args = positional(x)
                for pos in kindargs(x, helpers)
                    pos <= length(args) || continue
                    symbols!(out, args[pos])
                    k = args[pos]
                    !(name(x.args[1]) in errtypes) || k isa QuoteNode || argname(k) isa Symbol ||
                        (k isa Expr && k.head === :if) || push!(other, k)
                end
            elseif x.head === :(=) && x.args[1] === :kind
                symbols!(out, x.args[2])
            end
            foreach(a -> kinds!(out, other, a, helpers), x.args)
            return out
        end
        trees = [Meta.parseall(read(joinpath(root, f), String)) for (root, _, files) in
                 walkdir(joinpath(pkgdir(Allspice), "src")) for f in files if endswith(f, ".jl")]
        helpers = Dict{Any,Tuple}()
        while any(t -> helpers!(helpers, t), trees)   # until no function hands a kind on to a new one
        end
        thrown, other = Set{Symbol}(), Any[]
        foreach(t -> kinds!(thrown, other, t, helpers), trees)
        @test isempty(other)                     # a first argument of a form that is not read
        @test issorted(Allspice.ERROR_KINDS) && allunique(Allspice.ERROR_KINDS)
        @test sort!(collect(setdiff(thrown, Allspice.ERROR_KINDS))) == Symbol[]     # thrown, not registered
        @test sort!(collect(setdiff(Allspice.ERROR_KINDS, thrown))) == Symbol[]     # registered, never thrown
        @test_throws ArgumentError SpiceError(:nosuchkind, "a kind that is not registered")
    end

    @apart @testset "contexts are independent" begin
        kclear()                                         # default context: nothing loaded
        c440 = neo_context([Kernels.DE440S])
        c405 = neo_context([Kernels.DE405])
        @test ktotal() == 0
        @test_throws Allspice.SpiceError spkezr(301, epochs[1], "J2000", "NONE", 399)
        a = [query(c440, et) for et in epochs]
        b = [query(c405, et) for et in epochs]
        @test all(a[i] != b[i] for i in eachindex(epochs))
        # Loading into the default context changes neither explicit context.
        furnsh(Kernels.DE405)
        @test ktotal() == 1 && ktotal(c440) == 1 && ktotal(c405) == 1
        @test [query(c440, et) for et in epochs] == a
        @test [spkezr(301, et, "J2000", "NONE", 10) for et in epochs] == b
        # Changing one explicit context changes neither the other nor the default.
        furnsh(c440, Kernels.PCK11)
        @test throws_kind(() -> gdpool(c405, "BODY399_RADII"), :kernelvarnotfound)
        @test throws_kind(() -> gdpool("BODY399_RADII"), :kernelvarnotfound)
        kclear(c440)
        @test ktotal(c440) == 0 && ktotal(c405) == 1 && ktotal() == 1
        @test [query(c405, et) for et in epochs] == b
        @test [spkezr(301, et, "J2000", "NONE", 10) for et in epochs] == b
        kclear()
        @test [query(c405, et) for et in epochs] == b
        @test default_context() === default_context()
        @test default_context() !== c405
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
