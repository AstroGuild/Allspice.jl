# Group 10: orientation kernels (CK types 1 to 6) and CK frames (class 3) against CSPICE.
#
# The kernels are written by `test/ckfixtures.jl` with the CSPICE writers. Error
# measures: absolute difference per element for rotation matrices (elements are at most
# 1), and difference over `CK_RATE` for angular velocities and matrix derivatives, where
# `CK_RATE` is the size of the rotation rates of the test attitude (2e-3 rad/s).
#
# CSPICE keeps state between CK queries, and two defects of that state change its
# answers with the history of the calls (both are shown in the testset "CSPICE answers
# that depend on the call history"). `oracle_ck_reset` clears the state before each
# oracle query, so that the comparison is with the answer of a fresh CSPICE.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of NeoSpice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`NEOSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

isdefined(@__MODULE__, :CKFix) || include("ckfixtures.jl")

const CK_RATE = 2e-3
# Measured: 7e-16 (C-matrix), 2e-16 relative (angular velocity), 0 (clock time).
const CK_TOL_CMAT = 5e-15
const CK_TOL_AV = 5e-15
# Measured: 8e-16 (pxform), 4e-13 relative (sxform derivative block), 9e-16 (spkezr).
const CK_TOL_ROT = 5e-15
const CK_TOL_DROT = 1e-13
const CK_TOL_STATE = 5e-15

"One CK instrument of the fixture: frame name, CK ID, clock, data type, has angular velocity."
const CK_CASES = [
    (label="type 1", frame="TESTSC_T1", inst=-77001, sc=CKFix.SC_TDB, dtype=1, av=true),
    (label="type 1, no av", frame="TESTSC_T1_NOAV", inst=-77011, sc=CKFix.SC_TDB, dtype=1, av=false),
    (label="type 2", frame="TESTSC_T2", inst=-77002, sc=CKFix.SC_TDB, dtype=2, av=true),
    (label="type 3", frame="TESTSC_T3", inst=-77003, sc=CKFix.SC_TDB, dtype=3, av=true),
    (label="type 3, no av", frame="TESTSC_T3_NOAV", inst=-77013, sc=CKFix.SC_TDB, dtype=3, av=false),
    (label="type 3, base IAU_EARTH", frame="TESTSC_T3_EARTH", inst=-77030, sc=CKFix.SC_TDB, dtype=3, av=true),
    (label="type 3, base ECLIPJ2000", frame="TESTSC_OTHERCLOCK", inst=-99001, sc=CKFix.SC_TDB, dtype=3, av=true),
    (label="type 3, TDT clock", frame="TESTSC_TDT", inst=-78100, sc=CKFix.SC_TDT, dtype=3, av=true),
    (label="type 4", frame="TESTSC_T4", inst=-77004, sc=CKFix.SC_TDB, dtype=4, av=true),
    (label="type 5 subtype 0", frame="TESTSC_T5S0", inst=-77050, sc=CKFix.SC_TDB, dtype=5, av=true),
    (label="type 5 subtype 1", frame="TESTSC_T5S1", inst=-77051, sc=CKFix.SC_TDB, dtype=5, av=true),
    (label="type 5 subtype 2", frame="TESTSC_T5S2", inst=-77052, sc=CKFix.SC_TDB, dtype=5, av=true),
    (label="type 5 subtype 3", frame="TESTSC_T5S3", inst=-77053, sc=CKFix.SC_TDB, dtype=5, av=true),
    (label="type 6", frame="TESTSC_T6", inst=-77006, sc=CKFix.SC_TDB, dtype=6, av=true),
]

"Time tags of the segments of `case` in ticks of its clock."
ck_tags(case) = case.sc == CKFix.SC_TDT ? [3000.0 * (10 + i) for i in 1:80] : CKFix.ticks_tdb()

"Instruments of the type 5 segments of the fixture."
const CK5_INSTS = (-77050, -77051, -77052, -77053)

"Clock value strictly inside the first interpolation interval of the type 5 segments."
const CK5_RESET_TICK = CKFix.ticks_tdb()[2] + 1.0

"""
    oracle_ck_reset(inst)

Clear the state that CSPICE keeps for the instrument `inst` between CK queries: a
search far outside the coverage replaces the saved segment of `ckbsr`, and a query in a
gap of the type 6 segment with tolerance zero makes `ckr06` forget its saved
mini-segment. Postcondition: both queries found nothing.

For a type 5 segment a query inside the first interpolation interval follows (it finds
pointing), because `ckr05_` keeps the interval of the last request and keeps it wrong
in one case. A request in a gap that is nearer to the next interval is moved to the
start of that interval (ckr05.c, lines 1187 and 1188: `start = nstart`,
`nstart = nnstrt`), and the routine then saves `nnstrt` as the start after the next
one (line 1300), which is now equal to the saved next start. If the next request
reuses the saved interval (line 1017) and lies in the following gap, nearer to its end,
the window of epochs is empty (line 1257) and CSPICE finds no pointing, although the
same request finds pointing after any other request. A request inside the first
interval saves three correct starts.
"""
function oracle_ck_reset(inst::Integer)
    r = SPICE.ckgp(inst, 1e11, 0.0, "J2000")
    @assert r === nothing
    if inst == -77006
        r = SPICE.ckgp(inst, CKFix.ticks_tdb()[125], 0.0, "J2000")
        @assert r === nothing
    elseif inst in CK5_INSTS
        r = SPICE.ckgp(inst, CK5_RESET_TICK, 0.0, "J2000")
        @assert r !== nothing
    end
    return nothing
end

"Reset CSPICE for every instrument that a frame chain of the fixture can reach."
oracle_frame_reset(inst::Integer) = (oracle_ck_reset(inst); oracle_ck_reset(-77003))

"Result of `f()`, or `nothing` if it throws the error type `E`."
function try_or_nothing(f, ::Type{E}) where {E}
    try
        return f()
    catch e
        e isa E || rethrow()
        return nothing
    end
end
neo_try(f) = try_or_nothing(f, NeoSpice.SpiceError)
oracle_try(f) = try_or_nothing(f, SPICE.SpiceError)

"Random clock value of `case`: inside or just outside the coverage, or exactly a time tag."
function ck_random_tick(rng, case, tags)
    rand(rng) < 0.2 && return tags[rand(rng, eachindex(tags))]
    margin = 0.02 * (tags[end] - tags[1])
    return tags[1] - margin + (tags[end] - tags[1] + 2margin) * rand(rng)
end

"""
    ck_interior_tick(rng, case, tags, range) -> Float64

Clock value strictly inside one piece of the interpolation of `case`, between the tags
`i` and `i + 1` with `i` in `range`, at least a quarter of the spacing from both: for
type 2 inside a stored interval and for type 4 at least 5 s inside a Chebyshev interval.
"""
function ck_interior_tick(rng, case, tags, range)
    while true
        i = rand(rng, range)
        case.dtype == 2 && iseven(i) && continue
        t = tags[i] + (0.25 + 0.35 * rand(rng)) * (tags[i+1] - tags[i])
        if case.dtype == 4
            m = 5 * CKFix.TICKS_PER_SEC
            any(iv -> iv[1] + m <= t <= iv[2] - m, CKFix.type4_intervals()) || continue
        end
        return t
    end
end

"`ckmeta_c` of CSPICE (SPICE.jl has no wrapper of it)."
function oracle_ckmeta(ckid::Integer, meta::AbstractString)
    out = Ref{Cint}(0)
    ccall((:ckmeta_c, SPICE.libcspice), Cvoid, (Cint, Cstring, Ref{Cint}), ckid, meta, out)
    SPICE.handleerror()
    return Int(out[])
end

alloc_ckgp(ks, inst, t, tol, ref) = @allocated ckgp(ks, inst, t, tol, ref)
alloc_ckgpav(ks, inst, t, tol, ref) = @allocated ckgpav(ks, inst, t, tol, ref)
alloc_ck_pxform(ctx, from, to, et) = @allocated pxform(ctx, from, to, et)
alloc_ck_sxform(ctx, from, to, et) = @allocated sxform(ctx, from, to, et)
alloc_ck_pxder(ctx, from, to, et, ::Val{N}) where {N} = @allocated pxder(ctx, from, to, et, Val(N))
alloc_ck_spkezr(ctx, et, ref) = @allocated spkezr(ctx, 301, et, ref, "NONE", 399)

group("10. CK kernels and CK frames against CSPICE", Kernels.LSK, Kernels.PCK11,
    Kernels.DE440S) do
    rng = MersenneTwister(SEED + 10)
    generated = CKFix.all_kernels()
    set = [Kernels.LSK, Kernels.PCK11, Kernels.DE440S, generated...]
    oracle_load(set)
    ctx = neo_context(set)
    ks = NeoSpice.kernelset(ctx)
    ckpaths = CKFix.ck_files()

    @apart @testset "CK files are recognised and listed" begin
        for p in ckpaths
            @test NeoSpice.kernel_arch(p) === :CK
            daf = NeoSpice.open_daf(p)
            @test daf.kind === :CK
            @test isempty(daf.segments) && !isempty(daf.cksegments)
            @test all(s -> 1 <= s.baddr <= s.eaddr <= length(daf.data) && s.t0 <= s.t1,
                daf.cksegments)
            @test ckobj(p) == Int[i for i in SPICE.ckobj(p)]
        end
        @test sort(collect(keys(ks.ck))) == sort([c.inst for c in CK_CASES])
        # Priority order: later file first, and inside a file the later segment first.
        @test [s.avflag for s in ks.ck[Int32(-77003)]] == [true, false, true]
        @test [s.t1 - s.t0 for s in ks.ck[Int32(-77003)]] ==
              sort([s.t1 - s.t0 for s in ks.ck[Int32(-77003)]])
        @test throws_kind(() -> ckobj(Kernels.DE440S), :invalidfiletype)
        @test throws_kind(() -> ckobj("/no/such/file.bc"), :nofile)
    end

    @apart @testset "ckmeta" begin
        for id in (-77001, -77003, -77053, -78100, -99001, -77099, -5, 0, 12, -1000, -999,
            -2147483000), item in ("SCLK", "SPK", " sclk ")
            @test ckmeta(ctx, id, item) == oracle_ckmeta(id, item)
        end
        @test ckmeta(ctx, -99001, "SCLK") == CKFix.SC_TDB
        @test ckmeta(ctx, -77003, "SCLK") == -77           # the default rule
        @test throws_kind(() -> ckmeta(ctx, -77001, "CLOCK"), :unknownckmeta)
    end

    @apart @testset "ckgp and ckgpav" begin
        for case in CK_CASES
            tags = ck_tags(case)
            dcm = Float64[]; dav = Float64[]; dclk = Float64[]
            nfound = 0; nmiss = 0
            agree_found = true
            spacing = tags[2] - tags[1]
            for tol in (0.0, 0.3 * spacing, 4.0 * spacing), _ in 1:400
                t = ck_random_tick(rng, case, tags)
                for ref in ("J2000", "ECLIPJ2000", "IAU_EARTH")
                    oracle_ck_reset(case.inst)
                    r = SPICE.ckgp(case.inst, t, tol, ref)
                    cmat, clkout, found = ckgp(ctx, case.inst, t, tol, ref)
                    agree_found &= found == (r !== nothing)
                    if found && r !== nothing
                        nfound += 1
                        append!(dcm, abs.(cmat .- r[1]))
                        push!(dclk, abs(clkout - r[2]))
                    else
                        nmiss += 1
                    end
                    oracle_ck_reset(case.inst)
                    r = SPICE.ckgpav(case.inst, t, tol, ref)
                    cmat, av, clkout, found = ckgpav(ctx, case.inst, t, tol, ref)
                    agree_found &= found == (r !== nothing)
                    if found && r !== nothing
                        append!(dcm, abs.(cmat .- r[1]))
                        append!(dav, abs.(av .- r[2]) ./ CK_RATE)
                        push!(dclk, abs(clkout - r[3]))
                    end
                end
            end
            @test agree_found
            @test nfound > 300 && nmiss > 20
            a = report_quantiles("ckgp/ckgpav $(case.label): C-matrix", dcm)
            @test a.max <= CK_TOL_CMAT
            @test maximum(dclk) == 0.0
            if case.av
                b = report_quantiles("ckgpav $(case.label): av / $CK_RATE", dav)
                @test b.max <= CK_TOL_AV
            else
                @test isempty(dav)          # no segment with angular velocity: never found
            end
        end
    end

    @apart @testset "ckgp at the pointing instances and with a wrong tolerance" begin
        for case in CK_CASES
            tags = ck_tags(case)
            ok = true
            for i in eachindex(tags)
                oracle_ck_reset(case.inst)
                r = SPICE.ckgp(case.inst, tags[i], 0.0, "J2000")
                cmat, clkout, found = ckgp(ctx, case.inst, tags[i], 0.0, "J2000")
                ok &= found == (r !== nothing)
                found && r !== nothing &&
                    (ok &= maximum(abs.(cmat .- r[1])) <= CK_TOL_CMAT && clkout == r[2])
            end
            @test ok
        end
        # Fixed case of the saved interval of `ckr05_` (see `oracle_ck_reset`): two
        # requests in the gaps before two successive intervals, each nearer to the start
        # of the interval. CSPICE finds no pointing for the second one after the first
        # one, and finds it after a reset; NeoSpice keeps no state and finds it.
        tk = CKFix.ticks_tdb()
        ia, ib = CKFix.START_INDEX[end-1], CKFix.START_INDEX[end]
        ta = tk[ia] - 0.25 * (tk[ia] - tk[ia-1])
        tb = tk[ib] - 0.25 * (tk[ib] - tk[ib-1])
        wide = 4.0 * (tk[2] - tk[1])
        for inst in CK5_INSTS
            oracle_ck_reset(inst)
            @test SPICE.ckgp(inst, ta, wide, "J2000") !== nothing
            @test SPICE.ckgp(inst, tb, wide, "J2000") === nothing      # the CSPICE defect
            oracle_ck_reset(inst)
            r = SPICE.ckgp(inst, tb, wide, "J2000")
            cmat, clkout, found = ckgp(ctx, inst, tb, wide, "J2000")
            @test r !== nothing && found && clkout == r[2] == tk[ib]
            @test maximum(abs.(cmat .- r[1])) <= CK_TOL_CMAT
        end
        @test throws_kind(() -> ckgp(ctx, -77003, 1e6, -1.0, "J2000"), :valueoutofrange)
        @test throws_kind(() -> ckgpav(ctx, -77003, 1e6, -1.0, "J2000"), :valueoutofrange)
        @test throws_kind(() -> ckgp(ctx, -77003, 1e6, 0.0, "NO_SUCH_FRAME"), :unknownframe)
        @test ckgp(ctx, -12345, 1e6, 10.0, "J2000")[3] == false
        @test ckgp(ctx, -77003, 1e11, 10.0, "J2000") == (zero(SMatrix{3,3,Float64}), 0.0, false)
    end

    @apart @testset "priority, tolerance and the angular-velocity rule" begin
        tk = CKFix.ticks_tdb()
        inst = -77003
        base, noav, withav = ks.ck[Int32(inst)][3], ks.ck[Int32(inst)][2], ks.ck[Int32(inst)][1]
        data(seg) = ks.daf[seg.file]
        # Inside the two later segments: `ckgp` takes the one with the highest priority.
        t = tk[150] + 100.0
        @test ckgp(ctx, inst, t, 0.0, 1)[1] == NeoSpice._ck_eval(data(withav), withav, t, 0.0, false)[2]
        # In the later segment without angular velocity only: `ckgp` takes it, `ckgpav`
        # passes over it and takes the first loaded segment.
        t = tk[125] + 100.0
        @test ckgp(ctx, inst, t, 0.0, 1)[1] == NeoSpice._ck_eval(data(noav), noav, t, 0.0, false)[2]
        @test ckgpav(ctx, inst, t, 0.0, 1)[1] == NeoSpice._ck_eval(data(base), base, t, 0.0, true)[2]
        @test ckgp(ctx, inst, t, 0.0, 1)[1] != ckgpav(ctx, inst, t, 0.0, 1)[1]
        # Just outside a segment of higher priority, with a tolerance that reaches it: its
        # end point is returned, and the clock time of the answer is that end point.
        t = tk[140] - 50.0
        cmat, av, clkout, found = ckgpav(ctx, inst, t, 100.0, 1)
        @test found && clkout == tk[140]
        @test cmat == NeoSpice._ck_eval(data(withav), withav, tk[140], 0.0, true)[2]
        oracle_ck_reset(inst)
        r = SPICE.ckgpav(inst, t, 100.0, "J2000")
        @test r[3] == clkout && maximum(abs.(r[1] .- cmat)) <= CK_TOL_CMAT
        # `pxform` follows the rule of `ckgp`, `sxform` the rule of `ckgpav`.
        t = tk[125] + 100.0
        et = sct2e(ctx, CKFix.SC_TDB, t)
        R = pxform(ctx, "J2000", "TESTSC_T3", et)
        X = sxform(ctx, "J2000", "TESTSC_T3", et)
        @test maximum(abs.(R .- ckgp(ctx, inst, t, 0.0, 1)[1])) <= 1e-13
        @test maximum(abs.(X[SOneTo(3), SOneTo(3)] .- ckgpav(ctx, inst, t, 0.0, 1)[1])) <= 1e-13
        @test maximum(abs.(R .- X[SOneTo(3), SOneTo(3)])) > 0.1
    end

    @apart @testset "pxform, sxform and spkezr with CK frames" begin
        others = ("J2000", "IAU_EARTH", CKFix.TK_FRAME, "TESTSC_T5S1")
        for case in CK_CASES
            tags = ck_tags(case)
            dpx = Float64[]; dsx = Float64[]; dds = Float64[]; dst = Float64[]
            agree_avail = true
            structure_ok = true
            navail = 0
            navail_off = 0      # rotations that exist at a tick that is not a time tag
            nticks = nsample(100, 300)
            for _ in 1:nticks
                t = ck_random_tick(rng, case, tags)
                0 <= t || continue
                ontag = t in tags
                et = sct2e(ctx, case.sc, t)
                for other in others
                    other == case.frame && continue
                    oracle_frame_reset(case.inst)
                    Rc = oracle_try(() -> SPICE.pxform(other, case.frame, et))
                    R = neo_try(() -> pxform(ctx, other, case.frame, et))
                    agree_avail &= (R === nothing) == (Rc === nothing)
                    if R !== nothing && Rc !== nothing
                        navail += 1
                        navail_off += !ontag
                        append!(dpx, abs.(R .- Rc))
                        Rb = pxform(ctx, case.frame, other, et)
                        structure_ok &= maximum(abs.(Rb .- R')) <= 1e-14
                    end
                    oracle_frame_reset(case.inst)
                    Xc = oracle_try(() -> SPICE.sxform(other, case.frame, et))
                    X = neo_try(() -> sxform(ctx, other, case.frame, et))
                    agree_avail &= (X === nothing) == (Xc === nothing)
                    if X !== nothing && Xc !== nothing
                        append!(dsx, abs.(X[SOneTo(3), SOneTo(3)] .- Xc[1:3, 1:3]))
                        append!(dds, abs.(X[SVector(4, 5, 6), SOneTo(3)] .- Xc[4:6, 1:3]) ./ CK_RATE)
                        structure_ok &= X[SVector(4, 5, 6), SVector(4, 5, 6)] == X[SOneTo(3), SOneTo(3)]
                        structure_ok &= all(iszero, X[SOneTo(3), SVector(4, 5, 6)])
                    end
                end
                oracle_frame_reset(case.inst)
                sc = oracle_try(() -> SPICE.spkezr("MOON", et, case.frame, "NONE", "EARTH"))
                s = neo_try(() -> spkezr(ctx, "MOON", et, case.frame, "NONE", "EARTH"))
                agree_avail &= (s === nothing) == (sc === nothing)
                if s !== nothing && sc !== nothing
                    r = norm(sc[1][1:3])
                    append!(dst, abs.(s[1][SOneTo(3)] .- sc[1][1:3]) ./ r)
                    append!(dst, abs.(s[1][SVector(4, 5, 6)] .- sc[1][4:6]) ./ (r * CK_RATE))
                    structure_ok &= abs(s[2] - sc[2]) <= 1e-14 * abs(sc[2])
                    oracle_frame_reset(case.inst)
                    pc = SPICE.spkpos("MOON", et, case.frame, "NONE", "EARTH")
                    p = spkpos(ctx, "MOON", et, case.frame, "NONE", "EARTH")
                    append!(dst, abs.(p[1] .- pc[1]) ./ r)
                end
            end
            @test agree_avail
            @test structure_ok
            if case.dtype == 1
                # Tolerance zero: a type 1 frame exists at the pointing instances only.
                # `navail` counts each tick once for each of the other frames, and one
                # tick in five is a time tag, so there is no fixed bound on `navail`;
                # the property is that no rotation exists away from a time tag.
                @test navail_off == 0
                continue
            end
            @test navail > nticks
            a = report_quantiles("pxform $(case.label)", dpx)
            @test a.max <= CK_TOL_ROT
            if case.av
                b = report_quantiles("sxform dR / $CK_RATE $(case.label)", dds)
                c = report_quantiles("spkezr, relative, $(case.label)", dst)
                @test maximum(dsx) <= CK_TOL_ROT
                @test b.max <= CK_TOL_DROT
                @test c.max <= CK_TOL_STATE
            else
                @test isempty(dds)
            end
        end
    end

    @apart @testset "type 1 frames at their pointing instances" begin
        for case in CK_CASES[1:2]
            tags = ck_tags(case)
            n = 0
            ok = true
            for i in eachindex(tags)
                et = sct2e(ctx, case.sc, tags[i])
                oracle_frame_reset(case.inst)
                Rc = oracle_try(() -> SPICE.pxform("J2000", case.frame, et))
                R = neo_try(() -> pxform(ctx, "J2000", case.frame, et))
                ok &= (R === nothing) == (Rc === nothing)
                R === nothing && continue
                Rc === nothing && continue
                n += 1
                ok &= maximum(abs.(R .- Rc)) <= CK_TOL_ROT
                oracle_frame_reset(case.inst)
                Xc = oracle_try(() -> SPICE.sxform("J2000", case.frame, et))
                X = neo_try(() -> sxform(ctx, "J2000", case.frame, et))
                ok &= (X === nothing) == (Xc === nothing) == !case.av
                X === nothing || Xc === nothing ||
                    (ok &= maximum(abs.(X .- Xc)) <= CK_TOL_ROT)
                # The exact derivative of a stored instance is zero.
                ok &= all(iszero, pxder(ctx, "J2000", case.frame, et, Val(1))[2])
            end
            @test ok
            @test n > 100       # the clock conversion returns most of the tags exactly
        end
    end

    @apart @testset "derivatives of the interpolated rotation, and the CSPICE angular velocity" begin
        measured("  exact derivative (pxder) against CSPICE angular velocity (sxform), |ΔdR| / $CK_RATE:")
        for case in CK_CASES
            (case.dtype == 1 || !case.av) && continue
            tags = ck_tags(case)
            tdt = case.sc == CKFix.SC_TDT
            # Two ranges of tags: in the first one the clock runs at its nominal rate
            # (first coefficient record), in the second one 1e-5 faster (second record).
            ranges = tdt ? (2:28, 31:78) : (2:35, 181:225)
            h = 0.05
            fd1 = Float64[]; fd2 = Float64[]
            dav = (Float64[], Float64[])
            consistent = true
            for (k, range) in enumerate(ranges), _ in 1:60
                t = ck_interior_tick(rng, case, tags, range)
                et = sct2e(ctx, case.sc, t)
                J = pxder(ctx, "J2000", case.frame, et, Val(2))
                Jp = pxder(ctx, "J2000", case.frame, et + h, Val(1))
                Jm = pxder(ctx, "J2000", case.frame, et - h, Val(1))
                consistent &= maximum(abs.(J[1] .- pxform(ctx, "J2000", case.frame, et))) <= 1e-14
                consistent &= maximum(abs.(pxder(ctx, "J2000", case.frame, et, Val(1))[2] .- J[2])) <= 1e-14 * CK_RATE
                append!(fd1, abs.(J[2] .- (Jp[1] .- Jm[1]) ./ 2h) ./ CK_RATE)
                append!(fd2, abs.(J[3] .- (Jp[2] .- Jm[2]) ./ 2h) ./ CK_RATE^2)
                X = sxform(ctx, "J2000", case.frame, et)
                append!(dav[k], abs.(J[2] .- X[SVector(4, 5, 6), SOneTo(3)]) ./ CK_RATE)
                # The rotation stays a rotation to first order: Rᵀ dR is antisymmetric.
                w = J[1]' * J[2]
                consistent &= maximum(abs.(w .+ w')) <= 1e-12 * CK_RATE * 1e3
            end
            @test consistent
            a = report_quantiles("pxder order 1 against finite difference, $(case.label)", fd1)
            b = report_quantiles("pxder order 2 against finite difference, $(case.label)", fd2)
            @test a.max <= 1e-6
            @test b.max <= 1e-5
            measured(@sprintf("    %-26s nominal clock rate: max %.2e   clock 1e-5 fast: max %.2e",
                case.label, maximum(dav[1]), maximum(dav[2])))
            # The two derivatives describe the same motion: they differ by the
            # interpolation error of the data type, never by more than 1 percent here.
            @test max(maximum(dav[1]), maximum(dav[2])) <= 1e-2
            if case.dtype == 2 || case.label in ("type 5 subtype 0", "type 5 subtype 1")
                # Here the CSPICE angular velocity is the derivative of the interpolated
                # rotation with respect to the nominal clock, so the two agree where the
                # clock runs at the nominal rate, and differ by the rate offset elsewhere.
                @test maximum(dav[1]) <= 1e-9
                @test 1e-6 <= maximum(dav[2]) <= 2e-5
            end
        end
    end

    @apart @testset "no allocation, inferred types" begin
        tk = CKFix.ticks_tdb()
        for case in CK_CASES
            tags = ck_tags(case)
            t = case.dtype == 1 ? tags[17] : ck_interior_tick(rng, case, tags, 2:28)
            @test CKFix.warm(alloc_ckgp, ks, case.inst, t, 0.0, 1) == 0
            @test CKFix.warm(alloc_ckgp, ks, case.inst, t, 5000.0, 17) == 0
            @test CKFix.warm(alloc_ckgpav, ks, case.inst, t, 5000.0, 17) == 0
            @test (@inferred ckgp(ks, case.inst, t, 0.0, 1)) isa Tuple{SMatrix{3,3,Float64,9},Float64,Bool}
            @test (@inferred ckgpav(ks, case.inst, t, 0.0, 1)) isa
                  Tuple{SMatrix{3,3,Float64,9},SVector{3,Float64},Float64,Bool}
            case.dtype == 1 && continue
            et = sct2e(ctx, case.sc, t)
            id = Int(namfrm(ctx, case.frame))
            iau = Int(namfrm(ctx, "IAU_EARTH"))
            @test CKFix.warm(alloc_ck_pxform, ctx, 1, id, et) == 0
            @test CKFix.warm(alloc_ck_pxform, ctx, id, iau, et) == 0
            @test CKFix.warm(alloc_ck_pxder, ctx, 1, id, et, Val(1)) == 0
            @test CKFix.warm(alloc_ck_pxder, ctx, iau, id, et, Val(3)) == 0
            @test (@inferred pxform(ctx, 1, id, et)) isa SMatrix{3,3,Float64,9}
            @test (@inferred pxder(ctx, 1, id, et, Val(2))) isa NTuple{3,SMatrix{3,3,Float64,9}}
            @test (@inferred NeoSpice.frame_rotation(ks, Int32(1), Int32(id),
                NeoSpice.taylor_var(et, Val(2)))) isa SMatrix{3,3,NeoSpice.Taylor{2,Float64},9}
            case.av || continue
            @test CKFix.warm(alloc_ck_sxform, ctx, 1, id, et) == 0
            @test CKFix.warm(alloc_ck_sxform, ctx, iau, id, et) == 0
            @test CKFix.warm(alloc_ck_spkezr, ctx, et, id) == 0
            @test (@inferred sxform(ctx, 1, id, et)) isa SMatrix{6,6,Float64,36}
            @test (@inferred spkezr(ctx, 301, et, id, "NONE", 399)) isa Tuple{SVector{6,Float64},Float64}
        end
        # A fixed-offset frame on top of a CK frame, and two CK frames.
        et = sct2e(ctx, CKFix.SC_TDB, tk[20] + 1000.0)
        cam = Int(namfrm(ctx, CKFix.TK_FRAME))
        t5 = Int(namfrm(ctx, "TESTSC_T5S1"))
        @test CKFix.warm(alloc_ck_pxform, ctx, cam, t5, et) == 0
        @test CKFix.warm(alloc_ck_sxform, ctx, cam, t5, et) == 0
        @test CKFix.warm(alloc_ck_pxder, ctx, cam, t5, et, Val(2)) == 0
    end

    @apart @testset "ckcov" begin
        intervals(cell) = [Tuple(SPICE.wnfetd(cell, i)) for i in 1:SPICE.wncard(cell)]
        nwindows = 0
        for p in ckpaths, id in ckobj(p), needav in (false, true),
            level in ("SEGMENT", "INTERVAL"), tol in (0.0, 700.0), timsys in ("SCLK", "TDB")
            c = intervals(SPICE.ckcov(p, id, needav, level, tol, timsys))
            n = ckcov(ctx, p, id, needav, level, tol, timsys)
            same = length(c) == length(n) &&
                   all(maximum(abs.(collect(a) .- collect(b))) <= 1e-6 for (a, b) in zip(c, n))
            same || @info "ckcov differs" p id needav level tol timsys length(c) length(n)
            @test same
            nwindows += !isempty(n)
        end
        @test nwindows > 100
        p = ckpaths[1]
        @test ckcov(ctx, p, -77003, false, "segment", 0.0, "sclk") ==
              [(CKFix.ticks_tdb()[1], CKFix.ticks_tdb()[end])]
        @test length(ckcov(ctx, p, -77002, false, "INTERVAL", 0.0, "SCLK")) ==
              length(CKFix.type2_intervals())
        @test isempty(ckcov(ctx, p, -77011, true, "INTERVAL", 0.0, "SCLK"))
        @test isempty(ckcov(ctx, p, -4242, false, "INTERVAL", 0.0, "SCLK"))
        @test throws_kind(() -> ckcov(ctx, p, -77003, false, "RECORD", 0.0, "SCLK"), :invalidoption)
        @test throws_kind(() -> ckcov(ctx, p, -77003, false, "SEGMENT", 0.0, "UTC"), :unsupported)
        @test throws_kind(() -> ckcov(ctx, p, -77003, false, "SEGMENT", -1.0, "SCLK"), :valueoutofrange)
        @test throws_kind(() -> ckcov(ctx, Kernels.DE440S, 399, false, "SEGMENT", 0.0, "SCLK"),
            :invalidfiletype)
    end

    @apart @testset "frames without data" begin
        tk = CKFix.ticks_tdb()
        et = sct2e(ctx, CKFix.SC_TDB, tk[20] + 1000.0)
        # A CK frame without any segment, an epoch outside the coverage, a frame whose
        # segments have no angular velocity, and an epoch outside the clock.
        @test throws_kind(() -> pxform(ctx, "J2000", "TESTSC_NODATA", et), :noframeconnect)
        @test_throws SPICE.SpiceError SPICE.pxform("J2000", "TESTSC_NODATA", et)
        @test throws_kind(() -> pxform(ctx, "J2000", "TESTSC_T3", et + 5e4), :noframeconnect)
        @test_throws SPICE.SpiceError SPICE.pxform("J2000", "TESTSC_T3", et + 5e4)
        @test pxform(ctx, "J2000", "TESTSC_T3_NOAV", et) isa SMatrix
        @test throws_kind(() -> sxform(ctx, "J2000", "TESTSC_T3_NOAV", et), :noframeconnect)
        @test throws_kind(() -> spkezr(ctx, "MOON", et, "TESTSC_T3_NOAV", "NONE", "EARTH"),
            :noframeconnect)
        @test spkpos(ctx, "MOON", et, "TESTSC_T3_NOAV", "NONE", "EARTH")[1] isa SVector
        @test throws_kind(() -> pxform(ctx, "J2000", "TESTSC_T3", -1e6), :valueoutofrange)
        @test_throws SPICE.SpiceError SPICE.pxform("J2000", "TESTSC_T3", -1e6)
        # The CK frame to itself and to the frame fixed to it need no CK data.
        @test pxform(ctx, "TESTSC_NODATA", "TESTSC_NODATA", et) == one(SMatrix{3,3,Float64})
        Rc = SPICE.pxform(CKFix.TK_FRAME, "TESTSC_T3", et + 5e4)
        @test maximum(abs.(pxform(ctx, CKFix.TK_FRAME, "TESTSC_T3", et + 5e4) .- Rc)) <= CK_TOL_ROT
        Xc = SPICE.sxform(CKFix.TK_FRAME, "TESTSC_T3", et + 5e4)
        @test maximum(abs.(sxform(ctx, CKFix.TK_FRAME, "TESTSC_T3", et + 5e4) .- Xc)) <= CK_TOL_ROT
        # Without the SCLK kernel, and without the CK files, the frame has no connection.
        noclock = neo_context([Kernels.LSK, CKFix.frame_kernel(), ckpaths...])
        @test throws_kind(() -> pxform(noclock, "J2000", "TESTSC_T3", et), :noframeconnect)
        @test ckgp(noclock, -77003, tk[20], 0.0, "J2000")[3]
        @test throws_kind(() -> ckgp(noclock, -77003, tk[20], 0.0, "IAU_EARTH"), :kernelvarnotfound)
        nock = neo_context([Kernels.LSK, CKFix.sclk_kernel(), CKFix.frame_kernel()])
        @test throws_kind(() -> pxform(nock, "J2000", "TESTSC_T3", et), :noframeconnect)
        @test throws_kind(() -> sxform(nock, "J2000", "TESTSC_T3", et), :noframeconnect)
        # Unloading the CK file of higher priority gives the first profile back.
        t = tk[150] + 100.0
        e150 = sct2e(ctx, CKFix.SC_TDB, t)
        two = neo_context(set)
        before = pxform(two, "J2000", "TESTSC_T3", e150)
        unload(two, ckpaths[3])
        after = pxform(two, "J2000", "TESTSC_T3", e150)
        @test maximum(abs.(before .- after)) > 0.1
        @test after == pxform(two, "J2000", "TESTSC_T3_NOAV", e150)
    end

    @apart @testset "CSPICE answers that depend on the call history" begin
        tk = CKFix.ticks_tdb()
        # 1. `ckr06` keeps the mini-segment of the last call. When a request in a gap is
        #    moved to the next mini-segment, it rounds the clock rate of that
        #    mini-segment to a whole number (here 0) and keeps it, so the next request
        #    in that mini-segment is interpolated with wrong derivatives. NeoSpice has
        #    no state and gives the answer of a fresh CSPICE.
        t = tk[192] + 0.4 * (tk[193] - tk[192])
        oracle_ck_reset(-77006)
        fresh = SPICE.ckgpav(-77006, t, 300.0, "J2000")
        SPICE.ckgpav(-77006, tk[130] - 100.0, 300.0, "J2000")
        stale = SPICE.ckgpav(-77006, t, 300.0, "J2000")
        ours = ckgpav(ctx, -77006, t, 300.0, "J2000")
        @test maximum(abs.(ours[1] .- fresh[1])) <= CK_TOL_CMAT
        @test ours == ckgpav(ctx, -77006, t, 300.0, "J2000")
        measured(@sprintf("  CSPICE type 6, same request after a gap request: C-matrix changes by %.2e",
            maximum(abs.(stale[1] .- fresh[1]))))
        # 2. `ckbsr` keeps the segment of the last search. After a search that needs
        #    angular velocity, a search that does not need it takes that saved segment
        #    although a segment of higher priority (without angular velocity) covers
        #    the request. NeoSpice applies the priority rule of the CK Required Reading.
        t = tk[125] + 100.0
        oracle_ck_reset(-77003)
        fresh = SPICE.ckgp(-77003, t, 0.0, "J2000")
        SPICE.ckgpav(-77003, t, 0.0, "J2000")
        stale = SPICE.ckgp(-77003, t, 0.0, "J2000")
        ours = ckgp(ctx, -77003, t, 0.0, "J2000")
        @test maximum(abs.(ours[1] .- fresh[1])) <= CK_TOL_CMAT
        measured(@sprintf("  CSPICE ckgp after ckgpav, same request: C-matrix changes by %.2e",
            maximum(abs.(stale[1] .- fresh[1]))))
        oracle_ck_reset(-77003)
        oracle_ck_reset(-77006)
    end

    @apart @testset "concurrent readers" begin
        if Threads.nthreads() == 1
            skip_notice("10. concurrent CK readers", "one thread only")
        else
            jobs = Tuple{String,Float64}[]
            for case in CK_CASES, _ in 1:40
                (case.dtype == 1 || !case.av) && continue
                tags = ck_tags(case)
                push!(jobs, (case.frame, sct2e(ctx, case.sc, ck_interior_tick(rng, case, tags, 2:28))))
            end
            serial = [(pxform(ctx, "IAU_EARTH", f, et), sxform(ctx, "J2000", f, et)) for (f, et) in jobs]
            results = Vector{Bool}(undef, 16)
            @sync for k in 1:16
                Threads.@spawn begin
                    good = true
                    for rep in 1:20, (i, (f, et)) in enumerate(jobs)
                        good &= pxform(ctx, "IAU_EARTH", f, et) == serial[i][1]
                        good &= sxform(ctx, "J2000", f, et) == serial[i][2]
                    end
                    results[k] = good
                end
            end
            @test all(results)
        end
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
