# Group 17: dynamic frames (class 5) and switch frames (class 6) against CSPICE.
#
# The frame kernels are written by `test/dynfixtures.jl`: every frame family, every kind
# of defining vector, frozen frames, frames with an inertial rotation state, dynamic
# frames that are defined relative to dynamic frames, and switch frames.
#
# `pxform` is compared element by element. The derivative block of `sxform` and the
# velocity of `spkezr` are compared relative to their size. CSPICE takes the derivative
# of a velocity vector by a central difference over about one second, so its derivative
# block has a rounding noise of about 1e-16 divided by the rotation rate of the frame;
# NeoSpice makes the same difference, and the two results agree to that noise only. The
# thresholds below are therefore per frame, each about ten times the measured value.

isdefined(@__MODULE__, :CKFix) || include("ckfixtures.jl")
isdefined(@__MODULE__, :DynFix) || include("dynfixtures.jl")

"""
    df_oracle_kind(call, from, to, et) -> Symbol

Short message of the error that CSPICE signals for `pxform_c` (`call = :pxform_c`) or
`sxform_c` from frame `from` to frame `to` at `et`, as a lower-case symbol without
`SPICE(...)`; `:none` without an error. The C functions are called directly, because
SPICE.jl keeps the long message only.
"""
function df_oracle_kind(call::Symbol, from::AbstractString, to::AbstractString, et::Float64)
    out = zeros(36)
    if call === :pxform_c
        ccall((:pxform_c, SPICE.libcspice), Cvoid, (Cstring, Cstring, Cdouble, Ptr{Float64}),
            from, to, et, out)
    else
        ccall((:sxform_c, SPICE.libcspice), Cvoid, (Cstring, Cstring, Cdouble, Ptr{Float64}),
            from, to, et, out)
    end
    ccall((:failed_c, SPICE.libcspice), Cint, ()) == 0 && return :none
    buf = zeros(UInt8, 64)
    ccall((:getmsg_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Ptr{UInt8}), "SHORT", 64, buf)
    ccall((:reset_c, SPICE.libcspice), Cvoid, ())
    short = unsafe_string(pointer(buf))
    m = match(r"SPICE\(([A-Z0-9_]+)\)", short)
    return m === nothing ? :unknown : Symbol(lowercase(m.captures[1]))
end

"Kind of the `SpiceError` that NeoSpice throws for `f()`; `:none` without an error."
function df_neo_kind(f)
    try
        f()
    catch e
        e isa NeoSpice.SpiceError || rethrow()
        return e.kind
    end
    return :none
end

"""
    df_rel(a, b, floor=0.0)

Largest element of `|a - b|`, divided by the largest element of `|b|` or by `floor` if
that is larger. If both arrays are below 1e-17 the result is zero: the derivative block
of a frame that does not rotate is the rounding noise of a cancellation in both engines.
"""
function df_rel(a, b, floor=0.0)
    s = maximum(abs, b)
    (s <= 1e-17 && maximum(abs, a) <= 1e-17) && return 0.0
    return maximum(abs.(a .- b)) / max(s, floor, floatmin())
end

df_alloc_px(ks, from, to, et) = @allocated pxform(ks, from, to, et)
df_alloc_sx(ks, from, to, et) = @allocated sxform(ks, from, to, et)
df_alloc_der(ks, from, to, et, ::Val{N}) where {N} = @allocated pxder(ks, from, to, et, Val(N))
df_alloc_ezr(ks, et, ref, ab) = @allocated spkezr(ks, 301, et, ref, ab, 399)
df_alloc_pos(ks, et, ref, ab) = @allocated spkpos(ks, 301, et, ref, ab, 399)

"Five-point central difference of the matrix function `f` at `t` with the step `h`."
df_fd(f, t, h) = (8 .* (f(t + h) .- f(t - h)) .- (f(t + 2h) .- f(t - 2h))) ./ (12h)

"""
Thresholds per frame, `(pxform, sxform derivative block, spkezr position, spkezr
velocity)`: the default, and the frames whose derivative comes from a numerical
difference in CSPICE.
"""
const DF_DEFAULT_TOL = (5e-14, 2e-12, 5e-14, 2e-12)
const DF_TOL = Dict(
    "DT_GSE" => (5e-14, 5e-9, 5e-14, 5e-10),
    "DT_GSE_LTS" => (5e-14, 5e-9, 5e-14, 5e-10),
    "DT_GSE_CN" => (5e-14, 5e-9, 5e-14, 5e-10),
    "DT_NADIR" => (5e-14, 5e-9, 5e-14, 5e-10),
    "DT_NADIR_LT" => (5e-14, 5e-9, 5e-14, 5e-10),
    "DT_VEL_PCK" => (5e-14, 5e-9, 5e-14, 5e-7),
    "DT_CONST_ABCORR" => (5e-14, 5e-9, 5e-14, 5e-10),
    "DT_EULER_ON_GSE" => (5e-14, 5e-9, 5e-14, 5e-10),
    "DT_SWITCH_TIME" => (5e-14, 5e-9, 5e-14, 5e-10),
    "DT_SWITCH_NEST" => (5e-14, 5e-9, 5e-14, 5e-10),
    "DT_PRODUCT" => (5e-14, 5e-9, 5e-14, 5e-10),
    # A velocity that is given in a dynamic frame: CSPICE differentiates a quantity
    # that already holds a numerical derivative.
    "DT_TWOVEC_IN_DYN" => (5e-10, 5e-4, 5e-10, 5e-4),
)

group("17. dynamic frames and switch frames against CSPICE", Kernels.LSK, Kernels.PCK11,
    Kernels.DE440S) do
    rng = MersenneTwister(SEED + 17)
    set = [Kernels.LSK, Kernels.PCK11, Kernels.DE440S, DynFix.frame_kernel()]
    oracle_load(set)
    ctx = neo_context(set)
    ks = NeoSpice.kernelset(ctx)
    span = 5.5e8
    nepoch = 150

    @testset "frame definitions are read" begin
        for (name, family, center) in DynFix.FRAMES
            id = namfrm(ctx, name)
            @test id == SPICE.namfrm(name) != 0
            @test frinfo(ctx, id) == SPICE.frinfo(id)
            @test frinfo(ctx, id)[1] == center
            class = frinfo(ctx, id)[2]
            @test class in (5, 6)
            @test haskey(class == 5 ? ks.dyn : ks.switch, Int32(id))
        end
        @test length(ks.dyn) + length(ks.switch) == length(DynFix.FRAMES)
        @test all(d -> d.errkind === :none, values(ks.dyn))
        @test all(d -> d.errkind === :none, values(ks.switch))
        # A snapshot without a dynamic frame has no inner snapshot.
        @test NeoSpice.kernelset(neo_context([Kernels.LSK, Kernels.PCK11])).inner === nothing
        @test ks.inner !== nothing && ks.inner.ks.inner !== nothing
        @test ks.inner.ks.inner.ks.inner === nothing
    end

    rows = String[]
    @testset "pxform, sxform and spkezr for each frame" begin
        for (name, family, center) in DynFix.FRAMES
            tol = get(DF_TOL, name, DF_DEFAULT_TOL)
            dp = 0.0; ds = 0.0; dr = 0.0; dv = 0.0; dra = 0.0; dva = 0.0
            for et in random_epochs(rng, nepoch, -span, span)
                a = pxform(ctx, "J2000", name, et)
                b = SPICE.pxform("J2000", name, et)
                dp = max(dp, maximum(abs.(a .- b)))
                a6 = sxform(ctx, "J2000", name, et)
                b6 = SPICE.sxform("J2000", name, et)
                dp = max(dp, maximum(abs.(a6[1:3, 1:3] .- b6[1:3, 1:3])))
                ds = max(ds, df_rel(a6[4:6, 1:3], b6[4:6, 1:3]))
                for (targ, obs) in (("MOON", "EARTH"), ("SUN", "MOON"))
                    s, lt = spkezr(ctx, targ, et, name, "NONE", obs)
                    t, lu = SPICE.spkezr(targ, et, name, "NONE", obs)
                    dr = max(dr, norm(s[1:3] .- t[1:3]) / norm(t[1:3]))
                    dv = max(dv, norm(s[4:6] .- t[4:6]) / norm(t[4:6]))
                    s, lt = spkezr(ctx, targ, et, name, "LT+S", obs)
                    t, lu = SPICE.spkezr(targ, et, name, "LT+S", obs)
                    dra = max(dra, norm(s[1:3] .- t[1:3]) / norm(t[1:3]))
                    dva = max(dva, norm(s[4:6] .- t[4:6]) / norm(t[4:6]))
                    @test abs(lt - lu) <= 1e-12 * lu
                end
            end
            push!(rows, @sprintf("  %-17s px=%.1e sx=%.1e ezr NONE r=%.1e v=%.1e  LT+S r=%.1e v=%.1e  %s",
                name, dp, ds, dr, dv, dra, dva, family))
            ok = dp <= tol[1] && ds <= tol[2] && max(dr, dra) <= tol[3] && max(dv, dva) <= tol[4]
            ok || @info "above the threshold" name tol dp ds dr dv dra dva
            @test dp <= tol[1]
            @test ds <= tol[2]
            @test max(dr, dra) <= tol[3]
            @test max(dv, dva) <= tol[4]
        end
        measured("  largest difference from CSPICE, $nepoch epochs per frame (px: element of the " *
                 "matrix; sx: derivative block, relative; ezr: relative)")
        foreach(measured, rows)
    end

    @testset "other corrections, spkpos, and frame to frame" begin
        d = Float64[]; dv = Float64[]; dq = Float64[]; dx = Float64[]
        names = first.(DynFix.FRAMES)
        for ab in ("LT", "CN", "CN+S", "XLT", "XLT+S", "XCN", "XCN+S"), name in names
            name == "DT_TWOVEC_IN_DYN" && continue
            for et in random_epochs(rng, 6, -span, span)
                s, lt = spkezr(ctx, 301, et, name, ab, 10)
                t, lu = SPICE.spkezr("301", et, name, ab, "10")
                push!(d, norm(s[1:3] .- t[1:3]) / norm(t[1:3]))
                push!(dv, norm(s[4:6] .- t[4:6]) / norm(t[4:6]))
                p, lp = spkpos(ctx, 301, et, name, ab, 10)
                q, lq = SPICE.spkpos("301", et, name, ab, "10")
                push!(dq, norm(p .- q) / norm(q))
            end
        end
        for _ in 1:400
            a, b = rand(rng, names), rand(rng, names)
            et = (2 * rand(rng) - 1) * 2.9e8
            push!(dx, maximum(abs.(pxform(ctx, a, b, et) .- SPICE.pxform(a, b, et))))
        end
        s1 = report_quantiles("spkezr, 7 other corrections: position, relative", d)
        s2 = report_quantiles("spkezr, 7 other corrections: velocity, relative", dv)
        s3 = report_quantiles("spkpos, 7 other corrections: position, relative", dq)
        s4 = report_quantiles("pxform between two dynamic or switch frames", dx)
        @test s1.max <= 5e-13
        @test s2.max <= 5e-7
        @test s3.max <= 5e-13
        @test s4.max <= 5e-10
        @test s4.q99 <= 5e-13
    end

    @testset "pxder against finite differences and sxform" begin
        lines = String[]
        for (name, family, center) in DynFix.FRAMES
            id = namfrm(ctx, name)
            def = get(ks.dyn, Int32(id), nothing)
            switch = def === nothing
            frozen = !switch && def.frozen
            inertial = !switch && def.inertial
            loose = name == "DT_TWOVEC_IN_DYN"
            e1 = 0.0; e2 = 0.0; es = 0.0
            for et in record_safe_epochs(rng, 12; span=2.8e8)
                # Not near an epoch at which a switch frame changes its base.
                switch && any(e -> abs(et - e) < 200, DynFix.SWITCH_EDGES) && continue
                R, dR, ddR = pxder(ctx, "J2000", name, et, Val(2))
                # Order 0 is the rotation of `pxform`, computed with series arithmetic.
                @test maximum(abs.(R .- pxform(ctx, "J2000", name, et))) <= (loose ? 1e-9 : 1e-13)
                if frozen
                    @test iszero(dR) && iszero(ddR)
                    @test iszero(sxform(ctx, "J2000", name, et)[4:6, 1:3])
                    continue
                end
                h = 16.0
                fd1 = df_fd(t -> pxform(ctx, "J2000", name, t), et, h)
                fd2 = df_fd(t -> pxder(ctx, "J2000", name, t, Val(1))[2], et, h)
                # The floors are the rounding noise of the differences: a frame that
                # turns slower than this is compared in absolute terms.
                e1 = max(e1, df_rel(dR, fd1, 1e-10))
                e2 = max(e2, df_rel(ddR, fd2, 1e-15))
                sx = sxform(ctx, "J2000", name, et)[4:6, 1:3]
                if inertial
                    @test iszero(sx)
                    @test !iszero(dR)
                else
                    es = max(es, df_rel(sx, dR, 1e-13))
                end
            end
            frozen && continue
            push!(lines, @sprintf("  %-17s order 1 vs FD=%.1e  order 2 vs FD=%.1e  sxform vs order 1=%.1e",
                name, e1, e2, es))
            @test e1 <= (loose ? 1e-3 : 2e-6)
            @test e2 <= (loose ? 1e-1 : 2e-5)
            # The derivative block of CSPICE is not the derivative of its rotation for
            # every definition (see the docstring of `sxform`); the difference is small.
            @test es <= (loose ? 1e-3 : 5e-5)
        end
        measured("  pxder of the rotation from J2000 to each frame (relative to the largest element)")
        foreach(measured, lines)
        # Orders 0 to 3 are consistent with each other.
        for name in ("DT_GSE", "DT_ETE_OF_DATE", "DT_EULER", "DT_GSM", "DT_SWITCH_ALL", "DT_NADIR_LT")
            et = record_safe_epochs(rng, 1; span=2.8e8)[1]
            j3 = pxder(ctx, "J2000", name, et, Val(3))
            j1 = pxder(ctx, "J2000", name, et, Val(1))
            j0 = pxder(ctx, "J2000", name, et, Val(0))
            @test maximum(abs.(j0[1] .- j3[1])) <= 1e-13 && maximum(abs.(j1[1] .- j3[1])) <= 1e-13
            @test df_rel(j1[2], j3[2]) <= 1e-12
        end
    end

    @testset "frozen frames and inertial rotation state" begin
        for (frozen, live, t0) in (("DT_GSE_FROZEN", "DT_GSE", "2015-MAR-03 01:02:03 TDB"),
            ("DT_ETE_FROZEN", "DT_ETE_OF_DATE", "1995-JUN-15 00:00:00 TDB"))
            et0 = str2et(ctx, t0)
            fixed = pxform(ctx, "J2000", live, et0)
            for et in random_epochs(rng, 20, -span, span)
                @test pxform(ctx, "J2000", frozen, et) == fixed
            end
        end
        for (inert, live) in (("DT_ETE_INERTIAL", "DT_ETE_OF_DATE"), ("DT_GSE_INERTIAL", "DT_GSE"))
            for et in random_epochs(rng, 20, -span, span)
                @test pxform(ctx, "J2000", inert, et) == pxform(ctx, "J2000", live, et)
                @test iszero(sxform(ctx, "J2000", inert, et)[4:6, 1:3])
                @test iszero(SPICE.sxform("J2000", inert, et)[4:6, 1:3])
                @test pxder(ctx, "J2000", inert, et, Val(1))[2] == pxder(ctx, "J2000", live, et, Val(1))[2]
            end
        end
    end

    @testset "switch frames at the ends of the intervals" begin
        worst = 0.0
        n = 0
        for edge in DynFix.SWITCH_EDGES, delta in (-1.0, -1e-3, 0.0, 1e-3, 1.0),
            name in ("DT_SWITCH_TIME", "DT_SWITCH_NEST")
            et = edge + delta
            want = df_oracle_kind(:pxform_c, "J2000", name, et)
            got = df_neo_kind(() -> pxform(ctx, "J2000", name, et))
            want == got || @info "switch frame" name et want got
            @test got == want
            got === :none || continue
            worst = max(worst, maximum(abs.(pxform(ctx, "J2000", name, et) .- SPICE.pxform("J2000", name, et))))
            n += 1
        end
        measured(@sprintf("  %-58s n=%-7d max=%.3e", "switch frames at interval ends: pxform element", n, worst))
        @test worst <= 5e-14
        @test n >= 40
        # The base that is listed last has priority where the intervals overlap.
        @test pxform(ctx, "J2000", "DT_SWITCH_TIME", 2.0e8) == pxform(ctx, "J2000", "DT_EULER", 2.0e8)
        @test pxform(ctx, "J2000", "DT_SWITCH_TIME", 0.0) == pxform(ctx, "J2000", "ECLIPJ2000", 0.0)
        @test pxform(ctx, "J2000", "DT_SWITCH_TIME", -1.5e8) == pxform(ctx, "J2000", "DT_GSE", -1.5e8)
        @test pxform(ctx, "J2000", "DT_SWITCH_TIME", -5.0e8) == pxform(ctx, "J2000", "IAU_EARTH", -5.0e8)
        @test pxform(ctx, "J2000", "DT_SWITCH_ALL", 1.0e7) == pxform(ctx, "J2000", "DT_ETE_OF_DATE", 1.0e7)
        @test pxform(ctx, "J2000", "DT_SWITCH_NEST", 4.0e8) == pxform(ctx, "J2000", "IAU_MOON", 4.0e8)
        @test throws_kind(() -> pxform(ctx, "J2000", "DT_SWITCH_TIME", 7.0e8), :noframeconnect)
        @test throws_kind(() -> sxform(ctx, "J2000", "DT_SWITCH_TIME", -7.0e8), :noframeconnect)
        @test throws_kind(() -> spkezr(ctx, 301, 7.0e8, "DT_SWITCH_TIME", "NONE", 399), :noframeconnect)
    end

    @testset "no allocation and inferred types" begin
        lts = Aberration("LT+S"); none = Aberration("NONE"); cn = Aberration("XCN+S")
        et = 1.2345e8
        nalloc = 0
        for (name, family, center) in DynFix.FRAMES
            id = namfrm(ctx, name)
            pxform(ks, 1, id, et); sxform(ks, 1, id, et); pxder(ks, 1, id, et, Val(2))
            pxder(ks, id, 17, et, Val(1))
            for ab in (lts, none, cn)
                spkezr(ks, 301, et, id, ab, 399); spkpos(ks, 301, et, id, ab, 399)
            end
            a = (df_alloc_px(ks, 1, id, et), df_alloc_px(ks, id, 10013, et), df_alloc_sx(ks, 1, id, et),
                df_alloc_der(ks, 1, id, et, Val(2)), df_alloc_der(ks, id, 17, et, Val(1)),
                df_alloc_ezr(ks, et, id, lts), df_alloc_ezr(ks, et, id, none), df_alloc_ezr(ks, et, id, cn),
                df_alloc_pos(ks, et, id, lts), df_alloc_pos(ks, et, id, none))
            all(==(0), a) || @info "allocation" name a
            @test all(==(0), a)
            nalloc += length(a)
            @test (@inferred pxform(ks, 1, id, et)) isa SMatrix{3,3,Float64,9}
            @test (@inferred sxform(ks, 1, id, et)) isa SMatrix{6,6,Float64,36}
            @test (@inferred pxder(ks, 1, id, et, Val(2))) isa NTuple{3,SMatrix{3,3,Float64,9}}
            @test (@inferred spkezr(ks, 301, et, id, lts, 399)) isa Tuple{SVector{6,Float64},Float64}
            @test (@inferred spkezr(ks, 301, et, id, none, 399)) isa Tuple{SVector{6,Float64},Float64}
            @test (@inferred spkpos(ks, 301, et, id, lts, 399)) isa Tuple{SVector{3,Float64},Float64}
        end
        println("  ", nalloc, " calls with zero allocation")
        # The rotation is generic over the number type of the epoch.
        tv = NeoSpice.taylor_var(et, Val(2))
        id = namfrm(ctx, "DT_GSE_LTS")
        m = @inferred NeoSpice.frame_rotation(ks, Int32(1), Int32(id), tv)
        @test m isa SMatrix{3,3,NeoSpice.Taylor{2,Float64},9}
        @test maximum(abs.(map(x -> x.v, m) .- pxform(ks, 1, id, et))) <= 1e-13
    end

    if Threads.nthreads() == 1
        skip_notice("17. threads", "one thread only")
    else
        @testset "threads: the same answers from many tasks" begin
            names = first.(DynFix.FRAMES)
            jobs = [(rand(rng, names), (2 * rand(rng) - 1) * 2.9e8) for _ in 1:400]
            serial = [(pxform(ctx, "J2000", n, et), sxform(ctx, "J2000", n, et),
                spkezr(ctx, 301, et, n, "LT+S", 399)[1]) for (n, et) in jobs]
            # Order 3 was not used before: the tasks race to make its function pointer.
            results = Vector{Any}(undef, length(jobs))
            der = Vector{Any}(undef, length(jobs))
            Threads.@threads for i in eachindex(jobs)
                n, et = jobs[i]
                der[i] = pxder(ctx, "J2000", n, et, Val(3))
                results[i] = (pxform(ctx, "J2000", n, et), sxform(ctx, "J2000", n, et),
                    spkezr(ctx, 301, et, n, "LT+S", 399)[1])
            end
            @test all(i -> results[i] == serial[i], eachindex(jobs))
            @test all(i -> der[i] == pxder(ctx, "J2000", jobs[i][1], jobs[i][2], Val(3)), eachindex(jobs))
        end
    end

    @testset "definitions that are not valid, and the recursion limit" begin
        set2 = [set..., DynFix.error_kernel()]
        oracle_load(set2)
        ctx2 = neo_context(set2)     # the load itself does not throw
        et = 1.0e8
        lines = String[]
        for (name, kind) in DynFix.ERROR_FRAMES
            @test namfrm(ctx2, name) == SPICE.namfrm(name) != 0
            want = df_oracle_kind(:pxform_c, "J2000", name, et)
            got = df_neo_kind(() -> pxform(ctx2, "J2000", name, et))
            push!(lines, @sprintf("  %-17s NeoSpice :%-20s CSPICE :%s", name, got, want))
            got == want || @info "error kind" name got want
            @test got == want
            @test got == kind
            @test df_neo_kind(() -> sxform(ctx2, "J2000", name, et)) == kind
            @test df_neo_kind(() -> spkezr(ctx2, 301, et, name, "LT+S", 399)) == kind
            if kind === :none
                @test maximum(abs.(pxform(ctx2, "J2000", name, et) .- SPICE.pxform("J2000", name, et))) <= 5e-10
            end
        end
        foreach(println, lines)
        # The valid frames of the same context are not affected.
        @test maximum(abs.(pxform(ctx2, "J2000", "DT_GSE", et) .- SPICE.pxform("J2000", "DT_GSE", et))) <= 5e-14
        # A class 5 frame without any definition, and a class above 6.
        extra = joinpath(mktempdir(), "empty.tf")
        write(extra, "KPL/FK\n\\begindata\n" * DynFix.head("DT_EMPTY", 1850131, 5) *
                     DynFix.head("DT_CLASS7", 1850132, 7) * "\\begintext\n")
        SPICE.furnsh(extra); furnsh(ctx2, extra)
        @test df_neo_kind(() -> pxform(ctx2, "J2000", "DT_EMPTY", et)) !== :none
        @test df_oracle_kind(:pxform_c, "J2000", "DT_EMPTY", et) !== :none
        k7 = df_neo_kind(() -> pxform(ctx2, "J2000", "DT_CLASS7", et))
        println("  class 5 frame without a definition: NeoSpice :",
            df_neo_kind(() -> pxform(ctx2, "J2000", "DT_EMPTY", et)), "; frame of class 7: NeoSpice :", k7)
        @test k7 !== :none
        @test df_oracle_kind(:pxform_c, "J2000", "DT_CLASS7", et) !== :none
    end

    @testset "definitions follow furnsh, unload and pool edits" begin
        c = neo_context([Kernels.LSK, Kernels.PCK11, Kernels.DE440S])
        et = 2.0e8
        @test namfrm(c, "DT_GSE") == 0
        furnsh(c, DynFix.frame_kernel())
        a = pxform(c, "J2000", "DT_EULER", et)
        @test a == pxform(ctx, "J2000", "DT_EULER", et)
        # An edit of a definition gives a new snapshot with the new definition.
        pdpool(c, "FRAME_1850013_ANGLE_2_COEFFS", [0.2])
        SPICE.kclear(); foreach(SPICE.furnsh, set); SPICE.pdpool("FRAME_1850013_ANGLE_2_COEFFS", [0.2])
        b = pxform(c, "J2000", "DT_EULER", et)
        @test b != a
        @test maximum(abs.(b .- SPICE.pxform("J2000", "DT_EULER", et))) <= 5e-14
        unload(c, DynFix.frame_kernel())
        @test namfrm(c, "DT_EULER") == 0
        @test NeoSpice.kernelset(c).inner === nothing
        @test throws_kind(() -> pxform(c, "J2000", "DT_EULER", et), :unknownframe) ||
              throws_kind(() -> pxform(c, "J2000", "DT_EULER", et), :noframe)
    end

    @testset "switch frames over CK frames" begin
        set3 = [Kernels.LSK, Kernels.PCK11, Kernels.DE440S, CKFix.all_kernels()...,
            DynFix.switch_ck_kernel()]
        oracle_load(set3)
        c = neo_context(set3)
        tk = CKFix.ticks_tdb()
        t0 = sct2e(c, CKFix.SC_TDB, tk[1]) - 300.0
        t1 = sct2e(c, CKFix.SC_TDB, tk[end]) + 300.0
        # CSPICE keeps the segment of the last CK search of an instrument, and its answer
        # can then depend on the earlier calls (see `test_ck.jl`). A search far outside
        # the coverage clears that state before each query of the oracle.
        reset() = foreach(i -> (@assert SPICE.ckgp(i, 1e11, 0.0, "J2000") === nothing),
            (-77001, -77002, -77003, -77013))
        for name in ("DT_SWITCH_CK", "DT_SWITCH_CK_NOAV")
            dp = Float64[]; ds = Float64[]
            bases = name == "DT_SWITCH_CK" ? ("IAU_EARTH", "TESTSC_T3", "TESTSC_T2") :
                    ("ECLIPJ2000", "TESTSC_T1", "TESTSC_T3_NOAV")
            used = zeros(Int, 3)
            nerr = 0
            for et in random_epochs(rng, 500, t0, t1)
                reset()
                want = df_oracle_kind(:pxform_c, "J2000", name, et)
                got = df_neo_kind(() -> pxform(c, "J2000", name, et))
                @test got == want
                got === :none || (nerr += 1; continue)
                a = pxform(c, "J2000", name, et)
                reset()
                push!(dp, maximum(abs.(a .- SPICE.pxform("J2000", name, et))))
                for k in 3:-1:1
                    df_neo_kind(() -> pxform(c, "J2000", bases[k], et)) === :none || continue
                    a == pxform(c, "J2000", bases[k], et) && (used[k] += 1; break)
                end
                reset()
                wants = df_oracle_kind(:sxform_c, "J2000", name, et)
                gots = df_neo_kind(() -> sxform(c, "J2000", name, et))
                @test gots == wants
                gots === :none || continue
                reset()
                push!(ds, maximum(abs.(sxform(c, "J2000", name, et) .- SPICE.sxform("J2000", name, et))))
            end
            println("  ", name, ": bases used (lowest priority first) ", used, ", epochs with an error: ", nerr)
            s1 = report_quantiles("$name: pxform element", dp)
            s2 = report_quantiles("$name: sxform element", ds)
            @test s1.max <= 1e-12
            @test s2.max <= 1e-12
            @test count(>(0), used) >= 2
        end
    end
end
