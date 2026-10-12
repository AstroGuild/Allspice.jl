# Geometry finder, visibility in a field of view: `gftfov` against `gftfov_c`, `gfrfov`
# against `gfrfov_c`, and `gffove`.
#
# The seven instruments of the geometry tests (circle, ellipse, rectangle and polygon,
# given by corners and by angles), a camera frame that is fixed in J2000 (the Moon
# crosses the field of view in some hours) or fixed to the Earth (the Sun or a star
# crosses it once a day, so a window has several intervals), targets as ellipsoids and
# as points with the nine corrections, and rays with the three corrections of a ray
# (`NONE`, `S`, `XS`), in an inertial frame and in a rotating frame.
#
# Tolerance: the convergence tolerance of every search is the default of CSPICE,
# 1e-6 s. The number of intervals must be equal and each end point must agree within
# 2e-6 s, twice the convergence tolerance: the two libraries make the same bisection of
# a state that is true or false, so the brackets differ only if the state differs at a
# sample, which is possible only within rounding of an entry or of an exit.
#
# Behaviour specification: test/features/gf_fov.feature

isdefined(@__MODULE__, :GFO) || include("gforacle.jl")
isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")

group("gf fov: gftfov, gfrfov, gffove", Kernels.FULL_SET_FILES...) do
    ctx = geo_setup()
    ks = NeoSpice.kernelset(ctx)
    rng = MersenneTwister(SEED + 6800)
    day = 86400.0
    hour = 3600.0
    found = Dict{String,Int}()
    start_epoch() = -3.0e8 + 9.0e8 * rand(rng)
    "Make the camera frame an offset of `frame` (J2000 or IAU_EARTH) in the two libraries."
    camera_base!(frame) = both_pcpool(ctx, "TKFRAME_-999100_RELATIVE", frame)
    "A direction at a random small angle (up to `ang` radians) from `dir`."
    function near(dir, ang)
        z = normalize(collect(Float64, dir))
        h = abs(z[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
        x = normalize(cross(h, z))
        y = cross(z, x)
        a, p = ang * rand(rng), 2pi * rand(rng)
        return cos(a) * z + sin(a) * (cos(p) * x + sin(p) * y)
    end
    function check(label, neo, ora)
        w = neo()
        ok = gf_compare(label, w, ora())
        @test ok
        found[label] = get(found, label, 0) + wncard(w)
        return ok, w
    end

    @testset "gftfov: the Moon in a camera fixed in J2000" begin
        camera_base!("J2000")
        n = 0
        for inst in GEO_INSTRUMENTS, corr in GEO_CORRS, (tshape, tframe) in (("ELLIPSOID", "IAU_MOON"),
                ("POINT", " "))
            t0 = start_epoch()
            dir = spkpos(ctx, "MOON", t0, "J2000", corr, "EARTH")[1]
            point_camera!(ctx, near(dir, 0.012), 2pi * rand(rng))
            cn = rand(rng) < 0.7 ? Window([t0 - day, t0 + day]) : gf_random_window(rng, 3, t0 - day, t0 + day)
            step = rand(rng, (300.0, 600.0))
            ok, w = check("gftfov Moon, J2000 camera, " * lowercase(tshape),
                () -> gftfov(ctx, inst, "MOON", tshape, tframe, corr, "EARTH", step, cn),
                () -> GFO.gftfov(inst, "MOON", tshape, tframe, corr, "EARTH", step, collect(cn)))
            ok || @info "gftfov case" inst corr tshape t0 step cn
            @test all(t -> wnelmd(cn, t), w)
            n += 1
        end
        @test n == 7 * 9 * 2
        @test all(found[k] >= 40 for k in keys(found))
    end

    @testset "gftfov: the Sun in a camera fixed to the Earth" begin
        camera_base!("IAU_EARTH")
        for inst in GEO_INSTRUMENTS, corr in GEO_CORRS, (tshape, tframe) in (("ELLIPSOID", "IAU_SUN"),
                ("POINT", " "))
            t0 = start_epoch()
            dir = spkpos(ctx, "SUN", t0, "IAU_EARTH", corr, "EARTH")[1]
            point_camera!(ctx, near(dir, 0.01), 2pi * rand(rng))
            cn = Window([t0 - 3.5 * day, t0 + 3.5 * day])
            ok, w = check("gftfov Sun, Earth-fixed camera, " * lowercase(tshape),
                () -> gftfov(ctx, inst, "SUN", tshape, tframe, corr, "EARTH", 120.0, cn),
                () -> GFO.gftfov(inst, "SUN", tshape, tframe, corr, "EARTH", 120.0, collect(cn)))
            ok || @info "gftfov case" inst corr tshape t0 cn
        end
        @test found["gftfov Sun, Earth-fixed camera, ellipsoid"] >= 3 * 63
        @test found["gftfov Sun, Earth-fixed camera, point"] >= 2 * 63
    end

    @testset "gfrfov: a ray fixed in J2000, a camera fixed to the Earth" begin
        camera_base!("IAU_EARTH")
        for inst in GEO_INSTRUMENTS, corr in ("NONE", "S", "XS"), _ in 1:3
            t0 = start_epoch()
            ray = near([1.0, 0.0, 0.0], pi) .* (0.5 + 3.0 * rand(rng))
            dir = pxform(ctx, "J2000", "IAU_EARTH", t0) * ray
            point_camera!(ctx, near(dir, 0.01), 2pi * rand(rng))
            cn = Window([t0 - 2.5 * day, t0 + 2.5 * day])
            ok, w = check("gfrfov J2000 ray, Earth-fixed camera",
                () -> gfrfov(ctx, inst, ray, "J2000", corr, "EARTH", 60.0, cn),
                () -> GFO.gfrfov(inst, ray, "J2000", corr, "EARTH", 60.0, collect(cn)))
            ok || @info "gfrfov case" inst corr ray t0
        end
        @test found["gfrfov J2000 ray, Earth-fixed camera"] >= 3 * 63
    end

    @testset "gfrfov: a ray fixed to the Earth, a camera fixed in J2000" begin
        camera_base!("J2000")
        for inst in GEO_INSTRUMENTS, corr in ("NONE", "S", "XS"), _ in 1:2
            t0 = start_epoch()
            ray = near([0.0, 1.0, 0.0], 1.0)
            dir = pxform(ctx, "IAU_EARTH", "J2000", t0) * ray
            point_camera!(ctx, near(dir, 0.01), 2pi * rand(rng))
            cn = rand(rng) < 0.5 ? Window([t0 - 2.5 * day, t0 + 2.5 * day]) :
                 gf_random_window(rng, 4, t0 - 2.5 * day, t0 + 2.5 * day)
            ok, w = check("gfrfov Earth-fixed ray, J2000 camera",
                () -> gfrfov(ctx, inst, ray, "IAU_EARTH", corr, "EARTH", 60.0, cn),
                () -> GFO.gfrfov(inst, ray, "IAU_EARTH", corr, "EARTH", 60.0, collect(cn)))
            ok || @info "gfrfov case" inst corr ray t0 cn
            @test all(t -> wnelmd(cn, t), w)
        end
        @test found["gfrfov Earth-fixed ray, J2000 camera"] >= 42
        measured(@sprintf("  gf fov: intervals found: %s", join(["$k: $v" for (k, v) in
            sort!(collect(found); by=first)], "; ")))
    end

    @testset "gffove, forms and threads" begin
        camera_base!("IAU_EARTH")
        t0 = SPICE.str2et("2021 MAR 01")
        dir = spkpos(ctx, "SUN", t0, "IAU_EARTH", "LT+S", "EARTH")[1]
        point_camera!(ctx, dir, 0.3)
        cn = Window([t0 - 3.5 * day, t0 + 3.5 * day])
        a = ("GEO_RECT", "SUN", "ELLIPSOID", "IAU_SUN", "LT+S", "EARTH")
        base = gftfov(ctx, a..., 120.0, cn)
        @test wncard(base) >= 5
        ks2 = NeoSpice.kernelset(ctx)            # the kernel set after the changes of the pool
        @test gftfov(ks2, a..., 120.0, cn) == base
        @test gftfov(ks2, -999003, 10, "ellipsoid", namfrm(ctx, "IAU_SUN"), Aberration("LT+S"), 399, 120.0,
            cn) == base
        @test gftfov(ctx, a..., 120.0, cn; threads=true) == base
        many = gf_random_window(rng, 12, t0 - 3.5 * day, t0 + 3.5 * day)
        @test gftfov(ctx, a..., 120.0, many; threads=true) == gftfov(ctx, a..., 120.0, many)
        @test isempty(gftfov(ctx, a..., 120.0, Window()))
        @test gffove(ctx, "GEO_RECT", "ELLIPSOID", [0.0, 0.0, 0.0], "SUN", "IAU_SUN", "LT+S", "EARTH", 1.0e-6,
            120.0, gfrefn, cn) == base
        @test gffove(ks2, "GEO_RECT", "ELLIPSOID", [0.0, 0.0, 0.0], "SUN", "IAU_SUN", "LT+S", "EARTH", 1.0e-6,
            et -> 120.0, gfrefn, cn; threads=true) == base
        quarter(t1, t2, s1, s2) = t1 + 0.25 * (t2 - t1)
        other = gffove(ctx, "GEO_RECT", "ELLIPSOID", [0.0, 0.0, 0.0], "SUN", "IAU_SUN", "LT+S", "EARTH",
            1.0e-6, 120.0, quarter, cn)
        @test length(other) == length(base) && maximum(abs.(collect(other) .- collect(base))) <= 2.0e-6
        loose = gftfov(ctx, a..., 120.0, cn; tol=1.0e-1)
        @test length(loose) == length(base) && 1.0e-6 < maximum(abs.(collect(loose) .- collect(base))) <= 0.2
        gfstol(ctx, 1.0e-1)
        @test gftfov(ctx, a..., 120.0, cn) == loose
        gfstol(ctx, 1.0e-6)
        # a ray
        ray = pxform(ctx, "IAU_EARTH", "J2000", t0) * dir
        rbase = gfrfov(ctx, "GEO_RECT", ray, "J2000", "S", "EARTH", 60.0, cn)
        @test wncard(rbase) >= 5
        @test gfrfov(ks2, "GEO_RECT", SVector{3}(ray), "J2000", "S", "EARTH", 60.0, cn) == rbase
        @test gfrfov(ctx, "GEO_RECT", ray, "J2000", "S", "EARTH", 60.0, cn; threads=true) == rbase
        @test gffove(ctx, "GEO_RECT", "RAY", ray, " ", "J2000", "S", "EARTH", 1.0e-6, 60.0, gfrefn, cn) == rbase
        log = String[]
        gftfov(ctx, a..., 120.0, cn; report=(e, x...) -> e === :init && push!(log, x[2]))
        @test log == ["Target visibility search "]
        @test_throws InterruptException gfrfov(ctx, "GEO_RECT", ray, "J2000", "S", "EARTH", 60.0, cn;
            bail=() -> true)
        serial = [gftfov(ctx, inst, "SUN", "ELLIPSOID", "IAU_SUN", corr, "EARTH", 120.0, cn)
                  for inst in GEO_INSTRUMENTS, corr in ("NONE", "LT+S", "XCN")]
        tasks = [Threads.@spawn gftfov(ctx, inst, "SUN", "ELLIPSOID", "IAU_SUN", corr, "EARTH", 120.0, cn)
                 for inst in GEO_INSTRUMENTS, corr in ("NONE", "LT+S", "XCN")]
        @test all(fetch(t) == s for (t, s) in zip(tasks, serial))
    end

    @testset "errors, with the kinds of CSPICE" begin
        camera_base!("J2000")
        point_camera!(ctx, [0.0, 0.0, 1.0], 0.0)
        cn = Window([0.0, 2 * day])
        c = collect(cn)
        both(kind, neo, ora) = (@test gf_kind(neo) === kind; @test GFO.kind(ora) === kind)
        # inst, target, tshape, tframe, abcorr, obsrvr, step
        good = ("GEO_RECT", "MOON", "ELLIPSOID", "IAU_MOON", "LT+S", "EARTH", 600.0)
        with(p...) = (d = Dict(p); ntuple(k -> get(d, k, good[k]), length(good)))
        both(:none, () -> gftfov(ctx, good..., cn), () -> GFO.gftfov(good..., c))
        for (kind, args) in (
            (:idcodenotfound, with(1 => "NOSUCHINST")),
            (:idcodenotfound, with(2 => "NOSUCHBODY")),
            (:idcodenotfound, with(6 => "NOSUCHBODY")),
            (:bodiesnotdistinct, with(6 => "MOON")),
            (:invalidshape, with(3 => "SPHERE")),
            (:invalidoption, with(3 => "RAY")),
            (:invalidframe, with(4 => " ")),
            (:invalidframe, with(4 => "NOSUCHFRAME")),
            (:invalidframe, with(4 => "IAU_EARTH")),
            (:invalidoption, with(5 => "LT+X")),
            (:fovtoowide, with(1 => "GEO_FLAT")),
            (:invalidstep, with(7 => 0.0)),
            (:invalidstep, with(7 => -1.0)))
            both(kind, () -> gftfov(ctx, args..., cn), () -> GFO.gftfov(args..., c))
        end
        # the frame of a point is not read; case and blanks
        both(:none, () -> gftfov(ctx, with(3 => "POINT", 4 => "NOSUCHFRAME")..., cn),
            () -> GFO.gftfov(with(3 => "POINT", 4 => "NOSUCHFRAME")..., c))
        both(:none, () -> gftfov(ctx, "geo_rect", "moon", " ellipsoid", "iau_moon", "lt + s", "earth", 600.0, cn),
            () -> GFO.gftfov("geo_rect", "moon", " ellipsoid", "iau_moon", "lt + s", "earth", 600.0, c))
        # a body that has no field of view
        both(:framemissing, () -> gftfov(ctx, with(1 => "MARS")..., cn),
            () -> GFO.gftfov(with(1 => "MARS")..., c))
        # inst, raydir, rframe, abcorr, obsrvr, step
        rgood = ("GEO_RECT", [0.0, 0.0, 1.0], "J2000", "S", "EARTH", 600.0)
        rwith(p...) = (d = Dict(p); ntuple(k -> get(d, k, rgood[k]), length(rgood)))
        both(:none, () -> gfrfov(ctx, rgood..., cn), () -> GFO.gfrfov(rgood..., c))
        @test gfrfov(ctx, rgood..., cn) == cn        # the ray is the axis of the camera
        for (kind, args) in (
            (:idcodenotfound, rwith(1 => "NOSUCHINST")),
            (:zerovector, rwith(2 => [0.0, 0.0, 0.0])),
            (:invalidframe, rwith(3 => "NOSUCHFRAME")),
            (:invalidframe, rwith(3 => " ")),
            (:invalidoption, rwith(4 => "LT")),
            (:invalidoption, rwith(4 => "LT+S")),
            (:invalidoption, rwith(4 => "XCN+S")),
            (:invalidoption, rwith(4 => "S+X")),
            (:idcodenotfound, rwith(5 => "NOSUCHBODY")),
            (:fovtoowide, rwith(1 => "GEO_FLAT")),
            (:invalidstep, rwith(6 => 0.0)))
            both(kind, () -> gfrfov(ctx, args..., cn), () -> GFO.gfrfov(args..., c))
        end
        @test throws_kind(() -> gfrfov(ctx, rwith(2 => [0.0, 1.0])..., cn), :badarg)
        @test throws_kind(() -> gffove(ctx, "GEO_RECT", "RAY", [0.0, 0.0, 1.0], " ", "J2000", "S", "EARTH", 0.0,
                600.0, gfrefn, cn), :invalidtolerance)
    end

    @testset "search time against CSPICE (reported, not asserted)" begin
        camera_base!("IAU_EARTH")
        t0 = SPICE.str2et("2020 JAN 01")
        dir = spkpos(ctx, "SUN", t0 + 15 * day, "IAU_EARTH", "LT+S", "EARTH")[1]
        point_camera!(ctx, dir, 0.0)
        ray = pxform(ctx, "IAU_EARTH", "J2000", t0) * dir
        cn = Window([t0, t0 + 30 * day])
        c = collect(cn)
        cases2 = (
            ("gftfov Sun as an ellipsoid, rectangle, 30 days, step 120 s, LT+S",
                () -> gftfov(ctx, "GEO_RECT", "SUN", "ELLIPSOID", "IAU_SUN", "LT+S", "EARTH", 120.0, cn),
                () -> GFO.gftfov("GEO_RECT", "SUN", "ELLIPSOID", "IAU_SUN", "LT+S", "EARTH", 120.0, c)),
            ("gftfov Sun as an ellipsoid, ellipse, 30 days, step 120 s, CN",
                () -> gftfov(ctx, "GEO_ELLIPSE", "SUN", "ELLIPSOID", "IAU_SUN", "CN", "EARTH", 120.0, cn),
                () -> GFO.gftfov("GEO_ELLIPSE", "SUN", "ELLIPSOID", "IAU_SUN", "CN", "EARTH", 120.0, c)),
            ("gftfov Sun as a point, polygon, 30 days, step 120 s, NONE",
                () -> gftfov(ctx, "GEO_POLY", "SUN", "POINT", " ", "NONE", "EARTH", 120.0, cn),
                () -> GFO.gftfov("GEO_POLY", "SUN", "POINT", " ", "NONE", "EARTH", 120.0, c)),
            ("gfrfov J2000 ray, circle, 30 days, step 60 s, S",
                () -> gfrfov(ctx, "GEO_CIRCLE", ray, "J2000", "S", "EARTH", 60.0, cn),
                () -> GFO.gfrfov("GEO_CIRCLE", ray, "J2000", "S", "EARTH", 60.0, c)),
        )
        for (name, neo, ora) in cases2
            @test gf_compare("timing", neo(), ora())
            tn, to = gf_time(neo), gf_time(ora)
            measured(@sprintf("  %-66s NeoSpice %8.2f ms  CSPICE %8.2f ms  ratio %.2f", name,
                1e3 * tn, 1e3 * to, to / tn))
        end
    end

    camera_base!("J2000")
    point_camera!(ctx, [0.0, 0.0, 1.0], 0.0)
    gf_report("gftfov")
    gf_report("gfrfov")
end
