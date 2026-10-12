# Geometry finder, coordinates of a surface point: `gfsubc` against `gfsubc_c` and
# `gfsntc` against `gfsntc_c`, and the matching vector definitions of `gfevnt`.
#
# Sub-observer point: three target, frame and observer triples, the two methods, every
# coordinate system and coordinate, the seven relations and the two absolute extrema
# with an adjust value, a random aberration correction among the nine of `spkezr`, and
# confinement windows of one and of several intervals.
#
# Surface intercept: a ray that follows the target (a two-vector dynamic frame, so the
# intercept exists on the whole window), a ray that meets the target for a part of each
# month only (the existence window has several intervals), and a ray that is fixed in
# J2000 and that the target crosses in about an hour (the frame of the ray is not
# centered on the observer, the other branch of `zzgfssin`).
#
# Coordinates that are constant to rounding are not compared: the altitude of a surface
# point (zero) and the radius of a point of the Moon (a sphere). A search on such a
# coordinate follows the rounding noise in each library, so the two windows have no
# relation. The test checks only that the two libraries accept the search.
#
# A ray that crosses the limb: with a converged correction (`CN`, `XCN` and the `+S`
# forms) the flag "intercept found" of `sincpt` is not monotone within about 0.3 m of
# the limb, in CSPICE and in NeoSpice, because the light-time loop stops after a
# number of iterations that depends on rounding. An end point of the existence window
# then differs by some milliseconds between the two libraries. The searches on such a
# ray use the other corrections, and the existence window is compared for the nine
# corrections with a tolerance of 5e-2 s for the converged ones. A ray that stays near
# the limb amplifies the rounding of the ray in the intercept point; the end points of
# the searches on that ray agree within 2e-5 s.
#
# Tolerance: the convergence tolerance of every search is the default of CSPICE,
# 1e-6 s. The number of intervals must be equal and each end point must agree within
# 2e-6 s, twice the convergence tolerance (the bracket of an event is not longer than
# the tolerance in each library, and the two brackets differ by one bisection step at
# most when a state differs by rounding at a sample). A root or an extremum where the
# coordinate is almost stationary is accepted with a larger difference, up to 1 s, if
# the coordinate has the same value at the two epochs within 1e-12 of its scale (the
# rounding noise of a coordinate of a surface point is about 1e-13 of its scale).
#
# Behaviour specification: test/features/gf_surface.feature

isdefined(@__MODULE__, :GFO) || include("gforacle.jl")
isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")

const GFS_TRACK_PATH = Ref("")

"""
    gfs_track_kernel() -> String

Write (once per process) a frame kernel with the dynamic frame `GF_TRACK`, centered on
the Moon: +X is the geometric position of the Earth as seen from the Moon and +Y is in
the half plane of the velocity of the Earth relative to the Moon. A ray that is fixed
in this frame near +X follows the Earth. Return the path of the file.
"""
function gfs_track_kernel()
    isfile(GFS_TRACK_PATH[]) && return GFS_TRACK_PATH[]
    path = joinpath(mktempdir(), "gf_track.tf")
    write(path, """
    KPL/FK

    Dynamic frame of the geometry finder tests (written by the tests).

    \\begindata

       FRAME_GF_TRACK                 = -998100
       FRAME_-998100_NAME             = 'GF_TRACK'
       FRAME_-998100_CLASS            = 5
       FRAME_-998100_CLASS_ID         = -998100
       FRAME_-998100_CENTER           = 301
       FRAME_-998100_RELATIVE         = 'J2000'
       FRAME_-998100_DEF_STYLE        = 'PARAMETERIZED'
       FRAME_-998100_FAMILY           = 'TWO-VECTOR'
       FRAME_-998100_PRI_AXIS         = 'X'
       FRAME_-998100_PRI_VECTOR_DEF   = 'OBSERVER_TARGET_POSITION'
       FRAME_-998100_PRI_OBSERVER     = 'MOON'
       FRAME_-998100_PRI_TARGET       = 'EARTH'
       FRAME_-998100_PRI_ABCORR       = 'NONE'
       FRAME_-998100_SEC_AXIS         = 'Y'
       FRAME_-998100_SEC_VECTOR_DEF   = 'OBSERVER_TARGET_VELOCITY'
       FRAME_-998100_SEC_OBSERVER     = 'MOON'
       FRAME_-998100_SEC_TARGET       = 'EARTH'
       FRAME_-998100_SEC_ABCORR       = 'NONE'
       FRAME_-998100_SEC_FRAME        = 'J2000'

    \\begintext
    """)
    GFS_TRACK_PATH[] = path
    @assert isfile(path)
    return path
end

"""
    gfs_reference(rng, q, w, rel, withadjust) -> (refval, adjust)

Reference value and adjust value for a search of the coordinate `q` on the window `w`
on which its vector exists: as `gf_reference`, with the samples taken inside the
intervals of `w`. Precondition: `w` is not empty.
"""
function gfs_reference(rng::AbstractRNG, q, w::Window, rel::AbstractString, withadjust::Bool)
    @assert wncard(w) >= 1
    vals = Float64[]
    for i in 1:wncard(w), x in range(0.02, 0.98; length=24)
        push!(vals, NeoSpice._gf_value(q, w[2i-1] + x * (w[2i] - w[2i-1])))
    end
    lo, hi = quantile_of(vals, 0.25), quantile_of(vals, 0.75)
    refval = rel in ("<", "=", ">") ? lo + (hi - lo) * rand(rng) : 0.0
    adjust = withadjust ? (0.1 + 0.2 * rand(rng)) * (maximum(vals) - minimum(vals)) : 0.0
    return refval, adjust
end

"Largest absolute value of the coordinate `q` at 20 epochs of the first interval of `w`."
gfs_scale(q, w::Window) = maximum(abs(NeoSpice._gf_value(q, w[1] + x * (w[2] - w[1])))
                                  for x in range(0.0, 1.0; length=20))

"Tolerance of the value of the coordinate `q` in `gf_compare`: 1e-12 of its scale on `w`."
gfs_vtol(q, w::Window) = 1.0e-12 * max(1.0, gfs_scale(q, w))

"""
    gfs_constant(target, sys, crd) -> Bool

True if the coordinate of a point of the surface of `target` is constant to rounding:
an altitude, or the distance from the center of the Moon (whose radii are equal).
"""
gfs_constant(target, sys, crd) = crd == "ALTITUDE" ||
    (target == "MOON" && sys != "CYLINDRICAL" && crd in ("RADIUS", "RANGE"))

"The corrections whose light-time loop has a fixed number of passes."
const GFS_FIXED_CORRS = ("NONE", "LT", "LT+S", "XLT", "XLT+S")

group("gf surface: gfsubc, gfsntc", Kernels.FULL_SET_FILES...) do
    paths = [geo_kernels()..., gfs_track_kernel()]
    oracle_load(paths)
    ctx = neo_context(paths)
    ks = NeoSpice.kernelset(ctx)
    rng = MersenneTwister(SEED + 6500)
    day = 86400.0
    hour = 3600.0
    methods = ("Near point: ellipsoid", "Intercept: ellipsoid")
    start_epoch() = -3.0e8 + 9.0e8 * rand(rng)
    function confinement(span)
        t0 = start_epoch()
        return rand(rng) < 0.7 ? Window([t0, t0 + span]) : gf_random_window(rng, 3, t0, t0 + span)
    end
    subq(targ, fix, meth, corr, obs, sys, crd) = NeoSpice._gf_surface_coordinate(ks, false, true,
        targ, fix, meth, corr, obs, 0, [0.0, 0.0, 0.0], sys, crd)
    sntq(targ, fix, corr, obs, dref, dvec, sys, crd) = NeoSpice._gf_surface_coordinate(ks, true,
        true, targ, fix, "Ellipsoid", corr, obs, dref, dvec, sys, crd)

    @testset "gfsubc: every system, coordinate and relation against CSPICE" begin
        # (target, frame, observer, steps, span)
        cases = (("EARTH", "IAU_EARTH", "MOON", (2 * hour, 3 * hour), 4 * day),
            ("MOON", "IAU_MOON", "EARTH", (day, day / 2), 90 * day),
            ("EARTH", "IAU_EARTH", "SUN", (3 * hour,), 3 * day))
        n = 0
        nconst = 0
        for (targ, fix, obs, steps, span) in cases, (sys, crds) in GF_SYSTEMS, crd in crds,
            (rel, wadj) in GF_VARIANTS

            corr = rand(rng, GF_ABCORRS)
            meth = rand(rng, methods)
            step = rand(rng, steps)
            cn = confinement(span)
            q = subq(targ, fix, meth, corr, obs, sys, crd)
            refval, adjust = gf_reference(rng, q, cn, rel, wadj)
            w = gfsubc(ctx, targ, fix, meth, corr, obs, sys, crd, rel, refval, adjust, step, cn)
            ref = GFO.gfsubc(targ, fix, meth, corr, obs, sys, crd, rel, refval, adjust, step, collect(cn))
            @test all(t -> wnelmd(cn, t), w)
            if gfs_constant(targ, sys, crd)
                nconst += 1
                continue
            end
            angular = (sys, crd) in GF_ANGULAR
            label = (angular ? "gfsubc longitude/RA " : "gfsubc other ") * gf_family(rel, adjust)
            ok = gf_compare(label, w, ref; value=t -> NeoSpice._gf_value(q, t),
                vtol=gfs_vtol(q, cn))
            ok || @info "gfsubc case" targ fix obs meth corr sys crd rel refval adjust step cn
            @test ok
            n += 1
        end
        @test n + nconst == 3 * 21 * 9 && nconst == (2 * 3 + 3) * 9
        measured(@sprintf("  gfsubc: %d searches compared, %d on a constant coordinate not compared", n, nconst))
    end

    @testset "gfsubc: a root is a root, without CSPICE" begin
        # The coordinate crosses the reference value at each end point of an equality
        # search (radius of the sub-Moon point on the Earth, which has no branch cut).
        for meth in methods, corr in ("NONE", "LT+S", "XCN")
            cn = Window([0.0, 60 * day])
            q = subq("EARTH", "IAU_EARTH", meth, corr, "MOON", "LATITUDINAL", "LATITUDE")
            refval = 0.2
            w = gfsubc(ctx, "EARTH", "IAU_EARTH", meth, corr, "MOON", "LATITUDINAL", "LATITUDE", "=",
                refval, 0.0, 6 * hour, cn)
            @test wncard(w) >= 2
            for i in 1:wncard(w)
                t = w[2i-1]
                @test w[2i] == t
                @test (NeoSpice._gf_value(q, t - 1.0) - refval) * (NeoSpice._gf_value(q, t + 1.0) - refval) < 0
            end
        end
    end

    @testset "gfsntc: a ray that follows the target" begin
        # The intercept exists on the whole window; every system and coordinate.
        n = 0
        nconst = 0
        for (sys, crds) in GF_SYSTEMS, crd in crds, (rel, wadj) in GF_VARIANTS
            corr = rand(rng, GF_ABCORRS)
            dvec = [1.0, 0.012 * (rand(rng) - 0.5), 0.012 * (rand(rng) - 0.5)]
            step = rand(rng, (2 * hour, 3 * hour))
            t0 = start_epoch()
            cn = Window([t0, t0 + 3 * day])
            q = sntq("EARTH", "IAU_EARTH", corr, "MOON", "GF_TRACK", dvec, sys, crd)
            exist = NeoSpice._gf_vecexist(q.vec, cn, NeoSpice._gf_run(step, 1.0e-6))
            @test collect(exist) ≈ [cn[1] + 1.0 + 1.0e-6, cn[2] - 1.0 - 1.0e-6]
            refval, adjust = gfs_reference(rng, q, exist, rel, wadj)
            w = gfsntc(ctx, "EARTH", "IAU_EARTH", "Ellipsoid", corr, "MOON", "GF_TRACK", dvec, sys,
                crd, rel, refval, adjust, step, cn)
            ref = GFO.gfsntc("EARTH", "IAU_EARTH", "Ellipsoid", corr, "MOON", "GF_TRACK", dvec, sys,
                crd, rel, refval, adjust, step, collect(cn))
            @test all(t -> wnelmd(exist, t), w)
            if gfs_constant("EARTH", sys, crd)
                nconst += 1
                continue
            end
            angular = (sys, crd) in GF_ANGULAR
            label = (angular ? "gfsntc tracking ray, longitude/RA " : "gfsntc tracking ray, other ") *
                    gf_family(rel, adjust)
            ok = gf_compare(label, w, ref; value=t -> NeoSpice._gf_value(q, t),
                vtol=gfs_vtol(q, exist))
            ok || @info "gfsntc case" corr dvec sys crd rel refval adjust step cn
            @test ok
            n += 1
        end
        @test n + nconst == 21 * 9 && nconst == 2 * 9
    end

    @testset "gfsntc: a ray that meets the target for a part of each month" begin
        # The ray passes 0.0166 rad from the center of the Earth, so it meets the Earth
        # only while the Moon is nearer than about 384,000 km. The intercept is always
        # near the limb, where a rounding error of the ray is amplified in the point, so
        # the two libraries find a root at epochs that differ by more than the
        # convergence tolerance: the tolerance of an end point is 2e-5 s here.
        dvec = [1.0, 0.0166, 0.0]
        nexist = 0
        for (sys, crd) in (("LATITUDINAL", "LATITUDE"), ("RECTANGULAR", "Z"), ("RECTANGULAR", "X"),
                ("GEODETIC", "LONGITUDE")), (rel, wadj) in GF_VARIANTS
            corr = rand(rng, GFS_FIXED_CORRS)
            t0 = start_epoch()
            cn = Window([t0, t0 + 60 * day])
            step = 3 * hour
            q = sntq("EARTH", "IAU_EARTH", corr, "MOON", "GF_TRACK", dvec, sys, crd)
            exist = NeoSpice._gf_vecexist(q.vec, cn, NeoSpice._gf_run(step, 1.0e-6))
            @test wncard(exist) >= 2
            nexist = max(nexist, wncard(exist))
            refval, adjust = gfs_reference(rng, q, exist, rel, wadj)
            w = gfsntc(ctx, "EARTH", "IAU_EARTH", "ELLIPSOID", corr, "MOON", "GF_TRACK", dvec, sys,
                crd, rel, refval, adjust, step, cn)
            ref = GFO.gfsntc("EARTH", "IAU_EARTH", "ELLIPSOID", corr, "MOON", "GF_TRACK", dvec, sys,
                crd, rel, refval, adjust, step, collect(cn))
            ok = gf_compare("gfsntc grazing ray " * gf_family(rel, adjust), w, ref; atol=2.0e-5,
                value=t -> NeoSpice._gf_value(q, t), vtol=gfs_vtol(q, exist))
            ok || @info "gfsntc grazing case" corr sys crd rel refval adjust cn
            @test ok
            @test all(t -> wnelmd(exist, t), w)
        end
        measured(@sprintf("  gfsntc grazing ray: largest number of intervals of the existence window %d", nexist))
    end

    @testset "gfsntc: a ray fixed in J2000 that the target crosses" begin
        n = 0
        for (sys, crds) in GF_SYSTEMS, crd in crds, rel in ("<", "=", "LOCMAX", "ABSMIN")
            gfs_constant("MOON", sys, crd) && continue
            corr = rand(rng, GFS_FIXED_CORRS)
            t0 = start_epoch()
            # the direction of the Moon at t0, as the search will see it
            dvec = collect(spkpos(ctx, "MOON", t0, "J2000", corr, "EARTH")[1])
            cn = Window([t0 - 6 * hour, t0 + 6 * hour])
            step = 300.0
            q = sntq("MOON", "IAU_MOON", corr, "EARTH", "J2000", dvec, sys, crd)
            exist = NeoSpice._gf_vecexist(q.vec, cn, NeoSpice._gf_run(step, 1.0e-6))
            @test wncard(exist) == 1 && cn[1] < exist[1] < t0 < exist[2] < cn[2]
            refval, adjust = gfs_reference(rng, q, exist, rel, false)
            w = gfsntc(ctx, "MOON", "IAU_MOON", "ellipsoid", corr, "EARTH", "J2000", dvec, sys, crd,
                rel, refval, adjust, step, cn)
            ref = GFO.gfsntc("MOON", "IAU_MOON", "ellipsoid", corr, "EARTH", "J2000", dvec, sys, crd,
                rel, refval, adjust, step, collect(cn))
            ok = gf_compare("gfsntc crossing ray " * gf_family(rel, adjust), w, ref;
                value=t -> NeoSpice._gf_value(q, t), vtol=gfs_vtol(q, exist))
            ok || @info "gfsntc crossing case" corr sys crd rel refval t0
            @test ok
            n += 1
        end
        @test n == (21 - 5) * 4
    end

    @testset "gfsntc: the existence window, the nine corrections" begin
        # A coordinate that is less than the reference value wherever the intercept
        # exists: the result is the existence window, contracted as in CSPICE.
        worst = Dict(false => 0.0, true => 0.0)
        for corr in GF_ABCORRS, k in 1:4
            converged = !(corr in GFS_FIXED_CORRS)
            atol = converged ? 5.0e-2 : 2.0e-6
            t0 = start_epoch()
            if isodd(k)
                a = ("EARTH", "IAU_EARTH", "Ellipsoid", corr, "MOON", "GF_TRACK", [1.0, 0.0166, 0.0],
                    "RECTANGULAR", "X", "<", 1.0e9, 0.0, 3 * hour)
                cn = Window([t0, t0 + 60 * day])
            else
                dvec = collect(spkpos(ctx, "MOON", t0, "J2000", corr, "EARTH")[1])
                a = ("MOON", "IAU_MOON", "Ellipsoid", corr, "EARTH", "J2000", dvec, "RECTANGULAR", "X",
                    "<", 1.0e9, 0.0, 300.0)
                cn = Window([t0 - 6 * hour, t0 + 6 * hour])
            end
            w = gfsntc(ctx, a..., cn)
            ref = GFO.gfsntc(a..., collect(cn))
            @test wncard(w) >= 1
            ok = gf_compare("gfsntc existence window, " * (converged ? "converged" : "fixed passes"),
                w, ref; atol)
            ok || @info "gfsntc existence case" corr k t0
            @test ok
            length(w) == length(ref) && (worst[converged] = max(worst[converged],
                maximum(abs.(collect(w) .- ref))))
        end
        measured(@sprintf("  gfsntc existence window: largest end-point difference %.3e s (fixed passes), %.3e s (converged)",
            worst[false], worst[true]))
    end

    @testset "gfevnt: the surface vector definitions" begin
        cn = Window([0.0, 30 * day])
        names = ["TARGET", "OBSERVER", "ABCORR", "COORDINATE SYSTEM", "COORDINATE", "REFERENCE FRAME",
            "VECTOR DEFINITION", "METHOD", "DVEC", "DREF"]
        ev(vals, dvec; op="=", refval=0.1) = gfevnt(ctx, 6 * hour, gfrefn, "COORDINATE", names, vals,
            dvec, Int[], Bool[], op, refval, 1.0e-6, 0.0, cn)
        sub = gfsubc(ctx, "EARTH", "IAU_EARTH", "Near point: ellipsoid", "LT+S", "MOON", "GEODETIC",
            "LATITUDE", "=", 0.1, 0.0, 6 * hour, cn)
        @test !isempty(sub)
        @test ev(["EARTH", "MOON", "LT+S", "GEODETIC", "LATITUDE", "IAU_EARTH", "SUB-OBSERVER POINT",
                "Near point: ellipsoid", " ", " "], Float64[]) == sub
        @test ev([" earth", "moon", "lt+s", "geodetic", "latitude", "iau_earth", "sub-observer point",
                "near point: ellipsoid", " ", " "], [0.0, 0.0, 0.0]) == sub
        dvec = [1.0, 0.004, -0.003]
        snt = gfsntc(ctx, "EARTH", "IAU_EARTH", "Ellipsoid", "CN", "MOON", "GF_TRACK", dvec,
            "LATITUDINAL", "LONGITUDE", "<", 1.0, 0.0, 6 * hour, cn)
        @test !isempty(snt)
        @test ev(["EARTH", "MOON", "CN", "LATITUDINAL", "LONGITUDE", "IAU_EARTH", "SURFACE INTERCEPT POINT",
                "Ellipsoid", " ", "GF_TRACK"], dvec; op="<", refval=1.0) == snt
        @test throws_kind(() -> ev(["EARTH", "MOON", "CN", "LATITUDINAL", "LONGITUDE", "IAU_EARTH",
                "SURFACE INTERCEPT POINT", "Ellipsoid", " ", "GF_TRACK"], [1.0, 0.0]), :invalidcount)
        # CSPICE checks the center of the frame only for the upper-case spelling of the
        # vector definition; with another spelling the geometry call makes the check.
        @test throws_kind(() -> ev(["EARTH", "MOON", "NONE", "LATITUDINAL", "LATITUDE", "IAU_MOON",
                "SUB-OBSERVER POINT", "Near point: ellipsoid", " ", " "], Float64[]), :invalidframe)
        @test throws_kind(() -> ev(["EARTH", "MOON", "NONE", "LATITUDINAL", "LATITUDE", "IAU_MOON",
                "sub-observer point", "Near point: ellipsoid", " ", " "], Float64[]), :invalidframe)
    end

    @testset "forms, keywords and threads" begin
        cn = Window([0.0, 4 * day])
        many = gf_random_window(rng, 12, 0.0, 12 * day)
        a = ("EARTH", "IAU_EARTH", "Intercept: ellipsoid", "LT+S", "MOON", "LATITUDINAL", "LONGITUDE",
            "=", 1.0, 0.0, 3 * hour)
        base = gfsubc(ctx, a..., cn)
        @test wncard(base) >= 3
        @test gfsubc(ks, a..., cn) == base
        @test gfsubc(ks, 399, namfrm(ctx, "IAU_EARTH"), a[3], Aberration("LT+S"), 301, a[6:end]..., cn) == base
        @test gfsubc(ctx, a..., cn; threads=true) == base
        @test gfsubc(ctx, a..., many; threads=true) == gfsubc(ctx, a..., many)
        @test isempty(gfsubc(ctx, a..., Window()))
        # a looser tolerance moves the roots by less than that tolerance
        loose = gfsubc(ctx, a..., cn; tol=1.0e-2)
        @test length(loose) == length(base) && maximum(abs.(collect(loose) .- collect(base))) < 2.0e-2
        # the tolerance of the context
        gfstol(ctx, 1.0e-2)
        @test gfsubc(ctx, a..., cn) == loose
        gfstol(ctx, 1.0e-6)
        dvec = SVector(1.0, 0.0166, 0.0)
        b = ("EARTH", "IAU_EARTH", "Ellipsoid", "NONE", "MOON", "GF_TRACK", dvec, "RECTANGULAR", "Z",
            ">", 0.0, 0.0, 3 * hour)
        long = Window([0.0, 60 * day])
        sbase = gfsntc(ctx, b..., long)
        @test wncard(sbase) >= 2
        @test gfsntc(ks, b..., long) == sbase
        @test gfsntc(ctx, b..., long; threads=true) == sbase
        @test isempty(gfsntc(ctx, b..., Window()))
        # the callbacks: the existence pass is reported, and the search can be stopped
        log = String[]
        gfsntc(ctx, b..., long; report=(e, x...) -> e === :init && push!(log, x[2]))
        @test log[1] == "Intercept existence pass 1 of 1" && length(log) >= 2
        @test_throws InterruptException gfsntc(ctx, b..., cn; bail=() -> true)
        @test_throws InterruptException gfsubc(ctx, a..., cn; bail=() -> true)
    end

    @testset "errors, with the kinds of CSPICE" begin
        cn = Window([0.0, 2 * day])
        c = collect(cn)
        both(kind, neo, ora) = (@test gf_kind(neo) === kind; @test GFO.kind(ora) === kind)
        NP = "Near point: ellipsoid"
        # gfsubc: the arguments in the order target, fixref, method, abcorr, obsrvr,
        # crdsys, coord, relate, refval, adjust, step
        good = ("EARTH", "IAU_EARTH", NP, "NONE", "MOON", "LATITUDINAL", "LATITUDE", "=", 0.1, 0.0, 3 * hour)
        with(i, v) = ntuple(k -> k == i ? v : good[k], length(good))
        both(:none, () -> gfsubc(ctx, good..., cn), () -> GFO.gfsubc(good..., c))
        for (kind, args) in (
            (:idcodenotfound, with(1, "NOSUCHBODY")),
            (:idcodenotfound, with(5, "NOSUCHBODY")),
            (:bodiesnotdistinct, with(5, "399")),
            (:noframe, with(2, "NOSUCHFRAME")),
            (:invalidframe, with(2, "IAU_MOON")),
            (:notsupported, with(3, "NEAR POINT/ELLIPSOID")),
            (:notsupported, with(3, "Ellipsoid")),
            (:notsupported, with(3, "NADIR/DSK/UNPRIORITIZED")),
            (:invalidoption, with(4, "LT+X")),
            (:notsupported, with(6, "POLAR")),
            (:notsupported, with(7, "ALTITUDE")),
            (:notrecognized, with(8, "==")),
            (:valueoutofrange, with(10, -1.0)),
            (:invalidvalue, with(10, 1.0)),
            (:invalidstep, with(11, 0.0)),
            (:invalidstep, with(11, -1.0)))
            both(kind, () -> gfsubc(ctx, args..., cn), () -> GFO.gfsubc(args..., c))
        end
        # the method is read at the first evaluation: an empty window gives no error
        both(:none, () -> gfsubc(ctx, with(3, "NEAR POINT/ELLIPSOID")..., Window()),
            () -> GFO.gfsubc(with(3, "NEAR POINT/ELLIPSOID")..., Float64[]))
        # case and blanks of the method, the system and the coordinate
        both(:none, () -> gfsubc(ctx, "earth", " iau_earth", " near point:ellipsoid", "lt + s", "moon",
                "lati tudinal", " latitude ", " =", 0.1, 0.0, 3 * hour, cn),
            () -> GFO.gfsubc("earth", " iau_earth", " near point:ellipsoid", "lt + s", "moon",
                "lati tudinal", " latitude ", " =", 0.1, 0.0, 3 * hour, c))
        # a frame that is not centered on the target, and unequal equatorial radii
        both(:invalidframe, () -> gfsubc(ctx, "SUN", "J2000", NP, "NONE", "EARTH", "LATITUDINAL",
                "LATITUDE", "=", 0.1, 0.0, day, cn),
            () -> GFO.gfsubc("SUN", "J2000", NP, "NONE", "EARTH", "LATITUDINAL", "LATITUDE", "=", 0.1,
                0.0, day, c))
        radii = gdpool(ctx, "BODY399_RADII")
        both_pdpool(ctx, "BODY399_RADII", [6378.0, 6379.0, 6356.0])
        both(:notsupported, () -> gfsubc(ctx, with(6, "GEODETIC")..., cn),
            () -> GFO.gfsubc(with(6, "GEODETIC")..., c))
        both(:none, () -> gfsubc(ctx, good..., cn), () -> GFO.gfsubc(good..., c))
        both_pdpool(ctx, "BODY399_RADII", radii)

        # gfsntc: target, fixref, method, abcorr, obsrvr, dref, dvec, crdsys, coord,
        # relate, refval, adjust, step
        sgood = ("EARTH", "IAU_EARTH", "Ellipsoid", "NONE", "MOON", "GF_TRACK", [1.0, 0.004, 0.003],
            "LATITUDINAL", "LATITUDE", "=", 0.1, 0.0, 3 * hour)
        swith(i, v) = ntuple(k -> k == i ? v : sgood[k], length(sgood))
        both(:none, () -> gfsntc(ctx, sgood..., cn), () -> GFO.gfsntc(sgood..., c))
        for (kind, args) in (
            (:idcodenotfound, swith(1, "NOSUCHBODY")),
            (:idcodenotfound, swith(5, "NOSUCHBODY")),
            (:bodiesnotdistinct, swith(5, "EARTH")),
            (:noframe, swith(2, "NOSUCHFRAME")),
            (:invalidframe, swith(2, "IAU_MOON")),
            (:notsupported, swith(3, "DSK/UNPRIORITIZED")),
            (:notsupported, swith(3, "Intercept: ellipsoid")),
            (:invalidoption, swith(4, "S")),
            (:noframe, swith(6, "NOSUCHFRAME")),
            (:zerovector, swith(7, [0.0, 0.0, 0.0])),
            (:notsupported, swith(8, "POLAR")),
            (:notsupported, swith(9, "RANGE")),
            (:notrecognized, swith(10, "LOCALMAX")),
            (:valueoutofrange, swith(12, -1.0)),
            (:invalidvalue, swith(12, 1.0)),
            (:invalidstep, swith(13, 0.0)))
            both(kind, () -> gfsntc(ctx, args..., cn), () -> GFO.gfsntc(args..., c))
        end
        @test throws_kind(() -> gfsntc(ctx, swith(7, [1.0, 0.0])..., cn), :badarg)
        # a ray that never meets the target: an empty window, no error
        both(:none, () -> gfsntc(ctx, swith(7, [-1.0, 0.0, 0.0])..., cn),
            () -> GFO.gfsntc(swith(7, [-1.0, 0.0, 0.0])..., c))
        @test isempty(gfsntc(ctx, swith(7, [-1.0, 0.0, 0.0])..., cn))
    end

    @testset "search time against CSPICE (reported, not asserted)" begin
        t0 = SPICE.str2et("2020 JAN 01")
        cases2 = (
            ("gfsubc sub-Moon latitude on the Earth, 1 year, '=', LT+S", 365 * day,
                c -> gfsubc(ctx, "EARTH", "IAU_EARTH", methods[1], "LT+S", "MOON", "GEODETIC", "LATITUDE",
                    "=", 0.2, 0.0, 6 * hour, c),
                c -> GFO.gfsubc("EARTH", "IAU_EARTH", methods[1], "LT+S", "MOON", "GEODETIC", "LATITUDE",
                    "=", 0.2, 0.0, 6 * hour, c)),
            ("gfsubc sub-Moon longitude on the Earth, 30 days, '<', NONE", 30 * day,
                c -> gfsubc(ctx, "EARTH", "IAU_EARTH", methods[2], "NONE", "MOON", "LATITUDINAL", "LONGITUDE",
                    "<", 1.0, 0.0, 3 * hour, c),
                c -> GFO.gfsubc("EARTH", "IAU_EARTH", methods[2], "NONE", "MOON", "LATITUDINAL", "LONGITUDE",
                    "<", 1.0, 0.0, 3 * hour, c)),
            ("gfsntc tracking ray, latitude, 30 days, LOCMAX, CN+S", 30 * day,
                c -> gfsntc(ctx, "EARTH", "IAU_EARTH", "Ellipsoid", "CN+S", "MOON", "GF_TRACK",
                    [1.0, 0.004, 0.003], "LATITUDINAL", "LATITUDE", "LOCMAX", 0.0, 0.0, 3 * hour, c),
                c -> GFO.gfsntc("EARTH", "IAU_EARTH", "Ellipsoid", "CN+S", "MOON", "GF_TRACK",
                    [1.0, 0.004, 0.003], "LATITUDINAL", "LATITUDE", "LOCMAX", 0.0, 0.0, 3 * hour, c)),
        )
        for (name, span, neo, ora) in cases2
            cn = Window([t0, t0 + span])
            c = collect(cn)
            @test gf_compare("timing", neo(cn), ora(c))
            tn, to = gf_time(() -> neo(cn)), gf_time(() -> ora(c))
            measured(@sprintf("  %-62s NeoSpice %8.2f ms  CSPICE %8.2f ms  ratio %.2f", name,
                1e3 * tn, 1e3 * to, to / tn))
        end
    end

    gf_report("gfsubc")
    gf_report("gfsntc")
end
