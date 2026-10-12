# Geometry group: instrument fields of view (`getfov`, `getfvn`, `fovtrg`, `fovray`)
# against CSPICE.
#
# Feature file: features/geometry_fov.feature
#
# The instruments are those of `geo_instrument_kernel`: a circle, an ellipse, a
# rectangle and a polygon given by corners, a circle, an ellipse and a rectangle given
# by angles, a wide polygon whose axis comes from the convex hull search, and a polygon
# that is too wide. The frame of the instruments is a fixed offset from J2000; each case
# sets its matrix in the two kernel pools, so that the target is near the edge of the
# field of view.
#
# The outputs of `fovtrg` and `fovray` are flags, so the comparison is exact.
#
# The default run has 500 random cases of `fovtrg` and `fovray`; the exhaustive run
# (`NEOSPICE_TEST_FULL=1`) has 2000.

isdefined(@__MODULE__, :GEOLIB) || include("geomfixtures.jl")

const FOV_SCRATCH = -999050
const FOV_VARS = ["FOV_FRAME", "FOV_SHAPE", "BORESIGHT", "FOV_CLASS_SPEC", "FOV_BOUNDARY",
    "FOV_BOUNDARY_CORNERS", "FOV_REF_VECTOR", "FOV_REF_ANGLE", "FOV_CROSS_ANGLE", "FOV_ANGLE_UNITS"]

"Define the scratch instrument in the two kernel pools: only the variables in `vars` exist."
function fov_define!(ctx, vars::Pair...)
    for v in FOV_VARS
        both_dvpool(ctx, "INS$(FOV_SCRATCH)_$v")
    end
    for (k, v) in vars
        v isa AbstractString ? both_pcpool(ctx, "INS$(FOV_SCRATCH)_$k", v) : both_pdpool(ctx, "INS$(FOV_SCRATCH)_$k", v)
    end
    return nothing
end

"Compare `getfov` for the scratch instrument: the same error kind, or the same values."
function fov_same(ctx, room=100)
    kc = c_kind(() -> g_getfov(FOV_SCRATCH, room))
    kn = neo_kind(() -> getfov(ctx, FOV_SCRATCH, room))
    kc == kn || (println("  getfov kinds: CSPICE ", repr(kc), ", NeoSpice ", repr(kn)); return false, kc)
    kc === nothing || return true, kc
    a = getfov(ctx, FOV_SCRATCH, room)
    b = g_getfov(FOV_SCRATCH, room)
    ok = a[1] == b[1] && a[2] == b[2] && grel(a[3], b[3]) < 1e-15 && length(a[4]) == size(b[4], 2) &&
         all(grel(a[4][i], b[4][:, i]) < 1e-14 for i in eachindex(a[4]))
    return ok, nothing
end

@gf_light group("geom_fov: getfov, getfvn, fovtrg, fovray", Kernels.FULL_SET_FILES...) do
    ctx = geo_setup()
    rng = MersenneTwister(SEED + 107)
    allinst = [GEO_INSTRUMENTS..., "GEO_WIDE", "GEO_FLAT"]

    @testset "getfov and getfvn: the instruments of the kernel" begin
        worst = 0.0
        for (k, name) in enumerate(allinst)
            id = -999000 - k
            a = getfov(ctx, id)
            b = g_getfov(id)
            @test a[1] == b[1] && a[2] == b[2] == "GEO_CAM"
            @test length(a[4]) == size(b[4], 2)
            worst = max(worst, grel(a[3], b[3]), maximum(grel(a[4][i], b[4][:, i]) for i in eachindex(a[4])))
            c = getfvn(ctx, name)
            @test c == a
            d = g_getfvn(name)
            @test c[1] == d[1] && length(c[4]) == size(d[4], 2)
            @test a isa Tuple{String,String,SVector{3,Float64},Vector{SVector{3,Float64}}}
        end
        measured(@sprintf("  %-64s max=%.3e", "getfov: boresight and boundary vectors", worst))
        @test worst < 1e-15
        @test getfov(ctx, -999003)[1] == "RECTANGLE" && length(getfov(ctx, -999003)[4]) == 4
        @test getfov(ctx, -999002)[1] == "ELLIPSE"                 # the kernel has 'ellipse'
        @test same_kind(() -> g_getfvn("NO_INSTRUMENT"), () -> getfvn(ctx, "NO_INSTRUMENT"); expect=:idcodenotfound)
        @test getfvn(ctx, "-999001") == getfov(ctx, -999001)
    end

    @testset "getfov: every form of the specification and every error" begin
        F = "FOV_FRAME" => "GEO_CAM"
        B = "BORESIGHT" => [0.0, 0.0, 1.0]
        A = "FOV_CLASS_SPEC" => "ANGLES"
        R = "FOV_REF_VECTOR" => [1.0, 0.0, 0.0]
        U = "FOV_ANGLE_UNITS" => "DEGREES"
        sq = [1.0, 1.0, 9.0, -1.0, 1.0, 9.0, -1.0, -1.0, 9.0, 1.0, -1.0, 9.0]
        cases = [
            (:framemissing, 100, ["FOV_SHAPE" => "CIRCLE", B]),
            (:shapemissing, 100, [F, B]),
            (:shapenotsupported, 100, [F, "FOV_SHAPE" => "TRIANGLE", B]),
            (:shapenotsupported, 100, [F, "FOV_SHAPE" => "CIRC", B]),
            (:boresightmissing, 100, [F, "FOV_SHAPE" => "CIRCLE"]),
            (:badboresightspec, 100, [F, "FOV_SHAPE" => "CIRCLE", "BORESIGHT" => [0.0, 1.0]]),
            (:badboresightspec, 100, [F, "FOV_SHAPE" => "CIRCLE", "BORESIGHT" => "Z"]),
            (:zeroboresight, 100, [F, "FOV_SHAPE" => "CIRCLE", "BORESIGHT" => zeros(3)]),
            (:boundarymissing, 100, [F, "FOV_SHAPE" => "CIRCLE", B]),
            (:badboundary, 100, [F, "FOV_SHAPE" => "CIRCLE", B, "FOV_BOUNDARY" => [0.1, 0.0, 1.0, 0.0]]),
            (:badboundary, 100, [F, "FOV_SHAPE" => "CIRCLE", B, "FOV_BOUNDARY" => sq[1:6]]),
            (:badboundary, 100, [F, "FOV_SHAPE" => "ELLIPSE", B, "FOV_BOUNDARY" => sq[1:3]]),
            (:badboundary, 100, [F, "FOV_SHAPE" => "RECTANGLE", B, "FOV_BOUNDARY" => sq[1:9]]),
            (:badboundary, 100, [F, "FOV_SHAPE" => "POLYGON", B, "FOV_BOUNDARY" => sq[1:6]]),
            (:boundarytoobig, 3, [F, "FOV_SHAPE" => "RECTANGLE", B, "FOV_BOUNDARY" => sq]),
            (:unsupportedspec, 100, [F, "FOV_SHAPE" => "CIRCLE", B, "FOV_CLASS_SPEC" => "FOO", "FOV_BOUNDARY" => sq[1:3]]),
            (:shapenotsupported, 100, [F, "FOV_SHAPE" => "POLYGON", B, A, R, U, "FOV_REF_ANGLE" => [1.0]]),
            (:refvectormissing, 100, [F, "FOV_SHAPE" => "CIRCLE", B, A, U, "FOV_REF_ANGLE" => [1.0]]),
            (:badrefvectorspec, 100, [F, "FOV_SHAPE" => "CIRCLE", B, A, "FOV_REF_VECTOR" => [1.0, 0.0], U, "FOV_REF_ANGLE" => [1.0]]),
            (:badrefvectorspec, 100, [F, "FOV_SHAPE" => "CIRCLE", B, A, "FOV_REF_VECTOR" => "X", U, "FOV_REF_ANGLE" => [1.0]]),
            (:badrefvectorspec, 100, [F, "FOV_SHAPE" => "CIRCLE", B, A, "FOV_REF_VECTOR" => [0.0, 0.0, 2.0], U, "FOV_REF_ANGLE" => [1.0]]),
            (:refanglemissing, 100, [F, "FOV_SHAPE" => "CIRCLE", B, A, R, U]),
            (:unitsmissing, 100, [F, "FOV_SHAPE" => "CIRCLE", B, A, R, "FOV_REF_ANGLE" => [1.0]]),
            (:crossanglemissing, 100, [F, "FOV_SHAPE" => "ELLIPSE", B, A, R, U, "FOV_REF_ANGLE" => [1.0]]),
            (:crossanglemissing, 100, [F, "FOV_SHAPE" => "RECTANGLE", B, A, R, U, "FOV_REF_ANGLE" => [1.0]]),
            (:boundarytoobig, 0, [F, "FOV_SHAPE" => "CIRCLE", B, A, R, U, "FOV_REF_ANGLE" => [1.0]]),
            (:boundarytoobig, 1, [F, "FOV_SHAPE" => "ELLIPSE", B, A, R, U, "FOV_REF_ANGLE" => [1.0], "FOV_CROSS_ANGLE" => [2.0]]),
            (:boundarytoobig, 3, [F, "FOV_SHAPE" => "RECTANGLE", B, A, R, U, "FOV_REF_ANGLE" => [1.0], "FOV_CROSS_ANGLE" => [2.0]]),
            (:badboundary, 100, [F, "FOV_SHAPE" => "RECTANGLE", B, A, R, U, "FOV_REF_ANGLE" => [90.0], "FOV_CROSS_ANGLE" => [2.0]]),
            (:badboundary, 100, [F, "FOV_SHAPE" => "RECTANGLE", B, A, R, U, "FOV_REF_ANGLE" => [1.0], "FOV_CROSS_ANGLE" => [95.0]]),
            # Valid definitions.
            (nothing, 100, [F, "FOV_SHAPE" => "  circle", B, "FOV_BOUNDARY" => sq[1:3]]),
            (nothing, 100, [F, "FOV_SHAPE" => "ELLIPSE", B, "FOV_BOUNDARY_CORNERS" => sq[1:6]]),
            (nothing, 4, [F, "FOV_SHAPE" => "RECTANGLE", B, "FOV_CLASS_SPEC" => "corners", "FOV_BOUNDARY" => sq]),
            (nothing, 100, [F, "FOV_SHAPE" => "POLYGON", B, "FOV_BOUNDARY" => vcat(sq, [0.0, -3.0, 9.0])]),
            (nothing, 1, [F, "FOV_SHAPE" => "CIRCLE", "BORESIGHT" => [0.3, -0.2, 2.0], A, "FOV_REF_VECTOR" => [1.0, 2.0, 3.0], U, "FOV_REF_ANGLE" => [12.0]]),
            (nothing, 2, [F, "FOV_SHAPE" => "ELLIPSE", "BORESIGHT" => [0.3, -0.2, 2.0], A, "FOV_REF_VECTOR" => [1.0, 2.0, 3.0],
                "FOV_ANGLE_UNITS" => "ARCSECONDS", "FOV_REF_ANGLE" => [7200.0], "FOV_CROSS_ANGLE" => [1800.0]]),
            (nothing, 4, [F, "FOV_SHAPE" => "RECTANGLE", "BORESIGHT" => [0.3, -0.2, 2.0], "FOV_CLASS_SPEC" => "Angles", "FOV_REF_VECTOR" => [1.0, 2.0, 3.0],
                "FOV_ANGLE_UNITS" => "RADIANS", "FOV_REF_ANGLE" => [0.7], "FOV_CROSS_ANGLE" => [0.2]])]
        kinds_seen = Set{Symbol}()
        for (expect, room, vars) in cases
            fov_define!(ctx, vars...)
            ok, kind = fov_same(ctx, room)
            @test ok
            @test kind == expect
            kind === nothing || push!(kinds_seen, kind)
        end
        @test length(kinds_seen) == 15
        # Units that `convrt` does not know.
        fov_define!(ctx, F, "FOV_SHAPE" => "CIRCLE", B, A, R, "FOV_ANGLE_UNITS" => "FURLONGS", "FOV_REF_ANGLE" => [1.0])
        @test c_kind(() -> g_getfov(FOV_SCRATCH)) !== nothing
        @test neo_kind(() -> getfov(ctx, FOV_SCRATCH)) !== nothing
        fov_define!(ctx)
    end

    @testset "fovtrg and fovray, 2000 random cases" begin
        nbad = 0
        ntrg = 0
        ntrgvis = 0
        nray = 0
        nrayvis = 0
        nwide = 0
        ncases = nsample(500, 2000)
        for k in 1:ncases
            inst = rand(rng, allinst)
            et = geo_epoch(rng)
            t = rand(rng, ["SUN", "MOON", "VENUS", "MERCURY"])
            cor = rand(rng, GEO_CORRS)
            r, _ = SPICE.spkpos(t, et, "J2000", cor, "EARTH")
            ang = asin(SPICE.bodvrd(t, "RADII", 3)[1] / norm(r))
            # The offset of the boresight from the target: up to the size of the field of
            # view plus the target (the wide polygon is 80 degrees across).
            off = (0.035 * rand(rng) + 1.3 * ang * rand(rng)) * (inst == "GEO_WIDE" ? 30 : 1)
            zdir = geo_offset_dir(rng, r, off)
            point_camera!(ctx, zdir, 2pi * rand(rng))
            sh = rand(rng, ["ELLIPSOID", "POINT", "ellipsoid", " Point"])
            fr = occursin("oint", sh) && rand(rng, Bool) ? " " : "IAU_" * t
            if inst == "GEO_FLAT"
                nwide += same_kind(() -> g_fovtrg(inst, t, sh, fr, cor, "EARTH", et), () -> fovtrg(ctx, inst, t, sh, fr, cor, "EARTH", et); expect=:fovtoowide)
                continue
            end
            a = fovtrg(ctx, inst, t, sh, fr, cor, "EARTH", et)
            b = g_fovtrg(inst, t, sh, fr, cor, "EARTH", et)
            ntrg += 1
            ntrgvis += b
            a == b || (nbad += 1; nbad <= 5 && println("  fovtrg differs: ", (inst, t, sh, cor, et, off), " NeoSpice ", a, ", CSPICE ", b))
            if k % 20 == 0
                @test fovtrg(ctx, bodn2c(ctx, inst), GEO_IDS[t], sh, "IAU_" * t, Aberration(cor), 399, et) == a
            end
            # A ray near the edge of the field of view, in one of several frames.
            rc = rand(rng, ["NONE", "S", "XS", "none", " xs"])
            rf = rand(rng, ["J2000", "IAU_EARTH", "GEO_CAM", "ECLIPJ2000", "GEO_SITE"])
            dj = geo_offset_dir(rng, zdir, off) .* (0.2 + 4rand(rng))
            d = SPICE.pxform("J2000", rf, et) * dj
            a = fovray(ctx, inst, d, rf, rc, "EARTH", et)
            b = g_fovray(inst, d, rf, rc, "EARTH", et)
            nray += 1
            nrayvis += b
            a == b || (nbad += 1; nbad <= 5 && println("  fovray differs: ", (inst, rf, rc, et, off), " NeoSpice ", a, ", CSPICE ", b))
        end
        measured("  fovtrg: $ntrg cases, $ntrgvis visible; fovray: $nray cases, $nrayvis visible; differences: $nbad; too wide: $nwide")
        @test nbad == 0
        @test 0.2 * ntrg < ntrgvis < 0.8 * ntrg
        @test 0.2 * nray < nrayvis < 0.8 * nray
        @test nwide > ncases ÷ 20
    end

    @testset "error conditions, kinds of CSPICE" begin
        et = 1.0e8
        point_camera!(ctx, [0.0, 0.0, 1.0], 0.0)
        C = "GEO_CIRCLE"
        for (inst, t, sh, fr, cor, obs) in (
                (C, "MOON", "RAY", "IAU_MOON", "LT", "EARTH"),
                (C, "MOON", "ELLIPSOID", "IAU_MOON", "LT", "NOBODY"),
                (C, "MOON", "SPHERE", "IAU_MOON", "LT", "EARTH"),
                (C, "MOON", "DSK/UNPRIORITIZED", "IAU_MOON", "LT", "EARTH"),
                (C, "NOBODY", "ELLIPSOID", "IAU_MOON", "LT", "EARTH"),
                (C, "EARTH", "ELLIPSOID", "IAU_EARTH", "LT", "EARTH"),
                (C, "EARTH", "POINT", " ", "LT", "EARTH"),
                (C, "MOON", "ELLIPSOID", " ", "LT", "EARTH"),
                (C, "MOON", "ELLIPSOID", "NOFRAME", "LT", "EARTH"),
                (C, "MOON", "ELLIPSOID", "IAU_EARTH", "LT", "EARTH"),
                (C, "MOON", "ELLIPSOID", "IAU_MOON", "FOO", "EARTH"),
                (C, "MOON", "POINT", " ", "S", "EARTH"),
                (C, "EARTH BARYCENTER", "ELLIPSOID", "J2000", "LT", "MOON"),
                ("NO_INSTRUMENT", "MOON", "ELLIPSOID", "IAU_MOON", "LT", "EARTH"),
                ("MOON", "SUN", "ELLIPSOID", "IAU_SUN", "LT", "EARTH"),
                ("GEO_FLAT", "MOON", "POINT", " ", "LT", "EARTH"))
            @test same_kind(() -> g_fovtrg(inst, t, sh, fr, cor, obs, et), () -> fovtrg(ctx, inst, t, sh, fr, cor, obs, et))
        end
        z = [0.0, 0.0, 1.0]
        for (inst, d, rf, cor, obs) in (
                (C, z, "J2000", "NONE", "NOBODY"),
                (C, z, " ", "NONE", "EARTH"),
                (C, z, "NOFRAME", "NONE", "EARTH"),
                (C, z, "J2000", "LT", "EARTH"),
                (C, z, "J2000", "CN+S", "EARTH"),
                (C, z, "J2000", "XLT", "EARTH"),
                (C, z, "J2000", "FOO", "EARTH"),
                (C, zeros(3), "J2000", "NONE", "EARTH"),
                ("NO_INSTRUMENT", z, "J2000", "NONE", "EARTH"),
                ("GEO_FLAT", z, "J2000", "NONE", "EARTH"))
            @test same_kind(() -> g_fovray(inst, d, rf, cor, obs, et), () -> fovray(ctx, inst, d, rf, cor, obs, et))
        end
        @test neo_kind(() -> fovray(ctx, C, [1.0, 0.0], "J2000", "NONE", "EARTH", et)) == :badarg
        # An instrument frame that is not known, and degenerate polygons.
        B = "BORESIGHT" => [0.0, 0.0, 1.0]
        fov_define!(ctx, "FOV_FRAME" => "NOFRAME", "FOV_SHAPE" => "CIRCLE", B, "FOV_BOUNDARY" => [0.1, 0.0, 1.0])
        both_pcpool(ctx, "NAIF_BODY_NAME", "GEO_SCRATCH")
        both_pdpool(ctx, "NAIF_BODY_CODE", [FOV_SCRATCH])
        S = "GEO_SCRATCH"
        @test c_kind(() -> g_fovray(S, z, "J2000", "NONE", "EARTH", et)) !== nothing
        @test neo_kind(() -> fovray(ctx, S, z, "J2000", "NONE", "EARTH", et)) == :unknownframe
        F = "FOV_FRAME" => "GEO_CAM"
        # A star with three convex and three concave corners, 165 degrees across: no side
        # of the polygon is a side of its convex hull, so the hull search takes pairs of
        # corners that are not neighbours.
        dir(x, y) = [sind(x) * cosd(y), sind(y), cosd(x) * cosd(y)]
        star = vcat(dir(-85, 0), dir(-2.5, 3), dir(80, 12), dir(75, 0), dir(80, -12), dir(-2.5, -3))
        for bnd in (star, [1.0, 0.0, 9.0, 2.0, 0.0, 18.0, 0.0, 1.0, 9.0],              # two parallel vectors
                [1.0, 0.0, 0.1, -1.0, 0.0, 0.1, 0.0, 1.0, -0.1, 0.0, -1.0, -0.1],          # not in one half space
                [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0],            # in a plane
                [1.0, 0.0, 0.02, 0.0, 1.0, 0.02, -1.0, 0.0, 0.02, 0.0, -1.0, 0.02])        # wide, but valid
            fov_define!(ctx, F, "FOV_SHAPE" => "POLYGON", B, "FOV_BOUNDARY" => bnd)
            kc = c_kind(() -> g_fovray(S, z, "GEO_CAM", "NONE", "EARTH", et))
            kn = neo_kind(() -> fovray(ctx, S, z, "GEO_CAM", "NONE", "EARTH", et))
            @test kc == kn
            if kc === nothing
                for d in (z, [0.5, 0.0, 1.0], [-0.9, 0.0, 0.3], [0.6, 0.15, 0.7], [0.0, 0.03, 1.0], [0.0, 0.08, 1.0])
                    @test fovray(ctx, S, d, "GEO_CAM", "NONE", "EARTH", et) == g_fovray(S, d, "GEO_CAM", "NONE", "EARTH", et)
                end
            end
            kc2 = c_kind(() -> g_fovtrg(S, "MOON", "ELLIPSOID", "IAU_MOON", "LT", "EARTH", et))
            kn2 = neo_kind(() -> fovtrg(ctx, S, "MOON", "ELLIPSOID", "IAU_MOON", "LT", "EARTH", et))
            @test kc2 == kn2
        end
        fov_define!(ctx)
    end
end
