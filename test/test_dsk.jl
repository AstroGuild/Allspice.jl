# Group 12: DSK type 2 (triangular plates).
#
# Feature file: features/dsk_type2.feature
#
# CSPICE is the oracle. No DSK file is needed on the machine: the group makes plate
# models (a tessellated ellipsoid, rough random surfaces, cubes, patches of a surface),
# writes them with the CSPICE writer and with the Allspice writer, and compares:
#
#   (a) the spatial index, the bounds of the third coordinate and the files of the two
#       writers (every integer and every double; the bytes outside of record padding);
#   (b) every access call on the files, with both engines;
#   (c) `dskx02` on single segments for random rays;
#   (d) `dskxsi` and `dskxv` on several loaded files: surface lists, coverage in time,
#       three body-fixed frames, a frame whose centre is not the body, rays that hit, miss,
#       graze plate edges and vertices, and start inside;
#   (e) the plate utilities, the surface name calls, `latsrf` and `srfnrm`;
#   (f) the tolerances, kernel management, and every input check;
#   (g) allocation, inference, time per intercept, and use from many threads.
#
# Kernels (frames, body constants, surface names, one SPK segment) are written by the
# group. CSPICE is called from the main task only.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

isdefined(@__MODULE__, :DSKLIB) || include("dskfixtures.jl")

const DSK_FULL = (-pi, pi, -pi / 2, pi / 2)

"Text kernel with the frames, the Mars constants and the surface names of the group."
function dsk_text_kernel(dir)
    path = joinpath(dir, "dsk_test.tf")
    write(path, """
KPL/FK
\\begindata
BODY499_POLE_RA  = ( 317.269202  -0.10927547  0. )
BODY499_POLE_DEC = (  54.432516  -0.05827105  0. )
BODY499_PM       = ( 176.049863 350.891982443297 0. )
BODY499_RADII    = ( 3.0 2.5 2.0 )
BODY401_RADII    = ( 1.3 1.1 0.9 )

FRAME_MARS_TILT          = 1499001
FRAME_1499001_NAME       = 'MARS_TILT'
FRAME_1499001_CLASS      = 4
FRAME_1499001_CLASS_ID   = 1499001
FRAME_1499001_CENTER     = 499
TKFRAME_1499001_RELATIVE = 'IAU_MARS'
TKFRAME_1499001_SPEC     = 'ANGLES'
TKFRAME_1499001_UNITS    = 'DEGREES'
TKFRAME_1499001_AXES     = ( 3 1 2 )
TKFRAME_1499001_ANGLES   = ( 20 -35 50 )

FRAME_MARS_PCK2          = 1499002
FRAME_1499002_NAME       = 'MARS_PCK2'
FRAME_1499002_CLASS      = 2
FRAME_1499002_CLASS_ID   = 499
FRAME_1499002_CENTER     = 499

FRAME_MARS_OFF           = 1499003
FRAME_1499003_NAME       = 'MARS_OFF'
FRAME_1499003_CLASS      = 4
FRAME_1499003_CLASS_ID   = 1499003
FRAME_1499003_CENTER     = 4
TKFRAME_1499003_RELATIVE = 'IAU_MARS'
TKFRAME_1499003_SPEC     = 'ANGLES'
TKFRAME_1499003_UNITS    = 'DEGREES'
TKFRAME_1499003_AXES     = ( 1 2 3 )
TKFRAME_1499003_ANGLES   = ( 10 0 0 )

FRAME_PHOBOS_FIX         = 1401001
FRAME_1401001_NAME       = 'PHOBOS_FIX'
FRAME_1401001_CLASS      = 4
FRAME_1401001_CLASS_ID   = 1401001
FRAME_1401001_CENTER     = 401
TKFRAME_1401001_RELATIVE = 'J2000'
TKFRAME_1401001_SPEC     = 'ANGLES'
TKFRAME_1401001_UNITS    = 'DEGREES'
TKFRAME_1401001_AXES     = ( 3 2 1 )
TKFRAME_1401001_ANGLES   = ( 30 40 50 )

NAIF_SURFACE_NAME = ( 'Ellipsoid One', 'ROUGH', 'tilted  cube', 'rough', 'old name', 'Old Name', 'p shape' )
NAIF_SURFACE_CODE = ( 1, 2, 3, 22, 5, 6, 1 )
NAIF_SURFACE_BODY = ( 499, 499, 499, 499, 499, 499, 401 )
\\begintext
""")
    return path
end

"SPK file with Mars (499) at a fixed place relative to the Mars barycentre (4)."
function dsk_offset_spk(dir)
    path = joinpath(dir, "mars_offset.bsp")
    isfile(path) && rm(path)
    h = spkopn(path, "offset", 0)
    states = repeat([0.31, -0.22, 0.13, 0.0, 0.0, 0.0], 1, 4)
    spkw08(h, 499, 4, "J2000", -1.5e9, 1.5e9, "mars offset", 1, states, -2.0e9, 1.5e9)
    spkcls(h)
    return path
end

"The segments of the three test files."
function dsk_specs(rng)
    ev, ep = ellipsoid_shape(3.0, 2.5, 2.0, 60, 30)
    rv, rp = rough_shape(rng, 3.2, 0.08, 72, 36)
    r2v, r2p = rough_shape(rng, 3.6, 0.05, 48, 24)
    wv, wp = patch_shape(r2v, r2p, -pi, 0.0, -pi / 2, pi / 2)
    xv, xp = patch_shape(r2v, r2p, 0.0, pi, -pi / 2, pi / 2)
    s4v, s4p = ellipsoid_shape(4.0, 4.0, 3.9, 36, 18)
    s5v, s5p = ellipsoid_shape(4.5, 4.4, 4.3, 36, 18)
    e4v, e4p = ellipsoid_shape(2.8, 2.6, 2.3, 40, 20)     # crosses the ellipsoid of surface 1
    cv, cp_ = cube_shape(1.0)
    c2v, c2p = cube_shape(1.5)
    pv, pp = ellipsoid_shape(1.3, 1.1, 0.9, 30, 16)
    a = [DskSpec(499, 1, 2, "IAU_MARS", 1, Float64[], DSK_FULL, -1e9, 1e9, 5.0, 4, ev, ep),
        DskSpec(499, 2, 2, "IAU_MARS", 1, Float64[], DSK_FULL, -1e9, 1e9, 3.0, 3, rv, rp),
        DskSpec(499, 3, 2, "MARS_TILT", 3, Float64[], (-1.0, 1.0, -1.0, 1.0), -1e9, 1e9, 1.0, 1, cv, cp_)]
    b = [DskSpec(499, 4, 1, "MARS_PCK2", 4, [3.0, 1 / 3], DSK_FULL, -1e9, 1e9, 5.0, 4, e4v, e4p),
        DskSpec(499, 5, 2, "IAU_MARS", 1, Float64[], (-pi, 0.0, -pi / 2, pi / 2), -1e9, 1e9, 4.0, 2, wv, wp),
        DskSpec(499, 5, 2, "IAU_MARS", 1, Float64[], (0.0, pi, -pi / 2, pi / 2), -1e9, 1e9, 4.0, 2, xv, xp),
        DskSpec(499, 6, 1, "IAU_MARS", 1, Float64[], DSK_FULL, 0.0, 100.0, 4.0, 3, s4v, s4p),
        DskSpec(499, 6, 1, "IAU_MARS", 1, Float64[], DSK_FULL, 50.0, 200.0, 4.0, 3, s5v, s5p),
        DskSpec(499, 7, 2, "MARS_OFF", 3, Float64[], (-1.5, 1.5, -1.5, 1.5), -1e9, 1e9, 1.0, 1, c2v, c2p)]
    c = [DskSpec(401, 1, 1, "PHOBOS_FIX", 4, [1.3, 0.2], DSK_FULL, -1e9, 1e9, 4.0, 3, pv, pp)]
    return a, b, c
end

"""
Random rays for a body of size about 3: `kind` 1 from outside toward the body, 2 from
outside in a random direction, 3 from inside, 4 from outside toward a vertex, an edge
point or the centroid of a plate of `shape` (a grazing ray), 5 with the vertex on a
plate vertex.
"""
function dsk_ray(rng, kind::Int, shape)
    verts, plates = shape
    if kind == 1
        v = V3(randn(rng, 3)) * (4 + 6 * rand(rng))
        return v, V3(randn(rng, 3)) * 2.2 - v
    elseif kind == 2
        return V3(randn(rng, 3)) * (4 + 6 * rand(rng)), V3(randn(rng, 3))
    elseif kind == 3
        return V3(randn(rng, 3)) * 0.8, V3(randn(rng, 3))
    elseif kind == 4
        p = plates[rand(rng, 1:length(plates))]
        a, b, c = verts[p[1]], verts[p[2]], verts[p[3]]
        w = rand(rng, 1:4)
        t = rand(rng)
        target = w == 1 ? a : w == 2 ? (a + b) / 2 : w == 3 ? a + t * (c - a) : (a + b + c) / 3
        v = target * (1.5 + 2 * rand(rng)) + V3(randn(rng, 3)) * (w == 4 ? 3.0 : 0.3)
        return v, target - v
    end
    v = verts[rand(rng, 1:length(verts))]
    return v, rand(rng, Bool) ? -v + V3(randn(rng, 3)) * 0.1 : V3(randn(rng, 3))
end

alloc_dskxsi(ctx, srf, et, fixfid, v, d) = @allocated dskxsi(ctx, false, 499, srf, et, fixfid, v, d)
alloc_dskxv!(x, f, ctx, srf, et, fixfid, v, d) = @allocated dskxv!(x, f, ctx, false, 499, srf, et, fixfid, v, d)
function alloc_dskx02(ks, seg, v, d)
    file = ks.dskfiles[seg.file]
    return @allocated Allspice._dsk_x02(file.ints, file.dps, seg, Allspice._DSK_DEFAULT_TOLS, v, d)
end
dsk_warm(f, args...) = (f(args...); f(args...))

"""
Compare `dskxsi` of the two engines for one ray. Returns (found by CSPICE, the two
agree on found and on the segment, the plates agree, the CSPICE plate is that of
another segment, |dxpt|, the segment that wins has the frame of the ray).
"""
function compare_xsi(ctx, ks, target, srf, et, fixref, v, d)
    cf, cx, cp, chh, cdla, cdsk = c_dskxsi(target, srf, et, fixref, v, d)
    r = dskxsi(ctx, false, target, srf, et, fixref, v, d)
    cf == (r !== nothing) || return cf, false, false, false, NaN, true, NaN
    cf || return false, true, true, false, 0.0, true, 0.0
    xpt, seg, plid = r
    same = dla_tuple(cdla) == dla_tuple(seg.dladsc) && dsk_tuple(cdsk) == dsk_tuple(seg.dskdsc) &&
           SPICE.dashfn(chh) == ks.dskfiles[seg.file].path
    stale = false
    fixfid = Allspice.namfrm(ctx, fixref)
    inframe = seg.dskdsc.frmcde == fixfid
    if same && plid != cp
        # CSPICE returns the plate of the last searched segment that has an intercept,
        # which is not always the segment with the closest one. The plate of
        # Allspice must be the one that `dskx02` of CSPICE finds in the segment that
        # both engines return, for the ray in the frame of that segment.
        # CSPICE searches a segment in the frame of the ray from the vertex moved onto
        # the outer bounding sphere, unless a segment in another frame takes part.
        sv, sd = v, d
        body = Int32(seg.dskdsc.center)
        segs = ks.dsk[body]
        multfr = any(Allspice._dsk_eligible(g, body, srf, et) && g.dskdsc.frmcde != fixfid for g in segs)
        if !inframe
            _, sv, sd = Allspice._dsk_to_segment(ks, seg.dskdsc, Int32(fixfid), et, v, d)
        elseif !multfr
            maxrad = Allspice._dsk_sph(ks, segs, body, srf)[2] * 1.01
            norm(v) > maxrad && (sv = Allspice._dsk_surfpt(v, d, maxrad, maxrad, maxrad)[1])
        end
        stale = c_dskx02(chh, cdla, sv, sd)[2] == plid
    end
    dx = cx - xpt
    across = norm(dx - d * (dot(dx, d) / dot(d, d)))
    return true, same, plid == cp, stale, norm(dx), inframe, across
end

"Part of group 12: (a) spatial index, bounds and files of the two writers."
function dsk_part_1(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    nident = 0
    for (name, specs) in ((:a, specs_a), (:b, specs_b), (:c, specs_c))
        # The file comparison below holds the spatial index of every segment (with the
        # vertex-plate list). The direct comparison of `dskmi2`, with and without that
        # list and in the matrix form, runs on the first segment of each file in the
        # default run and on every segment in the exhaustive run.
        for s in (FULL ? specs : specs[1:1]), makvtl in (true, false)
            z = index_sizes(s)
            cd, ci = c_dskmi2(s; makvtl)
            nd, ni = dskmi2(s.verts, s.plates, s.finscl, s.corscl, z.worksz, z.voxpsz,
                z.voxlsz, makvtl, index_spxisz(s))
            @test bits_equal(cd, nd)
            @test ci == ni
            # The matrix form gives the same index.
            md, mi = dskmi2(reduce(hcat, s.verts), reduce(hcat, s.plates), s.finscl, s.corscl,
                z.worksz, z.voxpsz, z.voxlsz, makvtl, index_spxisz(s))
            @test bits_equal(md, nd) && mi == ni
        end
        # The coordinate system is the branch of `dskrb2`; the shapes of file `a` (ellipsoid,
        # rough surface, cube) are the default sample.
        for s in (FULL || name === :a ? specs : specs[1:1]), (corsys, corpar) in ((1, zeros(10)), (3, zeros(10)),
                (4, [3.0, 1 / 3, 0, 0, 0, 0, 0, 0, 0, 0]), (4, [2.0, -0.4, 0, 0, 0, 0, 0, 0, 0, 0]))
            c3 = c_dskrb2(s.verts, s.plates, corsys, corpar)
            n3 = dskrb2(s.verts, s.plates, corsys, corpar)
            @test c3 === n3
        end
        cpath = c_write_dsk(joinpath(dir, "c_$name.bds"), specs; ncomch=300, comments)
        npath = neo_write_dsk(joinpath(dir, "n_$name.bds"), specs; ncomch=300, comments,
            ctx=wctx)
        files[name] = (cpath, npath)
        ci, cd, _ = c_das_data(cpath)
        ni, nd, _ = neo_das_data(npath)
        @test ci == ni
        @test bits_equal(cd, nd)
        @test c_das_data(npath)[1] == ci                    # CSPICE reads the Allspice file
        @test bits_equal(c_das_data(npath)[2], cd)
        @test neo_das_data(cpath)[1] == ci                  # and the reverse
        bc = read(cpath)
        bn = read(npath)
        @test length(bc) == length(bn)
        # The bytes can differ only in the padding after the last integer and the last
        # double: CSPICE writes there what its record buffers held before.
        h = dasopr(npath)
        _, _, ncomr, _, free, lastla, _, _ = dashfs(h)
        dascls(h)
        dprec = cld(lastla[2], 128)
        first_data = (2 + ncomr) * 1024
        pad = falses(length(bn))
        pad[first_data+lastla[2]*8+1:first_data+dprec*1024] .= true
        pad[first_data+dprec*1024+lastla[3]*4+1:end] .= true
        @test all(bc[.!pad] .== bn[.!pad])
        @test all(==(0x00), bn[pad])
        nident += bc == bn
    end
    measured("  DSK writer: $(nident) of 3 files are equal to the CSPICE files byte " *
             "for byte; the others differ only in record padding")
    # A file without the vertex-plate list.
    c0 = c_write_dsk(joinpath(dir, "c_novtl.bds"), specs_a[1:2]; makvtl=false)
    n0 = neo_write_dsk(joinpath(dir, "n_novtl.bds"), specs_a[1:2]; makvtl=false, ctx=wctx)
    @test c_das_data(c0)[1] == neo_das_data(n0)[1]
    @test bits_equal(c_das_data(c0)[2], neo_das_data(n0)[2])
    return nothing
end

"Part of group 12: (b) access calls."
function dsk_part_2(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    # The two files of a pair hold the same integers and doubles (part a), so the default
    # run reads the CSPICE file of each pair and the Allspice file of the first pair only.
    for name in (:a, :b, :c), path in (FULL || name === :a ? files[name] : files[name][1:1])
        @test dskobj(path) == collect(SPICE.dskobj!(SPICE.SpiceIntCell(10), path))
        for body in (499, 401, 5)
            @test dsksrf(path, body) == collect(SPICE.dsksrf!(SPICE.SpiceIntCell(10), path, body))
        end
        h = dasopr(path)
        ch = c_dasopr(path)
        d = dlabfs(h)
        cd = SPICE.dlabfs(ch)
        while d !== nothing
            @test dsk_tuple(dskgd(h, d)) == dsk_tuple(SPICE.dskgd(ch, cd))
            nb = dskb02(h, d)
            cb = c_dskb02(ch, cd)
            @test nb[1:3] == cb[1:3] && nb[8:11] == cb[8:11]
            @test nb[4] == cb[4] && nb[5] === cb[5] && nb[6] == cb[6] && nb[7] == cb[7]
            nv, np = dskz02(h, d)
            @test (nv, np) == SPICE.dskz02(ch, cd)
            @test dskv02(h, d, 1, nv) == c_dskv02(ch, cd, 1, nv)
            @test dskp02(h, d, 1, np) == c_dskp02(ch, cd, 1, np)
            for _ in 1:20
                a = rand(rng, 1:nv); r = rand(rng, 1:50)
                @test dskv02(h, d, a, r) == c_dskv02(ch, cd, a, r)
                a = rand(rng, 1:np); r = rand(rng, 1:50)
                @test dskp02(h, d, a, r) == c_dskp02(ch, cd, a, r)
            end
            # Every plate in the exhaustive run, every tenth plate in the default run.
            @test all(dskn02(h, d, k) == c_dskn02(ch, cd, k) for k in 1:(FULL ? 1 : 10):np)
            for item in 1:14
                full = dski02(h, d, item, 1, 10_000_000)
                @test full == c_dski02(ch, cd, item, 1, 2_000_000)
                n = length(full)
                for _ in 1:(n > 1 ? 5 : 0)
                    a = rand(rng, 1:n); r = rand(rng, 1:100)
                    @test dski02(h, d, item, a, r) == c_dski02(ch, cd, item, a, r)
                end
            end
            for item in 15:19
                full = dskd02(h, d, item, 1, 10_000_000)
                @test bits_equal(full, c_dskd02(ch, cd, item, 1, 2_000_000))
                n = length(full)
                for _ in 1:(n > 1 ? 5 : 0)
                    a = rand(rng, 1:n); r = rand(rng, 1:100)
                    @test bits_equal(dskd02(h, d, item, a, r), c_dskd02(ch, cd, item, a, r))
                end
            end
            d = dlafns(h, d)
            cd = SPICE.dlafns(ch, cd)
        end
        c_dascls(ch)
        dascls(h)
    end
    return nothing
end

"Part of group 12: (c) dskx02 on single segments."
function dsk_part_3(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    nray = 0; nhit = 0; nbad = 0; maxerr = 0.0
    for name in (:a, :b, :c)
        specs = name === :a ? specs_a : name === :b ? specs_b : specs_c
        path = files[name][1]
        h = dasopr(path)
        ch = c_dasopr(path)
        d = dlabfs(h)
        cd = SPICE.dlabfs(ch)
        for s in specs
            scale = maximum(norm, s.verts) / 3
            for k in 1:nsample(100, 2500)
                v, dirv = dsk_ray(rng, rand(rng, 1:5), (s.verts, s.plates))
                kind_scale = k % 2 == 0 ? scale : 1.0
                v *= kind_scale
                cf, cp, cx = c_dskx02(ch, cd, v, dirv)
                r = dskx02(h, d, v, dirv)
                nray += 1
                if cf != (r !== nothing) || (cf && cp != r[1])
                    nbad += 1
                elseif cf
                    nhit += 1
                    maxerr = max(maxerr, norm(cx - r[2]))
                end
            end
            d = dlafns(h, d)
            cd = SPICE.dlafns(ch, cd)
        end
        c_dascls(ch)
        dascls(h)
    end
    measured(@sprintf("  dskx02: %d rays on 10 segments, %d intercepts, %d differences in found or plate, max |dxpt| = %.3e km",
        nray, nhit, nbad, maxerr))
    @test nbad == 0
    @test maxerr == 0.0
    @test nhit > nsample(250, 8_000)
    return nothing
end

"Part of group 12: (d) dskxsi on the loaded files."
function dsk_part_4(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    surface_lists = ((), [1], [2], [3], [4], [5], [6], [7], [1, 2], [2, 5], [3, 7], [5, 6],
        [7, 6, 5, 4, 3, 2, 1])
    rows = Tuple{String,Int,Int,Int,Int,Int,Float64,Float64}[]
    function run_case(label, nfull, fixref, kinds, shapes, epochs, tol)
        n = nsample(nfull ÷ 25, nfull)
        nhit = 0; nfound = 0; nplate = 0; nstale = 0; maxerr = 0.0; maxrot = 0.0; maxacr = 0.0
        for k in 1:n
            srf = surface_lists[rand(rng, 1:length(surface_lists))]
            et = epochs[rand(rng, 1:length(epochs))]
            v, d = dsk_ray(rng, kinds[rand(rng, 1:length(kinds))], shapes[rand(rng, 1:length(shapes))])
            cf, same, plate, stale, err, inframe, across = compare_xsi(ctx, ks, "MARS", srf, et, fixref, v, d)
            nhit += cf
            nfound += !same
            nstale += stale
            nplate += same && cf && !plate && !stale
            same && cf && inframe && (maxerr = max(maxerr, err))
            same && cf && !inframe && (maxrot = max(maxrot, err))
            same && cf && !inframe && (maxacr = max(maxacr, across))
        end
        push!(rows, (label, n, nhit, nfound, nplate, nstale, maxerr, maxrot))
        @test nfound == 0
        @test nplate == 0
        @test maxerr == 0.0          # the segment that wins has the frame of the ray
        # The ray is rotated into the frame of the segment, so the two rays differ by
        # rounding. Across the ray the intercepts then agree to `tol`. Along the ray
        # the difference grows with the reciprocal of the angle between the ray and
        # the plate, which has no lower bound for a grazing ray.
        # Measured on 65,000 rays: 9e-15 km across the ray, 8e-12 km along it.
        @test maxacr <= tol / 100
        @test maxrot <= 100 * tol
        @test nhit > n ÷ 50          # CSPICE finds no intercept from a vertex below the inner radius bound
    end
    epochs = (0.0, 75.0, 150.0, 1.0e5)
    run_case("IAU_MARS, toward the body", 15_000, "IAU_MARS", (1,), (shape_e,), epochs, 1e-11)
    run_case("IAU_MARS, random direction", 8_000, "IAU_MARS", (2,), (shape_e,), epochs, 1e-11)
    run_case("IAU_MARS, from inside", 8_000, "IAU_MARS", (3,), (shape_e,), epochs, 1e-11)
    run_case("IAU_MARS, vertices and edges", 12_000, "IAU_MARS", (4,), (shape_e, shape_r, shape_w), epochs, 1e-11)
    run_case("IAU_MARS, vertex on the surface", 6_000, "IAU_MARS", (5,), (shape_e, shape_r, shape_w), epochs, 1e-11)
    run_case("MARS_TILT, mixed", 8_000, "MARS_TILT", (1, 2, 3), (shape_e,), epochs, 1e-11)
    run_case("MARS_PCK2, mixed", 8_000, "MARS_PCK2", (1, 2, 3), (shape_e,), epochs, 1e-11)
    measured("  dskxsi against CSPICE (13 surface lists, 4 epochs, 9 segments of body 499):")
    measured("    case                              rays    hits  found/segment diff  plate diff  CSPICE stale plate  max |dxpt| km (same frame)  (other frame)")
    for r in rows
        measured(@sprintf("    %-32s %6d  %6d  %18d  %10d  %18d  %.3e  %.3e", r...))
    end
    @test sum(r -> r[6], rows) > 0   # the defect of CSPICE is seen; see the report

    # The other body, and the longitude seam of the two halves of surface 5.
    nbad = 0
    for k in 1:nsample(300, 3000)
        v = V3(randn(rng, 3)) * 5
        d = V3(randn(rng, 3)) - v
        cf, same, plate, stale, err, _ = compare_xsi(ctx, ks, "PHOBOS", (), 10.0, "PHOBOS_FIX", v, d)
        nbad += !(same && (plate || stale) && err <= 1e-11)
    end
    @test nbad == 0
    nbad = 0; nhit = 0
    nseam = nsample(400, 4000)
    for k in 1:nseam
        lat = (rand(rng) - 0.5) * 3
        lon = rand(rng, (0.0, pi, -pi)) + (rand(rng, Bool) ? 0.0 : 1e-9 * randn(rng))
        target = 3.6 * V3(cos(lon) * cos(lat), sin(lon) * cos(lat), sin(lat))
        v = target * 2 + V3(0.0, 0.0, randn(rng))
        cf, same, plate, stale, err, _ = compare_xsi(ctx, ks, "MARS", [5], 0.0, "IAU_MARS", v, target - v)
        nhit += cf
        nbad += !(same && (plate || stale) && err == 0.0)
    end
    @test nbad == 0 && nhit > 3 * nseam ÷ 4
    return nothing
end

"Part of group 12: (d) dskxv, argument forms and coverage in time."
function dsk_part_5(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    n = nsample(1000, 5000)
    vtx = V3[]; dirs = V3[]
    for k in 1:n
        v, d = dsk_ray(rng, rand(rng, 1:4), shape_r)
        push!(vtx, v); push!(dirs, d)
    end
    for (srf, et, fixref) in (((), 0.0, "IAU_MARS"), ([2, 5], 75.0, "IAU_MARS"),
            ([6], 75.0, "IAU_MARS"), ([3, 7], 0.0, "MARS_TILT"))
        cx, cf = c_dskxv("MARS", srf, et, fixref, vtx, dirs)
        nx, nf = dskxv(ctx, false, "MARS", srf, et, fixref, vtx, dirs)
        @test nf == cf
        @test maximum(norm(cx[i] - nx[i]) for i in 1:n if cf[i]) <= (isempty(srf) || fixref != "IAU_MARS" ? 1e-11 : 0.0)
        @test all(iszero(nx[i]) for i in 1:n if !nf[i])
        # The matrix form, vectors of plain vectors, and integer IDs give the same.
        mx, mf = dskxv(ctx, false, 499, srf, et, Allspice.namfrm(ctx, fixref),
            reduce(hcat, vtx), reduce(hcat, dirs))
        @test mx == nx && mf == nf
        px, pf = dskxv(ctx, false, "mars", collect(srf), et, lowercase(fixref),
            Vector.(vtx[1:50]), Vector.(dirs[1:50]))
        @test px == nx[1:50] && pf == nf[1:50]
    end
    # Coverage in time of surface 6: the first segment only, both, the second only, none.
    v = V3(10.0, 1.0, 2.0)
    r(et) = dskxsi(ctx, false, "MARS", [6], et, "IAU_MARS", v, -v)
    @test r(0.0)[2].dskdsc.stop == 100.0 && r(25.0)[2].dskdsc.stop == 100.0
    @test r(75.0)[2].dskdsc.stop == 200.0 && r(100.0)[2].dskdsc.stop == 200.0
    @test r(150.0)[2].dskdsc.stop == 200.0 && r(200.0) !== nothing
    @test r(-1e-9) === nothing && r(200.0 + 1e-9) === nothing && r(1e5) === nothing
    for et in (0.0, 25.0, 50.0, 75.0, 100.0, 150.0, 200.0, 201.0, -1.0)
        cf, cx = c_dskxsi("MARS", [6], et, "IAU_MARS", v, -v)
        @test cf == (r(et) !== nothing) && (!cf || cx == r(et)[1])
    end
    # The default context and a KernelSet give what the explicit context gives.
    @test dskxsi(ks, false, 499, (), 0.0, "IAU_MARS", v, -v) == dskxsi(ctx, false, 499, (), 0.0, "IAU_MARS", v, -v)
    kclear()
    foreach(furnsh, loadlist)
    @test dskxsi(false, 499, (), 0.0, "IAU_MARS", v, -v) == dskxsi(ctx, false, 499, (), 0.0, "IAU_MARS", v, -v)
    @test dskxv(false, 499, (), 0.0, "IAU_MARS", [v], [-v]) == dskxv(ctx, false, 499, (), 0.0, "IAU_MARS", [v], [-v])
    xpt, seg, plid = dskxsi(false, 499, [2], 0.0, "IAU_MARS", v, -v)
    @test dskgd(seg).surfce == 2 && dskz02(seg) == (length(shape_r[1]), length(shape_r[2]))
    @test dskb02(seg)[1:2] == dskz02(seg)
    @test dskp02(seg, plid, 1) == dskp02(ctx, seg, plid, 1) == [shape_r[2][plid]]
    @test dskv02(seg, 3, 2) == dskv02(ks, seg, 3, 2) == shape_r[1][3:4]
    @test dskn02(seg, plid) == dskn02(ctx, seg, plid) == dskn02(ks, seg, plid)
    @test dski02(seg, 2, 1, 1) == dski02(ctx, seg, 2, 1, 1) == [length(shape_r[2])]
    @test dskd02(seg, 19, 1, 3) == dskd02(ks, seg, 19, 1, 3) == collect(shape_r[1][1])
    @test dskx02(seg, v, -v) == dskx02(ctx, seg, v, -v) == dskx02(ks, seg, v, -v)
    # `dskxsi` moves a far vertex onto a sphere around the body first, so the point can
    # differ from that of `dskx02` in the last bits.
    @test dskx02(seg, v, -v)[1] == plid && norm(dskx02(seg, v, -v)[2] - xpt) <= 1e-13
    @test pltnp(xpt, dskv02(seg, shape_r[2][plid][1], 1)[1], dskv02(seg, shape_r[2][plid][2], 1)[1],
        dskv02(seg, shape_r[2][plid][3], 1)[1])[2] <= 1e-12
    kclear()
    return nothing
end

"Part of group 12: (e) plate utilities."
function dsk_part_6(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    maxd = 0.0
    for _ in 1:nsample(200, 2000)
        p = V3(randn(rng, 3)) * 3
        a, b, c = V3(randn(rng, 3)), V3(randn(rng, 3)), V3(randn(rng, 3))
        rand(rng) < 0.05 && (c = a + rand(rng) * (b - a))      # a degenerate plate
        rand(rng) < 0.02 && (b = a; c = a)
        np_, nd = pltnp(p, a, b, c)
        cp_, cd = SPICE.pltnp(collect(p), collect(a), collect(b), collect(c))
        maxd = max(maxd, norm(np_ - cp_), abs(nd - cd))
        @test pltnrm(a, b, c) == SPICE.pltnrm(collect(a), collect(b), collect(c))
        delta = randn(rng)
        @test collect(pltexp([a, b, c], delta)) == SPICE.pltexp([collect(a), collect(b), collect(c)], delta)
        @test pltexp(hcat(a, b, c), delta) == pltexp((a, b, c) |> collect, delta)
    end
    # Fixed cases: three vertices exactly on one line, with the longest edge first, second
    # and third. A random plate reaches these branches about once in a thousand.
    for (a, b, c) in ((V3(0.0, 0.0, 0.0), V3(2.0, 0.0, 0.0), V3(1.0, 0.0, 0.0)),
            (V3(1.0, 0.0, 0.0), V3(0.0, 0.0, 0.0), V3(2.0, 0.0, 0.0)),
            (V3(2.0, 0.0, 0.0), V3(1.0, 0.0, 0.0), V3(0.0, 0.0, 0.0))),
        p in (V3(0.5, 1.0, 0.0), V3(-1.0, 2.0, 3.0), V3(5.0, 0.0, 1.0))

        np_, nd = pltnp(p, a, b, c)
        cp_, cd = SPICE.pltnp(collect(p), collect(a), collect(b), collect(c))
        maxd = max(maxd, norm(np_ - cp_), abs(nd - cd))
    end
    @test maxd == 0.0
    for s in (specs_a..., specs_b..., specs_c...)
        @test pltar(s.verts, s.plates) === c_pltar(s.verts, s.plates)
        @test pltar(reduce(hcat, s.verts), reduce(hcat, s.plates)) === pltar(s.verts, s.plates)
        length(s.plates) >= 4 && @test pltvol(s.verts, s.plates) === c_pltvol(s.verts, s.plates)
    end
    ev, ep = shape_e
    @test isapprox(pltvol(ev, ep), 4 / 3 * pi * 3.0 * 2.5 * 2.0; rtol=2e-2)
    @test pltvol(cube_shape(1.0)...) ≈ 8.0 && pltar(cube_shape(1.0)...) ≈ 24.0
    @test pltar(ev, P3[]) == 0.0
    @test neo_kind(() -> pltar(ev[1:2], [P3(1, 2, 1)])) === :toofewvertices
    @test neo_kind(() -> pltar(ev, [P3(1, 2, 0)])) === :indexoutofrange
    @test neo_kind(() -> pltvol(ev[1:3], ep[1:4])) === :toofewvertices
    @test neo_kind(() -> pltvol(ev, ep[1:3])) === :toofewplates
    @test neo_kind(() -> pltvol(ev, [P3(1, 2, 3), P3(1, 2, 3), P3(1, 2, 3), P3(1, 2, length(ev) + 1)])) === :indexoutofrange
    @test neo_kind(() -> pltnrm([1.0, 2.0], ev[1], ev[2])) === :badarg
    @test neo_kind(() -> pltexp([ev[1], ev[2]], 0.1)) === :badarg
    return nothing
end

"Part of group 12: (e) surface names and codes."
function dsk_part_7(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    for (code, body) in ((1, 499), (2, 499), (3, 499), (22, 499), (5, 499), (6, 499), (1, 401),
            (7, 499), (1, 5), (-3, 499))
        @test srfc2s(ctx, code, body) == SPICE.srfc2s(code, body)
    end
    for (code, body) in ((1, "MARS"), (22, "499"), (1, "phobos"), (6, "Mars"), (5, "MARS"),
            (1, "NO SUCH BODY"), (9, "MARS"))
        @test srfcss(ctx, code, body) == SPICE.srfcss(code, body)
    end
    for name in ("Ellipsoid One", "ellipsoid   one", " ROUGH ", "rough", "tilted cube",
            "old name", "OLD NAME", "p shape", "12", " -7 ", "nothing", " ", "1.5", "-1.5", "2e1", "x1")
        for body in (499, 401, 5)
            @test srfscc(ctx, name, body) == SPICE.srfscc(name, body)
        end
        for body in ("MARS", "401", "NO SUCH BODY")
            @test srfs2c(ctx, name, body) == SPICE.srfs2c(name, body)
        end
    end
    # An empty string is an error of the C interface of CSPICE.
    c_empty = c_kind(() -> ccall((:srfscc_c, DSKLIB), Cvoid, (Cstring, Cint, Ref{Cint}, Ref{Cint}),
        "", 499, Ref{Cint}(0), Ref{Cint}(0)))
    @test neo_kind(() -> srfscc(ctx, "", 499)) === :emptystring === c_empty
    @test neo_kind(() -> srfs2c(ctx, "", "MARS")) === :emptystring
    @test neo_kind(() -> srfs2c(ctx, "ROUGH", "")) === :emptystring
    @test neo_kind(() -> srfcss(ctx, 1, "")) === :emptystring
    # A context without the variables has no names; a bad map is an error.
    empty = SpiceContext()
    @test srfc2s(empty, 1, 499) == ("1", false) && srfscc(empty, "ROUGH", 499) === nothing
    @test srfscc(empty, "17", 499) == 17
    function bad_map(text)
        p = joinpath(dir, "badmap_$(hash(text)).tk")
        write(p, "KPL/FK\n\\begindata\n" * text * "\n\\begintext\n")
        c = SpiceContext()
        furnsh(c, p)
        SPICE.kclear()
        SPICE.furnsh(p)
        ck = c_kind(() -> ccall((:srfscc_c, DSKLIB), Cvoid, (Cstring, Cint, Ref{Cint}, Ref{Cint}),
            "X", 499, Ref{Cint}(0), Ref{Cint}(0)))
        return neo_kind(() -> srfscc(c, "X", 499)), neo_kind(() -> srfc2s(c, 1, 499)), ck
    end
    for (text, kind) in (("NAIF_SURFACE_NAME = ( 'A' )", :badsurfacemap),
            ("NAIF_SURFACE_NAME = ( 'A' )\nNAIF_SURFACE_CODE = ( 1 )\nNAIF_SURFACE_BODY = ( 'B' )", :badvariabletype),
            ("NAIF_SURFACE_NAME = ( 'A', 'B' )\nNAIF_SURFACE_CODE = ( 1 )\nNAIF_SURFACE_BODY = ( 499, 499 )", :arraysizemismatch),
            ("NAIF_SURFACE_NAME = ( 'A', ' ' )\nNAIF_SURFACE_CODE = ( 1, 2 )\nNAIF_SURFACE_BODY = ( 499, 499 )", :blanknameassigned))
        @test bad_map(text) == (kind, kind, kind)
    end
    oracle_load(loadlist)
    return nothing
end

"Part of group 12: (e) latsrf and srfnrm."
function dsk_part_8(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    lonlat = [[2pi * rand(rng) - pi, pi * rand(rng) - pi / 2] for _ in 1:nsample(100, 400)]
    push!(lonlat, [0.0, pi / 2], [0.0, -pi / 2], [pi, 0.0], [-pi, 0.0])
    for (method, target, fixref, tol) in (("ELLIPSOID", "MARS", "IAU_MARS", 1e-13),
            ("ellipsoid", "499", "MARS_TILT", 1e-13),
            ("DSK/UNPRIORITIZED", "MARS", "IAU_MARS", 1e-13),
            ("DSK/UNPRIORITIZED/SURFACES = 1", "MARS", "IAU_MARS", 1e-13),
            ("dsk/unprioritized/surfaces = \"Ellipsoid One\", rough", "MARS", "IAU_MARS", 1e-13),
            ("UNPRIORITIZED/DSK/SURFACES=2,5", "MARS", "IAU_MARS", 1e-13),
            ("DSK / SURFACES = \"tilted cube\" / UNPRIORITIZED", "MARS", "MARS_TILT", 1e-11),
            ("DSK/UNPRIORITIZED/SURFACES = 7, 3", "MARS", "MARS_PCK2", 1e-11),
            ("DSK/UNPRIORITIZED/SURFACES = \"p shape\"", "PHOBOS", "PHOBOS_FIX", 1e-11))
        for et in (0.0, 75.0)
            cpts = c_latsrf(method, target, et, fixref, lonlat)
            npts = latsrf(ctx, method, target, et, fixref, lonlat)
            @test maximum(norm(cpts[i] - npts[i]) for i in eachindex(cpts)) <= tol
            @test latsrf(ctx, method, target, et, fixref, reduce(hcat, lonlat)) == npts
            cn = c_srfnrm(method, target, et, fixref, cpts)
            nn = srfnrm(ctx, method, target, et, fixref, cpts)
            @test maximum(norm(cn[i] - nn[i]) for i in eachindex(cn)) <= max(tol, 0.0) * 1e3
            @test all(abs(norm(v) - 1) <= 4eps() for v in nn)
            @test srfnrm(ctx, method, target, et, fixref, reduce(hcat, cpts)) == nn
        end
    end
    # Error kinds of the method string, the target and the frame, as CSPICE.
    c_latsrf_kind(method, target, fixref) = c_kind(() -> ccall((:latsrf_c, DSKLIB), Cvoid,
        (Cstring, Cstring, Float64, Cstring, Cint, Ptr{Float64}, Ptr{Float64}), method, target,
        0.0, fixref, 1, [0.1, 0.2], zeros(3)))
    c_srfnrm_kind(method, target, fixref, p) = c_kind(() -> ccall((:srfnrm_c, DSKLIB), Cvoid,
        (Cstring, Cstring, Float64, Cstring, Cint, Ptr{Float64}, Ptr{Float64}), method, target,
        0.0, fixref, 1, collect(p), zeros(3)))
    for (method, target, fixref) in (("DSK", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/SURFACES", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/SURFACES =", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/SURFACES = 1,", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/SURFACES = 1 2", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/SURFACES = \"\"", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/SURFACES = \"open", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/SURFACES = nosuchname", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/UNPRIORITIZED", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/SURFACES = 1/SURFACES = 2", "MARS", "IAU_MARS"),
            ("DSK//UNPRIORITIZED", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/BOGUS", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/NADIR", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/TANGENT", "MARS", "IAU_MARS"),
            ("DSK/ELLIPSOID", "MARS", "IAU_MARS"), ("", "MARS", "IAU_MARS"),
            ("SPHERE", "MARS", "IAU_MARS"), ("ELLIPSOID/NADIR", "MARS", "IAU_MARS"),
            ("ELLIPSOID/NEAR POINT", "MARS", "IAU_MARS"), ("NEAR POINT: ELLIPSOID", "MARS", "IAU_MARS"),
            ("ELLIPSOID/TANGENT", "MARS", "IAU_MARS"), ("ELLIPSOID/BOGUS", "MARS", "IAU_MARS"),
            ("ELLIPSOID ELLIPSOID", "MARS", "IAU_MARS"), ("/ELLIPSOID", "MARS", "IAU_MARS"),
            ("ELLIPSOID", "NO SUCH BODY", "IAU_MARS"), ("ELLIPSOID", "MARS", "NO_SUCH_FRAME"),
            ("ELLIPSOID", "MARS", "IAU_EARTH"), ("ELLIPSOID", "MARS", "J2000"),
            ("DSK/UNPRIORITIZED", "EARTH", "IAU_EARTH"),
            ("DSK/UNPRIORITIZED/SURFACES = 99", "MARS", "IAU_MARS"),
            ("DSK/UNPRIORITIZED/SURFACES = 3", "MARS", "IAU_MARS"))
        nk = neo_kind(() -> latsrf(ctx, method, target, 0.0, fixref, [[0.1, 0.2]]))
        ck = c_latsrf_kind(method, target, fixref)
        @test (method, nk) == (method, ck)
    end
    # (regression: a pair with three components threw `ArgumentError`, not a `SpiceError`)
    @test neo_kind(() -> latsrf(ctx, "ELLIPSOID", "MARS", 0.0, "IAU_MARS", [[0.1, 0.2, 0.3]])) === :badarg
    on = c_latsrf("ELLIPSOID", "MARS", 0.0, "IAU_MARS", [[0.3, 0.4]])[1]
    for (method, p) in (("ELLIPSOID", on * 1.01), ("ELLIPSOID", on * (1 + 1e-9)),
            ("DSK/UNPRIORITIZED/SURFACES = 1", on * 1.5),
            ("DSK/UNPRIORITIZED/SURFACES = 1", on * 0.99),
            ("DSK/UNPRIORITIZED/SURFACES = 3", V3(5.0, 0.0, 0.0)))
        nk = neo_kind(() -> srfnrm(ctx, method, "MARS", 0.0, "IAU_MARS", [p]))
        @test (method, nk) == (method, c_srfnrm_kind(method, "MARS", "IAU_MARS", p))
    end
    @test neo_kind(() -> latsrf(ctx, "ELLIPSOID", "EARTH", 0.0, "IAU_EARTH", [[0.1, 0.2]])) !== nothing
    return nothing
end

"Part of group 12: (f) tolerances."
function dsk_part_9(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    for k in 1:6
        @test dskgtl(ctx, k) === c_dskgtl(k)
    end
    @test neo_kind(() -> dskgtl(ctx, 0)) === :indexoutofrange === c_kind(() -> ccall((:dskgtl_c, DSKLIB), Cvoid, (Cint, Ref{Float64}), 0, Ref(0.0)))
    @test neo_kind(() -> dskgtl(ctx, 7)) === :indexoutofrange
    @test neo_kind(() -> dskstl(ctx, 0, 1.0)) === :indexoutofrange
    @test neo_kind(() -> dskstl(ctx, 5, 1.0)) === :immutablevalue === c_kind(() -> ccall((:dskstl_c, DSKLIB), Cvoid, (Cint, Float64), 5, 1.0))
    @test neo_kind(() -> dskstl(ctx, 6, 1.0)) === :immutablevalue
    # A large plate expansion and a large segment margin change the intercepts in the
    # same way in both engines; another context keeps the default values.
    other = neo_context(loadlist)
    vtx = V3[]; dirs = V3[]
    for k in 1:nsample(1000, 3000)
        v, d = dsk_ray(rng, 4, shape_w)
        push!(vtx, v); push!(dirs, d)
    end
    base = dskxv(ctx, false, "MARS", [5], 0.0, "IAU_MARS", vtx, dirs)
    try
        for (k, val) in ((1, 1e-2), (2, 1e-3), (4, 1e-3))
            c_dskstl(k, val)
            dskstl(ctx, k, val)
            @test dskgtl(ctx, k) === val === c_dskgtl(k)
        end
        cx, cf = c_dskxv("MARS", [5], 0.0, "IAU_MARS", vtx, dirs)
        nx, nf = dskxv(ctx, false, "MARS", [5], 0.0, "IAU_MARS", vtx, dirs)
        @test nf == cf && all(cx[i] == nx[i] for i in eachindex(cx) if cf[i])
        @test (nx, nf) != base
        @test dskxv(other, false, "MARS", [5], 0.0, "IAU_MARS", vtx, dirs) == base
        @test dskgtl(other, 1) == 1e-10 && dskgtl(other, 2) == 1e-8
        on = c_latsrf("ELLIPSOID", "MARS", 0.0, "IAU_MARS", [[0.3, 0.4]])[1] * (1 + 1e-5)
        @test srfnrm(ctx, "ELLIPSOID", "MARS", 0.0, "IAU_MARS", [on]) == c_srfnrm("ELLIPSOID", "MARS", 0.0, "IAU_MARS", [on])
        @test neo_kind(() -> srfnrm(other, "ELLIPSOID", "MARS", 0.0, "IAU_MARS", [on])) === :pointnotonsurface
    finally
        for (k, val) in ((1, 1e-10), (2, 1e-8), (4, 1e-7))
            c_dskstl(k, val)
            dskstl(ctx, k, val)
        end
    end
    @test dskxv(ctx, false, "MARS", [5], 0.0, "IAU_MARS", vtx, dirs) == base
    return nothing
end

"Part of group 12: (f) kernel management."
function dsk_part_10(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    @test ktotal(ctx, "DSK") == 3 == SPICE.ktotal("DSK")
    @test ktotal(ctx, "ALL") == 5 == SPICE.ktotal("ALL")
    for i in 1:3
        n = kdata(ctx, i, "DSK")
        c = SPICE.kdata(i, "DSK")
        @test n[1:3] == c[1:3] && n[2] == "DSK"
    end
    @test kdata(ctx, 4, "DSK") === nothing
    @test kinfo(ctx, files[:b][1])[1] == "DSK" == SPICE.kinfo(files[:b][1])[1]
    @test getfat(files[:a][2]) == ("DAS", "DSK")
    c2 = neo_context(loadlist)
    v = V3(10.0, 1.0, 2.0)
    xpt, seg, plid = dskxsi(c2, false, 499, [5], 0.0, "IAU_MARS", v, -v)
    unload(c2, files[:b][1])
    @test ktotal(c2, "DSK") == 2
    @test neo_kind(() -> dskxsi(c2, false, 499, [5], 0.0, "IAU_MARS", v, -v)) === :dskdatanotfound
    @test neo_kind(() -> dskv02(c2, seg, 1, 1)) === :nosuchhandle
    @test dskxsi(c2, false, 499, [1], 0.0, "IAU_MARS", v, -v) == dskxsi(ctx, false, 499, [1], 0.0, "IAU_MARS", v, -v)
    # A later file has priority only through the order of the search: the same file
    # loaded twice gives the same answers.
    furnsh(c2, files[:b][2])
    @test dskxsi(c2, false, 499, [5], 0.0, "IAU_MARS", v, -v)[1] == xpt
    unload(c2, files[:a][1])
    unload(c2, files[:b][2])
    @test neo_kind(() -> dskxsi(c2, false, 499, (), 0.0, "IAU_MARS", v, -v)) === :nodsksegments
    kclear(c2)
    @test ktotal(c2, "DSK") == 0
    return nothing
end

"Part of group 12: (f) input checks of the intercept and the access calls."
function dsk_part_11(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    v = V3(10.0, 1.0, 2.0)
    c_xsi_kind(target, srf, fixref, vtx, d; pri=false) = c_kind(() -> c_dskxsi(target, srf, 0.0, fixref, vtx, d; pri))
    for (target, srf, fixref, vtx, d, pri) in (("MARS", (), "IAU_MARS", v, -v, true),
            ("NO SUCH BODY", (), "IAU_MARS", v, -v, false),
            ("MARS", (), "NO_SUCH_FRAME", v, -v, false),
            ("MARS", (), "IAU_EARTH", v, -v, false), ("MARS", (), "J2000", v, -v, false),
            ("MARS", [99], "IAU_MARS", v, -v, false), ("MARS", [1, 99], "IAU_MARS", v, -v, false),
            ("EARTH", (), "IAU_EARTH", v, -v, false),
            ("MARS", (), "IAU_MARS", v, V3(0.0, 0.0, 0.0), false),
            ("MARS", [1], "IAU_MARS", V3(0.0, 0.0, 0.0), V3(0.0, 0.0, 0.0), false))
        nk = neo_kind(() -> dskxsi(ctx, pri, target, srf, 0.0, fixref, vtx, d))
        @test (target, fixref, nk) == (target, fixref, c_xsi_kind(target, srf, fixref, vtx, d; pri))
    end
    @test neo_kind(() -> dskxv(ctx, true, "MARS", (), 0.0, "IAU_MARS", [v], [-v])) === :badpriorityspec
    @test neo_kind(() -> dskxv(ctx, false, "MARS", (), 0.0, "IAU_MARS", V3[], V3[])) === :invalidcount
    @test neo_kind(() -> dskxv(ctx, false, "MARS", (), 0.0, "IAU_MARS", [v, v], [-v])) === :badarg
    @test neo_kind(() -> dskxsi(ctx, false, "MARS", (), 0.0, "IAU_MARS", [1.0, 2.0], -v)) === :badarg

    path = files[:a][1]
    h = dasopr(path)
    ch = c_dasopr(path)
    d = dlabfs(h)
    cd = SPICE.dlabfs(ch)
    nv, np = dskz02(h, d)
    c_i02_kind(item, start, room) = c_kind(() -> ccall((:dski02_c, DSKLIB), Cvoid, (Cint,
            Ref{SPICE.DLADescr}, Cint, Cint, Cint, Ref{Cint}, Ptr{Cint}), ch, cd, item, start - 1,
        room, Ref{Cint}(0), zeros(Cint, 10)))
    c_d02_kind(item, start, room) = c_kind(() -> ccall((:dskd02_c, DSKLIB), Cvoid, (Cint,
            Ref{SPICE.DLADescr}, Cint, Cint, Cint, Ref{Cint}, Ptr{Float64}), ch, cd, item, start - 1,
        room, Ref{Cint}(0), zeros(10)))
    for (item, start, room) in ((1, 1, 0), (9, 1, -1), (0, 1, 1), (15, 1, 1), (9, 0, 1),
            (9, 3 * np + 1, 1), (4, 4, 1), (14, 100000, 1), (5, 0, 1))
        @test (item, start, neo_kind(() -> dski02(h, d, item, start, room))) == (item, start, c_i02_kind(item, start, room))
    end
    # Item 1 (the number of vertices) has one element. For a start index above 1 the
    # answer of CSPICE is not stable: the same call gives SPICE(INDEXOUTOFRANGE),
    # SPICE(DASNOSUCHADDRESS) or no error, and which one depends on the calls before it
    # (three equal calls in a row gave two of them). Allspice always gives the first one,
    # and CSPICE is not asked: it is no oracle here.
    @test neo_kind(() -> dski02(h, d, 1, 2, 1)) === :indexoutofrange
    for (item, start, room) in ((19, 1, 0), (14, 1, 1), (20, 1, 1), (19, 0, 1), (19, 3 * nv + 1, 1),
            (15, 25, 1), (16, 7, 1), (17, 4, 1), (18, 2, 1))
        @test (item, start, neo_kind(() -> dskd02(h, d, item, start, room))) == (item, start, c_d02_kind(item, start, room))
    end
    @test neo_kind(() -> dskv02(h, d, 0, 1)) === :indexoutofrange
    @test neo_kind(() -> dskv02(h, d, nv + 1, 1)) === :indexoutofrange
    @test neo_kind(() -> dskp02(h, d, np + 1, 1)) === :indexoutofrange
    @test neo_kind(() -> dskp02(h, d, 1, 0)) === :valueoutofrange
    @test neo_kind(() -> dskn02(h, d, 0)) === :indexoutofrange
    @test neo_kind(() -> dskn02(h, d, np + 1)) === :indexoutofrange
    @test neo_kind(() -> dskx02(h, d, v, V3(0.0, 0.0, 0.0))) === :rayiszerovector
    @test dskv02(h, d, nv, 5) == [shape_e[1][nv]] && dskp02(h, d, np, 5) == [shape_e[2][np]]
    c_dascls(ch)
    dascls(h)
    @test neo_kind(() -> dskz02(h, d)) === :dasnosuchhandle
    @test neo_kind(() -> dskgd(h, d)) === :dasnosuchhandle
    # Files that are not DSK files.
    gen = joinpath(dir, "notdsk.das")
    w = dasonw(gen, "TEST", "x", 0)
    dasadi(w, Cint[1000000, -1, -1])
    dascls(w)
    c_obj_kind(p) = c_kind(() -> (cell = SPICE.SpiceIntCell(10); ccall((:dskobj_c, DSKLIB), Cvoid,
        (Cstring, Ref{SPICE.Cell{Cint}}), p, cell.cell)))
    @test neo_kind(() -> dskobj(gen)) === :invalidfiletype === c_obj_kind(gen)
    @test neo_kind(() -> dskobj(tk)) === :invalidarchtype === c_obj_kind(tk)
    @test neo_kind(() -> dskobj(spk)) === :invalidarchtype === c_obj_kind(spk)
    @test neo_kind(() -> dsksrf(gen, 499)) === :invalidfiletype
    @test neo_kind(() -> dskobj(joinpath(dir, "missing.bds"))) === :filenotfound
    c = SpiceContext()
    @test neo_kind(() -> furnsh(c, gen)) !== nothing && ktotal(c) == 0
    return nothing
end

"Part of group 12: (f) input checks of the writer."
function dsk_part_12(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    s = specs_a[1]
    z = index_sizes(s)
    spaixd, spaixi = c_dskmi2(s)
    mn3, mx3 = c_dskrb2(s.verts, s.plates, 1, zeros(10))
    nmi2(; finscl=s.finscl, corscl=s.corscl, worksz=z.worksz, voxpsz=z.voxpsz, voxlsz=z.voxlsz,
        spxisz=index_spxisz(s), verts=s.verts, plates=s.plates) = neo_kind(() -> dskmi2(verts,
        plates, finscl, corscl, worksz, voxpsz, voxlsz, true, spxisz))
    cmi2(; finscl=s.finscl, corscl=s.corscl, worksz=z.worksz, voxpsz=z.voxpsz, voxlsz=z.voxlsz,
        spxisz=index_spxisz(s)) = c_kind(() -> c_dskmi2(s; finscl, corscl,
        sizes=(; worksz, voxpsz, voxlsz), spxisz))
    for kw in ((finscl=0.0,), (finscl=-1.0,), (corscl=0,), (worksz=10,), (voxpsz=0,), (voxlsz=10,),
            (spxisz=1000,), (corscl=7,), (corscl=100,), (finscl=1e-3,), (voxpsz=100,),
            (worksz=length(s.plates) + 1, voxpsz=10), (voxlsz=length(s.plates) + 1,),
            (finscl=0.3, corscl=1))
        @test (kw, nmi2(; kw...)) == (kw, cmi2(; kw...))
    end
    @test nmi2(verts=s.verts[1:2], plates=[P3(1, 2, 1)]) === :badvertexcount
    @test nmi2(plates=P3[]) === :badplatecount
    @test nmi2(plates=[P3(1, 2, length(s.verts) + 1)]) === :badvertexindex
    @test neo_kind(() -> dskrb2(s.verts, s.plates, 2, zeros(10))) === :notsupported
    @test neo_kind(() -> dskrb2(s.verts, s.plates, 4, [-1.0, 0.1])) === :valueoutofrange
    @test neo_kind(() -> dskrb2(s.verts, s.plates, 4, [1.0, 1.0])) === :valueoutofrange

    function w02_kinds(; bounds=s.bounds, first=s.first, last=s.last, dclass=s.dclass,
            corsys=s.corsys, frame=s.frame, corpar=zeros(10), plates=s.plates, mn=mn3, mx=mx3,
            ixd=spaixd, ixi=spaixi)
        pc = joinpath(dir, "wc.bds")
        pn = joinpath(dir, "wn.bds")
        rm(pc; force=true)
        rm(pn; force=true)
        ch = c_dskopn(pc, "x", 0)
        ck = c_kind(() -> c_dskw02(ch, s, mn, mx, ixd, Cint.(ixi); bounds, first, last, dclass,
            corsys, frame, corpar, plates))
        c_dskcls(ch)
        nh = dskopn(pn, "x", 0)
        nk = neo_kind(() -> dskw02(wctx, nh, s.center, s.surfid, dclass, frame, corsys, corpar,
            bounds..., mn, mx, first, last, s.verts, plates, ixd, ixi))
        dskcls(nh)
        # After a failed call the file has no segment in both engines.
        same = c_das_data(pc)[1] == neo_das_data(pn)[1]
        return nk, ck, same
    end
    tp = 2pi
    cases = Any[(frame="NO_SUCH_FRAME",), (first=1.0, last=1.0), (first=2.0, last=1.0),
        (bounds=(0.5, 0.5, -1.0, 1.0),), (bounds=(-7.0, 1.0, -1.0, 1.0),),
        (bounds=(tp, 1.0, -1.0, 1.0),), (bounds=(0.0, 7.0, -1.0, 1.0),),
        (bounds=(0.0, -tp, -1.0, 1.0),), (bounds=(-tp, tp, -1.0, 1.0),),
        (bounds=(-1.0, 1.0, -2.0, 1.0),), (bounds=(-1.0, 1.0, pi / 2, pi / 2),),
        (bounds=(-1.0, 1.0, -1.0, 2.0),), (bounds=(-1.0, 1.0, 1.0, -1.0),),
        (bounds=(-1.0, 1.0, -1.0, -pi / 2),), (bounds=(1.0, -1.0, -1.0, 1.0),),
        (bounds=(0.0, tp, 0.0, 1.0),), (bounds=(-tp, 0.0, -pi / 2, 0.0),),
        (mn=-1.0,), (mx=0.0, mn=0.0), (mn=5.0, mx=4.0), (mn=3.0, mx=3.0),
        (corsys=4, corpar=[0.0, 0.1, 0, 0, 0, 0, 0, 0, 0, 0]),
        (corsys=4, corpar=[3.0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0]),
        (corsys=4, corpar=[3.0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0], mn=-10.0),
        (corsys=4, corpar=[3.0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0], mn=-1.0, mx=1.0),
        (corsys=3, bounds=(1.0, -1.0, -1.0, 1.0)), (corsys=3, bounds=(-1.0, 1.0, 1.0, 1.0)),
        (corsys=3, mn=1.0, mx=1.0), (corsys=3, bounds=(-9.0, 9.0, -9.0, 9.0), mn=-9.0, mx=9.0),
        (corsys=2,), (corsys=0,), (dclass=0,), (dclass=3,),
        (plates=vcat(s.plates[1:end-1], [P3(1, 2, 0)]),),
        (plates=vcat(s.plates[1:end-1], [P3(1, 2, length(s.verts) + 1)]),),
        (ixi=(x = copy(spaixi); x[1] = 0; x),), (ixi=(x = copy(spaixi); x[4] = 0; x),),
        (ixi=(x = copy(spaixi); x[4] = 5; x),), (ixi=(x = copy(spaixi); x[4] = 64; x),),
        (ixi=(x = copy(spaixi); x[1:3] .= 1000; x),), (ixi=(x = copy(spaixi); x[1:3] .= 400; x[4] = 1; x),)]
    nerr = 0
    for kw in cases
        nk, ck, same = w02_kinds(; kw...)
        @test (kw, nk) == (kw, ck)
        @test same
        nerr += nk !== nothing
    end
    @test nerr >= 30
    h = dasopr(files[:a][1])
    @test neo_kind(() -> dskw02(wctx, h, 499, 1, 1, "IAU_MARS", 1, zeros(10), s.bounds..., mn3, mx3,
        0.0, 1.0, s.verts, s.plates, spaixd, spaixi)) === :dasinvalidaccess
    dascls(h)
    @test neo_kind(() -> dskopn(files[:a][1], "x", 0)) === :fileopenfailed
    @test neo_kind(() -> dskopn(joinpath(dir, "neg.bds"), "x", -1)) === :badrecordcount
    return nothing
end

"Part of group 12: (g) allocation and inference."
function dsk_part_13(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    iau = Int(Allspice.namfrm(ctx, "IAU_MARS"))
    tilt = Int(Allspice.namfrm(ctx, "MARS_TILT"))
    pck2 = Int(Allspice.namfrm(ctx, "MARS_PCK2"))
    v = V3(10.0, 1.0, 2.0)
    miss = V3(0.0, 0.0, 1.0)
    srf = [1, 2, 3, 4, 5, 6, 7]
    for fixfid in (iau, tilt, pck2), (lst, et) in (((), 0.0), (srf, 75.0), ((2, 5), 0.0)),
        (vtx, d) in ((v, -v), (v, miss), (V3(0.1, 0.2, 0.3), miss))
        @test dsk_warm(alloc_dskxsi, ctx, lst, et, fixfid, vtx, d) == 0
    end
    n = 256
    vtx = [dsk_ray(rng, 1, shape_e)[1] for _ in 1:n]
    dirs = [-x + V3(randn(rng, 3)) for x in vtx]
    xpt = Vector{V3}(undef, n)
    fnd = Vector{Bool}(undef, n)
    for fixfid in (iau, tilt)
        @test dsk_warm(alloc_dskxv!, xpt, fnd, ctx, srf, 75.0, fixfid, vtx, dirs) == 0
    end
    seg = dskxsi(ctx, false, 499, [2], 0.0, iau, v, -v)[2]
    @test dsk_warm(alloc_dskx02, ks, seg, v, -v) == 0
    @test dsk_warm(alloc_dskx02, ks, seg, v, miss) == 0
    R = Tuple{V3,DSKSegment,Int}
    @test (@inferred Union{Nothing,R} dskxsi(ctx, false, 499, srf, 0.0, iau, v, -v)) isa R
    @test (@inferred Union{Nothing,R} dskxsi(ctx, false, "MARS", (), 0.0, "IAU_MARS", v, -v)) isa R
    @test (@inferred Union{Nothing,R} dskxsi(ctx, false, 499, (), 0.0, iau, v, miss)) === nothing
    @test (@inferred Allspice._dsk_bux(ks, Allspice._DSK_DEFAULT_TOLS, Int32(499), srf, 0.0,
        Int32(iau), v, -v)) isa Tuple{V3,Int,Int,Bool}
    file = ks.dskfiles[seg.file]
    @test (@inferred Allspice._dsk_x02(file.ints, file.dps, seg, Allspice._DSK_DEFAULT_TOLS,
        v, -v)) isa Tuple{Int,V3,Bool}
    @test (@inferred dskxv!(xpt, fnd, ctx, false, 499, srf, 0.0, iau, vtx, dirs)) isa Tuple{Vector{V3},Vector{Bool}}
    @test (@inferred dskxv(ctx, false, 499, srf, 0.0, iau, vtx, dirs)) isa Tuple{Vector{V3},Vector{Bool}}
    @test (@inferred dskgtl(ctx, 1)) isa Float64
    @test (@inferred pltnp(v, vtx[1], vtx[2], vtx[3])) isa Tuple{V3,Float64}
    @test (@inferred pltnrm(vtx[1], vtx[2], vtx[3])) isa V3
    @test isbitstype(DSKSegment) && isbitstype(DLADescr) && isbitstype(DSKDescr)
    return nothing
end

"Part of group 12: (g) time per intercept."
function dsk_part_14(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    iau = Int(Allspice.namfrm(ctx, "IAU_MARS"))
    n = nsample(500, 20_000)
    vtx = [dsk_ray(rng, 1, shape_e)[1] for _ in 1:n]
    dirs = [V3(randn(rng, 3)) * 2.2 - x for x in vtx]
    xpt = Vector{V3}(undef, n)
    fnd = Vector{Bool}(undef, n)
    fv = flat(vtx); fd = flat(dirs); cx = zeros(3, n); cf = zeros(Cint, n)
    function c_time(srf, fixref)
        lst = collect(Cint, srf)
        best = Inf
        for _ in 1:3
            t = @elapsed ccall((:dskxv_c, DSKLIB), Cvoid, (Cint, Cstring, Cint, Ptr{Cint}, Float64,
                    Cstring, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cint}), 0, "MARS",
                length(lst), lst, 0.0, fixref, n, fv, fd, cx, cf)
            best = min(best, t)
        end
        c_check()
        return best / n * 1e9
    end
    function n_time(srf, fixfid)
        best = Inf
        for _ in 1:3
            t = @elapsed dskxv!(xpt, fnd, ctx, false, 499, srf, 0.0, fixfid, vtx, dirs)
            best = min(best, t)
        end
        return best / n * 1e9
    end
    measured("  time per intercept (ns), best of 3 runs of $n rays, one thread:")
    measured("    case                                      Allspice    CSPICE   ratio")
    for (label, srf, fixref) in (("one segment, 3480 plates (surface 1)", [1], "IAU_MARS"),
            ("one segment, 5040 plates (surface 2)", [2], "IAU_MARS"),
            ("9 segments, all surfaces", Int[], "IAU_MARS"),
            ("9 segments, other frame (MARS_TILT)", Int[], "MARS_TILT"))
        tn = n_time(srf, Int(Allspice.namfrm(ctx, fixref)))
        tc = c_time(srf, fixref)
        measured(@sprintf("    %-40s %9.0f %9.0f %7.2f", label, tn, tc, tn / tc))
        @test tn < 20 * tc           # a coarse bound; the measured values are in the report
    end
    return nothing
end

"Part of group 12: (g) many threads, while another task loads and unloads a DSK."
function dsk_part_15(env)
    (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments) = env
    (; loadlist, ctx, ks, shape_e, shape_r, shape_w) = env
    if Threads.nthreads() < 2
        skip_notice("12. DSK threads", "Threads.nthreads() == 1")
    else
        iau = Int(Allspice.namfrm(ctx, "IAU_MARS"))
        base = [tk, spk, files[:a][1]]
        extra = files[:b][1]
        n = 2000
        vtx = [dsk_ray(rng, 1, shape_e)[1] for _ in 1:n]
        dirs = [V3(randn(rng, 3)) * 2.2 - x for x in vtx]
        ctxA = neo_context(base)
        ctxB = neo_context([base..., extra])
        q(c, i) = dskxsi(c, false, 499, (), 60.0, iau, vtx[i], dirs[i])
        key(r) = r === nothing ? nothing : (r[1], r[2].dskdsc.surfce, r[3])
        refA = [key(q(ctxA, i)) for i in 1:n]
        refB = [key(q(ctxB, i)) for i in 1:n]
        @test count(refA .!= refB) > n ÷ 2      # the two snapshots give different answers
        shared = neo_context(base)
        stop = Threads.Atomic{Bool}(false)
        bad = Threads.Atomic{Int}(0)
        seenA = Threads.Atomic{Int}(0)
        seenB = Threads.Atomic{Int}(0)
        total = Threads.Atomic{Int}(0)
        writer = Threads.@spawn begin
            cycles = 0
            while !stop[]
                furnsh(shared, extra)
                yield()
                unload(shared, extra)
                cycles += 1
            end
            cycles
        end
        readers = [Threads.@spawn begin
            for rep in 1:nsample(8, 40), i in 1:n
                r = try
                    key(q(shared, i))
                catch
                    Threads.atomic_add!(bad, 1)
                    continue
                end
                if r == refA[i]
                    Threads.atomic_add!(seenA, 1)
                elseif r == refB[i]
                    Threads.atomic_add!(seenB, 1)
                else
                    Threads.atomic_add!(bad, 1)
                end
                Threads.atomic_add!(total, 1)
            end
        end for _ in 1:max(2, Threads.nthreads() - 1)]
        foreach(wait, readers)
        stop[] = true
        cycles = fetch(writer)
        measured("  DSK threads: $(total[]) intercepts on $(length(readers)) tasks during " *
                 "$cycles load/unload cycles; $(seenA[]) saw the set without the file, " *
                 "$(seenB[]) the set with it, $(bad[]) saw neither")
        @test bad[] == 0
        @test cycles > 0 && seenA[] > 0
        # A pure parallel loop gives the serial results.
        par = Vector{Any}(undef, n)
        Threads.@threads for i in 1:n
            par[i] = key(q(ctxB, i))
        end
        @test par == refB
    end
    return nothing
end

group("12. DSK type 2 against CSPICE") do
    rng = MersenneTwister(SEED + 12)
    dir = mktempdir()
    tk = dsk_text_kernel(dir)
    spk = dsk_offset_spk(dir)
    specs_a, specs_b, specs_c = dsk_specs(rng)
    files = Dict{Symbol,NTuple{2,String}}()
    oracle_load([tk, spk])               # the writers need the frame names
    wctx = neo_context([tk, spk])
    comments = ["DSK test file", "second comment line"]
    env = (; rng, dir, tk, spk, specs_a, specs_b, specs_c, files, wctx, comments)

    @apart @testset "(a) spatial index, bounds and files of the two writers" begin
        dsk_part_1(env)
    end

    @apart @testset "(b) access calls" begin
        dsk_part_2(env)
    end

    @apart @testset "(c) dskx02 on single segments" begin
        dsk_part_3(env)
    end

    # All files in both engines.
    loadlist = [tk, spk, files[:a][1], files[:b][1], files[:c][1]]
    oracle_load(loadlist)
    ctx = neo_context(loadlist)
    ks = Allspice.kernelset(ctx)
    shape_e = (specs_a[1].verts, specs_a[1].plates)
    shape_r = (specs_a[2].verts, specs_a[2].plates)
    shape_w = (specs_b[2].verts, specs_b[2].plates)

    env = (; env..., loadlist, ctx, ks, shape_e, shape_r, shape_w)

    @apart @testset "(d) dskxsi on the loaded files" begin
        dsk_part_4(env)
    end

    @apart @testset "(d) dskxv, argument forms and coverage in time" begin
        dsk_part_5(env)
    end

    @apart @testset "(e) plate utilities" begin
        dsk_part_6(env)
    end

    oracle_load(loadlist)
    @apart @testset "(e) surface names and codes" begin
        dsk_part_7(env)
    end

    @apart @testset "(e) latsrf and srfnrm" begin
        dsk_part_8(env)
    end

    @apart @testset "(f) tolerances" begin
        dsk_part_9(env)
    end

    @apart @testset "(f) kernel management" begin
        dsk_part_10(env)
    end

    @apart @testset "(f) input checks of the intercept and the access calls" begin
        dsk_part_11(env)
    end

    @apart @testset "(f) input checks of the writer" begin
        dsk_part_12(env)
    end

    oracle_load(loadlist)
    @apart @testset "(g) allocation and inference" begin
        dsk_part_13(env)
    end

    @apart @testset "(g) time per intercept" begin
        dsk_part_14(env)
    end

    @apart @testset "(g) many threads, while another task loads and unloads a DSK" begin
        dsk_part_15(env)
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
