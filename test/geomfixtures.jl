# Fixtures of the geometry groups: direct calls into CSPICE, an instrument kernel that
# the tests write, DSK files for the Moon that the tests write with the NeoSpice DSK
# writer, and synthetic bodies for the occultation cases.
#
# CSPICE is called from the main task only. Each direct call ends with `c_check`, which
# turns a CSPICE error into a Julia error and resets the CSPICE error state; `c_kind`
# gives the kind of that error.

isdefined(@__MODULE__, :DSKLIB) || include("dskfixtures.jl")

const GEOLIB = SPICE.libcspice
const GV = Vector{Float64}

# ---------------------------------------------------------------------------------------
# Direct calls into CSPICE
# ---------------------------------------------------------------------------------------

function g_subpnt(method, target, et, fixref, abcorr, obsrvr)
    sp = zeros(3); te = Ref(0.0); sv = zeros(3)
    ccall((:subpnt_c, GEOLIB), Cvoid, (Cstring, Cstring, Float64, Cstring, Cstring, Cstring,
            Ptr{Float64}, Ref{Float64}, Ptr{Float64}), method, target, et, fixref, abcorr, obsrvr, sp, te, sv)
    c_check()
    return sp, te[], sv
end

function g_subslr(method, target, et, fixref, abcorr, obsrvr)
    sp = zeros(3); te = Ref(0.0); sv = zeros(3)
    ccall((:subslr_c, GEOLIB), Cvoid, (Cstring, Cstring, Float64, Cstring, Cstring, Cstring,
            Ptr{Float64}, Ref{Float64}, Ptr{Float64}), method, target, et, fixref, abcorr, obsrvr, sp, te, sv)
    c_check()
    return sp, te[], sv
end

"`sincpt_c`: `(spoint, trgepc, srfvec)` or `nothing`."
function g_sincpt(method, target, et, fixref, abcorr, obsrvr, dref, dvec)
    sp = zeros(3); te = Ref(0.0); sv = zeros(3); found = Ref{Cint}(0)
    ccall((:sincpt_c, GEOLIB), Cvoid, (Cstring, Cstring, Float64, Cstring, Cstring, Cstring, Cstring,
            Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Cint}), method, target, et,
        fixref, abcorr, obsrvr, dref, collect(Float64, dvec), sp, te, sv, found)
    c_check()
    return found[] != 0 ? (sp, te[], sv) : nothing
end

function g_ilumin(method, target, et, fixref, abcorr, obsrvr, spoint)
    te = Ref(0.0); sv = zeros(3); a = Ref(0.0); b = Ref(0.0); c = Ref(0.0)
    ccall((:ilumin_c, GEOLIB), Cvoid, (Cstring, Cstring, Float64, Cstring, Cstring, Cstring,
            Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        method, target, et, fixref, abcorr, obsrvr, collect(Float64, spoint), te, sv, a, b, c)
    c_check()
    return te[], sv, a[], b[], c[]
end

function g_illumg(method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
    te = Ref(0.0); sv = zeros(3); a = Ref(0.0); b = Ref(0.0); c = Ref(0.0)
    ccall((:illumg_c, GEOLIB), Cvoid, (Cstring, Cstring, Cstring, Float64, Cstring, Cstring, Cstring,
            Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        method, target, ilusrc, et, fixref, abcorr, obsrvr, collect(Float64, spoint), te, sv, a, b, c)
    c_check()
    return te[], sv, a[], b[], c[]
end

function g_illumf(method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
    te = Ref(0.0); sv = zeros(3); a = Ref(0.0); b = Ref(0.0); c = Ref(0.0)
    vis = Ref{Cint}(0); lit = Ref{Cint}(0)
    ccall((:illumf_c, GEOLIB), Cvoid, (Cstring, Cstring, Cstring, Float64, Cstring, Cstring, Cstring,
            Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
            Ref{Cint}, Ref{Cint}), method, target, ilusrc, et, fixref, abcorr, obsrvr,
        collect(Float64, spoint), te, sv, a, b, c, vis, lit)
    c_check()
    return te[], sv, a[], b[], c[], vis[] != 0, lit[] != 0
end

function g_phaseq(et, target, illmn, obsrvr, abcorr)
    r = ccall((:phaseq_c, GEOLIB), Float64, (Float64, Cstring, Cstring, Cstring, Cstring), et, target,
        illmn, obsrvr, abcorr)
    c_check()
    return r
end

function g_occult(t1, s1, f1, t2, s2, f2, abcorr, obsrvr, et)
    id = Ref{Cint}(0)
    ccall((:occult_c, GEOLIB), Cvoid, (Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring,
            Cstring, Float64, Ref{Cint}), t1, s1, f1, t2, s2, f2, abcorr, obsrvr, et, id)
    c_check()
    return Int(id[])
end

"`limbpt_c`: `(npts, points, epochs, tangts)` with the points as the columns of a matrix."
function g_limbpt(method, target, et, fixref, abcorr, corloc, obsrvr, refvec, rolstp, ncuts,
        schstp, soltol, maxn)
    m = max(maxn, 1)
    npts = zeros(Cint, m); pts = zeros(3, m); ep = zeros(m); tg = zeros(3, m)
    ccall((:limbpt_c, GEOLIB), Cvoid, (Cstring, Cstring, Float64, Cstring, Cstring, Cstring, Cstring,
            Ptr{Float64}, Float64, Cint, Float64, Float64, Cint, Ptr{Cint}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}), method, target, et, fixref, abcorr, corloc, obsrvr,
        collect(Float64, refvec), rolstp, ncuts, schstp, soltol, maxn, npts, pts, ep, tg)
    c_check()
    n = sum(npts[1:ncuts])
    return Int.(npts[1:ncuts]), pts[:, 1:n], ep[1:n], tg[:, 1:n]
end

"`termpt_c`: `(npts, points, epochs, trmvcs)` with the points as the columns of a matrix."
function g_termpt(method, ilusrc, target, et, fixref, abcorr, corloc, obsrvr, refvec, rolstp,
        ncuts, schstp, soltol, maxn)
    m = max(maxn, 1)
    npts = zeros(Cint, m); pts = zeros(3, m); ep = zeros(m); tg = zeros(3, m)
    ccall((:termpt_c, GEOLIB), Cvoid, (Cstring, Cstring, Cstring, Float64, Cstring, Cstring, Cstring,
            Cstring, Ptr{Float64}, Float64, Cint, Float64, Float64, Cint, Ptr{Cint}, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}), method, ilusrc, target, et, fixref, abcorr, corloc, obsrvr,
        collect(Float64, refvec), rolstp, ncuts, schstp, soltol, maxn, npts, pts, ep, tg)
    c_check()
    n = sum(npts[1:max(ncuts, 0)])
    return Int.(npts[1:max(ncuts, 0)]), pts[:, 1:n], ep[1:n], tg[:, 1:n]
end

function g_tangpt(method, target, et, fixref, abcorr, corloc, obsrvr, dref, dvec)
    tp = zeros(3); alt = Ref(0.0); rg = Ref(0.0); sp = zeros(3); te = Ref(0.0); sv = zeros(3)
    ccall((:tangpt_c, GEOLIB), Cvoid, (Cstring, Cstring, Float64, Cstring, Cstring, Cstring, Cstring,
            Cstring, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Float64}, Ptr{Float64},
            Ref{Float64}, Ptr{Float64}), method, target, et, fixref, abcorr, corloc, obsrvr, dref,
        collect(Float64, dvec), tp, alt, rg, sp, te, sv)
    c_check()
    return tp, alt[], rg[], sp, te[], sv
end

function g_fovtrg(inst, target, tshape, tframe, abcorr, obsrvr, et)
    vis = Ref{Cint}(0)
    ccall((:fovtrg_c, GEOLIB), Cvoid, (Cstring, Cstring, Cstring, Cstring, Cstring, Cstring,
            Ref{Float64}, Ref{Cint}), inst, target, tshape, tframe, abcorr, obsrvr, et, vis)
    c_check()
    return vis[] != 0
end

function g_fovray(inst, raydir, rframe, abcorr, obsrvr, et)
    vis = Ref{Cint}(0)
    ccall((:fovray_c, GEOLIB), Cvoid, (Cstring, Ptr{Float64}, Cstring, Cstring, Cstring,
            Ref{Float64}, Ref{Cint}), inst, collect(Float64, raydir), rframe, abcorr, obsrvr, et, vis)
    c_check()
    return vis[] != 0
end

"`getfov_c`: `(shape, frame, bsight, bounds)` with the boundary vectors as the columns of a matrix."
function g_getfov(instid, room=100)
    shape = zeros(UInt8, 64); frame = zeros(UInt8, 64); bs = zeros(3); n = Ref{Cint}(0)
    bounds = zeros(3, max(room, 1))
    ccall((:getfov_c, GEOLIB), Cvoid, (Cint, Cint, Cint, Cint, Ptr{UInt8}, Ptr{UInt8}, Ptr{Float64},
            Ref{Cint}, Ptr{Float64}), instid, room, 64, 64, shape, frame, bs, n, bounds)
    c_check()
    return unsafe_string(pointer(shape)), unsafe_string(pointer(frame)), bs, bounds[:, 1:n[]]
end

function g_getfvn(inst, room=100)
    shape = zeros(UInt8, 64); frame = zeros(UInt8, 64); bs = zeros(3); n = Ref{Cint}(0)
    bounds = zeros(3, max(room, 1))
    ccall((:getfvn_c, GEOLIB), Cvoid, (Cstring, Cint, Cint, Cint, Ptr{UInt8}, Ptr{UInt8}, Ptr{Float64},
            Ref{Cint}, Ptr{Float64}), inst, room, 64, 64, shape, frame, bs, n, bounds)
    c_check()
    return unsafe_string(pointer(shape)), unsafe_string(pointer(frame)), bs, bounds[:, 1:n[]]
end

function g_srfrec(body, lon, lat)
    r = zeros(3)
    ccall((:srfrec_c, GEOLIB), Cvoid, (Cint, Float64, Float64, Ptr{Float64}), body, lon, lat, r)
    c_check()
    return r
end

function g_spkcpo(target, et, outref, refloc, abcorr, obspos, obsctr, obsref)
    st = zeros(6); lt = Ref(0.0)
    ccall((:spkcpo_c, GEOLIB), Cvoid, (Cstring, Float64, Cstring, Cstring, Cstring, Ptr{Float64},
            Cstring, Cstring, Ptr{Float64}, Ref{Float64}), target, et, outref, refloc, abcorr,
        collect(Float64, obspos), obsctr, obsref, st, lt)
    c_check()
    return st, lt[]
end

function g_spkcpt(trgpos, trgctr, trgref, et, outref, refloc, abcorr, obsrvr)
    st = zeros(6); lt = Ref(0.0)
    ccall((:spkcpt_c, GEOLIB), Cvoid, (Ptr{Float64}, Cstring, Cstring, Float64, Cstring, Cstring,
            Cstring, Cstring, Ptr{Float64}, Ref{Float64}), collect(Float64, trgpos), trgctr, trgref,
        et, outref, refloc, abcorr, obsrvr, st, lt)
    c_check()
    return st, lt[]
end

function g_spkcvo(target, et, outref, refloc, abcorr, obssta, obsepc, obsctr, obsref)
    st = zeros(6); lt = Ref(0.0)
    ccall((:spkcvo_c, GEOLIB), Cvoid, (Cstring, Float64, Cstring, Cstring, Cstring, Ptr{Float64},
            Float64, Cstring, Cstring, Ptr{Float64}, Ref{Float64}), target, et, outref, refloc,
        abcorr, collect(Float64, obssta), obsepc, obsctr, obsref, st, lt)
    c_check()
    return st, lt[]
end

function g_spkcvt(trgsta, trgepc, trgctr, trgref, et, outref, refloc, abcorr, obsrvr)
    st = zeros(6); lt = Ref(0.0)
    ccall((:spkcvt_c, GEOLIB), Cvoid, (Ptr{Float64}, Float64, Cstring, Cstring, Float64, Cstring,
            Cstring, Cstring, Cstring, Ptr{Float64}, Ref{Float64}), collect(Float64, trgsta), trgepc,
        trgctr, trgref, et, outref, refloc, abcorr, obsrvr, st, lt)
    c_check()
    return st, lt[]
end

function g_azlcpo(method, target, et, abcorr, azccw, elplsz, obspos, obsctr, obsref)
    st = zeros(6); lt = Ref(0.0)
    ccall((:azlcpo_c, GEOLIB), Cvoid, (Cstring, Cstring, Float64, Cstring, Cint, Cint, Ptr{Float64},
            Cstring, Cstring, Ptr{Float64}, Ref{Float64}), method, target, et, abcorr, azccw, elplsz,
        collect(Float64, obspos), obsctr, obsref, st, lt)
    c_check()
    return st, lt[]
end

function g_azlrec(range, az, el, azccw, elplsz)
    r = zeros(3)
    ccall((:azlrec_c, GEOLIB), Cvoid, (Float64, Float64, Float64, Cint, Cint, Ptr{Float64}), range,
        az, el, azccw, elplsz, r)
    c_check()
    return r
end

function g_recazl(rectan, azccw, elplsz)
    a = Ref(0.0); b = Ref(0.0); c = Ref(0.0)
    ccall((:recazl_c, GEOLIB), Cvoid, (Ptr{Float64}, Cint, Cint, Ref{Float64}, Ref{Float64},
            Ref{Float64}), collect(Float64, rectan), azccw, elplsz, a, b, c)
    c_check()
    return a[], b[], c[]
end

"`dazldr_c`: the Jacobian as a Julia matrix (row `i`, column `j`)."
function g_dazldr(x, y, z, azccw, elplsz)
    j = zeros(3, 3)
    ccall((:dazldr_c, GEOLIB), Cvoid, (Float64, Float64, Float64, Cint, Cint, Ptr{Float64}), x, y, z,
        azccw, elplsz, j)
    c_check()
    return permutedims(j)
end

function g_drdazl(range, az, el, azccw, elplsz)
    j = zeros(3, 3)
    ccall((:drdazl_c, GEOLIB), Cvoid, (Float64, Float64, Float64, Cint, Cint, Ptr{Float64}), range,
        az, el, azccw, elplsz, j)
    c_check()
    return permutedims(j)
end

"`pdpool_c` and `pdpool` of NeoSpice with the same values."
function both_pdpool(ctx, name, values)
    v = collect(Float64, values)
    ccall((:pdpool_c, GEOLIB), Cvoid, (Cstring, Cint, Ptr{Float64}), name, length(v), v)
    c_check()
    pdpool(ctx, name, v)
    return nothing
end

"`pcpool_c` and `pcpool` of NeoSpice with one string."
function both_pcpool(ctx, name, value::AbstractString)
    buf = zeros(UInt8, 81)
    copyto!(buf, codeunits(value))
    ccall((:pcpool_c, GEOLIB), Cvoid, (Cstring, Cint, Cint, Ptr{UInt8}), name, 1, 81, buf)
    c_check()
    pcpool(ctx, name, [String(value)])
    return nothing
end

"Remove the pool variable `name` from CSPICE and from the NeoSpice context."
function both_dvpool(ctx, name)
    ccall((:dvpool_c, GEOLIB), Cvoid, (Cstring,), name)
    c_check()
    dvpool(ctx, name)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Comparison helpers
# ---------------------------------------------------------------------------------------

"Relative difference of two vectors, `|a - b| / |b|`."
grel(a, b) = norm(collect(a) .- collect(b)) / max(norm(collect(b)), floatmin())

const GEO_CORRS = ["NONE", "LT", "LT+S", "CN", "CN+S", "XLT", "XLT+S", "XCN", "XCN+S"]
const GEO_BODIES = ["SUN", "MERCURY", "VENUS", "EARTH", "MOON"]
const GEO_IDS = Dict("SUN" => 10, "MERCURY" => 199, "VENUS" => 299, "EARTH" => 399, "MOON" => 301)
"Random epochs of the geometry groups: 1990 to 2010 (ET seconds), inside all kernels."
geo_epoch(rng) = (rand(rng) - 0.5) * 6.0e8

"The error kinds of CSPICE and of NeoSpice for the same call."
kinds(cf, nf) = (c_kind(cf), neo_kind(nf))

# ---------------------------------------------------------------------------------------
# Instrument kernel
# ---------------------------------------------------------------------------------------

const GEO_IK_PATH = Ref("")

"""
    geo_instrument_kernel() -> String

Write (once per process) a text kernel with the frame `GEO_CAM` (class 4, fixed offset
from J2000, center Earth; the tests change its matrix with `pdpool`), the frame
`GEO_SITE` (class 4, fixed offset from `IAU_EARTH`) and the instruments

| ID      | name         | shape     | specification                    |
|---------|--------------|-----------|----------------------------------|
| -999001 | GEO_CIRCLE   | CIRCLE    | corners                          |
| -999002 | GEO_ELLIPSE  | ELLIPSE   | corners                          |
| -999003 | GEO_RECT     | RECTANGLE | corners                          |
| -999004 | GEO_POLY     | POLYGON   | corners, 6 vertices, not convex  |
| -999005 | GEO_CIRCLE_A | CIRCLE    | angles, degrees                  |
| -999006 | GEO_ELLIPS_A | ELLIPSE   | angles, radians                  |
| -999007 | GEO_RECT_A   | RECTANGLE | angles, arc minutes              |
| -999008 | GEO_WIDE     | POLYGON   | corners, one vertex more than 90 degrees from the mean axis |
| -999009 | GEO_FLAT     | POLYGON   | corners, 90 degrees from the axis (too wide) |

and return its path.
"""
function geo_instrument_kernel()
    isfile(GEO_IK_PATH[]) && return GEO_IK_PATH[]
    path = joinpath(mktempdir(), "geo_instruments.ti")
    write(path, """
    KPL/IK

    Instruments and frames of the geometry tests (written by the tests).

    \\begindata

       FRAME_GEO_CAM                  = -999100
       FRAME_-999100_NAME             = 'GEO_CAM'
       FRAME_-999100_CLASS            = 4
       FRAME_-999100_CLASS_ID         = -999100
       FRAME_-999100_CENTER           = 399
       TKFRAME_-999100_SPEC           = 'MATRIX'
       TKFRAME_-999100_RELATIVE       = 'J2000'
       TKFRAME_-999100_MATRIX         = ( 1 0 0  0 1 0  0 0 1 )

       FRAME_GEO_SITE                 = -999200
       FRAME_-999200_NAME             = 'GEO_SITE'
       FRAME_-999200_CLASS            = 4
       FRAME_-999200_CLASS_ID         = -999200
       FRAME_-999200_CENTER           = 399
       TKFRAME_-999200_SPEC           = 'ANGLES'
       TKFRAME_-999200_RELATIVE       = 'IAU_EARTH'
       TKFRAME_-999200_ANGLES         = ( 0.7 -0.4 1.1 )
       TKFRAME_-999200_AXES           = ( 3 2 1 )
       TKFRAME_-999200_UNITS          = 'RADIANS'

       NAIF_BODY_NAME += ( 'GEO_CIRCLE', 'GEO_ELLIPSE', 'GEO_RECT', 'GEO_POLY',
                           'GEO_CIRCLE_A', 'GEO_ELLIPS_A', 'GEO_RECT_A', 'GEO_WIDE',
                           'GEO_FLAT' )
       NAIF_BODY_CODE += ( -999001, -999002, -999003, -999004,
                           -999005, -999006, -999007, -999008, -999009 )

       INS-999001_FOV_SHAPE           = 'CIRCLE'
       INS-999001_FOV_FRAME           = 'GEO_CAM'
       INS-999001_BORESIGHT           = ( 0.0 0.0 2.0 )
       INS-999001_FOV_BOUNDARY_CORNERS = ( 0.02 0.0 1.0 )

       INS-999002_FOV_SHAPE           = 'ellipse'
       INS-999002_FOV_FRAME           = 'GEO_CAM'
       INS-999002_BORESIGHT           = ( 0.0 0.0 1.0 )
       INS-999002_FOV_BOUNDARY        = ( 0.03 0.0 1.0   0.0 0.015 1.0 )

       INS-999003_FOV_SHAPE           = 'RECTANGLE'
       INS-999003_FOV_FRAME           = 'GEO_CAM'
       INS-999003_BORESIGHT           = ( 0.0 0.0 1.0 )
       INS-999003_FOV_CLASS_SPEC      = 'CORNERS'
       INS-999003_FOV_BOUNDARY_CORNERS = (  0.03  0.02 1.0   -0.03  0.02 1.0
                                           -0.03 -0.02 1.0    0.03 -0.02 1.0 )

       INS-999004_FOV_SHAPE           = 'POLYGON'
       INS-999004_FOV_FRAME           = 'GEO_CAM'
       INS-999004_BORESIGHT           = ( 0.001 0.002 1.0 )
       INS-999004_FOV_BOUNDARY_CORNERS = (  0.04  0.00 1.0    0.01  0.01 1.0
                                            0.00  0.04 1.0   -0.03  0.01 1.0
                                           -0.02 -0.03 1.0    0.02 -0.02 1.0 )

       INS-999005_FOV_SHAPE           = 'CIRCLE'
       INS-999005_FOV_FRAME           = 'GEO_CAM'
       INS-999005_BORESIGHT           = ( 0.0 0.0 1.0 )
       INS-999005_FOV_CLASS_SPEC      = 'ANGLES'
       INS-999005_FOV_REF_VECTOR      = ( 1.0 0.0 0.0 )
       INS-999005_FOV_REF_ANGLE       = 1.5
       INS-999005_FOV_ANGLE_UNITS     = 'DEGREES'

       INS-999006_FOV_SHAPE           = 'ELLIPSE'
       INS-999006_FOV_FRAME           = 'GEO_CAM'
       INS-999006_BORESIGHT           = ( 0.0 0.1 3.0 )
       INS-999006_FOV_CLASS_SPEC      = 'angles'
       INS-999006_FOV_REF_VECTOR      = ( 1.0 0.3 0.0 )
       INS-999006_FOV_REF_ANGLE       = 0.03
       INS-999006_FOV_CROSS_ANGLE     = 0.012
       INS-999006_FOV_ANGLE_UNITS     = 'RADIANS'

       INS-999007_FOV_SHAPE           = 'RECTANGLE'
       INS-999007_FOV_FRAME           = 'GEO_CAM'
       INS-999007_BORESIGHT           = ( 0.0 0.0 1.0 )
       INS-999007_FOV_CLASS_SPEC      = 'ANGLES'
       INS-999007_FOV_REF_VECTOR      = ( 0.2 1.0 0.5 )
       INS-999007_FOV_REF_ANGLE       = 90.0
       INS-999007_FOV_CROSS_ANGLE     = 55.0
       INS-999007_FOV_ANGLE_UNITS     = 'ARCMINUTES'

       INS-999008_FOV_SHAPE           = 'POLYGON'
       INS-999008_FOV_FRAME           = 'GEO_CAM'
       INS-999008_BORESIGHT           = ( 0.0 0.0 1.0 )
       INS-999008_FOV_BOUNDARY        = ( -1.0  0.0  0.05    1.0  0.1  0.3
                                           1.0  0.0  0.5     1.0 -0.1  0.3 )

       INS-999009_FOV_SHAPE           = 'POLYGON'
       INS-999009_FOV_FRAME           = 'GEO_CAM'
       INS-999009_BORESIGHT           = ( 0.0 0.0 1.0 )
       INS-999009_FOV_BOUNDARY        = (  1.0  0.0  1.0E-8    0.0  1.0  1.0E-8
                                          -1.0  0.0  1.0E-8    0.0 -1.0  1.0E-8 )

    \\begintext
    """)
    GEO_IK_PATH[] = path
    return path
end

const GEO_INSTRUMENTS = ["GEO_CIRCLE", "GEO_ELLIPSE", "GEO_RECT", "GEO_POLY", "GEO_CIRCLE_A",
    "GEO_ELLIPS_A", "GEO_RECT_A"]

"""
    point_camera!(ctx, zdir, twist)

Set the matrix of the frame `GEO_CAM` in CSPICE and in `ctx`, so that its +Z axis is
along `zdir` (J2000) and its +X axis is turned by `twist` radians about +Z.
"""
function point_camera!(ctx, zdir, twist)
    z = normalize(collect(Float64, zdir))
    h = abs(z[3]) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
    x0 = normalize(cross(h, z))
    y0 = cross(z, x0)
    x = cos(twist) * x0 + sin(twist) * y0
    y = cross(z, x)
    # TKFRAME matrix: the columns are the axes of GEO_CAM in J2000, stored by column.
    both_pdpool(ctx, "TKFRAME_-999100_MATRIX", vcat(x, y, z))
    return nothing
end

# ---------------------------------------------------------------------------------------
# DSK files for the Moon
# ---------------------------------------------------------------------------------------

const GEO_DSK_PATH = Ref("")

"""
    geo_moon_dsk() -> String

Write (once per process), with the DSK writer of NeoSpice, a DSK file for the Moon in
the frame `IAU_MOON`: surface 1 is a tessellated ellipsoid with the radii of the Moon
(120 by 60 grid), surface 2 a rough closed surface (radius 1737.4 km, ±1.5 percent,
90 by 45 grid). Return its path.
"""
function geo_moon_dsk()
    isfile(GEO_DSK_PATH[]) && return GEO_DSK_PATH[]
    rng = MersenneTwister(SEED + 77)
    ev, ep = ellipsoid_shape(1737.4, 1737.4, 1737.4, 120, 60)
    rv, rp = rough_shape(rng, 1737.4, 0.015, 90, 45)
    full = (-pi, pi, -pi / 2, pi / 2)
    specs = [DskSpec(301, 1, 2, "IAU_MOON", 1, Float64[], full, -1e10, 1e10, 3.0, 4, ev, ep),
        DskSpec(301, 2, 2, "IAU_MOON", 1, Float64[], full, -1e10, 1e10, 3.0, 4, rv, rp)]
    path = joinpath(mktempdir(), "geo_moon.bds")
    neo_write_dsk(path, specs; ctx=neo_context(Kernels.full_set()))
    GEO_DSK_PATH[] = path
    return path
end

"The kernel list of the geometry groups: the standard set, the instrument kernel and, with `dsk`, the Moon DSK."
geo_kernels(; dsk::Bool=false) = dsk ? [Kernels.full_set()..., geo_instrument_kernel(), geo_moon_dsk()] :
                                 [Kernels.full_set()..., geo_instrument_kernel()]

"Load the geometry kernels in CSPICE and in a new NeoSpice context, which is returned."
function geo_setup(; dsk::Bool=false)
    paths = geo_kernels(; dsk)
    oracle_load(paths)
    return neo_context(paths)
end
