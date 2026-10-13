# Fixtures of the DAS and DSK groups: direct calls into CSPICE, plate models that the
# tests generate (no DSK file is needed on the machine), and file comparison.
#
# CSPICE is called from the main task only. Each direct call ends with `c_check`, which
# turns a CSPICE error into a Julia error and resets the CSPICE error state.

const DSKLIB = SPICE.libcspice

"Short message of the pending CSPICE error as a lower-case symbol, or `nothing`; the error state is reset."
function c_pending_kind()
    ccall((:failed_c, DSKLIB), Cint, ()) != 0 || return nothing
    buf = zeros(UInt8, 1841)
    ccall((:getmsg_c, DSKLIB), Cvoid, (Cstring, Cint, Ptr{UInt8}), "SHORT", 1841, buf)
    short = unsafe_string(pointer(buf))
    ccall((:reset_c, DSKLIB), Cvoid, ())
    return Symbol(lowercase(replace(short, "SPICE(" => "", ")" => "")))
end

"Throw if CSPICE has a pending error."
function c_check()
    kind = c_pending_kind()
    kind === nothing || error("CSPICE error: $kind")
    return nothing
end

"The kind of the CSPICE error that `f()` raises, or `nothing`. `f` must use direct calls without `c_check`."
function c_kind(f)
    try
        f()
    catch e
        e isa ErrorException || rethrow()
        m = match(r"CSPICE error: (\w+)", e.msg)
        m === nothing && rethrow()
        return Symbol(m.captures[1])
    end
    return c_pending_kind()
end

"The kind of the `SpiceError` that `f()` throws, or `nothing`."
function neo_kind(f)
    try
        f()
    catch e
        e isa SpiceError || rethrow()
        return e.kind
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Plate models
# ---------------------------------------------------------------------------------------

const V3 = SVector{3,Float64}
const P3 = SVector{3,Int}

flat(v::Vector{V3}) = collect(reinterpret(Float64, v))
flat(p::Vector{P3}) = Cint[x for q in p for x in q]

# SPICE.jl converts a list of vectors with `hcat(list...)`, which compiles a new method
# for each list length; the calls below pass flat arrays.
"`pltar_c` for vertices `verts` and plates `plates`."
c_pltar(verts, plates) = (r = ccall((:pltar_c, DSKLIB), Float64, (Cint, Ptr{Float64}, Cint,
        Ptr{Cint}), length(verts), flat(verts), length(plates), flat(plates)); c_check(); r)

"`pltvol_c` for vertices `verts` and plates `plates`."
c_pltvol(verts, plates) = (r = ccall((:pltvol_c, DSKLIB), Float64, (Cint, Ptr{Float64}, Cint,
        Ptr{Cint}), length(verts), flat(verts), length(plates), flat(plates)); c_check(); r)

"`latsrf_c`: the surface points for a list of (longitude, latitude) pairs."
function c_latsrf(method, target, et, fixref, lonlat)
    n = length(lonlat)
    out = zeros(3, n)
    ccall((:latsrf_c, DSKLIB), Cvoid, (Cstring, Cstring, Float64, Cstring, Cint, Ptr{Float64},
            Ptr{Float64}), method, target, et, fixref, n, Float64[x for q in lonlat for x in q], out)
    c_check()
    return [V3(out[1, i], out[2, i], out[3, i]) for i in 1:n]
end

"`srfnrm_c`: the outward unit normals at a list of surface points."
function c_srfnrm(method, target, et, fixref, pts)
    n = length(pts)
    out = zeros(3, n)
    ccall((:srfnrm_c, DSKLIB), Cvoid, (Cstring, Cstring, Float64, Cstring, Cint, Ptr{Float64},
            Ptr{Float64}), method, target, et, fixref, n, Float64[x for q in pts for x in q], out)
    c_check()
    return [V3(out[1, i], out[2, i], out[3, i]) for i in 1:n]
end

"""
Closed surface r(lon, lat) on a grid of `nlon` longitudes and `nlat` latitude bands with
the two poles; the plates have outward normals.
"""
function latlon_shape(radius, nlon::Int, nlat::Int)
    verts = V3[]
    push!(verts, V3(0.0, 0.0, radius(0.0, pi / 2)))
    for i in 1:nlat-1, j in 0:nlon-1
        lat = pi / 2 - i * pi / nlat
        lon = j * 2pi / nlon
        r = radius(lon, lat)
        push!(verts, V3(r * cos(lon) * cos(lat), r * sin(lon) * cos(lat), r * sin(lat)))
    end
    push!(verts, V3(0.0, 0.0, -radius(0.0, -pi / 2)))
    plates = P3[]
    idx(i, j) = 1 + (i - 1) * nlon + mod(j, nlon) + 1
    south = length(verts)
    for j in 0:nlon-1
        push!(plates, P3(1, idx(1, j), idx(1, j + 1)))
    end
    for i in 1:nlat-2, j in 0:nlon-1
        push!(plates, P3(idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)))
        push!(plates, P3(idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)))
    end
    for j in 0:nlon-1
        push!(plates, P3(south, idx(nlat - 1, j + 1), idx(nlat - 1, j)))
    end
    return verts, plates
end

"Tessellated ellipsoid with the semi-axes `a`, `b`, `c`."
ellipsoid_shape(a, b, c, nlon, nlat) = latlon_shape(nlon, nlat) do lon, lat
    u = V3(cos(lon) * cos(lat), sin(lon) * cos(lat), sin(lat))
    1 / sqrt((u[1] / a)^2 + (u[2] / b)^2 + (u[3] / c)^2)
end

"Rough closed surface: each grid radius is `r0` times a random factor in [1 - amp, 1 + amp]."
function rough_shape(rng::AbstractRNG, r0, amp, nlon, nlat)
    tab = Dict{Tuple{Float64,Float64},Float64}()
    latlon_shape(nlon, nlat) do lon, lat
        get!(() -> r0 * (1 + amp * (2rand(rng) - 1)), tab, (lon, lat))
    end
end

"Cube with the half edge `h`, 12 plates."
function cube_shape(h)
    verts = [V3(x, y, z) for z in (-h, h) for y in (-h, h) for x in (-h, h)]
    plates = [P3(1, 3, 2), P3(2, 3, 4), P3(5, 6, 7), P3(6, 8, 7), P3(1, 2, 5), P3(2, 6, 5),
        P3(3, 7, 4), P3(4, 7, 8), P3(1, 5, 3), P3(3, 5, 7), P3(2, 4, 6), P3(4, 8, 6)]
    return verts, plates
end

"""
The plates of `(verts, plates)` whose centroid has a longitude in `[lon1, lon2]` and a
latitude in `[lat1, lat2]`, with the vertices that they use, renumbered.
"""
function patch_shape(verts::Vector{V3}, plates::Vector{P3}, lon1, lon2, lat1, lat2)
    keep = P3[]
    for p in plates
        c = (verts[p[1]] + verts[p[2]] + verts[p[3]]) / 3
        lon = atan(c[2], c[1])
        lat = atan(c[3], hypot(c[1], c[2]))
        (lon1 <= lon <= lon2 && lat1 <= lat <= lat2) && push!(keep, p)
    end
    used = sort!(unique!([k for p in keep for k in p]))
    renum = Dict(k => i for (i, k) in enumerate(used))
    return verts[used], [P3(renum[p[1]], renum[p[2]], renum[p[3]]) for p in keep]
end

"One type 2 segment to write: the arguments of `dskw02` without the bounds of the third coordinate."
struct DskSpec
    center::Int
    surfid::Int
    dclass::Int
    frame::String
    corsys::Int
    corpar::Vector{Float64}
    bounds::NTuple{4,Float64}
    first::Float64
    last::Float64
    finscl::Float64
    corscl::Int
    verts::Vector{V3}
    plates::Vector{P3}
end
corpar10(s::DskSpec) = (c = zeros(10); c[1:length(s.corpar)] = s.corpar; c)
index_sizes(s::DskSpec) = (voxpsz=300_000, worksz=max(60 * length(s.plates) + 10_000, 300_000),
    voxlsz=60 * length(s.plates) + 10_000)
index_spxisz(s::DskSpec) = 100_007 + index_sizes(s).voxpsz + index_sizes(s).voxlsz +
                           2 * length(s.verts) + 3 * length(s.plates) + 100

# ---------------------------------------------------------------------------------------
# Direct calls into CSPICE
# ---------------------------------------------------------------------------------------

function c_dskmi2(s::DskSpec; makvtl::Bool=true, finscl=s.finscl, corscl=s.corscl,
        sizes=index_sizes(s), spxisz=index_spxisz(s))
    work = zeros(Cint, 2, max(sizes.worksz, 1))
    spaixd = zeros(10)
    spaixi = zeros(Cint, max(spxisz, 1))
    ccall((:dskmi2_c, DSKLIB), Cvoid, (Cint, Ptr{Float64}, Cint, Ptr{Cint}, Float64, Cint, Cint,
            Cint, Cint, Cint, Cint, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}), length(s.verts),
        flat(s.verts), length(s.plates), flat(s.plates), finscl, corscl, sizes.worksz,
        sizes.voxpsz, sizes.voxlsz, makvtl, spxisz, work, spaixd, spaixi)
    c_check()
    return spaixd, spaixi
end

function c_dskrb2(verts, plates, corsys, corpar)
    mn = Ref(0.0)
    mx = Ref(0.0)
    ccall((:dskrb2_c, DSKLIB), Cvoid, (Cint, Ptr{Float64}, Cint, Ptr{Cint}, Cint, Ptr{Float64},
            Ref{Float64}, Ref{Float64}), length(verts), flat(verts), length(plates),
        flat(plates), corsys, corpar, mn, mx)
    c_check()
    return mn[], mx[]
end

function c_dskopn(path, ifname, ncomch)
    h = Ref{Cint}(0)
    ccall((:dskopn_c, DSKLIB), Cvoid, (Cstring, Cstring, Cint, Ref{Cint}), path, ifname, ncomch, h)
    c_check()
    return h[]
end

function c_dskw02(h, s::DskSpec, mn3, mx3, spaixd, spaixi; bounds=s.bounds, first=s.first,
        last=s.last, dclass=s.dclass, corsys=s.corsys, frame=s.frame, corpar=corpar10(s),
        plates=s.plates)
    ccall((:dskw02_c, DSKLIB), Cvoid, (Cint, Cint, Cint, Cint, Cstring, Cint, Ptr{Float64},
            Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Cint,
            Ptr{Float64}, Cint, Ptr{Cint}, Ptr{Float64}, Ptr{Cint}), h, s.center, s.surfid,
        dclass, frame, corsys, corpar, bounds[1], bounds[2], bounds[3], bounds[4], mn3, mx3,
        first, last, length(s.verts), flat(s.verts), length(plates), flat(plates), spaixd,
        spaixi)
    c_check()
end

function c_dskcls(h, optmiz::Bool=true)
    ccall((:dskcls_c, DSKLIB), Cvoid, (Cint, Cint), h, optmiz)
    c_check()
end

"Write the segments `specs` to the new DSK file `path` with the CSPICE writer."
function c_write_dsk(path, specs; ifname="Allspice test DSK", ncomch=0, makvtl=true,
        comments=String[])
    isfile(path) && rm(path)
    h = c_dskopn(path, ifname, ncomch)
    for s in specs
        spaixd, spaixi = c_dskmi2(s; makvtl)
        mn3, mx3 = c_dskrb2(s.verts, s.plates, s.corsys, corpar10(s))
        c_dskw02(h, s, mn3, mx3, spaixd, spaixi)
    end
    isempty(comments) || (SPICE.dasac(h, comments))
    c_dskcls(h)
    return path
end

"Write the segments `specs` to the new DSK file `path` with the Allspice writer."
function neo_write_dsk(path, specs; ifname="Allspice test DSK", ncomch=0, makvtl=true,
        comments=String[], ctx=default_context())
    isfile(path) && rm(path)
    h = dskopn(path, ifname, ncomch)
    for s in specs
        z = index_sizes(s)
        spaixd, spaixi = dskmi2(s.verts, s.plates, s.finscl, s.corscl, z.worksz, z.voxpsz,
            z.voxlsz, makvtl, index_spxisz(s))
        mn3, mx3 = dskrb2(s.verts, s.plates, s.corsys, s.corpar)
        dskw02(ctx, h, s.center, s.surfid, s.dclass, s.frame, s.corsys, s.corpar, s.bounds...,
            mn3, mx3, s.first, s.last, s.verts, s.plates, spaixd, spaixi)
    end
    isempty(comments) || dasac(h, comments)
    dskcls(h, true)
    return path
end

c_dasopr(path) = (h = Ref{Cint}(0); ccall((:dasopr_c, DSKLIB), Cvoid, (Cstring, Ref{Cint}), path, h); c_check(); h[])
c_dasopw(path) = (h = Ref{Cint}(0); ccall((:dasopw_c, DSKLIB), Cvoid, (Cstring, Ref{Cint}), path, h); c_check(); h[])
c_dascls(h) = (ccall((:dascls_c, DSKLIB), Cvoid, (Cint,), h); c_check())

function c_daslla(h)
    c = Ref{Cint}(0); d = Ref{Cint}(0); i = Ref{Cint}(0)
    ccall((:daslla_, DSKLIB), Cint, (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}), h, c, d, i)
    c_check()
    return Int(c[]), Int(d[]), Int(i[])
end

function c_dashfs(h)
    r = [Ref{Cint}(0) for _ in 1:5]
    lastla = zeros(Cint, 3); lastrc = zeros(Cint, 3); lastwd = zeros(Cint, 3)
    ccall((:dashfs_, DSKLIB), Cint, (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint},
            Ref{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), h, r[1], r[2], r[3], r[4], r[5],
        lastla, lastrc, lastwd)
    c_check()
    return Int(r[1][]), Int(r[2][]), Int(r[3][]), Int(r[4][]), Int(r[5][]),
        Tuple(Int.(lastla)), Tuple(Int.(lastrc)), Tuple(Int.(lastwd))
end

function c_dasrdi(h, first, last)
    out = zeros(Cint, max(last - first + 1, 0))
    ccall((:dasrdi_c, DSKLIB), Cvoid, (Cint, Cint, Cint, Ptr{Cint}), h, first, last, out)
    c_check()
    return out
end

function c_dasrdd(h, first, last)
    out = zeros(Float64, max(last - first + 1, 0))
    ccall((:dasrdd_c, DSKLIB), Cvoid, (Cint, Cint, Cint, Ptr{Float64}), h, first, last, out)
    c_check()
    return out
end

function c_dasrdc(h, first, last)
    n = last - first + 1
    out = fill(UInt8(' '), n)
    ccall((:dasrdc_, DSKLIB), Cint, (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint},
            Ptr{UInt8}, Cint), h, first, last, 1, n, out, n)
    c_check()
    return String(out)
end

"All integers and all doubles of the DAS file `path`, read with CSPICE."
function c_das_data(path)
    h = c_dasopr(path)
    lastc, lastd, lasti = c_daslla(h)
    ints = lasti > 0 ? c_dasrdi(h, 1, lasti) : Cint[]
    dps = lastd > 0 ? c_dasrdd(h, 1, lastd) : Float64[]
    chars = lastc > 0 ? c_dasrdc(h, 1, lastc) : ""
    c_dascls(h)
    return ints, dps, chars
end

"All integers, doubles and characters of the DAS file `path`, read with Allspice."
function neo_das_data(path)
    h = dasopr(path)
    lastc, lastd, lasti = daslla(h)
    out = (Cint.(dasrdi(h, 1, lasti)), dasrdd(h, 1, lastd), dasrdc(h, 1, lastc))
    dascls(h)
    return out
end

bits_equal(a::AbstractVector{Float64}, b::AbstractVector{Float64}) =
    length(a) == length(b) && all(reinterpret(UInt64, a) .== reinterpret(UInt64, b))

function c_dskx02(h, dla, vertex, raydir)
    plid = Ref{Cint}(0); xpt = zeros(3); found = Ref{Cint}(0)
    ccall((:dskx02_c, DSKLIB), Cvoid, (Cint, Ref{SPICE.DLADescr}, Ptr{Float64}, Ptr{Float64},
            Ref{Cint}, Ptr{Float64}, Ref{Cint}), h, dla, collect(Float64, vertex),
        collect(Float64, raydir), plid, xpt, found)
    c_check()
    return found[] != 0, Int(plid[]), V3(xpt)
end

function c_dskxsi(target, srflst, et, fixref, vertex, raydir; pri=false)
    xpt = zeros(3); h = Ref{Cint}(0); dla = Ref{SPICE.DLADescr}(); dsk = Ref{SPICE.DSKDescr}()
    dc = zeros(1); ic = zeros(Cint, 1); found = Ref{Cint}(0)
    ccall((:dskxsi_c, DSKLIB), Cvoid, (Cint, Cstring, Cint, Ptr{Cint}, Float64, Cstring,
            Ptr{Float64}, Ptr{Float64}, Cint, Cint, Ptr{Float64}, Ref{Cint},
            Ref{SPICE.DLADescr}, Ref{SPICE.DSKDescr}, Ptr{Float64}, Ptr{Cint}, Ref{Cint}),
        pri, target, length(srflst), collect(Cint, srflst), et, fixref, collect(Float64, vertex),
        collect(Float64, raydir), 1, 1, xpt, h, dla, dsk, dc, ic, found)
    c_check()
    return found[] != 0, V3(xpt), Int(ic[1]), Int(h[]), dla[], dsk[]
end

function c_dskxv(target, srflst, et, fixref, vtx::Vector{V3}, dir::Vector{V3}; pri=false)
    n = length(vtx)
    xpt = zeros(3, max(n, 1)); fnd = zeros(Cint, max(n, 1))
    ccall((:dskxv_c, DSKLIB), Cvoid, (Cint, Cstring, Cint, Ptr{Cint}, Float64, Cstring, Cint,
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Cint}), pri, target, length(srflst),
        collect(Cint, srflst), et, fixref, n, flat(vtx), flat(dir), xpt, fnd)
    c_check()
    return [V3(xpt[:, i]) for i in 1:n], fnd[1:n] .!= 0
end

function c_dski02(h, dla, item, start, room)
    n = Ref{Cint}(0); out = zeros(Cint, max(room, 1))
    ccall((:dski02_c, DSKLIB), Cvoid, (Cint, Ref{SPICE.DLADescr}, Cint, Cint, Cint, Ref{Cint},
            Ptr{Cint}), h, dla, item, start - 1, room, n, out)
    c_check()
    return Int.(out[1:n[]])
end

function c_dskd02(h, dla, item, start, room)
    n = Ref{Cint}(0); out = zeros(Float64, max(room, 1))
    ccall((:dskd02_c, DSKLIB), Cvoid, (Cint, Ref{SPICE.DLADescr}, Cint, Cint, Cint, Ref{Cint},
            Ptr{Float64}), h, dla, item, start - 1, room, n, out)
    c_check()
    return out[1:n[]]
end

function c_dskb02(h, dla)
    i = [Ref{Cint}(0) for _ in 1:7]
    vtxbds = zeros(2, 3); voxsiz = Ref(0.0); voxori = zeros(3); vgrext = zeros(Cint, 3)
    ccall((:dskb02_c, DSKLIB), Cvoid, (Cint, Ref{SPICE.DLADescr}, Ref{Cint}, Ref{Cint}, Ref{Cint},
            Ptr{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint},
            Ref{Cint}), h, dla, i[1], i[2], i[3], vtxbds, voxsiz, voxori, vgrext, i[4], i[5],
        i[6], i[7])
    c_check()
    # dskb02_c fills the bounds as [3][2] in C order: (min, max) for x, then y, then z.
    return Int(i[1][]), Int(i[2][]), Int(i[3][]), vtxbds, voxsiz[], voxori, Int.(vgrext),
        Int(i[4][]), Int(i[5][]), Int(i[6][]), Int(i[7][])
end

c_dskgtl(k) = (v = Ref(0.0); ccall((:dskgtl_c, DSKLIB), Cvoid, (Cint, Ref{Float64}), k, v); c_check(); v[])
c_dskstl(k, v) = (ccall((:dskstl_c, DSKLIB), Cvoid, (Cint, Float64), k, v); c_check())

"The fields of a DLA descriptor of either engine as a tuple."
dla_tuple(d) = (Int(d.bwdptr), Int(d.fwdptr), Int(d.ibase), Int(d.isize), Int(d.dbase),
    Int(d.dsize), Int(d.cbase), Int(d.csize))
"The fields of a DSK descriptor of either engine as a tuple."
dsk_tuple(d) = (Int(d.surfce), Int(d.center), Int(d.dclass), Int(d.dtype), Int(d.frmcde),
    Int(d.corsys), Tuple(d.corpar), d.co1min, d.co1max, d.co2min, d.co2max, d.co3min,
    d.co3max, d.start, d.stop)

# SPICE.jl passes `start - 1` to `dskv02_c` and `dskp02_c`, which take a number that starts
# at 1, and it does not check for errors there; the calls are made directly.
function c_dskv02(h, dla, start, room)
    n = Ref{Cint}(0); out = zeros(3, max(room, 1))
    ccall((:dskv02_c, DSKLIB), Cvoid, (Cint, Ref{SPICE.DLADescr}, Cint, Cint, Ref{Cint},
            Ptr{Float64}), h, dla, start, room, n, out)
    c_check()
    return [V3(out[:, i]) for i in 1:n[]]
end

function c_dskp02(h, dla, start, room)
    n = Ref{Cint}(0); out = zeros(Cint, 3, max(room, 1))
    ccall((:dskp02_c, DSKLIB), Cvoid, (Cint, Ref{SPICE.DLADescr}, Cint, Cint, Ref{Cint},
            Ptr{Cint}), h, dla, start, room, n, out)
    c_check()
    return [P3(out[:, i]) for i in 1:n[]]
end

function c_dskn02(h, dla, plid)
    out = zeros(3)
    ccall((:dskn02_c, DSKLIB), Cvoid, (Cint, Ref{SPICE.DLADescr}, Cint, Ptr{Float64}), h, dla, plid, out)
    c_check()
    return V3(out)
end
