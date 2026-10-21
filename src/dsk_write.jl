# The DSK type 2 writer: `dskopn`, the spatial index (`dskmi2`), the bounds of the third
# coordinate (`dskrb2`), the segment writer (`dskw02`) and `dskcls`.
#
# The spatial index has the layout of CSPICE. `spaixd` holds the vertex bounds (1:6),
# the voxel grid origin (7:9) and the voxel size (10). `spaixi` holds the voxel grid
# extents (1:3), the coarse voxel scale (4), the sizes of the voxel-plate pointer array
# (5), of the voxel-plate list (6) and of the vertex-plate list (7), the coarse voxel grid
# (8:100007), and then the voxel-plate pointers, the voxel-plate list, the vertex-plate
# pointers and the vertex-plate list, one directly after the other.

const _DSK_IXIFIX = _DSK_MAXCGR + 7       # size of the fixed part of the integer index

"""
    dskopn(fname, ifname, ncomch) -> DASHandle

Create the DSK file `fname` and open it for writing, as `dskopn_c`.

`ifname` is the internal file name and `ncomch` is the number of comment characters to
make room for. The new file has no segment. Segments that [`dskw02`](@ref) adds are in
memory, and the file on disk is not a complete DSK file, until [`dskcls`](@ref) runs. An
existing file is never replaced.

Throws `SpiceError(:emptystring)` if `fname` is empty, `:badrecordcount` if `ncomch < 0`,
`:fileopenfailed` if `fname` exists or cannot be created, and the other error kinds of
[`dasonw`](@ref).

See also [`dskw02`](@ref) for an example, [`dasac`](@ref) to add comments.
"""
dskopn(fname::AbstractString, ifname::AbstractString, ncomch::Integer) =
    dlaopn(fname, "DSK", ifname, ncomch)

"""
    dskcls(handle, optmiz=true) -> Nothing

Close the DSK file of `handle`, as `dskcls_c`. A file that is open for writing is written
to disk by this call.

`optmiz` has no effect. The file is always written segregated (see
[`DASHandle`](@ref)), which is what `optmiz = true` gives in CSPICE.

Throws `SpiceError(:daswritefail)` if the file cannot be written; the handle then stays
open.
"""
dskcls(h::DASHandle, optmiz::Bool=true) = dascls(h)

# ---------------------------------------------------------------------------------------
# Linked cells (zzinilnk, zzaddlnk, zztrvlnk, zzuntngl)
# ---------------------------------------------------------------------------------------

"`zzinilnk`: all `maxp` pointers of `a[poff+1:...]` become -1 and the cell count 0."
function _dsk_inilnk!(a::Vector{Int}, poff::Int, maxp::Int, maxc::Int)
    maxp < 1 && _dsk_throw(:valueoutofrange, "the pointer array size MAXP = $maxp; it must be >= 1")
    maxc < maxp && _dsk_throw(:valueoutofrange, "Cell array size MAXC = $maxc; size must " *
        "be at least as large as pointer array size $maxp.")
    fill!(view(a, poff+1:poff+maxp), -1)
    return 0
end

"`zzaddlnk`: add the value `bval` to the list of the pointer `aval`; returns the new cell count."
@inline function _dsk_addlnk!(a::Vector{Int}, poff::Int, aval::Int, bval::Int, maxa::Int,
        cval::Vector{Int32}, cnext::Vector{Int32}, ncell::Int)
    (1 <= aval <= maxa) || _dsk_throw(:avaloutofrange,
        "the index AVAL is out of range: AVAL = $aval, MAXA = $maxa")
    ncell += 1
    ncell > length(cval) && _dsk_throw(:cellarraytoosmall,
        "NCELL is larger than the cell array size $(length(cval))")
    @inbounds cval[ncell] = bval
    @inbounds cnext[ncell] = a[poff+aval]
    @inbounds a[poff+aval] = ncell
    return ncell
end

"""
    _dsk_untngl!(a, poff, nptr, cval, cnext, loff, maxb) -> nout

`zzuntngl`: write the lists of the pointers `1:nptr` one after the other to
`a[loff+1:loff+maxb]`, each as its element count followed by its elements (the last
one that was added first), and make each pointer the position of its count, or -1 for
an empty list. Returns the number of elements that were written.
"""
function _dsk_untngl!(a::Vector{Int}, poff::Int, nptr::Int, cval::Vector{Int32},
        cnext::Vector{Int32}, loff::Int, maxb::Int)
    nptr > length(cval) && _dsk_throw(:barraytoosmall,
        "the input pointer array is larger than the cell array")
    room = maxb
    ptrdex = 0
    for aval in 1:nptr
        # invariant: a[loff+1:loff+ptrdex] holds the lists of the pointers 1 .. aval-1
        (ptrdex + 2 > maxb || room <= 0) && _dsk_throw(:barraytoosmall,
            "the output list array of size $maxb is too small")
        nb = 0
        a[loff+ptrdex+1] = 0
        a[loff+ptrdex+2] = 0
        ptr = a[poff+aval]
        while ptr != -1
            (ptr < -1 || ptr == 0 || ptr > length(cval)) && _dsk_throw(:pointeroutofrange,
                "the cell pointer $ptr is out of range")
            nb += 1
            nb > room && _dsk_throw(:barraytoosmall,
                "the output list array of size $maxb is too small")
            # CSPICE stores the count in the first element after the walk; an element
            # count of `room` then does not fit, which the check above reports.
            loff + ptrdex + 1 + nb <= length(a) || _dsk_throw(:barraytoosmall,
                "the output list array of size $maxb is too small")
            a[loff+ptrdex+1+nb] = cval[ptr]
            ptr = Int(cnext[ptr])
        end
        a[loff+ptrdex+1] = nb
        if nb > 0
            a[poff+aval] = ptrdex + 1
            room -= nb + 1
            ptrdex += 1 + nb
        else
            a[poff+aval] = -1
        end
    end
    return ptrdex
end

# ---------------------------------------------------------------------------------------
# The spatial index
# ---------------------------------------------------------------------------------------

"`zzgetvox`: the voxel (1-based) of the point `(x, y, z)`, or `nothing` if it is outside of the grid."
@inline function _dsk_getvox(voxsiz::Float64, voxori::NTuple{3,Float64}, nvox::NTuple{3,Int},
        xyz::NTuple{3,Float64})
    voxsiz <= 0.0 && _dsk_throw(:nonpositivevalue, "the voxel size must be positive")
    t1 = (xyz[1] - voxori[1]) / voxsiz
    (t1 < 0.0 || t1 > nvox[1]) && return nothing
    t2 = (xyz[2] - voxori[2]) / voxsiz
    (t2 < 0.0 || t2 > nvox[2]) && return nothing
    t3 = (xyz[3] - voxori[3]) / voxsiz
    (t3 < 0.0 || t3 > nvox[3]) && return nothing
    return (min(unsafe_trunc(Int, t1) + 1, nvox[1]), min(unsafe_trunc(Int, t2) + 1, nvox[2]),
        min(unsafe_trunc(Int, t3) + 1, nvox[3]))
end

"`d_nint`: the nearest integer, with a tie away from zero."
@inline _dsk_nint(x::Float64) = round(x, RoundNearestTiesAway)

"""
    _dsk_mkspin!(spaixd, spaixi, vf, pf, np, voxscl, cgscal, maxptr, maxvxl, cval, cnext) -> nvxtot

`zzmkspin`: the voxel grid and the voxel-plate lists of the plates `pf` (flat, 3 `np`
vertex numbers) on the vertices `vf` (flat). The results go to the places of the
spatial index: the extents to `spaixd[1:6]`, the origin to `spaixd[7:9]`, the voxel
size to `spaixd[10]`, the grid extents to `spaixi[1:3]`, the pointer count to
`spaixi[5]`, the list size to `spaixi[6]`, the coarse grid to `spaixi[8:...]`, the
pointers to `spaixi[100008:...]` and the list after `maxptr` pointers.
"""
function _dsk_mkspin!(spaixd::Vector{Float64}, spaixi::Vector{Int}, vf::Vector{Float64},
        pf::Vector{Int}, np::Int, voxscl::Float64, cgscal::Int, maxptr::Int, maxvxl::Int,
        cval::Vector{Int32}, cnext::Vector{Int32})
    (np < 1 || np > _DSK_MAXPLT) && _dsk_throw(:valueoutofrange,
        "the plate count NP = $np; the count must be in the range 1:$_DSK_MAXPLT")
    cgscal < 1 && _dsk_throw(:valueoutofrange, "the coarse voxel scale = $cgscal; it must be >= 1")
    avext = 0.0
    xmin = ymin = zmin = _DSK_DPMAX
    xmax = ymax = zmax = -_DSK_DPMAX
    for i in 1:np                        # invariant: the extents cover the plates 1 .. i-1
        k1 = 3 * pf[3i-2]
        k2 = 3 * pf[3i-1]
        k3 = 3 * pf[3i]
        bxmin = min(vf[k1-2], vf[k2-2], vf[k3-2])
        bxmax = max(vf[k1-2], vf[k2-2], vf[k3-2])
        bymin = min(vf[k1-1], vf[k2-1], vf[k3-1])
        bymax = max(vf[k1-1], vf[k2-1], vf[k3-1])
        bzmin = min(vf[k1], vf[k2], vf[k3])
        bzmax = max(vf[k1], vf[k2], vf[k3])
        xmin = min(xmin, bxmin)
        xmax = max(xmax, bxmax)
        ymin = min(ymin, bymin)
        ymax = max(ymax, bymax)
        zmin = min(zmin, bzmin)
        zmax = max(zmax, bzmax)
        avext = avext + abs(bxmax - bxmin) + abs(bymax - bymin) + abs(bzmax - bzmin)
    end
    extent = (xmin, xmax, ymin, ymax, zmin, zmax)
    for k in 1:6
        spaixd[k] = extent[k]
    end
    avext /= np * 3
    voxsiz = voxscl * avext
    mdltol = voxsiz * 0.001
    xe = (xmin - mdltol, xmax + mdltol, ymin - mdltol, ymax + mdltol, zmin - mdltol, zmax + mdltol)
    cvxsiz = voxsiz * cgscal
    xvmin = _dsk_nint(xmin / cvxsiz - 1.0)
    yvmin = _dsk_nint(ymin / cvxsiz - 1.0)
    zvmin = _dsk_nint(zmin / cvxsiz - 1.0)
    xvmax = _dsk_nint(xmax / cvxsiz + 1.0)
    yvmax = _dsk_nint(ymax / cvxsiz + 1.0)
    zvmax = _dsk_nint(zmax / cvxsiz + 1.0)
    voxori = (xvmin * cvxsiz, yvmin * cvxsiz, zvmin * cvxsiz)
    spans = (xvmax - xvmin, yvmax - yvmin, zvmax - zvmin)
    all(s -> isfinite(s) && s * cgscal <= _DSK_MAXVOX, spans) || _dsk_throw(:valueoutofrange,
        "the fine voxel grid is too large or not defined; check the fine voxel scale")
    nx = Int(_dsk_nint(spans[1])) * cgscal
    ny = Int(_dsk_nint(spans[2])) * cgscal
    nz = Int(_dsk_nint(spans[3])) * cgscal
    nvox = (nx, ny, nz)
    nvxtot = nx * ny * nz
    nvxtot > _DSK_MAXVOX && _dsk_throw(:valueoutofrange,
        "the fine voxel count NVXTOT = $nvxtot; the count must be in the range 1:$_DSK_MAXVOX")
    cgscal > Float64(nvxtot)^0.33333333333333331 && _dsk_throw(:valueoutofrange,
        "the coarse voxel scale = $cgscal; the scale must be in the range " *
        "1:NVXTOT**(1/3), where NVXTOT = $nvxtot")
    npcg = cgscal * (cgscal * cgscal)
    nvxtot % npcg == 0 || _dsk_throw(:incompatiblescale,
        "the total number of voxels, $nvxtot, is not a multiple of the coarse scale cubed")
    ncgflg = nvxtot ÷ npcg
    ncgflg > _DSK_MAXCGR && _dsk_throw(:coarsegridoverflow, "the number of coarse voxels " *
        "$ncgflg is larger than the maximum $_DSK_MAXCGR")
    spaixd[7], spaixd[8], spaixd[9] = voxori
    spaixd[10] = voxsiz
    spaixi[1], spaixi[2], spaixi[3] = nvox
    vxpoff = _DSK_IXIFIX
    ncell = _dsk_inilnk!(spaixi, vxpoff, maxptr, length(cval))
    cg1 = nx ÷ cgscal
    cg2 = ny ÷ cgscal
    fill!(view(spaixi, 8:7+ncgflg), 0)
    to = 1
    for i in 1:np
        # invariant: the plates 1 .. i-1 are in the lists of all voxels that they can touch
        k1 = 3 * pf[3i-2]
        k2 = 3 * pf[3i-1]
        k3 = 3 * pf[3i]
        bxmin = _dsk_brcktd(min(min(vf[k1-2], vf[k2-2]), vf[k3-2]) - mdltol, xe[1], xe[2])
        bxmax = _dsk_brcktd(max(max(vf[k1-2], vf[k2-2]), vf[k3-2]) + mdltol, xe[1], xe[2])
        bymin = _dsk_brcktd(min(min(vf[k1-1], vf[k2-1]), vf[k3-1]) - mdltol, xe[3], xe[4])
        bymax = _dsk_brcktd(max(max(vf[k1-1], vf[k2-1]), vf[k3-1]) + mdltol, xe[3], xe[4])
        bzmin = _dsk_brcktd(min(min(vf[k1], vf[k2]), vf[k3]) - mdltol, xe[5], xe[6])
        bzmax = _dsk_brcktd(max(max(vf[k1], vf[k2]), vf[k3]) + mdltol, xe[5], xe[6])
        gmin = _dsk_getvox(voxsiz, voxori, nvox, (bxmin, bymin, bzmin))
        gmax = _dsk_getvox(voxsiz, voxori, nvox, (bxmax, bymax, bzmax))
        (gmin === nothing || gmax === nothing) && _dsk_throw(:bug,
            "a plate of the model is outside of the voxel grid")
        for iz in gmin[3]:gmax[3], iy in gmin[2]:gmax[2], ix in gmin[1]:gmax[1]
            c1 = (ix - 1) ÷ cgscal + 1
            c2 = (iy - 1) ÷ cgscal + 1
            c3 = (iz - 1) ÷ cgscal + 1
            cgof1d = (iz - cgscal * (c3 - 1) - 1) * cgscal * cgscal +
                     (iy - cgscal * (c2 - 1) - 1) * cgscal + (ix - cgscal * (c1 - 1))
            cvid = c1 + cg1 * (c2 - 1 + (c3 - 1) * cg2)
            if spaixi[7+cvid] == 0
                spaixi[7+cvid] = to
                to += npcg
            end
            ixptr = spaixi[7+cvid] - 1 + cgof1d
            ncell = _dsk_addlnk!(spaixi, vxpoff, ixptr, i, maxptr, cval, cnext, ncell)
        end
    end
    nvxptr = to - 1
    spaixi[5] = nvxptr
    spaixi[6] = _dsk_untngl!(spaixi, vxpoff, nvxptr, cval, cnext, vxpoff + maxptr, maxvxl)
    return nvxtot
end

"""
    dskmi2(vrtces, plates, finscl, corscl, worksz, voxpsz, voxlsz, makvtl, spxisz)
        -> (spaixd, spaixi)

Return the spatial index of a type 2 segment, as `dskmi2_c`: the double part `spaixd` (10
elements) and the integer part `spaixi` (`spxisz` elements). [`dskw02`](@ref) takes both.

# Arguments
- `vrtces`, `plates`: the vertices (km) and the plates, each a 3 x n matrix or a vector of
  3-vectors. A plate is three vertex numbers, which start at 1.
- `finscl`: the edge length of a fine voxel, in units of the mean extent of the plates.
- `corscl`: the number of fine voxels along one edge of a coarse voxel.
- `worksz`: the size of the workspace. It must hold the pairs of a plate and a fine
  voxel that meet: `np` times the mean number of fine voxels that a plate meets. With
  `makvtl` it must also hold `nv + 3np`. The least value is `np + 1`.
- `voxpsz`: the size of the voxel-plate pointer array, at least `corscl^3` times the
  number of coarse voxels that are not empty.
- `voxlsz`: the size of the voxel-plate list, at least the first bound of `worksz`. The
  least value is `np + 1`.
- `makvtl`: `true` to make the vertex-plate list.
- `spxisz`: the size of `spaixi`, at least `100007 + voxpsz + voxlsz`, plus `2nv + 3np`
  if `makvtl` is `true`.

The index is equal to the one that CSPICE makes.

Throws `SpiceError`:
- `:badfinevoxelscale` unless `finscl > 0`, `:badcoarsevoxscale` unless `corscl >= 1`;
- `:badvertexcount` or `:badplatecount` if the number of vertices or of plates is out of
  range;
- `:workspacetoosmall` if `worksz < np + 1`, `:ptrarraytoosmall` if `voxpsz < 1`,
  `:platelisttoosmall` if `voxlsz < np + 1`, `:intindextoosmall` if `spxisz` is below
  its bound;
- `:cellarraytoosmall` if `worksz` cannot hold the pairs, `:barraytoosmall` if a list
  does not fit in `voxlsz` or in `worksz`;
- `:valueoutofrange`, `:incompatiblescale` or `:coarsegridoverflow` if the voxel grid that
  follows from the scales is not valid;
- `:badvertexindex` if a vertex number is not in `1:nv`, which CSPICE does not check
  here.

See also [`dskrb2`](@ref), and [`dskw02`](@ref) for an example.
"""
function dskmi2(vrtces, plates, finscl::Real, corscl::Integer, worksz::Integer,
        voxpsz::Integer, voxlsz::Integer, makvtl::Bool, spxisz::Integer)
    nv = _dsk_count3(vrtces)
    np = _dsk_count3(plates)
    finscl <= 0.0 && _dsk_throw(:badfinevoxelscale, "the fine voxel scale = $finscl; the scale must be positive")
    corscl < 1 && _dsk_throw(:badcoarsevoxscale, "the coarse voxel scale = $corscl; the scale must be >= 1")
    (nv < 3 || nv > _DSK_MAXVRT) && _dsk_throw(:badvertexcount,
        "the vertex count NV = $nv; the count must be in the range 3:$_DSK_MAXVRT")
    (np < 1 || np > _DSK_MAXPLT) && _dsk_throw(:badplatecount,
        "the plate count NP = $np; the count must be in the range 1:$_DSK_MAXPLT")
    worksz < np + 1 && _dsk_throw(:workspacetoosmall,
        "the workspace size = $worksz; the size is too small to hold all voxel-plate associations")
    voxpsz < 1 && _dsk_throw(:ptrarraytoosmall, "the voxel-plate pointer array size = $voxpsz; " *
        "the size is too small to be compatible with the coarse voxel scale $corscl")
    voxlsz < np + 1 && _dsk_throw(:platelisttoosmall,
        "the voxel-plate list array size = $voxlsz; the size is too small")
    reqsiz = _DSK_IXIFIX + Int(voxpsz) + Int(voxlsz)
    vtxlsz = makvtl ? nv + 3 * np : 0
    makvtl && (reqsiz += nv + vtxlsz)
    spxisz < reqsiz && _dsk_throw(:intindextoosmall,
        "the integer spatial index size = $spxisz; the required size is >= $reqsiz")
    vf = _dsk_flat(Float64, vrtces)
    pf = _dsk_flat(Int, plates)
    for k in eachindex(pf)               # invariant: the vertex numbers before k are valid
        (1 <= pf[k] <= nv) || _dsk_throw(:badvertexindex, "vertex index $(k - 3 * ((k - 1) ÷ 3)) " *
            "of plate $((k - 1) ÷ 3 + 1) was $(pf[k]); vertex indices must be in the range 1:$nv")
    end
    spaixd = zeros(Float64, 10)
    spaixi = zeros(Int, spxisz)
    cval = zeros(Int32, worksz)
    cnext = fill(Int32(-1), worksz)
    spaixi[4] = corscl
    vxlidx = _DSK_IXIFIX + 1 + Int(voxpsz)
    _dsk_mkspin!(spaixd, spaixi, vf, pf, np, Float64(finscl), Int(corscl), Int(voxpsz),
        Int(voxlsz), cval, cnext)
    # Move the voxel-plate list down so that it starts directly after the pointers.
    nshift = Int(voxpsz) - spaixi[5]
    for i in 1:spaixi[6]                 # invariant: the elements 1 .. i-1 of the list are moved
        j = vxlidx - 1 + i
        spaixi[j-nshift] = spaixi[j]
    end
    vxlidx -= nshift
    if makvtl
        vtpoff = vxlidx - 1 + spaixi[6]
        vtloff = vtpoff + nv
        worksz < 3 * np && _dsk_throw(:cellarraytoosmall,
            "the cell array size = $worksz; the size must be >= 3*NP, which is $(3 * np)")
        fill!(cnext, Int32(-1))
        ncell = _dsk_inilnk!(spaixi, vtpoff, nv, length(cval))
        for i in 1:np, j in 1:3
            ncell = _dsk_addlnk!(spaixi, vtpoff, pf[3*(i-1)+j], i, nv, cval, cnext, ncell)
        end
        spaixi[7] = _dsk_untngl!(spaixi, vtpoff, nv, cval, cnext, vtloff, vtxlsz)
    else
        spaixi[7] = 0
    end
    return spaixd, spaixi
end

# ---------------------------------------------------------------------------------------
# Altitude above a spheroid, for the altitude bounds of `dskrb2`
# ---------------------------------------------------------------------------------------

"`recgeo`: the altitude of `rectan` above the spheroid with the equatorial radius `re` and the flattening `f`."
function _dsk_altitude(rectan::Vec3{Float64}, re::Float64, f::Float64)
    re <= 0.0 && _dsk_throw(:valueoutofrange, "the equatorial radius was $re; it must be positive")
    f >= 1.0 && _dsk_throw(:valueoutofrange, "the flattening coefficient was $f; it must be less than 1")
    return nearpt(rectan, re, re, re - f * re)[2]
end

"""
    dskrb2(vrtces, plates, corsys, corpar) -> (mncor3, mxcor3)

Return the bounds of the third coordinate of the plates of a type 2 segment, as
`dskrb2_c`, in the unit of the vertices (km).

- `corsys = 1`, latitudinal: the radius.
- `corsys = 3`, rectangular: z.
- `corsys = 4`, planetodetic: the altitude above the spheroid with the equatorial radius
  `corpar[1]` and the flattening `corpar[2]`.

`vrtces` and `plates` are as for [`dskmi2`](@ref). For the planetodetic system the lower
bound is a safe bound, not the exact minimum, as in CSPICE.

Throws `SpiceError(:notsupported)` if `corsys` is not 1, 3 or 4, `:badvertexindex` if
a vertex number is not in `1:nv`, and for the planetodetic system `:valueoutofrange`
unless the radius is positive and the flattening is less than 1.
"""
function dskrb2(vrtces, plates, corsys::Integer, corpar)
    nv = _dsk_count3(vrtces)
    np = _dsk_count3(plates)
    for i in 1:np, k in _dsk_plate3(plates, i)
        (1 <= k <= nv) || _dsk_throw(:badvertexindex,
            "a vertex index of plate $i was $k; vertex indices must be in the range 1:$nv")
    end
    if corsys == 1
        mx = 0.0
        for i in 1:nv                    # invariant: mx is the largest norm of the vertices 1 .. i-1
            mx = max(_dsk_vnorm(_dsk_get3(vrtces, i)), mx)
        end
        mn = _DSK_DPMAX
        for i in 1:np                    # invariant: mn is the least distance of the plates 1 .. i-1
            i1, i2, i3 = _dsk_plate3(plates, i)
            _, dist = _dsk_pltnp(_DSK_ORIGIN, _dsk_get3(vrtces, i1), _dsk_get3(vrtces, i2),
                _dsk_get3(vrtces, i3))
            mn = min(dist, mn)
        end
        return mn, mx
    elseif corsys == 3
        mn = _DSK_DPMAX
        mx = -_DSK_DPMAX
        for i in 1:nv
            z = _dsk_get3(vrtces, i)[3]
            mn = min(mn, z)
            mx = max(mx, z)
        end
        return mn, mx
    elseif corsys == 4
        re = Float64(corpar[begin])
        f = Float64(corpar[begin+1])
        mx = -_DSK_DPMAX
        mn = _DSK_DPMAX
        for i in 1:nv
            mx = max(mx, _dsk_altitude(_dsk_get3(vrtces, i), re, f))
        end
        third = 0.33333333333333331
        for i in 1:np
            i1, i2, i3 = _dsk_plate3(plates, i)
            v1 = _dsk_get3(vrtces, i1)
            v2 = _dsk_get3(vrtces, i2)
            v3 = _dsk_get3(vrtces, i3)
            center = _dsk_vlcom3(third, v1, third, v2, third, v3)
            maxd = max(max(_dsk_vdist(v1, center), _dsk_vdist(v2, center)), _dsk_vdist(v3, center))
            mn = min(mn, _dsk_altitude(center, re, f) - maxd)
        end
        return mn, mx
    end
    _dsk_throw(:notsupported, "the coordinate system code $corsys is not recognized")
end

# ---------------------------------------------------------------------------------------
# The segment writer
# ---------------------------------------------------------------------------------------

"""
    dskw02([ctx,] handle, center, surfid, dclass, frame, corsys, corpar, mncor1, mxcor1,
           mncor2, mxcor2, mncor3, mxcor3, first, last, vrtces, plates, spaixd, spaixi)
        -> Nothing

Write a type 2 segment (triangular plates) to the DSK file of `handle`, as `dskw02_c`.

# Arguments
- `handle`: a [`DASHandle`](@ref) from [`dskopn`](@ref).
- `center`, `surfid`: the body ID and the surface ID.
- `dclass`: the data class, 1 for a surface that each ray from the center meets one time,
  2 for a general surface.
- `frame`: the name or the ID of the frame of the vertices. A name that is not built in is
  looked up in the kernel pool of `ctx` (default: the default context).
- `corsys`, `corpar`: the coordinate system of the bounds (1 latitudinal, 3 rectangular,
  4 planetodetic) and its parameters; for system 4 these are the equatorial radius (km)
  and the flattening. `corpar` can have fewer than 10 elements; the rest is 0.
- `mncor1`, `mxcor1`, `mncor2`, `mxcor2`, `mncor3`, `mxcor3`: the bounds of the three
  coordinates, radians for longitude and latitude and km for the others.
  [`dskrb2`](@ref) gives those of the third coordinate.
- `first`, `last`: the coverage, seconds past J2000 TDB.
- `vrtces`, `plates`: the vertices (km) and the plates, as for [`dskmi2`](@ref), which
  gives `spaixd` and `spaixi`.

The segment becomes the last one of the file. It is on disk after [`dskcls`](@ref). A
call that throws adds nothing to the file.

Throws `SpiceError`, with the kinds of CSPICE:
- `:frameidnotfound` if the frame is not known;
- `:timesoutoforder` unless `first < last`;
- `:valueoutofrange`, `:invalidlonextent`, `:boundsoutoforder`, `:zeroboundsextent` or
  `:degeneratesurface` if the bounds are out of range, out of order or of zero extent;
- `:notsupported` unless `corsys` is 1, 3 or 4 and `dclass` is 1 or 2;
- `:valueoutofrange` or `:incompatiblescale` if a count or a parameter of the spatial
  index is out of range;
- `:badvertexindex` if a vertex number is not in `1:nv`.

Throws `SpiceError(:badarg)` if `corpar` has more than 10 elements or `spaixd` or
`spaixi` is too short.

# Examples
Write a tetrahedron, then read it back from the file and as a loaded kernel:
```jldoctest
julia> vrtces = [[0.0, 0.0, 1.0], [1.0, 0.0, -0.5], [-0.5, 0.8, -0.5], [-0.5, -0.8, -0.5]];

julia> plates = [[1, 2, 3], [1, 3, 4], [1, 4, 2], [2, 4, 3]];

julia> spaixd, spaixi = dskmi2(vrtces, plates, 1.0, 1, 1000, 1000, 1000, true, 103000);

julia> mncor3, mxcor3 = dskrb2(vrtces, plates, 1, zeros(10));

julia> path = joinpath(mktempdir(), "tetrahedron.bds");

julia> h = dskopn(path, "tetrahedron", 0);

julia> dskw02(h, 499, 1, 2, "IAU_MARS", 1, zeros(10), -pi, pi, -pi / 2, pi / 2, mncor3,
              mxcor3, -1.0e9, 1.0e9, vrtces, plates, spaixd, spaixi)

julia> dskcls(h)

julia> h = dasopr(path);

julia> dladsc = dlabfs(h);

julia> dskz02(h, dladsc)
(4, 4)

julia> dskx02(h, dladsc, [0.0, 0.0, -5.0], [0.0, 0.0, 1.0])
(4, [0.0, 0.0, -0.5])

julia> dascls(h)

julia> ctx = SpiceContext();

julia> furnsh(ctx, path)

julia> xpt, segment, plid = dskxsi(ctx, false, 499, Int[], 0.0, "IAU_MARS", [0.0, 0.0, -5.0],
                                   [0.0, 0.0, 1.0]);

julia> xpt ≈ [0.0, 0.0, -0.5], plid
(true, 4)

julia> dskp02(ctx, segment, plid, 1)[1] == plates[plid]
true
```

See also [`dskopn`](@ref), [`dskmi2`](@ref), [`dskrb2`](@ref), [`dskcls`](@ref).
"""
function dskw02(ctx::SpiceContext, h::DASHandle, center::Integer, surfid::Integer,
        dclass::Integer, frame::NameOrID, corsys::Integer, corpar, mncor1::Real, mxcor1::Real,
        mncor2::Real, mxcor2::Real, mncor3::Real, mxcor3::Real, first::Real, last::Real,
        vrtces, plates, spaixd, spaixi)
    _das_check_write(h, "dskw02")
    frmcde = frame isa Integer ? Int(frame) : namfrm(kernelset(ctx), frame)
    frmcde == 0 && _dsk_throw(:frameidnotfound, "the frame name $frame could not be mapped to an ID code")
    last <= first && _dsk_throw(:timesoutoforder,
        "the segment stop time $last must be greater than the start time $first")
    mn1, mx1, mn2, mx2, mn3, mx3 = Float64(mncor1), Float64(mxcor1), Float64(mncor2),
        Float64(mxcor2), Float64(mncor3), Float64(mxcor3)
    length(corpar) <= 10 || throw(SpiceError(:badarg, "dskw02: CORPAR has more than 10 elements"))
    par = zeros(Float64, 10)
    copyto!(par, corpar)
    twopi_ = _DSK_TWOPI
    hpi = _DSK_HALFPI
    if corsys == 1 || corsys == 4
        mn1 == mx1 && _dsk_throw(:zeroboundsextent,
            "the minimum and maximum longitude were both $(rad2deg(mn1)) degrees")
        (mn1 < -twopi_ - 1e-12 || mn1 > twopi_ - 1e-12) && _dsk_throw(:valueoutofrange,
            "the minimum longitude $(rad2deg(mn1)) degrees is out of range")
        (mx1 < -twopi_ + 1e-12 || mx1 > twopi_ + 1e-12) && _dsk_throw(:valueoutofrange,
            "the maximum longitude $(rad2deg(mx1)) degrees is out of range")
        b1 = max(-twopi_, mn1)
        b2 = min(twopi_, mx1)
        (mx1 > mn1 + twopi_ + 1e-12 || mx1 < mn1 - twopi_ - 1e-12) && _dsk_throw(
            :invalidlonextent, "the longitude bounds $(rad2deg(mn1)):$(rad2deg(mx1)) " *
            "degrees are too far apart")
        if b2 > b1
            b2 = min(b2, b1 + twopi_)
        elseif b2 < b1
            b2 = max(b2, b1 - twopi_)
        end
        (b2 == b1 || b2 == b1 - twopi_) && _dsk_throw(:zeroboundsextent,
            "the longitude bounds have zero extent after adjustment")
        (mn2 < -hpi - 1e-12 || mn2 > hpi - 1e-12) && _dsk_throw(:valueoutofrange,
            "the minimum latitude $(rad2deg(mn2)) degrees is out of range")
        (mx2 < -hpi + 1e-12 || mx2 > hpi + 1e-12) && _dsk_throw(:valueoutofrange,
            "the maximum latitude $(rad2deg(mx2)) degrees is out of range")
        b3 = max(-hpi, mn2)
        b4 = min(hpi, mx2)
        mx2 < mn2 && _dsk_throw(:boundsoutoforder, "the latitude bounds are out of order")
        if corsys == 1
            mn3 < 0.0 && _dsk_throw(:valueoutofrange, "the radius lower bound must be non-negative but was $mn3")
            mx3 <= 0.0 && _dsk_throw(:valueoutofrange, "the radius upper bound must be strictly positive but was $mx3")
        else
            par[1] <= 0.0 && _dsk_throw(:valueoutofrange, "the equatorial radius was $(par[1]); it must be positive")
            par[2] >= 1.0 && _dsk_throw(:valueoutofrange, "the flattening coefficient was $(par[2]); it must be less than 1")
            a = par[1]
            b = a * (1 - par[2])
            altlim = max(-(a * a) / b, -(b * b) / a)
            mn3 <= altlim && _dsk_throw(:degeneratesurface,
                "the minimum altitude $mn3 is not above the limit $altlim")
        end
        mx3 < mn3 && _dsk_throw(:boundsoutoforder, "the bounds of the third coordinate are out of order")
        mx3 == mn3 && _dsk_throw(:zeroboundsextent, "the bounds of the third coordinate are equal")
    elseif corsys == 3
        (mx1 <= mn1 || mx2 <= mn2 || mx3 <= mn3) && _dsk_throw(:boundsoutoforder,
            "the rectangular coordinate bounds must be strictly increasing")
        b1, b2, b3, b4 = mn1, mx1, mn2, mx2
    else
        _dsk_throw(:notsupported, "the coordinate system code $corsys is not recognized")
    end
    (dclass < 1 || dclass > 2) && _dsk_throw(:notsupported, "the data class code $dclass is not recognized")
    nv = _dsk_count3(vrtces)
    np = _dsk_count3(plates)
    (nv < 1 || nv > _DSK_MAXVRT) && _dsk_throw(:valueoutofrange,
        "the vertex count NV = $nv; the count must be in the range 1:$_DSK_MAXVRT")
    (np < 1 || np > _DSK_MAXPLT) && _dsk_throw(:valueoutofrange,
        "the plate count NP = $np; the count must be in the range 1:$_DSK_MAXPLT")
    pf = _dsk_flat(Int32, plates)
    for k in eachindex(pf)               # invariant: the vertex numbers before k are valid
        (1 <= pf[k] <= nv) || _dsk_throw(:badvertexindex, "vertex index $(k - 3 * ((k - 1) ÷ 3)) " *
            "of plate $((k - 1) ÷ 3 + 1) was $(pf[k]); vertex indices must be in the range 1:$nv")
    end
    (length(spaixd) >= 10 && length(spaixi) >= _DSK_IXIFIX) || throw(SpiceError(:badarg, 
        "dskw02: SPAIXD needs 10 elements and SPAIXI at least $_DSK_IXIFIX"))
    ix0 = firstindex(spaixi) - 1
    vgrext = (Int(spaixi[ix0+1]), Int(spaixi[ix0+2]), Int(spaixi[ix0+3]))
    cgrscl = Int(spaixi[ix0+4])
    voxnpt = Int(spaixi[ix0+5])
    voxnpl = Int(spaixi[ix0+6])
    vtxnpl = Int(spaixi[ix0+7])
    for i in 1:3
        (vgrext[i] < 1 || vgrext[i] > _DSK_MAXVOX) && _dsk_throw(:valueoutofrange,
            "the voxel grid extents are $vgrext; all must be in the range 1:$_DSK_MAXVOX")
    end
    nvxtot = vgrext[1] * vgrext[2] * vgrext[3]
    nvxtot > _DSK_MAXVOX && _dsk_throw(:valueoutofrange,
        "the fine voxel count NVXTOT = $nvxtot; the count must be in the range 1:$_DSK_MAXVOX")
    (cgrscl < 1 || cgrscl > Float64(nvxtot)^0.33333333333333331) && _dsk_throw(:valueoutofrange,
        "the coarse voxel scale = $cgrscl; the scale must be in the range 1:NVXTOT**(1/3)")
    npcg = cgrscl * (cgrscl * cgrscl)
    nvxtot % npcg == 0 || _dsk_throw(:incompatiblescale,
        "the total number of voxels, $nvxtot, is not a multiple of the coarse scale cubed")
    ncgr = nvxtot ÷ npcg
    (ncgr < 1 || ncgr > _DSK_MAXCGR) && _dsk_throw(:valueoutofrange,
        "the coarse voxel count = $ncgr; the count must be in the range 1:$_DSK_MAXCGR")
    (voxnpt >= 0 && voxnpl >= 0 && vtxnpl >= 0 &&
     length(spaixi) >= _DSK_IXIFIX + voxnpt + voxnpl + nv + vtxnpl) || throw(SpiceError(:badarg, 
        "dskw02: SPAIXI is shorter than the sizes that it holds"))
    vf = _dsk_flat(Float64, vrtces)
    pvoxpt = ix0 + _DSK_IXIFIX
    pvoxpl = pvoxpt + voxnpt
    pvtxpt = pvoxpl + voxnpl
    pvtxpl = pvtxpt + nv

    dlabns(h)
    descr = zeros(Float64, 24)
    descr[1] = surfid
    descr[2] = center
    descr[3] = dclass
    descr[4] = 2.0
    descr[5] = frmcde
    descr[6] = corsys
    descr[7:16] = par
    descr[17], descr[18], descr[19], descr[20] = b1, b2, b3, b4
    descr[21] = mn3
    descr[22] = mx3
    descr[23] = first
    descr[24] = last
    dasadd(h, descr)
    dasadd(h, Float64[spaixd[begin+k] for k in 0:9])
    dasadd(h, vf)
    dasadi(h, Int32[nv, np, nvxtot, vgrext[1], vgrext[2], vgrext[3], cgrscl, voxnpt, voxnpl, vtxnpl])
    dasadi(h, pf)
    dasadi(h, Int32[spaixi[pvoxpt+k] for k in 1:voxnpt])
    dasadi(h, Int32[spaixi[pvoxpl+k] for k in 1:voxnpl])
    dasadi(h, Int32[spaixi[pvtxpt+k] for k in 1:nv])
    dasadi(h, Int32[spaixi[pvtxpl+k] for k in 1:vtxnpl])
    dasadi(h, Int32[spaixi[ix0+7+k] for k in 1:ncgr])
    dlaens(h)
    return nothing
end
dskw02(h::DASHandle, center::Integer, surfid::Integer, dclass::Integer, frame::NameOrID,
    corsys::Integer, corpar, mncor1::Real, mxcor1::Real, mncor2::Real, mxcor2::Real,
    mncor3::Real, mxcor3::Real, first::Real, last::Real, vrtces, plates, spaixd, spaixi) =
    dskw02(default_context(), h, center, surfid, dclass, frame, corsys, corpar, mncor1, mxcor1,
        mncor2, mxcor2, mncor3, mxcor3, first, last, vrtces, plates, spaixd, spaixi)
