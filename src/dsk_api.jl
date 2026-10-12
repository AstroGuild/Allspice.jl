# The public DSK calls: tolerances, access to the items of a type 2 segment, the
# ray-surface intercept, and the plate utilities.
#
# A segment is named in one of two ways. The file form takes a `DASHandle` and a
# `DLADescr`, as `(handle, dladsc)` of CSPICE. The loaded form takes a `DSKSegment`, which
# `dskxsi` returns in place of the triple (handle, DLA descriptor, DSK descriptor), with
# a `KernelSet`, a `SpiceContext`, or nothing (the default context).

# ---------------------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------------------

"""
    dsktolerances(ctx::SpiceContext) -> DSKTolerances

Contract: returns the current DSK tolerances of `ctx` with one atomic read; the result
never changes.
"""
dsktolerances(ctx::SpiceContext) = @atomic ctx.dsktols

"""
    dskgtl([ctx], keywrd) -> Float64

The DSK tolerance or margin with the keyword `keywrd`, as `dskgtl_c`: 1 plate expansion
fraction (`XFRACT`), 2 segment boundary margin (`SGREED`), 3 segment pad margin
(`SGPADM`), 4 surface point membership margin (`PTMEMM`), 5 angular rounding margin
(`ANGMRG`), 6 longitude alias margin (`LONALI`). The values belong to the context.
Throws `SpiceError(:indexoutofrange)` for another keyword.
"""
function dskgtl(ctx::SpiceContext, keywrd::Integer)
    (1 <= keywrd <= 6) || throw(SpiceError(:indexoutofrange,
        "dskgtl: the keyword $keywrd is not in the range 1:6"))
    return getfield(dsktolerances(ctx), Int(keywrd))::Float64
end
dskgtl(keywrd::Integer) = dskgtl(default_context(), keywrd)

"""
    dskstl([ctx], keywrd, dpval) -> Nothing

Set the DSK tolerance or margin with the keyword `keywrd` (see [`dskgtl`](@ref)) to
`dpval`, as `dskstl_c`. Precondition: `keywrd` in 1:4 (`SpiceError(:indexoutofrange)`
outside 1:6, `SpiceError(:immutablevalue)` for 5 and 6). As CSPICE, the value is not
checked. Postcondition: a new immutable snapshot is stored under the lock of the
context; a query that runs on another thread uses the old or the new set, never a
mixture.
"""
function dskstl(ctx::SpiceContext, keywrd::Integer, dpval::Real)
    (1 <= keywrd <= 6) || throw(SpiceError(:indexoutofrange,
        "dskstl: the keyword $keywrd is not in the range 1:6"))
    keywrd > 4 && throw(SpiceError(:immutablevalue,
        "dskstl: the parameter with the keyword $keywrd cannot be changed"))
    lock(ctx.lock) do
        t = dsktolerances(ctx)
        v = Float64(dpval)
        @atomic ctx.dsktols = DSKTolerances(keywrd == 1 ? v : t.xfract,
            keywrd == 2 ? v : t.sgreed, keywrd == 3 ? v : t.sgpadm,
            keywrd == 4 ? v : t.ptmemm, t.angmrg, t.lonali)
    end
    return nothing
end
dskstl(keywrd::Integer, dpval::Real) = dskstl(default_context(), keywrd, dpval)

# ---------------------------------------------------------------------------------------
# Input conversion
# ---------------------------------------------------------------------------------------

"A 3-vector as `Vec3{Float64}`; throws `ArgumentError` unless `v` has three elements."
@inline _dsk_vec3(v::Vec3{Float64}) = v
function _dsk_vec3(v)
    length(v) == 3 || throw(ArgumentError("a vector of 3 elements is required, got $(length(v))"))
    return Vec3{Float64}(v[begin], v[begin+1], v[begin+2])
end

"Number of 3-vectors in `a`: the columns of a 3 x n matrix or the elements of a vector of 3-vectors."
function _dsk_count3(a::AbstractMatrix)
    size(a, 1) == 3 || throw(ArgumentError("a 3 x n matrix is required, got $(size(a))"))
    return size(a, 2)
end
_dsk_count3(a::AbstractVector) = length(a)

"Element `i` of a list of 3-vectors (see `_dsk_count3`) as `Vec3{Float64}`."
@inline _dsk_get3(a::AbstractMatrix, i::Int) = Vec3{Float64}(a[1, i], a[2, i], a[3, i])
@inline _dsk_get3(a::AbstractVector, i::Int) = _dsk_vec3(a[begin+i-1])
@inline _dsk_get3(a::AbstractVector{Vec3{Float64}}, i::Int) = a[begin+i-1]

"The list of 3-vectors `a` as one flat `Vector{T}` of length 3n, in the order of the file."
function _dsk_flat(::Type{T}, a) where {T}
    n = _dsk_count3(a)
    out = Vector{T}(undef, 3n)
    for i in 1:n                         # invariant: out[1:3(i-1)] holds the elements 1 .. i-1
        if a isa AbstractMatrix
            out[3i-2] = a[1, i]
            out[3i-1] = a[2, i]
            out[3i] = a[3, i]
        else
            e = a[begin+i-1]
            length(e) == 3 || throw(ArgumentError("element $i does not have 3 components"))
            out[3i-2] = e[begin]
            out[3i-1] = e[begin+1]
            out[3i] = e[begin+2]
        end
    end
    return out
end

# ---------------------------------------------------------------------------------------
# Segments: the file form and the loaded form
# ---------------------------------------------------------------------------------------

"The arrays and the `DSKSegment` of the segment `d` of the open DSK file `h`."
function _dsk_open_segment(h::DASHandle, d::DLADescr, who::AbstractString)
    _das_check_open(h, who)
    return h.ints, h.dps, _dsk_segment(h.ints, h.dps, d, 0, h.path)
end

"""
    _dsk_file(ks, seg) -> DSKFile

The loaded file of `seg`. Throws `SpiceError(:nosuchhandle)` if `seg` is not a segment
of the set `ks` (for example when its file was unloaded).
"""
function _dsk_file(ks::KernelSet, seg::DSKSegment)
    if 1 <= seg.file <= length(ks.dskfiles)
        file = ks.dskfiles[seg.file]
        for s in file.segments           # invariant: no segment before `s` is `seg`
            s.dladsc == seg.dladsc && s.dskdsc == seg.dskdsc && return file
        end
    end
    throw(SpiceError(:nosuchhandle, "the DSK segment does not belong to a file of this " *
        "kernel set; its file was unloaded or the set is another one"))
end

"Throws `SpiceError(:typenotsupported)` unless the segment has the data type 2."
@inline function _dsk_need02(seg::DSKSegment, who::AbstractString)
    seg.dskdsc.dtype == 2 || throw(SpiceError(:typenotsupported,
        "$who: the DSK segment has the data type $(seg.dskdsc.dtype), not 2"))
    return nothing
end

"""
    dskgd(handle, dladsc) -> DSKDescr
    dskgd(segment::DSKSegment) -> DSKDescr

The DSK descriptor of a segment, as `dskgd_c`. Throws `SpiceError(:invalidformat)` if
the segment has fewer than 24 doubles.
"""
function dskgd(h::DASHandle, d::DLADescr)
    _das_check_open(h, "dskgd")
    return _dsk_descr(h.dps, d, h.path)
end
dskgd(seg::DSKSegment) = seg.dskdsc

"`dskb02` on a header: the items in the order of CSPICE."
_dsk_b02(hdr::DSK02Header) = (hdr.nv, hdr.np, hdr.nvxtot,
    [hdr.vtxbds[1] hdr.vtxbds[3] hdr.vtxbds[5]; hdr.vtxbds[2] hdr.vtxbds[4] hdr.vtxbds[6]],
    hdr.voxsiz, hdr.voxori, collect(hdr.vgrext), hdr.cgscal, hdr.vtxnpl, hdr.voxnpt, hdr.voxnpl)

"""
    dskb02(handle, dladsc) -> (nv, np, nvxtot, vtxbds, voxsiz, voxori, vgrext, cgscal, vtxnpl, voxnpt, voxnpl)
    dskb02(segment::DSKSegment) -> the same

The bookkeeping items of a type 2 segment, as `dskb02_c`: vertex count, plate count,
voxel count, vertex bounds (2 x 3: row 1 the lower bounds, row 2 the upper bounds of x,
y, z), voxel size, voxel grid origin, voxel grid extents, coarse voxel scale, and the
sizes of the vertex-plate list, of the voxel-plate pointer array and of the voxel-plate
list.
"""
function dskb02(h::DASHandle, d::DLADescr)
    _, _, seg = _dsk_open_segment(h, d, "dskb02")
    _dsk_need02(seg, "dskb02")
    return _dsk_b02(seg.hdr)
end
function dskb02(seg::DSKSegment)
    _dsk_need02(seg, "dskb02")
    return _dsk_b02(seg.hdr)
end

"""
    dskz02(handle, dladsc) -> (nv, np)
    dskz02(segment::DSKSegment) -> (nv, np)

Vertex count and plate count of a type 2 segment, as `dskz02_c`.
"""
function dskz02(h::DASHandle, d::DLADescr)
    _, _, seg = _dsk_open_segment(h, d, "dskz02")
    _dsk_need02(seg, "dskz02")
    return seg.hdr.nv, seg.hdr.np
end
function dskz02(seg::DSKSegment)
    _dsk_need02(seg, "dskz02")
    return seg.hdr.nv, seg.hdr.np
end

"Checks of `start` and `room` that the fetch calls share; returns the element count to fetch."
function _dsk_fetch_count(who::AbstractString, what::AbstractString, start::Integer,
        room::Integer, size::Integer)
    room <= 0 && throw(SpiceError(:valueoutofrange, "$who: ROOM was $room; it must be positive"))
    (1 <= start <= size) || throw(SpiceError(:indexoutofrange,
        "$who: START was $start; the valid range of the $what is 1:$size"))
    return min(Int(room), Int(size) - Int(start) + 1)
end

function _dsk_p02(ints::Vector{Int32}, seg::DSKSegment, start::Integer, room::Integer)
    _dsk_need02(seg, "dskp02")
    n = _dsk_fetch_count("dskp02", "plates", start, room, seg.hdr.np)
    return [SVector{3,Int}(_dsk_plate(ints, seg, start + i - 1)) for i in 1:n]
end

function _dsk_v02(dps::Vector{Float64}, seg::DSKSegment, start::Integer, room::Integer)
    _dsk_need02(seg, "dskv02")
    n = _dsk_fetch_count("dskv02", "vertices", start, room, seg.hdr.nv)
    return [_dsk_vertex(dps, seg, start + i - 1) for i in 1:n]
end

function _dsk_n02(ints::Vector{Int32}, dps::Vector{Float64}, seg::DSKSegment, plid::Integer)
    _dsk_need02(seg, "dskn02")
    i1, i2, i3 = _dsk_plate(ints, seg, plid)
    v1 = _dsk_vertex(dps, seg, i1)
    v2 = _dsk_vertex(dps, seg, i2)
    v3 = _dsk_vertex(dps, seg, i3)
    return _dsk_ucrss(_dsk_vsub(v2, v1), _dsk_vsub(v3, v1))
end

"""
    _dsk_i02(ints, seg, item, start, room) -> Vector{Int}

`dski02`: at most `room` integers of the item `item` (1:14) from element `start`.
"""
function _dsk_i02(ints::Vector{Int32}, seg::DSKSegment, item::Integer, start::Integer,
        room::Integer)
    room <= 0 && throw(SpiceError(:valueoutofrange, "dski02: ROOM was $room; it must be positive"))
    hdr = seg.hdr
    hdr.cgscal < 1 && throw(SpiceError(:valueoutofrange,
        "dski02: the coarse voxel grid scale is $(hdr.cgscal); it must be at least 1"))
    scalar = item == 1 ? hdr.nv : item == 2 ? hdr.np : item == 3 ? hdr.nvxtot :
             item == 5 ? hdr.cgscal : item == 6 ? hdr.voxnpt : item == 7 ? hdr.voxnpl :
             item == 8 ? hdr.vtxnpl : -1
    ibase = seg.dladsc.ibase
    if item in (1, 2, 3, 5, 6, 7, 8)
        # As CSPICE: the value is returned for START = 1; the size is not defined for
        # another START, which is an error.
        start == 1 && return Int[scalar]
        throw(SpiceError(:indexoutofrange, "dski02: START was $start; the item $item has one element"))
    end
    size, off = if item == 4
        3, ibase + 3
    elseif item == 9
        3 * hdr.np, _dsk_off_plates(seg)
    elseif item == 10
        hdr.voxnpt, _dsk_off_voxptr(seg)
    elseif item == 11
        hdr.voxnpl, _dsk_off_voxlst(seg)
    elseif item == 12
        hdr.nv, _dsk_off_vtxptr(seg)
    elseif item == 13
        hdr.vtxnpl, _dsk_off_vtxlst(seg)
    elseif item == 14
        _dsk_ncgr(hdr), _dsk_off_cgrptr(seg)
    else
        throw(SpiceError(:notsupported, "dski02: the keyword $item is not recognized"))
    end
    n = _dsk_fetch_count("dski02", "item $item", start, room, size)
    return Int[ints[off+start+i-1] for i in 1:n]
end

"""
    _dsk_d02(dps, seg, item, start, room) -> Vector{Float64}

`dskd02`: at most `room` doubles of the item `item` (15:19) from element `start`.
"""
function _dsk_d02(dps::Vector{Float64}, seg::DSKSegment, item::Integer, start::Integer,
        room::Integer)
    room <= 0 && throw(SpiceError(:valueoutofrange, "dskd02: ROOM was $room; it must be positive"))
    dbase = seg.dladsc.dbase
    size, off = if item == 19
        3 * seg.hdr.nv, dbase + 34
    elseif item == 15
        24, dbase
    elseif item == 16
        6, dbase + 24
    elseif item == 17
        3, dbase + 30
    elseif item == 18
        1, dbase + 33
    else
        throw(SpiceError(:notsupported, "dskd02: the keyword $item is not recognized"))
    end
    n = _dsk_fetch_count("dskd02", "item $item", start, room, size)
    return Float64[dps[off+start+i-1] for i in 1:n]
end

"""
    dskp02(handle, dladsc, start, room) -> Vector{SVector{3,Int}}
    dskp02([ctx | ks], segment::DSKSegment, start, room) -> the same

At most `room` plates of a type 2 segment from plate `start` (1 is the first), as
`dskp02_c`; a plate is three vertex numbers that start at 1. Throws `SpiceError`
(`:valueoutofrange` for `room <= 0`, `:indexoutofrange` for `start` outside `1:np`).
"""
function dskp02(h::DASHandle, d::DLADescr, start::Integer, room::Integer)
    ints, _, seg = _dsk_open_segment(h, d, "dskp02")
    return _dsk_p02(ints, seg, start, room)
end
dskp02(ks::KernelSet, seg::DSKSegment, start::Integer, room::Integer) =
    _dsk_p02(_dsk_file(ks, seg).ints, seg, start, room)

"""
    dskv02(handle, dladsc, start, room) -> Vector{SVector{3,Float64}}
    dskv02([ctx | ks], segment::DSKSegment, start, room) -> the same

At most `room` vertices of a type 2 segment from vertex `start` (1 is the first), as
`dskv02_c`. Throws `SpiceError` (`:valueoutofrange` for `room <= 0`, `:indexoutofrange`
for `start` outside `1:nv`).
"""
function dskv02(h::DASHandle, d::DLADescr, start::Integer, room::Integer)
    _, dps, seg = _dsk_open_segment(h, d, "dskv02")
    return _dsk_v02(dps, seg, start, room)
end
dskv02(ks::KernelSet, seg::DSKSegment, start::Integer, room::Integer) =
    _dsk_v02(_dsk_file(ks, seg).dps, seg, start, room)

"""
    dskn02(handle, dladsc, plid) -> SVector{3,Float64}
    dskn02([ctx | ks], segment::DSKSegment, plid) -> the same

The outward unit normal of plate `plid` of a type 2 segment, as `dskn02_c`: the unit
vector of `(v2 - v1) x (v3 - v1)`. Throws `SpiceError(:indexoutofrange)` for `plid`
outside `1:np`.
"""
function dskn02(h::DASHandle, d::DLADescr, plid::Integer)
    ints, dps, seg = _dsk_open_segment(h, d, "dskn02")
    return _dsk_n02(ints, dps, seg, plid)
end
function dskn02(ks::KernelSet, seg::DSKSegment, plid::Integer)
    file = _dsk_file(ks, seg)
    return _dsk_n02(file.ints, file.dps, seg, plid)
end

"""
    dski02(handle, dladsc, item, start, room) -> Vector{Int}
    dski02([ctx | ks], segment::DSKSegment, item, start, room) -> the same

At most `room` integers of the item `item` of a type 2 segment from element `start`
(1 is the first), as `dski02_c`. The items: 1 `NV`, 2 `NP`, 3 `NVXTOT`, 4 `VGREXT`,
5 `CGSCAL`, 6 `VOXNPT`, 7 `VOXNPL`, 8 `VTXNPL`, 9 plates, 10 voxel-plate pointers,
11 voxel-plate list, 12 vertex-plate pointers, 13 vertex-plate list, 14 coarse voxel
grid. Throws `SpiceError` (`:valueoutofrange` for `room <= 0`, `:notsupported` for
another item, `:indexoutofrange` for a bad `start`).
"""
function dski02(h::DASHandle, d::DLADescr, item::Integer, start::Integer, room::Integer)
    ints, _, seg = _dsk_open_segment(h, d, "dski02")
    return _dsk_i02(ints, seg, item, start, room)
end
dski02(ks::KernelSet, seg::DSKSegment, item::Integer, start::Integer, room::Integer) =
    _dsk_i02(_dsk_file(ks, seg).ints, seg, item, start, room)

"""
    dskd02(handle, dladsc, item, start, room) -> Vector{Float64}
    dskd02([ctx | ks], segment::DSKSegment, item, start, room) -> the same

At most `room` doubles of the item `item` of a type 2 segment from element `start`
(1 is the first), as `dskd02_c`. The items: 15 DSK descriptor, 16 vertex bounds,
17 voxel grid origin, 18 voxel size, 19 vertices. Throws `SpiceError`
(`:valueoutofrange` for `room <= 0`, `:notsupported` for another item,
`:indexoutofrange` for a bad `start`).
"""
function dskd02(h::DASHandle, d::DLADescr, item::Integer, start::Integer, room::Integer)
    _, dps, seg = _dsk_open_segment(h, d, "dskd02")
    return _dsk_d02(dps, seg, item, start, room)
end
dskd02(ks::KernelSet, seg::DSKSegment, item::Integer, start::Integer, room::Integer) =
    _dsk_d02(_dsk_file(ks, seg).dps, seg, item, start, room)

for f in (:dskp02, :dskv02)
    @eval begin
        $f(ctx::SpiceContext, seg::DSKSegment, start::Integer, room::Integer) =
            $f(kernelset(ctx), seg, start, room)
        $f(seg::DSKSegment, start::Integer, room::Integer) = $f(default_context(), seg, start, room)
    end
end
for f in (:dski02, :dskd02)
    @eval begin
        $f(ctx::SpiceContext, seg::DSKSegment, item::Integer, start::Integer, room::Integer) =
            $f(kernelset(ctx), seg, item, start, room)
        $f(seg::DSKSegment, item::Integer, start::Integer, room::Integer) =
            $f(default_context(), seg, item, start, room)
    end
end
dskn02(ctx::SpiceContext, seg::DSKSegment, plid::Integer) = dskn02(kernelset(ctx), seg, plid)
dskn02(seg::DSKSegment, plid::Integer) = dskn02(default_context(), seg, plid)

"""
    dskx02(handle, dladsc, vertex, raydir) -> Union{Nothing,Tuple{Int,SVector{3,Float64}}}
    dskx02([ctx | ks], segment::DSKSegment, vertex, raydir) -> the same

The plate number and the point where the ray from `vertex` along `raydir` first meets
the plates of one type 2 segment, as `dskx02_c`; `nothing` if there is no intercept.
The ray is in the frame of the segment, relative to the centre of that frame. The file
form uses the default tolerances, and the form with a context uses the tolerances of
the context. Throws `SpiceError(:rayiszerovector)` for a zero `raydir`.
"""
function dskx02(h::DASHandle, d::DLADescr, vertex, raydir)
    ints, dps, seg = _dsk_open_segment(h, d, "dskx02")
    _dsk_need02(seg, "dskx02")
    plid, xpt, found = _dsk_x02(ints, dps, seg, _DSK_DEFAULT_TOLS, _dsk_vec3(vertex),
        _dsk_vec3(raydir))
    return found ? (plid, xpt) : nothing
end
function _dsk_x02_loaded(ks::KernelSet, tols::DSKTolerances, seg::DSKSegment, vertex, raydir)
    file = _dsk_file(ks, seg)
    _dsk_need02(seg, "dskx02")
    plid, xpt, found = _dsk_x02(file.ints, file.dps, seg, tols, _dsk_vec3(vertex),
        _dsk_vec3(raydir))
    return found ? (plid, xpt) : nothing
end
dskx02(ks::KernelSet, seg::DSKSegment, vertex, raydir) =
    _dsk_x02_loaded(ks, _DSK_DEFAULT_TOLS, seg, vertex, raydir)
dskx02(ctx::SpiceContext, seg::DSKSegment, vertex, raydir) =
    _dsk_x02_loaded(kernelset(ctx), dsktolerances(ctx), seg, vertex, raydir)
dskx02(seg::DSKSegment, vertex, raydir) = dskx02(default_context(), seg, vertex, raydir)

# ---------------------------------------------------------------------------------------
# Bodies and surfaces of a file
# ---------------------------------------------------------------------------------------

"The descriptors of all segments of the DSK file `path`, after the checks of `dskobj`."
function _dsk_file_descrs(path::AbstractString, who::AbstractString)
    arch, kind = getfat(path)
    arch == "XFR" && throw(SpiceError(:invalidformat,
        "$who: the file $path has the architecture $arch; it must be a binary DSK file"))
    arch == "DAS" || throw(SpiceError(:invalidarchtype,
        "$who: the file $path has the architecture $arch; DAS is required"))
    kind == "DSK" || throw(SpiceError(:invalidfiletype,
        "$who: the file $path has the type $kind; DSK is required"))
    h = dasopr(path)
    try
        return DSKDescr[_dsk_descr(h.dps, d, h.path) for d in _dla_segments(h.ints, h.path)]
    finally
        dascls(h)
    end
end

"""
    dskobj(path) -> Vector{Int}

The ID codes of the bodies that the segments of the DSK file `path` belong to, sorted
and without repetition, as the set that `dskobj_c` returns. Throws `SpiceError`
(`:invalidformat` for a transfer file, `:invalidarchtype`, `:invalidfiletype`, and the
kinds of `getfat` and `dasopr`).
"""
dskobj(path::AbstractString) = sort!(unique!(Int[d.center for d in _dsk_file_descrs(path, "dskobj")]))

"""
    dsksrf(path, bodyid) -> Vector{Int}

The ID codes of the surfaces of the body `bodyid` in the DSK file `path`, sorted and
without repetition, as the set that `dsksrf_c` returns. Errors as [`dskobj`](@ref).
"""
dsksrf(path::AbstractString, bodyid::Integer) = sort!(unique!(Int[d.surfce
    for d in _dsk_file_descrs(path, "dsksrf") if d.center == bodyid]))

# ---------------------------------------------------------------------------------------
# Ray-surface intercept on the loaded segments
# ---------------------------------------------------------------------------------------

"The body code of the target (a name, or a string that holds an integer, or a code)."
@inline _dsk_target(::KernelSet, target::Integer) = Int32(target)
function _dsk_target(ks::KernelSet, target::AbstractString)
    code = bods2c(ks, target)
    code === nothing && throw(SpiceError(:idcodenotfound, "the target, '$target', is not " *
        "a recognized name for an ephemeris object"))
    return Int32(code)
end

"Centre of the frame `id`, or `nothing` if the frame is unknown. Does not allocate."
@inline function _dsk_frame_center(ks::KernelSet, id::Int32)
    f = get(BUILTIN_FRAMES, id, nothing)
    f === nothing || return f.center
    f = get(ks.frames, id, nothing)
    f === nothing || return f.center
    13000 <= id < 14000 && return Int32(399)
    return nothing
end

@noinline _dsk_bad_frame(kind::Symbol, fixref, target) = throw(SpiceError(kind,
    kind === :idcodenotfound ? "the reference frame $fixref is not recognized" :
    kind === :noframeinfo ? "attributes of the reference frame $fixref could not be determined" :
    "the reference frame $fixref is not centered at the target body $target"))

"""
    _dsk_fixfid(ks, fixref, trgcde, target) -> Int32

The ID of the body-fixed frame `fixref` (a name or an ID) after the checks of
`dskxsi`: the frame is known (`SpiceError(:idcodenotfound)`, `:noframeinfo`) and its
centre is the target (`:invalidframe`).
"""
@inline function _dsk_fixfid(ks::KernelSet, fixref, trgcde::Int32, target)
    fixfid = _dsk_frame_code(ks, fixref)
    fixfid == 0 && _dsk_bad_frame(:idcodenotfound, fixref, target)
    center = _dsk_frame_center(ks, fixfid)
    center === nothing && _dsk_bad_frame(:noframeinfo, fixref, target)
    center == trgcde || _dsk_bad_frame(:invalidframe, fixref, target)
    return fixfid
end
@inline _dsk_frame_code(::KernelSet, fixref::Integer) = Int32(fixref)
_dsk_frame_code(ks::KernelSet, fixref::AbstractString) = Int32(namfrm(ks, fixref))

@noinline _dsk_bad_pri() = throw(SpiceError(:badpriorityspec, "the PRI flag was set to " *
    "true; prioritized searches are not supported, as in this version of CSPICE"))

"""
    dskxsi([ctx | ks], pri, target, srflst, et, fixref, vertex, raydir)
        -> Union{Nothing,Tuple{SVector{3,Float64},DSKSegment,Int}}

The first point where the ray from `vertex` along `raydir` meets the surface that the
loaded DSK segments of `target` define, as `dskxsi_c`: `(xpt, segment, plid)`, or
`nothing` if there is no intercept.

- `pri` must be `false`: all segments have the same priority and the nearest intercept
  wins (`SpiceError(:badpriorityspec)` otherwise, as CSPICE).
- `target` is a body name or ID; `srflst` holds surface IDs (a vector or a tuple), and
  an empty list accepts every surface of the body.
- `et` is the epoch (seconds past J2000 TDB) that selects the segments by coverage and
  at which the frames are evaluated.
- `fixref` is the name or the ID of a body-fixed frame centred on the target; `vertex`
  and `raydir` are in that frame, `vertex` relative to the target centre, and so is
  `xpt`.

`segment` is a [`DSKSegment`](@ref): it replaces the handle, the DLA descriptor and the
DSK descriptor of CSPICE (`segment.dladsc`, `segment.dskdsc`), and `plid` replaces
`ic[1]`; `dc` of CSPICE has no content for type 2. The access calls take the segment.

Contract: the call reads one `KernelSet` snapshot and one tolerance snapshot and writes
nothing, so that many threads can call it while another thread loads kernels. With an
integer `target`, an integer `fixref` and `SVector` arguments it does not allocate.
Throws `SpiceError` (`:badpriorityspec`, `:idcodenotfound`, `:noframeinfo`,
`:invalidframe`, `:nodsksegments`, `:dskdatanotfound`, `:rayiszerovector`).
"""
function dskxsi(ks::KernelSet, pri::Bool, target::NameOrID, srflst, et::Real, fixref::NameOrID,
        vertex, raydir)
    return _dsk_xsi(ks, _DSK_DEFAULT_TOLS, pri, target, srflst, et, fixref, vertex, raydir)
end
function dskxsi(ctx::SpiceContext, pri::Bool, target::NameOrID, srflst, et::Real,
        fixref::NameOrID, vertex, raydir)
    return _dsk_xsi(kernelset(ctx), dsktolerances(ctx), pri, target, srflst, et, fixref,
        vertex, raydir)
end
dskxsi(pri::Bool, target::NameOrID, srflst, et::Real, fixref::NameOrID, vertex, raydir) =
    dskxsi(default_context(), pri, target, srflst, et, fixref, vertex, raydir)

@inline function _dsk_xsi(ks::KernelSet, tols::DSKTolerances, pri::Bool, target, srflst,
        et::Real, fixref, vertex, raydir)
    pri && _dsk_bad_pri()
    trgcde = _dsk_target(ks, target)
    fixfid = _dsk_fixfid(ks, fixref, trgcde, target)
    xpt, idx, plid, found = _dsk_bux(ks, tols, trgcde, srflst, Float64(et), fixfid,
        _dsk_vec3(vertex), _dsk_vec3(raydir))
    found || return nothing
    return xpt, (@inbounds ks.dsk[trgcde][idx]), plid
end

"""
    dskxv!(xptarr, fndarr, [ctx | ks], pri, target, srflst, et, fixref, vtxarr, dirarr) -> (xptarr, fndarr)

As [`dskxv`](@ref), with the results stored in `xptarr` (a vector of `SVector{3,Float64}`)
and `fndarr` (a vector of `Bool`). Precondition: both have one element for each ray.
With integer `target` and `fixref` and vectors of `SVector{3,Float64}` the call does
not allocate.
"""
function dskxv!(xptarr::AbstractVector{Vec3{Float64}}, fndarr::AbstractVector{Bool},
        ks::KernelSet, pri::Bool, target::NameOrID, srflst, et::Real, fixref::NameOrID, vtxarr,
        dirarr)
    return _dsk_xv!(xptarr, fndarr, ks, _DSK_DEFAULT_TOLS, pri, target, srflst, et, fixref,
        vtxarr, dirarr)
end
function dskxv!(xptarr::AbstractVector{Vec3{Float64}}, fndarr::AbstractVector{Bool},
        ctx::SpiceContext, pri::Bool, target::NameOrID, srflst, et::Real, fixref::NameOrID,
        vtxarr, dirarr)
    return _dsk_xv!(xptarr, fndarr, kernelset(ctx), dsktolerances(ctx), pri, target, srflst,
        et, fixref, vtxarr, dirarr)
end
dskxv!(xptarr::AbstractVector{Vec3{Float64}}, fndarr::AbstractVector{Bool}, pri::Bool,
    target::NameOrID, srflst, et::Real, fixref::NameOrID, vtxarr, dirarr) =
    dskxv!(xptarr, fndarr, default_context(), pri, target, srflst, et, fixref, vtxarr, dirarr)

function _dsk_xv!(xptarr::AbstractVector{Vec3{Float64}}, fndarr::AbstractVector{Bool},
        ks::KernelSet, tols::DSKTolerances, pri::Bool, target, srflst, et::Real, fixref,
        vtxarr, dirarr)
    pri && _dsk_bad_pri()
    nrays = _dsk_count3(vtxarr)
    nrays < 1 && throw(SpiceError(:invalidcount, "dskxv: the number of rays must be positive"))
    (_dsk_count3(dirarr) == nrays && length(xptarr) == nrays && length(fndarr) == nrays) ||
        throw(ArgumentError("dskxv: the ray arrays and the result arrays must have one element for each ray"))
    trgcde = _dsk_target(ks, target)
    fixfid = _dsk_fixfid(ks, fixref, trgcde, target)
    t = Float64(et)
    for i in 1:nrays                     # invariant: the results of the rays 1 .. i-1 are stored
        xpt, _, _, found = _dsk_bux(ks, tols, trgcde, srflst, t, fixfid, _dsk_get3(vtxarr, i),
            _dsk_get3(dirarr, i))
        xptarr[begin+i-1] = xpt
        fndarr[begin+i-1] = found
    end
    return xptarr, fndarr
end

"""
    dskxv([ctx | ks], pri, target, srflst, et, fixref, vtxarr, dirarr) -> (xptarr, fndarr)

The ray-surface intercepts of a set of rays on the loaded DSK segments of `target`, as
`dskxv_c`. `vtxarr` and `dirarr` are 3 x n matrices or vectors of n 3-vectors.
`xptarr[i]` is the intercept of ray `i` (an `SVector{3,Float64}`, zero if there is
none) and `fndarr[i]` tells whether it exists. The other arguments and the errors are
those of [`dskxsi`](@ref); `SpiceError(:invalidcount)` for no ray.
"""
function dskxv(ks::Union{KernelSet,SpiceContext}, pri::Bool, target::NameOrID, srflst, et::Real,
        fixref::NameOrID, vtxarr, dirarr)
    n = _dsk_count3(vtxarr)
    return dskxv!(Vector{Vec3{Float64}}(undef, n), Vector{Bool}(undef, n), ks, pri, target,
        srflst, et, fixref, vtxarr, dirarr)
end
dskxv(pri::Bool, target::NameOrID, srflst, et::Real, fixref::NameOrID, vtxarr, dirarr) =
    dskxv(default_context(), pri, target, srflst, et, fixref, vtxarr, dirarr)

# ---------------------------------------------------------------------------------------
# Plate utilities
# ---------------------------------------------------------------------------------------

"""
    pltexp(iverts, delta) -> SVector{3,SVector{3,Float64}}

Expand the plate with the three vertices `iverts` (three 3-vectors, or the columns of
a 3 x 3 matrix) about its centroid by the fraction `delta`, as `pltexp_c`. The result
is in the same plane and has the same orientation.
"""
function pltexp(iverts, delta::Real)
    _dsk_count3(iverts) == 3 || throw(ArgumentError("pltexp: three vertices are required"))
    a, b, c = _dsk_pltexp(_dsk_get3(iverts, 1), _dsk_get3(iverts, 2), _dsk_get3(iverts, 3),
        Float64(delta))
    return SVector{3,Vec3{Float64}}(a, b, c)
end

"""
    pltnrm(v1, v2, v3) -> SVector{3,Float64}

The outward normal `(v2 - v1) x (v3 - v2)` of the plate with the vertices `v1`, `v2`,
`v3`, as `pltnrm_c`. The result does not have unit length; it is zero for a degenerate
plate.
"""
pltnrm(v1, v2, v3) = _dsk_pltnrm(_dsk_vec3(v1), _dsk_vec3(v2), _dsk_vec3(v3))

"""
    pltnp(point, v1, v2, v3) -> (pnear, dist)

The point of the plate with the vertices `v1`, `v2`, `v3` that is nearest to `point`,
and its distance from `point`, as `pltnp_c`. A degenerate plate (a segment or a point)
is accepted.
"""
pltnp(point, v1, v2, v3) = _dsk_pltnp(_dsk_vec3(point), _dsk_vec3(v1), _dsk_vec3(v2), _dsk_vec3(v3))

"Throws `SpiceError(:indexoutofrange)` unless every vertex number of `plates` is in `1:nv`."
function _dsk_check_plates(who::AbstractString, plates, np::Int, nv::Int)
    for i in 1:np                        # invariant: the plates 1 .. i-1 are valid
        p = plates isa AbstractMatrix ? (plates[1, i], plates[2, i], plates[3, i]) : plates[begin+i-1]
        length(p) == 3 || throw(ArgumentError("$who: plate $i does not have 3 vertex numbers"))
        for j in 1:3
            (1 <= p[j] <= nv) || throw(SpiceError(:indexoutofrange, "$who: vertex number " *
                "$(p[j]) of plate $i is not in the range 1:$nv"))
        end
    end
    return nothing
end

"The three vertex numbers of plate `i` of a list of plates."
@inline _dsk_plate3(plates::AbstractMatrix, i::Int) = (Int(plates[1, i]), Int(plates[2, i]), Int(plates[3, i]))
@inline function _dsk_plate3(plates::AbstractVector, i::Int)
    p = plates[begin+i-1]
    return Int(p[begin]), Int(p[begin+1]), Int(p[begin+2])
end

"""
    pltar(vrtces, plates) -> Float64

The total area of the plates, as `pltar_c`. `vrtces` and `plates` are 3 x n matrices or
vectors of 3-vectors; vertex numbers start at 1. No plate gives 0. Throws `SpiceError`
(`:toofewvertices` for fewer than 3 vertices when there is a plate, `:indexoutofrange`
for a vertex number outside `1:nv`).
"""
function pltar(vrtces, plates)
    nv = _dsk_count3(vrtces)
    np = _dsk_count3(plates)
    np == 0 && return 0.0
    nv < 3 && throw(SpiceError(:toofewvertices,
        "pltar: at least 3 vertices are needed, but NV = $nv"))
    _dsk_check_plates("pltar", plates, np, nv)
    area = 0.0
    for i in 1:np                        # invariant: area is the sum over the plates 1 .. i-1
        i1, i2, i3 = _dsk_plate3(plates, i)
        v1 = _dsk_get3(vrtces, i1)
        v2 = _dsk_get3(vrtces, i2)
        v3 = _dsk_get3(vrtces, i3)
        area += _dsk_vnorm(_dsk_vcrss(_dsk_vsub(v2, v1), _dsk_vsub(v3, v2))) * 0.5
    end
    return area
end

"""
    pltvol(vrtces, plates) -> Float64

The volume of the region that the plates enclose, as `pltvol_c`. Precondition: the
plates form a closed surface with outward normals; the function does not check it.
Throws `SpiceError` (`:toofewvertices` for fewer than 4 vertices, `:toofewplates` for
fewer than 4 plates, `:indexoutofrange` for a vertex number outside `1:nv`).
"""
function pltvol(vrtces, plates)
    nv = _dsk_count3(vrtces)
    np = _dsk_count3(plates)
    nv < 4 && throw(SpiceError(:toofewvertices,
        "pltvol: at least 4 vertices are needed, but NV = $nv"))
    np < 4 && throw(SpiceError(:toofewplates, "pltvol: at least 4 plates are needed, but NP = $np"))
    _dsk_check_plates("pltvol", plates, np, nv)
    vol = 0.0
    for i in 1:np                        # invariant: vol is the sum over the plates 1 .. i-1
        i1, i2, i3 = _dsk_plate3(plates, i)
        a = _dsk_get3(vrtces, i1)
        b = _dsk_get3(vrtces, i2)
        c = _dsk_get3(vrtces, i3)
        det = a[1] * (b[2] * c[3] - c[2] * b[3]) - b[1] * (a[2] * c[3] - c[2] * a[3]) +
              c[1] * (a[2] * b[3] - b[2] * a[3])
        vol += det / 6.0
    end
    return vol
end
