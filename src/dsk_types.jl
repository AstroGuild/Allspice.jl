# Data types of the DLA and DSK layers. This file comes before `types.jl` because
# `Loaded` and `KernelSet` hold these types. Every type here is immutable. The segment table of a `KernelSet` is built here too,
# because the empty set of the default context needs it when the package loads.

"""
    DLADescr

Descriptor of one DLA segment (a segment of a DSK file), as `SpiceDLADescr` of CSPICE.

- `bwdptr`, `fwdptr`: addresses of the descriptors of the segment before and of the
  segment after in the list, `-1` at the ends.
- `ibase`, `isize`: base address and element count of the segment in the integer address
  space; `dbase`, `dsize` and `cbase`, `csize` are those of the double and of the
  character address space. A base address is the address before the first element.

[`dlabfs`](@ref), [`dlafns`](@ref), [`dlabbs`](@ref) and [`dlafps`](@ref) return a
descriptor.
"""
struct DLADescr
    bwdptr::Int
    fwdptr::Int
    ibase::Int
    isize::Int
    dbase::Int
    dsize::Int
    cbase::Int
    csize::Int
end

"""
    DSKDescr

Descriptor of one DSK segment, as `SpiceDSKDescr` of CSPICE.

- `surfce`, `center`: surface ID and body ID.
- `dclass`: data class, 1 for a single-valued surface and 2 for a general surface.
- `dtype`: data type, 2 for triangular plates.
- `frmcde`: frame ID.
- `corsys`: coordinate system, 1 latitudinal, 2 cylindrical, 3 rectangular,
  4 planetodetic.
- `corpar`: the ten coordinate system parameters; for system 4 the first two are the
  equatorial radius (km) and the flattening.
- `co1min`, `co1max`, `co2min`, `co2max`, `co3min`, `co3max`: bounds of the three
  coordinates (radians and km).
- `start`, `stop`: coverage, seconds past J2000 TDB.

[`dskgd`](@ref) returns a descriptor.
"""
struct DSKDescr
    surfce::Int
    center::Int
    dclass::Int
    dtype::Int
    frmcde::Int
    corsys::Int
    corpar::NTuple{10,Float64}
    co1min::Float64
    co1max::Float64
    co2min::Float64
    co2max::Float64
    co3min::Float64
    co3max::Float64
    start::Float64
    stop::Float64
end

"""
    DSK02Header

The fixed-size items of a type 2 DSK segment, as `dskb02` returns them: vertex and
plate counts, total voxel count, vertex bounds (xmin, xmax, ymin, ymax, zmin, zmax),
voxel size, voxel grid origin, voxel grid extents, coarse voxel scale, and the sizes
of the vertex-plate list, the voxel-plate pointer array and the voxel-plate list.
"""
struct DSK02Header
    nv::Int
    np::Int
    nvxtot::Int
    vtxbds::NTuple{6,Float64}
    voxsiz::Float64
    voxori::SVector{3,Float64}
    vgrext::NTuple{3,Int}
    cgscal::Int
    vtxnpl::Int
    voxnpt::Int
    voxnpl::Int
end

"""
    DSKSegment

One segment of a loaded DSK file. [`dskxsi`](@ref) returns it in place of the handle, the
DLA descriptor and the DSK descriptor that `dskxsi_c` returns.

The access calls (`dskgd`, `dskb02`, `dskz02`, `dskp02`, `dskv02`, `dskn02`, `dski02`,
`dskd02`, `dskx02`) take the segment. The calls that read data of the file also take the
`KernelSet` or the `SpiceContext` that holds the file, and throw
`SpiceError(:nosuchhandle)` if the file is not loaded there any more.

A caller reads these fields:
- `dladsc`, `dskdsc`: the [`DLADescr`](@ref) and the [`DSKDescr`](@ref) of the segment.
- `boxctr`, `boxrad`: center and radius (km) of the bounding sphere of the segment, in the
  frame of the segment. `boxrad` is NaN unless the coordinate system of the segment is
  latitudinal, rectangular or planetodetic.
"""
struct DSKSegment
    file::Int32
    dladsc::DLADescr
    dskdsc::DSKDescr
    boxctr::SVector{3,Float64}
    boxrad::Float64
    hdr::DSK02Header
end

"""
    DSKFile

One loaded DSK file. `ints[a]` and `dps[a]` are the integer and the double at DAS
address `a` (memory-mapped when the file is segregated and of native byte order, an
owned converted copy otherwise). `segments` are in file order, so a later segment comes
first in the search order. `kind` is always `:DSK`. Nothing in a `DSKFile` is modified
after `furnsh` has built it.
"""
struct DSKFile
    path::String
    kind::Symbol
    ints::Vector{Int32}
    dps::Vector{Float64}
    segments::Vector{DSKSegment}
end

"""
    DSKTolerances

Immutable snapshot of the six DSK tolerances of `dskgtl`/`dskstl`, in keyword order:
`xfract` (1, plate expansion fraction), `sgreed` (2, segment boundary margin),
`sgpadm` (3, segment pad margin), `ptmemm` (4, surface point membership margin),
`angmrg` (5, angular rounding margin, fixed), `lonali` (6, longitude alias margin,
fixed). Declared `mutable` with `const` fields only so that it lives behind one
pointer, as `TimeDefaults` does.
"""
mutable struct DSKTolerances
    const xfract::Float64
    const sgreed::Float64
    const sgpadm::Float64
    const ptmemm::Float64
    const angmrg::Float64
    const lonali::Float64
end

"The tolerances of CSPICE at start."
DSKTolerances() = DSKTolerances(1e-10, 1e-8, 1e-10, 1e-7, 1e-12, 1e-12)

_dsk_withfile(s::DSKSegment, file::Integer) = DSKSegment(Int32(file), s.dladsc, s.dskdsc, s.boxctr,
    s.boxrad, s.hdr)

"""
    _dsk_segment_table(files) -> Dict{Int32,Vector{DSKSegment}}

The segments of the loaded DSK `files` (load order) by center body, in the search
order of CSPICE: later file first, and inside one file later segment first. Each
segment gets the number of its file.
"""
function _dsk_segment_table(files::Vector{DSKFile})
    table = Dict{Int32,Vector{DSKSegment}}()
    for i in length(files):-1:1
        segs = files[i].segments
        for j in length(segs):-1:1
            seg = _dsk_withfile(segs[j], i)
            push!(get!(() -> DSKSegment[], table, Int32(seg.dskdsc.center)), seg)
        end
    end
    return table
end
