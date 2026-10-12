# The DLA layer (DAS linked array): a doubly linked list of segments in a DAS file.
#
# The first three integers of the file are the DLA format version (1000000) and the
# addresses of the first and of the last segment descriptor (-1 if the file has no
# segment). A segment descriptor is eight integers: the backward and forward pointers,
# and the base address and size of the segment in the integer, double and character
# address spaces.

const _DLA_VERSION = 1000000
const _DLA_NULPTR = -1
const _DLA_DSCSIZ = 8

"The DLA descriptor in `ints[addr:addr+7]`; throws `SpiceError(:dasnosuchaddress)` if that range is not in `ints`."
function _dla_descr(ints::Vector{Int32}, addr::Integer, path::AbstractString)
    (1 <= addr && addr + _DLA_DSCSIZ - 1 <= length(ints)) || throw(SpiceError(
        :dasnosuchaddress, "the DLA descriptor address $addr is not inside the " *
        "$(length(ints)) integers of $path"))
    @inbounds return DLADescr(ints[addr], ints[addr+1], ints[addr+2], ints[addr+3],
        ints[addr+4], ints[addr+5], ints[addr+6], ints[addr+7])
end

"Throws `SpiceError(:dasnosuchaddress)` unless `ints` holds the three DLA file pointers."
function _dla_check(ints::Vector{Int32}, path::AbstractString)
    length(ints) >= 3 || throw(SpiceError(:dasnosuchaddress,
        "the file $path has $(length(ints)) integers; a DLA file has at least 3"))
    return nothing
end

"""
    _dla_segments(ints, path) -> Vector{DLADescr}

All segment descriptors of the DLA file with the integer address space `ints`, in file
order (the forward list). Throws `SpiceError(:baddladescriptor)` if the list has a
cycle.
"""
function _dla_segments(ints::Vector{Int32}, path::AbstractString)
    _dla_check(ints, path)
    out = DLADescr[]
    addr = Int(ints[2])
    while addr != _DLA_NULPTR
        # invariant: `out` holds the descriptors before `addr`, each one visited once
        length(out) * _DLA_DSCSIZ <= length(ints) || throw(SpiceError(:baddladescriptor,
            "the forward list of DLA segments has a cycle in $path"))
        d = _dla_descr(ints, addr, path)
        push!(out, d)
        addr = d.fwdptr
    end
    return out
end

"""
    dlaopn(path, ftype, ifname, ncomch) -> DASHandle

Create the new DLA file `path` of type `ftype` with the internal file name `ifname`
and room for `ncomch` comment characters, and open it for write, as `dlaopn` of
SPICELIB. Preconditions: those of [`dasonw`](@ref), and `ncomch >= 0`
(`SpiceError(:badrecordcount)`). Postcondition: the file has no segment.
"""
function dlaopn(path::AbstractString, ftype::AbstractString, ifname::AbstractString,
        ncomch::Integer)
    ncomch < 0 && throw(SpiceError(:badrecordcount,
        "dlaopn: the requested number of comment characters must be non-negative but was $ncomch"))
    ncomr = ncomch > 0 ? (Int(ncomch) - 1) ÷ _DAS_RECLEN + 1 : 0
    h = dasonw(path, ftype, ifname, ncomr)
    dasadi(h, Int32[_DLA_VERSION, _DLA_NULPTR, _DLA_NULPTR])
    return h
end

"""
    dlabns(handle) -> Nothing

Begin a new segment at the end of the DLA file, as `dlabns` of SPICELIB: append its
descriptor and link it to the list. Precondition: write access. Postcondition: the
segment is the last one and has the size zero in every address space until
[`dlaens`](@ref) is called.
"""
function dlabns(h::DASHandle)
    _das_check_write(h, "dlabns")
    _dla_check(h.ints, h.path)
    first, last = Int(h.ints[2]), Int(h.ints[3])
    lastc, lastd, lasti = daslla(h)
    this = lasti + 1
    dasadi(h, Int32[last, _DLA_NULPTR, lasti + _DLA_DSCSIZ, _DLA_NULPTR, lastd,
        _DLA_NULPTR, lastc, _DLA_NULPTR])
    last != _DLA_NULPTR && (h.ints[last+1] = this)
    first == _DLA_NULPTR && (h.ints[2] = this)
    h.ints[3] = this
    return nothing
end

"""
    dlaens(handle) -> Nothing

End the segment that [`dlabns`](@ref) began, as `dlaens` of SPICELIB: store the sizes
of the segment in its descriptor. Precondition: write access, and a segment exists.
"""
function dlaens(h::DASHandle)
    _das_check_write(h, "dlaens")
    _dla_check(h.ints, h.path)
    this = Int(h.ints[3])
    d = _dla_descr(h.ints, this, h.path)
    lastc, lastd, lasti = daslla(h)
    h.ints[this+3] = lasti - d.ibase
    h.ints[this+5] = lastd - d.dbase
    h.ints[this+7] = lastc - d.cbase
    return nothing
end

"""
    dlabfs(handle) -> Union{DLADescr,Nothing}

Begin a forward segment search: the descriptor of the first segment of the DLA file,
or `nothing` if the file has no segment, as `dlabfs_c` (and SPICE.jl). The search has
no hidden state: continue with [`dlafns`](@ref) on the descriptor.
"""
function dlabfs(h::DASHandle)
    _das_check_open(h, "dlabfs")
    _dla_check(h.ints, h.path)
    this = Int(h.ints[2])
    return this == _DLA_NULPTR ? nothing : _dla_descr(h.ints, this, h.path)
end

"""
    dlafns(handle, dladsc) -> Union{DLADescr,Nothing}

The descriptor of the segment that follows `dladsc`, or `nothing` if `dladsc` is the
last segment, as `dlafns_c`.
"""
function dlafns(h::DASHandle, d::DLADescr)
    _das_check_open(h, "dlafns")
    return d.fwdptr == _DLA_NULPTR ? nothing : _dla_descr(h.ints, d.fwdptr, h.path)
end

"""
    dlabbs(handle) -> Union{DLADescr,Nothing}

Begin a backward segment search: the descriptor of the last segment of the DLA file,
or `nothing` if the file has no segment, as `dlabbs_c`.
"""
function dlabbs(h::DASHandle)
    _das_check_open(h, "dlabbs")
    _dla_check(h.ints, h.path)
    this = Int(h.ints[3])
    return this == _DLA_NULPTR ? nothing : _dla_descr(h.ints, this, h.path)
end

"""
    dlafps(handle, dladsc) -> Union{DLADescr,Nothing}

The descriptor of the segment that precedes `dladsc`, or `nothing` if `dladsc` is the
first segment, as `dlafps_c`.
"""
function dlafps(h::DASHandle, d::DLADescr)
    _das_check_open(h, "dlafps")
    return d.bwdptr == _DLA_NULPTR ? nothing : _dla_descr(h.ints, d.bwdptr, h.path)
end
