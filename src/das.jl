# The DAS container (direct access, segregated): the file architecture under DSK files.
#
# A DAS file is a sequence of 1024-byte records: the file record, reserved records,
# comment records, and then directory records and data records. It holds three address
# spaces, each numbered from 1: characters (1024 in a record), doubles (128 in a record)
# and 32-bit integers (256 in a record). A directory record lists runs ("clusters") of
# data records of one type; a file that was closed with `dascls` is segregated, which
# means that it has one directory record and then all character records, all double
# records and all integer records.
#
# CSPICE keeps a process-wide table of open DAS files and buffers their records. Here
# an opened file is an object, `DASHandle`, that the caller owns and that holds the
# three address spaces as plain arrays: memory-mapped for a read handle on a
# segregated file of native byte order, an owned copy in every other case. A write
# handle works on the arrays in memory, and `dascls` writes the complete, segregated
# file and moves it into place in one step, so a reader that has the old file mapped
# keeps a consistent image. There is no global state. One handle is not safe for
# concurrent use: a caller that shares a handle between tasks must hold a lock.

const _DAS_RECLEN = 1024
const _DAS_NW = (1024, 128, 256)          # words in a record of type 1 (char), 2 (dp), 3 (int)
const _DAS_NEXT = (2, 3, 1)
const _DAS_PREV = (3, 1, 2)
const _DAS_NATIVE_LITTLE = ENDIAN_BOM == 0x04030201
const _DAS_NATIVE_BFF = _DAS_NATIVE_LITTLE ? "LTL-IEEE" : "BIG-IEEE"

"""
    DASHandle

One DAS file opened with `dasopr`, `dasopw`, `dasonw` or `dskopn`. It replaces the
integer handle of CSPICE and holds what CSPICE keeps in global tables.

Invariants: `chars`, `dps` and `ints` are the three address spaces (`ints[a]` is the
integer at DAS address `a`); `length(comments) == 1024 * ncomr` and
`0 <= ncomc <= length(comments)`; the arrays of a handle that is not `writable` are
never modified (they can be memory-mapped). For a read handle `free`, `lastrc` and
`lastwd` describe the file as it is on disk; for a write handle they describe the
segregated file that `dascls` will write. A handle must not be used by two tasks at
the same time without a lock that the caller holds; different handles are independent.
"""
mutable struct DASHandle
    const path::String
    const writable::Bool
    isopen::Bool
    const idword::String
    ifname::String
    const nresvr::Int
    const nresvc::Int
    ncomr::Int
    ncomc::Int
    const resv::Vector{UInt8}        # the reserved records, kept as they are
    comments::Vector{UInt8}          # the comment records
    chars::Vector{UInt8}
    dps::Vector{Float64}
    ints::Vector{Int32}
    free::Int
    lastrc::NTuple{3,Int}
    lastwd::NTuple{3,Int}
end

function Base.show(io::IO, h::DASHandle)
    print(io, "DASHandle(", repr(h.path), ", ", h.writable ? "write" : "read",
        h.isopen ? "" : ", closed", ")")
end
Base.isopen(h::DASHandle) = h.isopen
Base.close(h::DASHandle) = dascls(h)

"Throws `SpiceError(:dasnosuchhandle)` unless `h` is open."
@inline function _das_check_open(h::DASHandle, who::AbstractString)
    h.isopen || throw(SpiceError(:dasnosuchhandle, "$who: the DAS file $(h.path) is closed"))
    return nothing
end

"Throws unless `h` is open for writing (`:dasnosuchhandle`, `:dasinvalidaccess`)."
@inline function _das_check_write(h::DASHandle, who::AbstractString)
    _das_check_open(h, who)
    h.writable || throw(SpiceError(:dasinvalidaccess,
        "$who: the DAS file $(h.path) is not open for write"))
    return nothing
end

# ---------------------------------------------------------------------------------------
# Reading a file image
# ---------------------------------------------------------------------------------------

"32-bit integer at 0-based byte offset `off` of `raw`; `swap` tells that the file has the other byte order."
@inline function _das_i32(raw::AbstractVector{UInt8}, off::Int, swap::Bool)
    v = UInt32(raw[off+1]) | (UInt32(raw[off+2]) << 8) | (UInt32(raw[off+3]) << 16) |
        (UInt32(raw[off+4]) << 24)
    v = _DAS_NATIVE_LITTLE ? v : bswap(v)       # `v` is now the stored word, native layout
    return Int(reinterpret(Int32, swap ? bswap(v) : v))
end

"""
    _das_probe(path, who) -> (raw_file_record, swap)

Check that `path` is an existing DAS file that this machine can read, and return its
file record and whether its byte order is the other one. Throws `SpiceError` with kind
`:blankfilename`, `:filenotfound`, `:filereadfailed`, `:unknownfilarc`,
`:filarchmismatch`, `:ftpxfererror` or `:unknownbff`. A file from before the
binary-file-format string (bytes 85 to 92 blank or zero) is taken as native.
"""
function _das_probe(path::AbstractString, who::AbstractString)
    isempty(strip(path)) && throw(SpiceError(:blankfilename, "$who: the file name is blank"))
    isfile(path) || throw(SpiceError(:filenotfound, "$who: the file $path does not exist"))
    raw = open(io -> read(io, _DAS_RECLEN), path, "r")
    length(raw) == _DAS_RECLEN || throw(SpiceError(:filereadfailed,
        "$who: $path is too short ($(length(raw)) bytes) to hold a DAS file record"))
    arch, _ = idw2at(String(raw[1:8]))
    arch == "?" && throw(SpiceError(:unknownfilarc,
        "$who: the file $path has an unidentified file architecture"))
    arch == "DAS" || throw(SpiceError(:filarchmismatch,
        "$who: the file $path has the architecture $arch, not DAS"))
    if raw[700:705] == _DAF_FTPSTR[1:6] && raw[700:727] != _DAF_FTPSTR
        throw(SpiceError(:ftpxfererror, "$who: the FTP validation string of $path is " *
            "damaged; the file was most likely moved in text mode"))
    end
    bff = rstrip(String(raw[85:92]), [' ', '\0'])
    little = if bff == "LTL-IEEE"
        true
    elseif bff == "BIG-IEEE"
        false
    elseif isempty(bff)
        _DAS_NATIVE_LITTLE
    else
        throw(SpiceError(:unknownbff,
            "$who: the binary file format $(repr(bff)) of $path is not known"))
    end
    return raw, little != _DAS_NATIVE_LITTLE
end

"""
    _das_extent(io, raw, T, clusters, n, swap, copy) -> Vector{T}

The first `n` words of one DAS address space whose data records are the `clusters`
(first record, record count), in address order. Contract: the clusters lie inside
`raw` and hold at least `n` words. The result is memory-mapped from `io` if the
clusters are adjacent, the byte order is native and `copy` is false; else it is an
owned array in native byte order.
"""
function _das_extent(io::IOStream, raw::Vector{UInt8}, ::Type{T},
        clusters::Vector{NTuple{2,Int}}, n::Int, swap::Bool, copy::Bool) where {T}
    n == 0 && return T[]
    adjacent = true
    for k in 2:length(clusters)          # invariant: clusters 1 .. k-1 are one run of records
        adjacent &= clusters[k-1][1] + clusters[k-1][2] == clusters[k][1]
    end
    if adjacent && !swap && !copy
        return Mmap.mmap(io, Vector{T}, (n,), (clusters[1][1] - 1) * _DAS_RECLEN)
    end
    out = Vector{T}(undef, n)
    nbytes = n * sizeof(T)
    done = 0
    GC.@preserve out raw begin
        dst = Ptr{UInt8}(pointer(out))
        for (first, count) in clusters   # invariant: `done` bytes of `out` are filled
            take = min(count * _DAS_RECLEN, nbytes - done)
            take <= 0 && break
            unsafe_copyto!(dst + done, pointer(raw, (first - 1) * _DAS_RECLEN + 1), take)
            done += take
        end
    end
    if swap && sizeof(T) > 1
        map!(bswap, out, out)
    end
    return out
end

"""
    _das_load(path, who, writable) -> DASHandle

Read the DAS file `path` and return an open handle on it. Contract: `path` is a DAS
file of either byte order whose directory records are consistent; the handle holds the
complete address spaces (owned copies if `writable`, else memory-mapped where the file
allows it). Throws the errors of `_das_probe`, and `SpiceError(:baddasdirectory)` for
a directory that points outside the file, has a cycle, or does not cover the last
addresses that it states.
"""
function _das_load(path::AbstractString, who::AbstractString, writable::Bool)
    frec, swap = _das_probe(path, who)
    idword = String(frec[1:8])
    ifname = String(rstrip(String(frec[9:68]), [' ', '\0']))
    nresvr = _das_i32(frec, 68, swap)
    nresvc = _das_i32(frec, 72, swap)
    ncomr = _das_i32(frec, 76, swap)
    ncomc = _das_i32(frec, 80, swap)
    bad(msg) = throw(SpiceError(:baddasdirectory, "$who: $msg in $path"))
    return open(path, "r") do io
        nrecords = Int(filesize(io) ÷ _DAS_RECLEN)
        raw = Mmap.mmap(io, Vector{UInt8}, (nrecords * _DAS_RECLEN,))
        (nresvr >= 0 && ncomr >= 0 && nresvr + ncomr + 2 <= nrecords) ||
            bad("the file record gives $nresvr reserved and $ncomr comment records, " *
                "but the file has $nrecords records")
        (0 <= ncomc <= ncomr * _DAS_RECLEN) ||
            bad("the comment area of $ncomr records cannot hold $ncomc characters")
        clusters = (NTuple{2,Int}[], NTuple{2,Int}[], NTuple{2,Int}[])
        lastla = [0, 0, 0]
        lastrc = [0, 0, 0]
        lastwd = [0, 0, 0]
        lastrec = 0                      # the last directory or data record that is in use
        rec = nresvr + ncomr + 2
        visited = 0
        while rec > 0
            # invariant: `visited` directory records were read, none twice while
            # `visited <= nrecords`; `clusters` holds their data records in address order
            visited += 1
            visited <= nrecords || bad("the chain of directory records has a cycle")
            rec <= nrecords || bad("directory record $rec is outside the file")
            base = (rec - 1) * _DAS_RECLEN
            word(k) = _das_i32(raw, base + 4 * (k - 1), swap)
            lastrec = max(lastrec, rec)
            for t in 1:3
                hi = word(2t + 2)
                hi > 0 && (lastla[t] = hi)
            end
            typ = word(9)
            if 1 <= typ <= 3
                prvtyp = _DAS_PREV[typ]
                datarec = rec + 1
                pos = 10
                while pos <= 256 && word(pos) != 0
                    # invariant: descriptors 10 .. pos-1 are in `clusters`
                    count = word(pos)
                    curtyp = count > 0 ? _DAS_NEXT[prvtyp] : _DAS_PREV[prvtyp]
                    count = abs(count)
                    datarec + count - 1 <= nrecords ||
                        bad("directory record $rec describes records past the end of the file")
                    push!(clusters[curtyp], (datarec, count))
                    lastrc[curtyp] = rec
                    lastwd[curtyp] = pos
                    datarec += count
                    lastrec = max(lastrec, datarec - 1)
                    prvtyp = curtyp
                    pos += 1
                end
            elseif typ != 0
                bad("directory record $rec has the data type $typ")
            end
            rec = word(2)
        end
        for t in 1:3
            room = sum(c -> c[2], clusters[t]; init=0) * _DAS_NW[t]
            lastla[t] <= room || bad("the last address $(lastla[t]) of type $t is not " *
                "inside the $room words of its data records")
        end
        chars = _das_extent(io, raw, UInt8, clusters[1], lastla[1], swap, writable)
        dps = _das_extent(io, raw, Float64, clusters[2], lastla[2], swap, writable)
        ints = _das_extent(io, raw, Int32, clusters[3], lastla[3], swap, writable)
        resv = raw[_DAS_RECLEN+1:_DAS_RECLEN*(1+nresvr)]
        comments = raw[_DAS_RECLEN*(1+nresvr)+1:_DAS_RECLEN*(1+nresvr+ncomr)]
        h = DASHandle(String(path), writable, true, idword, ifname, nresvr, nresvc, ncomr,
            ncomc, resv, comments, chars, dps, ints, lastrec + 1, Tuple(lastrc),
            Tuple(lastwd))
        writable && _das_layout!(h)
        h
    end
end

"Record counts (chars, doubles, integers) of the segregated file of `h`."
_das_nrecs(h::DASHandle) = (cld(length(h.chars), _DAS_NW[1]), cld(length(h.dps), _DAS_NW[2]),
    cld(length(h.ints), _DAS_NW[3]))

"""
    _das_layout!(h) -> h

Set `free`, `lastrc` and `lastwd` of the write handle `h` to the values of the
segregated file that `dascls` writes, as `dassdr` of CSPICE leaves them: one directory
record, in which the first type that has data is at word 10 and the next types follow.
"""
function _das_layout!(h::DASHandle)
    dirrec = h.nresvr + h.ncomr + 2
    lastrc = [0, 0, 0]
    lastwd = [0, 0, 0]
    pos = 10
    for (t, n) in enumerate((length(h.chars), length(h.dps), length(h.ints)))
        n > 0 || continue
        lastrc[t] = dirrec
        lastwd[t] = pos
        pos += 1
    end
    h.lastrc = Tuple(lastrc)
    h.lastwd = Tuple(lastwd)
    h.free = dirrec + sum(_das_nrecs(h)) + 1
    return h
end

# ---------------------------------------------------------------------------------------
# Writing a file image
# ---------------------------------------------------------------------------------------

"The 1024 bytes of the file record of `h`, as `zzdasnfr` of CSPICE writes them."
function _das_file_record(h::DASHandle)
    rec = zeros(UInt8, _DAS_RECLEN)
    fill!(view(rec, 1:68), UInt8(' '))
    idw = codeunits(h.idword)
    copyto!(rec, 1, idw, 1, min(8, length(idw)))
    ifn = codeunits(h.ifname)
    copyto!(rec, 9, ifn, 1, min(60, length(ifn)))
    ints = reinterpret(Int32, view(rec, 69:84))
    ints[1] = h.nresvr
    ints[2] = h.nresvc
    ints[3] = h.ncomr
    ints[4] = h.ncomc
    copyto!(rec, 85, codeunits(_DAS_NATIVE_BFF), 1, 8)
    copyto!(rec, 700, _DAF_FTPSTR, 1, length(_DAF_FTPSTR))
    return rec
end

"The directory record of the segregated file of `h`, as `dassdr` of CSPICE writes it."
function _das_directory_record(h::DASHandle)
    dir = zeros(Int32, 256)
    counts = _das_nrecs(h)
    last = (length(h.chars), length(h.dps), length(h.ints))
    pos = 9
    prvtyp = 0
    for t in 1:3                         # invariant: types before t are described
        last[t] > 0 || continue
        dir[2t+1] = 1
        dir[2t+2] = last[t]
        if pos == 9
            dir[9] = t
            dir[10] = counts[t]
            pos = 11
        else
            dir[pos] = t == _DAS_NEXT[prvtyp] ? counts[t] : -counts[t]
            pos += 1
        end
        prvtyp = t
    end
    return dir
end

"Write `data` to `io` and then `fillbyte` up to the next multiple of 1024 bytes."
function _das_write_padded(io::IO, data::Vector, fillbyte::UInt8)
    write(io, data)
    pad = mod(-sizeof(data), _DAS_RECLEN)
    pad > 0 && write(io, fill(fillbyte, pad))
    return nothing
end

"""
    _das_write_image(h) -> Nothing

Write the complete segregated file of the write handle `h` to a temporary file beside
`h.path` and move it into place. Postcondition: the file holds the file record, the
reserved records, the comment records, one directory record, and the character, double
and integer records; the unused words of a last record are blank (characters) or zero.
"""
function _das_write_image(h::DASHandle)
    tmp = string(h.path, ".", getpid(), ".", objectid(h), ".tmp")
    try
        open(tmp, "w") do io
            write(io, _das_file_record(h))
            write(io, h.resv)
            write(io, h.comments)
            write(io, _das_directory_record(h))
            _das_write_padded(io, h.chars, UInt8(' '))
            _das_write_padded(io, h.dps, 0x00)
            _das_write_padded(io, h.ints, 0x00)
        end
        mv(tmp, h.path; force=true)
    catch e
        rm(tmp; force=true)
        e isa SpiceError && rethrow()
        throw(SpiceError(:daswritefail,
            "cannot write the DAS file $(h.path): $(sprint(showerror, e))"))
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Open and close
# ---------------------------------------------------------------------------------------

"""
    dasopr(path) -> DASHandle

Open the DAS file `path` for read access, as `dasopr_c`. A file of the other byte
order and a file that is not segregated can be read. Contract: `path` names an
existing DAS file; the result is an open read handle. Throws `SpiceError`
(`:blankfilename`, `:filenotfound`, `:filereadfailed`, `:unknownfilarc`,
`:filarchmismatch`, `:ftpxfererror`, `:unknownbff`, `:baddasdirectory`). The handle
belongs to the caller: close it with `dascls`, and do not share it between tasks
without a lock.
"""
dasopr(path::AbstractString) = _das_load(path, "dasopr", false)

"""
    dasopw(path) -> DASHandle

Open the existing DAS file `path` for write access, as `dasopw_c`. The file is read
into memory; all changes are on disk when `dascls` returns. Contract and errors as
[`dasopr`](@ref); in addition the file must be writable (`:fileopenfailed`).
"""
function dasopw(path::AbstractString)
    h = _das_load(path, "dasopw", true)
    (uperm(path) & 0x02) != 0 || throw(SpiceError(:fileopenfailed,
        "dasopw: the file $path cannot be opened for write"))
    return h
end

"""
    dasonw(path, ftype, ifname, ncomr) -> DASHandle

Create the new DAS file `path` of type `ftype` (for example `"DSK"`) with the internal
file name `ifname` and `ncomr` comment records, and open it for write, as `dasonw` of
SPICELIB.

Preconditions (each throws `SpiceError`): `path` is not blank (`:blankfilename`) and
does not exist (`:fileopenfailed`); `ftype` is not blank (`:blankfiletype`) and
printable (`:illegalcharacter`); `ncomr >= 0` (`:invalidcount`). Postcondition: the
file exists and holds the file record, `ncomr` comment records and an empty directory
record; the handle has no data.
"""
function dasonw(path::AbstractString, ftype::AbstractString, ifname::AbstractString,
        ncomr::Integer)
    isempty(strip(path)) && throw(SpiceError(:blankfilename, "dasonw: the file name is blank"))
    isempty(strip(ftype, ' ')) && throw(SpiceError(:blankfiletype, "dasonw: the file type is blank"))
    ttype = String(strip(ftype, ' '))
    all(c -> 0x20 <= c <= 0x7e, codeunits(ttype)) || throw(SpiceError(:illegalcharacter,
        "dasonw: the file type contains nonprinting characters"))
    ncomr < 0 && throw(SpiceError(:invalidcount,
        "dasonw: the number of comment records allocated must be non-negative but was $ncomr"))
    ispath(path) && throw(SpiceError(:fileopenfailed,
        "dasonw: the file $path already exists; a new DAS file cannot replace it"))
    idword = rpad("DAS/" * ttype[1:min(4, lastindex(ttype))], 8)
    h = DASHandle(String(path), true, true, idword, String(rstrip(ifname)), 0, 0, Int(ncomr),
        0, UInt8[], zeros(UInt8, Int(ncomr) * _DAS_RECLEN), UInt8[], Float64[], Int32[],
        0, (0, 0, 0), (0, 0, 0))
    _das_layout!(h)
    try
        open(path, "w") do io
            write(io, _das_file_record(h))
            write(io, h.comments)
            write(io, zeros(Int32, 256))
        end
    catch e
        throw(SpiceError(:fileopenfailed,
            "dasonw: cannot create $path: $(sprint(showerror, e))"))
    end
    return h
end

"""
    dascls(handle) -> Nothing

Close the DAS file of `handle`, as `dascls_c`. For a write handle the complete file is
written in segregated form and moved into place; all data are on disk when the call
returns. Closing a handle that is already closed does nothing.
"""
function dascls(h::DASHandle)
    h.isopen || return nothing
    if h.writable
        _das_layout!(h)
        _das_write_image(h)
    end
    h.isopen = false
    return nothing
end

"""
    dasllc(handle) -> Nothing

Close the DAS file of `handle`, as `dasllc` of SPICELIB. CSPICE does not segregate the
file in this call; here the data of a write handle exist only in memory before the
close, so the call does what [`dascls`](@ref) does and the file is segregated.
"""
dasllc(h::DASHandle) = dascls(h)

# ---------------------------------------------------------------------------------------
# File information
# ---------------------------------------------------------------------------------------

"""
    dashfn(handle) -> String

Name of the DAS file of `handle`, as `dashfn_c`. Throws `SpiceError(:dasnosuchhandle)`
for a closed handle.
"""
function dashfn(h::DASHandle)
    _das_check_open(h, "dashfn")
    return h.path
end

"""
    dasrfr(handle) -> (idword, ifname, nresvr, nresvc, ncomr, ncomc)

Contents of the file record, as `dasrfr_c`: the ID word and the internal file name
(trailing blanks removed), the numbers of reserved records and characters, and the
numbers of comment records and characters.
"""
function dasrfr(h::DASHandle)
    _das_check_open(h, "dasrfr")
    return String(rstrip(h.idword)), h.ifname, h.nresvr, h.nresvc, h.ncomr, h.ncomc
end

"""
    daslla(handle) -> (lastc, lastd, lasti)

Last logical addresses in use of the character, double and integer address spaces, as
`daslla` of SPICELIB.
"""
function daslla(h::DASHandle)
    _das_check_open(h, "daslla")
    return length(h.chars), length(h.dps), length(h.ints)
end

"""
    dashfs(handle) -> (nresvr, nresvc, ncomr, ncomc, free, lastla, lastrc, lastwd)

File summary, as `dashfs` of SPICELIB: the counts of the file record, the first free
record number, and for the three data types (character, double, integer) the last
logical address, the record number of the last directory record that describes the
type, and the word of that record that holds its last cluster descriptor. For a write
handle the last three describe the segregated file that `dascls` will write.
"""
function dashfs(h::DASHandle)
    _das_check_open(h, "dashfs")
    h.writable && _das_layout!(h)
    return h.nresvr, h.nresvc, h.ncomr, h.ncomc, h.free, daslla(h), h.lastrc, h.lastwd
end

# ---------------------------------------------------------------------------------------
# Data: read, add, update
# ---------------------------------------------------------------------------------------

"Throws `SpiceError(:dasnosuchaddress)` unless `first:last` lies in `1:n`; requires `first <= last`."
@inline function _das_check_range(h::DASHandle, who::AbstractString, first::Integer,
        last::Integer, n::Int)
    (1 <= first && last <= n) || throw(SpiceError(:dasnosuchaddress,
        "$who: the address range $first:$last is not inside 1:$n in $(h.path)"))
    return nothing
end

"""
    dasrdd(handle, first, last) -> Vector{Float64}

The doubles at the DAS addresses `first:last`, as `dasrdd` of SPICELIB. The result is
empty if `first > last`. Throws `SpiceError(:dasnosuchaddress)` for an address outside
the file.
"""
function dasrdd(h::DASHandle, first::Integer, last::Integer)
    _das_check_open(h, "dasrdd")
    first > last && return Float64[]
    _das_check_range(h, "dasrdd", first, last, length(h.dps))
    return h.dps[first:last]
end

"""
    dasrdi(handle, first, last) -> Vector{Int}

The integers at the DAS addresses `first:last`, as `dasrdi` of SPICELIB. The result is
empty if `first > last`. Throws `SpiceError(:dasnosuchaddress)` for an address outside
the file.
"""
function dasrdi(h::DASHandle, first::Integer, last::Integer)
    _das_check_open(h, "dasrdi")
    first > last && return Int[]
    _das_check_range(h, "dasrdi", first, last, length(h.ints))
    return Int[h.ints[a] for a in first:last]
end

"""
    dasrdc(handle, first, last) -> String

The characters at the DAS addresses `first:last` as one string, as `dasrdc` of SPICELIB
with one substring. The result is empty if `first > last`. Throws
`SpiceError(:dasnosuchaddress)` for an address outside the file.
"""
function dasrdc(h::DASHandle, first::Integer, last::Integer)
    _das_check_open(h, "dasrdc")
    first > last && return ""
    _das_check_range(h, "dasrdc", first, last, length(h.chars))
    return String(h.chars[first:last])
end

"""
    dasadd(handle, data) -> Nothing

Append the doubles `data` to the double address space, as `dasadd` of SPICELIB.
Precondition: write access (`SpiceError(:dasinvalidaccess)`).
"""
function dasadd(h::DASHandle, data::AbstractVector{<:Real})
    _das_check_write(h, "dasadd")
    append!(h.dps, data)
    return nothing
end

"""
    dasadi(handle, data) -> Nothing

Append the integers `data` to the integer address space, as `dasadi` of SPICELIB.
Preconditions: write access (`SpiceError(:dasinvalidaccess)`); every value fits in 32
bits (`InexactError`).
"""
function dasadi(h::DASHandle, data::AbstractVector{<:Integer})
    _das_check_write(h, "dasadi")
    append!(h.ints, data)
    return nothing
end

"""
    dasadc(handle, data) -> Nothing

Append the characters of the string `data` to the character address space, as `dasadc`
of SPICELIB with one substring. Precondition: write access.
"""
function dasadc(h::DASHandle, data::AbstractString)
    _das_check_write(h, "dasadc")
    append!(h.chars, codeunits(data))
    return nothing
end

"""
    dasudd(handle, first, last, data) -> Nothing

Replace the doubles at the addresses `first:last` with `data`, as `dasudd` of SPICELIB.
Nothing happens if `first > last`. Preconditions: write access; the range is inside the
file (`SpiceError(:dasnosuchaddress)`); `length(data) >= last - first + 1`.
"""
function dasudd(h::DASHandle, first::Integer, last::Integer, data::AbstractVector{<:Real})
    _das_check_write(h, "dasudd")
    first > last && return nothing
    _das_check_range(h, "dasudd", first, last, length(h.dps))
    length(data) >= last - first + 1 || throw(ArgumentError("dasudd: data is too short"))
    copyto!(h.dps, first, data, firstindex(data), last - first + 1)
    return nothing
end

"""
    dasudi(handle, first, last, data) -> Nothing

Replace the integers at the addresses `first:last` with `data`, as `dasudi` of
SPICELIB. Nothing happens if `first > last`. Preconditions as [`dasudd`](@ref).
"""
function dasudi(h::DASHandle, first::Integer, last::Integer, data::AbstractVector{<:Integer})
    _das_check_write(h, "dasudi")
    first > last && return nothing
    _das_check_range(h, "dasudi", first, last, length(h.ints))
    length(data) >= last - first + 1 || throw(ArgumentError("dasudi: data is too short"))
    copyto!(h.ints, first, data, firstindex(data), last - first + 1)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Comment area
# ---------------------------------------------------------------------------------------

"""
    dasac(handle, lines) -> Nothing

Append the comment lines `lines` to the comment area of the DAS file, as `dasac_c`.
Trailing blanks of a line are not stored. Comment records are added if the area does
not have room.

Preconditions (each throws `SpiceError`): write access (`:dasinvalidaccess`); at least
one line (`:invalidargument`); only printable ASCII characters, codes 32 to 126
(`:illegalcharacter`). Postcondition: `dasec(handle)` gives the old lines and then
`lines` without trailing blanks.
"""
function dasac(h::DASHandle, lines::AbstractVector{<:AbstractString})
    _das_check_write(h, "dasac")
    n = length(lines)
    n <= 0 && throw(SpiceError(:invalidargument, "dasac: the number of comment lines " *
        "to be added to the binary DAS file $(h.path) was not positive: $n"))
    texts = Vector{Vector{UInt8}}(undef, n)
    nchars = 0
    for (i, line) in enumerate(lines)
        text = Vector{UInt8}(codeunits(rstrip(line, ' ')))
        for (j, c) in enumerate(text)
            (0x20 <= c <= 0x7e) || throw(SpiceError(:illegalcharacter, "dasac: a " *
                "nonprinting character was encountered in the comment buffer. Value: " *
                "$(Int(c)), line: $i, position: $j"))
        end
        texts[i] = text
        nchars += length(text) + 1
    end
    space = h.ncomr * _DAS_RECLEN - h.ncomc
    if nchars > space
        newrec = (nchars - space - 1) ÷ _DAS_RECLEN + 1
        append!(h.comments, zeros(UInt8, newrec * _DAS_RECLEN))
        h.ncomr += newrec
    end
    pos = h.ncomc                        # characters of the area that are in use
    # A record that the text enters at its first character is blank before the text
    # goes in, as the record buffer of `dasac` is.
    blank_from(p) = fill!(view(h.comments, p+1:p+_DAS_RECLEN), UInt8(' '))
    pos == 0 && h.ncomr > 0 && blank_from(0)
    for text in texts                    # invariant: `pos` characters are in the area
        for c in text
            pos % _DAS_RECLEN == 0 && pos > 0 && blank_from(pos)
            h.comments[pos+1] = c
            pos += 1
        end
        pos % _DAS_RECLEN == 0 && pos > 0 && blank_from(pos)
        h.comments[pos+1] = 0x00
        pos += 1
    end
    h.ncomc = pos
    return nothing
end

"""
    dasec(handle) -> Vector{String}

All lines of the comment area of the DAS file, in order, as repeated calls of
`dasec_c` give them (trailing blanks are not stored). The result is empty if the file
has no comments. Throws `SpiceError(:baddascommentarea)` if the last line is not
terminated inside the comment characters that the file record states.
"""
function dasec(h::DASHandle)
    _das_check_open(h, "dasec")
    lines = String[]
    h.ncomc == 0 && return lines
    h.comments[h.ncomc] == 0x00 || throw(SpiceError(:baddascommentarea,
        "dasec: the last comment line of $(h.path) is not terminated; the count of " *
        "comment characters is $(h.ncomc)"))
    start = 1
    for k in 1:h.ncomc                   # invariant: lines hold the text before `start`
        h.comments[k] == 0x00 || continue
        push!(lines, String(h.comments[start:k-1]))
        start = k + 1
    end
    return lines
end

"""
    dasdc(handle) -> Nothing

Delete the whole comment area of the DAS file, as `dasdc_c`: the comment records are
removed and the count of comment characters becomes zero. Precondition: write access.
"""
function dasdc(h::DASHandle)
    _das_check_write(h, "dasdc")
    h.comments = UInt8[]
    h.ncomr = 0
    h.ncomc = 0
    return nothing
end
