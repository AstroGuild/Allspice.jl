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
# keeps a consistent image. The only global state is the table of open files of
# `daf_access.jl`. One handle is not safe for concurrent use: a caller that shares a
# handle between tasks must hold a lock.

const _DAS_RECLEN = 1024
const _DAS_NW = (1024, 128, 256)          # words in a record of type 1 (char), 2 (dp), 3 (int)
const _DAS_NEXT = (2, 3, 1)
const _DAS_PREV = (3, 1, 2)
const _DAS_NATIVE_LITTLE = ENDIAN_BOM == 0x04030201
const _DAS_NATIVE_BFF = _DAS_NATIVE_LITTLE ? "LTL-IEEE" : "BIG-IEEE"

"""
    DASHandle

An open DAS file. It replaces the integer handle of CSPICE.

[`dasopr`](@ref), [`dasopw`](@ref), [`dasonw`](@ref), [`dasops`](@ref), [`dlaopn`](@ref),
[`dskopn`](@ref), [`ekopr`](@ref), [`ekopw`](@ref), [`ekopn`](@ref) and [`ekops`](@ref)
return a handle. Close it with [`dascls`](@ref), [`dasllc`](@ref), [`dskcls`](@ref),
[`ekcls`](@ref) or `close`, which all do the same; `isopen(handle)` tells whether it is
open.

The handle holds the three address spaces of the file in memory: characters, doubles and
32-bit integers, each numbered from 1. A handle that is open for writing changes this
memory only, and the file on disk changes when `dascls` runs. The data of a handle that is
open for reading are never modified.

A DAS file is segregated if all records of one data type are together.

Two handles are independent. One handle is not safe for concurrent use: tasks that share a
handle must hold a lock around every call. A file that an open handle has open, for read or
for write, cannot be opened for write by another handle (`SpiceError(:fileopenconflict)`),
as in CSPICE.
Unlike CSPICE, a file that is open for write can be opened for read, and a kernel that
[`furnsh`](@ref) loaded can be opened for write.

A call on a closed handle throws `SpiceError(:dasnosuchhandle)`; the close calls do
nothing instead. A write call on a handle that is open for reading throws
`SpiceError(:dasinvalidaccess)`. A call that throws changes neither the data of the handle
nor the file.
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

@inline function _das_check_open(h::DASHandle, who::AbstractString)
    h.isopen || throw(SpiceError(:dasnosuchhandle, "$who: the DAS file $(h.path) is closed"))
    return nothing
end

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
`:emptystring`, `:blankfilename`, `:filenotfound`, `:filereadfailed`, `:unknownfilarc`,
`:filarchmismatch`, `:ftpxfererror`, `:unknownbff` or `:unsupportedbff`. A file from
before the binary-file-format string (bytes 85 to 92 zero) is taken as native.
"""
function _das_probe(path::AbstractString, who::AbstractString)
    _chk_filename(who, path)
    isfile(path) || throw(SpiceError(:filenotfound, "$who: the file $path does not exist"))
    raw = open(io -> read(io, _DAS_RECLEN), path, "r")
    length(raw) == _DAS_RECLEN || throw(SpiceError(:filereadfailed,
        "$who: $path is too short ($(length(raw)) bytes) to hold a DAS file record"))
    arch, _ = idw2at(String(raw[1:8]))
    arch == "DAF" && throw(SpiceError(:filarchmismatch,
        "$who: the file $path is a DAF, not a DAS file"))
    arch == "DAS" || throw(SpiceError(:unknownfilarc,
        "$who: the file $path is neither a DAF nor a DAS file"))
    if raw[700:705] == _DAF_FTPSTR[1:6] && raw[700:727] != _DAF_FTPSTR
        throw(SpiceError(:ftpxfererror, "$who: the FTP validation string of $path is " *
            "damaged; the file was most likely moved in text mode"))
    end
    bff = String(raw[85:92])
    little = if bff == "LTL-IEEE"
        true
    elseif bff == "BIG-IEEE"
        false
    elseif bff in ("VAX-GFLT", "VAX-DFLT")
        throw(SpiceError(:unsupportedbff,
            "$who: the binary file format $bff of $path is not supported"))
    elseif all(iszero, raw[85:92])
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
allows it). Throws the errors of `_das_probe`, `SpiceError(:dasfilereadfailed)` for a
directory record past the end of the file, and `SpiceError(:baddasdirectory)` for a
directory that points to data records outside the file, has a cycle, or does not
cover the last addresses that it states.
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
        try
            (nresvr >= 0 && ncomr >= 0) ||
                bad("the file record gives $nresvr reserved and $ncomr comment records")
            nresvr + ncomr + 2 <= nrecords || throw(SpiceError(:dasfilereadfailed,
                "$who: the first directory record $(nresvr + ncomr + 2) of $path is past " *
                "the end of the file ($nrecords records)"))
            (0 <= ncomc <= ncomr * _DAS_RECLEN) ||
                bad("the comment area of $ncomr records cannot hold $ncomc characters")
            clusters = (NTuple{2,Int}[], NTuple{2,Int}[], NTuple{2,Int}[])
            lastla = [0, 0, 0]
            lastrc = [0, 0, 0]
            lastwd = [0, 0, 0]
            lastrec = 0                  # the last directory or data record in use
            rec = nresvr + ncomr + 2
            visited = 0
            while rec > 0
                # invariant: `visited` directory records were read, none twice while
                # `visited <= nrecords`; `clusters` holds their data records in address
                # order
                visited += 1
                visited <= nrecords || bad("the chain of directory records has a cycle")
                rec <= nrecords || throw(SpiceError(:dasfilereadfailed,
                    "$who: directory record $rec of $path is past the end of the file"))
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
                            bad("directory record $rec describes records past the end " *
                                "of the file")
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
                lastla[t] <= room || bad("the last address $(lastla[t]) of type $t is " *
                    "not inside the $room words of its data records")
            end
            chars = _das_extent(io, raw, UInt8, clusters[1], lastla[1], swap, writable)
            dps = _das_extent(io, raw, Float64, clusters[2], lastla[2], swap, writable)
            ints = _das_extent(io, raw, Int32, clusters[3], lastla[3], swap, writable)
            resv = raw[_DAS_RECLEN+1:_DAS_RECLEN*(1+nresvr)]
            comments = raw[_DAS_RECLEN*(1+nresvr)+1:_DAS_RECLEN*(1+nresvr+ncomr)]
            h = DASHandle(String(path), writable, true, idword, ifname, nresvr, nresvc,
                ncomr, ncomc, resv, comments, chars, dps, ints, lastrec + 1, Tuple(lastrc),
                Tuple(lastwd))
            writable && _das_layout!(h)
            h
        finally
            _unmap!(raw)                 # `h` has copies or maps of its own; see `_unmap!`
        end
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

Open the DAS file `path` for reading and return its handle, as `dasopr_c`.

The file can have either byte order, and it need not be segregated (see
[`DASHandle`](@ref)). Close the handle with [`dascls`](@ref).

Throws `SpiceError`:
- `:emptystring` if `path` is empty, `:blankfilename` if it holds only blanks;
- `:filenotfound` if the file does not exist;
- `:filereadfailed` if the file is shorter than one record (1024 bytes);
- `:filarchmismatch` if the file is a DAF, `:unknownfilarc` if it is neither a DAF
  nor a DAS file (a text kernel, for example);
- `:ftpxfererror` if the file was damaged by a transfer in text mode;
- `:unknownbff` if the binary file format (the number format that the file record
  names) is not known, `:unsupportedbff` if it is `VAX-GFLT` or `VAX-DFLT`;
- `:dasfilereadfailed` if a directory record is past the end of the file;
- `:baddasdirectory` if the directory records do not agree with the file.

See also [`dasopw`](@ref), [`dlabfs`](@ref), [`dasrdd`](@ref).
"""
dasopr(path::AbstractString) = _open_register!(_das_load(path, "dasopr", false), path, false, "dasopr")

"""
    dasopw(path) -> DASHandle

Open the existing DAS file `path` for writing, as `dasopw_c`.

The file is read into memory. The file on disk does not change until [`dascls`](@ref),
which writes the complete file.

No open handle of this process may have the file open, for read or for write, under any
path. After the call, [`dasopr`](@ref) can still open the file for read, where CSPICE
refuses; that handle keeps the data of the file as they were when it opened.

Throws the `SpiceError` kinds of [`dasopr`](@ref), `:fileopenconflict` if another handle
has the file open, and `:fileopenfailed` if the file is not writable.
"""
function dasopw(path::AbstractString)
    h = _das_load(path, "dasopw", true)
    (uperm(path) & 0x02) != 0 || throw(SpiceError(:fileopenfailed,
        "dasopw: the file $path cannot be opened for write"))
    return _open_register!(h, path, true, "dasopw")
end

"""
    dasonw(path, ftype, ifname, ncomr) -> DASHandle

Create the DAS file `path` and open it for writing, as `dasonw` of SPICELIB.

`ftype` is the file type (for example `"DSK"`; four characters are stored), `ifname` is
the internal file name and `ncomr` is the number of comment records, 1024 characters in
each.

The file is on disk when the call returns, with the file record, `ncomr` comment records
and an empty directory record. Data that are added afterwards are in memory until
[`dascls`](@ref) writes them. An existing file is never replaced.

Throws `SpiceError`:
- `:blankfilename` if `path` is blank;
- `:fileopenfailed` if `path` exists or the file cannot be created;
- `:blankfiletype` if `ftype` is blank, `:illegalcharacter` if it has a character that is
  not printable;
- `:invalidcount` if `ncomr < 0`.

See also [`dlaopn`](@ref), [`dskopn`](@ref), [`dasops`](@ref).
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
    return _open_register!(h, path, true, "dasonw"; check=false)
end

"""
    dasops() -> DASHandle

Open a scratch DAS file for writing, as `dasops` of SPICELIB: a DAS file that exists only
while its handle is open.

CSPICE makes a temporary file and deletes it at the close. Here the data are in memory, so
no file exists at any time and [`dascls`](@ref) writes nothing. The handle has no data, the
file name `""` and the file type `SCR`.

# Examples
```jldoctest
julia> h = dasops();

julia> dasadd(h, [1.0, 2.0, 3.0])

julia> dasadi(h, [7, 8])

julia> daslla(h)
(0, 3, 2)

julia> dasudd(h, 2, 2, [20.0])

julia> dasrdd(h, 1, 3)
3-element Vector{Float64}:
  1.0
 20.0
  3.0

julia> dascls(h)
```
"""
function dasops()
    h = DASHandle("", true, true, rpad("DAS/SCR", 8), "", 0, 0, 0, 0, UInt8[], UInt8[], UInt8[],
        Float64[], Int32[], 0, (0, 0, 0), (0, 0, 0))
    _das_layout!(h)
    return h
end

"""
    dascls(handle) -> Nothing

Close the DAS file of `handle`, as `dascls_c`.

A handle that is open for writing writes the complete file, in segregated form (see
[`DASHandle`](@ref)). The file goes to a temporary file beside it and is moved into place
in one step. All data are on disk when the call returns. A scratch file
([`dasops`](@ref)) is not written. The call does nothing if the handle is closed already.
After the call the handle keeps no memory map of the file, so that Windows can delete or
replace it.

Throws `SpiceError(:daswritefail)` if the file cannot be written; the handle then stays
open. On Windows the file cannot be replaced, and the call throws, while another handle
has it open for reading or a context has it loaded.
"""
function dascls(h::DASHandle)
    h.isopen || return nothing
    if h.writable && !isempty(h.path)       # a scratch file (`dasops`) has no path
        _das_layout!(h)
        _das_write_image(h)
    end
    h.isopen = false
    _open_unregister!(h)
    _das_release!(h)
    return nothing
end

"""
    _das_release!(h) -> h

Remove the memory maps of the address spaces of `h` and put empty arrays in their place,
so that a closed handle does not keep its file mapped (see `_unmap!`). No result of a
call on `h` refers to these arrays, because the read calls copy.
"""
function _das_release!(h::DASHandle)
    _unmap!(h.chars)
    _unmap!(h.dps)
    _unmap!(h.ints)
    h.chars = UInt8[]
    h.dps = Float64[]
    h.ints = Int32[]
    return h
end

"""
    dasllc(handle) -> Nothing

Close the DAS file of `handle`, as `dasllc` of SPICELIB.

The call does what [`dascls`](@ref) does. CSPICE does not segregate the file in `dasllc`;
here the data of a write handle are only in memory before the close, so the file is
written and it is segregated.
"""
dasllc(h::DASHandle) = dascls(h)

# ---------------------------------------------------------------------------------------
# File information
# ---------------------------------------------------------------------------------------

"""
    dashfn(handle) -> String

Return the name of the DAS file of `handle`, as `dashfn_c`. A scratch file
([`dasops`](@ref)) has the name `""`.
"""
function dashfn(h::DASHandle)
    _das_check_open(h, "dashfn")
    return h.path
end

"""
    dasrfr(handle) -> (idword, ifname, nresvr, nresvc, ncomr, ncomc)

Return the contents of the file record, as `dasrfr_c`: the ID word and the internal file
name (both without trailing blanks), the numbers of reserved records and of reserved
characters, and the numbers of comment records and of comment characters.

For a handle that is open for writing, the comment counts include the changes that are not
on disk yet.
"""
function dasrfr(h::DASHandle)
    _das_check_open(h, "dasrfr")
    return String(rstrip(h.idword)), h.ifname, h.nresvr, h.nresvc, h.ncomr, h.ncomc
end

"""
    daslla(handle) -> (lastc, lastd, lasti)

Return the last addresses in use of the character, the double and the integer address
space, as `daslla` of SPICELIB. An address space without data gives 0.
"""
function daslla(h::DASHandle)
    _das_check_open(h, "daslla")
    return length(h.chars), length(h.dps), length(h.ints)
end

"""
    dashfs(handle) -> (nresvr, nresvc, ncomr, ncomc, free, lastla, lastrc, lastwd)

Return the file summary, as `dashfs` of SPICELIB.

- `nresvr`, `nresvc`: the numbers of reserved records and of reserved characters.
- `ncomr`, `ncomc`: the numbers of comment records and of comment characters.
- `free`: the first free record number.
- `lastla`, `lastrc`, `lastwd`: tuples of three values, for the character, the double
  and the integer data.
  - `lastla`: the last address in use.
  - `lastrc`: the number of the last directory record that describes the data type.
  - `lastwd`: the word of that record that describes the last run of records of the
    type.

For a handle that is open for writing, `free`, `lastrc` and `lastwd` describe the
segregated file that [`dascls`](@ref) will write.
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

# The update routines of SPICELIB check both addresses before they look at the order,
# so an empty range outside the address space is an error there.
@inline function _das_check_update(h::DASHandle, who::AbstractString, first::Integer,
        last::Integer, n::Int)
    (1 <= first <= n && 1 <= last <= n) || throw(SpiceError(:invalidaddress,
        "$who: first was $first and last was $last; the valid range is 1:$n in $(h.path)"))
    return nothing
end

"""
    dasrdd(handle, first, last) -> Vector{Float64}

Return the doubles at the DAS addresses `first:last`, as `dasrdd` of SPICELIB. Addresses
start at 1. The result is empty if `first > last`.

Throws `SpiceError(:dasnosuchaddress)` if the range is not inside the address space.

See also [`dasrdi`](@ref), [`dasrdc`](@ref), [`daslla`](@ref).
"""
function dasrdd(h::DASHandle, first::Integer, last::Integer)
    _das_check_open(h, "dasrdd")
    first > last && return Float64[]
    _das_check_range(h, "dasrdd", first, last, length(h.dps))
    return h.dps[first:last]
end

"""
    dasrdi(handle, first, last) -> Vector{Int}

Return the integers at the DAS addresses `first:last`, as `dasrdi` of SPICELIB. Addresses
start at 1. The result is empty if `first > last`.

Throws `SpiceError(:dasnosuchaddress)` if the range is not inside the address space.
"""
function dasrdi(h::DASHandle, first::Integer, last::Integer)
    _das_check_open(h, "dasrdi")
    first > last && return Int[]
    _das_check_range(h, "dasrdi", first, last, length(h.ints))
    return Int[h.ints[a] for a in first:last]
end

"""
    dasrdc(handle, first, last) -> String

Return the characters at the DAS addresses `first:last` as one string. `dasrdc` of
SPICELIB fills substrings of an array of strings. Addresses start at 1. The result is
empty if `first > last`.

Throws `SpiceError(:dasnosuchaddress)` if the range is not inside the address space.
"""
function dasrdc(h::DASHandle, first::Integer, last::Integer)
    _das_check_open(h, "dasrdc")
    first > last && return ""
    _das_check_range(h, "dasrdc", first, last, length(h.chars))
    return String(h.chars[first:last])
end

"""
    dasadd(handle, data) -> Nothing

Append the doubles `data` to the double address space, as `dasadd` of SPICELIB. The file
on disk changes when [`dascls`](@ref) runs.

See also [`dasadi`](@ref), [`dasadc`](@ref), [`dasudd`](@ref), and [`dasops`](@ref) for
an example.
"""
function dasadd(h::DASHandle, data::AbstractVector{<:Real})
    _das_check_write(h, "dasadd")
    append!(h.dps, data)
    return nothing
end

"Throws `SpiceError(:valueoutofrange)` unless every element of `data` fits in an `Int32`, the integer of a DAS file."
function _das_check_int32(who::AbstractString, data)
    for v in data
        typemin(Int32) <= v <= typemax(Int32) || throw(SpiceError(:valueoutofrange,
            "$who: the value $v does not fit in the 32-bit integer of a DAS file"))
    end
    return nothing
end

"""
    dasadi(handle, data) -> Nothing

Append the integers `data` to the integer address space, as `dasadi` of SPICELIB. The file
on disk changes when [`dascls`](@ref) runs.

Throws `SpiceError(:valueoutofrange)` if a value does not fit in 32 bits.
"""
function dasadi(h::DASHandle, data::AbstractVector{<:Integer})
    _das_check_write(h, "dasadi")
    _das_check_int32("dasadi", data)
    append!(h.ints, data)
    return nothing
end

"""
    dasadc(handle, data) -> Nothing

Append the code units of the string `data` to the character address space. `dasadc` of
SPICELIB takes substrings of an array of strings. The file on disk changes when
[`dascls`](@ref) runs.
"""
function dasadc(h::DASHandle, data::AbstractString)
    _das_check_write(h, "dasadc")
    append!(h.chars, codeunits(data))
    return nothing
end

"""
    dasudd(handle, first, last, data) -> Nothing

Replace the doubles at the DAS addresses `first:last` with the first `last - first + 1`
elements of `data`, as `dasudd` of SPICELIB.

Throws `SpiceError(:invalidaddress)` unless `first` and `last` are both addresses of the
double address space, and `SpiceError(:badarg)` if `data` is too short. If they are
addresses and `first > last`, nothing changes.
"""
function dasudd(h::DASHandle, first::Integer, last::Integer, data::AbstractVector{<:Real})
    _das_check_write(h, "dasudd")
    _das_check_update(h, "dasudd", first, last, length(h.dps))
    first > last && return nothing
    length(data) >= last - first + 1 || throw(SpiceError(:badarg, "dasudd: data is too short"))
    copyto!(h.dps, first, data, firstindex(data), last - first + 1)
    return nothing
end

"""
    dasudi(handle, first, last, data) -> Nothing

Replace the integers at the DAS addresses `first:last` with the first `last - first + 1`
elements of `data`, as `dasudi` of SPICELIB.

Throws `SpiceError(:invalidaddress)` unless `first` and `last` are both addresses of the
integer address space, `SpiceError(:badarg)` if `data` is too short, and
`SpiceError(:valueoutofrange)` if a value does not fit in 32 bits. If they are addresses
and `first > last`, nothing changes.
"""
function dasudi(h::DASHandle, first::Integer, last::Integer, data::AbstractVector{<:Integer})
    _das_check_write(h, "dasudi")
    _das_check_update(h, "dasudi", first, last, length(h.ints))
    first > last && return nothing
    length(data) >= last - first + 1 || throw(SpiceError(:badarg, "dasudi: data is too short"))
    _das_check_int32("dasudi", view(data, firstindex(data):firstindex(data)+(last-first)))
    copyto!(h.ints, first, data, firstindex(data), last - first + 1)
    return nothing
end

"""
    dasudc(handle, first, last, data) -> Nothing

Replace the characters at the DAS addresses `first:last` with the first
`last - first + 1` code units of the string `data`, as `dasudc_c`. `dasudc_c` takes
substrings of an array of strings.

Throws `SpiceError(:invalidaddress)` unless `first` and `last` are both addresses of the
character address space, and `SpiceError(:badarg)` if `data` is too short. If they are
addresses and `first > last`, nothing changes.

See also [`dasadc`](@ref), [`dasrdc`](@ref).
"""
function dasudc(h::DASHandle, first::Integer, last::Integer, data::AbstractString)
    _das_check_write(h, "dasudc")
    _das_check_update(h, "dasudc", first, last, length(h.chars))
    first > last && return nothing
    units = codeunits(data)
    length(units) >= last - first + 1 || throw(SpiceError(:badarg, "dasudc: data is too short"))
    copyto!(h.chars, first, units, 1, last - first + 1)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Comment area
# ---------------------------------------------------------------------------------------

"""
    dasac(handle, lines) -> Nothing

Append the comment lines `lines` to the comment area of the DAS file, as `dasac_c`.

Trailing blanks of a line are not stored. Comment records are added if the area does not
have room. The file on disk changes when [`dascls`](@ref) runs.

Throws `SpiceError(:invalidargument)` if `lines` is empty, and `:illegalcharacter` if a
line has a character outside the printable ASCII codes 32 to 126.

See also [`dasec`](@ref), [`dasdc`](@ref).
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

Return all lines of the comment area of the DAS file, in order.

One call returns the whole comment area, where `dasec_c` returns one buffer of lines on
each call. The result is empty if the file has no comments.

Throws `SpiceError(:baddascommentarea)` if the comment area is damaged.
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

Delete the comment area of the DAS file, as `dasdc_c`: the file then has no comment
record and no comment character. The file on disk changes when [`dascls`](@ref) runs.
"""
function dasdc(h::DASHandle)
    _das_check_write(h, "dasdc")
    h.comments = UInt8[]
    h.ncomr = 0
    h.ncomc = 0
    return nothing
end
