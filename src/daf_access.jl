# Low-level access to DAF files: open for read, the array search, summaries, names,
# data, the comment area, the file record, and the file architecture (`getfat`).
#
# CSPICE keeps a process-wide table of open DAF files and one hidden "current array"
# per search. Here the opened file is an object, `DAFHandle`, that the caller owns, and
# the search state lives inside it. The only global state is the table of open files
# below, under its own lock, so two tasks can work on two handles at the same time. One
# handle is not safe for concurrent use: a caller that shares a handle between tasks
# must hold a lock around every call.
#
# Rule for kernels that are written: `furnsh` a file only after `dafcls`, `spkcls` or
# `ckcls` has returned. `furnsh` memory-maps the file and takes its segment
# descriptors once; the writers never shorten or truncate a file (a new file must not
# exist, and an existing file only grows), so a map that is already in a `KernelSet`
# stays valid, but it does not see segments that are added later. Unload a kernel
# before appending to it, and load it again after the close.

const _DAF_CRECLEN = 1000                 # characters of a comment record or name record
const _DAF_NATIVE_LITTLE = ENDIAN_BOM == 0x04030201
const _DAF_NATIVE_BFF = _DAF_NATIVE_LITTLE ? "LTL-IEEE" : "BIG-IEEE"
const _DAF_EOC = 0x04                     # end of the comment area
const _DAF_EOL = 0x00                     # end of a comment line

# The table of open files, shared by `DAFHandle` and `DASHandle`. It maps a handle that a
# public open call returned to the device and inode of its file and its access. As in
# the file table of CSPICE (`zzddhopn`), a file that a handle has open, for read or for
# write, cannot be opened for write again; unlike CSPICE, a file that is open for write
# can be opened for read. The keys are weak, so a handle that is no longer referenced
# leaves the table when the garbage collector frees it.
const _OPEN_FILES = WeakKeyDict{Any,Tuple{UInt64,UInt64,Bool}}()
const _OPEN_FILES_LOCK = ReentrantLock()

function _file_id(path::AbstractString)
    st = stat(path)
    return UInt64(st.device), UInt64(st.inode)
end

"""
    _open_register!(h, path, writable, who; check=writable) -> h

Enter the open handle `h` of the existing file `path` in the table of open files. With
`check`, throw `SpiceError(:fileopenconflict)` if another open handle has the same file
(the same device and inode, under any path); `h` is then not entered.
"""
function _open_register!(h, path::AbstractString, writable::Bool, who::AbstractString;
        check::Bool=writable)
    @assert isopen(h)
    dev, ino = _file_id(path)
    lock(_OPEN_FILES_LOCK) do
        if check
            for (other, (odev, oino, ow)) in _OPEN_FILES
                # invariant: no handle before `other` has the file open
                (odev == dev && oino == ino && other !== h && isopen(other)) &&
                    throw(SpiceError(:fileopenconflict, "$who: the file $path is open " *
                        "for $(ow ? "write" : "read") through another handle; close it first"))
            end
        end
        _OPEN_FILES[h] = (dev, ino, writable)
    end
    return h
end

_open_unregister!(h) = lock(() -> delete!(_OPEN_FILES, h), _OPEN_FILES_LOCK)

"""
    DAFHandle

An open DAF file. It replaces the integer handle of CSPICE.

[`dafopr`](@ref), [`dafopw`](@ref), [`dafonw`](@ref), [`spkopn`](@ref), [`spkopa`](@ref),
[`ckopn`](@ref) and [`pckopn`](@ref) return a handle. Close it with [`dafcls`](@ref),
[`spkcls`](@ref), [`ckcls`](@ref), [`pckcls`](@ref) or `close`; `isopen(handle)` tells
whether it is open.

The handle holds the state that CSPICE keeps in global tables: the array search
(`dafbfs`, `daffna`, `dafbbs`, `daffpa`) and the array that is being added (`dafbna`,
`dafada`, `dafena`). Two handles are independent, so two tasks can work on two handles at
the same time. One handle is not safe for concurrent use: tasks that share a handle must
hold a lock around every call. A file that an open handle has open, for read or for write,
cannot be opened for write by another handle (`SpiceError(:fileopenconflict)`), as in
CSPICE.
Unlike CSPICE, a file that is open for write can be opened for read, and a kernel that
[`furnsh`](@ref) loaded can be opened for write.

A call on a closed handle throws `SpiceError(:dafnosuchhandle)`; the close calls do
nothing instead. A write call on a handle that is open for reading throws
`SpiceError(:dafinvalidaccess)`.

A write call puts its bytes on disk before it returns, except [`spk14b`](@ref) and
[`spk14a`](@ref), whose data [`spk14e`](@ref) writes. A write call that throws for an
invalid input does not change the file.
"""
mutable struct DAFHandle
    const path::String
    const io::IOStream
    const writable::Bool
    const swap::Bool                 # the file has the other byte order (read only)
    const nd::Int
    const ni::Int
    isopen::Bool
    # --- array search --------------------------------------------------------------
    searching::Bool                  # `dafbfs` or `dafbbs` was called
    sthis::Int                       # record number of the buffered summary record, or 0
    snext::Int
    sprev::Int
    snsum::Int                       # summaries in the buffered record
    scurr::Int                       # index of the current array in the buffered record
    const sbuf::Vector{Float64}      # the buffered summary record, native byte order
    # --- array that is being added ---------------------------------------------------
    adding::Bool
    abegin::Int                      # first address of the new array
    afree::Int                       # next free address while the array grows
    asum::Vector{Float64}
    aname::Vector{UInt8}
end

function Base.show(io::IO, h::DAFHandle)
    print(io, "DAFHandle(", repr(h.path), ", ", h.writable ? "write" : "read",
        h.isopen ? "" : ", closed", ")")
end
Base.isopen(h::DAFHandle) = h.isopen
Base.close(h::DAFHandle) = dafcls(h)

"Doubles in one summary of a DAF with the format `nd`, `ni`."
_daf_sumsize(nd::Integer, ni::Integer) = Int(nd) + (Int(ni) + 1) ÷ 2

@inline function _daf_check_open(h::DAFHandle, who::AbstractString)
    h.isopen || throw(SpiceError(:dafnosuchhandle, "$who: the DAF $(h.path) is closed"))
    return nothing
end

@inline function _daf_check_write(h::DAFHandle, who::AbstractString)
    _daf_check_open(h, who)
    h.writable || throw(SpiceError(:dafinvalidaccess,
        "$who: the DAF $(h.path) is not open for write"))
    return nothing
end

# ---------------------------------------------------------------------------------------
# Record input
# ---------------------------------------------------------------------------------------

"""
    _daf_read_at!(h, buf, offset) -> Int

Read up to `length(buf)` bytes of the file of `h` from the 0-based byte `offset` into
`buf` and return the number of bytes read (less than `length(buf)` at the end of file).
"""
function _daf_read_at!(h::DAFHandle, buf::Vector{UInt8}, offset::Integer)
    offset >= filesize(h.io) && return 0
    seek(h.io, offset)
    return readbytes!(h.io, buf)
end

"The 1024 bytes of record `recno` (zero-filled past the end of file) and the count of bytes read."
function _daf_raw_record(h::DAFHandle, recno::Integer)
    raw = zeros(UInt8, _DAF_RECLEN)
    recno >= 1 || return raw, 0
    n = _daf_read_at!(h, raw, (Int(recno) - 1) * _DAF_RECLEN)
    return raw, n
end

"""
    _daf_drec(h, recno) -> (rec::Vector{Float64}, found::Bool)

Record `recno` as 128 doubles in native byte order (every double is converted if the
file has the other byte order; use `_daf_sumrec` for a summary record). `found` is
false, and `rec` is zero, if the record is not wholly in the file.
"""
function _daf_drec(h::DAFHandle, recno::Integer)
    raw, n = _daf_raw_record(h, recno)
    n == _DAF_RECLEN || return zeros(Float64, _DAF_RECWORDS), false
    rec = collect(reinterpret(Float64, raw))
    if h.swap
        map!(bswap, rec, rec)
    end
    return rec, true
end

"""
    _daf_sumrec(h, recno) -> (rec::Vector{Float64}, found::Bool)

Summary record `recno` in native byte order: the three control doubles and the `nd`
doubles of every summary are converted as doubles and the packed integers as 32-bit
integers, so that `dafus` unpacks the result on this machine.
"""
function _daf_sumrec(h::DAFHandle, recno::Integer)
    raw, n = _daf_raw_record(h, recno)
    n == _DAF_RECLEN || return zeros(Float64, _DAF_RECWORDS), false
    if h.swap
        ss = _daf_sumsize(h.nd, h.ni)
        for w in 0:_DAF_RECWORDS-1       # invariant: words before w are converted
            base = 8w
            k = w < 3 ? 0 : (w - 3) % ss         # position of the word inside its summary
            if w < 3 || k < h.nd
                reverse!(raw, base + 1, base + 8)
            else
                reverse!(raw, base + 1, base + 4)
                reverse!(raw, base + 5, base + 8)
            end
        end
    end
    return collect(reinterpret(Float64, raw)), true
end

"The first 1000 bytes of character record `recno` (blank-filled past the end of file) and whether all were read."
function _daf_crec(h::DAFHandle, recno::Integer)
    rec = fill(UInt8(' '), _DAF_CRECLEN)
    recno >= 1 || return rec, false
    n = _daf_read_at!(h, rec, (Int(recno) - 1) * _DAF_RECLEN)
    return rec, n == _DAF_CRECLEN
end

"Record number and word (both 1-based) of the DAF address `addr`, as `dafarw`."
_daf_arw(addr::Integer) = ((Int(addr) - 1) ÷ _DAF_RECWORDS + 1, (Int(addr) - 1) % _DAF_RECWORDS + 1)

"First address of record `recno`, as `dafrwa` with word 1."
_daf_rwa(recno::Integer) = (Int(recno) - 1) * _DAF_RECWORDS + 1

# ---------------------------------------------------------------------------------------
# File record and file architecture
# ---------------------------------------------------------------------------------------

"32-bit integer at 0-based byte offset `off` of the file record `raw` of `h`."
_daf_frint(raw::Vector{UInt8}, off::Int, swap::Bool) =
    Int(_daf_i32(raw, off, swap ? !_DAF_NATIVE_LITTLE : _DAF_NATIVE_LITTLE))

"""
    idw2at(idword) -> (arch, kind)

Architecture and type of a kernel from its ID word (the first eight characters of the
file), as `idw2at` of SPICELIB: `("DAF", "SPK")` for `DAF/SPK`, `("DAF", "?")` for the
old `NAIF/DAF`, `("KPL", "LSK")` for `KPL/LSK`, and `("?", "?")` if the word is not
recognized.
"""
function idw2at(idword::AbstractString)
    word = String(idword)
    slash = findfirst('/', word)
    (isempty(strip(word)) || slash === nothing) && return "?", "?"
    part1 = String(rstrip(word[1:prevind(word, slash)]))
    rest = word[nextind(word, slash):end]
    part2 = String(rstrip(rest[1:min(lastindex(rest), 8)]))
    kind = isempty(strip(part2)) ? "?" : part2
    if part1 in ("DAF", "DAS", "TXT", "KPL")
        return part1, kind
    elseif part1 == "ASC"
        return "TXT", kind
    elseif part1 == "NAIF"
        (part2 == "DAF" || part2 == "NIP") && return "DAF", "?"
        part2 == "DAS" && return "DAS", "PRE"
    end
    return "?", "?"
end

"""
    _daf_probe(path, who) -> (nd, ni, swap)

Check that `path` is an existing DAF file that this machine can read, and return its
summary format and whether its byte order is the other one. Throws `SpiceError` with
kind `:emptystring`, `:blankfilename`, `:filenotfound`, `:filereadfailed`,
`:unknownfilarc`, `:filarchmismatch`, `:ftpxfererror`, `:unknownbff`, `:unsupportedbff` or
`:dafinvalidparams`.
"""
function _daf_probe(path::AbstractString, who::AbstractString)
    _chk_filename(who, path)
    isfile(path) || throw(SpiceError(:filenotfound, "$who: the file $path does not exist"))
    raw, little = open(path, "r") do io
        raw = read(io, _DAF_RECLEN)
        length(raw) == _DAF_RECLEN || throw(SpiceError(:filereadfailed,
            "$who: $path is too short ($(length(raw)) bytes) to hold a DAF file record"))
        arch, _ = idw2at(String(raw[1:8]))
        arch == "DAS" && throw(SpiceError(:filarchmismatch,
            "$who: the file $path is a DAS file, not a DAF"))
        arch == "DAF" || throw(SpiceError(:unknownfilarc,
            "$who: the file $path is neither a DAF nor a DAS file"))
        if raw[700:705] == _DAF_FTPSTR[1:6] && raw[700:727] != _DAF_FTPSTR
            throw(SpiceError(:ftpxfererror, "$who: the FTP validation string of $path is " *
                "damaged; the file was most likely moved in text mode"))
        end
        raw, _daf_little(io, raw, path, who)
    end
    swap = little != _DAF_NATIVE_LITTLE
    nd = Int(_daf_i32(raw, 8, little))
    ni = Int(_daf_i32(raw, 12, little))
    (0 <= nd <= 124 && 2 <= ni <= 250 && _daf_sumsize(nd, ni) <= 125) || throw(SpiceError(
        :dafinvalidparams, "$who: ND = $nd, NI = $ni of $path are not a DAF summary format"))
    return nd, ni, swap
end

function _daf_handle(path::AbstractString, io::IOStream, writable::Bool, swap::Bool,
        nd::Int, ni::Int)
    return DAFHandle(String(path), io, writable, swap, nd, ni, true,
        false, 0, 0, 0, 0, 0, zeros(Float64, _DAF_RECWORDS),
        false, 0, 0, Float64[], UInt8[])
end

"""
    dafopr(path) -> DAFHandle

Open the DAF file `path` for reading and return its handle, as `dafopr_c`.

The file can have either byte order. No array search is begun. Close the handle with
[`dafcls`](@ref).

Throws `SpiceError`:
- `:emptystring` if `path` is empty, `:blankfilename` if it holds only blanks;
- `:filenotfound` if the file does not exist;
- `:filereadfailed` if the file is shorter than one record (1024 bytes);
- `:filarchmismatch` if the file is a DAS file, `:unknownfilarc` if it is neither
  a DAF nor a DAS file (a text kernel, for example);
- `:ftpxfererror` if the file was damaged by a transfer in text mode;
- `:unknownbff` if the binary file format (the number format that the file record
  names) is not known, `:unsupportedbff` if it is `VAX-GFLT` or `VAX-DFLT`;
- `:dafinvalidparams` if the summary format in the file record is not valid;
- `:fileopenfailed` if the system cannot open the file.

# Examples
Write an SPK file with one segment, then search it as a DAF:
```jldoctest
julia> path = joinpath(mktempdir(), "example.bsp");

julia> h = spkopn(path, "example SPK");

julia> cdata = [[1.0, 0.5, 2.0, 0.0, 3.0, 0.0]];

julia> spkw02(h, 499, 4, "J2000", 0.0, 100.0, "mars", 100.0, cdata, 0.0)

julia> spkcls(h)

julia> h = dafopr(path);

julia> dafbfs(h)

julia> daffna(h)
true

julia> dc, ic = dafus(dafgs(h), 2, 6)   # see dafgs for the meaning of the components
([0.0, 100.0], [499, 4, 1, 2, 385, 396])

julia> dafgn(h)
"mars"

julia> dafgda(h, ic[5], ic[5] + 1)      # midpoint and radius of the first record, s
2-element Vector{Float64}:
 50.0
 50.0

julia> daffna(h)
false

julia> dafcls(h)
```

See also [`dafopw`](@ref), [`dafbfs`](@ref), [`dafgda`](@ref), [`dafec`](@ref),
[`getfat`](@ref).
"""
function dafopr(path::AbstractString)
    nd, ni, swap = _daf_probe(path, "dafopr")
    io = try
        open(path, "r")
    catch e
        throw(SpiceError(:fileopenfailed, "dafopr: cannot open $path: $(sprint(showerror, e))"))
    end
    return _open_register!(_daf_handle(path, io, false, swap, nd, ni), path, false, "dafopr")
end

"""
    dafcls(handle) -> Nothing

Close the DAF file of `handle`, as `dafcls_c`.

The close writes no data: every write call puts its bytes on disk before it returns,
except [`spk14b`](@ref) and [`spk14a`](@ref) (see [`DAFHandle`](@ref)). An array that
`dafbna` began and `dafena` did not end has no summary in the file, so a search does not
find it. The call does nothing if the handle is closed already. After the call the
handle keeps the file neither open nor memory-mapped, so that Windows can delete or
replace it.
"""
function dafcls(h::DAFHandle)
    h.isopen || return nothing
    h.isopen = false
    h.searching = false
    h.adding = false
    _open_unregister!(h)
    close(h.io)
    return nothing
end

"""
    dafhsf(handle) -> (nd, ni)

Return the summary format of the DAF file, as `dafhsf` of SPICELIB: the number of double
precision components `nd` and of integer components `ni` of an array summary.
"""
function dafhsf(h::DAFHandle)
    _daf_check_open(h, "dafhsf")
    return h.nd, h.ni
end

"Raw file record of `h` and its fields `(fward, bward, free)`; throws `SpiceError(:daffrnotfound)`."
function _daf_filerecord(h::DAFHandle, who::AbstractString)
    raw, n = _daf_raw_record(h, 1)
    n == _DAF_RECLEN || throw(SpiceError(:daffrnotfound,
        "$who: the file record of $(h.path) cannot be read"))
    return raw, _daf_frint(raw, 76, h.swap), _daf_frint(raw, 80, h.swap), _daf_frint(raw, 84, h.swap)
end

"""
    dafrfr(handle) -> (nd, ni, ifname, fward, bward, free)

Return the contents of the file record, as `dafrfr_c`: the summary format `nd`, `ni`, the
internal file name without trailing blanks, the record numbers of the first and of the
last summary record, and the first free address, at which the next array begins.

Throws `SpiceError(:daffrnotfound)` if the file record cannot be read.
"""
function dafrfr(h::DAFHandle)
    _daf_check_open(h, "dafrfr")
    raw, fward, bward, free = _daf_filerecord(h, "dafrfr")
    ifname = String(rstrip(String(raw[17:76]), [' ', '\0']))
    return h.nd, h.ni, ifname, fward, bward, free
end

"""
    getfat(path) -> (arch, kind)

Return the file architecture and the kernel type of the file `path`, as `getfat_c`.

Results are, for example, `("DAF", "SPK")`, `("DAF", "CK")`, `("DAF", "PCK")`,
`("DAS", "DSK")` and `("KPL", "LSK")`. The result is `("?", "?")` for a file that is not
recognized; a text kernel without the `KPL/<type>` first line is one.

Throws `SpiceError(:emptystring)` if `path` is empty, `:blankfilename` if it holds only
blanks, `:filenotfound` if the file does not exist, and `:filereadfailed` if the file has
fewer than 12 bytes and no complete line.
"""
function getfat(path::AbstractString)
    _chk_filename("getfat", path)
    isfile(path) || throw(SpiceError(:filenotfound, "getfat: the kernel file $path does not exist"))
    head = open(io -> read(io, 12), path, "r")
    # A file of fewer than 12 bytes that has no complete line cannot be read by CSPICE.
    (length(head) < 12 && !(UInt8('\n') in head)) && throw(SpiceError(:filereadfailed,
        "getfat: attempt to read the first 12 bytes of the file $path failed"))
    for i in eachindex(head)
        (head[i] < 0x20 || head[i] > 0x7e) && (head[i] = UInt8(' '))
    end
    words = split(uppercase(String(head)))
    isempty(words) && return "?", "?"
    idword = String(words[1])
    idword == "DAFETF" && return "XFR", "DAF"
    idword == "DASETF" && return "XFR", "DAS"
    startswith(idword, "'NAIF/DAF'") && return "DEC", "DAF"
    idword == "NAIF/DAS" && return "DAS", "PRE"
    arch, kind = idw2at(idword[1:min(lastindex(idword), 8)])
    if arch == "DAF" && kind == "?"
        h = dafopr(path)
        try
            kind = _daf_ckspk(h)
        finally
            dafcls(h)
        end
    end
    return arch, kind
end

# `_daf_sizeok` is in `daf.jl`.

"Type `\"SPK\"`, `\"CK\"` or `\"?\"` of an old DAF without a type in its ID word, as `zzckspk`."
function _daf_ckspk(h::DAFHandle)
    (h.nd == 2 && h.ni == 6) || return "?"
    dafbfs(h)
    daffna(h) || return "?"
    _, ic = dafus(dafgs(h), 2, 6)
    angvel = ic[4]
    dtype = ic[3]
    angvel == 0 && return "CK"
    angvel > 1 && return "SPK"
    ic[2] == 0 && return "SPK"
    first = ic[5]
    last = ic[6]
    size = last - first + 1
    spkok, nspk = _daf_sizeok(size - 1, 72, 100, 0)
    spkok || return "CK"
    dafgda(h, last, last)[1] == nspk || return "CK"
    dtype == 2 || return "SPK"
    ck2ok, nck2 = _daf_sizeok(size, 10, 100, 1)
    ck2ok || return "SPK"
    if nck2 < 201
        frsttm = dafgda(h, last - nspk, last - nspk)[1]
        from = first + (nspk - 1) * 71
        return dafgda(h, from, from)[1] > frsttm ? "SPK" : "CK"
    end
    from = last - (nck2 - 1) ÷ 100
    times = dafgda(h, from, from + 1)
    return times[1] > times[2] ? "CK" : "SPK"
end

# ---------------------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------------------

"""
    dafps(dc, ic) -> Vector{Float64}

Pack the double precision components `dc` and the integer components `ic` of an array
summary into `length(dc) + (length(ic) + 1) ÷ 2` doubles, as `dafps_c`.

The integers are stored as 32-bit values, two in the space of one double, in native byte
order. [`dafus`](@ref) is the inverse.

Throws `SpiceError(:dafinvalidparams)` if an integer does not fit in 32 bits.

# Examples
```jldoctest
julia> packed = dafps([0.0, 100.0], [499, 4, 1, 2, 385, 396]);

julia> length(packed)
5

julia> dafus(packed, 2, 6)
([0.0, 100.0], [499, 4, 1, 2, 385, 396])
```
"""
function dafps(dc::AbstractVector{<:Real}, ic::AbstractVector{<:Integer})
    nd = length(dc)
    ni = length(ic)
    all(i -> typemin(Int32) <= i <= typemax(Int32), ic) || throw(SpiceError(:dafinvalidparams,
        "dafps: an integer component does not fit in 32 bits: $(collect(ic))"))
    sum = zeros(Float64, _daf_sumsize(nd, ni))
    for j in 1:nd
        sum[j] = Float64(dc[j])
    end
    ints = reinterpret(Int32, sum)
    for j in 1:ni                        # invariant: integers 1 .. j-1 are packed
        ints[2nd+j] = Int32(ic[j])
    end
    return sum
end

"""
    dafus(sum, nd, ni) -> (dc::Vector{Float64}, ic::Vector{Int})

Unpack the array summary `sum` into its `nd` double precision components and its `ni`
integer components, as `dafus_c`.

Throws `SpiceError(:dafinvalidparams)` if `sum` has fewer than `nd + (ni + 1) ÷ 2`
doubles.

See also [`dafps`](@ref), [`dafgs`](@ref).
"""
function dafus(sum::AbstractVector{<:Real}, nd::Integer, ni::Integer)
    (nd >= 0 && ni >= 0 && length(sum) >= _daf_sumsize(nd, ni)) || throw(SpiceError(
        :dafinvalidparams, "dafus: a summary of $(length(sum)) doubles does not hold " *
        "ND = $nd, NI = $ni"))
    s = Float64[Float64(x) for x in sum]
    dc = s[1:nd]
    ints = reinterpret(Int32, s)
    ic = Int[ints[2nd+j] for j in 1:ni]
    return dc, ic
end

# ---------------------------------------------------------------------------------------
# Array search
# ---------------------------------------------------------------------------------------

"Load summary record `recno` into the search buffer of `h`; throws `SpiceError(:recordnotfound)`."
function _daf_load_sumrec!(h::DAFHandle, recno::Int, who::AbstractString)
    rec, found = _daf_sumrec(h, recno)
    found || throw(SpiceError(:recordnotfound,
        "$who: summary record $recno is not in the file $(h.path)"))
    copyto!(h.sbuf, rec)
    maxsum = 125 ÷ _daf_sumsize(h.nd, h.ni)
    nsum = rec[3]
    (isfinite(rec[1]) && isfinite(rec[2]) && 0 <= nsum <= maxsum && nsum == trunc(nsum) &&
        rec[1] >= 0 && rec[2] >= 0) || throw(SpiceError(:baddaf,
        "$who: summary record $recno of $(h.path) has the control words $(rec[1:3])"))
    h.sthis = recno
    h.snext = Int(rec[1])
    h.sprev = Int(rec[2])
    h.snsum = Int(nsum)
    return nothing
end

"""
    dafbfs(handle) -> Nothing

Begin a forward search for the arrays of the DAF file, as `dafbfs_c`.

No array is current after the call; the first [`daffna`](@ref) moves to the first array.

See also [`dafbbs`](@ref), and [`dafopr`](@ref) for an example.
"""
function dafbfs(h::DAFHandle)
    _daf_check_open(h, "dafbfs")
    _, fward, _, _ = _daf_filerecord(h, "dafbfs")
    h.searching = true
    h.sthis = 0
    h.snext = fward
    h.sprev = 0
    h.snsum = 0
    h.scurr = 0
    return nothing
end

"""
    dafbbs(handle) -> Nothing

Begin a backward search for the arrays of the DAF file, as `dafbbs_c`.

No array is current after the call; the first [`daffpa`](@ref) moves to the last array.

See also [`dafbfs`](@ref).
"""
function dafbbs(h::DAFHandle)
    _daf_check_open(h, "dafbbs")
    _, _, bward, _ = _daf_filerecord(h, "dafbbs")
    h.searching = true
    h.sthis = 0
    h.snext = 0
    h.sprev = bward
    h.snsum = 0
    h.scurr = 1
    return nothing
end

function _daf_check_search(h::DAFHandle, who::AbstractString)
    _daf_check_open(h, who)
    h.searching || throw(SpiceError(:dafnosearch,
        "$who: no search was begun on the DAF $(h.path); call dafbfs or dafbbs first"))
    return nothing
end

"""
    daffna(handle) -> Bool

Move to the next array of the search on `handle`, as `daffna_c`. Return `true` if there
is one and `false` if the search is past the last array.

Throws `SpiceError(:dafnosearch)` if neither `dafbfs` nor `dafbbs` was called on the
handle, and `:recordnotfound` or `:baddaf` if the chain of summary records is damaged.

See also [`daffpa`](@ref), [`dafgs`](@ref), [`dafgn`](@ref).
"""
function daffna(h::DAFHandle)
    _daf_check_search(h, "daffna")
    h.scurr += 1
    steps = 0
    while h.scurr > h.snsum && h.snext != 0
        # invariant: every summary record before `snext` in the chain is exhausted
        steps += 1
        steps <= 1 + filesize(h.io) ÷ _DAF_RECLEN || throw(SpiceError(:baddaf,
            "daffna: the chain of summary records has a cycle in $(h.path)"))
        _daf_load_sumrec!(h, h.snext, "daffna")
        h.scurr = 1
    end
    if h.scurr > h.snsum
        h.scurr = h.snsum + 1
        return false
    end
    return true
end

"""
    daffpa(handle) -> Bool

Move to the previous array of the search on `handle`, as `daffpa_c`. Return `true` if
there is one and `false` if the search is before the first array.

Throws `SpiceError(:dafnosearch)` if neither `dafbfs` nor `dafbbs` was called on the
handle, and `:recordnotfound` or `:baddaf` if the chain of summary records is damaged.

See also [`daffna`](@ref).
"""
function daffpa(h::DAFHandle)
    _daf_check_search(h, "daffpa")
    h.scurr -= 1
    steps = 0
    while h.scurr <= 0 && h.sprev != 0
        # invariant: every summary record after `sprev` in the chain is exhausted
        steps += 1
        steps <= 1 + filesize(h.io) ÷ _DAF_RECLEN || throw(SpiceError(:baddaf,
            "daffpa: the chain of summary records has a cycle in $(h.path)"))
        _daf_load_sumrec!(h, h.sprev, "daffpa")
        h.scurr = h.snsum
    end
    if h.scurr <= 0
        h.scurr = 0
        return false
    end
    return true
end

function _daf_check_current(h::DAFHandle, who::AbstractString)
    _daf_check_search(h, who)
    (h.sthis != 0 && 1 <= h.scurr <= h.snsum) || throw(SpiceError(:nocurrentarray,
        "$who: no array is current in the search on $(h.path)"))
    return nothing
end

"""
    dafgs(handle) -> Vector{Float64}

Return the packed summary of the current array of the search on `handle`, as `dafgs_c`.

The result has `nd + (ni + 1) ÷ 2` doubles; unpack it with [`dafus`](@ref). An SPK file
and a CK file have `nd = 2` and `ni = 6`:

| File | `dc`                                        | `ic`                                                   |
|:-----|:--------------------------------------------|:-------------------------------------------------------|
| SPK  | first and last epoch, seconds past J2000 TDB | body, center, frame ID, data type, first and last address |
| CK   | first and last time, encoded SCLK ticks      | instrument, frame ID, data type, angular velocity flag, first and last address |

The addresses are those of the data of the array, for [`dafgda`](@ref).

Throws `SpiceError(:dafnosearch)` if no search was begun, and `:nocurrentarray` unless the
last `daffna` or `daffpa` returned `true`.

See also [`dafgn`](@ref), and [`dafopr`](@ref) for an example.
"""
function dafgs(h::DAFHandle)
    _daf_check_current(h, "dafgs")
    ss = _daf_sumsize(h.nd, h.ni)
    off = 3 + (h.scurr - 1) * ss
    return h.sbuf[off+1:off+ss]
end

"""
    dafgn(handle) -> String

Return the name of the current array of the search on `handle`, without trailing blanks,
as `dafgn_c`.

Throws `SpiceError(:dafnosearch)` if no search was begun, and `:nocurrentarray` unless the
last `daffna` or `daffpa` returned `true`.
"""
function dafgn(h::DAFHandle)
    _daf_check_current(h, "dafgn")
    namsiz = 8 * _daf_sumsize(h.nd, h.ni)
    names, _ = _daf_crec(h, h.sthis + 1)
    off = (h.scurr - 1) * namsiz
    return String(rstrip(String(names[off+1:off+namsiz]), [' ', '\0']))
end

# ---------------------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------------------

"""
    dafgda(handle, first, last) -> Vector{Float64}

Return the doubles at the DAF addresses `first:last`, as `dafgda_c`. Addresses start at 1,
as in CSPICE.

Two results differ from CSPICE:
- A range that is not inside the file throws. CSPICE returns zeros without an error.
- If the file ends inside its last record, the doubles of that record that the file
  holds are returned. CSPICE returns zeros for the whole record.

Throws `SpiceError(:dafnegaddr)` if `first < 1`, `:dafbeggtend` if `first > last`, and
`:filereadfailed` if the range is not inside the file.
"""
function dafgda(h::DAFHandle, first::Integer, last::Integer)
    _daf_check_open(h, "dafgda")
    first <= 0 && throw(SpiceError(:dafnegaddr, "dafgda: negative value for the first address: $first"))
    first > last && throw(SpiceError(:dafbeggtend,
        "dafgda: beginning address ($first) greater than ending address ($last)"))
    n = Int(last) - Int(first) + 1
    data = Vector{Float64}(undef, n)
    8 * Int(last) <= filesize(h.io) || throw(SpiceError(:filereadfailed,
        "dafgda: the addresses $first:$last are not inside the file $(h.path)"))
    seek(h.io, 8 * (Int(first) - 1))
    read!(h.io, data)
    if h.swap
        map!(bswap, data, data)
    end
    return data
end

"""
    dafgsr(handle, recno, first, last) -> Union{Vector{Float64},Nothing}

Return the doubles `first:last` of summary record `recno`, as `dafgsr_c`, or `nothing`
if the record is not in the file.

A record has 128 doubles, numbered from 1. `nothing` replaces the `found` flag of CSPICE.

Throws `SpiceError(:dafinvalidparams)` unless `1 <= first <= last <= 128`.
"""
function dafgsr(h::DAFHandle, recno::Integer, first::Integer, last::Integer)
    _daf_check_open(h, "dafgsr")
    (1 <= first <= last <= _DAF_RECWORDS) || throw(SpiceError(:dafinvalidparams,
        "dafgsr: the words $first:$last are not in 1:$(_DAF_RECWORDS)"))
    rec, found = _daf_sumrec(h, recno)
    found || return nothing
    return rec[first:last]
end

# ---------------------------------------------------------------------------------------
# Comment area
# ---------------------------------------------------------------------------------------

"""
    _daf_comment_extent(h, fward, who) -> (ncomc, ncomr, notused)

Number of comment characters before the end-of-comments mark, the number of the
reserved record that holds the mark (counted from 1 at record 2; 0 if the area is
empty) and the number of all-null reserved records after it, as `dafac` and `dafec`
find them. Throws `SpiceError(:badcommentarea)` if the mark is missing from a
comment area that is not empty.
"""
function _daf_comment_extent(h::DAFHandle, fward::Int, who::AbstractString)
    ncomr = fward - 2
    notused = 0
    while ncomr > 0                      # invariant: reserved records after ncomr are all null
        rec, _ = _daf_crec(h, ncomr + 1)
        eocpos = findfirst(==(_DAF_EOC), rec)
        eocpos !== nothing && return (ncomr - 1) * _DAF_CRECLEN + eocpos - 1, ncomr, notused
        all(==(_DAF_EOL), rec) || throw(SpiceError(:badcommentarea,
            "$who: the comment area of $(h.path) may be damaged; the end of the " *
            "comments could not be found"))
        ncomr -= 1
        notused += 1
    end
    return 0, 0, notused
end

"""
    dafec(handle) -> Vector{String}

Return all lines of the comment area of the DAF file, in order.

One call returns the whole comment area, where `dafec_c` returns one buffer of lines on
each call. The result is empty if the file has no comments.

Throws `SpiceError(:badcommentarea)` if the comment area is damaged.

See also [`dafac`](@ref), [`dafdc`](@ref).
"""
function dafec(h::DAFHandle)
    _daf_check_open(h, "dafec")
    _, fward, _, _ = _daf_filerecord(h, "dafec")
    ncomc, _, _ = _daf_comment_extent(h, fward, "dafec")
    lines = String[]
    ncomc == 0 && return lines
    chars = UInt8[]
    recno = 2
    while length(chars) < ncomc          # invariant: chars holds records 2 .. recno-1
        rec, _ = _daf_crec(h, recno)
        append!(chars, rec)
        recno += 1
    end
    resize!(chars, ncomc)
    chars[end] == _DAF_EOL || throw(SpiceError(:badcommentarea,
        "dafec: the last comment line of $(h.path) is not terminated"))
    start = 1
    for k in 1:ncomc                     # invariant: lines holds the lines that end before k
        if chars[k] == _DAF_EOL
            push!(lines, String(rstrip(String(chars[start:k-1]))))
            start = k + 1
        end
    end
    return lines
end
