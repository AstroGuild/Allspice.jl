# Low-level access to DAF files: open for read, the array search, summaries, names,
# data, the comment area, the file record, and the file architecture (`getfat`).
#
# CSPICE keeps a process-wide table of open DAF files and one hidden "current array"
# per search. Here the opened file is an object, `DAFHandle`, that the caller owns, and
# the search state lives inside it. There is no global state, so two tasks can work on
# two handles at the same time. One handle is not safe for concurrent use: a caller
# that shares a handle between tasks must hold a lock around every call.
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

"""
    DAFHandle

One DAF file opened for access with `dafopr`, `dafopw`, `dafonw`, `spkopn`, `spkopa`
or `ckopn`. It replaces the integer handle of CSPICE and holds what CSPICE keeps in
global tables: the open file, the state of the array search (`dafbfs`, `daffna`,
`dafbbs`, `daffpa`) and the state of the array that is being added (`dafbna`,
`dafada`, `dafena`).

Invariants: `nd`, `ni` are the summary format of the file; while `isopen`, `io` is
open for reading and, if `writable`, for writing; a writable handle has the native
byte order (`swap == false`). A handle must not be used by two tasks at the same time
without a lock that the caller holds; different handles are independent.
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

"Throws `SpiceError(:dafnosuchhandle)` unless `h` is open."
@inline function _daf_check_open(h::DAFHandle, who::AbstractString)
    h.isopen || throw(SpiceError(:dafnosuchhandle, "$who: the DAF $(h.path) is closed"))
    return nothing
end

"Throws unless `h` is open for writing (`:dafnosuchhandle`, `:dafinvalidaccess`)."
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
kind `:blankfilename`, `:filenotfound`, `:filereadfailed`, `:unknownfilarc`,
`:filarchmismatch`, `:ftpxfererror`, `:unknownbff` or `:dafinvalidparams`.
"""
function _daf_probe(path::AbstractString, who::AbstractString)
    isempty(strip(path)) && throw(SpiceError(:blankfilename, "$who: the file name is blank"))
    isfile(path) || throw(SpiceError(:filenotfound, "$who: the file $path does not exist"))
    raw = open(io -> read(io, _DAF_RECLEN), path, "r")
    length(raw) == _DAF_RECLEN || throw(SpiceError(:filereadfailed,
        "$who: $path is too short ($(length(raw)) bytes) to hold a DAF file record"))
    arch, _ = idw2at(String(raw[1:8]))
    arch == "?" && throw(SpiceError(:unknownfilarc,
        "$who: the file $path has an unidentified file architecture"))
    arch == "DAF" || throw(SpiceError(:filarchmismatch,
        "$who: the file $path has the architecture $arch, not DAF"))
    if raw[700:705] == _DAF_FTPSTR[1:6] && raw[700:727] != _DAF_FTPSTR
        throw(SpiceError(:ftpxfererror, "$who: the FTP validation string of $path is " *
            "damaged; the file was most likely moved in text mode"))
    end
    bff = rstrip(String(raw[89:96]))
    little = if bff == "LTL-IEEE"
        true
    elseif bff == "BIG-IEEE"
        false
    else
        # Files older than the binary-file-format string: ND tells the byte order.
        ndl = Int(_daf_i32(raw, 8, true))
        ndb = Int(_daf_i32(raw, 8, false))
        if all(b -> b == 0x00 || b == UInt8(' '), raw[89:96]) && (1 <= ndl <= 124) != (1 <= ndb <= 124)
            1 <= ndl <= 124
        else
            throw(SpiceError(:unknownbff,
                "$who: the binary file format $(repr(bff)) of $path is not known"))
        end
    end
    swap = little != _DAF_NATIVE_LITTLE
    nd = Int(_daf_i32(raw, 8, little))
    ni = Int(_daf_i32(raw, 12, little))
    (0 <= nd <= 124 && 2 <= ni <= 250 && _daf_sumsize(nd, ni) <= 125) || throw(SpiceError(
        :dafinvalidparams, "$who: ND = $nd, NI = $ni of $path are not a DAF summary format"))
    return nd, ni, swap
end

"New handle on the open stream `io` of the DAF `path`."
function _daf_handle(path::AbstractString, io::IOStream, writable::Bool, swap::Bool,
        nd::Int, ni::Int)
    return DAFHandle(String(path), io, writable, swap, nd, ni, true,
        false, 0, 0, 0, 0, 0, zeros(Float64, _DAF_RECWORDS),
        false, 0, 0, Float64[], UInt8[])
end

"""
    dafopr(path) -> DAFHandle

Open the DAF file `path` for read access, as `dafopr_c`. A file of the other byte
order can be read. Contract: `path` names an existing DAF; the result is an open
handle with no search begun. Throws `SpiceError` (`:blankfilename`, `:filenotfound`,
`:filereadfailed`, `:unknownfilarc`, `:filarchmismatch`, `:ftpxfererror`,
`:unknownbff`). The handle belongs to the caller: close it with `dafcls`, and do not
share it between tasks without a lock.
"""
function dafopr(path::AbstractString)
    nd, ni, swap = _daf_probe(path, "dafopr")
    io = try
        open(path, "r")
    catch e
        throw(SpiceError(:fileopenfailed, "dafopr: cannot open $path: $(sprint(showerror, e))"))
    end
    return _daf_handle(path, io, false, swap, nd, ni)
end

"""
    dafcls(handle) -> Nothing

Close the DAF of `handle`, as `dafcls_c`. All data are on disk when the call returns.
Closing a handle that is already closed does nothing. An array that was begun with
`dafbna` and not ended with `dafena` is not in the file.
"""
function dafcls(h::DAFHandle)
    h.isopen || return nothing
    h.isopen = false
    h.searching = false
    h.adding = false
    close(h.io)
    return nothing
end

"""
    dafhsf(handle) -> (nd, ni)

Summary format of the DAF: the number of double precision and of integer components of
a summary, as `dafhsf`.
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

Contents of the file record, as `dafrfr_c`: the summary format, the internal file
name (trailing blanks removed), the record numbers of the first and of the last
summary record, and the first free address.
"""
function dafrfr(h::DAFHandle)
    _daf_check_open(h, "dafrfr")
    raw, fward, bward, free = _daf_filerecord(h, "dafrfr")
    ifname = String(rstrip(String(raw[17:76]), [' ', '\0']))
    return h.nd, h.ni, ifname, fward, bward, free
end

"""
    getfat(path) -> (arch, kind)

File architecture and kernel type of the file `path`, as `getfat_c`: for example
`("DAF", "SPK")`, `("DAF", "CK")`, `("DAF", "PCK")`, `("KPL", "LSK")`, and
`("?", "?")` for a file that is not recognized (a text kernel without the `KPL/…`
first line is one). For an old `NAIF/DAF` file the type follows from the summary
format and the first array, as `zzckspk` decides it. Throws `SpiceError`
(`:blankfilename`, `:filenotfound`).
"""
function getfat(path::AbstractString)
    isempty(strip(path)) && throw(SpiceError(:blankfilename, "getfat: the file name is blank"))
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
summary into `length(dc) + (length(ic) + 1) ÷ 2` doubles, as `dafps_c`: the integers
are stored as 32-bit values, two in the space of one double, in native byte order.
Precondition: every integer fits in 32 bits, else `SpiceError(:dafinvalidparams)`.
Postcondition: `dafus(dafps(dc, ic), length(dc), length(ic)) == (dc, ic)`.
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

Unpack an array summary into its `nd` double precision and `ni` integer components,
as `dafus_c`. Precondition: `sum` holds at least `nd + (ni + 1) ÷ 2` doubles, else
`SpiceError(:dafinvalidparams)`.
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

Begin a forward search for arrays in the DAF, as `dafbfs_c`. Postcondition: no array is
current; the first `daffna(handle)` moves to the first array.
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

Begin a backward search for arrays in the DAF, as `dafbbs_c`. Postcondition: no array
is current; the first `daffpa(handle)` moves to the last array.
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

"Throws `SpiceError(:dafnosearch)` unless a search was begun on `h`."
function _daf_check_search(h::DAFHandle, who::AbstractString)
    _daf_check_open(h, who)
    h.searching || throw(SpiceError(:dafnosearch,
        "$who: no search was begun on the DAF $(h.path); call dafbfs or dafbbs first"))
    return nothing
end

"""
    daffna(handle) -> Bool

Move to the next array of the search on `handle` and return true, or return false if
there is no next array, as `daffna_c` (which has no handle argument and works on the
hidden current DAF). Precondition: `dafbfs` or `dafbbs` was called, else
`SpiceError(:dafnosearch)`.
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

Move to the previous array of the search on `handle` and return true, or return false
if there is none, as `daffpa_c`. Precondition: `dafbfs` or `dafbbs` was called, else
`SpiceError(:dafnosearch)`.
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

"Throws `SpiceError(:nocurrentarray)` unless the search on `h` has a current array."
function _daf_check_current(h::DAFHandle, who::AbstractString)
    _daf_check_search(h, who)
    (h.sthis != 0 && 1 <= h.scurr <= h.snsum) || throw(SpiceError(:nocurrentarray,
        "$who: no array is current in the search on $(h.path)"))
    return nothing
end

"""
    dafgs(handle) -> Vector{Float64}

Packed summary of the current array of the search on `handle`
(`nd + (ni + 1) ÷ 2` doubles), as `dafgs_c`; unpack it with `dafus`. Precondition: the
last `daffna` or `daffpa` returned true, else `SpiceError(:nocurrentarray)`.
"""
function dafgs(h::DAFHandle)
    _daf_check_current(h, "dafgs")
    ss = _daf_sumsize(h.nd, h.ni)
    off = 3 + (h.scurr - 1) * ss
    return h.sbuf[off+1:off+ss]
end

"""
    dafgn(handle) -> String

Name of the current array of the search on `handle`, trailing blanks removed, as
`dafgn_c`. Precondition: the last `daffna` or `daffpa` returned true, else
`SpiceError(:nocurrentarray)`.
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

The doubles at the DAF addresses `first:last`, as `dafgda_c`. Preconditions:
`first >= 1` (else `SpiceError(:dafnegaddr)`), `first <= last` (else
`SpiceError(:dafbeggtend)`), and the range is inside the file (else
`SpiceError(:filereadfailed)`; CSPICE returns zeros without an error in that case, and
also for all the words of a last record that the file holds in part, which this port
returns as they are).
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

Words `first:last` (1 to 128) of summary record `recno`, in native byte order, or
`nothing` if the record is not in the file, as `dafgsr_c`. Precondition:
`1 <= first <= last <= 128`, else `SpiceError(:dafinvalidparams)`.
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

All lines of the comment area of the DAF, in order, as repeated calls of `dafec_c`
give them (trailing blanks are not stored). The result is empty if the file has no
comments. Throws `SpiceError(:badcommentarea)` if the comment area is damaged.
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
