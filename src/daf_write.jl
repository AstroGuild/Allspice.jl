# Write access to DAF files: create a file, open one for update, add arrays, change
# summaries and names, and grow or delete the comment area. The routines follow
# `dafah.f` (`dafonw`, `dafopw`), `dafana.f` (`dafbna`, `dafada`, `dafena`), `daffa.f`
# (`dafws`, `dafrs`, `dafrn`), `dafac.f`, `dafarr.f`, `dafdc.f` and `dafrrr.f` of
# SPICELIB, so a file that these routines write has the bytes that CSPICE writes.
# The state of CSPICE that is global there is in the `DAFHandle` here.
#
# Files are written in the native byte order only. A file never shrinks: the bytes of
# an array are in place before its summary is, so a reader that opens the file between
# two calls sees a valid DAF without the array that is not yet ended.

# ---------------------------------------------------------------------------------------
# Record output
# ---------------------------------------------------------------------------------------

"""
    _daf_write_at!(h, bytes, offset) -> Nothing

Write `bytes` (bytes or doubles, as they are in memory) at the 0-based byte `offset` of
the file of `h` and flush. A gap between
the end of the file and `offset` is filled with zero bytes. Throws
`SpiceError(:dafwritefail)` if the system reports an error.
"""
function _daf_write_at!(h::DAFHandle, bytes::Union{Vector{UInt8},Vector{Float64}},
        offset::Integer)
    try
        size = filesize(h.io)
        if offset > size
            seek(h.io, size)
            write(h.io, zeros(UInt8, offset - size))
        end
        seek(h.io, offset)
        write(h.io, bytes)
        flush(h.io)
    catch e
        e isa SpiceError && rethrow()
        throw(SpiceError(:dafwritefail,
            "write to the DAF $(h.path) failed: $(sprint(showerror, e))"))
    end
    return nothing
end

"Write the 128 doubles `rec` as record `recno`, as `dafwdr`."
function _daf_write_drec!(h::DAFHandle, recno::Integer, rec::Vector{Float64})
    @assert length(rec) == _DAF_RECWORDS
    _daf_write_at!(h, rec, (Int(recno) - 1) * _DAF_RECLEN)
    return nothing
end

"Write the 1000 characters `rec` at the start of record `recno`, as `dafwcr` (bytes 1001 to 1024 are not touched)."
function _daf_write_crec!(h::DAFHandle, recno::Integer, rec::Vector{UInt8})
    @assert length(rec) == _DAF_CRECLEN
    _daf_write_at!(h, rec, (Int(recno) - 1) * _DAF_RECLEN)
    return nothing
end

"""
    _daf_write_data!(h, addr, data) -> Nothing

Write `data` at the DAF addresses `addr : addr + length(data) - 1`, as `dafwda`, and
fill the rest of the last record with zeros if the file ends inside it.
Postcondition: the file length is a whole number of records.
"""
function _daf_write_data!(h::DAFHandle, addr::Integer, data::Vector{Float64})
    isempty(data) && return nothing
    _daf_write_at!(h, data, 8 * (Int(addr) - 1))
    size = filesize(h.io)
    rest = mod(-size, _DAF_RECLEN)
    rest == 0 || _daf_write_at!(h, zeros(UInt8, rest), size)
    return nothing
end

"""
    _daf_write_filerecord!(h; ifname=nothing, fward, bward, free) -> Nothing

Update the internal file name (if given) and the three pointers of the file record
and leave all other bytes of the record as they are, as `dafwfr`.
"""
function _daf_write_filerecord!(h::DAFHandle; ifname::Union{Nothing,AbstractString}=nothing,
        fward::Integer, bward::Integer, free::Integer)
    (1 <= free <= typemax(Int32)) || throw(SpiceError(:dafwritefail,
        "the free address $free of $(h.path) is outside the 32-bit range of a DAF"))
    raw, _, _, _ = _daf_filerecord(h, "dafwfr")
    if ifname !== nothing
        raw[17:76] = _daf_padded(ifname, 60)
    end
    raw[77:88] = reinterpret(UInt8, Int32[fward, bward, free])
    _daf_write_at!(h, raw, 0)
    return nothing
end

"The code units of `s`, cut or padded with blanks to exactly `n` bytes (Fortran character assignment)."
function _daf_padded(s::AbstractString, n::Int)
    out = fill(UInt8(' '), n)
    cu = codeunits(String(s))
    copyto!(out, 1, cu, 1, min(n, length(cu)))
    return out
end

"Number of code units of `s` without trailing blanks, as `lastnb`."
function _daf_lastnb(s::AbstractString)
    cu = codeunits(String(s))
    n = length(cu)
    while n > 0 && cu[n] == UInt8(' ')
        n -= 1
    end
    return n
end

# ---------------------------------------------------------------------------------------
# Open for write
# ---------------------------------------------------------------------------------------

"""
    dafonw(path, ftype, nd, ni, ifname, resv) -> DAFHandle

Create the DAF file `path` and open it for writing, as `dafonw` of SPICELIB.

- `ftype`: the file type, for example `"SPK"`. Four characters are stored.
- `nd`, `ni`: the numbers of doubles and of integers in an array summary.
- `ifname`: the internal file name. 60 characters are stored.
- `resv`: the number of reserved records, which hold the comments.

The file is on disk when the call returns: the file record, `resv` comment records with an
empty comment area, one empty summary record and one name record, byte for byte as CSPICE
writes them. The byte order is the native one. An existing file is never replaced.

Throws `SpiceError`:
- `:blankfilename` if `path` is blank;
- `:fileopenfailed` if the file exists or cannot be created;
- `:blankfiletype` if `ftype` is blank, `:illegalcharacter` if it has a character that is
  not printable;
- `:dafinvalidparams` unless `0 <= nd <= 124`, `2 <= ni <= 250` and
  `nd + (ni + 1) ÷ 2 <= 125`;
- `:dafnoresv` if `resv < 0`.

See also [`spkopn`](@ref), [`ckopn`](@ref), [`dafbna`](@ref), [`dafcls`](@ref).
"""
function dafonw(path::AbstractString, ftype::AbstractString, nd::Integer, ni::Integer,
        ifname::AbstractString, resv::Integer)
    isempty(strip(path)) && throw(SpiceError(:blankfilename, "dafonw: the file name is blank"))
    isempty(strip(ftype, ' ')) && throw(SpiceError(:blankfiletype, "dafonw: the file type is blank"))
    ttype = String(strip(ftype, ' '))
    all(c -> 0x20 <= c <= 0x7e, codeunits(ttype)) || throw(SpiceError(:illegalcharacter,
        "dafonw: the file type contains nonprinting characters"))
    (0 <= nd <= 124) || throw(SpiceError(:dafinvalidparams, "dafonw: ND was $nd, should be in range [0,124]"))
    (2 <= ni <= 250) || throw(SpiceError(:dafinvalidparams, "dafonw: NI was $ni, should be in range [2,250]"))
    _daf_sumsize(nd, ni) <= 125 || throw(SpiceError(:dafinvalidparams,
        "dafonw: summary size was $(_daf_sumsize(nd, ni)), should not exceed 125"))
    resv < 0 && throw(SpiceError(:dafnoresv,
        "dafonw: an attempt was made to reserve a negative number ($resv) of records"))
    # Exclusive creation: the call fails if the path exists, so no file is ever truncated.
    io = try
        fd = Base.Filesystem.open(path, Base.Filesystem.JL_O_CREAT | Base.Filesystem.JL_O_EXCL |
            Base.Filesystem.JL_O_RDWR, 0o644)
        close(fd)
        open(path, "r+")
    catch e
        throw(SpiceError(:fileopenfailed,
            "dafonw: attempt to create the new file $path failed: $(sprint(showerror, e))"))
    end
    h = _daf_handle(path, io, true, false, Int(nd), Int(ni))
    fward = Int(resv) + 2
    free = _daf_rwa(fward + 2)
    filerec = zeros(UInt8, _DAF_RECLEN)
    filerec[1:8] = _daf_padded("DAF/" * String(_daf_padded(ttype, 4)), 8)
    filerec[9:16] = reinterpret(UInt8, Int32[nd, ni])
    filerec[17:76] = _daf_padded(ifname, 60)
    filerec[77:88] = reinterpret(UInt8, Int32[fward, fward, free])
    filerec[89:96] = codeunits(_DAF_NATIVE_BFF)
    filerec[700:727] = _DAF_FTPSTR
    _daf_write_at!(h, filerec, 0)
    # SPICELIB uses one character buffer for the comment records and for the first name
    # record: the name record is blank only if there is no reserved record, and else it
    # holds the nulls (and, for one reserved record, the end-of-comments mark) of the
    # last comment record that was written. The same bytes are written here.
    crec = resv > 0 ? zeros(UInt8, _DAF_CRECLEN) : fill(UInt8(' '), _DAF_CRECLEN)
    for recno in 2:resv+1                # invariant: comment records 2 .. recno-1 are written
        crec[1] = recno == 2 ? _DAF_EOC : 0x00
        _daf_write_crec!(h, recno, crec)
    end
    _daf_write_drec!(h, fward, zeros(Float64, _DAF_RECWORDS))
    _daf_write_crec!(h, fward + 1, crec)
    return h
end

"""
    dafopw(path) -> DAFHandle

Open the existing DAF file `path` for writing, as `dafopw_c`, to add arrays or comments.

The file must have the native byte order. Make sure that no other handle writes the same
file. Unload the file from every context before the call and load it again after the
close: a kernel that is loaded does not see arrays that are added later.

Throws the `SpiceError` kinds of [`dafopr`](@ref), `:unsupportedbff` if the file has the
other byte order, and `:fileopenfailed` if the file cannot be opened for update.
"""
function dafopw(path::AbstractString)
    nd, ni, swap = _daf_probe(path, "dafopw")
    swap && throw(SpiceError(:unsupportedbff,
        "dafopw: the file $path has a non-native binary file format; it can be read only"))
    io = try
        open(path, "r+")
    catch e
        throw(SpiceError(:fileopenfailed, "dafopw: cannot open $path for write: $(sprint(showerror, e))"))
    end
    return _daf_handle(path, io, true, false, nd, ni)
end

# ---------------------------------------------------------------------------------------
# Add a new array
# ---------------------------------------------------------------------------------------

"""
    dafbna(handle, sum, name) -> Nothing

Begin a new array in the DAF file of `handle`, as `dafbna_c`.

`sum` is the packed summary ([`dafps`](@ref)); [`dafena`](@ref) sets its two address
components. `name` is the array name, cut or padded to `8 * (nd + (ni + 1) ÷ 2)`
characters (40 for an SPK or a CK file). Add the data with [`dafada`](@ref) and end the
array with `dafena`. An array that was begun before and not ended is abandoned, as in
CSPICE.

Throws `SpiceError(:dafinvalidparams)` if `sum` has fewer doubles than a summary of the
file.

# Examples
```jldoctest
julia> h = dafonw(joinpath(mktempdir(), "example.daf"), "TEST", 2, 6, "example", 0);

julia> dafbna(h, dafps([0.0, 1.0], [1, 2, 3, 4, 0, 0]), "first array")

julia> dafada(h, [10.0, 20.0, 30.0])

julia> dafena(h)

julia> dafbfs(h)

julia> daffna(h)
true

julia> dafus(dafgs(h), 2, 6)
([0.0, 1.0], [1, 2, 3, 4, 385, 387])

julia> dafcls(h)
```
"""
function dafbna(h::DAFHandle, sum::AbstractVector{<:Real}, name::AbstractString)
    _daf_check_write(h, "dafbna")
    ss = _daf_sumsize(h.nd, h.ni)
    length(sum) >= ss || throw(SpiceError(:dafinvalidparams,
        "dafbna: the summary has $(length(sum)) doubles, the DAF needs $ss"))
    _, _, _, free = _daf_filerecord(h, "dafbna")
    h.asum = Float64[Float64(sum[j]) for j in 1:ss]
    h.aname = _daf_padded(name, 8ss)
    h.abegin = free
    h.afree = free
    h.adding = true
    return nothing
end

"""
    dafada(handle, data) -> Nothing

Append the doubles `data` to the array that [`dafbna`](@ref) began on `handle`, as
`dafada_c`. The data are written to the file at once.

Throws `SpiceError(:dafnewconflict)` if no array is begun on the handle, and
`:dafwritefail` if the array would pass the largest address of a DAF or the system reports
a write error.
"""
function dafada(h::DAFHandle, data::AbstractVector{<:Real})
    _daf_check_write(h, "dafada")
    h.adding || throw(SpiceError(:dafnewconflict,
        "dafada: an attempt was made to add data to an array that has not yet been " *
        "begun, in file $(h.path)"))
    n = length(data)
    n == 0 && return nothing
    h.afree + n - 1 <= typemax(Int32) - 2 * _DAF_RECWORDS || throw(SpiceError(:dafwritefail,
        "dafada: the array would pass the largest address of a DAF in $(h.path)"))
    _daf_write_data!(h, h.afree, data isa Vector{Float64} ? data : Float64[Float64(x) for x in data])
    h.afree += n
    return nothing
end

"""
    dafena(handle) -> Nothing

End the array that [`dafbna`](@ref) began on `handle`, as `dafena_c`.

The call writes the summary and the name of the array and the new free address of the
file. A search that is begun after the call finds the array. Before the call the file on
disk is a valid DAF without the array.

Throws `SpiceError(:dafnewconflict)` if no array is begun on the handle, and `:baddaf` if
the last summary record of the file is damaged.
"""
function dafena(h::DAFHandle)
    _daf_check_write(h, "dafena")
    h.adding || throw(SpiceError(:dafnewconflict,
        "dafena: an attempt was made to end an array that has not yet been begun, " *
        "in file $(h.path)"))
    nd, ni = h.nd, h.ni
    ss = _daf_sumsize(nd, ni)
    dc, ic = dafus(h.asum, nd, ni)
    ic[ni-1] = h.abegin
    ic[ni] = h.afree - 1
    sum = dafps(dc, ic)
    _, fward, bward, _ = _daf_filerecord(h, "dafena")
    sumrec, found = _daf_drec(h, bward)
    found || throw(SpiceError(:baddaf, "dafena: summary record $bward is not in $(h.path)"))
    namrec, _ = _daf_crec(h, bward + 1)
    narray = Int(sumrec[3])
    (0 <= narray < 125 ÷ ss) || throw(SpiceError(:baddaf,
        "dafena: summary record $bward of $(h.path) holds $narray summaries"))
    dloc = narray * ss + 4
    copyto!(sumrec, dloc, sum, 1, ss)
    namsiz = 8ss
    copyto!(namrec, narray * namsiz + 1, h.aname, 1, namsiz)
    sumrec[3] += 1.0
    narray += 1
    free = h.afree
    last = bward
    if narray < 125 ÷ ss
        _daf_write_drec!(h, bward, sumrec)
        _daf_write_crec!(h, bward + 1, namrec)
    else
        # The summary record is full: a new summary record and name record follow the data.
        next = _daf_arw(h.afree - 1)[1] + 1
        sumrec[1] = Float64(next)
        _daf_write_drec!(h, bward, sumrec)
        _daf_write_crec!(h, bward + 1, namrec)
        newrec = zeros(Float64, _DAF_RECWORDS)
        newrec[2] = Float64(bward)
        _daf_write_drec!(h, next, newrec)
        _daf_write_crec!(h, next + 1, fill(UInt8(' '), _DAF_CRECLEN))
        last = next
        free = _daf_rwa(next + 2)
    end
    _daf_write_filerecord!(h; fward=fward, bward=last, free=free)
    h.adding = false
    if h.searching && h.sthis == bward   # keep a search on this handle consistent
        copyto!(h.sbuf, sumrec)
        h.snsum = narray
        h.snext = Int(sumrec[1])
    end
    return nothing
end

function _daf_add_array!(h::DAFHandle, dc, ic, name::AbstractString, data::Vector{Float64})
    dafbna(h, dafps(dc, ic), name)
    dafada(h, data)
    dafena(h)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Change summaries and names
# ---------------------------------------------------------------------------------------

"""
    dafws(handle, sum) -> Nothing

Replace the packed summary of the current array of the search on `handle` with `sum`, as
`dafws_c`. The two address components are replaced too; [`dafrs`](@ref) keeps them.

Throws `SpiceError(:nocurrentarray)` if the search has no current array, and
`:dafinvalidparams` if `sum` has fewer doubles than a summary of the file.
"""
function dafws(h::DAFHandle, sum::AbstractVector{<:Real})
    _daf_check_write(h, "dafws")
    _daf_check_current(h, "dafws")
    ss = _daf_sumsize(h.nd, h.ni)
    length(sum) >= ss || throw(SpiceError(:dafinvalidparams,
        "dafws: the summary has $(length(sum)) doubles, the DAF needs $ss"))
    off = 3 + (h.scurr - 1) * ss
    for j in 1:ss
        h.sbuf[off+j] = Float64(sum[j])
    end
    _daf_write_drec!(h, h.sthis, copy(h.sbuf))
    return nothing
end

"""
    dafrs(handle, sum) -> Nothing

Replace the summary of the current array of the search on `handle` with `sum`, as
`dafrs_c`, but keep its two address components (the last two integers).

Throws `SpiceError(:nocurrentarray)` if the search has no current array, and
`:dafinvalidparams` if `sum` has fewer doubles than a summary of the file.
"""
function dafrs(h::DAFHandle, sum::AbstractVector{<:Real})
    _daf_check_write(h, "dafrs")
    _daf_check_current(h, "dafrs")
    _, exic = dafus(dafgs(h), h.nd, h.ni)
    newdc, newic = dafus(sum, h.nd, h.ni)
    newic[h.ni-1] = exic[h.ni-1]
    newic[h.ni] = exic[h.ni]
    dafws(h, dafps(newdc, newic))
    return nothing
end

"""
    dafrn(handle, name) -> Nothing

Replace the name of the current array of the search on `handle` with `name`, as
`dafrn_c`. `name` is cut or padded to the name size of the file.

Throws `SpiceError(:nocurrentarray)` if the search has no current array.
"""
function dafrn(h::DAFHandle, name::AbstractString)
    _daf_check_write(h, "dafrn")
    _daf_check_current(h, "dafrn")
    namsiz = 8 * _daf_sumsize(h.nd, h.ni)
    namrec, _ = _daf_crec(h, h.sthis + 1)
    copyto!(namrec, (h.scurr - 1) * namsiz + 1, _daf_padded(name, namsiz), 1, namsiz)
    _daf_write_crec!(h, h.sthis + 1, namrec)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Reserved records and comments
# ---------------------------------------------------------------------------------------

"Add `delta` to the two address components of every summary of `h`, as the end of `dafarr` and `dafrrr`."
function _daf_shift_addresses!(h::DAFHandle, delta::Int)
    dafbfs(h)
    while daffna(h)                      # invariant: every summary before the current one is shifted
        dc, ic = dafus(dafgs(h), h.nd, h.ni)
        ic[h.ni-1] += delta
        ic[h.ni] += delta
        dafws(h, dafps(dc, ic))
    end
    return nothing
end

"""
    _daf_add_reserved!(h, resv) -> Nothing

Add `resv` reserved records after the existing ones, as `dafarr`: every record from
the first summary record on moves down by `resv` records, the links of the summary
records, the pointers of the file record and the addresses in all summaries are
adjusted, and the new reserved records are filled with nulls. (Where SPICELIB writes
an uninitialized buffer to extend the file, this routine writes zeros.)
"""
function _daf_add_reserved!(h::DAFHandle, resv::Int)
    resv < 1 && return nothing
    _, fward, bward, free = _daf_filerecord(h, "dafarr")
    incr = resv * _DAF_RECWORDS
    lastrec = _daf_arw(free)[1]
    zero = zeros(Float64, _DAF_RECWORDS)
    for i in 1:resv
        _daf_write_drec!(h, lastrec + i, zero)
    end
    begblk = bward
    endblk = lastrec
    while begblk > 0
        # invariant: every record after `endblk` is at its new place
        for recno in endblk:-1:begblk+2
            drec, _ = _daf_drec(h, recno)
            _daf_write_drec!(h, recno + resv, drec)
        end
        crec, _ = _daf_crec(h, begblk + 1)
        _daf_write_crec!(h, begblk + 1 + resv, crec)
        drec, found = _daf_drec(h, begblk)
        found || throw(SpiceError(:baddaf, "dafarr: summary record $begblk is not in $(h.path)"))
        next = Int(drec[2])
        Int(drec[1]) > 0 && (drec[1] += resv)
        Int(drec[2]) > 0 && (drec[2] += resv)
        _daf_write_drec!(h, begblk + resv, drec)
        endblk = begblk - 1
        begblk = next
    end
    fward += resv
    bward += resv
    free += incr
    _daf_write_filerecord!(h; fward=fward, bward=bward, free=free)
    _daf_shift_addresses!(h, incr)
    nulls = zeros(UInt8, _DAF_CRECLEN)
    for recno in fward-resv:fward-1
        _daf_write_crec!(h, recno, nulls)
    end
    return nothing
end

"""
    dafac(handle, lines) -> Nothing

Append the comment lines `lines` to the comment area of the DAF file, as `dafac_c`.

Trailing blanks of a line are not stored. If the reserved records do not have room,
records are added and the data of the file move. Every address that was read from a
summary before the call is then stale, so read the summaries again.

Throws `SpiceError`:
- `:invalidargument` if `lines` is empty;
- `:illegalcharacter` if a line has a character outside the printable ASCII codes 32 to
  126;
- `:dafnewconflict` if an array is begun and not ended on the handle.

# Examples
```jldoctest
julia> h = dafonw(joinpath(mktempdir(), "example.daf"), "TEST", 2, 6, "example", 0);

julia> dafac(h, ["first line", "second line   "])

julia> dafec(h)
2-element Vector{String}:
 "first line"
 "second line"

julia> dafcls(h)
```

See also [`dafec`](@ref), [`dafdc`](@ref).
"""
function dafac(h::DAFHandle, lines::AbstractVector{<:AbstractString})
    _daf_check_write(h, "dafac")
    n = length(lines)
    n <= 0 && throw(SpiceError(:invalidargument, "dafac: the number of comment lines " *
        "to be added to the binary DAF file $(h.path) was not positive: $n"))
    h.adding && throw(SpiceError(:dafnewconflict,
        "dafac: an array is begun and not ended in $(h.path); call dafena first"))
    texts = Vector{Vector{UInt8}}(undef, n)
    nchars = 0
    for (i, line) in enumerate(lines)
        s = String(line)
        cu = codeunits(s)[1:_daf_lastnb(s)]
        bad = findfirst(c -> c > 0x7e || c < 0x20, cu)
        bad === nothing || throw(SpiceError(:illegalcharacter, "dafac: a nonprinting " *
            "character was encountered in the comment buffer. Value: $(Int(cu[bad]))"))
        texts[i] = cu
        nchars += length(cu)
    end
    nchars += n + 1
    _, fward, _, _ = _daf_filerecord(h, "dafac")
    ncomc, ncomr, notused = _daf_comment_extent(h, fward, "dafac")
    space = (ncomr + notused) * _DAF_CRECLEN - ncomc
    newrec = nchars > space ? (nchars - space - 1) ÷ _DAF_CRECLEN + 1 : 0
    _daf_add_reserved!(h, newrec)
    crec = fill(UInt8(' '), _DAF_CRECLEN)
    recno = 2
    curpos = 1
    if ncomc != 0
        rinuse = ncomc ÷ _DAF_CRECLEN + 1
        recno = rinuse + 1
        curpos = ncomc - (rinuse - 1) * _DAF_CRECLEN + 1
        crec, _ = _daf_crec(h, recno)
    end
    for text in texts
        # invariant: all earlier lines are in records 2 .. recno-1 and crec[1:curpos-1]
        for c in text
            if curpos > _DAF_CRECLEN
                _daf_write_crec!(h, recno, crec)
                recno += 1
                curpos = 1
                fill!(crec, UInt8(' '))
            end
            crec[curpos] = c
            curpos += 1
        end
        if curpos > _DAF_CRECLEN
            # As in SPICELIB, the buffer is not blanked here: the next record keeps the
            # characters of this one after the position that is written.
            _daf_write_crec!(h, recno, crec)
            recno += 1
            curpos = 1
        end
        crec[curpos] = _DAF_EOL
        curpos += 1
    end
    if curpos > _DAF_CRECLEN
        _daf_write_crec!(h, recno, crec)
        recno += 1
        curpos = 1
        fill!(crec, UInt8(' '))
    end
    crec[curpos] = _DAF_EOC
    _daf_write_crec!(h, recno, crec)
    return nothing
end

"""
    dafdc(handle) -> Nothing

Delete the comment area of the DAF file, as `dafdc_c`.

All reserved records are removed and the data of the file move up. Every address that was
read from a summary before the call is then stale. The file does not get shorter on disk,
as in CSPICE.

Throws `SpiceError(:dafnewconflict)` if an array is begun and not ended on the handle.
"""
function dafdc(h::DAFHandle)
    _daf_check_write(h, "dafdc")
    h.adding && throw(SpiceError(:dafnewconflict,
        "dafdc: an array is begun and not ended in $(h.path); call dafena first"))
    _, fward, bward, free = _daf_filerecord(h, "dafdc")
    remove = fward - 2
    remove < 1 && return nothing
    decr = remove * _DAF_RECWORDS
    begblk = fward
    while begblk > 0
        # invariant: every record before `begblk` is at its new place
        drec, found = _daf_drec(h, begblk)
        found || throw(SpiceError(:baddaf, "dafdc: summary record $begblk is not in $(h.path)"))
        next = Int(drec[1]) > 0 ? Int(drec[1]) : 0
        endblk = next > 0 ? next - 1 : _daf_arw(free)[1]
        Int(drec[1]) > 0 && (drec[1] -= remove)
        Int(drec[2]) > 0 && (drec[2] -= remove)
        _daf_write_drec!(h, begblk - remove, drec)
        crec, _ = _daf_crec(h, begblk + 1)
        _daf_write_crec!(h, begblk + 1 - remove, crec)
        for recno in begblk+2:endblk
            data, there = _daf_drec(h, recno)
            there && _daf_write_drec!(h, recno - remove, data)
        end
        begblk = next
    end
    _daf_write_filerecord!(h; fward=fward - remove, bward=bward - remove, free=free - decr)
    _daf_shift_addresses!(h, -decr)
    return nothing
end
