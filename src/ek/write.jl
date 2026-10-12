# Writing an EK file, part 3: the public calls.
#
# `ekopn` and `ekopw` give a write handle (a `DASHandle`); `ekbseg` starts a segment,
# `ekappr`, `ekinsr` and `ekdelr` add and delete records, `ekacec`, `ekaced`, `ekacei`
# add column entries and `ekucec`, `ekuced`, `ekucei` replace them. `ekifld`, `ekaclc`,
# `ekacld`, `ekacli` and `ekffld` are the fast writers, which store a segment column by
# column. `ekcls` (in `api.jl`) closes the file.
#
# Segment, record and element numbers start at 1 (the C interface of CSPICE counts from
# 0; the Fortran interface, which this file ports, counts from 1).
#
# CSPICE keeps one "shadow" mechanism for the backup of changed records (`ekshdw`,
# `zzekrbck`); in CSPICE it is not active (`ekshdw` always returns false), so its
# branches are not ported.

# ---------------------------------------------------------------------------------------
# Opening a file for write
# ---------------------------------------------------------------------------------------

"""
    ekopn(path, ifname, ncomch=0) -> DASHandle

Create the new EK file `path` with the internal file name `ifname` and room for
`ncomch` comment characters, and open it for write, as `ekopn_c`. Close the handle
with [`ekcls`](@ref); the data are on disk when `ekcls` returns.

Preconditions (each throws `SpiceError`): `ncomch >= 0` (`:invalidcount`); the names
are not empty (`:emptystring`); the file does not exist (`:fileopenfailed`), and the
other conditions of `dasonw`. Postcondition: the file is a paged EK with no segment.
"""
function ekopn(path::AbstractString, ifname::AbstractString, ncomch::Integer=0)
    (isempty(path) || isempty(ifname)) && throw(SpiceError(:emptystring,
        "ekopn: the file name and the internal file name must not be empty"))
    ncomch < 0 && throw(SpiceError(:invalidcount, "The number of reserved comment " *
        "characters must be non-negative but was $ncomch."))
    h = dasonw(path, "EK", ifname, (Int(ncomch) + 1023) ÷ 1024)
    _ekw_pgin(h)
    p, base = _ekw_pgan(h, 3)
    tree = _ekw_trit(h)
    h.ints[base+1] = tree
    return h
end

"""
    ekopw(path) -> DASHandle

Open the existing EK file `path` for write, as `ekopw_c`. The file is read into
memory; all changes are on disk when [`ekcls`](@ref) returns. Contract: `path` names a
DAS file that holds a paged EK of architecture 8, else `SpiceError` with the kinds of
`dasopw`, `:wrongarchitecture` or `:invalidformat`.
"""
function ekopw(path::AbstractString)
    isempty(path) && throw(SpiceError(:emptystring, "ekopw: the file name must not be empty"))
    h = dasopw(path)
    _ek_check(h, h.path, "ekopw")
    return h
end

# ---------------------------------------------------------------------------------------
# Column declarations
# ---------------------------------------------------------------------------------------

"A string as the bytes of a Fortran string: without trailing blanks."
@inline _ekw_bytes(s::AbstractString) = Vector{UInt8}(codeunits(rstrip(String(s), ' ')))
@inline _ekw_bytes(s::AbstractVector{UInt8}) = Vector{UInt8}(s)

"""
    _ekw_chckid(class, maxlen, id) -> Nothing

The check of `chckid`: the last non-blank character of `id` is at a position not after
`maxlen` (`:idstringtoolong`) and `id` has printable characters only
(`:nonprintablechars`).
"""
function _ekw_chckid(class::AbstractString, maxlen::Int, id::AbstractString)
    b = codeunits(id)
    l = something(findlast(!=(UInt8(' ')), b), 0)
    l > maxlen && throw(SpiceError(:idstringtoolong, "The $class '$id' is invalid; the last " *
        "non-blank character is located at position $l; the maximum allowed length is $maxlen."))
    i = findfirst(c -> c < 0x20 || c > 0x7e, b)
    i === nothing || throw(SpiceError(:nonprintablechars, "The $class is invalid; this " *
        "string contains a non-printing character (ICHAR = $(b[i])) at position $i."))
    return nothing
end

"""
    _ekw_isident(name) -> Bool

True if `name` (trailing blanks do not count) is an identifier as `lxidnt` finds it
with the default character sets: a letter, then letters, digits, `\$` and `_`.
"""
function _ekw_isident(name::AbstractString)
    b = codeunits(name)
    l = something(findlast(!=(UInt8(' ')), b), 0)
    l == 0 && return false
    letter(c) = UInt8('A') <= c <= UInt8('Z') || UInt8('a') <= c <= UInt8('z')
    letter(b[1]) || return false
    for k in 2:l                         # invariant: the characters before k are valid
        c = b[k]
        (letter(c) || UInt8('0') <= c <= UInt8('9') || c == UInt8('$') || c == UInt8('_')) ||
            return false
    end
    return true
end

@noinline _ekw_baddecl(msg::AbstractString, decl::AbstractString) =
    throw(SpiceError(:badcolumndecl, msg * " Declaration was: " * decl))

"The integer that `nparsi` reads from `token`, or `nothing` if the token is not a number."
function _ekw_nparsi(token::AbstractString)
    x = try
        _nparsd(token)
    catch e
        e isa SpiceError || rethrow()
        return nothing
    end
    (typemin(Int32) <= x <= typemax(Int32)) || return nothing
    return trunc(Int32, x)
end

"""
    _ekw_lparsm(list, delims, nmax) -> Vector{String}

Items of `list` that the characters of `delims` (no blank among them) separate, at
most `nmax`, as `lparsm`: blanks around an item are dropped, two delimiters in a row
give a blank item, and a list that ends with a delimiter has a last blank item. A
blank list has one blank item.
"""
function _ekw_lparsm(list::AbstractString, delims::AbstractString, nmax::Int)
    c = codeunits(rstrip(String(list), ' '))
    items = String[]
    eol = length(c)
    eol == 0 && return push!(items, "")
    isdelim(k) = Char(c[k]) in delims
    b = 1
    while b <= eol
        # invariant: `items` holds the items that end before position b
        while b <= eol && c[b] == UInt8(' ')
            b += 1
        end
        e = b
        while e <= eol && !isdelim(e)
            e += 1
        end
        push!(items, String(rstrip(String(c[b:e-1]), ' ')))
        length(items) < nmax || return items
        b = e + 1
    end
    isdelim(eol) && length(items) < nmax && push!(items, "")
    return items
end

"""
    _ekw_pdec(decl) -> Vector{Int32}

The column descriptor (11 integers) that the declaration `decl` gives, as `zzekpdec`.
A declaration is a list of `KEYWORD = value` items that commas separate, in any case:

- `DATATYPE` (required): `INTEGER`, `DOUBLE PRECISION`, `TIME`, `CHARACTER*(n)` or
  `CHARACTER*(*)`;
- `SIZE`: a number or `VARIABLE` (1 if absent);
- `INDEXED`, `NULLS_OK`, `FIXED_COUNT`: `TRUE` or `FALSE` (false if absent).

Postcondition: class, type, string length, size, index type (1 or -1) and null flag
(1 or -1) are set; name base, index pointer, position and metadata pointer are 0.
Throws `SpiceError(:badcolumndecl)` for a declaration that does not follow these
rules, and `:badcolumdecl` (the spelling of CSPICE) if `DATATYPE` is absent.
"""
function _ekw_pdec(decl::AbstractString)
    pardsc = zeros(Int32, _EK_CDSCSZ)
    tokens = String[uppercase(String(first(t, 32))) for t in _ekw_lparsm(decl, ",=", 20)]
    n = length(tokens)
    keys = ("DATATYPE", "SIZE", "INDEXED", "NULLS_OK", "FIXED_COUNT")
    attloc = zeros(Int, 5)
    for i in 1:n                         # invariant: `attloc` holds the last keyword places before i
        j = findfirst(==(tokens[i]), keys)
        j === nothing || (attloc[j] = i)
    end
    attloc[1] == 0 && throw(SpiceError(:badcolumdecl,
        "Required keyword DATATYPE was not found in column declaration $decl."))
    tokloc = attloc[1] + 1
    n < tokloc && _ekw_baddecl("Column data type specification did not follow \"DATATYPE\" " *
        "keyword.", decl)
    tok = tokens[tokloc]
    if tok == "INTEGER"
        pardsc[_EK_TYPIDX] = _EK_INT
        pardsc[_EK_LENIDX] = 1
    elseif _ek_eqstr(tok, "DOUBLE PRECISION")
        pardsc[_EK_TYPIDX] = _EK_DP
        pardsc[_EK_LENIDX] = 1
    elseif _ek_eqstr(tok, "TIME")
        pardsc[_EK_TYPIDX] = _EK_TIME
        pardsc[_EK_LENIDX] = 1
    elseif startswith(tok, "CHARACTER")
        pardsc[_EK_TYPIDX] = _EK_CHR
        tok = filter(c -> c != ' ' && c != '(' && c != ')', tok)
        b = codeunits(tok)
        (length(b) >= 10 && b[10] == UInt8('*')) || _ekw_baddecl("Required asterisk missing " *
            "from character column declaration:  $tok.", decl)
        if length(b) >= 11 && b[11] == UInt8('*')
            pardsc[_EK_LENIDX] = -1
        else
            len = _ekw_nparsi(String(b[11:end]))
            len === nothing && _ekw_baddecl("String length specification $tok didn't parse " *
                "as an integer.", decl)
            pardsc[_EK_LENIDX] = len
        end
    else
        _ekw_baddecl("Data type specification $tok is unrecognized.", decl)
    end
    if attloc[2] != 0
        tokloc = attloc[2] + 1
        n < tokloc && _ekw_baddecl("Column size specification did not follow \"SIZE\" " *
            "keyword.", decl)
        tok = tokens[tokloc]
        if tok == "VARIABLE"
            (pardsc[_EK_TYPIDX] == _EK_CHR && pardsc[_EK_LENIDX] == -1) && _ekw_baddecl(
                "Column size specification was VARIABLE for a CHARACTER*(*) column.", decl)
            pardsc[_EK_SIZIDX] = -1
        else
            size = _ekw_nparsi(tok)
            size === nothing && _ekw_baddecl("Column element size specification $tok " *
                "didn't parse as an integer.", decl)
            pardsc[_EK_SIZIDX] = size
        end
    else
        pardsc[_EK_SIZIDX] = 1
    end
    scalar = pardsc[_EK_SIZIDX] == 1
    dtype = pardsc[_EK_TYPIDX]
    pardsc[_EK_CLSIDX] = dtype == _EK_CHR ? (scalar ? 3 : 6) : dtype == _EK_INT ?
                         (scalar ? 1 : 4) : (scalar ? 2 : 5)
    # `flag(j, name)` reads the boolean that follows keyword j: true, false, or absent.
    function flag(j::Int, name::AbstractString)
        attloc[j] == 0 && return nothing
        loc = attloc[j] + 1
        n < loc && _ekw_baddecl("Boolean value did not follow \"$name\" keyword.", decl)
        tokens[loc] == "TRUE" && return true
        tokens[loc] == "FALSE" && return false
        _ekw_baddecl("Invalid token $(tokens[loc]) follows $name keyword.", decl)
    end
    pardsc[_EK_NFLIDX] = flag(4, "NULLS_OK") === true ? 1 : -1
    if flag(3, "INDEXED") === true
        scalar || _ekw_baddecl("Non-scalar columns cannot be indexed.", decl)
        pardsc[_EK_IXTIDX] = 1
    else
        pardsc[_EK_IXTIDX] = -1
    end
    if flag(5, "FIXED_COUNT") === true
        pardsc[_EK_CLSIDX] <= 3 || _ekw_baddecl("FIXED_COUNT attribute used in non-scalar " *
            "column declaration.", decl)
        pardsc[_EK_CLSIDX] += 6
    end
    return pardsc
end

"""
    _ekw_styp(cdscrs) -> Int

Segment type of a set of column descriptors (`zzekstyp`): 1 if all classes are 1 to 6,
2 if all are 7 to 9. Throws `SpiceError(:badattributes)` for a mixture.
"""
function _ekw_styp(cdscrs::Vector{Vector{Int32}})
    var = any(d -> 1 <= d[_EK_CLSIDX] <= 6, cdscrs)
    fixed = any(d -> 7 <= d[_EK_CLSIDX] <= 9, cdscrs)
    (var && !fixed) && return 1
    (fixed && !var) && return 2
    throw(SpiceError(:badattributes, "Column set contains a mixture of variable and " *
        "fixed-count columns.  Segments must contain all variable or all fixed count columns."))
end

# ---------------------------------------------------------------------------------------
# Segments
# ---------------------------------------------------------------------------------------

"""
    _ekw_bs(h, stype, tabnam, cnames, cdscrs) -> segno

Write the metadata of a new, empty segment of type `stype` (`zzekbs01`, `zzekbs02`):
adjacent integer pages with the segment descriptor and the column descriptors,
adjacent character pages with the table name and the column names (upper case), the
record tree (type 1 only), the three data page trees and, in a type 1 segment, an
empty index tree for every indexed column. The segment is appended to the segment tree
of the file. Mutates `cdscrs` (position, name base, index). Postcondition: `segno` is
the number of segments of the file.
"""
function _ekw_bs(h::DASHandle, stype::Int, tabnam::AbstractString,
        cnames::AbstractVector{<:AbstractString}, cdscrs::Vector{Vector{Int32}})
    ncols = length(cdscrs)
    nipage = (ncols * _EK_CDSCSZ + _EK_SDSCSZ + 255) ÷ 256
    p1, p1base = _ekw_pgan(h, 3)
    for _ in 2:nipage
        _ekw_pgan(h, 3)
    end
    ncpage = (ncols * _EK_CNAMSZ + 1087) ÷ 1024
    _, cbase = _ekw_pgan(h, 1)
    for _ in 2:ncpage
        _ekw_pgan(h, 1)
    end
    rec = stype == 1 ? _ekw_trit(h) : 0
    cpt = _ekw_trit(h)
    dpt = _ekw_trit(h)
    ipt = _ekw_trit(h)
    segdsc = zeros(Int32, _EK_SDSCSZ)
    segdsc[_EK_EKTIDX] = stype
    segdsc[_EK_SNOIDX] = _ek_nseg(h) + 1
    segdsc[_EK_IMDIDX] = p1base
    segdsc[_EK_TNMIDX] = cbase
    segdsc[_EK_NCIDX] = ncols
    segdsc[_EK_RTIDX] = rec
    segdsc[_EK_CPTIDX] = cpt
    segdsc[_EK_DPTIDX] = dpt
    segdsc[_EK_IPTIDX] = ipt
    segdsc[_EK_MFLIDX] = 1
    segdsc[_EK_IFLIDX] = -1
    segdsc[_EK_SHDIDX] = -1
    segdsc[_EK_LCWIDX] = _EK_CPSIZE
    segdsc[_EK_LCWIDX+1] = _EK_DPSIZE
    segdsc[_EK_LCWIDX+2] = _EK_IPSIZE
    segdsc[_EK_NMLIDX] = cbase + _EK_TNAMSZ
    tname = Vector{UInt8}(codeunits(uppercase(String(lstrip(tabnam, ' ')))))
    _ekw_putc!(h.chars, cbase + 1, tname, 1, _EK_TNAMSZ)
    for i in 1:ncols                     # invariant: the columns before i are in the file
        coldsc = cdscrs[i]
        coldsc[_EK_ORDIDX] = i
        nambas = cbase + _EK_TNAMSZ + (i - 1) * _EK_CNAMSZ
        coldsc[_EK_NAMIDX] = nambas
        cname = Vector{UInt8}(codeunits(uppercase(String(cnames[i]))))
        _ekw_putc!(h.chars, nambas + 1, cname, 1, _EK_CNAMSZ)
        if stype == 1 && coldsc[_EK_IXTIDX] != -1
            coldsc[_EK_IXTIDX] = 1       # `zzekcix1`
            coldsc[_EK_IXPIDX] = _ekw_trit(h)
        end
        copyto!(h.ints, p1base + _EK_SDSCSZ + (i - 1) * _EK_CDSCSZ + 1, coldsc, 1, _EK_CDSCSZ)
    end
    copyto!(h.ints, p1base + 1, segdsc, 1, _EK_SDSCSZ)
    return _ekw_trap(h, _ek_segtree(h), p1)
end

"""
    ekbseg(handle, tabnam, cnames, decls) -> segno

Start a new segment of the table `tabnam` in the EK file of the write handle, as
`ekbseg_c`, and return its number (from 1). `cnames` holds the column names and
`decls` the column declarations, one for each column: a list of `KEYWORD = value`
items with the keywords `DATATYPE` (`INTEGER`, `DOUBLE PRECISION`, `TIME`,
`CHARACTER*(n)`, `CHARACTER*(*)`), `SIZE` (a number or `VARIABLE`), `INDEXED`,
`NULLS_OK` and `FIXED_COUNT` (`TRUE` or `FALSE`).

Preconditions (each throws `SpiceError`): write handle (`:dasnosuchhandle`,
`:dasinvalidaccess`); names of at most 64 (table) and 32 (column) printable characters
(`:emptystring`, `:idstringtoolong`, `:nonprintablechars`) that start with a letter
and hold letters, digits, `\$` and `_` (`:invalidname`); 1 to 100 columns and as many
declarations (`:invalidcount`); valid declarations (`:badcolumndecl`, `:badcolumdecl`
for a missing data type); no mixture of fixed-count and other columns
(`:badattributes`). Postcondition: the segment exists and has no rows.
"""
function ekbseg(h::DASHandle, tabnam::AbstractString, cnames::AbstractVector{<:AbstractString},
        decls::AbstractVector{<:AbstractString})
    _ekw_check(h, "ekbseg")
    isempty(tabnam) && throw(SpiceError(:emptystring, "ekbseg: the table name is empty"))
    _ekw_chckid("EK table name", _EK_TNAMSZ, tabnam)
    _ekw_isident(tabnam) || throw(SpiceError(:invalidname,
        "Table name <$tabnam> violates syntax rules."))
    ncols = length(cnames)
    (1 <= ncols <= _EK_MXCLSG) || throw(SpiceError(:invalidcount,
        "Number of columns must be in range 1:$_EK_MXCLSG, was $ncols."))
    length(decls) == ncols || throw(SpiceError(:invalidcount,
        "ekbseg: $ncols column names but $(length(decls)) declarations"))
    for name in cnames
        _ekw_chckid("EK column name", _EK_CNAMSZ, name)
        _ekw_isident(name) || throw(SpiceError(:invalidname,
            "Column name <$name> violates syntax rules."))
    end
    cdscrs = Vector{Int32}[_ekw_pdec(d) for d in decls]
    stype = _ekw_styp(cdscrs)
    return _ekw_bs(h, stype, tabnam, cnames, cdscrs)
end

"""
    _ekw_segment(h, who, segno) -> (segdsc, mbase)

A copy of the descriptor of segment `segno` for a routine that changes records, and
the integer base of the segment metadata. Preconditions: write handle; valid segment
number (`:invalidindex`); segment of type 1, else `SpiceError(:invalidtype)` (CSPICE
has no such check and damages a type 2 segment).
"""
function _ekw_segment(h::DASHandle, who::AbstractString, segno::Integer)
    _ekw_check(h, who)
    _, mbase = _ek_mloc(h, segno)
    segdsc = h.ints[mbase+1:mbase+_EK_SDSCSZ]
    segdsc[_EK_EKTIDX] == 1 || throw(SpiceError(:invalidtype, "$who: segment $segno has " *
        "type $(segdsc[_EK_EKTIDX]); records can be changed in segments of type 1 only."))
    return segdsc, mbase
end

"""
    ekinsr(handle, segno, recno) -> Nothing

Add an empty record to segment `segno` so that it becomes row `recno`, as `ekinsr_c`;
the rows from `recno` on move up by one. Fill the record with [`ekacec`](@ref),
[`ekaced`](@ref), [`ekacei`](@ref): a segment that holds a record with an empty entry
cannot be read by all of CSPICE.

Preconditions (each throws `SpiceError`): write handle; valid segment number
(`:invalidindex`) of a segment that the fast writers did not make with fixed-count
columns (`:invalidtype`); `1 <= recno <= nrows + 1` (`:invalidindex`). Postcondition:
the segment has one row more.
"""
function ekinsr(h::DASHandle, segno::Integer, recno::Integer)
    segdsc, _ = _ekw_segment(h, "ekinsr", segno)
    _ekw_insr(h, segdsc, recno)
    return nothing
end

"The body of `ekinsr` for the descriptor copy `segdsc`."
function _ekw_insr(h::DASHandle, segdsc::Vector{Int32}, recno::Integer)
    ints = h.ints
    size = Int(segdsc[_EK_NCIDX]) + _EK_DPTBAS
    size > _EK_IPSIZE && _ekw_bug("Record pointer requires $size integer words; EK software " *
        "assumes size is <= $_EK_IPSIZE.")
    nrec = Int(segdsc[_EK_NRIDX])
    (1 <= recno <= nrec + 1) || throw(SpiceError(:invalidindex,
        "Record number = $recno; valid range is 1:$(nrec + 1)."))
    lastp = Int(segdsc[_ekw_lcp(3)])
    lastw = Int(segdsc[_ekw_lcw(3)])
    local recbas::Int
    if size <= _EK_IPSIZE - lastw
        base = _ek_ibase(lastp)
        recbas = base + lastw
        ints[base+256] += 1
        segdsc[_ekw_lcw(3)] += size
    else
        p, recbas = _ekw_aps(h, segdsc, 3)
        ints[recbas+256] = 1
        segdsc[_ekw_lcp(3)] = p
        segdsc[_ekw_lcw(3)] = size
    end
    ints[recbas+1] = 1
    fill!(view(ints, recbas+2:recbas+size), Int32(_EK_UNINIT))
    _ekw_trin(h, segdsc[_EK_RTIDX], recno, recbas)
    segdsc[_EK_NRIDX] += 1
    _ekw_put_segdsc(h, segdsc)
    return nothing
end

"""
    ekappr(handle, segno) -> recno

Add an empty record at the end of segment `segno`, as `ekappr_c`, and return its row
number (from 1). Contract as [`ekinsr`](@ref).
"""
function ekappr(h::DASHandle, segno::Integer)
    segdsc, _ = _ekw_segment(h, "ekappr", segno)
    recno = Int(segdsc[_EK_NRIDX]) + 1
    _ekw_insr(h, segdsc, recno)
    return recno
end

"""
    ekdelr(handle, segno, recno) -> Nothing

Delete row `recno` of segment `segno`, as `ekdelr_c`: its entries leave the indexes
and the data pages, and the rows after it move down by one. Preconditions (each throws
`SpiceError`): write handle; valid segment number (`:invalidindex`) of a type 1
segment (`:invalidtype`); `1 <= recno <= nrows` (`:invalidindex`). Postcondition: the
segment has one row less. Pages that become empty go to the free list of the file.

Rows of a segment that the fast writers made can be deleted too. CSPICE cannot do
that: it frees the record page of the row with all the records in it, it does not find
the row in the index of a character column with equal strings, and it takes links
from data pages that other array entries use. The record pages of a fast load are
therefore never freed here.
"""
function ekdelr(h::DASHandle, segno::Integer, recno::Integer)
    segdsc, mbase = _ekw_segment(h, "ekdelr", segno)
    ncols = Int(segdsc[_EK_NCIDX])
    nrec = Int(segdsc[_EK_NRIDX])
    (1 <= recno <= nrec) || throw(SpiceError(:invalidindex,
        "Record number = $recno; valid range is 1:$nrec."))
    recptr = _ek_trdp(h.ints, segdsc[_EK_RTIDX], recno)
    for i in 1:ncols                     # invariant: the entries of the columns before i are deleted
        dscbas = mbase + _EK_SDSCSZ + (i - 1) * _EK_CDSCSZ
        _ekw_de(h, segdsc, h.ints[dscbas+1:dscbas+_EK_CDSCSZ], recptr)
    end
    # A record page of the fast writers has no links: CSPICE frees it with the first
    # record that is deleted, and so loses the other records of the page. Here such a
    # page stays in the segment.
    fastrec = _ekw_fast_record(h, recptr)
    p, _ = _ekw_pgpg(3, recptr + 1)
    fastrec || _ekw_unlink(h, segdsc, 3, p)
    _ekw_trdl(h, segdsc[_EK_RTIDX], recno)
    segdsc[_EK_NRIDX] -= 1
    _ekw_put_segdsc(h, segdsc)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Column entries
# ---------------------------------------------------------------------------------------

"The 11 integers of the descriptor of the column `column` of a segment (`zzekcdsc`); `:bug` if unknown."
function _ekw_cdsc(h::DASHandle, segdsc::AbstractVector{Int32}, column::AbstractString)
    col = _ek_find_column(h, segdsc, column)
    b = Int(segdsc[_EK_IMDIDX]) + _EK_SDSCSZ + (Int(col.ord) - 1) * _EK_CDSCSZ
    return h.ints[b+1:b+_EK_CDSCSZ]
end

"Name of the data types `types` for an error message."
_ekw_typename(types::Tuple) = types[1] == _EK_CHR ? "character" : types[1] == _EK_INT ?
                              "integer" : "d.p. or time"

"""
    _ekw_locate(h, who, segno, recno, column, types) -> (segdsc, coldsc, recptr)

The common start of the calls that add or replace one entry: descriptor copies of the
segment and of the column, and the record pointer of row `recno`. Throws `SpiceError`:
`:dasnosuchhandle`; `:emptystring`; `:invalidindex` for a bad segment number; `:bug`
for an unknown column (as CSPICE); `:wrongdatatype` unless the data type code of the
column is in `types`; `:indexoutofrange` for a row that does not exist;
`:dasinvalidaccess` for a read handle.
"""
function _ekw_locate(h::DASHandle, who::AbstractString, segno::Integer, recno::Integer,
        column::AbstractString, types::Tuple)
    _ek_check_read(h, who)
    isempty(column) && throw(SpiceError(:emptystring, "$who: the column name is empty"))
    segdsc = Vector{Int32}(_ek_sdsc(h, segno))
    coldsc = _ekw_cdsc(h, segdsc, column)
    coldsc[_EK_TYPIDX] in types || throw(SpiceError(:wrongdatatype, "Column $column is of " *
        "type $(coldsc[_EK_TYPIDX]); $(uppercase(who)) only works with " *
        "$(_ekw_typename(types)) columns.  RECNO = $recno; SEGNO = $segno; EK = $(h.path)."))
    recptr = _ek_trdp(h.ints, segdsc[_EK_RTIDX], recno)
    _das_check_write(h, who)
    return segdsc, coldsc, recptr
end

@noinline _ekw_noclass(h::DASHandle, coldsc, what, column, recno, segno) = throw(SpiceError(
    :noclass, "Class $(coldsc[_EK_CLSIDX]) from input column descriptor is not a supported " *
    "$what class.  COLUMN = $column; RECNO = $recno; SEGNO = $segno; EK = $(h.path)."))

"The strings `vals` as byte vectors without trailing blanks."
_ekw_cvals(vals::AbstractVector) = Vector{UInt8}[_ekw_bytes(v) for v in vals]
"The first of the strings `vals` (a blank string if there is none) as bytes."
_ekw_cval1(vals::AbstractVector) = isempty(vals) ? UInt8[] : _ekw_bytes(first(vals))
"The first of the numbers `vals`, or zero if there is none."
_ekw_nval1(::Type{T}, vals::AbstractVector) where {T} = isempty(vals) ? zero(T) :
                                                         convert(T, first(vals))
"The numbers `vals` as a vector of `T`."
_ekw_nvals(::Type{T}, vals::AbstractVector) where {T} = vals isa Vector{T} ? vals : T[v for v in vals]

"""
    _ekw_entry_c(h, who, update, segno, recno, column, vals, isnull) -> Nothing

The body of `ekacec` (`update` false) and `ekucec` (`update` true).
"""
function _ekw_entry_c(h::DASHandle, who::AbstractString, update::Bool, segno::Integer,
        recno::Integer, column::AbstractString, vals::AbstractVector, isnull::Bool)
    segdsc, coldsc, recptr = _ekw_locate(h, who, segno, recno, column, (_EK_CHR,))
    class = coldsc[_EK_CLSIDX]
    if class == 3
        update && _ekw_de03(h, segdsc, coldsc, recptr)
        _ekw_ad03(h, segdsc, coldsc, recptr, _ekw_cval1(vals), isnull)
    elseif class == 6
        update && _ekw_de456(h, segdsc, coldsc, recptr)
        _ekw_ad06(h, segdsc, coldsc, recptr, _ekw_cvals(vals), isnull)
    else
        _ekw_noclass(h, coldsc, "character", column, recno, segno)
    end
    return nothing
end

"""
    _ekw_entry_n(h, who, update, T, types, segno, recno, column, vals, isnull) -> Nothing

The body of `ekaced`, `ekacei` (`update` false) and `ekuced`, `ekucei` (`update`
true); `T` is the element type of the column and `types` its data type codes.
"""
function _ekw_entry_n(h::DASHandle, who::AbstractString, update::Bool, ::Type{T}, types::Tuple,
        segno::Integer, recno::Integer, column::AbstractString, vals::AbstractVector,
        isnull::Bool) where {T}
    segdsc, coldsc, recptr = _ekw_locate(h, who, segno, recno, column, types)
    class = coldsc[_EK_CLSIDX]
    if class == 1 || class == 2
        val = _ekw_nval1(T, vals)
        if update
            _ekw_ue12(h, segdsc, coldsc, recptr, val, isnull)
        else
            _ekw_ad12(h, segdsc, coldsc, recptr, val, isnull)
        end
    elseif class == 4 || class == 5
        nvals = _ekw_nvals(T, vals)
        update && _ekw_de456(h, segdsc, coldsc, recptr)
        _ekw_ad45(h, segdsc, coldsc, recptr, nvals, isnull)
    else
        _ekw_noclass(h, coldsc, T === Int32 ? "integer" : "d.p.", column, recno, segno)
    end
    return nothing
end

"""
    ekacec(handle, segno, recno, column, cvals, isnull=false) -> Nothing

Add the entry of the character column `column` to row `recno` of segment `segno`, as
`ekacec_c`. `cvals` holds the strings of the entry (one for a scalar column); trailing
blanks do not count, and a string longer than the string length of the column is cut.
With `isnull` the entry is null and `cvals` is not used.

Preconditions (each throws `SpiceError`): write handle (`:dasnosuchhandle`,
`:dasinvalidaccess`); valid segment number (`:invalidindex`); the column exists
(`:bug`, as CSPICE) and has character type (`:wrongdatatype`); the row exists
(`:indexoutofrange`); a null entry only in a column that allows nulls
(`:badattribute`); for an array column at least one string, and as many as a fixed
entry size states (`:invalidcount`). Postcondition: the entry reads back as `cvals`
(cut to the string length), and the index of an indexed column holds the row.
"""
ekacec(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString,
    cvals::AbstractVector, isnull::Bool=false) =
    _ekw_entry_c(h, "ekacec", false, segno, recno, column, cvals, isnull)

"""
    ekaced(handle, segno, recno, column, dvals, isnull=false) -> Nothing

Add the entry of the double precision or time column `column` to row `recno` of
segment `segno`, as `ekaced_c`. Contract as [`ekacec`](@ref); in addition the entry of
a scalar column must not have a value yet (`:nonemptyentry`).
"""
ekaced(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString,
    dvals::AbstractVector, isnull::Bool=false) = _ekw_entry_n(h, "ekaced", false, Float64,
    (_EK_DP, _EK_TIME), segno, recno, column, dvals, isnull)

"""
    ekacei(handle, segno, recno, column, ivals, isnull=false) -> Nothing

Add the entry of the integer column `column` to row `recno` of segment `segno`, as
`ekacei_c`. Contract as [`ekaced`](@ref); the values must be 32-bit integers.
"""
ekacei(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString,
    ivals::AbstractVector, isnull::Bool=false) = _ekw_entry_n(h, "ekacei", false, Int32,
    (_EK_INT,), segno, recno, column, ivals, isnull)

"""
    ekucec(handle, segno, recno, column, cvals, isnull=false) -> Nothing

Replace the entry of the character column `column` in row `recno` of segment `segno`,
as `ekucec_c`. The new entry can have another size than the old one, and it can be
null. Contract as [`ekacec`](@ref). Postcondition: the entry reads back as `cvals`;
the data pages that only the old entry used are free.
"""
ekucec(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString,
    cvals::AbstractVector, isnull::Bool=false) =
    _ekw_entry_c(h, "ekucec", true, segno, recno, column, cvals, isnull)

"""
    ekuced(handle, segno, recno, column, dvals, isnull=false) -> Nothing

Replace the entry of the double precision or time column `column` in row `recno` of
segment `segno`, as `ekuced_c`. Contract as [`ekucec`](@ref).
"""
ekuced(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString,
    dvals::AbstractVector, isnull::Bool=false) = _ekw_entry_n(h, "ekuced", true, Float64,
    (_EK_DP, _EK_TIME), segno, recno, column, dvals, isnull)

"""
    ekucei(handle, segno, recno, column, ivals, isnull=false) -> Nothing

Replace the entry of the integer column `column` in row `recno` of segment `segno`, as
`ekucei_c`. Contract as [`ekucec`](@ref).
"""
ekucei(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString,
    ivals::AbstractVector, isnull::Bool=false) = _ekw_entry_n(h, "ekucei", true, Int32,
    (_EK_INT,), segno, recno, column, ivals, isnull)

# ---------------------------------------------------------------------------------------
# Fast writers
# ---------------------------------------------------------------------------------------

"""
    EKFastLoad <: AbstractVector{Int32}

The state of one fast load, which [`ekifld`](@ref) starts and [`ekffld`](@ref) ends.
As a vector it holds the record pointers of the rows of the new segment (the `rcptrs`
of CSPICE; empty for a segment with fixed-count columns). It also holds the handle,
the segment number and, for every column that was stored, the column position and the
data pointers of its rows (CSPICE keeps these in one global stack, so it can do one
fast load at a time; here every load has its own state).
"""
struct EKFastLoad <: AbstractVector{Int32}
    handle::DASHandle
    segno::Int
    ptrs::Vector{Int32}
    stack::Vector{Int32}
end
Base.size(r::EKFastLoad) = size(r.ptrs)
Base.getindex(r::EKFastLoad, i::Int) = r.ptrs[i]
Base.IndexStyle(::Type{EKFastLoad}) = IndexLinear()

"""
    ekifld(handle, tabnam, nrows, cnames, decls) -> (segno, rcptrs::EKFastLoad)

Start a new segment of `nrows` rows that is written column by column, as `ekifld_c`.
Store every column once with [`ekaclc`](@ref), [`ekacld`](@ref) or [`ekacli`](@ref)
and then call [`ekffld`](@ref); between `ekifld` and `ekffld` no other call must
change the file. Columns that are declared with `FIXED_COUNT = TRUE` give a segment
with a fixed layout (all columns of the segment must then have that attribute).

Preconditions: `nrows >= 1` (`:invalidcount`) and those of [`ekbseg`](@ref).
Postcondition: the segment exists with its metadata; it must not be read before
`ekffld` has returned.
"""
function ekifld(h::DASHandle, tabnam::AbstractString, nrows::Integer,
        cnames::AbstractVector{<:AbstractString}, decls::AbstractVector{<:AbstractString})
    nrows < 1 && throw(SpiceError(:invalidcount, "Number of rows must be > 0, was $nrows. "))
    segno = ekbseg(h, tabnam, cnames, decls)
    _, mbase = _ek_mloc(h, segno)
    h.ints[mbase+_EK_NRIDX] = nrows
    segdsc = h.ints[mbase+1:mbase+_EK_SDSCSZ]
    ncols = Int(segdsc[_EK_NCIDX])
    ptrs = Int32[]
    if segdsc[_EK_EKTIDX] == 1           # `zzekif01`: the pages of the records
        rpsize = ncols + _EK_DPTBAS
        nrp = _EK_IPSIZE ÷ rpsize
        npage = (Int(nrows) + nrp - 1) ÷ nrp
        remain = Int(nrows)
        sizehint!(ptrs, remain)
        for _ in 1:npage                 # invariant: `ptrs` holds the pointers of the pages before
            _, pbase = _ekw_aps(h, segdsc, 3)
            nr = min(nrp, remain)
            for j in 1:nr
                push!(ptrs, pbase + (j - 1) * rpsize)
            end
            remain -= nr
        end
    else                                 # `zzekif02`: the metadata of the columns
        _, base = _ekw_aps(h, segdsc, 3)
        offset = base
        for i in 1:ncols
            dscbas = mbase + _EK_SDSCSZ + (i - 1) * _EK_CDSCSZ
            h.ints[dscbas+_EK_METIDX] = offset
            offset += 2
        end
    end
    return segno, EKFastLoad(h, Int(segno), ptrs, Int32[])
end

"""
    _ekw_fast_locate(h, who, segno, column, rcptrs, types) -> (segdsc, coldsc)

The common start of `ekaclc`, `ekacld` and `ekacli`. Throws `SpiceError`:
`:dasnosuchhandle`, `:dasinvalidaccess`, `:emptystring`, `:invalidindex`,
`:wrongsegment` if `rcptrs` is the state of another fast load, `:bug` for an unknown
column, `:wrongdatatype`.
"""
function _ekw_fast_locate(h::DASHandle, who::AbstractString, segno::Integer,
        column::AbstractString, rcptrs::EKFastLoad, types::Tuple)
    _ekw_check(h, who)
    isempty(column) && throw(SpiceError(:emptystring, "$who: the column name is empty"))
    segdsc = Vector{Int32}(_ek_sdsc(h, segno))
    (rcptrs.handle === h && rcptrs.segno == segno) || throw(SpiceError(:wrongsegment,
        "$who: the record pointers belong to the fast load of segment $(rcptrs.segno) of " *
        "$(rcptrs.handle.path), not to segment $segno of $(h.path)."))
    coldsc = _ekw_cdsc(h, segdsc, column)
    coldsc[_EK_TYPIDX] in types || throw(SpiceError(:wrongdatatype, "Column $column is of " *
        "type $(coldsc[_EK_TYPIDX]); $(uppercase(who)) only works with " *
        "$(_ekw_typename(types)) columns."))
    return segdsc, coldsc
end

"""
    _ekw_fast_counts(who, segdsc, coldsc, nvals, entszs, nlflgs) -> Nothing

Check that the inputs of a fast writer hold what the column needs: a null flag for
every row of a column that allows nulls, an entry size for every row of a column with
entries of variable size, and the values of all rows (a null row takes the room of
its entry size). Throws `SpiceError(:invalidcount)`. CSPICE cannot make this check.
"""
function _ekw_fast_counts(who::AbstractString, segdsc::Vector{Int32}, coldsc::Vector{Int32},
        nvals::Int, entszs::AbstractVector, nlflgs::AbstractVector)
    nrows = Int(segdsc[_EK_NRIDX])
    size = Int(coldsc[_EK_SIZIDX])
    (coldsc[_EK_NFLIDX] != -1 && length(nlflgs) < nrows) && throw(SpiceError(:invalidcount,
        "$who: $(length(nlflgs)) null flags for $nrows rows"))
    need = nrows
    if coldsc[_EK_CLSIDX] in (4, 5, 6)
        (size == -1 && length(entszs) < nrows) && throw(SpiceError(:invalidcount,
            "$who: $(length(entszs)) entry sizes for $nrows rows"))
        need = size == -1 ? sum(Int(entszs[i]) for i in 1:nrows) : nrows * size
    end
    nvals < need && throw(SpiceError(:invalidcount,
        "$who: $nvals values, but the $nrows rows of the column need $need"))
    return nothing
end

@noinline _ekw_notsupported(coldsc, column) = throw(SpiceError(:notsupported, "Unsupported " *
    "column class code $(coldsc[_EK_CLSIDX]) found in descriptor for column $column."))

"The null flags as a vector of `Bool`."
_ekw_flags(nlflgs::AbstractVector) = nlflgs isa AbstractVector{Bool} ? nlflgs :
                                     Bool[f != 0 for f in nlflgs]

"""
    ekaclc(handle, segno, column, cvals, entszs, nlflgs, rcptrs) -> Nothing

Store the whole character column `column` of the segment that [`ekifld`](@ref)
started, as `ekaclc_c`. `cvals` holds the strings of all rows in row order (for an
array column the strings of row 1, then those of row 2, and so on); `entszs[r]` is the
number of strings of row `r` (used only for entries of variable size); `nlflgs[r]`
tells that row `r` is null (used only if the column allows nulls). A null row also
takes the room of its entry in `cvals`. `rcptrs` is the result of `ekifld`.

Preconditions (each throws `SpiceError`): write handle; valid segment number
(`:invalidindex`); `rcptrs` belongs to this segment (`:wrongsegment`); the column
exists (`:bug`) and has character type (`:wrongdatatype`); the inputs hold all rows
(`:invalidcount`). A fixed-count column needs a string length of 1 to 1014
(`:invalidsize`).
"""
function ekaclc(h::DASHandle, segno::Integer, column::AbstractString, cvals::AbstractVector,
        entszs::AbstractVector, nlflgs::AbstractVector, rcptrs::EKFastLoad)
    segdsc, coldsc = _ekw_fast_locate(h, "ekaclc", segno, column, rcptrs, (_EK_CHR,))
    class = coldsc[_EK_CLSIDX]
    (class == 3 || class == 6 || class == 9) || _ekw_notsupported(coldsc, column)
    _ekw_fast_counts("ekaclc", segdsc, coldsc, length(cvals), entszs, nlflgs)
    vals = _ekw_cvals(cvals)
    flags = _ekw_flags(nlflgs)
    if class == 3
        _ekw_ac03(h, segdsc, coldsc, vals, flags, rcptrs, rcptrs.stack)
    elseif class == 6
        _ekw_ac06(h, segdsc, coldsc, vals, entszs, flags, rcptrs.stack)
    else
        _ekw_ac09(h, segdsc, coldsc, vals, flags)
    end
    return nothing
end

"The body of `ekacld` and `ekacli`; `T` is the element type and `types` the data type codes."
function _ekw_fast_n(h::DASHandle, who::AbstractString, ::Type{T}, types::Tuple, segno::Integer,
        column::AbstractString, nvals::AbstractVector, entszs::AbstractVector,
        nlflgs::AbstractVector, rcptrs::EKFastLoad) where {T}
    segdsc, coldsc = _ekw_fast_locate(h, who, segno, column, rcptrs, types)
    class = coldsc[_EK_CLSIDX]
    kind = class == 1 || class == 2 ? 1 : class == 4 || class == 5 ? 2 :
           class == 7 || class == 8 ? 3 : 0
    kind == 0 && _ekw_notsupported(coldsc, column)
    _ekw_fast_counts(who, segdsc, coldsc, length(nvals), entszs, nlflgs)
    vals = _ekw_nvals(T, nvals)
    flags = _ekw_flags(nlflgs)
    if kind == 1
        _ekw_ac12(h, segdsc, coldsc, vals, flags, rcptrs, rcptrs.stack)
    elseif kind == 2
        _ekw_ac45(h, segdsc, coldsc, vals, entszs, flags, rcptrs.stack)
    else
        _ekw_ac78(h, segdsc, coldsc, vals, flags)
    end
    return nothing
end

"""
    ekacld(handle, segno, column, dvals, entszs, nlflgs, rcptrs) -> Nothing

Store the whole double precision or time column `column` of the segment that
[`ekifld`](@ref) started, as `ekacld_c`. Contract as [`ekaclc`](@ref).
"""
ekacld(h::DASHandle, segno::Integer, column::AbstractString, dvals::AbstractVector,
    entszs::AbstractVector, nlflgs::AbstractVector, rcptrs::EKFastLoad) =
    _ekw_fast_n(h, "ekacld", Float64, (_EK_DP, _EK_TIME), segno, column, dvals, entszs, nlflgs,
        rcptrs)

"""
    ekacli(handle, segno, column, ivals, entszs, nlflgs, rcptrs) -> Nothing

Store the whole integer column `column` of the segment that [`ekifld`](@ref) started,
as `ekacli_c`. Contract as [`ekaclc`](@ref).
"""
ekacli(h::DASHandle, segno::Integer, column::AbstractString, ivals::AbstractVector,
    entszs::AbstractVector, nlflgs::AbstractVector, rcptrs::EKFastLoad) =
    _ekw_fast_n(h, "ekacli", Int32, (_EK_INT,), segno, column, ivals, entszs, nlflgs, rcptrs)

"""
    _ekw_ff01(h, segno, rcptrs) -> Nothing

End the fast load of a type 1 segment (`zzekff01`): the records get their status and
the data pointers of all columns (from the stack of the load), the record tree is
made, and the segment descriptor gets the tree and marks the last data pages as full.
Precondition: every column was stored once, else `SpiceError(:invalidaddress)`.
"""
function _ekw_ff01(h::DASHandle, segno::Integer, rcptrs::EKFastLoad)
    _, mbase = _ek_mloc(h, segno)
    ints = h.ints
    nrows = Int(ints[mbase+_EK_NRIDX])
    ncols = Int(ints[mbase+_EK_NCIDX])
    stack = rcptrs.stack
    length(stack) >= ncols * (nrows + 1) || throw(SpiceError(:invalidaddress,
        "ekffld: $(length(stack) ÷ (nrows + 1)) of the $ncols columns of segment $segno " *
        "were stored"))
    colord = zeros(Int, ncols)
    for i in 1:ncols                     # invariant: `colord` knows the columns of the blocks before i
        colidx = Int(stack[(i-1)*(nrows+1)+1])
        (1 <= colidx <= ncols) || _ekw_bug("ekffld: column position $colidx on the stack")
        colord[colidx] = i
    end
    all(>(0), colord) || throw(SpiceError(:invalidaddress,
        "ekffld: a column of segment $segno was stored twice and another one not at all"))
    rpsize = ncols + _EK_DPTBAS
    for (r, recptr) in enumerate(rcptrs.ptrs)
        # invariant: the records of the rows before r are complete
        ints[recptr+1] = 1
        ints[recptr+2] = -1
        for col in 1:ncols
            ints[recptr+_EK_DPTBAS+col] = stack[(colord[col]-1)*(nrows+1)+1+r]
        end
    end
    # CSPICE writes whole pages and so clears the words after the last record of a page.
    nrp = _EK_IPSIZE ÷ rpsize
    for first in 1:nrp:nrows
        nr = min(nrp, nrows - first + 1)
        pbase = Int(rcptrs.ptrs[first])
        fill!(view(ints, pbase+nr*rpsize+1:pbase+_EK_IPSIZE), Int32(0))
    end
    tree = _ekw_trit(h)
    _ekw_tr1s(h, tree, rcptrs.ptrs)
    ints = h.ints
    ints[mbase+_EK_RTIDX] = tree
    ints[mbase+_EK_NRIDX] = nrows
    ints[mbase+_EK_LCWIDX] = _EK_CPSIZE
    ints[mbase+_EK_LCWIDX+1] = _EK_DPSIZE
    ints[mbase+_EK_LCWIDX+2] = _EK_IPSIZE
    return nothing
end

"""
    ekffld(handle, segno, rcptrs) -> Nothing

End the fast load of segment `segno` that [`ekifld`](@ref) started, as `ekffld_c`.
Preconditions (each throws `SpiceError`): write handle; valid segment number
(`:invalidindex`); `rcptrs` is the result of the `ekifld` call for this segment
(`:wrongsegment`); every column was stored once (`:invalidaddress`). Postcondition:
the segment is complete and can be read; rows can be added to it later, unless it has
fixed-count columns.
"""
function ekffld(h::DASHandle, segno::Integer, rcptrs::EKFastLoad)
    _ekw_check(h, "ekffld")
    segdsc = _ek_sdsc(h, segno)
    stype = segdsc[_EK_EKTIDX]
    (rcptrs.handle === h && rcptrs.segno == segno) || throw(SpiceError(:wrongsegment,
        "Attempt to finish fast load of wrong segment.  Input segment number is $segno; " *
        "stacked segment number is $(rcptrs.segno).  Input file is $(h.path); stacked " *
        "file is $(rcptrs.handle.path)."))
    if stype == 1
        _ekw_ff01(h, segno, rcptrs)
    elseif stype != 2
        _ekw_bug("Segment type $stype is not currently supported.")
    end
    return nothing
end
