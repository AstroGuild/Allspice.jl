# Writing an EK file, part 2: data pages of a segment, column entries and indexes.
#
# The nine column classes are listed in `columns.jl`. This file ports, for the classes
# 1 to 6 (segment type 1), the routines that add, update and delete one column entry
# (`zzekad0n`, `zzekue0n`, `zzekde0n`) and that keep the index tree of an indexed
# column in order (`zzekiii1`, `zzekiid1`, `zzekiic1`, `zzekixdl`), and for all classes
# the routines of the fast writers that store a whole column (`zzekac0n`).
#
# Conventions of this file:
#
#   - `segdsc` is a copy of the 24 integers of the segment descriptor and `coldsc` a copy
#     of the 11 integers of a column descriptor, as in CSPICE. A routine changes
#     `segdsc` and writes it back to the file where CSPICE does.
#   - The type code `t` is 1 for character, 2 for double and 3 for integer pages.
#   - A character value is a byte vector without trailing blanks; it stands for the
#     Fortran string that continues with blanks.
#
# The order of the page allocations and of the words written is that of CSPICE.

# Positions of the last data page and of its last word in use, by page type.
@inline _ekw_lcp(t::Int) = _EK_LCPIDX + t - 1
@inline _ekw_lcw(t::Int) = _EK_LCWIDX + t - 1

@inline _ekw_ntype(coldsc::Vector{Int32}) = coldsc[_EK_TYPIDX] == _EK_INT ? 3 : 2

@inline function _ekw_put_segdsc(h::DASHandle, segdsc::Vector{Int32})
    copyto!(h.ints, Int(segdsc[_EK_IMDIDX]) + 1, segdsc, 1, _EK_SDSCSZ)
    return nothing
end

@inline function _ekw_putc!(dst::Vector{UInt8}, at::Int, val::AbstractVector{UInt8}, pos::Int,
        n::Int)
    m = clamp(length(val) - pos + 1, 0, n)
    m > 0 && copyto!(dst, at, val, pos, m)
    m < n && fill!(view(dst, at+m:at+n-1), UInt8(' '))
    return nothing
end

@inline function _ekw_enc!(buf::Vector{UInt8}, at::Int, val::Integer)
    n = Int(val)
    for k in 4:-1:0                      # invariant: digits after k are stored
        buf[at+k] = (n % 128) % UInt8
        n ÷= 128
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Data pages of a segment
# ---------------------------------------------------------------------------------------

"""
    _ekw_aps(h, segdsc, t) -> (p, base)

Allocate a data page of type `t` for the segment (`zzekaps` with NEW = false): the page
gets the link count 0 and no forward pointer, and its number is appended to the data
page tree of its type.
"""
function _ekw_aps(h::DASHandle, segdsc::Vector{Int32}, t::Int)
    p, base = _ekw_pgal(h, t)
    _ekw_slnk(h, t, p, 0)
    _ekw_sfwd(h, t, p, 0)
    _ekw_trap(h, segdsc[_EK_CPTIDX+t-1], p)
    return p, base
end

"""
    _ekw_acps(h, segdsc, t, n) -> (p, base)

Allocate `n` new adjacent data pages of type `t` for the segment (`zzekacps`); `p` is
the first one. CSPICE appends one page also for `n < 1`; so does this function.
"""
function _ekw_acps(h::DASHandle, segdsc::Vector{Int32}, t::Int, n::Int)
    p, base = _ekw_pgan(h, t)
    for _ in 2:n
        _ekw_pgan(h, t)
    end
    for i in 1:n
        _ekw_slnk(h, t, p + i - 1, 0)
        _ekw_sfwd(h, t, p + i - 1, 0)
    end
    tree = segdsc[_EK_CPTIDX+t-1]
    for i in 1:n                         # invariant: pages p .. p+i-2 are in the tree
        _ekw_trap(h, tree, p + i - 1)
    end
    return p, base
end

"""
    _ekw_dps(h, segdsc, t, p) -> Nothing

Free the data page `p` of type `t` of the segment (`zzekdps`): the page goes to the
free list and leaves the data page tree. If it is the last data page of its type, the
segment descriptor (copy and file) marks that page as full, so the next entry starts
a new page.
"""
function _ekw_dps(h::DASHandle, segdsc::Vector{Int32}, t::Int, p::Integer)
    _ekw_pgfr(h, t, p)
    if segdsc[_ekw_lcp(t)] == p
        h.ints[Int(segdsc[_EK_IMDIDX])+_ekw_lcw(t)] = _ekw_DSIZE[t]
        segdsc[_ekw_lcw(t)] = _ekw_DSIZE[t]
    end
    tree = segdsc[_EK_CPTIDX+t-1]
    _ekw_trdl(h, tree, _ekw_trls(h, tree, p))
    return nothing
end

@inline function _ekw_unlink(h::DASHandle, segdsc::Vector{Int32}, t::Int, p::Integer)
    nlinks = _ekw_glnk(h, t, p)
    if nlinks > 1
        _ekw_slnk(h, t, p, nlinks - 1)
    else
        _ekw_dps(h, segdsc, t, p)
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Comparison of a column entry with a key (`zzekscmp`)
# ---------------------------------------------------------------------------------------

"The decoded descriptor of `coldsc`, without the name, for the readers of `columns.jl`."
@inline _ekw_col(coldsc::AbstractVector{Int32}) = _EKColumn("", coldsc[_EK_CLSIDX],
    coldsc[_EK_TYPIDX], coldsc[_EK_LENIDX], coldsc[_EK_SIZIDX], coldsc[_EK_IXTIDX],
    coldsc[_EK_IXPIDX], coldsc[_EK_NFLIDX] == 1, coldsc[_EK_ORDIDX], coldsc[_EK_METIDX])

_ekw_cnam(s, coldsc::AbstractVector{Int32}) =
    _ek_rstr(s.chars, Int(coldsc[_EK_NAMIDX]) + 1, Int(coldsc[_EK_NAMIDX]) + _EK_CNAMSZ)

"""
    _ekw_cmpc(a, b) -> Int

Order of two strings as Fortran compares them: byte by byte, the shorter one continued
with blanks. Returns -1, 0 or 1.
"""
function _ekw_cmpc(a::AbstractVector{UInt8}, b::AbstractVector{UInt8})
    la = length(a)
    lb = length(b)
    for k in 1:max(la, lb)               # invariant: the characters before k are equal
        ca = k <= la ? a[k] : UInt8(' ')
        cb = k <= lb ? b[k] : UInt8(' ')
        ca == cb || return ca < cb ? -1 : 1
    end
    return 0
end

"The first `min(length, 1024)` bytes of a character key, as `zzekscmp` uses them."
@inline _ekw_key1024(key::AbstractVector{UInt8}) = view(key, 1:min(length(key), 1024))

"""
    _ekw_rel(s, coldsc, row, key, null) -> Int

Relation of the first element of the entry of the column `coldsc` in the row `row` (a
record pointer; a row number for the classes 7 to 9) to the key, as `zzekscmp` finds
it: -1 less, 0 equal, 1 greater. A null entry is less than every value and equal to a
null key (`null` true). `key` is a byte vector for a character column and a number for
the other types; an integer column is compared with a double key as doubles. Of a
string, the first 1024 characters count. Throws `SpiceError`: `:invaliddatatype` for
an unknown column type, `:bug` for a key of the wrong type, and the errors of the
readers (`:uninitializedvalue`, `:uninitialized`).
"""
function _ekw_rel(s, coldsc::AbstractVector{Int32}, row::Integer, key, null::Bool)
    col = _ekw_col(coldsc)
    dtype = col.dtype
    if dtype == _EK_CHR
        cv, enull, _ = _ek_elt_c(s, col, row, 1)
        enull && return null ? 0 : -1
        null && return 1
        key isa AbstractVector{UInt8} || _ekw_bug("Column type is $dtype; value type is not character.")
        elt = codeunits(cv)
        return _ekw_cmpc(view(elt, 1:min(length(elt), 1024)), _ekw_key1024(key))
    elseif dtype == _EK_DP || dtype == _EK_TIME
        dv, enull, _ = _ek_elt_d(s, col, row, 1)
        enull && return null ? 0 : -1
        null && return 1
        key isa Real || _ekw_bug("Column type is $dtype; value type is not numeric.")
        numval = Float64(key)
        return dv < numval ? -1 : dv > numval ? 1 : 0
    elseif dtype == _EK_INT
        iv, enull, _ = _ek_elt_i(s, col, row, 1)
        enull && return null ? 0 : -1
        null && return 1
        key isa Real || _ekw_bug("Column type is $dtype; value type is not numeric.")
        numval = Float64(key)
        return iv < numval ? -1 : iv > numval ? 1 : 0
    end
    throw(SpiceError(:invaliddatatype, "Data type code $dtype not recognized."))
end

# ---------------------------------------------------------------------------------------
# Indexes
# ---------------------------------------------------------------------------------------

"""
    _ekw_ixlk(s, coldsc, key) -> recptr

The record pointer at position `key` of the index of the column `coldsc` (`zzekixlk`):
from the index tree (index type 1) or from the order vector (type 2, row numbers).
Throws `SpiceError(:invalidtype)` for another index type.
"""
function _ekw_ixlk(s, coldsc::AbstractVector{Int32}, key::Integer)
    itype = coldsc[_EK_IXTIDX]
    itype == 1 && return _ek_trdp(s.ints, coldsc[_EK_IXPIDX], key)
    itype == 2 || throw(SpiceError(:invalidtype, "The index type $itype is not supported."))
    q = (Int(key) - 1) ÷ _EK_IPSIZE
    r = Int(key) - q * _EK_IPSIZE
    return Int(s.ints[Int(coldsc[_EK_IXPIDX])+q*256+r])
end

function _ekw_index1_check(s, coldsc::AbstractVector{Int32})
    itype = coldsc[_EK_IXTIDX]
    itype == -1 && throw(SpiceError(:notindexed, "Column $(_ekw_cnam(s, coldsc)) is not indexed."))
    itype == 1 || throw(SpiceError(:invalidtype,
        "Column $(_ekw_cnam(s, coldsc)) has index type $itype."))
    return nothing
end

"""
    _ekw_ler(s, segdsc, coldsc, key, recptr, null) -> (prvidx, prvptr)

Position in the index tree of the column `coldsc` of the last entry that is not after
the pair (`key`, `recptr`), and the record pointer at that position; `(0, 0)` if no
entry qualifies (`zzeklerc`, `zzeklerd`, `zzekleri` with `zzekerc1`, `zzekerd1`,
`zzekeri1`). The index is in the order of the entry values, null entries first, and
entries with equal values in the order of their record pointers. Precondition: index
of type 1 (`:notindexed`, `:invalidtype`) that holds no more entries than the segment
has rows (`:sizemismatch`).
"""
function _ekw_ler(s, segdsc::AbstractVector{Int32}, coldsc::AbstractVector{Int32}, key,
        recptr::Integer, null::Bool)
    _ekw_index1_check(s, coldsc)
    ints = s.ints
    tree = coldsc[_EK_IXPIDX]
    tsize = _ek_trsz(ints, tree)
    nrec = Int(segdsc[_EK_NRIDX])
    tsize > nrec && throw(SpiceError(:sizemismatch,
        "Index size = $tsize but column contains $nrec records."))
    tsize == 0 && return 0, 0
    lo = 1
    hi = tsize
    begptr = _ek_trdp(ints, tree, lo)
    endptr = _ek_trdp(ints, tree, hi)
    r = _ekw_rel(s, coldsc, begptr, key, null)
    (r > 0 || (r == 0 && recptr < begptr)) && return 0, 0
    r = _ekw_rel(s, coldsc, endptr, key, null)
    (r < 0 || (r == 0 && endptr <= recptr)) && return tsize, endptr
    while hi > lo + 1
        # invariant: the entry at `lo` is not after the pair, the entry at `hi` is
        mid = (lo + hi) ÷ 2
        midptr = _ek_trdp(ints, tree, mid)
        r = _ekw_rel(s, coldsc, midptr, key, null)
        if r < 0 || (r == 0 && midptr <= recptr)
            lo = mid
        else
            hi = mid
        end
    end
    return lo, _ek_trdp(ints, tree, lo)
end

"""
    _ekw_ixin(h, segdsc, coldsc, key, recptr, null) -> Nothing

Insert the record `recptr`, whose entry in the indexed column `coldsc` is `key` (or
null), into the index of the column (`zzekiii1`, `zzekiid1`, `zzekiic1`).
Postcondition: the index is in order and has one entry more.
"""
function _ekw_ixin(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, key,
        recptr::Integer, null::Bool)
    prvidx, _ = _ekw_ler(h, segdsc, coldsc, key, recptr, null)
    _ekw_trin(h, coldsc[_EK_IXPIDX], prvidx + 1, recptr)
    return nothing
end

"The value of the scalar entry of `coldsc` in the record `recptr` as an index key: `(key, isnull)`."
function _ekw_own_key(s, coldsc::AbstractVector{Int32}, recptr::Integer)
    col = _ekw_col(coldsc)
    dtype = col.dtype
    if dtype == _EK_CHR
        cv, isnull, _ = _ek_elt_c(s, col, recptr, 1)
        bytes = Vector{UInt8}(codeunits(cv))
        length(bytes) > 1024 && resize!(bytes, 1024)
        return bytes, isnull
    elseif dtype == _EK_DP || dtype == _EK_TIME
        dv, isnull, _ = _ek_elt_d(s, col, recptr, 1)
        return dv, isnull
    elseif dtype == _EK_INT
        iv, isnull, _ = _ek_elt_i(s, col, recptr, 1)
        return iv, isnull
    end
    throw(SpiceError(:itemnotfound, "File = $(s.path); COLIDX = $(col.ord). Unrecognized " *
        "data type code $dtype found in descriptor."))
end

"""
    _ekw_frx(h, segdsc, coldsc, recptr) -> pos

Position of the record `recptr` in the index of the column `coldsc` (`zzekfrx`). The
entry of the record must still hold the value that the index was built with.

CSPICE looks only at the position that the order "by value, then by record pointer"
gives, and it signals ITEMNOTFOUND if another record is there. Two CSPICE writers do
not keep that order among equal values: the fast writer of a character column
(`zzekordc`, a sort that is not stable) and `zzekad03`, which puts a string that is
longer than the string length of the column at the place of the whole string. So this
function also looks at the other entries with an equal value, and it finds every
record that the index holds. Throws `SpiceError(:itemnotfound)` if the index does not
hold the record.
"""
function _ekw_frx(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Integer)
    key, isnull = _ekw_own_key(h, coldsc, recptr)
    pos, prvptr = _ekw_ler(h, segdsc, coldsc, key, recptr, isnull)
    prvptr == recptr && return pos
    ints = h.ints
    tree = coldsc[_EK_IXPIDX]
    tsize = _ek_trsz(ints, tree)
    # `pos` is in the run of the entries with an equal value, or just before that run.
    k = pos
    while k >= 1                         # invariant: the entries k+1 .. pos are not the record
        ptr = _ek_trdp(ints, tree, k)
        ptr == recptr && return k
        _ekw_rel(h, coldsc, ptr, key, isnull) == 0 || break
        k -= 1
    end
    k = pos + 1
    while k <= tsize                     # invariant: the entries pos+1 .. k-1 are not the record
        ptr = _ek_trdp(ints, tree, k)
        ptr == recptr && return k
        _ekw_rel(h, coldsc, ptr, key, isnull) == 0 || break
        k += 1
    end
    throw(SpiceError(:itemnotfound, "File = $(h.path); COLIDX = $(coldsc[_EK_ORDIDX]).  The " *
        "index does not hold the record.  This probably indicates a corrupted file or a " *
        "bug in the EK code."))
end

"""
    _ekw_ixdl(h, segdsc, coldsc, recptr) -> Nothing

Delete the record `recptr` from the index of the column `coldsc` (`zzekixdl`).
Precondition: the column has an index (`:bug`) of type 1 (`:invalidtype`), and the
entry of the record is not changed yet.
"""
function _ekw_ixdl(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Integer)
    idxtyp = coldsc[_EK_IXTIDX]
    idxtyp == -1 && _ekw_bug("Column was not indexed. File = $(h.path); COLIDX = " *
        "$(coldsc[_EK_ORDIDX]).")
    idxtyp == 1 || throw(SpiceError(:invalidtype, "The index type $idxtyp is not supported."))
    pos = _ekw_frx(h, segdsc, coldsc, recptr)
    _ekw_trdl(h, coldsc[_EK_IXPIDX], pos)
    return nothing
end

"""
    _ekw_llx(s, segdsc, coldsc, key, orequal) -> (prvloc, prvptr)

Position in the index of the column `coldsc` of the last entry that is less than `key`
(`orequal` false: `zzeklltc`, `zzeklltd`, `zzekllti`) or not greater than `key`
(`orequal` true: `zzekllec`, `zzeklled`, `zzekllei`), and the record pointer at that
position; `(0, 0)` if no entry qualifies. Works for both index types. Precondition:
the column has an index (`:notindexed`).
"""
function _ekw_llx(s, segdsc::AbstractVector{Int32}, coldsc::AbstractVector{Int32}, key,
        orequal::Bool)
    coldsc[_EK_IXTIDX] == -1 && throw(SpiceError(:notindexed,
        "Column $(_ekw_cnam(s, coldsc)) is not indexed."))
    nrows = Int(segdsc[_EK_NRIDX])
    nrows == 0 && return 0, 0
    # `before(ptr)` tells that the entry of the record is on the wanted side of the key.
    before(ptr) = (r = _ekw_rel(s, coldsc, ptr, key, false); orequal ? r <= 0 : r < 0)
    lo = 1
    hi = nrows
    before(_ekw_ixlk(s, coldsc, lo)) || return 0, 0
    endptr = _ekw_ixlk(s, coldsc, hi)
    before(endptr) && return nrows, endptr
    while hi > lo + 1                    # invariant: `before` holds at `lo` and not at `hi`
        mid = (lo + hi) ÷ 2
        if before(_ekw_ixlk(s, coldsc, mid))
            lo = mid
        else
            hi = mid
        end
    end
    return lo, _ekw_ixlk(s, coldsc, lo)
end

"The key `val` in the type that the column `coldsc` compares with (`zzekillt`, `zzekille`)."
function _ekw_index_key(coldsc::AbstractVector{Int32}, val)
    coltyp = coldsc[_EK_TYPIDX]
    if coltyp == _EK_CHR
        return val isa AbstractString ? Vector{UInt8}(codeunits(rstrip(val, ' '))) : val
    elseif coltyp == _EK_DP || coltyp == _EK_TIME
        return Float64(val)
    elseif coltyp == _EK_INT
        return val isa Integer ? Int(val) : round(Int, val, RoundNearestTiesAway)
    end
    throw(SpiceError(:invalidsize, "The data type $coltyp is not supported."))
end

"""
    _ekw_illt(s, segdsc, coldsc, nrows, val) -> Int

Number of entries of the index of the column `coldsc` that are less than `val`
(`zzekillt`); `val` is a string for a character column and a number else (a double is
rounded for an integer column). Precondition: `nrows >= 1`, else `:invalidsize`.
"""
function _ekw_illt(s, segdsc::AbstractVector{Int32}, coldsc::AbstractVector{Int32},
        nrows::Integer, val)
    nrows < 1 && throw(SpiceError(:invalidsize, "Number of rows must be positive; was $nrows."))
    return _ekw_llx(s, segdsc, coldsc, _ekw_index_key(coldsc, val), false)[1]
end

"As [`_ekw_illt`](@ref) for the entries that are not greater than `val` (`zzekille`)."
function _ekw_ille(s, segdsc::AbstractVector{Int32}, coldsc::AbstractVector{Int32},
        nrows::Integer, val)
    nrows < 1 && throw(SpiceError(:invalidsize, "Number of rows must be positive; was $nrows."))
    return _ekw_llx(s, segdsc, coldsc, _ekw_index_key(coldsc, val), true)[1]
end

# ---------------------------------------------------------------------------------------
# Add one entry (classes 1 to 6)
# ---------------------------------------------------------------------------------------

"""
    _ekw_entry_checks(segdsc, coldsc, isnull) -> colidx

The first checks of every `zzekad0n`: the ordinal position of the column is inside the
segment (`:invalidindex`) and a null entry is given only to a column that allows nulls
(`:badattribute`). Returns the ordinal position.
"""
@inline function _ekw_entry_checks(segdsc::Vector{Int32}, coldsc::Vector{Int32}, isnull::Bool)
    ncols = Int(segdsc[_EK_NCIDX])
    colidx = Int(coldsc[_EK_ORDIDX])
    (1 <= colidx <= ncols) || throw(SpiceError(:invalidindex,
        "Column index = $colidx; valid range is 1:$ncols."))
    (isnull && coldsc[_EK_NFLIDX] != 1) && throw(SpiceError(:badattribute, "Column having " *
        "index $colidx in segment $(segdsc[_EK_SNOIDX]) does not allow nulls, but a null " *
        "value was supplied for the element."))
    return colidx
end

"""
    _ekw_empty_check(h, segdsc, colidx, ptrloc) -> Nothing

Throws `SpiceError(:nonemptyentry)` unless the entry whose data pointer is at `ptrloc`
was never written. `zzekad01` and `zzekad02` make this check; `zzekad03` to `zzekad06`
do not, and a second write then leaves the pages of the first entry linked and puts the
record into the index of the column twice.
"""
@inline function _ekw_empty_check(h::DASHandle, segdsc::Vector{Int32}, colidx::Int, ptrloc::Int)
    datptr = Int(h.ints[ptrloc])
    (datptr == _EK_UNINIT || datptr == _EK_NOBACK) || throw(SpiceError(:nonemptyentry,
        "Column having index $colidx in segment $(segdsc[_EK_SNOIDX]) has non-empty " *
        "element in the record."))
    return nothing
end

"The count checks of `zzekad04`, `zzekad05` and `zzekad06` (`:invalidcount`)."
@inline function _ekw_count_checks(segdsc::Vector{Int32}, coldsc::Vector{Int32}, nvals::Int)
    nvals < 1 && throw(SpiceError(:invalidcount, "COLIDX = $(coldsc[_EK_ORDIDX]);  segment " *
        "= $(segdsc[_EK_SNOIDX]); NVALS = $nvals;  NVALS must be positive "))
    size = coldsc[_EK_SIZIDX]
    (size != -1 && nvals != size) && throw(SpiceError(:invalidcount, "COLIDX = " *
        "$(coldsc[_EK_ORDIDX]);  segment = $(segdsc[_EK_SNOIDX]); NVALS = $nvals; declared " *
        "entry size = $size.  Sizes must match."))
    return nothing
end

@inline function _ekw_putn(h::DASHandle, t::Int, addr::Int, val::Real)
    if t == 2
        h.dps[addr] = val
    else
        h.ints[addr] = val
    end
    return nothing
end

@inline function _ekw_indexed(segdsc::Vector{Int32}, coldsc::Vector{Int32})
    itype = coldsc[_EK_IXTIDX]
    itype == -1 && return false
    itype == 1 || throw(SpiceError(:invalidtype, "Column having index $(coldsc[_EK_ORDIDX]) " *
        "in segment $(segdsc[_EK_SNOIDX]) has index type $itype."))
    return true
end

"""
    _ekw_ad12(h, segdsc, coldsc, recptr, val, isnull) -> Nothing

Add the entry `val` (or null) of a scalar integer, double or time column (class 1 or
2) to the record at `recptr` (`zzekad01`, `zzekad02`). The value takes the next free
word of the last data page of its type, or the first word of a new page; an indexed
column gets the record in its index. Preconditions (each throws `SpiceError`): valid
column position (`:invalidindex`), null allowed (`:badattribute`), entry not written
yet (`:nonemptyentry`). Postcondition: the segment descriptor in the file is `segdsc`.
"""
function _ekw_ad12(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Int,
        val::Real, isnull::Bool)
    t = _ekw_ntype(coldsc)
    ints = h.ints
    colidx = _ekw_entry_checks(segdsc, coldsc, isnull)
    ptrloc = recptr + _EK_DPTBAS + colidx
    _ekw_empty_check(h, segdsc, colidx, ptrloc)
    datptr = Int(ints[ptrloc])
    if isnull
        ints[ptrloc] = _EK_NULL
    else
        lastw = Int(segdsc[_ekw_lcw(t)])
        if lastw < _ekw_DSIZE[t]
            p = Int(segdsc[_ekw_lcp(t)])
            datptr = _ekw_base(t, p) + lastw + 1
            ints[ptrloc] = datptr
            _ekw_putn(h, t, datptr, val)
            _ekw_slnk(h, t, p, _ekw_glnk(h, t, p) + 1)
            segdsc[_ekw_lcw(t)] = lastw + 1
        else
            p, pbase = _ekw_aps(h, segdsc, t)
            _ekw_putn(h, t, pbase + 1, val)
            _ekw_slnk(h, t, p, 1)
            segdsc[_ekw_lcp(t)] = p
            segdsc[_ekw_lcw(t)] = 1
            h.ints[ptrloc] = pbase + 1
        end
    end
    _ekw_put_segdsc(h, segdsc)
    _ekw_indexed(segdsc, coldsc) && _ekw_ixin(h, segdsc, coldsc, val, recptr, isnull)
    return nothing
end

"""
    _ekw_ad03(h, segdsc, coldsc, recptr, cval, isnull) -> Nothing

Add the string `cval` (or null) of a scalar character column (class 3) to the record at
`recptr` (`zzekad03`). The file takes a five-character count and the characters
without trailing blanks (at least one; at most the string length of the column); the
string starts on the last character data page if that has room for the count and one
character, and continues on new pages. The index compares the whole of `cval`, as in
CSPICE, also if the column stores less of it. Preconditions: as `_ekw_ad12`; CSPICE
does not check that the entry is empty (see `_ekw_empty_check`).
"""
function _ekw_ad03(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Int,
        cval::Vector{UInt8}, isnull::Bool)
    ints = h.ints
    colidx = _ekw_entry_checks(segdsc, coldsc, isnull)
    lnb = max(1, length(cval))
    strlen = coldsc[_EK_LENIDX] != -1 ? min(Int(coldsc[_EK_LENIDX]), lnb) : lnb
    ptrloc = recptr + _EK_DPTBAS + colidx
    _ekw_empty_check(h, segdsc, colidx, ptrloc)
    if isnull
        ints[ptrloc] = _EK_NULL
    else
        n = strlen
        pos = 1
        lastw = Int(segdsc[_ekw_lcw(1)])
        p = Int(segdsc[_ekw_lcp(1)])
        pcount = 0
        pbase = 0
        while n > 0
            # invariant: the characters before `pos` are in the file; `lastw` characters
            # of page `p` are in use
            if lastw < 1009
                pcount += 1
                local datptr::Int
                if pcount == 1
                    pbase = _ek_cbase(p)
                    datptr = pbase + lastw + 1
                    h.ints[ptrloc] = datptr
                    _ekw_sei(h, datptr, strlen)
                    datptr += 5
                    lastw += 5
                else
                    datptr = pbase + 1
                end
                nwrite = min(_EK_CPSIZE - lastw, n)
                _ekw_putc!(h.chars, datptr, cval, pos, nwrite)
                n -= nwrite
                pos += nwrite
                _ekw_slnk(h, 1, p, _ekw_glnk(h, 1, p) + 1)
                lastw += nwrite
                segdsc[_ekw_lcw(1)] = lastw
            else
                p2, pbase = _ekw_aps(h, segdsc, 1)
                pcount > 0 && _ekw_sfwd(h, 1, p, p2)
                p = p2
                lastw = 0
                segdsc[_ekw_lcp(1)] = p
                segdsc[_ekw_lcw(1)] = lastw
                _ekw_slnk(h, 1, p, 0)
            end
        end
    end
    _ekw_put_segdsc(h, segdsc)
    _ekw_indexed(segdsc, coldsc) && _ekw_ixin(h, segdsc, coldsc, cval, recptr, isnull)
    return nothing
end

"""
    _ekw_ad45(h, segdsc, coldsc, recptr, vals, isnull) -> Nothing

Add the array `vals` (or null) of an integer, double or time array column (class 4 or
5) to the record at `recptr` (`zzekad04`, `zzekad05`). The file takes the element
count and then the elements; the entry starts on the last data page of its type if
that has two free words, and continues on new pages that the forward pointers link.
Preconditions: valid column position (`:invalidindex`), null allowed
(`:badattribute`), `length(vals) >= 1` and equal to a fixed entry size
(`:invalidcount`), entry not written yet (`:nonemptyentry`).
"""
function _ekw_ad45(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Int,
        vals::AbstractVector{<:Real}, isnull::Bool)
    t = _ekw_ntype(coldsc)
    colidx = _ekw_entry_checks(segdsc, coldsc, isnull)
    nvals = length(vals)
    _ekw_count_checks(segdsc, coldsc, nvals)
    ptrloc = recptr + _EK_DPTBAS + colidx
    _ekw_empty_check(h, segdsc, colidx, ptrloc)
    if isnull
        h.ints[ptrloc] = _EK_NULL
    else
        dsize = _ekw_DSIZE[t]
        lastw = Int(segdsc[_ekw_lcw(t)])
        room = dsize - lastw
        remain = nvals
        start = firstindex(vals)
        fstpag = true
        p = 0
        while remain > 0
            # invariant: `nvals - remain` elements are in the file; `room` words of the
            # last data page are free
            if room >= 2
                p = Int(segdsc[_ekw_lcp(t)])
                datptr = _ekw_base(t, p) + lastw + 1
                if fstpag
                    h.ints[ptrloc] = datptr
                    _ekw_putn(h, t, datptr, nvals)
                    room -= 1
                    datptr += 1
                end
                nwrite = min(remain, room)
                if t == 2
                    copyto!(h.dps, datptr, vals, start, nwrite)
                else
                    copyto!(h.ints, datptr, vals, start, nwrite)
                end
                remain -= nwrite
                room -= nwrite
                start += nwrite
                _ekw_slnk(h, t, p, _ekw_glnk(h, t, p) + 1)
                if fstpag
                    segdsc[_ekw_lcw(t)] = lastw + 1 + nwrite
                    fstpag = false
                else
                    segdsc[_ekw_lcw(t)] = lastw + nwrite
                end
            else
                p2, _ = _ekw_aps(h, segdsc, t)
                fstpag || _ekw_sfwd(h, t, p, p2)
                p = p2
                lastw = 0
                segdsc[_ekw_lcp(t)] = p
                segdsc[_ekw_lcw(t)] = lastw
                room = dsize
                _ekw_slnk(h, t, p, 0)
            end
        end
    end
    _ekw_put_segdsc(h, segdsc)
    return nothing
end

"""
    _ekw_ad06(h, segdsc, coldsc, recptr, cvals, isnull) -> Nothing

Add the array of strings `cvals` (or null) of a character array column (class 6) to the
record at `recptr` (`zzekad06`). The file takes a five-character element count and
then every string with exactly the string length of the column (cut or continued with
blanks). Preconditions as `_ekw_ad45`.
"""
function _ekw_ad06(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Int,
        cvals::Vector{Vector{UInt8}}, isnull::Bool)
    colidx = _ekw_entry_checks(segdsc, coldsc, isnull)
    nvals = length(cvals)
    _ekw_count_checks(segdsc, coldsc, nvals)
    ptrloc = recptr + _EK_DPTBAS + colidx
    _ekw_empty_check(h, segdsc, colidx, ptrloc)
    if isnull
        h.ints[ptrloc] = _EK_NULL
    else
        strlen = Int(coldsc[_EK_LENIDX])
        lastw = Int(segdsc[_ekw_lcw(1)])
        room = _EK_CPSIZE - lastw
        fstpag = true
        p = 0
        datptr = 0
        if lastw < _EK_CPSIZE
            p = Int(segdsc[_ekw_lcp(1)])
            datptr = _ek_cbase(p) + lastw + 1
        end
        for cval in cvals
            pos = 0
            remain = strlen
            while remain > 0
                # invariant: `pos` characters of the string are in the file; the next
                # one goes to `datptr`, where `room` characters of the page are free
                if room >= (fstpag ? 6 : 1)
                    if fstpag
                        h.ints[ptrloc] = datptr
                        _ekw_sei(h, datptr, nvals)
                        room -= 5
                        datptr += 5
                        _ekw_slnk(h, 1, p, _ekw_glnk(h, 1, p) + 1)
                    end
                    nwrite = min(remain, room)
                    _ekw_putc!(h.chars, datptr, cval, pos + 1, nwrite)
                    pos += nwrite
                    datptr += nwrite
                    remain -= nwrite
                    room -= nwrite
                    lastw += fstpag ? 5 + nwrite : nwrite
                    fstpag = false
                    segdsc[_ekw_lcw(1)] = lastw
                else
                    p2, pbase = _ekw_aps(h, segdsc, 1)
                    fstpag || _ekw_sfwd(h, 1, p, p2)
                    p = p2
                    lastw = 0
                    segdsc[_ekw_lcp(1)] = p
                    segdsc[_ekw_lcw(1)] = lastw
                    room = _EK_CPSIZE
                    datptr = pbase + 1
                    _ekw_slnk(h, 1, p, fstpag ? 0 : 1)
                end
            end
        end
    end
    _ekw_put_segdsc(h, segdsc)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Delete one entry (classes 1 to 6)
# ---------------------------------------------------------------------------------------

@noinline _ekw_badptr(h::DASHandle, segdsc, coldsc) = _ekw_bug("Data pointer is corrupted. " *
    "SEGNO = $(segdsc[_EK_SNOIDX]); COLIDX =  $(coldsc[_EK_ORDIDX]); EK = $(h.path)")

"""
    _ekw_de12(h, segdsc, coldsc, recptr) -> Nothing

Delete the entry of a scalar integer, double or time column (class 1 or 2) from the
record at `recptr` (`zzekde01`, `zzekde02`): the record leaves the index of the column,
the data page loses a link (and is freed with its last link), and the data pointer
becomes "not initialized". An entry that was never written stays as it is.
"""
function _ekw_de12(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Int)
    t = _ekw_ntype(coldsc)
    ptrloc = recptr + _EK_DPTBAS + Int(coldsc[_EK_ORDIDX])
    datptr = Int(h.ints[ptrloc])
    if datptr > 0
        coldsc[_EK_IXTIDX] != -1 && _ekw_ixdl(h, segdsc, coldsc, recptr)
        p, _ = _ekw_pgpg(t, datptr)
        _ekw_unlink(h, segdsc, t, p)
        h.ints[ptrloc] = _EK_UNINIT
    elseif datptr == _EK_NULL
        coldsc[_EK_IXTIDX] != -1 && _ekw_ixdl(h, segdsc, coldsc, recptr)
        h.ints[ptrloc] = _EK_UNINIT
    elseif datptr != _EK_UNINIT
        _ekw_badptr(h, segdsc, coldsc)
    end
    return nothing
end

"""
    _ekw_unlink_chain(h, segdsc, t, datptr, nwords) -> Nothing

Take one link from every data page that CSPICE visits for an entry of `nwords` data
words whose count is at `datptr` (`zzekde03` to `zzekde06`). CSPICE takes the words
of the first page as `base + size - datptr`; for a character entry that is four more
than the page holds after the count, so the link of a last page that holds four
characters or less of the entry is not taken. The same count is used here.

With `fastrec` only the first page loses a link, and no page is freed. That is for the
array entries of a record that the fast writers made: `zzekac04` to `zzekac06` give a
page one link for every entry that starts on it and none for an entry that continues
on it. So the link that CSPICE takes from the other pages belongs to another entry,
and a page whose links are all gone can still hold the end of an entry that started
on the page before. Such pages stay in the segment.
"""
function _ekw_unlink_chain(h::DASHandle, segdsc::Vector{Int32}, t::Int, datptr::Int,
        nwords::Int, fastrec::Bool=false)
    p, base = _ekw_pgpg(t, datptr)
    if fastrec
        nlinks = _ekw_glnk(h, t, p)
        nlinks > 1 && _ekw_slnk(h, t, p, nlinks - 1)
        return nothing
    end
    next = _ekw_gfwd(h, t, p)
    _ekw_unlink(h, segdsc, t, p)
    nseen = min(nwords, base + _ekw_DSIZE[t] - datptr)
    while nseen < nwords                 # invariant: the pages that hold `nseen` words lost a link
        p = next
        next = _ekw_gfwd(h, t, p)
        _ekw_unlink(h, segdsc, t, p)
        nseen = min(nwords, nseen + _ekw_DSIZE[t])
    end
    return nothing
end

"""
    _ekw_fast_record(h, recptr) -> Bool

True if the record at `recptr` is in a record page that the fast writers made. Such a
page has the link count 0 (`zzekff01` writes it so), and a page that `ekinsr` made has
a link for every record that it holds, so it never has the count 0 while a record is
in it.
"""
@inline _ekw_fast_record(h::DASHandle, recptr::Integer) =
    _ekw_glnk(h, 3, _ekw_pgpg(3, Int(recptr) + 1)[1]) == 0

"""
    _ekw_de03(h, segdsc, coldsc, recptr) -> Nothing

Delete the entry of a scalar character column (class 3) from the record at `recptr`
(`zzekde03`): index, links of the pages of the string, data pointer.
"""
function _ekw_de03(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Int)
    ptrloc = recptr + _EK_DPTBAS + Int(coldsc[_EK_ORDIDX])
    datptr = Int(h.ints[ptrloc])
    if datptr > 0
        coldsc[_EK_IXTIDX] != -1 && _ekw_ixdl(h, segdsc, coldsc, recptr)
        nchars = _ek_gei(h.chars, datptr)
        h.ints[ptrloc] = _EK_UNINIT
        _ekw_unlink_chain(h, segdsc, 1, datptr, nchars)
    elseif datptr == _EK_NULL
        coldsc[_EK_IXTIDX] != -1 && _ekw_ixdl(h, segdsc, coldsc, recptr)
        h.ints[ptrloc] = _EK_UNINIT
    elseif datptr != _EK_UNINIT
        _ekw_badptr(h, segdsc, coldsc)
    end
    return nothing
end

"""
    _ekw_de456(h, segdsc, coldsc, recptr) -> Nothing

Delete the entry of an array column (class 4, 5 or 6) from the record at `recptr`
(`zzekde04`, `zzekde05`, `zzekde06`): links of the pages of the entry, data pointer.
CSPICE sets the status word of the record to 2 ("updated") in these three routines;
so does this function.

For a record that the fast writers made (its record page has the link count 0, see
[`_ekw_fast_record`](@ref)) the first page of the entry loses a link and no page is
freed. CSPICE takes a link from every page of the entry and frees a page without
links; after a fast load that frees pages that other entries still use.
"""
function _ekw_de456(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Int)
    class = coldsc[_EK_CLSIDX]
    t = class == 4 ? 3 : class == 5 ? 2 : 1
    ptrloc = recptr + _EK_DPTBAS + Int(coldsc[_EK_ORDIDX])
    datptr = Int(h.ints[ptrloc])
    if datptr > 0
        nwords = t == 3 ? Int(h.ints[datptr]) : t == 2 ? round(Int, h.dps[datptr]) :
                 Int(coldsc[_EK_LENIDX]) * _ek_gei(h.chars, datptr)
        h.ints[ptrloc] = _EK_UNINIT
        _ekw_unlink_chain(h, segdsc, t, datptr, nwords, _ekw_fast_record(h, recptr))
    elseif datptr == _EK_NULL
        h.ints[ptrloc] = _EK_UNINIT
    elseif datptr != _EK_UNINIT
        _ekw_badptr(h, segdsc, coldsc)
    end
    h.ints[recptr+1] = 2
    return nothing
end

function _ekw_de(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Int)
    class = coldsc[_EK_CLSIDX]
    if class == 1 || class == 2
        _ekw_de12(h, segdsc, coldsc, recptr)
    elseif class == 3
        _ekw_de03(h, segdsc, coldsc, recptr)
    elseif 4 <= class <= 6
        _ekw_de456(h, segdsc, coldsc, recptr)
    else
        throw(SpiceError(:noclass, "Class $class from input column descriptor is not " *
            "supported.  COLUMN = $(_ekw_cnam(h, coldsc)); SEGNO = $(segdsc[_EK_SNOIDX]); " *
            "EK = $(h.path)."))
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Update one entry (classes 1 to 6)
# ---------------------------------------------------------------------------------------

"""
    _ekw_ue12(h, segdsc, coldsc, recptr, val, isnull) -> Nothing

Replace the entry of a scalar integer, double or time column (class 1 or 2) in the
record at `recptr` with `val` or null (`zzekue01`, `zzekue02`). A value replaces a
value in place. As in CSPICE, a value that becomes null gives up its link of the data
page but the page is not freed, and that case does not check that the column allows
nulls.
"""
function _ekw_ue12(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32}, recptr::Int,
        val::Real, isnull::Bool)
    t = _ekw_ntype(coldsc)
    ptrloc = recptr + _EK_DPTBAS + Int(coldsc[_EK_ORDIDX])
    datptr = Int(h.ints[ptrloc])
    if datptr > 0
        if _ekw_indexed(segdsc, coldsc)
            _ekw_ixdl(h, segdsc, coldsc, recptr)
            _ekw_ixin(h, segdsc, coldsc, val, recptr, isnull)
        end
        if isnull
            p, _ = _ekw_pgpg(t, datptr)
            _ekw_slnk(h, t, p, _ekw_glnk(h, t, p) - 1)
            h.ints[ptrloc] = _EK_NULL
        else
            _ekw_putn(h, t, datptr, val)
        end
    elseif datptr == _EK_NULL
        if !isnull
            _ekw_indexed(segdsc, coldsc) && _ekw_ixdl(h, segdsc, coldsc, recptr)
            h.ints[ptrloc] = _EK_UNINIT
            _ekw_ad12(h, segdsc, coldsc, recptr, val, isnull)
        end
    elseif datptr == _EK_UNINIT || datptr == _EK_NOBACK
        _ekw_ad12(h, segdsc, coldsc, recptr, val, isnull)
    else
        _ekw_badptr(h, segdsc, coldsc)
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Order vectors of the fast writers
# ---------------------------------------------------------------------------------------

"""
    _ekw_order(n, nullok, nlflgs, inorder) -> Vector{Int32}

The order vector of `n` column entries, as `zzekordc`, `zzekordd` and `zzekordi` make
it with a Shell sort: `iorder[k]` is the number of the entry at rank `k`. Null entries
(if `nullok`) come first, in the order of their numbers; `inorder(a, b)` tells for two
entries that are not null that `a` stays before `b`. The sequence of the swaps is that
of CSPICE, so the result is the same also where the sort is not stable.
"""
function _ekw_order(n::Int, nullok::Bool, nlflgs::AbstractVector{Bool}, inorder::F) where {F}
    iorder = collect(Int32(1):Int32(n))
    gap = n ÷ 2
    while gap > 0
        for i in gap+1:n
            j = i - gap
            while j > 0
                # invariant: the entries at i, i - gap, ... down to j + gap are in order
                jg = j + gap
                a = iorder[j]
                b = iorder[jg]
                ok = if !nullok
                    inorder(a, b)
                else
                    na = nlflgs[a]
                    nb = nlflgs[b]
                    (na && !nb) || (na && nb && a < b) || (!(na || nb) && inorder(a, b))
                end
                if ok
                    j = 0
                else
                    iorder[j] = b
                    iorder[jg] = a
                end
                j -= gap
            end
        end
        gap ÷= 2
    end
    return iorder
end

"Order vector of numbers (`zzekordd`, `zzekordi`): by value, then by entry number."
_ekw_order_n(vals::AbstractVector, n::Int, nullok::Bool, nlflgs) = _ekw_order(n, nullok, nlflgs,
    (a, b) -> vals[a] < vals[b] || (vals[a] == vals[b] && a < b))

"Order vector of strings (`zzekordc`): equal strings are not swapped."
_ekw_order_c(vals::Vector{Vector{UInt8}}, n::Int, nullok::Bool, nlflgs) =
    _ekw_order(n, nullok, nlflgs, (a, b) -> _ekw_cmpc(vals[a], vals[b]) <= 0)

# ---------------------------------------------------------------------------------------
# Fast writers: one whole column (classes 1 to 9)
# ---------------------------------------------------------------------------------------

@inline _ekw_nnull(nlflgs, nrows::Int) = count(i -> nlflgs[i], 1:nrows)

"""
    _ekw_fast_index(h, segdsc, coldsc, iorder, rcptrs) -> Nothing

Make the index tree of a fast-loaded column of class 1, 2 or 3 from its order vector
(the end of `zzekac01`, `zzekac02`, `zzekac03`): the tree holds the record pointers in
the order of the entries, and its root page goes into the column descriptor.
"""
function _ekw_fast_index(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32},
        iorder::Vector{Int32}, rcptrs::AbstractVector{Int32})
    for i in eachindex(iorder)
        iorder[i] = rcptrs[iorder[i]]
    end
    tree = _ekw_trit(h)
    _ekw_tr1s(h, tree, iorder)
    dscbas = Int(segdsc[_EK_IMDIDX]) + _EK_SDSCSZ + (Int(coldsc[_EK_ORDIDX]) - 1) * _EK_CDSCSZ
    h.ints[dscbas+_EK_IXPIDX] = tree
    return nothing
end

@inline function _ekw_fast_indexed(h::DASHandle, coldsc::Vector{Int32})
    idxtyp = coldsc[_EK_IXTIDX]
    idxtyp == -1 && return false
    idxtyp == 1 || throw(SpiceError(:unrecognizedtype, "Index type code $idxtyp found in " *
        "descriptor for column $(_ekw_cnam(h, coldsc)).  Code should be 1."))
    return true
end

"""
    _ekw_ac12(h, segdsc, coldsc, vals, nlflgs, rcptrs, stack) -> Nothing

Store a whole scalar integer, double or time column (class 1 or 2) of a segment in
fast load (`zzekac01`, `zzekac02`). The values that are not null fill new data pages
in row order; each page gets the number of its values as link count. The column
position and then the data pointer of every row (null: -2) are pushed on `stack`. An
indexed column gets its index tree.
"""
function _ekw_ac12(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32},
        vals::AbstractVector{<:Real}, nlflgs::AbstractVector{Bool},
        rcptrs::AbstractVector{Int32}, stack::Vector{Int32})
    t = _ekw_ntype(coldsc)
    dsize = _ekw_DSIZE[t]
    nullok = coldsc[_EK_NFLIDX] != -1
    indexd = _ekw_fast_indexed(h, coldsc)
    push!(stack, coldsc[_EK_ORDIDX])
    nrows = Int(segdsc[_EK_NRIDX])
    ndata = nullok ? nrows - _ekw_nnull(nlflgs, nrows) : nrows
    page = t == 2 ? zeros(Float64, 128) : zeros(Int32, 256)
    p = 0
    pbase = 0
    ndata > 0 && ((p, pbase) = _ekw_aps(h, segdsc, t))
    dst = t == 2 ? h.dps : h.ints
    to = 1
    nwrite = 0
    n = 0
    for from in 1:nrows
        # invariant: `nwrite` values are stored; `n` of them are in the page buffer
        if nullok && nlflgs[from]
            push!(stack, _EK_NULL)
        else
            push!(stack, to + pbase)
            page[to] = vals[from]
            to += 1
            nwrite += 1
            n += 1
        end
        if n == dsize || (nwrite == ndata && ndata != 0)
            page[dsize+2] = n
            copyto!(dst, pbase + 1, page, 1, length(page))
            if nwrite < ndata
                p, pbase = _ekw_aps(h, segdsc, t)
                dst = t == 2 ? h.dps : h.ints
                fill!(page, 0)
                n = 0
                to = 1
            end
        end
    end
    indexd && _ekw_fast_index(h, segdsc, coldsc, _ekw_order_n(vals, nrows, nullok, nlflgs), rcptrs)
    return nothing
end

@inline function _ekw_blank_page!(page::Vector{UInt8})
    fill!(page, UInt8(' '))
    _ekw_enc!(page, 1020, 0)
    return nothing
end

@inline _ekw_inc_link!(page::Vector{UInt8}) = _ekw_enc!(page, 1020, _ek_gei(page, 1020) + 1)

"""
    _ekw_ac03(h, segdsc, coldsc, cvals, nlflgs, rcptrs, stack) -> Nothing

Store a whole scalar character column (class 3) of a segment in fast load
(`zzekac03`). Every string that is not null is a five-character count and its
characters; a string starts on a page only if six characters are free, and a long
string continues on linked pages. CSPICE writes whole page buffers here, so the
forward pointer of a page that no string leaves holds blanks; so it does here.
"""
function _ekw_ac03(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32},
        cvals::Vector{Vector{UInt8}}, nlflgs::AbstractVector{Bool},
        rcptrs::AbstractVector{Int32}, stack::Vector{Int32})
    colwid = Int(coldsc[_EK_LENIDX])
    nullok = coldsc[_EK_NFLIDX] != -1
    indexd = _ekw_fast_indexed(h, coldsc)
    push!(stack, coldsc[_EK_ORDIDX])
    nrows = Int(segdsc[_EK_NRIDX])
    ndata = nullok ? nrows - _ekw_nnull(nlflgs, nrows) : nrows
    page = Vector{UInt8}(undef, 1024)
    p = 0
    pbase = 0
    if ndata > 0
        p, pbase = _ekw_aps(h, segdsc, 1)
        _ekw_blank_page!(page)
    end
    to = 1
    nwrite = 0
    n = 0
    for from in 1:nrows
        # invariant: `nwrite` strings are stored; `n` characters of the page buffer are in use
        if nullok && nlflgs[from]
            push!(stack, _EK_NULL)
        else
            cval = cvals[from]
            lnb = max(1, length(cval))
            strlen = colwid != -1 ? min(lnb, colwid) : lnb
            push!(stack, to + pbase)
            pos = 1
            _ekw_enc!(page, to, strlen)
            n += 5
            to = n + 1
            nchars = strlen
            while nchars > 0
                room = _EK_CPSIZE - n
                if nchars <= room
                    _ekw_putc!(page, to, cval, pos, nchars)
                    n += nchars
                    to = n + 1
                    nchars = 0
                    _ekw_inc_link!(page)
                else
                    _ekw_putc!(page, to, cval, pos, room)
                    pos += room
                    nchars -= room
                    _ekw_inc_link!(page)
                    p2, pbase2 = _ekw_aps(h, segdsc, 1)
                    _ekw_enc!(page, 1015, p2)
                    copyto!(h.chars, pbase + 1, page, 1, 1024)
                    p = p2
                    pbase = pbase2
                    _ekw_blank_page!(page)
                    n = 0
                    to = 1
                end
            end
            nwrite += 1
        end
        if nwrite < ndata
            if _EK_CPSIZE - n < 6
                copyto!(h.chars, pbase + 1, page, 1, 1024)
                p, pbase = _ekw_aps(h, segdsc, 1)
                _ekw_blank_page!(page)
                n = 0
                to = 1
            end
        elseif n > 0
            copyto!(h.chars, pbase + 1, page, 1, 1024)
            n = 0
        end
    end
    indexd && _ekw_fast_index(h, segdsc, coldsc, _ekw_order_c(cvals, nrows, nullok, nlflgs), rcptrs)
    return nothing
end

@inline _ekw_cursiz(size::Int, entszs, row::Int) = size != -1 ? size : Int(entszs[row])

"""
    _ekw_ac45(h, segdsc, coldsc, vals, entszs, nlflgs, stack) -> Nothing

Store a whole integer, double or time array column (class 4 or 5) of a segment in fast
load (`zzekac04`, `zzekac05`). `vals` holds the elements of all rows in row order,
also the (unused) elements of the null rows. An entry is its element count and its
elements; it starts on a page only if two words are free there, and it continues on
linked pages. The link count of a page is the number of entries that start on it.
"""
function _ekw_ac45(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32},
        vals::AbstractVector{<:Real}, entszs::AbstractVector{<:Integer},
        nlflgs::AbstractVector{Bool}, stack::Vector{Int32})
    t = _ekw_ntype(coldsc)
    dsize = _ekw_DSIZE[t]
    size = Int(coldsc[_EK_SIZIDX])
    nullok = coldsc[_EK_NFLIDX] != -1
    push!(stack, coldsc[_EK_ORDIDX])
    nrows = Int(segdsc[_EK_NRIDX])
    ndata = 0
    for i in 1:nrows
        (nullok && nlflgs[i]) || (ndata += _ekw_cursiz(size, entszs, i))
    end
    page = t == 2 ? zeros(Float64, 128) : zeros(Int32, 256)
    p = 0
    pbase = 0
    ndata > 0 && ((p, pbase) = _ekw_aps(h, segdsc, t))
    remain = ndata
    from = 1
    to = 1
    row = 1
    nelt = 1
    n = 0
    nlink = 0
    cursiz = 0
    while row <= nrows
        # invariant: the entries of the rows before `row` and `nelt - 1` elements of
        # this row are stored; `n` words of the page buffer are in use
        newreq = false
        cntinu = false
        if nullok && nlflgs[row]
            from += _ekw_cursiz(size, entszs, row)
            push!(stack, _EK_NULL)
            row += 1
            nelt = 1
        else
            if nelt == 1
                cursiz = _ekw_cursiz(size, entszs, row)
                push!(stack, to + pbase)
                page[to] = cursiz
                to += 1
                n += 1
                nlink += 1
            end
            page[to] = vals[from]
            to += 1
            n += 1
            from += 1
            remain -= 1
            cntinu = nelt < cursiz && n == dsize
            if nelt == cursiz
                remain > 0 && (newreq = n > dsize - 2)
                nelt = 1
                row += 1
            else
                nelt += 1
            end
        end
        if cntinu || newreq || (row > nrows && ndata > 0)
            page[dsize+2] = nlink
            copyto!(t == 2 ? h.dps : h.ints, pbase + 1, page, 1, length(page))
            if remain > 0
                p2, pbase = _ekw_aps(h, segdsc, t)
                fill!(page, 0)
                n = 0
                nlink = 0
                to = 1
                cntinu && _ekw_sfwd(h, t, p, p2)
                p = p2
            end
        end
    end
    return nothing
end

"""
    _ekw_ac06(h, segdsc, coldsc, cvals, entszs, nlflgs, stack) -> Nothing

Store a whole character array column (class 6) of a segment in fast load
(`zzekac06`). `cvals` holds the strings of all rows in row order, also the (unused)
strings of the null rows. An entry is a five-character element count and its strings,
each with exactly the string length of the column; it starts on a page only if six
characters are free there.

`zzekac06` loses the last character of a string when the declared length of its input
strings is less than the string length of the column and a page ends just before that
character. A string here has no declared length (it continues with blanks), so that
case does not exist.
"""
function _ekw_ac06(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32},
        cvals::Vector{Vector{UInt8}}, entszs::AbstractVector{<:Integer},
        nlflgs::AbstractVector{Bool}, stack::Vector{Int32})
    size = Int(coldsc[_EK_SIZIDX])
    strlen = Int(coldsc[_EK_LENIDX])
    nullok = coldsc[_EK_NFLIDX] != -1
    push!(stack, coldsc[_EK_ORDIDX])
    nrows = Int(segdsc[_EK_NRIDX])
    ndata = 0
    for i in 1:nrows
        (nullok && nlflgs[i]) || (ndata += strlen * _ekw_cursiz(size, entszs, i))
    end
    page = fill(UInt8(' '), 1024)
    p = 0
    pbase = 0
    ndata > 0 && ((p, pbase) = _ekw_aps(h, segdsc, 1))
    remain = ndata
    from = 1
    to = 1
    row = 1
    cp = 1
    n = 0
    nchars = 0
    nlink = 0
    curchr = 0
    newent = true
    while row <= nrows
        # invariant: the entries of the rows before `row` and `nchars` characters of
        # this row are stored; `n` characters of the page buffer are in use
        newreq = false
        cntinu = false
        if nullok && nlflgs[row]
            from += _ekw_cursiz(size, entszs, row)
            push!(stack, _EK_NULL)
            row += 1
            newent = true
        else
            if newent
                cursiz = _ekw_cursiz(size, entszs, row)
                curchr = cursiz * strlen
                nchars = 0
                cp = 1
                push!(stack, to + pbase)
                _ekw_enc!(page, to, cursiz)
                to += 5
                n += 5
                nlink += 1
                newent = false
            end
            nwrite = min(curchr - nchars, _EK_CPSIZE - n)
            nw = nwrite
            while nw > 0
                l = min(strlen - cp + 1, nw)
                _ekw_putc!(page, to, cvals[from], cp, l)
                cp += l
                nw -= l
                to += l
                if cp > strlen
                    from += 1
                    cp = 1
                end
            end
            n += nwrite
            remain -= nwrite
            nchars += nwrite
            cntinu = nchars < curchr && n == _EK_CPSIZE
            if nchars == curchr
                remain > 0 && (newreq = n > 1008)
                row += 1
                newent = true
            end
        end
        if cntinu || newreq || (row > nrows && ndata > 0)
            _ekw_enc!(page, 1020, nlink)
            copyto!(h.chars, pbase + 1, page, 1, 1024)
            if remain > 0
                p2, pbase = _ekw_aps(h, segdsc, 1)
                fill!(page, UInt8(' '))
                n = 0
                nlink = 0
                to = 1
                cntinu && _ekw_sfwd(h, 1, p, p2)
                p = p2
            end
        end
    end
    return nothing
end

"""
    _ekw_wpai(h, segdsc, ivals) -> base

Store the integers `ivals` in new adjacent integer data pages of the segment, 254 in a
page, and return the base address of the first page (`zzekwpai`). The page buffer is
not cleared between pages, as in CSPICE, so the unused words of the last page repeat
words of the page before.
"""
function _ekw_wpai(h::DASHandle, segdsc::Vector{Int32}, ivals::AbstractVector{<:Integer})
    nvals = length(ivals)
    p, base = _ekw_acps(h, segdsc, 3, (nvals + 253) ÷ 254)
    page = zeros(Int32, 256)
    to = 1
    for from in 1:nvals                  # invariant: the pages before `p` are written
        page[to] = ivals[from]
        to += 1
        if to > 254 || from == nvals
            page[256] = to - 1
            copyto!(h.ints, _ek_ibase(p) + 1, page, 1, 256)
            p += 1
            to = 1
        end
    end
    return base
end

"""
    _ekw_wpal(h, segdsc, lvals, nvals) -> base

Store the first `nvals` flags of `lvals` as the characters `T` and `F` in new adjacent
character data pages of the segment, 1014 in a page, and return the base address of
the first page (`zzekwpal`). The page buffer is not cleared between pages.
"""
function _ekw_wpal(h::DASHandle, segdsc::Vector{Int32}, lvals::AbstractVector{Bool}, nvals::Int)
    p, base = _ekw_acps(h, segdsc, 1, (nvals + 1013) ÷ 1014)
    page = fill(UInt8(' '), 1024)
    to = 1
    for from in 1:nvals                  # invariant: the pages before `p` are written
        page[to] = lvals[from] ? UInt8('T') : UInt8('F')
        to += 1
        if to > 1014 || from == nvals
            copyto!(h.chars, _ek_cbase(p) + 1, page, 1, 1024)
            _ekw_slnk(h, 1, p, to - 1)
            p += 1
            to = 1
        end
    end
    return base
end

"""
    _ekw_wpac(h, segdsc, l, cvals) -> base

Store the strings `cvals`, each with exactly `l` characters, in new adjacent character
data pages of the segment and return the base address of the first page (`zzekwpac`).
Precondition: `1 <= l <= 1014`, else `SpiceError(:invalidsize)`.
"""
function _ekw_wpac(h::DASHandle, segdsc::Vector{Int32}, l::Int, cvals::Vector{Vector{UInt8}})
    (1 <= l <= _EK_CPSIZE) || throw(SpiceError(:invalidsize,
        "String length $l is just plain wrong."))
    nvals = length(cvals)
    spp = _EK_CPSIZE ÷ l
    p, base = _ekw_acps(h, segdsc, 1, (nvals + spp - 1) ÷ spp)
    page = fill(UInt8(' '), 1024)
    to = 1
    for from in 1:nvals                  # invariant: the pages before `p` are written
        _ekw_putc!(page, to, cvals[from], 1, l)
        to += l
        if to > _EK_CPSIZE - l + 1 || from == nvals
            copyto!(h.chars, _ek_cbase(p) + 1, page, 1, 1024)
            _ekw_slnk(h, 1, p, (to - l) ÷ l)
            p += 1
            to = 1
        end
    end
    return base
end

"""
    _ekw_fixed_finish(h, segdsc, coldsc, datbas, iorder, nlflgs) -> Nothing

The common end of `zzekac07`, `zzekac08` and `zzekac09`: the base of the data goes
into the metadata of the column; an indexed column gets its order vector (row numbers,
index type 2) and a column that allows nulls its page array of null flags.
"""
function _ekw_fixed_finish(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32},
        datbas::Int, iorder, nlflgs::AbstractVector{Bool})
    nrows = Int(segdsc[_EK_NRIDX])
    cmbase = Int(coldsc[_EK_METIDX])
    h.ints[cmbase+1] = datbas
    if iorder !== nothing
        idxbas = _ekw_wpai(h, segdsc, iorder)
        dscbas = Int(segdsc[_EK_IMDIDX]) + _EK_SDSCSZ + (Int(coldsc[_EK_ORDIDX]) - 1) * _EK_CDSCSZ
        h.ints[dscbas+_EK_IXPIDX] = idxbas
        h.ints[dscbas+_EK_IXTIDX] = 2
    end
    if coldsc[_EK_NFLIDX] != -1
        h.ints[cmbase+2] = _ekw_wpal(h, segdsc, nlflgs, nrows)
    end
    return nothing
end

"""
    _ekw_ac78(h, segdsc, coldsc, vals, nlflgs) -> Nothing

Store a whole fixed-count integer, double or time column (class 7 or 8) of a type 2
segment (`zzekac07`, `zzekac08`): row `r` is word `r` of adjacent data pages. The slot
of a null row keeps what the page buffer held, as in CSPICE (zero, or the value of the
row one page before).
"""
function _ekw_ac78(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32},
        vals::AbstractVector{<:Real}, nlflgs::AbstractVector{Bool})
    t = _ekw_ntype(coldsc)
    dsize = _ekw_DSIZE[t]
    nullok = coldsc[_EK_NFLIDX] != -1
    indexd = coldsc[_EK_IXTIDX] != -1
    nrows = Int(segdsc[_EK_NRIDX])
    p, datbas = _ekw_acps(h, segdsc, t, (nrows + dsize - 1) ÷ dsize)
    page = t == 2 ? zeros(Float64, 128) : zeros(Int32, 256)
    dst = t == 2 ? h.dps : h.ints
    to = 1
    for from in 1:nrows                  # invariant: the pages before `p` are written
        (!nullok || !nlflgs[from]) && (page[to] = vals[from])
        to += 1
        if to > dsize || from == nrows
            page[dsize+2] = to - 1
            copyto!(dst, _ekw_base(t, p) + 1, page, 1, length(page))
            p += 1
            to = 1
        end
    end
    iorder = indexd ? _ekw_order_n(vals, nrows, nullok, nlflgs) : nothing
    _ekw_fixed_finish(h, segdsc, coldsc, datbas, iorder, nlflgs)
    return nothing
end

"""
    _ekw_ac09(h, segdsc, coldsc, cvals, nlflgs) -> Nothing

Store a whole fixed-count character column (class 9) of a type 2 segment
(`zzekac09`): every row takes exactly the string length `l` of the column, and a page
holds `1014 ÷ l` rows. Precondition: `1 <= l <= 1014`, else `:invalidsize` (CSPICE
divides by zero for `l = 0`, and it also refuses a column whose string length is
greater than the declared length of the input strings, which a string here does not
have). CSPICE writes whole page buffers, so the forward pointers of these pages hold
blanks.
"""
function _ekw_ac09(h::DASHandle, segdsc::Vector{Int32}, coldsc::Vector{Int32},
        cvals::Vector{Vector{UInt8}}, nlflgs::AbstractVector{Bool})
    l = Int(coldsc[_EK_LENIDX])
    nullok = coldsc[_EK_NFLIDX] != -1
    indexd = coldsc[_EK_IXTIDX] != -1
    l == -1 && _ekw_bug("Column $(_ekw_cnam(h, coldsc)) has variable string length; class " *
        "9 supports fixed-length strings only.")
    (1 <= l <= _EK_CPSIZE) || throw(SpiceError(:invalidsize,
        "String length $l is just plain wrong."))
    spp = _EK_CPSIZE ÷ l
    nrows = Int(segdsc[_EK_NRIDX])
    p, datbas = _ekw_acps(h, segdsc, 1, (nrows + spp - 1) ÷ spp)
    page = fill(UInt8(' '), 1024)
    to = 1
    for from in 1:nrows                  # invariant: the pages before `p` are written
        (!nullok || !nlflgs[from]) && _ekw_putc!(page, to, cvals[from], 1, l)
        to += l
        if to > _EK_CPSIZE - l + 1 || from == nrows
            copyto!(h.chars, _ek_cbase(p) + 1, page, 1, 1024)
            _ekw_slnk(h, 1, p, (to - l) ÷ l)
            p += 1
            to = 1
        end
    end
    iorder = indexd ? _ekw_order_c(cvals, nrows, nullok, nlflgs) : nothing
    _ekw_fixed_finish(h, segdsc, coldsc, datbas, iorder, nlflgs)
    return nothing
end
