# Reading column entries of an EK segment, for the nine column classes:
#
#   class   type          entry                        segment type
#   1       INT           one value                    1
#   2       DP, TIME      one value                    1
#   3       CHR           one string                   1
#   4       INT           array                        1
#   5       DP, TIME      array                        1
#   6       CHR           array of strings             1
#   7       INT           one value, fixed layout      2 (fast writers)
#   8       DP, TIME      one value, fixed layout      2
#   9       CHR           one string, fixed layout     2
#
# In a segment of type 1 a row is a record: an array of integers at the record pointer
# `recptr`, with the status at `recptr + 1`, a companion pointer at `recptr + 2`, and
# then one data pointer for each column. A data pointer is the DAS address of the entry
# (of the type of the column), or a negative code: null, uninitialized. An entry that
# does not fit in the rest of its data page continues at the start of the page that the
# forward pointer of the page names; the words of an entry form a "stream".
#
# In a segment of type 2 the entries of a column are stored in row order in adjacent
# pages, and the null flags (one character for a row) in adjacent character pages.
#
# The functions take a store (fields `chars`, `dps`, `ints`), a decoded column
# descriptor and a record pointer (for classes 7 to 9 the row number). They read only.

const _EK_DPTBAS = 2              # the data pointers of a record follow two integers

"Data pointer of the column `col` in the record at `recptr`."
@inline _ek_datptr(s, col::_EKColumn, recptr::Integer) = Int(s.ints[Int(recptr)+_EK_DPTBAS+col.ord])

"Throws the error of CSPICE for a data pointer that is not an address and not null."
@noinline function _ek_badptr(col::_EKColumn, datptr::Int, recptr::Integer)
    if datptr == _EK_UNINIT || datptr == _EK_NOBACK
        kind = (col.cclass == 3 || col.cclass == 6) ? :uninitialized : :uninitializedvalue
        throw(SpiceError(kind, "Attempted to read uninitialized column entry.  " *
            "COLUMN = $(col.name); record pointer = $recptr"))
    end
    throw(SpiceError(:bug, "Data pointer is corrupted. COLUMN = $(col.name); " *
        "record pointer = $recptr"))
end

"Throws `SpiceError(:noclass)` for a column class that the data type does not have."
@noinline _ek_noclass(col::_EKColumn, what::AbstractString) = throw(SpiceError(:noclass,
    "Class $(col.cclass) from input column descriptor is not a supported $what class.  " *
    "COLUMN = $(col.name)."))

# ---------------------------------------------------------------------------------------
# Streams: words of one entry that continue over linked pages
# ---------------------------------------------------------------------------------------

"""
    _ek_stream_i!(out, ints, addr, first, n) -> out

Copy the elements `first:first+n-1` of the integer array entry whose element count is
at the address `addr` into `out[1:n]` (`zzekrd04`). Precondition: the entry has at
least `first + n - 1` elements.
"""
function _ek_stream_i!(out::AbstractVector, ints::AbstractVector{Int32}, addr::Int, first::Int, n::Int)
    base = _ek_ibase((addr - 1) ÷ 256)
    off = addr - base + first            # position of the element in the page, if it is there
    done = 0
    while done < n
        # invariant: `done` elements are copied; the next one is at stream offset `off` of the page
        while off > _EK_IPSIZE
            base = _ek_ibase(ints[base+255])
            off -= _EK_IPSIZE
        end
        take = min(n - done, _EK_IPSIZE - off + 1)
        for k in 0:take-1
            out[done+k+1] = ints[base+off+k]
        end
        done += take
        off += take
    end
    return out
end

"As `_ek_stream_i!` for a double array entry (`zzekrd05`); pointers are stored as doubles."
function _ek_stream_d!(out::AbstractVector, dps::AbstractVector{Float64}, addr::Int, first::Int, n::Int)
    base = _ek_dbase((addr + 127) ÷ 128)
    off = addr - base + first
    done = 0
    while done < n
        # invariant: as in `_ek_stream_i!`
        while off > _EK_DPSIZE
            base = _ek_dbase(round(Int, dps[base+127]))
            off -= _EK_DPSIZE
        end
        take = min(n - done, _EK_DPSIZE - off + 1)
        for k in 0:take-1
            out[done+k+1] = dps[base+off+k]
        end
        done += take
        off += take
    end
    return out
end

"""
    _ek_stream_c(chars, addr, skip, n) -> String

The `n` characters of a character entry that follow the first `skip` characters after
the five-character count at the address `addr` (`zzekrd03`, `zzekrd06`).
"""
function _ek_stream_c(chars::AbstractVector{UInt8}, addr::Int, skip::Int, n::Int)
    out = Vector{UInt8}(undef, n)
    base = _ek_cbase((addr + 1023) ÷ 1024)
    off = addr - base + 5 + skip         # stream offset of the first wanted character
    done = 0
    while done < n
        # invariant: as in `_ek_stream_i!`
        while off > _EK_CPSIZE
            base = _ek_cbase(_ek_gei(chars, base + 1015))
            off -= _EK_CPSIZE
        end
        take = min(n - done, _EK_CPSIZE - off + 1)
        copyto!(out, done + 1, chars, base + off, take)
        done += take
        off += take
    end
    return String(out)
end

# ---------------------------------------------------------------------------------------
# Entry size
# ---------------------------------------------------------------------------------------

"""
    _ek_nelt(s, col, recptr) -> Int

Number of elements of the entry of `col` in the record at `recptr`, as `zzekesiz`: 1 for
a scalar column, the fixed size of an array column, and for an array column of
variable size the stored count (1 for an entry that is null or not initialized).
"""
function _ek_nelt(s, col::_EKColumn, recptr::Integer)
    c = col.cclass
    (c <= 3 || c >= 7) && return 1
    col.size != -1 && return Int(col.size)
    datptr = _ek_datptr(s, col, recptr)
    datptr < 1 && return 1
    c == 4 && return Int(s.ints[datptr])
    c == 5 && return round(Int, s.dps[datptr])
    return _ek_gei(s.chars, datptr)
end

"True if row `recno` of the class 7 to 9 column `col` is null."
@inline function _ek_fixed_null(s, col::_EKColumn, recno::Int)
    col.nullok || return false
    nflbas = Int(s.ints[col.metptr+2])
    q = (recno - 1) ÷ _EK_CPSIZE
    r = recno - q * _EK_CPSIZE
    return s.chars[nflbas+r+q*1024] == UInt8('T')
end

# ---------------------------------------------------------------------------------------
# One element
# ---------------------------------------------------------------------------------------

"""
    _ek_elt_i(s, col, recptr, elt) -> (value, isnull, found)

Element `elt` of the entry of the integer column `col` (class 1, 4 or 7) in the record
at `recptr`, as `zzekrsi`. `found` is false if the entry has no element `elt`; a null
entry gives `(0, true, true)`. Throws `SpiceError` (`:uninitializedvalue`, `:bug`,
`:noclass`).
"""
function _ek_elt_i(s, col::_EKColumn, recptr::Integer, elt::Integer=1)
    c = col.cclass
    if c == 7
        recno = Int(recptr)
        _ek_fixed_null(s, col, recno) && return 0, true, true
        q = (recno - 1) ÷ _EK_IPSIZE
        r = recno - q * _EK_IPSIZE
        return Int(s.ints[Int(s.ints[col.metptr+1])+r+q*256]), false, true
    end
    (c == 1 || c == 4) || _ek_noclass(col, "integer")
    datptr = _ek_datptr(s, col, recptr)
    datptr == _EK_NULL && return 0, true, true
    datptr > 0 || _ek_badptr(col, datptr, recptr)
    c == 1 && return Int(s.ints[datptr]), false, true
    nelt = Int(s.ints[datptr])
    (1 <= elt <= nelt) || return 0, false, false
    v = Ref{Int32}(0)
    _ek_stream_i!(_EKRefVec(v), s.ints, datptr, Int(elt), 1)
    return Int(v[]), false, true
end

"A one-element vector on a `Ref`, so that a stream read of one element does not allocate."
struct _EKRefVec{T} <: AbstractVector{T}
    r::Base.RefValue{T}
end
Base.size(::_EKRefVec) = (1,)
Base.setindex!(v::_EKRefVec, x, ::Int) = (v.r[] = x)
Base.getindex(v::_EKRefVec, ::Int) = v.r[]

"""
    _ek_elt_d(s, col, recptr, elt) -> (value, isnull, found)

Element `elt` of the entry of the double or time column `col` (class 2, 5 or 8), as
`zzekrsd`. Results and errors as [`_ek_elt_i`](@ref).
"""
function _ek_elt_d(s, col::_EKColumn, recptr::Integer, elt::Integer=1)
    c = col.cclass
    if c == 8
        recno = Int(recptr)
        _ek_fixed_null(s, col, recno) && return 0.0, true, true
        q = (recno - 1) ÷ _EK_DPSIZE
        r = recno - q * _EK_DPSIZE
        return s.dps[Int(s.ints[col.metptr+1])+r+q*128], false, true
    end
    (c == 2 || c == 5) || _ek_noclass(col, "double precision")
    datptr = _ek_datptr(s, col, recptr)
    datptr == _EK_NULL && return 0.0, true, true
    datptr > 0 || _ek_badptr(col, datptr, recptr)
    c == 2 && return s.dps[datptr], false, true
    nelt = round(Int, s.dps[datptr])
    (1 <= elt <= nelt) || return 0.0, false, false
    v = Ref{Float64}(0.0)
    _ek_stream_d!(_EKRefVec(v), s.dps, datptr, Int(elt), 1)
    return v[], false, true
end

"""
    _ek_elt_c(s, col, recptr, elt) -> (value, isnull, found)

Element `elt` of the entry of the character column `col` (class 3, 6 or 9), as
`zzekrsc`: the string with all characters that the file holds for it (the stored
length for class 3, the string length of the column for classes 6 and 9, so trailing
blanks can be present). Results and errors as [`_ek_elt_i`](@ref), with
`:uninitialized` for an entry that was never written.
"""
function _ek_elt_c(s, col::_EKColumn, recptr::Integer, elt::Integer=1)
    c = col.cclass
    if c == 9
        recno = Int(recptr)
        _ek_fixed_null(s, col, recno) && return "", true, true
        l = Int(col.strlen)
        spp = _EK_CPSIZE ÷ l
        q = (recno - 1) ÷ spp
        r = recno - q * spp
        a = Int(s.ints[col.metptr+1]) + q * 1024 + (r - 1) * l + 1
        return String(s.chars[a:a+l-1]), false, true
    end
    (c == 3 || c == 6) || _ek_noclass(col, "character")
    datptr = _ek_datptr(s, col, recptr)
    datptr == _EK_NULL && return "", true, true
    datptr > 0 || _ek_badptr(col, datptr, recptr)
    n = _ek_gei(s.chars, datptr)
    c == 3 && return _ek_stream_c(s.chars, datptr, 0, n), false, true
    (1 <= elt <= n) || return "", false, false
    l = Int(col.strlen)
    return _ek_stream_c(s.chars, datptr, l * (Int(elt) - 1), l), false, true
end

# ---------------------------------------------------------------------------------------
# Whole entries
# ---------------------------------------------------------------------------------------

"""
    _ek_entry_i(s, col, recptr) -> (values, isnull)

All elements of the entry of the integer column `col` in the record at `recptr`
(`ekrcei` after the descriptor lookup). A null entry gives an empty vector and `true`.
"""
function _ek_entry_i(s, col::_EKColumn, recptr::Integer)
    if col.cclass == 4
        datptr = _ek_datptr(s, col, recptr)
        datptr == _EK_NULL && return Int[], true
        datptr > 0 || _ek_badptr(col, datptr, recptr)
        n = Int(s.ints[datptr])
        return _ek_stream_i!(Vector{Int}(undef, n), s.ints, datptr, 1, n), false
    end
    v, isnull, _ = _ek_elt_i(s, col, recptr, 1)
    return isnull ? Int[] : Int[v], isnull
end

"As [`_ek_entry_i`](@ref) for a double or time column (`ekrced`)."
function _ek_entry_d(s, col::_EKColumn, recptr::Integer)
    if col.cclass == 5
        datptr = _ek_datptr(s, col, recptr)
        datptr == _EK_NULL && return Float64[], true
        datptr > 0 || _ek_badptr(col, datptr, recptr)
        n = round(Int, s.dps[datptr])
        return _ek_stream_d!(Vector{Float64}(undef, n), s.dps, datptr, 1, n), false
    end
    v, isnull, _ = _ek_elt_d(s, col, recptr, 1)
    return isnull ? Float64[] : Float64[v], isnull
end

"As [`_ek_entry_i`](@ref) for a character column (`ekrcec`); the strings are as stored."
function _ek_entry_c(s, col::_EKColumn, recptr::Integer)
    if col.cclass == 6
        datptr = _ek_datptr(s, col, recptr)
        datptr == _EK_NULL && return String[], true
        datptr > 0 || _ek_badptr(col, datptr, recptr)
        n = _ek_gei(s.chars, datptr)
        l = Int(col.strlen)
        return String[_ek_stream_c(s.chars, datptr, l * (k - 1), l) for k in 1:n], false
    end
    v, isnull, _ = _ek_elt_c(s, col, recptr, 1)
    return isnull ? String[] : String[v], isnull
end
