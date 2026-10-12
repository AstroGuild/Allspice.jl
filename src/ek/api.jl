# Public EK calls: direct access to an EK file through an owned handle (`ekopr`,
# `eknseg`, `ekssum`, `ekrcec`, `ekrced`, `ekrcei`), and the calls on loaded EK data
# (`ekntab`, `ektnam`, `ekccnt`, `ekcii`, `ekpsel`, `ekfind`, `eknelt`, `ekgc`, `ekgd`,
# `ekgi`).
#
# CSPICE keeps the table of loaded EK files and the result of the last `ekfind` in
# saved variables. Here the loaded EK data belong to the `KernelSet` snapshot, and
# `ekfind` returns an `EKResult` object that the fetch calls take as first argument.

# ---------------------------------------------------------------------------------------
# Direct access through a handle
# ---------------------------------------------------------------------------------------

"""
    ekopr(path) -> DASHandle

Open the EK file `path` for read access, as `ekopr_c`. The handle is a `DASHandle` that
the caller owns; close it with [`ekcls`](@ref). Contract: `path` names a DAS file that
holds a paged EK of architecture 8. Throws `SpiceError` with the kinds of `dasopr`,
and `:wrongarchitecture` or `:invalidformat` for a DAS file that is not such an EK.
"""
function ekopr(path::AbstractString)
    h = _das_load(path, "ekopr", false)
    _ek_check(h, h.path, "ekopr")
    return h
end

"""
    ekcls(handle) -> Nothing

Close the EK file of `handle`, as `ekcls_c`. For a write handle all data are on disk
when the call returns. Closing a closed handle does nothing.
"""
ekcls(h::DASHandle) = dascls(h)

"Checks of every direct access call: open handle that holds a paged EK."
@inline function _ek_check_read(h::DASHandle, who::AbstractString)
    _das_check_open(h, who)
    _ek_check(h, h.path, who)
    return nothing
end

"""
    eknseg(handle) -> Int

Number of segments of the EK file of `handle`, as `eknseg_c`. Throws
`SpiceError(:dasnosuchhandle)` for a closed handle and `:wrongarchitecture` or
`:invalidformat` for a file that is not a paged EK.
"""
function eknseg(h::DASHandle)
    _ek_check_read(h, "eknseg")
    return _ek_nseg(h)
end

"""
    ekssum(handle, segno) -> EKSegmentSummary

Summary of segment `segno` (from 1) of the EK file of `handle`, as `ekssum_c`: table
name, number of rows, and the name and attributes of every column. Precondition:
`1 <= segno <= eknseg(handle)`, else `SpiceError(:indexoutofrange)`.
"""
function ekssum(h::DASHandle, segno::Integer)
    _ek_check_read(h, "ekssum")
    nseg = _ek_nseg(h)
    (1 <= segno <= nseg) || throw(SpiceError(:indexoutofrange,
        "Segment index was $segno; valid range is 1:$nseg"))
    _, base = _ek_mloc(h, segno)
    segdsc = view(h.ints, base+1:base+_EK_SDSCSZ)
    ncols = Int(segdsc[_EK_NCIDX])
    tnb = Int(segdsc[_EK_TNMIDX])
    cols = _EKColumn[_ek_column(h, base, i) for i in 1:ncols]
    return EKSegmentSummary(_ek_rstr(h.chars, tnb + 1, tnb + _EK_TNAMSZ),
        Int(segdsc[_EK_NRIDX]), ncols, String[c.name for c in cols],
        EKColumnAttr[_ek_attr(c) for c in cols])
end

"""
    _ek_locate(h, who, segno, recno, column, types) -> (col, recptr)

The common start of `ekrcec`, `ekrced` and `ekrcei`: the descriptor of `column` in
segment `segno` and the record pointer of row `recno`. Throws `SpiceError`:
`:invalidindex` for a bad segment number, `:bug` for an unknown column (as CSPICE),
`:wrongdatatype` unless the data type code of the column is in `types`,
`:indexoutofrange` for a row number outside `1:nrows`.
"""
function _ek_locate(h::DASHandle, who::AbstractString, segno::Integer, recno::Integer,
        column::AbstractString, types::Tuple)
    _ek_check_read(h, who)
    segdsc = _ek_sdsc(h, segno)
    col = _ek_find_column(h, segdsc, column)
    col.dtype in types || throw(SpiceError(:wrongdatatype, "Column $(col.name) is of type " *
        "$(col.dtype); $(uppercase(who)) only works with " *
        "$(types[1] == _EK_CHR ? "character" : types[1] == _EK_INT ? "integer" : "d.p. or TIME") " *
        "columns.  RECNO = $recno; SEGNO = $segno; EK = $(h.path)."))
    nrows = Int(segdsc[_EK_NRIDX])
    (1 <= recno <= nrows) || throw(SpiceError(:indexoutofrange,
        "$who: record number = $recno; valid range = 1:$nrows. SEGNO = $segno; EK = $(h.path)"))
    return col, _ek_rplk(h, segdsc, recno)
end

"""
    ekrcec(handle, segno, recno, column) -> (cvals::Vector{String}, isnull::Bool)

Read the entry of the character column `column` in row `recno` of segment `segno`, as
`ekrcec_c`: all strings of the entry, without trailing blanks (the C interface removes
them too). A null entry gives an empty vector and `true`. Preconditions (each throws
`SpiceError`): valid segment number (`:invalidindex`); the column exists (`:bug`, as
CSPICE) and has character type (`:wrongdatatype`); `1 <= recno <= nrows`
(`:indexoutofrange`); the entry was written (`:uninitialized`).
"""
function ekrcec(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString)
    col, recptr = _ek_locate(h, "ekrcec", segno, recno, column, (_EK_CHR,))
    vals, isnull = _ek_entry_c(h, col, recptr)
    return String[String(rstrip(v, ' ')) for v in vals], isnull
end

"""
    ekrced(handle, segno, recno, column) -> (dvals::Vector{Float64}, isnull::Bool)

Read the entry of the double precision or time column `column` in row `recno` of
segment `segno`, as `ekrced_c`. A null entry gives an empty vector and `true`.
Preconditions as [`ekrcec`](@ref), with `:uninitializedvalue` for an entry that was
never written.
"""
function ekrced(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString)
    col, recptr = _ek_locate(h, "ekrced", segno, recno, column, (_EK_DP, _EK_TIME))
    return _ek_entry_d(h, col, recptr)
end

"""
    ekrcei(handle, segno, recno, column) -> (ivals::Vector{Int}, isnull::Bool)

Read the entry of the integer column `column` in row `recno` of segment `segno`, as
`ekrcei_c`. A null entry gives an empty vector and `true`. Preconditions as
[`ekrced`](@ref).
"""
function ekrcei(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString)
    col, recptr = _ek_locate(h, "ekrcei", segno, recno, column, (_EK_INT,))
    return _ek_entry_i(h, col, recptr)
end
