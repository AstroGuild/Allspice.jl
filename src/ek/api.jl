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

# ---------------------------------------------------------------------------------------
# The catalogue of loaded tables
# ---------------------------------------------------------------------------------------

"""
    ekntab([ctx]) -> Int
    ekntab(ks::KernelSet) -> Int

Number of loaded EK tables, as `ekntab_c`. Postcondition: the result is 0 if no EK file
is loaded.
"""
ekntab(ks::KernelSet) = length(ks.ektables)
ekntab(ctx::SpiceContext) = ekntab(kernelset(ctx))
ekntab() = ekntab(default_context())

"""
    eklef([ctx,] path) -> String

Load the EK file at `path` for queries, as `eklef_c`, and return the handle for `ekuef`
(the path itself, since a context has no file handles). The call is `furnsh` for a file
that is an EK. Throws `SpiceError(:nofile)` if the file does not exist and
`SpiceError(:wrongarchitecture)` if it is not an EK; the context is then unchanged.
"""
function eklef(ctx::SpiceContext, path::AbstractString)::String
    p = String(path)
    isfile(p) || throw(SpiceError(:nofile, "eklef: kernel file not found: $p"))
    kernel_arch(p) === :EK || throw(SpiceError(:wrongarchitecture,
        "eklef: the file $p is not an EK file"))
    furnsh(ctx, p)
    return p
end
eklef(path::AbstractString) = eklef(default_context(), path)

"""
    ekuef([ctx,] handle) -> Nothing

Unload the EK file that `eklef` loaded with the result `handle`, as `ekuef_c`. A handle
that is not loaded changes nothing, as in CSPICE.
"""
ekuef(ctx::SpiceContext, handle::AbstractString) = unload(ctx, handle)
ekuef(handle::AbstractString) = ekuef(default_context(), handle)

"""
    ektnam([ctx,] n) -> String

Name of loaded table number `n` (from 1; `ektnam_c` counts from 0), in the order in
which the tables were loaded. Preconditions: an EK file is loaded
(`SpiceError(:noloadedfiles)`) and `1 <= n <= ekntab()` (`:invalidindex`).
"""
function ektnam(ks::KernelSet, n::Integer)
    isempty(ks.ektables) && throw(SpiceError(:noloadedfiles, "No E-kernels are currently loaded."))
    (1 <= n <= length(ks.ektables)) || throw(SpiceError(:invalidindex,
        "The index $n does not correspond to a loaded table."))
    return ks.ektables[n].name
end
ektnam(ctx::SpiceContext, n::Integer) = ektnam(kernelset(ctx), n)
ektnam(n::Integer) = ektnam(default_context(), n)

"The loaded table `table` of `ks` (name without regard to case); `:tablenotloaded` if absent."
function _ek_table(ks::KernelSet, table::AbstractString)
    name = uppercase(strip(table, ' '))
    k = findfirst(t -> t.name == name, ks.ektables)
    k === nothing && throw(SpiceError(:tablenotloaded, "The table $table is not currently loaded."))
    return ks.ektables[k]
end

"""
    ekccnt([ctx,] table) -> Int

Number of columns of the loaded table `table`, as `ekccnt_c`. Precondition: the table
is loaded, else `SpiceError(:tablenotloaded)`.
"""
ekccnt(ks::KernelSet, table::AbstractString) = length(_ek_table(ks, table).cnames)
ekccnt(ctx::SpiceContext, table::AbstractString) = ekccnt(kernelset(ctx), table)
ekccnt(table::AbstractString) = ekccnt(default_context(), table)

"""
    ekcii([ctx,] table, cindex) -> (column::String, attdsc::EKColumnAttr)

Name and attributes of column number `cindex` (from 1; `ekcii_c` counts from 0) of the
loaded table `table`, as `ekcii_c`. Preconditions: the table is loaded
(`SpiceError(:tablenotloaded)`) and `1 <= cindex <= ekccnt(table)` (`:invalidindex`).
The string length of a column that does not have character type is 0 here; CSPICE
gives the value that the file holds, which is 1.
"""
function ekcii(ks::KernelSet, table::AbstractString, cindex::Integer)
    t = _ek_table(ks, table)
    (1 <= cindex <= length(t.cnames)) || throw(SpiceError(:invalidindex, "Column indices for " *
        "table $(t.name) range from 1 to $(length(t.cnames)); requested index was $cindex."))
    return t.cnames[cindex], t.attrs[cindex]
end
ekcii(ctx::SpiceContext, table::AbstractString, cindex::Integer) = ekcii(kernelset(ctx), table, cindex)
ekcii(table::AbstractString, cindex::Integer) = ekcii(default_context(), table, cindex)

# ---------------------------------------------------------------------------------------
# Queries
# ---------------------------------------------------------------------------------------

"""
    ekpsel([ctx,] query) -> NamedTuple

Analyse the SELECT clause of `query`, as `ekpsel_c`. Returns `(n, xbegs, xends, xtypes,
xclass, tabs, cols, error, errmsg)`: the number of SELECT items, and for each item the
first and last character of its expression in `query` (from 1; `ekpsel_c` counts from
0), its data type (`:CHR`, `:DP`, `:INT`, `:TIME`), its class (`:COL`), and the names of
its table and column. A query that is not valid gives `error == true`, the message in
`errmsg`, and `n == 0`; the call does not throw for it.
"""
ekpsel(ks::KernelSet, td::TimeDefaults, query::AbstractString) =
    _ek_psel(ks.ektables, ks, td, query)
ekpsel(ctx::SpiceContext, query::AbstractString) =
    ekpsel(kernelset(ctx), timedefaults(ctx), query)
ekpsel(query::AbstractString) = ekpsel(default_context(), query)

"""
    ekfind([ctx,] query) -> (result::EKResult, error::Bool, errmsg::String)
    ekfind(ks::KernelSet, td::TimeDefaults, query)

Find the rows of the loaded EK tables that satisfy `query`, as `ekfind_c`. The rows are
in an [`EKResult`](@ref); `result.nmrows` is the output `nmrows` of CSPICE, and
[`eknelt`](@ref), [`ekgc`](@ref), [`ekgd`](@ref), [`ekgi`](@ref) read the result. A
query that is not valid gives `error == true`, the message of CSPICE in `errmsg` and
an empty result; the call does not throw for it. The call reads one snapshot of the
loaded kernels and changes nothing, so many tasks can run queries at the same time,
also while another task loads or unloads kernels. Postcondition: `error` is false or
`result.nmrows == 0`.
"""
function ekfind(ks::KernelSet, td::TimeDefaults, query::AbstractString)
    q, error, errmsg = _ek_encode_query(ks.ektables, ks, td, query)
    error && return _ek_empty_result(), true, errmsg
    return _ek_search(ks, q), false, ""
end
ekfind(ctx::SpiceContext, query::AbstractString) =
    ekfind(kernelset(ctx), timedefaults(ctx), query)
ekfind(query::AbstractString) = ekfind(default_context(), query)

"Segment, column and record pointer of row `row` and SELECT column `selidx` of `res`."
@inline function _ek_result_entry(res::EKResult, selidx::Integer, row::Integer)
    (1 <= row <= res.nmrows) || throw(SpiceError(:invalidindex, "Row indices for query " *
        "result range from 1 to $(res.nmrows); requested row index was $row."))
    (1 <= selidx <= length(res.seltab)) || throw(SpiceError(:invalidindex,
        "The SELECT column index $selidx is out of the valid range 1:$(length(res.seltab))"))
    t = res.seltab[selidx]
    qs = res.segs[res.rowseg[t, row]]
    return qs.file, qs.cols[res.selcol[selidx]], Int(res.rowptr[t, row])
end

"""
    eknelt(result, selidx, row) -> Int

Number of elements of the entry in row `row` and SELECT column `selidx` of `result`
(both from 1; `eknelt_c` counts from 0), as `eknelt_c`. A null entry has 1 element.
Preconditions: `1 <= row <= result.nmrows` and `selidx` names a SELECT column, else
`SpiceError(:invalidindex)`.
"""
function eknelt(res::EKResult, selidx::Integer, row::Integer)
    f, col, recptr = _ek_result_entry(res, selidx, row)
    return _ek_nelt(f, col, recptr)
end

"Checks of the fetch calls: positive element index and a column of one of the `types`."
@inline function _ek_fetch_check(col::_EKColumn, elment::Integer, types)
    elment >= 1 || throw(SpiceError(:invalidindex, "ELMENT must be positive but was $elment."))
    col.dtype in types || throw(SpiceError(:invalidtype,
        "Column $(col.name) has data type $(_EKQ_TYPSTR[col.dtype])."))
    return nothing
end

"""
    ekgc(result, selidx, row, elment=1) -> (cdata::String, null::Bool, found::Bool)

Element `elment` of the character entry in row `row` and SELECT column `selidx` of
`result` (all from 1; `ekgc_c` counts from 0), without trailing blanks, as `ekgc_c`.
`null` tells that the entry is null (then `cdata` is empty), `found` that the entry
has an element `elment`. Preconditions (each throws `SpiceError`): valid `row` and
`selidx` and `elment >= 1` (`:invalidindex`); the column has character type
(`:invalidtype`).
"""
function ekgc(res::EKResult, selidx::Integer, row::Integer, elment::Integer=1)
    f, col, recptr = _ek_result_entry(res, selidx, row)
    _ek_fetch_check(col, elment, (_EK_CHR,))
    v, isnull, found = _ek_elt_c(f, col, recptr, elment)
    return String(rstrip(v, ' ')), isnull, found
end

"""
    ekgd(result, selidx, row, elment=1) -> (ddata::Float64, null::Bool, found::Bool)

Element `elment` of the double precision or time entry in row `row` and SELECT column
`selidx` of `result`, as `ekgd_c`; a time is in seconds past J2000 TDB. Results and
preconditions as [`ekgc`](@ref).
"""
function ekgd(res::EKResult, selidx::Integer, row::Integer, elment::Integer=1)
    f, col, recptr = _ek_result_entry(res, selidx, row)
    _ek_fetch_check(col, elment, (_EK_DP, _EK_TIME))
    return _ek_elt_d(f, col, recptr, elment)
end

"""
    ekgi(result, selidx, row, elment=1) -> (idata::Int, null::Bool, found::Bool)

Element `elment` of the integer entry in row `row` and SELECT column `selidx` of
`result`, as `ekgi_c`. Results and preconditions as [`ekgc`](@ref).
"""
function ekgi(res::EKResult, selidx::Integer, row::Integer, elment::Integer=1)
    f, col, recptr = _ek_result_entry(res, selidx, row)
    _ek_fetch_check(col, elment, (_EK_INT,))
    return _ek_elt_i(f, col, recptr, elment)
end

"""
    ekcolumn(result, selidx) -> Vector

All entries of the SELECT column `selidx` of `result`, one item for each row: `missing`
for a null entry, the value for a scalar column, and the vector of the values for an
array column. Strings are without trailing blanks. Not a CSPICE call: it replaces the
loop over `eknelt` and the fetch call. Precondition: `selidx` names a SELECT column,
else `SpiceError(:invalidindex)`.
"""
function ekcolumn(res::EKResult, selidx::Integer)
    (1 <= selidx <= length(res.seltab)) || throw(SpiceError(:invalidindex,
        "The SELECT column index $selidx is out of the valid range 1:$(length(res.seltab))"))
    out = Vector{Any}(undef, res.nmrows)
    for row in 1:res.nmrows              # invariant: out[1:row-1] is set
        f, col, recptr = _ek_result_entry(res, selidx, row)
        vals, isnull = col.dtype == _EK_CHR ? _ek_entry_c(f, col, recptr) :
                       col.dtype == _EK_INT ? _ek_entry_i(f, col, recptr) :
                       _ek_entry_d(f, col, recptr)
        col.dtype == _EK_CHR && (vals = String[String(rstrip(v, ' ')) for v in vals])
        out[row] = isnull ? missing : (col.cclass <= 3 || col.cclass >= 7) ? vals[1] : vals
    end
    return identity.(out)
end
