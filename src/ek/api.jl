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

Open the EK file `path` for read access and return its handle, as `ekopr_c`.

The caller owns the handle; close it with [`ekcls`](@ref). The handle serves the calls
that read one file ([`eknseg`](@ref), [`ekssum`](@ref), [`ekrcec`](@ref),
[`ekrced`](@ref), [`ekrcei`](@ref)). To run queries, load the file with
[`eklef`](@ref) instead.

Throws `SpiceError(:filenotfound)` if the file does not exist, the other kinds of
`dasopr` if it is not a DAS file, and `SpiceError(:wrongarchitecture)` or
`SpiceError(:invalidformat)` if the DAS file is not a paged EK (architecture 8).

See also [`ekopw`](@ref), [`ekopn`](@ref).
"""
function ekopr(path::AbstractString)
    h = _das_load(path, "ekopr", false)
    _ek_check(h, h.path, "ekopr")
    return h
end

"""
    ekcls(handle) -> Nothing

Close the EK file of `handle`, as `ekcls_c`.

For a write handle all data are on disk when the call returns. Closing a closed handle
does nothing.

See also [`ekopr`](@ref), [`ekopw`](@ref), [`ekopn`](@ref).
"""
ekcls(h::DASHandle) = dascls(h)

@inline function _ek_check_read(h::DASHandle, who::AbstractString)
    _das_check_open(h, who)
    _ek_check(h, h.path, who)
    return nothing
end

"""
    eknseg(handle) -> Int

Return the number of segments of the EK file of `handle`, as `eknseg_c`.

Throws `SpiceError(:dasnosuchhandle)` if the handle is closed, and
`SpiceError(:wrongarchitecture)` or `SpiceError(:invalidformat)` if the file is not a
paged EK.
"""
function eknseg(h::DASHandle)
    _ek_check_read(h, "eknseg")
    return _ek_nseg(h)
end

"""
    ekssum(handle, segno) -> EKSegmentSummary

Return the summary of segment `segno` of the EK file of `handle`, as `ekssum_c`: the
table name, the number of rows, and the name and attributes of every column.

`segno` counts from 1; `ekssum_c` counts from 0.

Throws `SpiceError(:indexoutofrange)` unless `1 <= segno <= eknseg(handle)`, and the
errors of [`eknseg`](@ref) for the handle.

See also [`EKSegmentSummary`](@ref), [`EKColumnAttr`](@ref).
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

Return all strings of the entry of the character column `column` in row `recno` of
segment `segno`, without trailing blanks, as `ekrcec_c`.

`segno` and `recno` count from 1; `ekrcec_c` counts from 0. The column name is matched
without regard to case. The entry of a scalar column gives a vector of one string. A
null entry gives an empty vector and `isnull == true`; the C call returns the element
count and the flag in separate outputs.

Throws `SpiceError(:invalidindex)` if `segno` is not a segment of the file,
`SpiceError(:bug)` if the segment has no column `column` (as CSPICE),
`SpiceError(:wrongdatatype)` if the column does not have character type,
`SpiceError(:indexoutofrange)` unless `1 <= recno <= nrows`, and
`SpiceError(:uninitialized)` if the entry was never written. A closed handle gives the
errors of [`eknseg`](@ref).

See also [`ekrced`](@ref), [`ekrcei`](@ref), [`ekfind`](@ref).
"""
function ekrcec(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString)
    col, recptr = _ek_locate(h, "ekrcec", segno, recno, column, (_EK_CHR,))
    vals, isnull = _ek_entry_c(h, col, recptr)
    return String[String(rstrip(v, ' ')) for v in vals], isnull
end

"""
    ekrced(handle, segno, recno, column) -> (dvals::Vector{Float64}, isnull::Bool)

Return all values of the entry of the double precision or time column `column` in row
`recno` of segment `segno`, as `ekrced_c`. A time is in seconds past J2000 TDB.

Indices and results are as for [`ekrcec`](@ref): `segno` and `recno` count from 1, and
a null entry gives an empty vector and `isnull == true`.

Throws the errors of [`ekrcec`](@ref), with `SpiceError(:wrongdatatype)` if the column
has character or integer type and `SpiceError(:uninitializedvalue)` if the entry was
never written.
"""
function ekrced(h::DASHandle, segno::Integer, recno::Integer, column::AbstractString)
    col, recptr = _ek_locate(h, "ekrced", segno, recno, column, (_EK_DP, _EK_TIME))
    return _ek_entry_d(h, col, recptr)
end

"""
    ekrcei(handle, segno, recno, column) -> (ivals::Vector{Int}, isnull::Bool)

Return all values of the entry of the integer column `column` in row `recno` of
segment `segno`, as `ekrcei_c`.

Indices and results are as for [`ekrcec`](@ref): `segno` and `recno` count from 1, and
a null entry gives an empty vector and `isnull == true`.

Throws the errors of [`ekrcec`](@ref), with `SpiceError(:wrongdatatype)` if the column
does not have integer type and `SpiceError(:uninitializedvalue)` if the entry was never
written.
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

Return the number of loaded EK tables, as `ekntab_c`; 0 if no EK file is loaded.

Segments of the same table name, in one file or in several, count as one table.

See also [`ektnam`](@ref), [`eklef`](@ref).
"""
ekntab(ks::KernelSet) = length(ks.ektables)
ekntab(ctx::SpiceContext) = ekntab(kernelset(ctx))
ekntab() = ekntab(default_context())

"""
    eklef([ctx,] path) -> String

Load the EK file at `path` for queries, as `eklef_c`, and return the handle for
[`ekuef`](@ref).

The handle is the path itself, not an integer as in CSPICE, because a context has no
file handles. The call has the effect of [`furnsh`](@ref) on a file that is an EK.

Throws `SpiceError(:nofile)` if the file does not exist and
`SpiceError(:wrongarchitecture)` if it is not an EK; the context is then unchanged.
Throws `SpiceError(:ekcolnummismatch)`, `SpiceError(:ekmissingcolumn)` or
`SpiceError(:badattributes)` if a segment of the file does not have the columns of a
loaded table of the same name.

See also [`ekfind`](@ref), [`ekntab`](@ref).
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

Unload the EK file that [`eklef`](@ref) loaded with the result `handle`, as `ekuef_c`.

A handle that is not loaded changes nothing, as in CSPICE. An [`EKResult`](@ref) of an
earlier [`ekfind`](@ref) stays readable after the call.
"""
ekuef(ctx::SpiceContext, handle::AbstractString) = unload(ctx, handle)
ekuef(handle::AbstractString) = ekuef(default_context(), handle)

"""
    ektnam([ctx,] n) -> String

Return the name, in upper case, of loaded table number `n`, as `ektnam_c`.

`n` counts from 1; `ektnam_c` counts from 0. The tables are in the order in which they
were loaded.

Throws `SpiceError(:noloadedfiles)` if no EK file is loaded and
`SpiceError(:invalidindex)` unless `1 <= n <= ekntab()`.
"""
function ektnam(ks::KernelSet, n::Integer)
    isempty(ks.ektables) && throw(SpiceError(:noloadedfiles, "No E-kernels are currently loaded."))
    (1 <= n <= length(ks.ektables)) || throw(SpiceError(:invalidindex,
        "The index $n does not correspond to a loaded table."))
    return ks.ektables[n].name
end
ektnam(ctx::SpiceContext, n::Integer) = ektnam(kernelset(ctx), n)
ektnam(n::Integer) = ektnam(default_context(), n)

function _ek_table(ks::KernelSet, table::AbstractString)
    name = uppercase(strip(table, ' '))
    k = findfirst(t -> t.name == name, ks.ektables)
    k === nothing && throw(SpiceError(:tablenotloaded, "The table $table is not currently loaded."))
    return ks.ektables[k]
end

"""
    ekccnt([ctx,] table) -> Int

Return the number of columns of the loaded table `table`, as `ekccnt_c`.

The table name is matched without regard to case.

Throws `SpiceError(:tablenotloaded)` if the table is not loaded.

See also [`ekcii`](@ref).
"""
ekccnt(ks::KernelSet, table::AbstractString) = length(_ek_table(ks, table).cnames)
ekccnt(ctx::SpiceContext, table::AbstractString) = ekccnt(kernelset(ctx), table)
ekccnt(table::AbstractString) = ekccnt(default_context(), table)

"""
    ekcii([ctx,] table, cindex) -> (column::String, attdsc::EKColumnAttr)

Return the name and the attributes of column number `cindex` of the loaded table
`table`, as `ekcii_c`.

`cindex` counts from 1; `ekcii_c` counts from 0. The columns are in the order of the
first loaded segment of the table. The string length `attdsc.strlen` of a column that
does not have character type is 0; CSPICE gives the value that the file holds, which
is 1.

Throws `SpiceError(:tablenotloaded)` if the table is not loaded and
`SpiceError(:invalidindex)` unless `1 <= cindex <= ekccnt(table)`.

See also [`EKColumnAttr`](@ref), [`ekccnt`](@ref).
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

Return the description of the SELECT clause of `query`, as `ekpsel_c`, in the named
tuple `(n, xbegs, xends, xtypes, xclass, tabs, cols, error, errmsg)`.

- `n`: the number of SELECT items. The vectors have one element for each item.
- `xbegs`, `xends`: first and last byte of the item in `query`, from 1 (`ekpsel_c`
  counts from 0). A qualified name runs from the qualifier to the column name.
- `xtypes`: the data type, `:CHR`, `:DP`, `:INT` or `:TIME`.
- `xclass`: the class of the item, always `:COL` (a column).
- `tabs`, `cols`: the table name (never the alias) and the column name, in upper case.
- `error`, `errmsg`: the error flag and the message of CSPICE for the query.

A query that is not valid gives `error == true`, `n == 0` and empty vectors; the call
does not throw for it. The query is checked against the loaded tables, as in
[`ekfind`](@ref).
"""
ekpsel(ks::KernelSet, td::TimeDefaults, query::AbstractString) =
    _ek_psel(ks.ektables, ks, td, query)
ekpsel(ctx::SpiceContext, query::AbstractString) =
    ekpsel(kernelset(ctx), timedefaults(ctx), query)
ekpsel(query::AbstractString) = ekpsel(default_context(), query)

"""
    ekfind([ctx,] query) -> (result::EKResult, error::Bool, errmsg::String)

Find the rows of the loaded EK tables that satisfy `query`, as `ekfind_c`, and return
them in an [`EKResult`](@ref) with the error flag and the error message of the query.

The query has the form `SELECT <column>, ... FROM <table> [<alias>], ... [WHERE
<constraints>] [ORDER BY <column> [ASC|DESC], ...]`, with the language of CSPICE.

CSPICE keeps the rows of the last query in global state. Here the rows are in
`result`, which the fetch calls take as first argument: `result.nmrows` (or
`length(result)`) is the output `nmrows` of `ekfind_c`, [`ekgc`](@ref),
[`ekgd`](@ref) and [`ekgi`](@ref) read one element, [`eknelt`](@ref) gives the size
of an entry, and [`ekcolumn`](@ref) reads a whole SELECT column. Rows, SELECT columns
and elements count from 1; CSPICE counts from 0.

A query that is not valid gives `error == true`, the message of CSPICE in `errmsg`
and a result without rows; the call does not throw for it. Positions in a message are
byte positions in `query`, from 1. A time value in a query is converted with the
leap-second and clock kernels of the context, and a failure of that conversion is
also reported in `errmsg`.

The call reads one snapshot of the loaded kernels and changes nothing, so many tasks
can run queries at the same time, also while another task loads or unloads kernels.

The rows differ from those of CSPICE N0067 only where CSPICE has a defect: rows that
several OR branches of the WHERE clause match (CSPICE can return such a row twice, or
drop other rows), a double precision value compared with an indexed integer column
(`ID < 17.4` loses the rows with `ID = 17` in CSPICE), and a number such as `2.756`,
which the number parser of CSPICE reads one unit in the last place off.

Throws `SpiceError(:celltoosmall)` if the WHERE clause expands to more constraints
than an encoded query holds (1000 without other clauses),
`SpiceError(:toomanyconjunctions)` if it expands to more than 1000 OR branches, and
`SpiceError(:buffertoosmall)` if the converted times do not fit in the buffer of 100
numbers.

# Examples
```jldoctest
julia> path = joinpath(mktempdir(), "moons.bes");

julia> h = ekopn(path, "moons");

julia> seg = ekbseg(h, "MOONS", ["NAME", "PLANET", "RADIUS"],
           ["DATATYPE = CHARACTER*(16)", "DATATYPE = INTEGER",
            "DATATYPE = DOUBLE PRECISION, NULLS_OK = TRUE"]);

julia> for (name, planet, radius) in (("PHOBOS", 4, 11.1), ("DEIMOS", 4, nothing), ("IO", 5, 1821.6))
           row = ekappr(h, seg)
           ekacec(h, seg, row, "NAME", [name])
           ekacei(h, seg, row, "PLANET", [planet])
           if radius === nothing
               ekaced(h, seg, row, "RADIUS", Float64[], true)   # null entry
           else
               ekaced(h, seg, row, "RADIUS", [radius])
           end
       end

julia> ekcls(h)

julia> ctx = SpiceContext();

julia> eklef(ctx, path);

julia> result, failed, errmsg = ekfind(ctx, "SELECT NAME, RADIUS FROM MOONS WHERE PLANET = 4 ORDER BY NAME");

julia> failed
false

julia> result
EKResult(2 rows, columns MOONS.NAME, MOONS.RADIUS)

julia> ekgc(result, 1, 2)      # SELECT column 1, row 2: (value, null, found)
("PHOBOS", false, true)

julia> ekgd(result, 2, 1)      # the radius of DEIMOS is null
(0.0, true, true)

julia> ekcolumn(result, 2)
2-element Vector{Union{Missing, Float64}}:
   missing
 11.1

julia> ekfind(ctx, "SELECT NAME FROM MOONS WHERE")[2:3]
(true, "Empty WHERE clause.")
```

See also [`eklef`](@ref), [`ekpsel`](@ref), [`ekrcec`](@ref).
"""
function ekfind(ks::KernelSet, td::TimeDefaults, query::AbstractString)
    q, error, errmsg = _ek_encode_query(ks.ektables, ks, td, query)
    error && return _ek_empty_result(), true, errmsg
    return _ek_search(ks, q), false, ""
end
ekfind(ctx::SpiceContext, query::AbstractString) =
    ekfind(kernelset(ctx), timedefaults(ctx), query)
ekfind(query::AbstractString) = ekfind(default_context(), query)

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

Return the number of elements of the entry in row `row` and SELECT column `selidx` of
`result`, as `eknelt_c`.

`selidx` and `row` count from 1; `eknelt_c` counts from 0. The result is 1 for a
scalar column and the entry size for an array column of fixed size, also when the
entry is null. A null entry of an array column of variable size has 1 element.

Throws `SpiceError(:invalidindex)` unless `1 <= row <= result.nmrows` and `selidx` is
the position of a SELECT column.

See also [`ekfind`](@ref), [`ekcolumn`](@ref).
"""
function eknelt(res::EKResult, selidx::Integer, row::Integer)
    f, col, recptr = _ek_result_entry(res, selidx, row)
    return _ek_nelt(f, col, recptr)
end

@inline function _ek_fetch_check(col::_EKColumn, elment::Integer, types)
    elment >= 1 || throw(SpiceError(:invalidindex, "ELMENT must be positive but was $elment."))
    col.dtype in types || throw(SpiceError(:invalidtype,
        "Column $(col.name) has data type $(_EKQ_TYPSTR[col.dtype])."))
    return nothing
end

"""
    ekgc(result, selidx, row, elment=1) -> (cdata::String, null::Bool, found::Bool)

Return element `elment` of the character entry in row `row` and SELECT column `selidx`
of `result`, without trailing blanks, as `ekgc_c`.

`selidx`, `row` and `elment` count from 1; `ekgc_c` counts from 0. `null` is true if
the entry is null; `cdata` is then empty and `found` is true. `found` is false if an
array entry has fewer than `elment` elements. For a scalar column `elment` is not
used, and `found` is true.

Throws `SpiceError(:invalidindex)` unless `1 <= row <= result.nmrows`, `selidx` is the
position of a SELECT column and `elment >= 1`, and `SpiceError(:invalidtype)` if the
column does not have character type.

See [`ekfind`](@ref) for an example. See also [`ekgd`](@ref), [`ekgi`](@ref),
[`eknelt`](@ref), [`ekcolumn`](@ref).
"""
function ekgc(res::EKResult, selidx::Integer, row::Integer, elment::Integer=1)
    f, col, recptr = _ek_result_entry(res, selidx, row)
    _ek_fetch_check(col, elment, (_EK_CHR,))
    v, isnull, found = _ek_elt_c(f, col, recptr, elment)
    return String(rstrip(v, ' ')), isnull, found
end

"""
    ekgd(result, selidx, row, elment=1) -> (ddata::Float64, null::Bool, found::Bool)

Return element `elment` of the double precision or time entry in row `row` and SELECT
column `selidx` of `result`, as `ekgd_c`. A time is in seconds past J2000 TDB.

Indices, `null` and `found` are as for [`ekgc`](@ref); `ddata` is 0.0 for a null entry
and for an element that is not found.

Throws `SpiceError(:invalidindex)` for a bad index, as [`ekgc`](@ref), and
`SpiceError(:invalidtype)` if the column has character or integer type.
"""
function ekgd(res::EKResult, selidx::Integer, row::Integer, elment::Integer=1)
    f, col, recptr = _ek_result_entry(res, selidx, row)
    _ek_fetch_check(col, elment, (_EK_DP, _EK_TIME))
    return _ek_elt_d(f, col, recptr, elment)
end

"""
    ekgi(result, selidx, row, elment=1) -> (idata::Int, null::Bool, found::Bool)

Return element `elment` of the integer entry in row `row` and SELECT column `selidx` of
`result`, as `ekgi_c`.

Indices, `null` and `found` are as for [`ekgc`](@ref); `idata` is 0 for a null entry
and for an element that is not found.

Throws `SpiceError(:invalidindex)` for a bad index, as [`ekgc`](@ref), and
`SpiceError(:invalidtype)` if the column does not have integer type.
"""
function ekgi(res::EKResult, selidx::Integer, row::Integer, elment::Integer=1)
    f, col, recptr = _ek_result_entry(res, selidx, row)
    _ek_fetch_check(col, elment, (_EK_INT,))
    return _ek_elt_i(f, col, recptr, elment)
end

"""
    ekcolumn(result, selidx) -> Vector

Return all entries of the SELECT column `selidx` (from 1) of `result`, one item for
each row, in the order of the result.

An item is `missing` for a null entry, the value for a scalar column, and the vector
of the values for an array column. Strings are without trailing blanks; a time is in
seconds past J2000 TDB. The element type is the narrowest that holds the items, for
example `Vector{Union{Missing, Float64}}`. CSPICE has no such call: it replaces the
loop over [`eknelt`](@ref) and [`ekgc`](@ref), [`ekgd`](@ref) or [`ekgi`](@ref).

Throws `SpiceError(:invalidindex)` unless `selidx` is the position of a SELECT column,
and `SpiceError(:uninitialized)` or `SpiceError(:uninitializedvalue)` if an entry was
never written.

See [`ekfind`](@ref) for an example.
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
