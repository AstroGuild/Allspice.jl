# Data types of the E-kernel (EK) subsystem that the shared types of `types.jl` use.
#
# An EK file is a DAS file with the ID word `DAS/EK`. It holds segments; a segment holds
# the rows of one table (several segments, in one file or in several files, can hold
# rows of the same table). A row has one entry for each column; an entry is one value or
# an array of values of the column data type, or null.

"Codes of the EK data types, as `ektype.inc`: character, double, integer, time."
const _EK_CHR = 1
const _EK_DP = 2
const _EK_INT = 3
const _EK_TIME = 4

"Names of the EK data types in code order; `EKColumnAttr.dtype` holds one of them."
const _EK_TYPE_SYMBOLS = (:CHR, :DP, :INT, :TIME)

@inline _ek_dtype_code(s::Symbol) =
    s === :CHR ? 1 : s === :DP ? 2 : s === :INT ? 3 : s === :TIME ? 4 : 0

"""
    EKColumnAttr

Attributes of one EK column, as the attribute descriptor of `ekcii_c` and `ekssum_c`.
[`ekcii`](@ref) and [`ekssum`](@ref) return it.

- `cclass`: the column class of CSPICE, 1 to 9. It is the code of the storage layout
  in the file; a caller needs it only to compare with CSPICE. Classes 7 to 9 are the
  columns with `FIXED_COUNT = TRUE` (see [`ekbseg`](@ref));
- `dtype`: `:CHR`, `:DP`, `:INT` or `:TIME`;
- `strlen`: the string length of a character column, -1 for strings of variable
  length, and 0 for a column of another type (CSPICE gives 1 there);
- `size`: the number of values of an entry, -1 for entries of variable size;
- `indexd`: true if the column has an index;
- `nullok`: true if the column can hold null entries.
"""
struct EKColumnAttr
    cclass::Int
    dtype::Symbol
    strlen::Int
    size::Int
    indexd::Bool
    nullok::Bool
end

"""
    EKSegmentSummary

Summary of one EK segment, the result of [`ekssum`](@ref).

- `tabnam`: the name of the table of the segment;
- `nrows`, `ncols`: the numbers of rows and columns;
- `cnames`: the column names, in the order of the segment;
- `cdescrs`: the [`EKColumnAttr`](@ref) of each column, in the same order.

`cnames` and `cdescrs` have `ncols` elements.
"""
struct EKSegmentSummary
    tabnam::String
    nrows::Int
    ncols::Int
    cnames::Vector{String}
    cdescrs::Vector{EKColumnAttr}
end

"""
    _EKColumn

One decoded column descriptor of a segment (the 11 integers of `ekcoldsc.inc`) with the
column name: class, data type code, string length (-1 variable), entry size (-1
variable), index type (-1 no index, 1 the index tree of a class 1 to 6 column, 2 the
order vector of a class 7 to 9 column), index pointer, null flag, ordinal position in
the segment, and the pointer to the integer metadata of a class 7 to 9 column.
"""
struct _EKColumn
    name::String
    cclass::Int32
    dtype::Int32
    strlen::Int32
    size::Int32
    ixtype::Int32
    ixptr::Int32
    nullok::Bool
    ord::Int32
    metptr::Int32
end

"""
    _EKSegment

One segment of a loaded EK file. `stype` is the segment type (1: rows are records
found through the record tree; 2: fixed layout of the fast writers). `recptrs[r]` is
the record pointer of row `r` (for type 2 the row number itself). `segdsc` holds the 24
integers of the segment descriptor. Never modified after it was built.
"""
struct _EKSegment
    segno::Int32
    stype::Int32
    table::String
    nrows::Int32
    columns::Vector{_EKColumn}
    recptrs::Vector{Int32}
    segdsc::Vector{Int32}
end

"""
    EKFile

One loaded EK file: the three DAS address spaces (`chars[a]`, `dps[a]`, `ints[a]` are
the character, the double and the integer at DAS address `a`; memory-mapped where the
file allows it) and the decoded segments in file order. `kind` is always `:EK`.
Nothing in an `EKFile` is modified after `furnsh` has built it.
"""
struct EKFile
    path::String
    kind::Symbol
    chars::Vector{UInt8}
    dps::Vector{Float64}
    ints::Vector{Int32}
    segments::Vector{_EKSegment}
end

"""
    EKTable

One loaded EK table: all loaded segments that have the same table name. This is the
type of the entries of the table catalogue of a context. No call returns it and a
caller does not build it: read the catalogue with [`ekntab`](@ref), [`ektnam`](@ref),
[`ekccnt`](@ref) and [`ekcii`](@ref).

All segments of a table have columns with the same names and attributes; a file whose
segment does not fit is rejected at load time, as in CSPICE (see [`eklef`](@ref)).
"""
struct EKTable
    name::String
    cnames::Vector{String}
    attrs::Vector{EKColumnAttr}
    segs::Vector{NTuple{2,Int32}}
    cols::Vector{Vector{_EKColumn}}
end

# A catalogue entry without column descriptors.
EKTable(name, cnames, attrs, segs) = EKTable(name, cnames, attrs, segs, Vector{_EKColumn}[])

_ek_attr(c::_EKColumn) = EKColumnAttr(Int(c.cclass), _EK_TYPE_SYMBOLS[c.dtype],
    c.dtype == _EK_CHR ? Int(c.strlen) : 0, Int(c.size), c.ixtype != -1, c.nullok)

_ek_same_attributes(a::_EKColumn, b::_EKColumn) = a.cclass == b.cclass && a.dtype == b.dtype &&
    a.strlen == b.strlen && a.size == b.size && (a.ixtype != -1) == (b.ixtype != -1) &&
    a.nullok == b.nullok

"""
    _ek_build_tables(files) -> Vector{EKTable}

The table catalogue of the loaded EK `files` (load order), as `eklef` builds it: tables
in the order of their first segment, the segments of a table in load order, the columns
of a table in the order of its first segment. Throws `SpiceError` as `eklef` if a
segment does not fit a table that an earlier segment defined: `:ekcolnummismatch`
(number of columns), `:ekmissingcolumn`, `:badattributes` (class, type, string length,
size, index or null attribute of a column differs).
"""
function _ek_build_tables(files::Vector{EKFile})
    tables = EKTable[]
    for (fno, file) in enumerate(files), seg in file.segments
        k = findfirst(t -> t.name == seg.table, tables)
        if k === nothing
            push!(tables, EKTable(seg.table, String[c.name for c in seg.columns],
                EKColumnAttr[_ek_attr(c) for c in seg.columns],
                NTuple{2,Int32}[(Int32(fno), seg.segno)], Vector{_EKColumn}[seg.columns]))
            continue
        end
        t = tables[k]
        length(seg.columns) == length(t.cnames) || throw(SpiceError(:ekcolnummismatch,
            "The EK file $(file.path) could not be loaded; the number of columns " *
            "($(length(seg.columns))) in segment $(seg.segno) does not match the number of " *
            "columns ($(length(t.cnames))) in the parent table $(t.name)."))
        cols = Vector{_EKColumn}(undef, length(t.cnames))
        for (i, name) in enumerate(t.cnames)   # invariant: cols[1:i-1] match the table
            j = findfirst(c -> c.name == name, seg.columns)
            j === nothing && throw(SpiceError(:ekmissingcolumn, "The EK file $(file.path) " *
                "could not be loaded; the column $name in already loaded table $(t.name) " *
                "is not present in segment $(seg.segno) in the EK file."))
            _ek_same_attributes(seg.columns[j], t.cols[1][i]) || throw(SpiceError(:badattributes,
                "EK file $(file.path) contains a column whose attributes conflict with a " *
                "loaded column.  The offending column name is $name; the column is in " *
                "segment $(seg.segno) of the file."))
            cols[i] = seg.columns[j]
        end
        push!(t.segs, (Int32(fno), seg.segno))
        push!(t.cols, cols)
    end
    return tables
end
