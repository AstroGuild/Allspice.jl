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

"Code 1 to 4 of the data type symbol `s` (`:CHR`, `:DP`, `:INT`, `:TIME`), or 0."
@inline _ek_dtype_code(s::Symbol) =
    s === :CHR ? 1 : s === :DP ? 2 : s === :INT ? 3 : s === :TIME ? 4 : 0

"""
    EKColumnAttr

Attributes of one EK column, as the attribute descriptor of `ekcii_c` and `ekssum_c`:

- `cclass`: the column class, 1 to 9 (the storage layout; classes 7 to 9 are the
  columns of a segment that the fast writers made);
- `dtype`: `:CHR`, `:DP`, `:INT` or `:TIME`;
- `strlen`: the string length of a character column, -1 for strings of variable
  length, and 0 for a column of another type;
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

Summary of one EK segment, as `ekssum_c` returns it: the table name `tabnam`, the
numbers of rows and columns, and for each column (in the order of the segment) its
name in `cnames` and its attributes in `cdescrs`. Invariant: `length(cnames) ==
length(cdescrs) == ncols`.
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

One loaded table of a `KernelSet`: its name (upper case), the names (upper case) and
attributes of its columns, and its segments as (file number, segment number) pairs in
load order, where the file number indexes `KernelSet.ekfiles`. In `attrs` the fields
`indexd` and `nullok` are those of the first loaded segment, as in CSPICE; `cclass` is
that of the first segment too. Never modified after it was built.
"""
struct EKTable
    name::String
    cnames::Vector{String}
    attrs::Vector{EKColumnAttr}
    segs::Vector{NTuple{2,Int32}}
end
