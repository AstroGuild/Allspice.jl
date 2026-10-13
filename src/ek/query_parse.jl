# The query parser of the E-kernel (EK) subsystem: the part of `ekfind_c` that runs
# before the search, and `ekpsel_c`.
#
# A query has the form
#
#     SELECT <column>, ... FROM <table> [<alias>], ... [WHERE <constraints>]
#                                                       [ORDER BY <column> [ASC|DESC], ...]
#
# and CSPICE turns it into an "encoded query" in six steps, which this file ports in
# the same order and with the same rules:
#
# 1. `zzekscan`: the scanner splits the string into tokens (keywords, identifiers,
#    integers, doubles, quoted strings and special characters);
# 2. `zzekpars`: the parser reads the FROM clause, then the WHERE clause, then the
#    ORDER BY clause and at last the SELECT clause, and it finds each clause by the
#    first occurrence of its keyword;
# 3. `zzeknrml`: the WHERE clause becomes a disjunction of conjunctions of relations
#    ("normalisation"): BETWEEN and NOT BETWEEN expand to two conjunctions, NOT is
#    pushed down to the relations by the rules of De Morgan, and AND is distributed
#    over OR;
# 4. `zzeknres`: name resolution gives each table and each column its position;
# 5. `zzektres`: time resolution converts the strings that are compared with a column
#    of type TIME into ephemeris time (TDB seconds past J2000);
# 6. `zzeksemc`: the semantic check tests the data types and the column sizes.
#
# The encoded query of CSPICE is an integer cell, a character buffer and a double
# buffer. This port keeps that layout internally (`_EKQEnc`), because the layout fixes
# the limits and the order of all items, and it then decodes the result into the plain
# structure `_EKQuery`. Positions in a query are 1-based byte positions. Nothing here
# has global state: every function is a pure function of its arguments.
#
# Known differences from CSPICE are listed at `_ek_encode_query`.

# ---------------------------------------------------------------------------------
# Limits and codes
# ---------------------------------------------------------------------------------

"Largest number of characters of a query (`MAXQRY`)."
const _EKQ_MAXQRY = 2000
"Largest number of tokens of a query (`MAXTOK`)."
const _EKQ_MAXTOK = 500
"Largest number of numeric tokens of a query (`MAXNUM`), and the size of the double buffer."
const _EKQ_MAXNUM = 100
"Largest number of tables in a FROM clause (`MAXTAB`)."
const _EKQ_MAXTAB = 10
"Largest number of order-by columns (`MAXORD`)."
const _EKQ_MAXORD = 10
"Largest number of SELECT columns (`MAXSEL`)."
const _EKQ_MAXSEL = 50
"Largest number of constraints of a normalised WHERE clause (`MAXCON`)."
const _EKQ_MAXCON = 1000
"Size of the integer cell of an encoded query (`EQIMIN`)."
const _EKQ_ISIZE = 27869
"Size of the character buffer of an encoded query (`MAXCLN`)."
const _EKQ_CSIZE = 2000
"Length of a table name or an alias in name resolution (`TNAMSZ`)."
const _EKQ_TNAMSZ = 64
"Length of the time string that time resolution converts (`TIMLEN` of `zzektres`)."
const _EKQ_TIMLEN = 32
"Size of the relation, conjunction and disjunction tables of the normaliser."
const _EKQ_NPOOL = 5000
"Largest number of conjunctions that the normaliser can count (`MAXCNJ`)."
const _EKQ_MAXCNJ = 1000

# Token codes of `ektokn.inc`.
const _EKQ_TKKEY = 1        # keyword, or a comparison operator written with special characters
const _EKQ_TKID = 2         # identifier
const _EKQ_TKINT = 3        # integer
const _EKQ_TKDP = 4         # double
const _EKQ_TKQSTR = 5       # quoted string
const _EKQ_TKLPAR = 6
const _EKQ_TKRPAR = 7
const _EKQ_TKCOMA = 8
const _EKQ_TKDOT = 9
const _EKQ_TKSTAR = 10
const _EKQ_TKEOQ = 11       # end of query (the parser makes it)

"Keywords of the query language in code order (`ekkeyw.inc`): the code is the position."
const _EKQ_KEYWORDS = ("ALL", "AND", "ASC", "AVG", "BETWEEN", "BY", "COUNT", "DESC",
    "DISTINCT", "EQ", "FROM", "GE", "GROUP", "GT", "HAVING", "IS", "LE", "LIKE", "LT",
    "MAX", "MIN", "NE", "NOT", "NULL", "OR", "ORDER", "SELECT", "SUM", "WHERE")
const _EKQ_KWAND = 2
const _EKQ_KWASC = 3
const _EKQ_KWBETW = 5
const _EKQ_KWBY = 6
const _EKQ_KWDESC = 8
const _EKQ_KWEQ = 10
const _EKQ_KWFROM = 11
const _EKQ_KWGE = 12
const _EKQ_KWGT = 14
const _EKQ_KWIS = 16
const _EKQ_KWLE = 17
const _EKQ_KWLIKE = 18
const _EKQ_KWLT = 19
const _EKQ_KWNE = 22
const _EKQ_KWNOT = 23
const _EKQ_KWNULL = 24
const _EKQ_KWOR = 25
const _EKQ_KWORDR = 26
const _EKQ_KWSEL = 27
const _EKQ_KWWHER = 29

"Special character strings of the scanner with their token codes and keyword codes."
const _EKQ_SPECIAL = (("!=", 1, 22), ("^=", 1, 22), ("<>", 1, 22), ("<=", 1, 17),
    (">=", 1, 12), ("<", 1, 19), (">", 1, 14), ("=", 1, 10), ("(", 6, 0), (")", 7, 0),
    (",", 8, 0), (".", 9, 0), ("*", 10, 0))

# Operator codes of `ekopcd.inc`.
const _EKQ_EQ = 1
const _EKQ_GE = 2
const _EKQ_GT = 3
const _EKQ_LE = 4
const _EKQ_LT = 5
const _EKQ_NE = 6
const _EKQ_LIKE = 7
const _EKQ_UNLIKE = 8
const _EKQ_ISNULL = 9
const _EKQ_NOTNUL = 10
"Complement of the operator codes 1 to 8: EQ/NE, GE/LT, GT/LE, LIKE/UNLIKE."
const _EKQ_CMPNEG = (6, 5, 4, 3, 2, 1, 8, 7)

"Constraint type: column operator column (`EQCOL`)."
const _EKQ_EQCOL = 1
"Constraint type: column operator value (`EQVAL`)."
const _EKQ_EQVAL = 2
"Order sense: ascending (`EQASND`)."
const _EKQ_ASCEND = 0
"Order sense: descending (`EQDSND`)."
const _EKQ_DESCEND = 1

"Names of the data types in the messages of the semantic check."
const _EKQ_TYPSTR = ("CHARACTER", "DOUBLE PRECISION", "INTEGER", "TIME")

# Meta-token codes of the normaliser.
const _EKQ_MRPAR = -1       # right parenthesis
const _EKQ_MLPAR = -2       # left parenthesis
const _EKQ_MBETW = -3       # BETWEEN
const _EKQ_MNOTBT = -4      # NOT BETWEEN
const _EKQ_MBEXPR = -5      # the pair of bounds of a BETWEEN
const _EKQ_MNAME = -6       # column name with its table descriptor
const _EKQ_MIDENT = -7      # identifier
const _EKQ_MVALUE = -8      # value
const _EKQ_MDOT = -9
const _EKQ_MAND = -10
const _EKQ_MOR = -11
const _EKQ_MNOT = -12
const _EKQ_MEXPR = -13      # a logical expression in normal form

# ---------------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------------

"""
    _EKQCol

A column reference of a query: the qualifier `table` as written (a table name or an
alias, upper case, `""` if none), the position `tabidx` of the table in the FROM
clause, the column name `name` (upper case), its ordinal position `colidx` in the
table, and the first and last character `lxbeg`, `lxend` of its lexemes in the query
(1-based; for a qualified name from the qualifier to the column name), by the rule of
`zzekqsel`. All fields are 0 or `""` in the unused right side of a constraint.
"""
struct _EKQCol
    table::String
    tabidx::Int
    name::String
    colidx::Int
    lxbeg::Int
    lxend::Int
end

"""
    _EKQConstraint

One relation of the normalised WHERE clause. `cnstyp` is `_EKQ_EQCOL` (1, column
operator column) or `_EKQ_EQVAL` (2, column operator value); `opcode` is one of the
codes 1 to 10 (EQ, GE, GT, LE, LT, NE, LIKE, UNLIKE, ISNULL, NOTNUL). For a value,
`dtype` is its type code after time resolution (1 CHR, 2 DP, 3 INT, 4 TIME; 0 for
ISNULL and NOTNUL and for a column on the right side), `cval` the string of type 1 as
written between the quotes (blanks kept), `dval` the stored double of types 2, 3 and 4
(for 4 the ephemeris time), and `ival` the integer of type 3 as `zzekqcon` gives it.
`dval` is the double that the number parser of CSPICE (`nparsd`) gives for the text of
the query, which can be one unit in the last place off; `dexact` is the correctly
rounded double of the same text, and the public search compares with it.
"""
struct _EKQConstraint
    cnstyp::Int
    lhs::_EKQCol
    opcode::Int
    rhs::_EKQCol
    dtype::Int
    cval::String
    dval::Float64
    ival::Int
    dexact::Float64
end

# Without `dexact` the correctly rounded value is the value of CSPICE.
_EKQConstraint(cnstyp, lhs, opcode, rhs, dtype, cval, dval, ival) =
    _EKQConstraint(cnstyp, lhs, opcode, rhs, dtype, cval, dval, ival, dval)

"""
    _EKQuery

A parsed, resolved and checked query. `tables` and `aliases` are those of the FROM
clause (upper case, cut to 64 characters, alias `""` if none) and `tabcat[i]` is the
index of table `i` in the catalogue of loaded tables. The WHERE clause is the
disjunction of `length(conjsizes)` conjunctions; conjunction `k` holds `conjsizes[k]`
constraints, and `constraints` lists them conjunction after conjunction. Invariant:
`sum(conjsizes) == length(constraints)`. A query without a WHERE clause has no
conjunction and no constraint, and it selects all rows. `orderby` holds each order-by
column with its sense (0 ascending, 1 descending).
"""
struct _EKQuery
    tables::Vector{String}
    aliases::Vector{String}
    tabcat::Vector{Int}
    conjsizes::Vector{Int}
    constraints::Vector{_EKQConstraint}
    orderby::Vector{Tuple{_EKQCol,Int}}
    select::Vector{_EKQCol}
end

# ---------------------------------------------------------------------------------
# Messages: REPMC, REPMI
# ---------------------------------------------------------------------------------

"""
    _ek_qrep(msg, value) -> String

Contract: returns `msg` with its first `#` replaced by `value`, as `repmc_` and
`repmi_`: a string loses its leading and trailing blanks and a blank string becomes
one blank. `msg` without `#` is returned as it is.
"""
function _ek_qrep(msg::String, value::Union{AbstractString,Integer})::String
    r = findfirst(==(UInt8('#')), codeunits(msg))
    r === nothing && return msg
    v = value isa Integer ? string(value) : String(strip(value, ' '))
    isempty(v) && (v = " ")
    c = codeunits(msg)
    return string(String(c[1:r-1]), v, String(c[r+1:end]))
end

"Contract: `1 <= b`; returns the bytes `b:e` of `qb` as a string (`\"\"` if `e < b`), cut at the end of `qb`."
_ek_qlex(qb::AbstractVector{UInt8}, b::Int, e::Int) = String(qb[b:min(e, length(qb))])

"Contract: returns the message `template` with the lexeme `lxb:lxe` of the query at its first `#` and then `lxb` at the next one, or the other order if `locfirst`."
function _ek_qlexmsg(template::String, qb, lxb::Int, lxe::Int, locfirst::Bool)
    locfirst && return _ek_qrep(_ek_qrep(template, lxb), _ek_qlex(qb, lxb, lxe))
    return _ek_qrep(_ek_qrep(template, _ek_qlex(qb, lxb, lxe)), lxb)
end

@noinline _ek_qsignal(kind::Symbol, msg::String) = throw(SpiceError(kind, msg))

# ---------------------------------------------------------------------------------
# The encoded query: ZZEKQINI, ZZEKINQC, ZZEKINQN
# ---------------------------------------------------------------------------------

"""
    _EKQEnc

The encoded query of CSPICE. `i[k]` is element `k` of the integer cell `EQRYI` (the
vector grows with the cardinality of the cell, at most `_EKQ_ISIZE`), `c` the character
buffer `EQRYC` and `d` the double buffer `EQRYD`. The first 19 integers are the header:
element 2 architecture, 3 initialised, 4 parsed, 5 names resolved, 6 times resolved,
7 semantically checked, 8 number of tables, 9 of constraints, 10 of conjunctions, 11 of
order-by columns, 12 of SELECT columns, 13 size of `d`, 14 first free element of `d`,
15 size of `c`, 16 first free character of `c`. After the header come 12 integers for
each table (name and alias descriptor), 26 for each constraint, one size for each
conjunction, 13 for each order-by column and 12 for each SELECT column. A descriptor is
six integers: type, first and last character of the lexeme, pointer to the value (begin
and end in `c`, or the element of `d`), and the table or column index.
"""
struct _EKQEnc
    i::Vector{Int32}
    c::Vector{UInt8}
    d::Vector{Float64}
    dx::Vector{Float64}     # the correctly rounded value of each element of `d`
end

"Contract: returns a new encoded query in the state that `zzekqini_` gives."
function _ek_qini()
    i = zeros(Int32, 19)
    i[2] = 1
    i[3] = 1
    i[4] = -1
    i[5] = -1
    i[6] = -1
    i[13] = _EKQ_MAXNUM
    i[14] = 1
    i[15] = _EKQ_CSIZE
    i[16] = 1
    sizehint!(i, 256)
    return _EKQEnc(i, fill(UInt8(' '), _EKQ_CSIZE), zeros(Float64, _EKQ_MAXNUM),
        zeros(Float64, _EKQ_MAXNUM))
end

# The functions of f2c index the cell from its control area: C element `k` is the cell
# element `k - 5`. The two accessors below keep the indices of the C source. An element
# past the cardinality reads as 0 (CSPICE reads memory that it never set there, which
# can only happen after its read of one token past a WHERE clause that ends with NOT,
# zzeknrml.c line 1426; see `tokat` in `_ek_qnrml!`).
@inline _ek_qget(q::_EKQEnc, k::Int) = k - 5 <= length(q.i) ? Int(@inbounds q.i[k-5]) : 0
@inline _ek_qput!(q::_EKQEnc, k::Int, v::Integer) = (@inbounds q.i[k-5] = Int32(v); nothing)

"Contract: appends `x` to the integer cell as `appndi_`; throws `SpiceError(:celltoosmall)` if the cell is full."
@inline function _ek_qappnd!(q::_EKQEnc, x::Integer)
    length(q.i) < _EKQ_ISIZE || _ek_qsignal(:celltoosmall,
        "The cell cannot accommodate the addition of the element $x")
    push!(q.i, Int32(x))
    return nothing
end

"Contract: appends the six integers of the descriptor `d`."
@inline function _ek_qappnd!(q::_EKQEnc, d::NTuple{6,Int})
    for x in d
        _ek_qappnd!(q, x)
    end
    return nothing
end

"The cleared descriptor with type 1 that stands for a missing table qualifier or alias."
const _EKQ_NODESCR = (1, 0, 0, 0, 0, 0)

"""
    _ek_qinqc!(q, value, lxb, lxe) -> NTuple{6,Int}

Contract: `value` is not empty. Stores the bytes `value` in the character buffer of `q`
and returns the descriptor of the string with the lexeme bounds `lxb`, `lxe`
(`zzekinqc_`). Throws `SpiceError(:buffertoosmall)` if the buffer is full.
"""
function _ek_qinqc!(q::_EKQEnc, value::AbstractVector{UInt8}, lxb::Int, lxe::Int)
    n = length(value)
    n < 1 && _ek_qsignal(:invalidcount, "Length of string value was $n; must be > 0.")
    free = Int(q.i[16])
    room = Int(q.i[15]) - free + 1
    n > room && _ek_qsignal(:buffertoosmall, "Out of room in character portion of " *
        "encoded query; only $room elements were available; $n are needed.")
    copyto!(q.c, free, value, firstindex(value), n)
    q.i[16] = free + n
    return (1, lxb, lxe, free, free + n - 1, 0)
end

"""
    _ek_qinqn!(q, value, type, lxb, lxe) -> NTuple{6,Int}

Contract: stores the number `value` of type code `type` in the double buffer of `q` and
returns its descriptor (`zzekinqn_`); `exact` is the correctly rounded value of the same
text. Throws `SpiceError(:buffertoosmall)` if the
buffer is full.
"""
function _ek_qinqn!(q::_EKQEnc, value::Float64, type::Int, lxb::Int, lxe::Int,
        exact::Float64=value)
    free = Int(q.i[14])
    size = Int(q.i[13])
    size - free + 1 <= 0 && _ek_qsignal(:buffertoosmall, "Out of room in numeric " *
        "portion of encoded query; only $size elements were available.")
    q.d[free] = value
    q.dx[free] = exact
    q.i[14] = free + 1
    return (type, lxb, lxe, free, 0, 0)
end

"Contract: `1 <= b <= e` are bounds in the character buffer; returns the string `b:e`."
_ek_qstr(q::_EKQEnc, b::Int, e::Int) = String(q.c[b:e])

# ---------------------------------------------------------------------------------
# Scanner: ZZEKSCAN with LXQSTR, PARSQS, LX4NUM, LXIDNT, BEINT
# ---------------------------------------------------------------------------------

"""
    _EKQTokens

The tokens of a query: `tokens[k]` is the token code, `values[k]` the keyword code of a
keyword, the index into `numvls` of a number, or the index into `chbegs`/`chends` of an
identifier or a quoted string (their text is `chrbuf[chbegs[v]:chends[v]]`; identifiers
are in upper case), and `lxbegs[k]`, `lxends[k]` are the first and last character of the
lexeme in the query. Invariant: all vectors of tokens have length `n`.
"""
struct _EKQTokens
    n::Int
    tokens::Vector{Int}
    values::Vector{Int}
    lxbegs::Vector{Int}
    lxends::Vector{Int}
    numvls::Vector{Float64}
    numexact::Vector{Float64}
    chrbuf::Vector{UInt8}
    chbegs::Vector{Int}
    chends::Vector{Int}
end

@inline _ek_qisdigit(c::UInt8) = UInt8('0') <= c <= UInt8('9')
@inline _ek_qisalpha(c::UInt8) =
    (UInt8('A') <= c <= UInt8('Z')) || (UInt8('a') <= c <= UInt8('z'))
@inline _ek_qisidtail(c::UInt8) =
    _ek_qisalpha(c) || _ek_qisdigit(c) || c == UInt8('$') || c == UInt8('_')
@inline _ek_qupper(c::UInt8) = (UInt8('a') <= c <= UInt8('z')) ? c - 0x20 : c

"Contract: returns `(last, nchar)` of the unsigned integer that starts at `s[first]`, as `lx4uns_`; `nchar == 0` and `last == first - 1` if there is none."
function _ek_qlx4uns(s::Vector{UInt8}, first::Int)
    l = length(s)
    last = first - 1
    (first < 1 || first > l) && return last, 0
    i = first
    while i <= l && _ek_qisdigit(@inbounds s[i])    # invariant: s[first:last] are digits
        last += 1
        i += 1
    end
    return last, last + 1 - first
end

"Contract: as `_ek_qlx4uns` for a signed integer (`lx4sgn_`)."
function _ek_qlx4sgn(s::Vector{UInt8}, first::Int)
    (first < 1 || first > length(s)) && return first - 1, 0
    c = s[first]
    if c == UInt8('+') || c == UInt8('-')
        last, n = _ek_qlx4uns(s, first + 1)
        return n == 0 ? (first - 1, 0) : (last, n + 1)
    end
    return _ek_qlx4uns(s, first)
end

"Contract: as `_ek_qlx4uns` for a decimal number (`lx4dec_`); `nchar` is only exact for the forms that start with a digit or a sign and a digit, as in CSPICE."
function _ek_qlx4dec(s::Vector{UInt8}, first::Int)
    l = length(s)
    (first < 1 || first > l) && return first - 1, 0
    next = first + 1
    i = s[first]
    j = next < l ? s[next] : UInt8(' ')
    if i == UInt8('.')
        last, n = _ek_qlx4uns(s, first + 1)
        return n == 0 ? (first - 1, 0) : (last, n + 1)
    elseif (i == UInt8('-') || i == UInt8('+')) && j == UInt8('.')
        last, n = _ek_qlx4uns(s, next + 1)
        return n == 0 ? (first - 1, 0) : (last, n + 1)
    end
    last, n = _ek_qlx4sgn(s, first)
    if n > 0 && last < l && s[last+1] == UInt8('.')
        last, _ = _ek_qlx4uns(s, last + 2)
        n = last + 1 - first
    end
    return last, n
end

"Contract: as `_ek_qlx4uns` for a number with an optional exponent `E`, `D`, `e` or `d` (`lx4num_`)."
function _ek_qlx4num(s::Vector{UInt8}, first::Int)
    l = length(s)
    (first < 1 || first > l) && return first - 1, 0
    last, n = _ek_qlx4dec(s, first)
    if n > 0 && last < l
        c = s[last+1]
        if c == UInt8('e') || c == UInt8('E') || c == UInt8('D') || c == UInt8('d')
            temp, m = _ek_qlx4sgn(s, last + 2)
            if m > 0
                last = temp
                n = last + 1 - first
            end
        end
    end
    return last, n
end

"""
    _ek_qlxqstr(s, qchar, first) -> (last, nchar)

Contract: finds the quoted string that starts at `s[first]` with the quote character
`qchar`, as `lxqstr_`: the string ends at the first quote that is not doubled, and a
quote without a partner ends it at the last complete pair. `nchar == 0` if `s[first]`
is not `qchar` or no second quote follows.
"""
function _ek_qlxqstr(s::Vector{UInt8}, qchar::UInt8, first::Int)
    l = length(s)
    (first < 1 || first > l - 1 || qchar == UInt8(' ') || s[first] != qchar) &&
        return first - 1, 0
    loc = 0
    for k in first+1:l
        if s[k] == qchar
            loc = k - first
            break
        end
    end
    loc == 0 && return first - 1, 0
    last = first + loc
    even = true
    pos = last + 1
    while pos <= l      # invariant: s[first:last] is a complete quoted string
        if s[pos] == qchar
            even = !even
            even && (last = pos)
        elseif even
            return last, last - first + 1
        end
        pos += 1
    end
    return last, last - first + 1
end

"""
    _ek_qparsqs(s, qchar, outlen) -> (value, error)

Contract: `s` is a quoted string token with the quote character `qchar`. Returns the
bytes between the quotes with each doubled quote as one, or the error message of
`parsqs_` (`""` if none); `outlen` is the room for the value.
"""
function _ek_qparsqs(s::AbstractVector{UInt8}, qchar::UInt8, outlen::Int)
    value = UInt8[]
    first = findfirst(!=(UInt8(' ')), s)
    first === nothing && return value, "Blank input string"
    last = findlast(!=(UInt8(' ')), s)
    s[first] != qchar && return value, "String token does not start with quote character"
    s[last] != qchar && return value, "String token does not end with quote character"
    first == last - 1 && return value, "Null (zero length) string token"
    opos = 1
    ipos = first + 1
    while ipos <= last - 1 && opos <= outlen
        chr = s[ipos]
        if chr != qchar
            push!(value, chr)
            ipos += 1
            opos += 1
        elseif ipos == last - 1
            return value, "Quote character is unmatched or else string ends without " *
                          "final quote; take your pick"
        elseif s[ipos+1] != qchar
            return value, "Interior quote character is not doubled"
        else
            push!(value, chr)
            opos += 1
            ipos += 2
        end
    end
    ipos < last - 1 && return value, "Output string too short, truncated on right"
    return value, ""
end

"""
    _ek_qexact(lexeme, x) -> Float64

Contract: `lexeme` is a numeric token that `nparsd` reads as `x`. Returns the correctly
rounded double of the token, or `x` if the token has a form that Julia does not read
(the result then is the value of CSPICE). The two values differ by a few units in the
last place at most, so a result farther from `x` is not accepted.
"""
function _ek_qexact(lexeme::AbstractVector{UInt8}, x::Float64)::Float64
    text = replace(String(copy(lexeme)), r"[dD]" => "e", ' ' => "", ',' => "")
    v = tryparse(Float64, text)
    (v === nothing || !isfinite(v)) && return x
    return abs(v - x) <= 8 * eps(abs(x)) ? v : x
end

"""
    _ek_qnparsd(s) -> (value, error, ptr)

Contract: `s` is a numeric token of `_ek_qlx4num`. Returns its value with an empty
`error`, or the message and the pointer of `nparsd_` for a number that is too large
for a double (the only error that such a token can have).
"""
function _ek_qnparsd(s::AbstractVector{UInt8})
    try
        return _nparsd(s), "", 0
    catch err
        err isa SpiceError || rethrow()
    end
    # The pointer of NPARSD: the exponent digit at which the exponent, with the count
    # of the mantissa digits that did not fit, passes 308; or one past the end.
    intval = decval = ecount = expval = 0.0
    sgn2 = 1.0
    mode = 1                    # 1 integer part, 2 fraction, 3 exponent
    zeroi = false
    ptr = length(s) + 1
    for (i, c) in enumerate(s)
        if _ek_qisdigit(c)
            v = Float64(c - UInt8('0'))
            if mode == 1
                intval < 1e15 ? (intval = intval * 10.0 + v) : (ecount += 1.0)
            elseif mode == 2
                if zeroi && decval < 1e15
                    decval = decval * 10.0 + v
                    ecount -= 1.0
                end
            elseif expval + ecount > 308.0
                ptr = i
                break
            elseif expval + ecount >= -309.0
                expval = expval * 10.0 + sgn2 * v
            end
        elseif c == UInt8('.')
            mode = 2
            zeroi = intval == 0.0
        elseif c == UInt8('-') && mode == 3
            sgn2 = -1.0
        elseif c != UInt8('+') && c != UInt8('-')
            mode = 3
        end
    end
    return 0.0, "The number represented by the input string is too large to be " *
                "stored as a double precision number.", ptr
end

"""
    _ek_qscan(qb, maxntk, maxnum, buflen) -> (tokens, error)

Contract: `qb` holds the bytes of a query of at most `_EKQ_MAXQRY` characters without
trailing blanks. Returns the tokens of the query (`zzekscan_`) and `""`, or `nothing`
and the error message of the scanner; `maxntk` and `maxnum` are the largest numbers of
tokens and of numeric tokens and `buflen` the room for identifiers and strings.
"""
function _ek_qscan(qb::AbstractVector{UInt8}, maxntk::Int, maxnum::Int, buflen::Int)
    l = max(1, something(findlast(!=(UInt8(' ')), qb), 0))
    # TQUERY of ZZEKSCAN: 2000 characters, padded with blanks. The lexers see all of it.
    tq = fill(UInt8(' '), _EKQ_MAXQRY)
    copyto!(tq, 1, qb, 1, min(l, length(qb), _EKQ_MAXQRY))
    tokens = Int[]
    values = Int[]
    lxbegs = Int[]
    lxends = Int[]
    numvls = Float64[]
    numexact = Float64[]
    chrbuf = UInt8[]
    chbegs = Int[]
    chends = Int[]
    cptr = 1
    function token!(tok, val, b, e)
        push!(tokens, tok)
        push!(values, val)
        push!(lxbegs, b)
        push!(lxends, e)
        return nothing
    end
    while true
        # invariant: the tokens found so far cover the query up to lxends[end], and
        # length(chrbuf) is CHCARD, the number of characters stored
        ntoken = length(tokens)
        ntoken > 0 && (cptr = lxends[ntoken] + 1)
        cptr > l && break
        while tq[cptr] == UInt8(' ') && cptr < l
            cptr += 1
        end
        tq[cptr] == UInt8(' ') && break
        chr = tq[cptr]
        if ntoken == maxntk
            return nothing, _ek_qrep(_ek_qrep("Maximum allowed number of tokens is #; " *
                "at least # tokens are present in QUERY.", maxntk), maxntk + 1)
        end
        kind = 8
        if chr == UInt8('\'') || chr == UInt8('"')
            kind = 4
        elseif chr == UInt8('.')
            _, nchars = _ek_qlx4num(tq, cptr)
            kind = nchars > 0 ? 6 : 8
        elseif _ek_qisdigit(chr) || chr == UInt8('+') || chr == UInt8('-')
            kind = 6
        elseif _ek_qisalpha(chr)
            kind = 7
        end
        chcard = length(chrbuf)
        if kind == 4                                    # quoted string
            last, nchars = _ek_qlxqstr(tq, chr, cptr)
            nchars == 0 &&
                return nothing, _ek_qrep("Invalid quoted string at location #.", cptr)
            room = buflen - chcard
            if nchars > room
                msg = "Insufficient space to store quoted string at location #; # " *
                      "chars needed; only # are available."
                return nothing, _ek_qrep(_ek_qrep(_ek_qrep(msg, cptr), nchars), room)
            end
            value, perr = _ek_qparsqs(view(tq, cptr:cptr+nchars-1), chr, room)
            if !isempty(perr)
                msg = "Error occurred while parsing quoted string token at location #:"
                return nothing, string(_ek_qrep(msg, cptr), "  ", perr)
            end
            append!(chrbuf, value)
            push!(chbegs, chcard + 1)
            push!(chends, chcard + length(value))
            token!(_EKQ_TKQSTR, length(chbegs), cptr, last)
        elseif kind == 6                                # number
            last, nchars = _ek_qlx4num(tq, cptr)
            nchars == 0 &&
                return nothing, _ek_qrep("Invalid numeric token at location #.", cptr)
            if maxnum - length(numvls) < 1
                msg = "Insufficient space to store value of number at location #; # " *
                      "elements are available in the NUMVLS array; # are required."
                return nothing, _ek_qrep(_ek_qrep(_ek_qrep(msg, cptr), maxnum), maxnum + 1)
            end
            lexeme = view(tq, cptr:last)
            x, nerr, ptr = _ek_qnparsd(lexeme)
            if !isempty(nerr)
                msg = "Error found in numeric token at location #:"
                return nothing, string(_ek_qrep(msg, cptr + ptr - 1), "  ", nerr)
            end
            push!(numvls, x)
            push!(numexact, _ek_qexact(lexeme, x))
            # BEINT: an optional sign and then digits only.
            b = (lexeme[1] == UInt8('+') || lexeme[1] == UInt8('-')) ? 2 : 1
            isint = b <= length(lexeme) && all(_ek_qisdigit, view(lexeme, b:length(lexeme)))
            token!(isint ? _EKQ_TKINT : _EKQ_TKDP, length(numvls), cptr, last)
        elseif kind == 7                                # identifier or keyword
            last = cptr
            while last < l && _ek_qisidtail(tq[last+1])
                last += 1
            end
            nchars = last - cptr + 1
            room = buflen - chcard
            if nchars > room
                msg = "Insufficient space to store identifier string at location #; # " *
                      "chars needed; only # are available."
                return nothing, _ek_qrep(_ek_qrep(_ek_qrep(msg, cptr), nchars), room)
            end
            name = map(_ek_qupper, view(tq, cptr:last))
            kw = nchars <= 8 ? findfirst(==(String(copy(name))), _EKQ_KEYWORDS) : nothing
            if kw !== nothing
                token!(_EKQ_TKKEY, kw, cptr, last)
            else
                append!(chrbuf, name)
                push!(chbegs, chcard + 1)
                push!(chends, chcard + nchars)
                token!(_EKQ_TKID, length(chbegs), cptr, last)
            end
        else                                            # special character
            found = false
            for n in min(2, l - cptr + 1):-1:1, (str, tok, val) in _EKQ_SPECIAL
                # A candidate of n characters equals a special string with blank padding.
                c = codeunits(str)
                same = tq[cptr] == c[1] &&
                       (n == 1 ? length(c) == 1 :
                        tq[cptr+1] == (length(c) == 2 ? c[2] : UInt8(' ')))
                if same
                    token!(tok, val, cptr, cptr - 1 + length(c))
                    found = true
                    break
                end
            end
            if !found
                msg = _ek_qrep("Invalid character found at location #.", cptr)
                if 33 <= chr <= 126
                    msg = _ek_qrep(msg * "  <character> = '#'", String([chr]))
                else
                    msg = _ek_qrep(msg * "  ICHAR(<character>) = #", Int(chr))
                end
                return nothing, msg
            end
        end
    end
    return _EKQTokens(length(tokens), tokens, values, lxbegs, lxends, numvls, numexact, chrbuf,
        chbegs, chends), ""
end

"Contract: returns `(loc, found)` of the first token with the code `tokid` (and the keyword code `kwcode` for a keyword), as `zzektloc_`; `loc == n + 1` if there is none."
function _ek_qtloc(tk::_EKQTokens, tokid::Int, kwcode::Int)
    for loc in 1:tk.n
        if tk.tokens[loc] == tokid && (tokid != _EKQ_TKKEY || tk.values[loc] == kwcode)
            return loc, true
        end
    end
    return tk.n + 1, false
end

"Contract: `1 <= k <= tk.n` is an identifier or a quoted string token; returns its text."
@inline function _ek_qtext(tk::_EKQTokens, k::Int)
    v = tk.values[k]
    return view(tk.chrbuf, tk.chbegs[v]:tk.chends[v])
end

# ---------------------------------------------------------------------------------
# Doubly linked list pools: LNKINI, LNKAN, LNKILA, LNKILB, LNKFSL, LNKHL, LNKTL, ...
# ---------------------------------------------------------------------------------

"""
    _EKQPool

A pool of doubly linked lists as SPICELIB keeps it: node `n` has the forward pointer
`fwd[n]` and the backward pointer `bck[n]`. In a list the backward pointer of the head
is minus the tail and the forward pointer of the tail is minus the head; a free node
has the backward pointer 0 and the free nodes form a list through `fwd` from `free`.
Invariant: `nfree` is the number of free nodes.
"""
mutable struct _EKQPool
    const fwd::Vector{Int32}
    const bck::Vector{Int32}
    nfree::Int
    free::Int
end

"Contract: `size >= 1`; returns a pool of `size` free nodes (`lnkini_`)."
function _ek_qlnkini(size::Int)
    fwd = Int32.(2:size+1)
    fwd[size] = 0
    return _EKQPool(fwd, zeros(Int32, size), size, 1)
end

@noinline _ek_qlnkbad(node, p::_EKQPool) = _ek_qsignal(:invalidnode,
    "NODE was $node; valid range is 1 to $(length(p.fwd)).")
@noinline _ek_qlnkunalloc(node) = _ek_qsignal(:unallocatednode,
    "Node $node of the linked list pool is not allocated.")

"Contract: throws `SpiceError` unless `node` is an allocated node of `p`."
@inline function _ek_qlnkchk(node::Int, p::_EKQPool)
    (node < 1 || node > length(p.fwd)) && _ek_qlnkbad(node, p)
    p.bck[node] == 0 && _ek_qlnkunalloc(node)
    return nothing
end

"Contract: allocates a node of `p` as a list of one node and returns it (`lnkan_`)."
function _ek_qlnkan!(p::_EKQPool)
    p.nfree == 0 && _ek_qsignal(:nofreenodes,
        "There are no free nodes left for allocating in the supplied linked list pool.")
    new = p.free
    p.free = p.fwd[new]
    p.nfree -= 1
    p.fwd[new] = -new
    p.bck[new] = -new
    return new
end

"Contract: returns the forward pointer of `node` (`lnknxt_`): the next node, or minus the head at the tail."
@inline _ek_qlnknxt(node::Int, p::_EKQPool) = (_ek_qlnkchk(node, p); Int(p.fwd[node]))
"Contract: returns the backward pointer of `node` (`lnkprv_`): the previous node, or minus the tail at the head."
@inline _ek_qlnkprv(node::Int, p::_EKQPool) = (_ek_qlnkchk(node, p); Int(p.bck[node]))

"Contract: returns the head of the list of `node`, or 0 if `node < 1` (`lnkhl_`)."
function _ek_qlnkhl(node::Int, p::_EKQPool)
    node < 1 && return 0
    _ek_qlnkchk(node, p)
    while p.bck[node] > 0
        node = Int(p.bck[node])
    end
    return node
end

"Contract: returns the tail of the list of `node`, or 0 if `node < 1` (`lnktl_`)."
function _ek_qlnktl(node::Int, p::_EKQPool)
    node < 1 && return 0
    _ek_qlnkchk(node, p)
    while p.fwd[node] > 0
        node = Int(p.fwd[node])
    end
    return node
end

"Contract: inserts the whole list of the node `list` after the node `prev`; nothing if `prev <= 0` (`lnkila_`)."
function _ek_qlnkila!(prev::Int, list::Int, p::_EKQPool)
    prev <= 0 && return nothing
    _ek_qlnkchk(prev, p)
    _ek_qlnkchk(list, p)
    head = _ek_qlnkhl(list, p)
    tail = -Int(p.bck[head])
    next = Int(p.fwd[prev])
    p.fwd[prev] = head
    p.bck[head] = prev
    if next > 0
        p.bck[next] = tail
    else
        p.bck[-next] = -tail
    end
    p.fwd[tail] = next
    return nothing
end

"Contract: inserts the whole list of the node `list` before the node `next`; nothing if `next <= 0` (`lnkilb_`)."
function _ek_qlnkilb!(list::Int, next::Int, p::_EKQPool)
    next <= 0 && return nothing
    _ek_qlnkchk(next, p)
    _ek_qlnkchk(list, p)
    head = _ek_qlnkhl(list, p)
    tail = -Int(p.bck[head])
    prev = Int(p.bck[next])
    p.fwd[tail] = next
    p.bck[next] = tail
    if prev > 0
        p.fwd[prev] = head
    else
        p.fwd[-prev] = -head
    end
    p.bck[head] = prev
    return nothing
end

"Contract: `tail` follows `head` in one list; frees the nodes from `head` to `tail` and joins their neighbours (`lnkfsl_`)."
function _ek_qlnkfsl!(head::Int, tail::Int, p::_EKQPool)
    _ek_qlnkchk(head, p)
    _ek_qlnkchk(tail, p)
    count = 1
    node = head
    while node != tail && node > 0
        count += 1
        node = Int(p.fwd[node])
    end
    node != tail && _ek_qsignal(:invalidsublist,
        "Node $tail cannot be found by forward traversal, starting at node $head.")
    prev = Int(p.bck[head])
    next = Int(p.fwd[tail])
    if prev > 0
        p.fwd[prev] = next
        next <= 0 && (p.bck[-next] = -prev)
    end
    if next > 0
        p.bck[next] = prev
        prev <= 0 && (p.fwd[-prev] = -next)
    end
    node = head
    while node != next
        p.bck[node] = 0
        node = Int(p.fwd[node])
    end
    p.fwd[tail] = p.free
    p.free = head
    p.nfree += count
    return nothing
end

# ---------------------------------------------------------------------------------
# WHERE clause: ZZEKNRML
# ---------------------------------------------------------------------------------

"""
    _ek_qnrml!(q, qb, tk) -> Union{Nothing,String}

Contract: `tk` holds the tokens of the query `qb` and has a WHERE keyword, and `q`
holds the FROM clause. Appends the constraints of the WHERE clause in normal form (a
disjunction of conjunctions) and the sizes of the conjunctions to `q`, as
`zzeknrml_`, and returns `nothing`; or returns the error message of the normaliser.

The WHERE clause runs from the token after the first WHERE to the token before the
first FROM, ORDER or SELECT that follows it. Its tokens become "meta-tokens" in a
linked list, which a stack machine reduces: an identifier or `table.column` becomes a
name, `name operator value-or-name` a relation (an expression with one conjunction of
one relation), `value AND value` the bounds of a BETWEEN, `NOT expression` the negated
expression, `expression AND expression` the distributed product, `expression OR
expression` the joined lists. AND binds stronger than OR, and NOT stronger than AND.
"""
function _ek_qnrml!(q::_EKQEnc, qb::AbstractVector{UInt8}, tk::_EKQTokens)
    ntoken = tk.n
    tokens, values = tk.tokens, tk.values
    ntoken > _EKQ_MAXTOK &&
        return _ek_qrep("Too many tokens in query; max allowed is #.", _EKQ_MAXTOK)
    whrbeg, fnd = _ek_qtloc(tk, _EKQ_TKKEY, _EKQ_KWWHER)
    whrbeg += 1
    fnd || return "Missing WHERE keyword."
    whrend = ntoken
    for kw in (_EKQ_KWFROM, _EKQ_KWORDR, _EKQ_KWSEL)
        endloc, found = _ek_qtloc(tk, _EKQ_TKKEY, kw)
        found && endloc < whrend && endloc > whrbeg && (whrend = endloc - 1)
    end
    whrsiz = whrend - whrbeg + 1
    whrsiz == 0 && return "Empty WHERE clause."

    rlpool = _ek_qlnkini(_EKQ_NPOOL)        # relations
    cjpool = _ek_qlnkini(_EKQ_NPOOL)        # conjunctions: lists of relation pointers
    djpool = _ek_qlnkini(_EKQ_NPOOL)        # disjunctions: lists of conjunction pointers
    mtpool = _ek_qlnkini(_EKQ_MAXTOK)       # meta-tokens
    dspool = _ek_qlnkini(_EKQ_NPOOL)        # descriptors
    mtcode = zeros(Int, _EKQ_MAXTOK)
    mtexpr = zeros(Int, _EKQ_MAXTOK)
    dsc = zeros(Int, 7, 2 * _EKQ_MAXTOK + 2)    # descriptor and its meta code, per node
    rels = zeros(Int, 3, _EKQ_NPOOL)            # left name, operator, right name or value
    cjptrs = zeros(Int, _EKQ_NPOOL)
    djptrs = zeros(Int, _EKQ_NPOOL)
    mstart = zeros(Int, _EKQ_MAXTOK)
    popcnd = zeros(Int, _EKQ_MAXTOK)
    tokat(k) = k <= ntoken ? tokens[k] : 0      # CSPICE reads past the last token here
    valat(k) = k <= ntoken ? values[k] : 0
    setdsc!(node, d, code) = (for m in 1:6
            dsc[m, node] = d[m]
        end; dsc[7, node] = code; nothing)

    # Meta-tokens of the WHERE clause.
    nmeta = 0
    tail = 0
    i = whrbeg
    while i <= whrend       # invariant: the meta-tokens of tokens whrbeg:i-1 are in the list
        node = _ek_qlnkan!(mtpool)
        _ek_qlnkila!(tail, node, mtpool)
        tail = node
        mtexpr[node] = i
        t = tokens[i]
        if t == _EKQ_TKLPAR
            mtcode[node] = _EKQ_MLPAR
        elseif t == _EKQ_TKRPAR
            mtcode[node] = _EKQ_MRPAR
        elseif t == _EKQ_TKINT || t == _EKQ_TKDP
            mtcode[node] = _EKQ_MVALUE
            mtexpr[node] = _ek_qlnkan!(dspool)
            d = _ek_qinqn!(q, tk.numvls[values[i]], t == _EKQ_TKINT ? 3 : 2,
                tk.lxbegs[i], tk.lxends[i], tk.numexact[values[i]])
            setdsc!(mtexpr[node], d, _EKQ_MVALUE)
        elseif t == _EKQ_TKQSTR || t == _EKQ_TKID
            code = t == _EKQ_TKQSTR ? _EKQ_MVALUE : _EKQ_MIDENT
            mtcode[node] = code
            mtexpr[node] = _ek_qlnkan!(dspool)
            d = _ek_qinqc!(q, _ek_qtext(tk, i), tk.lxbegs[i], tk.lxends[i])
            setdsc!(mtexpr[node], d, code)
        elseif t == _EKQ_TKDOT
            mtcode[node] = _EKQ_MDOT
        elseif t == _EKQ_TKKEY
            v = values[i]
            if v == _EKQ_KWAND
                mtcode[node] = _EKQ_MAND
            elseif v == _EKQ_KWOR
                mtcode[node] = _EKQ_MOR
            elseif v == _EKQ_KWNOT
                if tokat(i + 1) == _EKQ_TKKEY && valat(i + 1) == _EKQ_KWLIKE
                    mtcode[node] = _EKQ_UNLIKE
                    i += 1
                elseif tokat(i + 1) == _EKQ_TKKEY && valat(i + 1) == _EKQ_KWBETW
                    mtcode[node] = _EKQ_MNOTBT
                    i += 1
                else
                    mtcode[node] = _EKQ_MNOT
                end
            elseif v == _EKQ_KWEQ
                mtcode[node] = _EKQ_EQ
            elseif v == _EKQ_KWGE
                mtcode[node] = _EKQ_GE
            elseif v == _EKQ_KWGT
                mtcode[node] = _EKQ_GT
            elseif v == _EKQ_KWLE
                mtcode[node] = _EKQ_LE
            elseif v == _EKQ_KWLT
                mtcode[node] = _EKQ_LT
            elseif v == _EKQ_KWNE
                mtcode[node] = _EKQ_NE
            elseif v == _EKQ_KWLIKE
                mtcode[node] = _EKQ_LIKE
            elseif v == _EKQ_KWBETW
                mtcode[node] = _EKQ_MBETW
            elseif v == _EKQ_KWIS
                # IS is EQ and IS NOT is NE; a NULL on the right makes them ISNULL, NOTNUL.
                if i < whrend && tokens[i+1] == _EKQ_TKKEY && values[i+1] == _EKQ_KWNOT
                    mtcode[node] = _EKQ_NE
                    i += 1
                else
                    mtcode[node] = _EKQ_EQ
                end
            elseif v == _EKQ_KWNULL
                mtcode[node] = _EKQ_MVALUE
                mtexpr[node] = 0
            else
                return _ek_qlexmsg("Unexpected keyword # found at location #.", qb,
                    tk.lxbegs[i], tk.lxends[i], false)
            end
        else
            return _ek_qlexmsg("Unexpected token # found at location #.", qb,
                tk.lxbegs[i], tk.lxends[i], false)
        end
        i += 1
        nmeta += 1
    end

    # Remove the parentheses around a single value.
    metahd = _ek_qlnkhl(tail, mtpool)
    node = metahd
    while node > 0
        if mtcode[node] == _EKQ_MNAME || mtcode[node] == _EKQ_MVALUE
            prev = _ek_qlnkprv(node, mtpool)
            next = _ek_qlnknxt(node, mtpool)
            if prev > 0 && next > 0 && mtcode[prev] == _EKQ_MLPAR &&
               mtcode[next] == _EKQ_MRPAR
                _ek_qlnkfsl!(prev, prev, mtpool)
                _ek_qlnkfsl!(next, next, mtpool)
                metahd = _ek_qlnkhl(node, mtpool)
                nmeta -= 2
            else
                node = next
            end
        else
            node = _ek_qlnknxt(node, mtpool)
        end
    end
    whrsiz < 3 && return "Too few tokens in WHERE clause."

    full_rel = "Relation table is full."
    full_cj = "Conjunction table is full."
    full_dj = "Disjunction table is full."
    level = 1
    mstart[1] = metahd
    popcnd[1] = 0
    state = 2                   # 0 push, 1 pop, 2 parse, 3 reduced, 4 done
    first = second = third = fourth = 0
    start = retcnd = 0
    newdj = 0
    while state != 4
        if state == 2
            first = mstart[level]
            second = first > 0 ? _ek_qlnknxt(first, mtpool) : 0
            third = second > 0 ? _ek_qlnknxt(second, mtpool) : 0
            fourth = third > 0 ? _ek_qlnknxt(third, mtpool) : 0
            first <= 0 && return "WHERE clause ran out of tokens unexpectedly.  This " *
                                 "may be due to an extra left parenthesis."
            c1 = mtcode[first]
            c2 = second > 0 ? mtcode[second] : 0
            c3 = third > 0 ? mtcode[third] : 0
            if c1 == _EKQ_MIDENT
                # An identifier: the table of `table.column`, or a column without table.
                if third > 0 && c2 == _EKQ_MDOT && c3 == _EKQ_MIDENT
                    _ek_qlnkila!(mtexpr[first], mtexpr[third], dspool)
                    mtcode[first] = _EKQ_MNAME
                    _ek_qlnkfsl!(second, third, mtpool)
                    nmeta -= 2
                else
                    tabptr = _ek_qlnkan!(dspool)
                    setdsc!(tabptr, _EKQ_NODESCR, _EKQ_MIDENT)
                    _ek_qlnkila!(tabptr, mtexpr[first], dspool)
                    mtexpr[first] = tabptr
                    mtcode[first] = _EKQ_MNAME
                end
                state = 3
            elseif c1 == _EKQ_MVALUE
                third <= 0 && return "Tokens were missing from comparison relation."
                (mtexpr[first] == 0 || mtexpr[third] == 0) && return "NULL values are " *
                    "not allowed in BETWEEN or NOT BETWEEN clauses."
                if c3 == _EKQ_MIDENT
                    start = third
                    retcnd = 1
                    state = 0
                elseif c2 == _EKQ_MAND && (c3 == _EKQ_MNAME || c3 == _EKQ_MVALUE)
                    _ek_qlnkilb!(mtexpr[first], mtexpr[third], dspool)
                    mtcode[first] = _EKQ_MBEXPR
                    _ek_qlnkfsl!(second, third, mtpool)
                    nmeta -= 2
                    state = 3
                elseif c2 > 0
                    return "Token sequence must be followed by a value."
                else
                    return "Token must be followed by the AND operator."
                end
            elseif c1 == _EKQ_MNAME
                third <= 0 && return "Tokens were missing from comparison relation."
                if c3 == _EKQ_MIDENT
                    start = third
                    retcnd = 1
                    state = 0
                elseif c2 == _EKQ_MAND && (c3 == _EKQ_MNAME || c3 == _EKQ_MVALUE)
                    (mtexpr[first] == 0 || mtexpr[third] == 0) && return "NULL values " *
                        "are not allowed in BETWEEN or NOT BETWEEN clauses."
                    _ek_qlnkilb!(mtexpr[first], mtexpr[third], dspool)
                    mtcode[first] = _EKQ_MBEXPR
                    _ek_qlnkfsl!(second, third, mtpool)
                    nmeta -= 2
                    state = 3
                elseif c2 > 0 && (c3 == _EKQ_MNAME || c3 == _EKQ_MVALUE)
                    # A relation: an expression of one conjunction of one relation.
                    rlpool.nfree < 1 && return full_rel
                    newrel = _ek_qlnkan!(rlpool)
                    rels[1, newrel] = mtexpr[first]
                    rels[2, newrel] = c2
                    rels[3, newrel] = mtexpr[third]
                    _ek_qlnkfsl!(second, third, mtpool)
                    nmeta -= 2
                    cjpool.nfree < 1 && return full_cj
                    newcj = _ek_qlnkan!(cjpool)
                    cjptrs[newcj] = newrel
                    djpool.nfree < 1 && return full_dj
                    newdj = _ek_qlnkan!(djpool)
                    djptrs[newdj] = newcj
                    mtcode[first] = _EKQ_MEXPR
                    mtexpr[first] = newdj
                    state = 3
                elseif c2 == _EKQ_MBETW || c2 == _EKQ_MNOTBT
                    if c3 == _EKQ_MBEXPR
                        # BETWEEN lo AND hi:      (>= lo AND <= hi) OR (<= lo AND >= hi)
                        # NOT BETWEEN lo AND hi:  (<  lo AND <  hi) OR (>  lo AND >  hi)
                        rel = zeros(Int, 4)
                        cj = zeros(Int, 4)
                        dj = zeros(Int, 2)
                        for m in 1:4
                            k = mtexpr[third]
                            if m == 2 || m == 4
                                skip = dsc[7, k] == _EKQ_MIDENT ? 1 : 0
                                for _ in 1:skip+1
                                    k = _ek_qlnknxt(k, dspool)
                                end
                            end
                            rlpool.nfree < 1 && return full_rel
                            rel[m] = _ek_qlnkan!(rlpool)
                            rels[1, rel[m]] = mtexpr[first]
                            rels[2, rel[m]] = if c2 == _EKQ_MBETW
                                (m == 1 || m == 4) ? _EKQ_GE : _EKQ_LE
                            else
                                m <= 2 ? _EKQ_LT : _EKQ_GT
                            end
                            rels[3, rel[m]] = k
                            cjpool.nfree < 1 && return full_cj
                            cj[m] = _ek_qlnkan!(cjpool)
                            cjptrs[cj[m]] = rel[m]
                        end
                        _ek_qlnkila!(cj[1], cj[2], cjpool)
                        _ek_qlnkila!(cj[3], cj[4], cjpool)
                        for m in 1:2
                            djpool.nfree < 1 && return full_dj
                            dj[m] = _ek_qlnkan!(djpool)
                            djptrs[dj[m]] = cj[2*m-1]
                        end
                        _ek_qlnkila!(dj[1], dj[2], djpool)
                        mtcode[first] = _EKQ_MEXPR
                        mtexpr[first] = dj[1]
                        _ek_qlnkfsl!(second, third, mtpool)
                        nmeta -= 2
                        state = 3
                    elseif c3 == _EKQ_MNAME || c3 == _EKQ_MVALUE
                        start = third
                        retcnd = 1
                        state = 0
                    else
                        return "Token following BETWEEN operator is invalid."
                    end
                elseif c2 > 0
                    return "Token sequence must be followed by a value."
                else
                    return "Token must be followed by a comparison operator."
                end
            elseif c1 == _EKQ_MEXPR
                if second <= 0
                    (level > 1 || nmeta > 1) && return "More tokens expected."
                    state = 4
                elseif c2 == _EKQ_MRPAR
                    (level > 1 && popcnd[level] == 2) ||
                        return "Unexpected right parenthesis found."
                    state = 1
                elseif third <= 0
                    return "More tokens expected."
                elseif c3 == _EKQ_MIDENT || c3 == _EKQ_MNAME || c3 == _EKQ_MNOT
                    start = third
                    retcnd = 1
                    state = 0
                elseif c3 == _EKQ_MLPAR
                    start = fourth
                    retcnd = 2
                    state = 0
                elseif c2 == _EKQ_MOR && c3 == _EKQ_MEXPR
                    # OR waits if an AND follows the right expression.
                    if fourth > 0 && mtcode[fourth] == _EKQ_MAND
                        retcnd = 1
                        start = third
                        state = 0
                    else
                        _ek_qlnkilb!(mtexpr[first], mtexpr[third], djpool)
                        _ek_qlnkfsl!(second, third, mtpool)
                        nmeta -= 2
                        state = 3
                    end
                elseif c2 == _EKQ_MAND && c3 == _EKQ_MEXPR
                    # Distribute AND over OR: one new conjunction for each pair of a
                    # conjunction on the left and one on the right, left index slowest.
                    dj1 = mtexpr[first]
                    djtail = 0
                    while dj1 > 0
                        dj2 = mtexpr[third]
                        while dj2 > 0
                            djpool.nfree < 1 && return full_dj
                            newdj = _ek_qlnkan!(djpool)
                            head1 = head2 = 0
                            for side in 1:2
                                cjn = djptrs[side == 1 ? dj1 : dj2]
                                tl = 0
                                while cjn > 0   # invariant: the copy ends at node tl
                                    cjpool.nfree < 1 && return full_cj
                                    newcj = _ek_qlnkan!(cjpool)
                                    _ek_qlnkila!(tl, newcj, cjpool)
                                    tl = newcj
                                    cjptrs[tl] = cjptrs[cjn]
                                    cjn = _ek_qlnknxt(cjn, cjpool)
                                end
                                side == 1 ? (head1 = _ek_qlnkhl(tl, cjpool)) :
                                (head2 = _ek_qlnkhl(tl, cjpool))
                            end
                            _ek_qlnkilb!(head1, head2, cjpool)
                            djptrs[newdj] = head1
                            _ek_qlnkila!(djtail, newdj, djpool)
                            djtail = newdj
                            dj2 = _ek_qlnknxt(dj2, djpool)
                        end
                        dj1 = _ek_qlnknxt(dj1, djpool)
                    end
                    # Free the two old expressions.
                    _ek_qlnkilb!(mtexpr[first], mtexpr[third], djpool)
                    djnode = mtexpr[first]
                    while djnode > 0
                        cjnode = djptrs[djnode]
                        _ek_qlnkfsl!(cjnode, _ek_qlnktl(cjnode, cjpool), cjpool)
                        djnode = _ek_qlnknxt(djnode, djpool)
                    end
                    _ek_qlnkfsl!(mtexpr[first], _ek_qlnktl(mtexpr[first], djpool), djpool)
                    mtexpr[first] = _ek_qlnkhl(newdj, djpool)
                    _ek_qlnkfsl!(second, third, mtpool)
                    nmeta -= 2
                    state = 3
                else
                    return "Unexpected token found following valid expression."
                end
            elseif c1 == _EKQ_MNOT
                if second <= 0
                    return "Tokens were missing from logical expression."
                elseif c2 == _EKQ_MLPAR
                    start = third
                    retcnd = 2
                    state = 0
                elseif c2 == _EKQ_MNOT || c2 == _EKQ_MIDENT || c2 == _EKQ_MNAME
                    start = second
                    retcnd = 1
                    state = 0
                elseif c2 == _EKQ_MEXPR
                    # De Morgan: the negation of a disjunction of conjunctions is the
                    # product over the conjunctions of the disjunction of their negated
                    # relations. Start with the first conjunction, one relation each.
                    djnode = mtexpr[second]
                    cjnode = djptrs[djnode]
                    tl = 0
                    while cjnode > 0
                        djpool.nfree < 1 && return full_dj
                        newdj = _ek_qlnkan!(djpool)
                        cjpool.nfree < 1 && return full_cj
                        newcj = _ek_qlnkan!(cjpool)
                        cjptrs[newcj] = cjptrs[cjnode]
                        djptrs[newdj] = newcj
                        _ek_qlnkila!(tl, newdj, djpool)
                        tl = newdj
                        cjnode = _ek_qlnknxt(cjnode, cjpool)
                    end
                    exprhd = _ek_qlnkhl(tl, djpool)
                    djnode = _ek_qlnknxt(djnode, djpool)
                    while djnode > 0
                        # invariant: exprhd is the product for the conjunctions before djnode
                        dj1 = exprhd
                        djtail = 0
                        while dj1 > 0
                            cj2 = djptrs[djnode]
                            while cj2 > 0
                                cjnode = djptrs[dj1]
                                tl = 0
                                while cjnode > 0
                                    cjpool.nfree < 1 && return full_cj
                                    newcj = _ek_qlnkan!(cjpool)
                                    _ek_qlnkila!(tl, newcj, cjpool)
                                    cjptrs[newcj] = cjptrs[cjnode]
                                    tl = newcj
                                    cjnode = _ek_qlnknxt(cjnode, cjpool)
                                end
                                cj1 = _ek_qlnkhl(tl, cjpool)
                                cjpool.nfree < 1 && return full_cj
                                newcj = _ek_qlnkan!(cjpool)
                                cjptrs[newcj] = cjptrs[cj2]
                                djpool.nfree < 1 && return full_dj
                                newdj = _ek_qlnkan!(djpool)
                                _ek_qlnkilb!(cj1, newcj, cjpool)
                                djptrs[newdj] = cj1
                                _ek_qlnkila!(djtail, newdj, djpool)
                                djtail = newdj
                                cj2 = _ek_qlnknxt(cj2, cjpool)
                            end
                            dj1 = _ek_qlnknxt(dj1, djpool)
                        end
                        dj1 = exprhd
                        while dj1 > 0
                            cj1 = djptrs[dj1]
                            _ek_qlnkfsl!(cj1, _ek_qlnktl(cj1, cjpool), cjpool)
                            dj1 = _ek_qlnknxt(dj1, djpool)
                        end
                        _ek_qlnkfsl!(exprhd, _ek_qlnktl(exprhd, djpool), djpool)
                        exprhd = _ek_qlnkhl(djtail, djpool)
                        djnode = _ek_qlnknxt(djnode, djpool)
                    end
                    # Negate each relation of the old expression once.
                    relset = Set{Int}()
                    djnode = mtexpr[second]
                    while djnode > 0
                        cjnode = djptrs[djnode]
                        while cjnode > 0
                            push!(relset, cjptrs[cjnode])
                            cjnode = _ek_qlnknxt(cjnode, cjpool)
                        end
                        djnode = _ek_qlnknxt(djnode, djpool)
                    end
                    for relptr in relset
                        rels[2, relptr] = _EKQ_CMPNEG[rels[2, relptr]]
                    end
                    mtexpr[first] = exprhd
                    mtcode[first] = _EKQ_MEXPR
                    _ek_qlnkfsl!(second, second, mtpool)
                    nmeta -= 1
                    state = 3
                else
                    return "Token following NOT operator was invalid."
                end
            elseif c1 == _EKQ_MLPAR
                retcnd = 2
                start = second
                state = 0
            else
                return "Unexpected token found."
            end
        elseif state == 3
            if popcnd[level] == 1
                state = 1
            else
                mstart[level] = first
                state = 2
            end
        elseif state == 0
            level += 1
            level > _EKQ_MAXTOK && return "Stack is full"
            mstart[level] = start
            popcnd[level] = retcnd
            state = 2
        else
            level > 1 || return "Syntax error:  badly formed WHERE clause."
            if popcnd[level] == 2
                # Remove the parentheses around the expression.
                prev = _ek_qlnkprv(first, mtpool)
                mstart[level-1] == prev && (mstart[level-1] = first)
                metahd == prev && (metahd = first)
                _ek_qlnkfsl!(prev, prev, mtpool)
                _ek_qlnkfsl!(second, second, mtpool)
                nmeta -= 2
            end
            level -= 1
            state = 2
        end
    end

    # Write the constraints, conjunction after conjunction, and then the sizes.
    # CSPICE counts the conjunctions in an array of 1000 elements while it writes them
    # and goes on after a full cell, so it stops the program at conjunction 1001.
    nconj = 0
    djnode = mtexpr[first]
    while djnode > 0
        nconj += 1
        djnode = _ek_qlnknxt(djnode, djpool)
    end
    nconj > _EKQ_MAXCNJ && _ek_qsignal(:toomanyconjunctions,
        "The normalised WHERE clause has $nconj conjunctions; the maximum is $_EKQ_MAXCNJ.")
    sizes = Int[]
    nrels = 0
    djnode = mtexpr[first]
    while djnode > 0
        push!(sizes, 0)
        cjnode = djptrs[djnode]
        while cjnode > 0
            nrels += 1
            sizes[end] += 1
            relptr = cjptrs[cjnode]
            tabptr = rels[1, relptr]
            op = rels[2, relptr]
            rhsptr = rels[3, relptr]
            desc(n) = ntuple(m -> dsc[m, n], 6)
            _ek_qappnd!(q, 0)
            k = length(q.i)
            _ek_qappnd!(q, desc(tabptr))
            _ek_qappnd!(q, desc(_ek_qlnknxt(tabptr, dspool)))
            if rhsptr == 0
                op = op == _EKQ_EQ ? _EKQ_ISNULL : op == _EKQ_NE ? _EKQ_NOTNUL :
                     return "NULL values can only be used with the operators \"IS " *
                            "NULL\", \"NOT NULL\", or equivalents."
                _ek_qappnd!(q, op)
                for _ in 1:12
                    _ek_qappnd!(q, 0)
                end
                q.i[k] = _EKQ_EQVAL
            else
                _ek_qappnd!(q, op)
                _ek_qappnd!(q, desc(rhsptr))
                if dsc[7, rhsptr] == _EKQ_MVALUE
                    for _ in 1:6
                        _ek_qappnd!(q, 0)
                    end
                    q.i[k] = _EKQ_EQVAL
                else
                    _ek_qappnd!(q, desc(_ek_qlnknxt(rhsptr, dspool)))
                    q.i[k] = _EKQ_EQCOL
                end
            end
            cjnode = _ek_qlnknxt(cjnode, cjpool)
        end
        djnode = _ek_qlnknxt(djnode, djpool)
    end
    q.i[9] = nrels
    q.i[10] = length(sizes)
    for s in sizes
        _ek_qappnd!(q, s)
    end
    return nothing
end

# ---------------------------------------------------------------------------------
# Parser: ZZEKPARS
# ---------------------------------------------------------------------------------

"""
    _ek_qpars!(q, qb, tk) -> Union{Nothing,String}

Contract: `tk` holds the tokens of the query `qb` and `q` is a new encoded query.
Fills `q` with the FROM clause, the normalised WHERE clause, the ORDER BY clause and
the SELECT clause, in this order, and marks it as parsed (`zzekpars_`); returns
`nothing`, or the error message of the parser.

State machine: the states are those of the C source (0 table, 1 after a table, 3 after
an alias, 4 WHERE, 5 ORDER, 6 to 10 order-by columns, 11 SELECT, 12 to 15 SELECT
columns, 16 end). Each pass of the loop takes the next token; the states 4, 5 and 11
jump to the first occurrence of their keyword.
"""
function _ek_qpars!(q::_EKQEnc, qb::AbstractVector{UInt8}, tk::_EKQTokens)
    ntoken = tk.n
    tokens, values = tk.tokens, tk.values
    toknum, fnd = _ek_qtloc(tk, _EKQ_TKKEY, _EKQ_KWFROM)
    fnd || return "Every query must contain a FROM clause. The FROM keyword was not found."
    state = 0
    ntabs = nsel = norder = 0
    errtyp = ""
    namdsc = _EKQ_NODESCR
    inqc(k) = _ek_qinqc!(q, _ek_qtext(tk, k), tk.lxbegs[k], tk.lxends[k])
    # Replace the last six integers of the cell (MOVEI to the end of the cell).
    function setlast!(d)
        n = length(q.i)
        for m in 1:6
            q.i[n-6+m] = d[m]
        end
    end
    endsfrom(v) = v == _EKQ_KWWHER || v == _EKQ_KWSEL || v == _EKQ_KWORDR
    endssel(v) = v == _EKQ_KWWHER || v == _EKQ_KWFROM || v == _EKQ_KWORDR
    endsord(v) = v == _EKQ_KWWHER || v == _EKQ_KWFROM || v == _EKQ_KWSEL
    while state != 16
        toknum += 1
        token = toknum > ntoken ? _EKQ_TKEOQ : tokens[toknum]
        value = toknum > ntoken ? 0 : values[toknum]
        if state == 0                               # a table name is expected
            if token == _EKQ_TKID
                ntabs += 1
                if ntabs > _EKQ_MAXTAB
                    errtyp = "TOO_MANY_TABLES"
                    state = 16
                else
                    _ek_qappnd!(q, inqc(toknum))
                    _ek_qappnd!(q, _EKQ_NODESCR)
                    q.i[8] = ntabs
                    state = 1
                end
            else
                errtyp = token == _EKQ_TKEOQ ? "MORE_TOKENS_EXP" : "TABLE_EXP"
                state = 16
            end
        elseif state == 1 || state == 3             # after a table name or an alias
            if token == _EKQ_TKEOQ
                state = 4
            elseif token == _EKQ_TKCOMA
                state = 0
            elseif token == _EKQ_TKID && state == 1
                setlast!(inqc(toknum))
                state = 3
            elseif token == _EKQ_TKKEY
                if endsfrom(value)
                    state = 4
                else
                    errtyp = "BAD_KEYWORD"
                    state = 16
                end
            else
                errtyp = state == 1 ? "ALIAS_EXP" : "COMMA_OR_KEY_EXP"
                state = 16
            end
        elseif state == 4                           # WHERE clause
            toknum, fnd = _ek_qtloc(tk, _EKQ_TKKEY, _EKQ_KWWHER)
            if fnd
                msg = _ek_qnrml!(q, qb, tk)
                msg === nothing || return msg
            end
            state = 5
        elseif state == 5                           # ORDER BY
            toknum, fnd = _ek_qtloc(tk, _EKQ_TKKEY, _EKQ_KWORDR)
            if !fnd
                state = 11
            elseif toknum < ntoken
                toknum += 1
                if tokens[toknum] == _EKQ_TKKEY && values[toknum] == _EKQ_KWBY
                    state = 6
                else
                    errtyp = "BY_EXPECTED"
                    state = 16
                end
            else
                errtyp = "BY_EXPECTED"
                state = 16
            end
        elseif state == 6 || state == 12            # an order-by or a SELECT column
            if token == _EKQ_TKID
                n = state == 6 ? (norder += 1) : (nsel += 1)
                if n > (state == 6 ? _EKQ_MAXORD : _EKQ_MAXSEL)
                    errtyp = state == 6 ? "TOO_MANY_ORD_COLS" : "TOO_MANY_SEL_COLS"
                    state = 16
                else
                    namdsc = inqc(toknum)
                    _ek_qappnd!(q, _EKQ_NODESCR)
                    q.i[state == 6 ? 11 : 12] = n
                    state += 2
                end
            else
                errtyp = token == _EKQ_TKEOQ ? "MORE_TOKENS_EXP" : "TABLE_OR_COLUMN_EXP"
                state = 16
            end
        elseif state == 8                           # after the first name, order-by
            if token == _EKQ_TKEOQ
                _ek_qappnd!(q, namdsc)
                _ek_qappnd!(q, _EKQ_ASCEND)
                state = 11
            elseif token == _EKQ_TKCOMA
                _ek_qappnd!(q, namdsc)
                _ek_qappnd!(q, _EKQ_ASCEND)
                state = 6
            elseif token == _EKQ_TKDOT
                state = 7
            elseif token == _EKQ_TKKEY
                _ek_qappnd!(q, namdsc)
                if value == _EKQ_KWASC
                    _ek_qappnd!(q, _EKQ_ASCEND)
                    state = 10
                elseif value == _EKQ_KWDESC
                    _ek_qappnd!(q, _EKQ_DESCEND)
                    state = 10
                elseif endsord(value)
                    _ek_qappnd!(q, _EKQ_ASCEND)
                    state = 11
                else
                    errtyp = "BAD_KEYWORD"
                    state = 16
                end
            else
                errtyp = "BAD_TOKEN"
                state = 16
            end
        elseif state == 7 || state == 13            # the column after `table.`
            if token == _EKQ_TKID
                setlast!(namdsc)
                _ek_qappnd!(q, inqc(toknum))
                state += 2
            else
                errtyp = token == _EKQ_TKEOQ ? "MORE_TOKENS_EXP" : "COLUMN_EXP"
                state = 16
            end
        elseif state == 9                           # after `table.column`, order-by
            if token == _EKQ_TKKEY
                if value == _EKQ_KWASC
                    _ek_qappnd!(q, _EKQ_ASCEND)
                    state = 10
                elseif value == _EKQ_KWDESC
                    _ek_qappnd!(q, _EKQ_DESCEND)
                    state = 10
                elseif endsord(value)
                    # CSPICE appends no sense here: the encoded query is one integer
                    # short and the SELECT columns that follow are shifted.
                    state = 11
                else
                    errtyp = "BAD_KEYWORD"
                    state = 16
                end
            elseif token == _EKQ_TKCOMA
                _ek_qappnd!(q, _EKQ_ASCEND)
                state = 6
            elseif token == _EKQ_TKEOQ
                _ek_qappnd!(q, _EKQ_ASCEND)
                state = 11
            else
                errtyp = "COMMA_OR_KEY_EXP"
                state = 16
            end
        elseif state == 10                          # after ASC or DESC
            if token == _EKQ_TKCOMA
                state = 6
            elseif token == _EKQ_TKEOQ
                state = 11
            elseif token == _EKQ_TKKEY
                if endsord(value)
                    state = 11
                else
                    errtyp = "BAD_KEYWORD"
                    state = 16
                end
            else
                errtyp = "COMMA_EXP"
                state = 16
            end
        elseif state == 11                          # SELECT clause
            toknum, fnd = _ek_qtloc(tk, _EKQ_TKKEY, _EKQ_KWSEL)
            if fnd
                state = 12
            else
                errtyp = "SELECT_NOT_FOUND"
                state = 16
            end
        elseif state == 14                          # after the first name, SELECT
            if token == _EKQ_TKEOQ
                _ek_qappnd!(q, namdsc)
                state = 16
            elseif token == _EKQ_TKCOMA
                _ek_qappnd!(q, namdsc)
                state = 12
            elseif token == _EKQ_TKDOT
                state = 13
            elseif token == _EKQ_TKKEY
                _ek_qappnd!(q, namdsc)
                endssel(value) || (errtyp = "BAD_KEYWORD")
                state = 16
            else
                errtyp = "BAD_TOKEN"
                state = 16
            end
        else                                        # state 15: after `table.column`
            if token == _EKQ_TKKEY
                endssel(value) || (errtyp = "BAD_KEYWORD")
                state = 16
            elseif token == _EKQ_TKCOMA
                state = 12
            elseif token == _EKQ_TKEOQ
                state = 16
            else
                errtyp = "COMMA_OR_KEY_EXP"
                state = 16
            end
        end
    end
    if errtyp == ""
        q.i[4] = 1
        return nothing
    end
    errtyp == "SELECT_NOT_FOUND" && return "Every query must contain a SELECT clause. " *
                                           "The SELECT keyword was not found."
    errtyp == "BY_EXPECTED" &&
        return "The BY keyword was not found following the ORDER keyword."
    errtyp == "MORE_TOKENS_EXP" && return "More tokens were expected in query."
    errtyp == "TOO_MANY_TABLES" && return _ek_qrep("Number of tables in \"FROM\" " *
        "clause exceeds allowed maximum of #.", _EKQ_MAXTAB)
    errtyp == "TOO_MANY_ORD_COLS" && return _ek_qrep("Number of order-by columns " *
        "exceeds allowed maximum of #.", _EKQ_MAXORD)
    errtyp == "TOO_MANY_SEL_COLS" && return _ek_qrep("Number of SELECT columns " *
        "exceeds allowed maximum of #.", _EKQ_MAXSEL)
    template = if errtyp == "BAD_KEYWORD"
        "Invalid keyword at location #. Actual token was: #"
    elseif errtyp == "TABLE_OR_COLUMN_EXP"
        "Table or column name expected at location #. Actual token was: #"
    elseif errtyp == "TABLE_EXP"
        "Table name expected at location #. Actual token was: #"
    elseif errtyp == "COLUMN_EXP"
        "Column name expected at location #. Actual token was: #"
    elseif errtyp == "ALIAS_EXP"
        "Table alias, comma, or keyword expected at location #. Actual token was: #"
    elseif errtyp == "COMMA_OR_KEY_EXP"
        "Comma or keyword expected at location #. Actual token was: #"
    elseif errtyp == "COMMA_EXP"
        "Comma expected at location #. Actual token was: #"
    else
        "Invalid token at location #. Token was: #"
    end
    return _ek_qlexmsg(template, qb, tk.lxbegs[toknum], tk.lxends[toknum], true)
end

# ---------------------------------------------------------------------------------
# Name resolution: ZZEKNRES, ZZEKCCHK
# ---------------------------------------------------------------------------------

"Contract: `1 <= n <= number of tables of q`; returns the name and the alias (`\"\"` if none) of table `n` of the FROM clause, cut to 64 characters (`zzekqtab_`)."
function _ek_qtab(q::_EKQEnc, n::Int)
    base = (n - 1) * 12 + 19
    tb, te = _ek_qget(q, base + 9), _ek_qget(q, base + 10)
    ab, ae = _ek_qget(q, base + 15), _ek_qget(q, base + 16)
    table = _ek_qstr(q, tb, min(te, tb + _EKQ_TNAMSZ - 1))
    alias = ab > 0 ? _ek_qstr(q, ab, min(ae, ab + _EKQ_TNAMSZ - 1)) : ""
    return table, alias
end

"""
    _ek_qcchk!(q, tables, tabcat, table, alias, base) -> Union{Nothing,String}

Contract: `base` is the C base of a pair of descriptors (qualifier, column) in `q`;
`table`, `alias` are the names of the FROM clause and `tabcat` their catalogue
indices. Writes the table index and the column index into the pair and returns
`nothing`, or returns the error message of `zzekcchk_`. A qualifier is looked up in
the table names first and then in the aliases; a column without a qualifier must be
in exactly one table of the FROM clause.
"""
function _ek_qcchk!(q::_EKQEnc, tables::Vector{EKTable}, tabcat::Vector{Int},
    table::Vector{String}, alias::Vector{String}, base::Int)
    tb, te = _ek_qget(q, base + 9), _ek_qget(q, base + 10)
    cb, ce = _ek_qget(q, base + 15), _ek_qget(q, base + 16)
    (cb <= 0 || ce <= 0) && _ek_qsignal(:bug, "Invalid string bounds $cb:$ce for " *
        "column.  Column name descriptor base is $base.")
    column = _ek_qstr(q, cb, ce)
    tabidx = colidx = 0
    if tb == 0
        nmatch = 0
        for i in eachindex(table), (j, name) in enumerate(tables[tabcat[i]].cnames)
            if name == column
                nmatch += 1
                colidx = j
                tabidx = i
            end
        end
        nmatch == 0 && return _ek_qrep("Column <#> is not present in any table in " *
            "FROM clause.", column)
        nmatch > 1 && return _ek_qrep("Column name <#> is ambiguous; a qualifying " *
            "table name or alias is required.", column)
    else
        qual = _ek_qstr(q, tb, te)
        tabidx = something(findfirst(==(qual), table), findfirst(==(qual), alias), 0)
        tabidx == 0 &&
            return _ek_qrep("Table name <#> is not present in FROM clause.", qual)
        colidx = something(findfirst(==(column), tables[tabcat[tabidx]].cnames), 0)
        colidx == 0 && return _ek_qrep(_ek_qrep("Column <#> does not exist in table " *
            "<#>.", column), qual)
    end
    _ek_qput!(q, base + 11, tabidx)
    _ek_qput!(q, base + 17, colidx)
    return nothing
end

"""
    _ek_qnres!(q, qb, tables) -> (tabcat, error)

Contract: `q` is a parsed query of the string `qb` and `tables` the catalogue of
loaded tables. Resolves the table and column names of `q` (`zzeknres_`): returns the
catalogue index of each table of the FROM clause and `nothing`, or an empty vector and
the error message. Rules, in the order of CSPICE: aliases are distinct; an alias is not
the name of another table of the FROM clause; each table is loaded; then the columns
of the constraints, of the SELECT clause and of the ORDER BY clause are resolved.
"""
function _ek_qnres!(q::_EKQEnc, qb::AbstractVector{UInt8}, tables::Vector{EKTable})
    q.i[4] == -1 && _ek_qsignal(:querynotparsed, "Encoded query has not been parsed.")
    ntab, ncns, ncnj, nord, nsel = (Int(q.i[k]) for k in (8, 9, 10, 11, 12))
    table = String[]
    alias = String[]
    for i in 1:ntab
        t, a = _ek_qtab(q, i)
        push!(table, t)
        push!(alias, a)
    end
    none = Int[]
    lexeme(base) = _ek_qlex(qb, _ek_qget(q, base + 7), _ek_qget(q, base + 8))
    for i in 1:ntab-1, j in i+1:ntab
        if alias[i] == alias[j] && alias[i] != ""
            return none, _ek_qrep("Non-distinct alias <#> was found.",
                lexeme(((j - 1) * 2 + 1) * 6 + 19))
        end
    end
    for i in 1:ntab
        j = alias[i] == "" ? nothing : findfirst(==(alias[i]), table)
        if j !== nothing && j != i
            return none, _ek_qrep("Alias <#> conflicts with table name.",
                lexeme(((i - 1) * 2 + 1) * 6 + 19))
        end
    end
    tabcat = Int[]
    for i in 1:ntab
        j = findfirst(t -> t.name == table[i], tables)
        j === nothing && return none, _ek_qrep("Table <#> is not currently loaded.",
            lexeme((i - 1) * 12 + 19))
        push!(tabcat, j)
    end
    for i in 1:ntab
        base = (i - 1) * 12 + 19
        _ek_qput!(q, base + 11, i)
        _ek_qput!(q, base + 17, i)
    end
    for i in 1:ncns
        base = ntab * 12 + 19 + (i - 1) * 26
        msg = _ek_qcchk!(q, tables, tabcat, table, alias, base + 1)
        msg === nothing || return none, msg
        if _ek_qget(q, base + 6) == _EKQ_EQCOL
            msg = _ek_qcchk!(q, tables, tabcat, table, alias, base + 14)
            msg === nothing || return none, msg
        end
    end
    for i in 1:nsel
        base = ntab * 12 + 19 + ncnj + ncns * 26 + nord * 13 + (i - 1) * 12
        msg = _ek_qcchk!(q, tables, tabcat, table, alias, base)
        msg === nothing || return none, msg
    end
    for i in 1:nord
        base = ntab * 12 + 19 + ncnj + ncns * 26 + (i - 1) * 13
        msg = _ek_qcchk!(q, tables, tabcat, table, alias, base)
        msg === nothing || return none, msg
    end
    q.i[5] = 1
    return tabcat, nothing
end

# ---------------------------------------------------------------------------------
# Time resolution: ZZEKTRES, ZZEKTCNV with SCN2ID, SCPARS, SCPS01
# ---------------------------------------------------------------------------------

"""
    _ek_scn2id(ks, clknam) -> Union{Int,Nothing}

Contract: `clknam` is the name of a spacecraft clock, `"<body name> SCLK"`. Returns
the ID code of the clock, which is the code of the body, or `nothing` if the body name
is unknown or missing (`scn2id_`): the last `SCLK` of the name is removed and the rest
is a body name. Reads the body names of `ks` only.
"""
function _ek_scn2id(ks::KernelSet, clknam::AbstractString)
    name = uppercase(String(first(clknam, 32)))
    r = findlast("SCLK", name)
    r === nothing || (name = string(SubString(name, 1, prevind(name, first(r))), "    ",
        SubString(name, last(r) + 1)))
    name = strip(name, ' ')
    isempty(name) && return nothing
    return bodn2c(ks, name)
end

"""
    _ek_qsclkfields(clkstr) -> Vector{String}

Contract: returns the fields of the clock string `clkstr` as `lparsm_` gives them with
the clock delimiters. `_lparsm` takes byte positions and so accepts ASCII text only:
each other byte is replaced by a placeholder for the split and put back afterwards, in
order, so that the fields hold the bytes of `clkstr` unchanged.
"""
function _ek_qsclkfields(clkstr::AbstractString)
    isascii(clkstr) && return String.(_lparsm(clkstr, _SCLK_DELIMS, _SCLK_MXNFLD))
    bytes = Vector{UInt8}(codeunits(String(clkstr)))
    held = filter(b -> b >= 0x80 || b == 0x01, bytes)
    safe = String(map(b -> b >= 0x80 ? 0x01 : b, bytes))
    k = 0
    fields = String[]
    for item in _lparsm(safe, _SCLK_DELIMS, _SCLK_MXNFLD)
        raw = Vector{UInt8}(codeunits(String(item)))
        for j in eachindex(raw)             # invariant: k placeholders are put back
            raw[j] == 0x01 && (k += 1; raw[j] = held[k])
        end
        push!(fields, String(raw))
    end
    return fields
end

"""
    _ek_qscps01(clock, clkstr) -> (ticks, error)

Contract: `clkstr` is a clock string without partition number. Returns its ticks and
`""`, or 0 and the error message of `scps01_`.
"""
function _ek_qscps01(clock::SclkType1, clkstr::AbstractString)
    isempty(strip(clkstr, ' ')) &&
        return 0.0, "Non partition part of the input clock string is blank."
    cmptks = _sclk_cmptks(clock)
    items = _ek_qsclkfields(clkstr)
    n = length(items)
    n > clock.nfield && return 0.0, _ek_qrep(_ek_qrep(_ek_qrep("Input clock string # " *
        "has # fields; maximum for this spacecraft clock is #.", clkstr), n),
        Int(clock.nfield))
    ticks = 0.0
    for (i, item) in enumerate(items)       # invariant: ticks = sum over the fields 1:i-1
        value = clock.offsets[i]
        if !isempty(strip(item, ' '))
            try
                isascii(item) || throw(SpiceError(:notanumber, "not a number"))
                value = _nparsd(strip(item, ' '))
            catch err
                err isa SpiceError || rethrow()
                return 0.0, _ek_qrep(_ek_qrep("Could not parse SCLK component # from " *
                    "# as a number.", item), clkstr)
            end
        end
        value -= clock.offsets[i]
        _dnint(value) < 0 && return 0.0, _ek_qrep(_ek_qrep(_ek_qrep("Component number " *
            "#, # in the SCLK string  # is invalid.", i), item), clkstr)
        ticks += value * cmptks[i]
    end
    return ticks, ""
end

"""
    _ek_qscpars(ks, sc, sclkch) -> (sclkdp, error)

Contract: returns the encoded clock value of the clock string `sclkch`
(`"partition/fields"`, the partition optional) of the clock `sc` and `""`, or 0 and the
error message of `scpars_`. Throws `SpiceError` if the clock `sc` is not loaded.
"""
function _ek_qscpars(ks::KernelSet, sc::Int, sclkch::AbstractString)
    clock = _sclk(ks, sc)
    pstart, pstop, ptotls = _sclk_partitions(clock)
    nparts = length(pstart)
    isempty(strip(sclkch, ' ')) && return 0.0, "Input spacecraft clock string is blank."
    pos = something(findfirst('/', sclkch), 0)
    badpart = "Unable to parse the partition number from SCLK string #."
    pos == 1 && return 0.0, _ek_qrep(badpart, sclkch)
    fields = sclkch[pos+1:end]
    ticks, pserr = _ek_qscps01(clock, fields)
    isempty(pserr) || return 0.0, pserr
    part = 1
    if pos > 1
        x = try
            _nparsd(SubString(sclkch, 1, prevind(sclkch, pos)))
        catch err
            err isa SpiceError || rethrow()
            return 0.0, _ek_qrep(badpart, sclkch)
        end
        abs(x) < 2.0^31 || return 0.0, _ek_qrep(badpart, sclkch)
        part = Int(_dnint(x))
        (part <= 0 || part > nparts) && return 0.0, _ek_qrep(_ek_qrep(_ek_qrep(
            "Partition number # taken from SCLK string # is not in acceptable range " *
            "1 to #.", part), sclkch), nparts)
        (ticks < pstart[part] || ticks > pstop[part]) && return 0.0, _ek_qrep(_ek_qrep(
            "SCLK count from # does not fall in the boundaries of partition number #.",
            fields), part)
    else
        while part <= nparts && (ticks < pstart[part] || ticks > pstop[part])
            part += 1                # invariant: partitions 1:part-1 do not hold `ticks`
        end
        part > nparts && return 0.0, _ek_qrep(_ek_qrep("SCLK count # does not fall in " *
            "the boundaries of any of the partitions for spacecraft #.", fields), sc)
    end
    return part > 1 ? ticks - pstart[part] + ptotls[part-1] : ticks - pstart[part], ""
end

"""
    _ek_qtcnv(ks, td, timstr) -> (et, error)

Contract: `timstr` holds the 32 bytes of a time value of a query. Returns its
ephemeris time and `""`, or 0 and the error message of `zzektcnv_`. The string is a
spacecraft clock string if it contains `SCLK` (`<clock name> SCLK <clock string>`),
and a time string of `str2et` with the defaults `td` in the other case. An error of
the time or clock system (no leap-second kernel, no SCLK kernel, a time that
`str2et` rejects after `tpartv` accepted it) gives the message "Unexpected SPICELIB
error ..."; CSPICE leaves its SPICE error signalled in this case as well.
"""
function _ek_qtcnv(ks::KernelSet, td::TimeDefaults, timstr::AbstractVector{UInt8})
    # CMPRSS (runs of blanks become one blank), LJUST and UCASE.
    loc = UInt8[]
    for c in timstr
        (c == UInt8(' ') && (isempty(loc) || loc[end] == UInt8(' '))) && continue
        push!(loc, _ek_qupper(c))
    end
    (!isempty(loc) && loc[end] == UInt8(' ')) && pop!(loc)
    locstr = String(loc)
    orig = String(copy(timstr))
    unexpected() = string("Unexpected SPICELIB error encountered while attempting to ",
        "parse the string <", rstrip(orig, ' '), ">")
    r = findlast("SCLK", locstr)
    if r !== nothing
        pos = first(r)
        clkid = _ek_scn2id(ks, locstr[1:pos+3])
        if clkid === nothing
            pos > 1 && return 0.0, _ek_qrep("Time conversion failed; SCLK type <#> " *
                "was not recognized.", String(timstr[1:pos-1]))
            return 0.0, "Time conversion failed; SCLK name was not supplied."
        end
        try
            sclkdp, msg = _ek_qscpars(ks, clkid, locstr[pos+4:end])
            isempty(msg) || return 0.0, string(_ek_qrep("The string <#> didn't parse " *
                "as a spacecraft clock string.", orig), "   ", msg)
            return Float64(sct2e(ks, clkid, sclkdp)), ""
        catch err
            err isa SpiceError || rethrow()
            return 0.0, unexpected()
        end
    end
    p = _tpartv(locstr, 32)
    p.ok || return 0.0, p.error
    try
        return str2et(ks, locstr, td), ""
    catch err
        err isa SpiceError || rethrow()
        return 0.0, unexpected()
    end
end

"""
    _ek_qtres!(q, tables, tabcat, ks, td) -> Union{Nothing,String}

Contract: the names of `q` are resolved. For each constraint that compares a column
of type TIME with a value, converts the value to ephemeris time, stores it in the
double buffer and replaces the value descriptor by one of type 4 (`zzektres_`).
Returns `nothing`, or the error message of the conversion. As in CSPICE the time
string is cut to 32 characters, and a numeric value is converted as a blank string,
which is an error.
"""
function _ek_qtres!(q::_EKQEnc, tables::Vector{EKTable}, tabcat::Vector{Int},
    ks::KernelSet, td::TimeDefaults)
    q.i[5] == -1 &&
        _ek_qsignal(:namesnotresolved, "Names are not resolved in encoded query.")
    ntab, ncns = Int(q.i[8]), Int(q.i[9])
    for i in 1:ncns
        base = ntab * 12 + 19 + (i - 1) * 26
        opcode = _ek_qget(q, base + 19)
        (_ek_qget(q, base + 6) == _EKQ_EQVAL && opcode != _EKQ_ISNULL &&
         opcode != _EKQ_NOTNUL) || continue
        attr = tables[tabcat[_ek_qget(q, base + 12)]].attrs[_ek_qget(q, base + 18)]
        attr.dtype === :TIME || continue
        lxb = _ek_qget(q, base + 21)
        sb, se = _ek_qget(q, base + 23), _ek_qget(q, base + 24)
        timstr = fill(UInt8(' '), _EKQ_TIMLEN)
        se >= sb && copyto!(timstr, 1, q.c, sb, min(se - sb + 1, _EKQ_TIMLEN))
        et, msg = _ek_qtcnv(ks, td, timstr)
        isempty(msg) || return msg
        descr = _ek_qinqn!(q, et, _EK_TIME, lxb, lxb)
        for m in 1:6
            _ek_qput!(q, base + 19 + m, descr[m])
        end
    end
    q.i[6] = 1
    return nothing
end

# ---------------------------------------------------------------------------------
# Semantic check: ZZEKSEMC
# ---------------------------------------------------------------------------------

"""
    _ek_qsemc!(q, qb, tables, tabcat) -> Union{Nothing,String}

Contract: the names and the times of `q` are resolved. Checks the constraints and the
order-by columns (`zzeksemc_`) and returns `nothing`, or the error message. Rules: a
column of a constraint or of the ORDER BY clause is scalar (size 1); LIKE and NOT LIKE
need a character column on the left; the two sides of a constraint have the same
type, or both are numeric (DP and INT mix). As in CSPICE, the check stops without an
error at the first constraint with IS NULL or IS NOT NULL: the constraints after it
and the order-by columns are not checked.

In the message for a non-scalar column CSPICE prints the query from the start of the
left column name to the end of the right side of the most recent constraint that had
one; before such a constraint that end is not defined in CSPICE, and this port prints
the column name alone.
"""
function _ek_qsemc!(q::_EKQEnc, qb::AbstractVector{UInt8}, tables::Vector{EKTable},
    tabcat::Vector{Int})
    q.i[5] == -1 &&
        _ek_qsignal(:unresolvednames, "Encoded query has not had names resolved.")
    q.i[6] == -1 &&
        _ek_qsignal(:unresolvedtimes, "Encoded query has not had time values resolved.")
    ntab, ncns, ncnj, nord = (Int(q.i[k]) for k in (8, 9, 10, 11))
    attr(tabidx, colidx) = tables[tabcat[tabidx]].attrs[colidx]
    lxe1 = 0            # LXE(2) of ZZEKSEMC: 0 while CSPICE has not set it
    nonscalar(what, lxb, lxe, size) = _ek_qrep(_ek_qrep("Non-scalar column <#> having " *
        "size # found in $what.", _ek_qlex(qb, lxb, lxe)), size)
    mismatch = "Data type mismatch: column <#> has data type #; # <#> has data type #."
    for i in 1:ncns
        base = ntab * 12 + 19 + (i - 1) * 26
        cnstyp = _ek_qget(q, base + 6)
        lxb0, lxe0 = _ek_qget(q, base + 14), _ek_qget(q, base + 15)
        lhs = attr(_ek_qget(q, base + 12), _ek_qget(q, base + 18))
        lhstyp = _ek_dtype_code(lhs.dtype)
        lhs.size != 1 && return nonscalar("query constraint", lxb0,
            lxe1 == 0 ? lxe0 : lxe1, lhs.size)
        opcode = _ek_qget(q, base + 19)
        if (opcode == _EKQ_LIKE || opcode == _EKQ_UNLIKE) && lhstyp != _EK_CHR
            return _ek_qrep(_ek_qrep("LIKE and NOT LIKE operators may be used only " *
                "with character columns.  Column <#> has type #.",
                _ek_qlex(qb, lxb0, lxe0)), _EKQ_TYPSTR[lhstyp])
        end
        rhstyp = 0
        what = "column"
        if cnstyp == _EKQ_EQCOL
            lxb1, lxe1 = _ek_qget(q, base + 27), _ek_qget(q, base + 28)
            rhs = attr(_ek_qget(q, base + 25), _ek_qget(q, base + 31))
            rhstyp = _ek_dtype_code(rhs.dtype)
            rhs.size != 1 && return nonscalar("query constraint", lxb0, lxe1, rhs.size)
        else
            (opcode == _EKQ_ISNULL || opcode == _EKQ_NOTNUL) && return nothing
            rhstyp = _ek_qget(q, base + 20)
            lxb1, lxe1 = _ek_qget(q, base + 21), _ek_qget(q, base + 22)
            what = "value"
        end
        if rhstyp != lhstyp && (lhstyp == _EK_TIME || lhstyp == _EK_CHR ||
                                rhstyp == _EK_TIME || rhstyp == _EK_CHR)
            msg = _ek_qrep(mismatch, _ek_qlex(qb, lxb0, lxe0))
            msg = _ek_qrep(msg, _EKQ_TYPSTR[lhstyp])
            msg = _ek_qrep(msg, what)
            msg = _ek_qrep(msg, _ek_qlex(qb, lxb1, lxe1))
            return _ek_qrep(msg, _EKQ_TYPSTR[rhstyp])
        end
    end
    for i in 1:nord
        base = ntab * 12 + 19 + ncnj + ncns * 26 + (i - 1) * 13
        a = attr(_ek_qget(q, base + 11), _ek_qget(q, base + 17))
        lxb0, lxe0 = _ek_qget(q, base + 13), _ek_qget(q, base + 14)
        a.size != 1 &&
            return nonscalar("order-by column", lxb0, lxe1 == 0 ? lxe0 : lxe1, a.size)
    end
    q.i[7] = 1
    return nothing
end

# ---------------------------------------------------------------------------------
# The pipeline and the decoded query: ZZEKENCD, ZZEKQCON, ZZEKQORD, ZZEKQSEL
# ---------------------------------------------------------------------------------

"""
    _ek_qencode(tables, ks, td, query) -> (q, tabcat, error, errmsg)

Contract: as `_ek_encode_query`, but returns the encoded query in the layout of CSPICE
(`zzekencd_`) with the catalogue indices of its tables; `q` is `nothing` on error.
"""
function _ek_qencode(tables::Vector{EKTable}, ks::KernelSet, td::TimeDefaults,
    query::AbstractString)
    qb = codeunits(String(query))
    l = something(findlast(!=(UInt8(' ')), qb), 0)
    none = Int[]
    l > _EKQ_MAXQRY && return nothing, none, true, _ek_qrep(_ek_qrep("Query has # " *
        "characters; the maximum allowed is #.", l), _EKQ_MAXQRY)
    tk, msg = _ek_qscan(qb, _EKQ_MAXTOK, _EKQ_MAXNUM, _EKQ_CSIZE)
    tk === nothing && return nothing, none, true, msg
    q = _ek_qini()
    msg = _ek_qpars!(q, qb, tk)
    msg === nothing || return nothing, none, true, msg
    tabcat, msg = _ek_qnres!(q, qb, tables)
    msg === nothing || return nothing, none, true, msg
    msg = _ek_qtres!(q, tables, tabcat, ks, td)
    msg === nothing || return nothing, none, true, msg
    msg = _ek_qsemc!(q, qb, tables, tabcat)
    msg === nothing || return nothing, none, true, msg
    return q, tabcat, false, ""
end

"Contract: returns the integer nearest to `x` as `i_dnnt` of f2c gives it on x86: the smallest 32-bit integer if the result does not fit."
function _ek_qnint(x::Float64)
    r = x >= 0.0 ? floor(x + 0.5) : -floor(0.5 - x)
    return (-2.0^31 <= r < 2.0^31) ? Int(r) : Int(typemin(Int32))
end

"Contract: `base` is the C base of a pair of descriptors (qualifier, column) of a resolved query; returns the column reference, with the lexeme bounds by the rule of `zzekqsel_`."
function _ek_qcol(q::_EKQEnc, base::Int)
    tb, te = _ek_qget(q, base + 9), _ek_qget(q, base + 10)
    cb, ce = _ek_qget(q, base + 15), _ek_qget(q, base + 16)
    return _EKQCol(tb > 0 ? _ek_qstr(q, tb, te) : "", _ek_qget(q, base + 11),
        _ek_qstr(q, cb, ce), _ek_qget(q, base + 17),
        _ek_qget(q, tb > 0 ? base + 7 : base + 13), _ek_qget(q, base + 14))
end

"The unused right side of a constraint that compares a column with a value."
const _EKQ_NOCOL = _EKQCol("", 0, "", 0, 0, 0)

"""
    _ek_qdecode(q, tabcat) -> _EKQuery

Contract: `q` is an encoded query that passed all steps and `tabcat` holds the
catalogue indices of its tables. Returns the query as plain data, with the values that
`zzekqtab_`, `zzekqcnj_`, `zzekqcon_`, `zzekqord_` and `zzekqsel_` return.
"""
function _ek_qdecode(q::_EKQEnc, tabcat::Vector{Int})
    ntab, ncns, ncnj, nord, nsel = (Int(q.i[k]) for k in (8, 9, 10, 11, 12))
    names = [_ek_qtab(q, n) for n in 1:ntab]
    cbase = ntab * 12 + 19
    constraints = _EKQConstraint[]
    for n in 1:ncns
        base = cbase + (n - 1) * 26
        cnstyp = _ek_qget(q, base + 6)
        opcode = _ek_qget(q, base + 19)
        lhs = _ek_qcol(q, base + 1)
        if cnstyp == _EKQ_EQCOL
            push!(constraints, _EKQConstraint(cnstyp, lhs, opcode, _ek_qcol(q, base + 14),
                0, "", 0.0, 0))
        elseif opcode == _EKQ_ISNULL || opcode == _EKQ_NOTNUL
            push!(constraints, _EKQConstraint(cnstyp, lhs, opcode, _EKQ_NOCOL, 0, "", 0.0, 0))
        else
            dtype = _ek_qget(q, base + 20)
            p = _ek_qget(q, base + 23)
            cval = dtype == _EK_CHR ? _ek_qstr(q, p, _ek_qget(q, base + 24)) : ""
            dval = dtype == _EK_CHR ? 0.0 : q.d[p]
            ival = dtype == _EK_INT ? _ek_qnint(dval) : 0
            push!(constraints,
                _EKQConstraint(cnstyp, lhs, opcode, _EKQ_NOCOL, dtype, cval, dval, ival,
                    dtype == _EK_CHR ? 0.0 : q.dx[p]))
        end
    end
    conjsizes = [Int(q.i[cbase+ncns*26+n]) for n in 1:ncnj]
    obase = cbase + ncnj + ncns * 26
    orderby = [(_ek_qcol(q, obase + (n - 1) * 13), _ek_qget(q, obase + (n - 1) * 13 + 18))
               for n in 1:nord]
    sbase = obase + nord * 13
    select = [_ek_qcol(q, sbase + (n - 1) * 12) for n in 1:nsel]
    return _EKQuery(first.(names), last.(names), copy(tabcat), conjsizes, constraints,
        orderby, select)
end

"""
    _ek_encode_query(tables, ks, td, query) -> (query, error, errmsg)

Contract: `tables` is the catalogue of loaded tables in load order, with names in
upper case. Runs the steps of `ekfind_` that come before the search (`zzekqini_`,
`zzekscan_`, `zzekpars_`, `zzeknres_`, `zzektres_`, `zzeksemc_`) on the string `query`
and returns `(q::_EKQuery, false, "")`, or `(nothing, true, errmsg)` with the error
message that `ekfind_c` returns in `errmsg` for a query that is not valid. Positions
in a message are 1-based byte positions in `query`. Time values are converted with
the leap-second and clock data of `ks` and the time defaults `td`. Pure function of
its arguments; does not modify them.

Throws `SpiceError` where CSPICE signals a SPICE error: `:celltoosmall` if the
normalised WHERE clause has more constraints than the encoded query can hold (1000
without other clauses), `:toomanyconjunctions` for more than 1000 conjunctions (CSPICE
stops the program there), `:buffertoosmall` if the converted times do not fit in the
buffer of 100 numbers, and `:bug` for a query in which a qualified order-by column
comes directly before FROM, WHERE or SELECT (CSPICE then writes an encoded query that
its own name resolution rejects).

Differences from CSPICE:
- A query of more than 2000 characters gives an error message of this port; CSPICE
  reads past its buffer.
- CSPICE reads one token past the end of a query that ends with NOT; this port takes
  that token as absent.
- The message for a non-scalar column, see `_ek_qsemc!`.
- An error of the time system during time resolution is reported through the error
  flag only, see `_ek_qtcnv`.
- `ekfind_c` signals `SPICE(EMPTYSTRING)` for a string of length zero; this function
  reports the missing FROM clause, as `ekfind_` does for a blank string.
- A message is not cut to the length of an output buffer.
"""
function _ek_encode_query(tables::Vector{EKTable}, ks::KernelSet, td::TimeDefaults,
    query::AbstractString)::Tuple{Union{Nothing,_EKQuery},Bool,String}
    q, tabcat, error, errmsg = _ek_qencode(tables, ks, td, query)
    error && return nothing, true, errmsg
    return _ek_qdecode(q, tabcat), false, ""
end

"""
    _ek_psel(tables, ks, td, query) -> NamedTuple

Contract: as `_ek_encode_query` for the arguments. Analyses the SELECT clause of
`query` as `ekpsel_c` and returns `(n, xbegs, xends, xtypes, xclass, tabs, cols, error,
errmsg)`: the number `n` of SELECT items and for each item the first and last
character of its expression in `query` (`xbegs`, `xends`; 1-based byte positions, one
more than the 0-based positions of `ekpsel_c`), its data type (`xtypes`: `:CHR`,
`:DP`, `:INT` or `:TIME`), its class (`xclass`: always `:COL`, a column), the name of
its table (`tabs`: the table name, never the alias) and of its column (`cols`). For a
query that is not valid, `error` is true, `errmsg` holds the message, `n` is 0 and the
vectors are empty.
"""
function _ek_psel(tables::Vector{EKTable}, ks::KernelSet, td::TimeDefaults,
    query::AbstractString)
    q, error, errmsg = _ek_encode_query(tables, ks, td, query)
    error && return (n=0, xbegs=Int[], xends=Int[], xtypes=Symbol[], xclass=Symbol[],
        tabs=String[], cols=String[], error=true, errmsg=errmsg)
    sel = q.select
    return (n=length(sel), xbegs=[c.lxbeg for c in sel], xends=[c.lxend for c in sel],
        xtypes=[tables[q.tabcat[c.tabidx]].attrs[c.colidx].dtype for c in sel],
        xclass=fill(:COL, length(sel)), tabs=[q.tables[c.tabidx] for c in sel],
        cols=[c.name for c in sel], error=false, errmsg="")
end

"""
    _ek_qpcol(tables, q, qcol) -> (tabidx, table, alias, colidx, column, error, errmsg)

Contract: `q` is a query that `_ek_encode_query` returned for the catalogue `tables`.
Parses the column name `qcol` (`<column>` or `<table>.<column>`, where the table is a
name or an alias of the FROM clause of `q`) as `zzekpcol_` and returns the position
and the name of its table in the FROM clause, the alias, and the position and the name
of the column; or `error == true` and the message.
"""
function _ek_qpcol(tables::Vector{EKTable}, q::_EKQuery, qcol::AbstractString)
    fail(msg) = (0, "", "", 0, "", true, msg)
    qb = codeunits(String(qcol))
    tk, msg = _ek_qscan(qb, 3, 0, 160)
    tk === nothing && return fail(msg)
    badcol = "Invalid column name; name should consist of an identifier."
    tmptab = ""
    column = ""
    if tk.n == 1
        tk.tokens[1] == _EKQ_TKID || return fail(badcol)
        column = uppercase(String(rstrip(qcol, ' ')))     # leading blanks stay, as in CSPICE
    elseif tk.n == 3
        tk.tokens[1] == _EKQ_TKID ||
            return fail("Invalid table name; name should consist of an identifier.")
        tk.tokens[2] == _EKQ_TKDOT || return fail("Invalid qualified column name; " *
            "table name should be followed by a period.")
        tk.tokens[3] == _EKQ_TKID || return fail(badcol)
        tmptab = String(first(String(_ek_qtext(tk, 1)), _EKQ_TNAMSZ))
        column = String(_ek_qtext(tk, 3))
    else
        return fail("Invalid tokens present in qualified column name. Valid syntax " *
                    "is <column> or <table>.<column>")
    end
    tabidx = colidx = 0
    if tk.n == 3
        tabidx = something(findfirst(==(tmptab), q.tables),
            findfirst(==(tmptab), q.aliases), 0)
        tabidx == 0 && return fail(_ek_qrep("Table name <#> does not match table or " *
            "alias from query.", tmptab))
        colidx = something(findfirst(==(column), tables[q.tabcat[tabidx]].cnames), 0)
        colidx == 0 && return fail(_ek_qrep(_ek_qrep("Column name <#> does not appear " *
            "in the qualifying table <#>.", column), tmptab))
    else
        nmatch = 0
        for i in eachindex(q.tables), (j, name) in enumerate(tables[q.tabcat[i]].cnames)
            if name == column
                nmatch += 1
                colidx = j
                tabidx = i
            end
        end
        nmatch == 0 && return fail(_ek_qrep("Column name <#> does not appear in any " *
            "table in FROM clause of query.", column))
        nmatch > 1 && return fail(_ek_qrep("Column name <#> is ambiguous without a " *
            "qualifying table name.", column))
    end
    return (tabidx, q.tables[tabidx], q.aliases[tabidx], colidx, column, false, "")
end
