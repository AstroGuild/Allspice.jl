# Group EK 2: loaded EK tables and queries (`ekntab`, `ektnam`, `ekccnt`, `ekcii`,
# `ekpsel`, `ekfind`, `eknelt`, `ekgc`, `ekgd`, `ekgi`).
#
# Feature file: features/ek_queries.feature
#
# Three references are used.
#
#   1. CSPICE, for every query of the corpus: the error flag, the error message, the row
#      count, and every value that the fetch calls return, in the order of CSPICE.
#      NeoSpice is run in its compatibility mode here (`_ek_search(...; compat=true)`),
#      which reproduces three defects of CSPICE that change results (see
#      `src/ek/query_run.jl`); it must agree for every query.
#   2. The model that the EK files were written from, for every generated query: the
#      rows are computed by brute force over the model with the comparison rules of the
#      EK query language. The public `ekfind` (correct mode) must give this set of rows.
#   3. The two modes against each other: where the public `ekfind` differs from CSPICE,
#      the query must be one that these defects can touch.
#
# The corpus has hand-written queries (every operator and clause, invalid queries) and
# generated queries on one, two and three tables.

# The code of the EK groups is long and most of it runs once, so Julia compiles the
# module of the tests without optimization while this file runs (5 s less for the five
# EK groups); the last line of the file gives the default level back. The level of
# NeoSpice does not change.
Base.Experimental.@optlevel 0

isdefined(@__MODULE__, :EKLIB) || include("ekfixtures.jl")

# ---------------------------------------------------------------------------------------
# The model as a query reference
# ---------------------------------------------------------------------------------------

"One table of the model for the generator: name, columns, and all rows of all loaded segments."
struct EKQTable
    name::String
    cols::Vector{EKColSpec}
    rows::Vector{Vector{Any}}
end

"""
One relation of a generated WHERE clause: `lhs` and `rhs` are (table position in the FROM
clause, column position) pairs; `rhs` is `nothing` for a comparison with the value
`lit` (and `lit2` for BETWEEN). `text` is the relation in the query language.
"""
struct EKQAtom
    lhs::Tuple{Int,Int}
    op::Symbol
    rhs::Union{Nothing,Tuple{Int,Int}}
    lit::Any
    lit2::Any
    text::String
end

const EKQ_NEGATED = Dict(:eq => :ne, :ne => :eq, :lt => :ge, :ge => :lt, :le => :gt, :gt => :le,
    :like => :unlike, :unlike => :like, :isnull => :notnull, :notnull => :isnull,
    :between => :notbetween, :notbetween => :between)

"Order of two strings as the EK does it: by character code, blank padded, 1024 characters."
function ekq_strcmp(a::AbstractString, b::AbstractString)
    ca = codeunits(a)[1:min(end, 1024)]
    cb = codeunits(b)[1:min(end, 1024)]
    n = max(length(ca), length(cb))
    pa = vcat(ca, fill(UInt8(' '), n - length(ca)))
    pb = vcat(cb, fill(UInt8(' '), n - length(cb)))
    return cmp(pa, pb)
end

ekq_cmp(a::AbstractString, b::AbstractString) = ekq_strcmp(a, b)
ekq_cmp(a::Real, b::Real) = cmp(Float64(a), Float64(b))

"The LIKE operator: `*` any characters, `%` one character, no regard to case or outer blanks."
function ekq_like(s::AbstractString, t::AbstractString)
    s = strip(String(codeunits(s)[1:min(end, 1024)]), ' ')
    t = strip(String(codeunits(t)[1:min(end, 1024)]), ' ')
    io = IOBuffer()
    for c in t
        print(io, c == '*' ? ".*" : c == '%' ? "." : (isletter(c) || isdigit(c)) ? string(c) : "\\" * c)
    end
    return occursin(Regex("\\A(?:" * String(take!(io)) * ")\\z", "is"), s)
end

ekq_relop(op::Symbol, rel::Int) = op === :eq ? rel == 0 : op === :ne ? rel != 0 : op === :lt ? rel < 0 :
                                  op === :le ? rel <= 0 : op === :gt ? rel > 0 : rel >= 0

"True if the relation `a` (with the operator `op`) holds for `rows`, one model row for each table."
function ekq_atom(a::EKQAtom, op::Symbol, rows)
    e = rows[a.lhs[1]][a.lhs[2]]
    lnull = e === nothing
    op === :isnull && return lnull
    op === :notnull && return !lnull
    if a.rhs !== nothing
        e2 = rows[a.rhs[1]][a.rhs[2]]
        rnull = e2 === nothing
        if op === :like || op === :unlike
            m = ekq_like(lnull ? "" : e[1], rnull ? "" : e2[1])
            return op === :like ? m : !m
        end
        rel = lnull && rnull ? 0 : lnull ? -1 : rnull ? 1 : ekq_cmp(e[1], e2[1])
        return ekq_relop(op, rel)
    end
    if op === :like || op === :unlike
        lnull && return false
        m = ekq_like(e[1], a.lit)
        return op === :like ? m : !m
    end
    rel = lnull ? -1 : ekq_cmp(e[1], a.lit)
    if op === :between || op === :notbetween
        # The order of the two bounds has no effect in the EK query language.
        lo, hi = ekq_cmp(a.lit, a.lit2) <= 0 ? (a.lit, a.lit2) : (a.lit2, a.lit)
        inside = !lnull && ekq_cmp(e[1], lo) >= 0 && ekq_cmp(e[1], hi) <= 0
        return op === :between ? inside : !inside
    end
    return ekq_relop(op, rel)
end

"""
Value of the WHERE tree `t` for `rows`. A tree is an `EKQAtom` or `(:and, l, r)`,
`(:or, l, r)`, `(:not, x)`. NOT is not a logical negation of the result: as in the EK
query language it replaces each operator by its opposite (a null entry is less than
every value, so `NOT (A < 1)` is `A >= 1`, which is false for a null entry).
"""
function ekq_eval(t, rows, neg::Bool=false)
    t isa EKQAtom && return ekq_atom(t, neg ? EKQ_NEGATED[t.op] : t.op, rows)
    t[1] === :not && return ekq_eval(t[2], rows, !neg)
    isand = (t[1] === :and) != neg
    return isand ? (ekq_eval(t[2], rows, neg) && ekq_eval(t[3], rows, neg)) :
           (ekq_eval(t[2], rows, neg) || ekq_eval(t[3], rows, neg))
end

"The query text of the WHERE tree `t`."
function ekq_text(rng, @nospecialize(t), parent::Symbol=:top)
    t isa EKQAtom && return t.text
    t[1] === :not && return "NOT (" * ekq_text(rng, t[2]) * ")"
    s = ekq_text(rng, t[2], t[1]) * (t[1] === :and ? " AND " : " OR ") * ekq_text(rng, t[3], t[1])
    return (parent === :and && t[1] === :or) || (parent !== :top && rand(rng) < 0.2) ? "(" * s * ")" : s
end

"A generated query: text, FROM tables, WHERE tree (or `nothing`), SELECT columns, and a note on its kind."
struct EKQGen
    text::String
    tabs::Vector{EKQTable}
    tree::Any
    select::Vector{Tuple{Int,Int}}
    kind::String
end

"An entry of the model as the fetch calls give it: `nothing`, or the values (strings without trailing blanks)."
ekq_entry(e) = e === nothing ? nothing : Any[v isa AbstractString ? String(rstrip(v, ' ')) : v for v in e]

"""
The rows of a query result in an order that does not depend on the engine: the `repr` of
each row, sorted. Two results hold the same rows in any order if and only if their keys
are equal. (`sort(rows; by=repr)` gives the same answer but calls `repr` at every
comparison, which took 33 s of the 61 s of this group.)
"""
ekq_rowkeys(rows) = sort!(String[repr(r) for r in rows])

"The rows that the generated query `g` must give, by brute force over the model (no order)."
function ekq_expected(g::EKQGen)
    out = Vector{Any}[]
    rows = Vector{Vector{Any}}(undef, length(g.tabs))
    function rec(t)
        if t > length(g.tabs)
            (g.tree === nothing || ekq_eval(g.tree, rows)) &&
                push!(out, Any[ekq_entry(rows[a][b]) for (a, b) in g.select])
            return
        end
        for r in g.tabs[t].rows
            rows[t] = r
            rec(t + 1)
        end
    end
    rec(1)
    return out
end

# ---------------------------------------------------------------------------------------
# The generator
# ---------------------------------------------------------------------------------------

const EKQ_OPTEXT = Dict(:eq => ["=", "EQ"], :ne => ["!=", "<>", "^=", "NE"], :lt => ["<", "LT"],
    :le => ["<=", "LE"], :gt => [">", "GT"], :ge => [">=", "GE"])
const EKQ_MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]

"A string literal of the query language for `s`, with one of the two quote characters."
function ekq_quote(rng, s::AbstractString)
    q = rand(rng, ['"', '\''])
    return string(q, replace(s, string(q) => string(q, q)), q)
end

"""
    ekq_number(text) -> Float64

The double that the EK query language assigns to the number `text`. Contract: `text` is
a number in the syntax of CSPICE `nparsd`; the result is the output of `nparsd_` of
CSPICE (the routine that the scanner of queries calls, zzekscan.c line 900) and not of
NeoSpice, so that the model stays independent of the code under test. `nparsd` does not
always give the correctly rounded double (it gives `nextfloat(2.756)` for "2.756"),
thus `parse(Float64, text)` is not the value that a query compares with.
"""
function ekq_number(text::AbstractString)
    x = Ref{Float64}(NaN)
    msg = fill(UInt8(' '), 160)
    ptr = Ref{Cint}(0)
    ccall((:nparsd_, EKLIB), Cint, (Ptr{UInt8}, Ref{Float64}, Ptr{UInt8}, Ref{Cint}, Cint, Cint),
        text, x, msg, ptr, ncodeunits(text), length(msg))
    @assert all(==(UInt8(' ')), msg) && !isnan(x[]) "nparsd accepts the number '$text'"
    return x[]
end

"A value of the scalar column `col` of `tab`: often one that the model holds. Returns (value, text)."
function ekq_literal(rng, ctx, tab::EKQTable, c::Int)
    col = tab.cols[c]
    e = rand(rng, tab.rows)[c]
    v = (e === nothing || rand(rng) < 0.25) ? nothing : e[1]
    if col.dtype === :INT
        v === nothing && (v = rand(rng, -50:250))
        r = rand(rng)
        r < 0.12 && return v + 0.5, string(v + 0.5)             # between two integers
        r < 0.18 && return v - 0.25, string(v - 0.25)
        r < 0.24 && return Float64(v), string(Float64(v))         # a double with an integer value
        return v, string(v)
    elseif col.dtype === :DP
        v === nothing && (v = round(500 * randn(rng); digits=1))
        rand(rng) < 0.15 && isinteger(v) && return v, string(Int(v))
        # The public search reads the text as the correctly rounded double, which is `v`.
        # CSPICE and the compatibility mode read it with `nparsd` (see `ekq_number`), which
        # can be one unit in the last place away, so they can find other rows.
        return v, string(v)
    elseif col.dtype === :TIME
        if v !== nothing && rand(rng) < 0.6
            text = timout(ctx, v, "YYYY MON DD HR:MN:SC.###### TDB ::TDB")
            return str2et(ctx, text), ekq_quote(rng, text)
        end
        y = rand(rng, 1991:2009); m = rand(rng, 1:12); d = rand(rng, 1:28)
        text = rand(rng, ["$y $(EKQ_MONTHS[m]) $d", "$(EKQ_MONTHS[m]) $d, $y 12:00:00 TDB",
            "$y-$(lpad(m, 2, '0'))-$(lpad(d, 2, '0'))T06:07:08", "$y $(EKQ_MONTHS[m]) $d 18:30:00.5"])
        return str2et(ctx, text), ekq_quote(rng, text)
    end
    v === nothing && (v = rand(rng, EK_WORDS))
    v = String(rstrip(String(v)[1:min(end, 60)], ' '))
    isempty(v) && (v = " ")                 # a string of length zero is not a valid token
    return v, ekq_quote(rng, v)
end

"A LIKE template made from a value of the character column `c` of `tab`."
function ekq_template(rng, tab::EKQTable, c::Int)
    e = rand(rng, tab.rows)[c]
    s = e === nothing ? rand(rng, EK_WORDS) : String(strip(String(e[1])[1:min(end, 30)], ' '))
    r = rand(rng)
    isempty(s) && return rand(rng, ["*", "%", "*a*", "%%"])
    r < 0.2 && return s[1:rand(rng, 1:length(s))] * "*"
    r < 0.4 && return "*" * s[rand(rng, 1:length(s)):end]
    r < 0.55 && return "*" * s[rand(rng, 1:length(s))] * "*"
    r < 0.7 && (i = rand(rng, 1:length(s)); return s[1:i-1] * "%" * s[i+1:end])
    r < 0.8 && return uppercase(s)
    r < 0.9 && return "*" * lowercase(s[end:end])
    return rand(rng, ["*", "%%%%%", "a*a", "*e*a*", "%*", "*%e", "A*", "*  ", "%eta"])
end

ekq_scalars(tab::EKQTable) = [i for (i, c) in enumerate(tab.cols) if c.size == 1]

"A relation between the scalar column `c` of FROM table `t` and a value."
function ekq_value_atom(rng, ctx, tabs, names, t::Int, c::Int)
    tab = tabs[t]
    col = tab.cols[c]
    lhs = names[t] * col.name
    r = rand(rng)
    if col.nullok && r < 0.1
        op = rand(rng, [:isnull, :notnull])
        return EKQAtom((t, c), op, nothing, nothing, nothing, lhs * (op === :isnull ? " IS NULL" : " IS NOT NULL"))
    elseif col.dtype === :CHR && r < 0.45
        op = rand(rng, [:like, :like, :unlike])
        tpl = ekq_template(rng, tab, c)
        return EKQAtom((t, c), op, nothing, tpl, nothing, lhs * (op === :like ? " LIKE " : " NOT LIKE ") * ekq_quote(rng, tpl))
    elseif r > 0.88
        a, ta = ekq_literal(rng, ctx, tab, c)
        b, tb = ekq_literal(rng, ctx, tab, c)
        ekq_cmp(a, b) > 0 && rand(rng) < 0.9 && ((a, ta, b, tb) = (b, tb, a, ta))
        op = rand(rng, [:between, :between, :notbetween])
        return EKQAtom((t, c), op, nothing, a, b, lhs * (op === :between ? " BETWEEN " : " NOT BETWEEN ") * ta * " AND " * tb)
    end
    op = rand(rng, [:eq, :eq, :ne, :lt, :le, :gt, :ge])
    v, text = ekq_literal(rng, ctx, tab, c)
    return EKQAtom((t, c), op, nothing, v, nothing, lhs * " " * rand(rng, EKQ_OPTEXT[op]) * " " * text)
end

"Scalar columns of `b` that the column `col` can be compared with."
ekq_partners(col::EKColSpec, b::EKQTable) = [i for i in ekq_scalars(b) if
    (col.dtype in (:INT, :DP) ? b.cols[i].dtype in (:INT, :DP) : b.cols[i].dtype === col.dtype)]

"A relation between a scalar column of table `t1` and one of table `t2`; `ops` are the operators to choose from."
function ekq_column_atom(rng, tabs, names, t1::Int, t2::Int, ops)
    c1 = rand(rng, ekq_scalars(tabs[t1]))
    col = tabs[t1].cols[c1]
    c2 = rand(rng, ekq_partners(col, tabs[t2]))
    op = rand(rng, ops)
    (op in (:like, :unlike) && col.dtype !== :CHR) && (op = :eq)
    optext = op === :like ? "LIKE" : op === :unlike ? "NOT LIKE" : rand(rng, EKQ_OPTEXT[op])
    return EKQAtom((t1, c1), op, (t2, c2), nothing, nothing,
        names[t1] * col.name * " " * optext * " " * names[t2] * tabs[t2].cols[c2].name)
end

"A random tree over the atoms `atoms` (all are used)."
function ekq_tree(rng, atoms::Vector)
    length(atoms) == 1 && return rand(rng) < 0.12 ? (:not, atoms[1]) : atoms[1]
    k = rand(rng, 1:length(atoms)-1)
    t = (rand(rng, [:and, :and, :or]), ekq_tree(rng, atoms[1:k]), ekq_tree(rng, atoms[k+1:end]))
    return rand(rng) < 0.1 ? (:not, t) : t
end

"""
SELECT, FROM, WHERE and ORDER BY text of a generated query. A WHERE tree is a nested
tuple, so each tree has its own type: the functions that only build the text do not
specialize on it (one compiled method for all trees). `ekq_eval` does specialize, since
it runs once for each combination of rows.
"""
function ekq_build(rng, tabs, aliases, names, @nospecialize(tree), kind; arrays::Bool=true)
    all = [(t, c) for t in eachindex(tabs) for c in eachindex(tabs[t].cols) if arrays || tabs[t].cols[c].size == 1]
    select = [rand(rng, all) for _ in 1:rand(rng, 1:4)]
    from = join([isempty(aliases[t]) ? tabs[t].name : tabs[t].name * " " * aliases[t] for t in eachindex(tabs)], ", ")
    text = "SELECT " * join([names[t] * tabs[t].cols[c].name for (t, c) in select], ", ") * " FROM " * from
    tree === nothing || (text *= " WHERE " * ekq_text(rng, tree))
    if rand(rng) < 0.55
        keys = String[]
        for _ in 1:rand(rng, 1:3)
            t = rand(rng, eachindex(tabs))
            c = rand(rng, ekq_scalars(tabs[t]))
            push!(keys, names[t] * tabs[t].cols[c].name * rand(rng, ["", "", " ASC", " DESC", " DESC"]))
        end
        text *= " ORDER BY " * join(keys, ", ")
    end
    return EKQGen(text, tabs, tree, select, kind)
end

"A generated query on one table with one to four relations."
function ekq_single(rng, ctx, tab::EKQTable)
    alias = rand(rng) < 0.3 ? rand(rng, ["A", "T1", "EV"]) : ""
    names = [isempty(alias) ? (rand(rng) < 0.2 ? tab.name * "." : "") : alias * "."]
    tabs = [tab]
    natoms = rand(rng, [0, 1, 1, 1, 2, 2, 2, 3, 3, 4])
    atoms = Any[]
    for _ in 1:natoms
        if rand(rng) < 0.12
            push!(atoms, ekq_column_atom(rng, tabs, names, 1, 1, [:eq, :ne, :lt, :le, :gt, :ge, :like, :unlike]))
        else
            push!(atoms, ekq_value_atom(rng, ctx, tabs, names, 1, rand(rng, ekq_scalars(tab))))
        end
    end
    tree = isempty(atoms) ? nothing : ekq_tree(rng, atoms)
    return ekq_build(rng, tabs, [alias], names, tree, "one table")
end

"""
A generated join of the tables `tabs`: each table after the first is tied to an earlier
one by a relation between columns, with more relations on values; sometimes the whole
clause is one alternative of an OR.
"""
function ekq_join(rng, ctx, tabs::Vector{EKQTable})
    aliases = ["T$i" for i in eachindex(tabs)]
    names = [a * "." for a in aliases]
    parts = Any[]
    for t in 2:length(tabs)
        ops = rand(rng) < 0.65 ? [:eq] : rand(rng, [[:lt, :le, :gt, :ge], [:ne], [:like, :unlike], [:eq]])
        a = rand(rng) < 0.5 ? ekq_column_atom(rng, tabs, names, t, rand(rng, 1:t-1), ops) :
            ekq_column_atom(rng, tabs, names, rand(rng, 1:t-1), t, ops)
        push!(parts, a)
        rand(rng) < 0.25 && push!(parts, ekq_column_atom(rng, tabs, names, t, rand(rng, 1:t-1), [:eq, :ne, :lt, :ge]))
    end
    # Relations on values keep the row sets (and the time of CSPICE) small.
    for t in eachindex(tabs)
        length(tabs[t].rows) > 60 || rand(rng) < 0.3 || continue
        c = findfirst(col -> col.dtype === :INT && col.indexed, tabs[t].cols)
        lo = rand(rng, 1:150)
        push!(parts, EKQAtom((t, c), :between, nothing, lo, lo + 40, names[t] * tabs[t].cols[c].name *
            " BETWEEN $lo AND $(lo + 40)"))
        rand(rng) < 0.3 && push!(parts, ekq_value_atom(rng, ctx, tabs, names, t, rand(rng, ekq_scalars(tabs[t]))))
    end
    tree = foldl((l, r) -> (:and, l, r), shuffle(rng, parts))
    if rand(rng) < 0.25
        extra = ekq_value_atom(rng, ctx, tabs, names, 1, rand(rng, ekq_scalars(tabs[1])))
        tree = (:and, tree, (:or, extra, ekq_value_atom(rng, ctx, tabs, names, length(tabs),
            rand(rng, ekq_scalars(tabs[end])))))
    end
    return ekq_build(rng, tabs, aliases, names, tree, "$(length(tabs)) tables"; arrays=rand(rng) < 0.3)
end

# ---------------------------------------------------------------------------------------
# Hand-written queries
# ---------------------------------------------------------------------------------------

"Valid and invalid queries on the tables EVENTS, GROUPS and LINKS of the standard model."
const EK_HAND_QUERIES = [
    # All rows, SELECT forms
    "SELECT ID FROM EVENTS",
    "select id, name from events",
    "SELECT EVENTS.ID, EVENTS.NAME, EVENTS.VAL, EVENTS.EPOCH FROM EVENTS",
    "SELECT E.ID, E.DESCR, E.LABEL FROM EVENTS E",
    "SELECT VEC, IVEC, TAGS, WORDS, TVEC, DVEC, IPAIR FROM EVENTS WHERE ID < 40",
    "SELECT GID, GNAME, WEIGHT, T0, NOTE, K FROM GROUPS",
    "SELECT EID, GID, TAG, SCORE, WHEN FROM LINKS",
    # Every comparison operator, both spellings, on every data type
    "SELECT ID FROM EVENTS WHERE ID = 17",
    "SELECT ID FROM EVENTS WHERE ID EQ 17",
    "SELECT ID FROM EVENTS WHERE ID != 17",
    "SELECT ID FROM EVENTS WHERE ID <> 17",
    "SELECT ID FROM EVENTS WHERE ID ^= 17",
    "SELECT ID FROM EVENTS WHERE ID NE 17",
    "SELECT ID FROM EVENTS WHERE ID < 17",
    "SELECT ID FROM EVENTS WHERE ID LT 17",
    "SELECT ID FROM EVENTS WHERE ID <= 17",
    "SELECT ID FROM EVENTS WHERE ID LE 17",
    "SELECT ID FROM EVENTS WHERE ID > 170",
    "SELECT ID FROM EVENTS WHERE ID GT 170",
    "SELECT ID FROM EVENTS WHERE ID >= 170",
    "SELECT ID FROM EVENTS WHERE ID GE 170",
    "SELECT ID, GRP FROM EVENTS WHERE GRP = 100",
    "SELECT ID, GRP FROM EVENTS WHERE GRP < -900",
    "SELECT ID, VAL FROM EVENTS WHERE VAL = 3",
    "SELECT ID, VAL FROM EVENTS WHERE VAL = 3.0",
    "SELECT ID, VAL FROM EVENTS WHERE VAL < -1.5E3",
    "SELECT ID, VAL FROM EVENTS WHERE VAL > 1.5D3",
    "SELECT ID, VAL FROM EVENTS WHERE VAL >= -.5 AND VAL <= +.5",
    "SELECT ID, X FROM EVENTS WHERE X > 2000.",
    "SELECT ID, X FROM EVENTS WHERE X = 0",
    "SELECT ID, NAME FROM EVENTS WHERE NAME = 'alpha'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME = \"Alpha\"",
    "SELECT ID, NAME FROM EVENTS WHERE NAME = 'ALPHA   '",
    "SELECT ID, NAME FROM EVENTS WHERE NAME = 'it''s'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME = \"it's\"",
    "SELECT ID, NAME FROM EVENTS WHERE NAME < 'B'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME >= 'p'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME != 'alpha'",
    "SELECT ID, LABEL FROM EVENTS WHERE LABEL = 'Beta gamma'",
    "SELECT ID, LABEL FROM EVENTS WHERE LABEL > 'zeta' ORDER BY LABEL",
    "SELECT ID, DESCR FROM EVENTS WHERE DESCR < 'b' AND ID < 60",
    "SELECT ID, EPOCH FROM EVENTS WHERE EPOCH < '1995 JAN 1'",
    "SELECT ID, EPOCH FROM EVENTS WHERE EPOCH > \"2005 JUL 4 12:00:00\"",
    "SELECT ID, EPOCH FROM EVENTS WHERE EPOCH = '2000 JAN 1 12:00:00 TDB'",
    "SELECT ID, EPOCH FROM EVENTS WHERE EPOCH >= '2000-01-01T00:00:00' AND EPOCH < '2001-01-01T00:00:00'",
    "SELECT ID, T2 FROM EVENTS WHERE T2 <= 'JD 2451545.0'",
    "SELECT ID, T2 FROM EVENTS WHERE T2 > 'March 3, 2003 11:12 PM'",
    "SELECT GID, T0 FROM GROUPS WHERE T0 < '1 JAN 1992'",
    # A double value against an indexed integer column (defect of CSPICE for fractions)
    "SELECT ID FROM EVENTS WHERE ID < 17.0",
    "SELECT ID FROM EVENTS WHERE ID < 17.4",
    "SELECT ID FROM EVENTS WHERE ID <= 17.6",
    "SELECT ID FROM EVENTS WHERE ID > 17.5",
    "SELECT ID FROM EVENTS WHERE ID >= 16.5",
    "SELECT ID FROM EVENTS WHERE ID = 17.5",
    "SELECT ID FROM EVENTS WHERE ID = 17.2",
    "SELECT GID FROM GROUPS WHERE GID < 30.7",
    "SELECT GID FROM GROUPS WHERE GID BETWEEN 10.5 AND 20.5",
    "SELECT ID, GRP FROM EVENTS WHERE GRP < 17.4",
    # An integer value against an indexed double column
    "SELECT ID, VAL FROM EVENTS WHERE VAL < 0",
    "SELECT ID, VAL FROM EVENTS WHERE VAL >= 10",
    "SELECT GID, WEIGHT FROM GROUPS WHERE WEIGHT = 5",
    "SELECT GID, WEIGHT FROM GROUPS WHERE WEIGHT BETWEEN -3 AND 3",
    # LIKE and NOT LIKE
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE 'a*'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE 'A*'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE '*a'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE '*e*'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE '%eta'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE '%%'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE '*'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE 'a%b'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE 'a*b'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE '*mma'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE 'lead'",
    "SELECT ID, NAME FROM EVENTS WHERE NAME LIKE '  ALPHA  '",
    "SELECT ID, NAME FROM EVENTS WHERE NAME NOT LIKE '*a*'",
    "SELECT ID, NAME FROM EVENTS WHERE NOT NAME LIKE '*a*'",
    "SELECT ID, LABEL FROM EVENTS WHERE LABEL LIKE '*gamma*delta*'",
    "SELECT ID, DESCR FROM EVENTS WHERE DESCR LIKE '*z'",
    "SELECT GID, GNAME FROM GROUPS WHERE GNAME LIKE 'mixed*'",
    "SELECT GID, NOTE FROM GROUPS WHERE NOTE LIKE '%%'",
    "SELECT ID FROM EVENTS WHERE NAME LIKE LABEL",
    "SELECT ID FROM EVENTS WHERE NAME NOT LIKE LABEL AND ID < 30",
    # Null values
    "SELECT ID, NAME FROM EVENTS WHERE NAME IS NULL",
    "SELECT ID, NAME FROM EVENTS WHERE NAME IS NOT NULL AND ID < 25",
    "SELECT ID, VAL FROM EVENTS WHERE VAL IS NULL ORDER BY ID",
    "SELECT ID, T2 FROM EVENTS WHERE T2 IS NULL AND GRP IS NULL",
    "SELECT ID FROM EVENTS WHERE NOT (VAL IS NULL) AND ID < 25",
    "SELECT ID FROM EVENTS WHERE NOT VAL IS NOT NULL",
    "SELECT ID, VAL FROM EVENTS WHERE VAL < -1E9",
    "SELECT ID, NAME FROM EVENTS WHERE NAME < ' '",
    "SELECT ID, NAME FROM EVENTS WHERE NAME != 'x' AND NAME NOT LIKE '*' ",
    "SELECT GID, K FROM GROUPS WHERE K IS NULL",
    "SELECT GID, GNAME, WEIGHT FROM GROUPS WHERE GNAME IS NULL OR WEIGHT IS NULL",
    "SELECT ID FROM EVENTS WHERE ID IS NULL",
    # BETWEEN
    "SELECT ID FROM EVENTS WHERE ID BETWEEN 10 AND 20",
    "SELECT ID FROM EVENTS WHERE ID NOT BETWEEN 3 AND 197",
    "SELECT ID FROM EVENTS WHERE ID BETWEEN 20 AND 10",
    "SELECT ID, VAL FROM EVENTS WHERE VAL BETWEEN -1.5 AND 1.5",
    "SELECT ID, NAME FROM EVENTS WHERE NAME BETWEEN 'a' AND 'c'",
    "SELECT ID, EPOCH FROM EVENTS WHERE EPOCH BETWEEN '1999 JAN 1' AND '2000 JAN 1'",
    "SELECT ID, EPOCH FROM EVENTS WHERE EPOCH NOT BETWEEN '1991 JAN 1' AND '2009 JAN 1'",
    "SELECT ID FROM EVENTS WHERE NOT (ID BETWEEN 5 AND 195)",
    "SELECT ID, X, VAL FROM EVENTS WHERE X BETWEEN VAL AND 100",
    # AND, OR, NOT, parentheses
    "SELECT ID FROM EVENTS WHERE ID < 10 OR ID > 190",
    "SELECT ID FROM EVENTS WHERE ID < 10 OR ID < 20",
    "SELECT ID FROM EVENTS WHERE ID < 20 OR ID < 10",
    "SELECT ID FROM EVENTS WHERE ID < 10 OR ID > 190 OR NAME LIKE 'a*'",
    "SELECT ID FROM EVENTS WHERE (ID < 10 OR ID > 190) AND (VAL > 0 OR VAL IS NULL)",
    "SELECT ID FROM EVENTS WHERE ID < 50 AND NOT (NAME = 'alpha' OR NAME = 'beta')",
    "SELECT ID FROM EVENTS WHERE NOT (ID >= 10 AND ID <= 190)",
    "SELECT ID FROM EVENTS WHERE NOT NOT ID < 10",
    "SELECT ID FROM EVENTS WHERE ((ID < 10))",
    "SELECT ID FROM EVENTS WHERE ID < 10 AND ID > 20",
    "SELECT ID FROM EVENTS WHERE ID < 100 AND ID < 50 AND ID < 20 AND ID >= 5",
    "SELECT ID, NAME, VAL FROM EVENTS WHERE ID < 30 AND NAME > 'b' AND VAL > 0",
    "SELECT ID, NAME, VAL FROM EVENTS WHERE NAME = 'alpha' AND VAL < 0",
    "SELECT ID FROM EVENTS WHERE NAME = 'alpha' OR NAME = 'ALPHA' OR NAME = 'Alpha' ORDER BY ID",
    "SELECT GID FROM GROUPS WHERE GID < 30 OR WEIGHT > 100.0",
    "SELECT GID FROM GROUPS WHERE GID < 30 OR GID < 60 OR GID < 90",
    "SELECT GID, GNAME FROM GROUPS WHERE GNAME = 'alpha' OR GNAME = 'beta' OR GID = 7",
    "SELECT GID FROM GROUPS WHERE GID > 100 OR GID > 50",
    # Relations between columns of one table
    "SELECT ID, GRP FROM EVENTS WHERE ID = GRP",
    "SELECT ID, GRP FROM EVENTS WHERE ID < GRP AND ID < 20",
    "SELECT ID, VAL, X FROM EVENTS WHERE VAL > X",
    "SELECT ID, VAL FROM EVENTS WHERE ID > VAL AND VAL > 0",
    "SELECT ID, EPOCH, T2 FROM EVENTS WHERE EPOCH < T2",
    "SELECT ID, NAME, LABEL FROM EVENTS WHERE NAME = LABEL",
    "SELECT ID, NAME, LABEL FROM EVENTS WHERE NAME < LABEL AND ID < 40",
    "SELECT GID, K FROM GROUPS WHERE GID >= K",
    # ORDER BY
    "SELECT ID FROM EVENTS ORDER BY ID",
    "SELECT ID FROM EVENTS ORDER BY ID DESC",
    "SELECT ID, NAME FROM EVENTS ORDER BY NAME, ID",
    "SELECT ID, NAME FROM EVENTS ORDER BY NAME DESC, ID ASC",
    "SELECT ID, VAL FROM EVENTS ORDER BY VAL",
    "SELECT ID, VAL FROM EVENTS ORDER BY VAL DESC",
    "SELECT ID, EPOCH FROM EVENTS WHERE ID < 60 ORDER BY EPOCH DESC",
    "SELECT ID, GRP, NAME FROM EVENTS ORDER BY GRP DESC, NAME DESC, ID DESC",
    "SELECT ID FROM EVENTS ORDER BY DESCR",
    "SELECT ID, LABEL FROM EVENTS WHERE ID > 150 ORDER BY LABEL DESC",
    "SELECT GID, GNAME FROM GROUPS ORDER BY GNAME, GID",
    "SELECT GID, WEIGHT FROM GROUPS ORDER BY WEIGHT DESC, GID DESC",
    "SELECT GID FROM GROUPS ORDER BY K, T0",
    "SELECT EID, TAG FROM LINKS ORDER BY TAG DESC, EID",
    "SELECT ID FROM EVENTS WHERE ID < 10 OR ID > 190 ORDER BY ID DESC",
    # Joins of two tables
    "SELECT E.ID, G.GID, G.GNAME FROM EVENTS E, GROUPS G WHERE E.GRP = G.GID AND E.ID < 40",
    "SELECT E.ID, G.GID FROM EVENTS E, GROUPS G WHERE G.GID = E.GRP AND E.ID < 40",
    "SELECT E.ID, G.GID FROM GROUPS G, EVENTS E WHERE E.GRP = G.GID AND E.ID < 40",
    "SELECT E.ID, G.GNAME FROM EVENTS E, GROUPS G WHERE E.GRP = G.GID ORDER BY G.GNAME, E.ID",
    "SELECT E.ID, G.GNAME FROM EVENTS E, GROUPS G WHERE E.ID = G.GID AND E.NAME = G.GNAME",
    "SELECT E.NAME, G.GNAME FROM EVENTS E, GROUPS G WHERE E.NAME = G.GNAME AND E.ID < 30 AND G.GID < 60",
    "SELECT E.ID, G.GID FROM EVENTS E, GROUPS G WHERE E.VAL = G.WEIGHT",
    "SELECT E.ID, G.GID FROM EVENTS E, GROUPS G WHERE E.VAL = G.GID AND E.VAL > 0",
    "SELECT E.ID, G.GID FROM EVENTS E, GROUPS G WHERE E.EPOCH = G.T0",
    "SELECT E.ID, G.GID FROM EVENTS E, GROUPS G WHERE E.T2 < G.T0 AND E.ID < 6 AND G.GID < 30",
    "SELECT E.ID, L.EID FROM EVENTS E, LINKS L WHERE E.ID < L.EID AND E.ID < 12",
    "SELECT E.ID, L.EID FROM EVENTS E, LINKS L WHERE E.ID <= L.EID AND E.ID > 150",
    "SELECT E.ID, L.EID FROM EVENTS E, LINKS L WHERE E.ID > L.EID AND E.ID < 40",
    "SELECT E.ID, L.EID FROM EVENTS E, LINKS L WHERE L.EID >= E.ID AND E.ID > 170",
    "SELECT E.ID, L.EID FROM EVENTS E, LINKS L WHERE E.ID != L.EID AND E.ID < 4",
    "SELECT E.ID, L.TAG FROM EVENTS E, LINKS L WHERE E.NAME LIKE L.TAG",
    "SELECT E.ID, L.TAG FROM EVENTS E, LINKS L WHERE E.NAME NOT LIKE L.TAG AND E.ID < 8",
    "SELECT E.ID, L.TAG FROM EVENTS E, LINKS L WHERE E.ID = L.EID",
    "SELECT E.ID, L.TAG FROM EVENTS E, LINKS L WHERE E.ID = L.EID AND E.GRP = L.GID",
    "SELECT E.ID, L.TAG FROM EVENTS E, LINKS L WHERE E.ID = L.EID OR E.GRP = L.GID AND E.ID < 50",
    "SELECT E.ID, L.EID FROM EVENTS E, LINKS L WHERE E.ID < 3 AND L.EID < 40",
    "SELECT L.EID, G.GID FROM LINKS L, GROUPS G WHERE L.SCORE < G.WEIGHT AND G.GID < 12",
    "SELECT L.EID, G.GID FROM LINKS L, GROUPS G WHERE L.SCORE >= G.WEIGHT AND G.GID < 12 ORDER BY G.GID DESC, L.EID",
    "SELECT L.EID, L.WHEN, G.T0 FROM LINKS L, GROUPS G WHERE L.WHEN > G.T0 AND G.GID < 9",
    "SELECT L.TAG, G.GNAME FROM LINKS L, GROUPS G WHERE L.TAG = G.GNAME",
    "SELECT L.TAG, G.GNAME FROM LINKS L, GROUPS G WHERE L.TAG <= G.GNAME AND G.GID < 9",
    "SELECT A.ID, B.ID FROM EVENTS A, EVENTS B WHERE A.ID = B.GRP AND A.ID < 25",
    "SELECT A.ID, B.ID, A.NAME FROM EVENTS A, EVENTS B WHERE A.NAME = B.NAME AND A.ID < B.ID AND A.ID < 12",
    "SELECT A.EID, B.EID FROM LINKS A, LINKS B WHERE A.EID < B.EID",
    "SELECT A.EID, B.EID FROM LINKS A, LINKS B",
    "SELECT A.EID, B.EID FROM LINKS A, LINKS B WHERE A.TAG = B.TAG ORDER BY A.EID, B.EID DESC",
    "SELECT A.GID, B.GID FROM GROUPS A, GROUPS B WHERE A.GID = B.K AND A.GID < 60",
    "SELECT A.GID, B.GID FROM GROUPS A, GROUPS B WHERE A.WEIGHT = B.WEIGHT AND A.GID < B.GID",
    "SELECT EVENTS.ID, GROUPS.GID FROM EVENTS, GROUPS WHERE EVENTS.GRP = GROUPS.GID AND GROUPS.GID < 20",
    "SELECT ID, GNAME FROM EVENTS, GROUPS WHERE GRP = EVENTS.ID AND EVENTS.ID = GROUPS.GID",
    "SELECT E.ID, G.GID FROM EVENTS E, GROUPS G WHERE (E.GRP = G.GID AND E.ID < 30) OR (E.ID = G.GID AND E.ID > 180)",
    # Joins of three tables
    "SELECT E.ID, L.TAG, G.GID FROM EVENTS E, LINKS L, GROUPS G WHERE E.ID = L.EID AND L.GID = G.GID",
    "SELECT E.ID, L.TAG, G.GNAME FROM EVENTS E, LINKS L, GROUPS G WHERE E.ID = L.EID AND L.GID = G.GID ORDER BY G.GNAME DESC, E.ID",
    "SELECT E.ID, L.TAG, G.GID FROM LINKS L, GROUPS G, EVENTS E WHERE E.ID = L.EID AND L.GID = G.GID AND E.VAL > 0",
    "SELECT E.ID, L.EID, G.GID FROM EVENTS E, LINKS L, GROUPS G WHERE E.ID = L.EID AND E.GRP = G.GID",
    "SELECT E.ID, L.EID, G.GID FROM EVENTS E, LINKS L, GROUPS G WHERE E.ID < L.EID AND L.GID > G.GID AND E.ID < 5 AND G.GID < 9",
    "SELECT E.ID, L.EID, G.GID FROM EVENTS E, LINKS L, GROUPS G WHERE E.ID = G.GID AND L.EID = G.GID",
    "SELECT A.EID, B.EID, C.EID FROM LINKS A, LINKS B, LINKS C WHERE A.EID < B.EID AND B.EID < C.EID AND A.EID < 60",
    "SELECT A.EID, B.EID, C.EID FROM LINKS A, LINKS B, LINKS C WHERE A.GID = B.GID AND B.GID = C.GID",
    "SELECT E.ID, L.EID, G.GID FROM EVENTS E, LINKS L, GROUPS G WHERE E.ID = L.EID AND L.GID = G.GID OR E.ID < 2 AND L.EID < 30 AND G.GID < 4",
    # Invalid queries: the error flag and the message must agree
    "   ",
    "SELECT",
    "SELECT ID",
    "SELECT ID FROM",
    "SELECT ID FROM NOTABLE",
    "SELECT FOO FROM EVENTS",
    "SELECT ID FROM EVENTS WHERE",
    "SELECT ID FROM EVENTS WHERE ID",
    "SELECT ID FROM EVENTS WHERE ID <",
    "SELECT ID FROM EVENTS WHERE ID < 'abc'",
    "SELECT ID FROM EVENTS WHERE NAME < 3",
    "SELECT ID FROM EVENTS WHERE NAME LIKE 3",
    "SELECT ID FROM EVENTS WHERE ID LIKE 'a*'",
    "SELECT ID FROM EVENTS WHERE EPOCH < 3.0",
    "SELECT ID FROM EVENTS WHERE EPOCH < 'not a time'",
    "SELECT ID FROM EVENTS ORDER BY",
    "SELECT ID FROM EVENTS ORDER ID",
    "SELECT ID FROM EVENTS WHERE (ID < 3",
    "SELECT ID FROM EVENTS WHERE ID < 3)",
    "SELECT ID FROM EVENTS WHERE ID < 3 AND",
    "SELECT ID FROM EVENTS WHERE ID BETWEEN 3",
    "SELECT ID FROM EVENTS WHERE ID BETWEEN 3 AND",
    "SELECT ID FROM EVENTS WHERE NAME = 'unfinished",
    "SELECT ID FROM EVENTS WHERE NAME = ''",
    "SELECT ID FROM EVENTS, GROUPS WHERE GID = EID",
    "SELECT GID FROM GROUPS, LINKS",
    "SELECT E.ID FROM EVENTS E, GROUPS E",
    "SELECT X.ID FROM EVENTS E",
    "SELECT ID, FROM EVENTS",
    "SELECT * FROM EVENTS",
    "SELECT COUNT(*) FROM EVENTS",
    "SELECT ID FROM EVENTS GROUP BY ID",
    "SELECT ID FROM EVENTS WHERE ID < 3 ORDER BY NAME DESC ASC",
    "SELECT ID FROM EVENTS WHERE ID = = 3",
    "SELECT ID FROM EVENTS WHERE NAME = LABEL = DESCR",
    "SELECT ID FROM EVENTS WHERE ID < 1E400",
    "SELECT ID FROM EVENTS WHERE ID # 3",
    "SELECT ID FROM EVENTS WHERE NAME IS 3",
    "SELECT ID FROM EVENTS WHERE E.ID < 3",
    "SELECT ID FROM EVENTS WHERE EVENTS.NOPE < 3",
    "FROM EVENTS SELECT ID",
    "SELECT ID FROM EVENTS WHERE EPOCH = NAME",
    "SELECT ID FROM EVENTS WHERE ID = EPOCH",
    "SELECT E.ID FROM EVENTS E, GROUPS G WHERE E.NAME = G.GID",
]

# ---------------------------------------------------------------------------------------
# Comparison of one query
# ---------------------------------------------------------------------------------------

"""
Run `query` in CSPICE and in NeoSpice (both modes). Returns a named tuple: `cerr`/`nerr`
the error flags, `msgsame`, `nc` the row count of CSPICE, `compat` true if the
compatibility mode gives the count and every value of CSPICE in the same order (or both
engines give the same error and message), `default` the same for the public `ekfind`,
`sameset` true if the public `ekfind` gives the rows of CSPICE in another order,
`nd` its row count, and `values` its table of values.
"""
function ek_compare_query(ctx, query::String)
    signal = c_kind(() -> c_ekfind(query))
    if signal !== nothing
        # CSPICE signals an error (it does not set the error flag): NeoSpice throws the same kind.
        same = neo_kind(() -> ekfind(ctx, query)) === signal
        return (cerr=true, nerr=true, msgsame=same, nc=0, compat=same, default=same, sameset=same, nd=0,
            values=Vector{Any}[], cmsg=string(signal), nmsg="", signal=signal)
    end
    nc, cerr, cmsg = c_ekfind(query)
    ks = NeoSpice.kernelset(ctx)
    eq, nerr, nmsg = NeoSpice._ek_encode_query(ks.ektables, ks, NeoSpice.timedefaults(ctx), query)
    res, derr, dmsg = ekfind(ctx, query)
    @assert derr == nerr && dmsg == nmsg
    if cerr || nerr
        same = cerr == nerr && cmsg == nmsg
        return (cerr=cerr, nerr=nerr, msgsame=cmsg == nmsg, nc=nc, compat=same, default=same, sameset=same,
            nd=res.nmrows, values=Vector{Any}[], cmsg=cmsg, nmsg=nmsg, signal=nothing)
    end
    rc = NeoSpice._ek_search(ks, eq; compat=true)
    cvals = c_ek_result(nc, rc.seltype)
    compat = rc.nmrows == nc && neo_ek_result(rc) == cvals
    dvals = neo_ek_result(res)
    default = dvals == cvals
    return (cerr=false, nerr=false, msgsame=true, nc=nc, compat=compat, default=default,
        sameset=default || ekq_rowkeys(dvals) == ekq_rowkeys(cvals), nd=res.nmrows, values=dvals, cmsg="",
        nmsg="", signal=nothing)
end

"True if one of the defects of CSPICE can change the result of the query `q`."
ek_defect_possible(q::String) = occursin(r"\bOR\b|NOT|\d\.\d*[1-9]|BETWEEN"i, q)

group("EK queries: catalogue, query corpus against CSPICE and against the model", Kernels.LSK) do
    model, paths = ek_standard_files(SEED + 1; scale=0.15)
    kernels = [Kernels.LSK, paths...]
    oracle_load(kernels)
    ctx = neo_context(kernels)
    allsegs = [s for (_, segs) in model for s in segs]
    qtable(name) = EKQTable(name, first(s.cols for s in allsegs if s.table == name),
        reduce(vcat, [s.rows for s in allsegs if s.table == name]))
    events, groups, links = qtable("EVENTS"), qtable("GROUPS"), qtable("LINKS")

    ek_part("kernel list: an EK file is a kernel of type EK") do
        @test ktotal(ctx, "EK") == SPICE.ktotal("EK") == 3
        @test ktotal(ctx, "ALL") == SPICE.ktotal("ALL") == 4
        for i in 1:3
            file, filtyp, source, _ = kdata(ctx, i, "EK")
            @test (file, filtyp, source) == SPICE.kdata(i, "EK")[1:3]
            @test file == paths[i] && filtyp == "EK"
            @test kinfo(ctx, file)[1] == "EK"
            @test getfat(file) == ("DAS", "EK") == SPICE.getfat(file)
        end
    end

    ek_part("catalogue: ekntab, ektnam, ekccnt, ekcii") do
        @test ekntab(ctx) == c_ekntab() == 3
        for n in 1:ekntab(ctx)
            table = ektnam(ctx, n)
            @test table == c_ektnam(n)
            @test ekccnt(ctx, table) == c_ekccnt(table)
            @test ekccnt(ctx, lowercase(table)) == c_ekccnt(lowercase(table))
            for c in 1:ekccnt(ctx, table)
                name, a = ekcii(ctx, table, c)
                cname, ca = c_ekcii(table, c)
                @test name == cname
                # CSPICE gives the string length 1 for a column that is not a string.
                @test (a.cclass, a.dtype, a.dtype === :CHR ? a.strlen : 1, a.size, a.indexd, a.nullok) == ca
            end
        end
        @test [ektnam(ctx, n) for n in 1:3] == ["EVENTS", "GROUPS", "LINKS"]
        for (f, cf) in ((() -> ektnam(ctx, 0), () -> ccall((:ektnam_c, EKLIB), Cvoid, (Cint, Cint, Ptr{UInt8}), -1, 65, zeros(UInt8, 65))),
            (() -> ektnam(ctx, 4), () -> ccall((:ektnam_c, EKLIB), Cvoid, (Cint, Cint, Ptr{UInt8}), 3, 65, zeros(UInt8, 65))),
            (() -> ekccnt(ctx, "NOTABLE"), () -> ccall((:ekccnt_c, EKLIB), Cvoid, (Cstring, Ref{Cint}), "NOTABLE", Ref{Cint}(0))),
            (() -> ekcii(ctx, "NOTABLE", 1), () -> ccall((:ekcii_c, EKLIB), Cvoid, (Cstring, Cint, Cint, Ptr{UInt8}, Ptr{Cint}), "NOTABLE", 0, 33, zeros(UInt8, 33), zeros(Cint, 6))),
            (() -> ekcii(ctx, "LINKS", 6), () -> ccall((:ekcii_c, EKLIB), Cvoid, (Cstring, Cint, Cint, Ptr{UInt8}, Ptr{Cint}), "LINKS", 5, 33, zeros(UInt8, 33), zeros(Cint, 6))),
            (() -> ekcii(ctx, "LINKS", 0), () -> ccall((:ekcii_c, EKLIB), Cvoid, (Cstring, Cint, Cint, Ptr{UInt8}, Ptr{Cint}), "LINKS", -1, 33, zeros(UInt8, 33), zeros(Cint, 6))))
            cf()
            kind = c_pending_kind()
            @test kind !== nothing
            @test neo_kind(f) == kind
        end
        empty = SpiceContext()
        @test ekntab(empty) == 0
        @test neo_kind(() -> ektnam(empty, 1)) == :noloadedfiles
        @test neo_kind(() -> ekccnt(empty, "EVENTS")) == :tablenotloaded
    end

    ek_part("ekpsel: the SELECT clause") do
        for q in ["SELECT ID, NAME, VAL, EPOCH, VEC FROM EVENTS",
            "select  e.id ,g.gname, G.T0   from events e, groups g where e.grp = g.gid order by e.id",
            "SELECT EVENTS.LABEL, IVEC FROM EVENTS WHERE ID < 3", "SELECT FOO FROM EVENTS", "SELECT ID FROM"]
            a = ekpsel(ctx, q)
            b = c_ekpsel(q)
            @test a == b
        end
    end

    rng = MersenneTwister(SEED + 2)
    generated = EKQGen[]
    # The model answers a join by brute force over all combinations of rows (1.8 million
    # for EVENTS, LINKS and GROUPS), which is most of the time of the generated queries.
    for _ in 1:nsample(60, 110)
        push!(generated, ekq_single(rng, ctx, events))
    end
    for _ in 1:nsample(25, 50)
        push!(generated, ekq_single(rng, ctx, groups))
    end
    for _ in 1:nsample(10, 20)
        push!(generated, ekq_single(rng, ctx, links))
    end
    for _ in 1:nsample(30, 70)
        push!(generated, ekq_join(rng, ctx, shuffle(rng, rand(rng, [[events, groups], [events, links],
            [links, groups], [events, events], [groups, groups], [links, links]]))))
    end
    for _ in 1:nsample(4, 30)
        push!(generated, ekq_join(rng, ctx, shuffle(rng, rand(rng, [[events, links, groups],
            [links, links, groups], [events, links, links]]))))
    end

    stats = Dict{String,Vector{Int}}()          # kind => [queries, compat equal, default equal, same rows, errors]
    examples = String[]
    function record(kind, q, r)
        s = get!(stats, kind, zeros(Int, 5))
        s[1] += 1
        s[2] += r.compat
        s[3] += r.default
        s[4] += r.sameset
        s[5] += r.cerr
        if !r.compat
            println("  differs from CSPICE in compatibility mode: ", q, "\n     CSPICE: ", r.nc, " rows",
                r.cerr ? " error " * r.cmsg : "", "; NeoSpice: ", r.nerr ? "error " * r.nmsg : "$(r.nd) rows in correct mode")
        elseif !r.default && length(examples) < 8
            push!(examples, "    $q\n      CSPICE $(r.nc) rows, NeoSpice $(r.nd) rows" *
                (r.sameset ? " (same rows, other order)" : ""))
        end
        return nothing
    end

    ek_part("hand-written queries against CSPICE") do
        for q in EK_HAND_QUERIES
            r = ek_compare_query(ctx, q)
            record(r.cerr ? "hand-written, invalid" : "hand-written, valid", q, r)
            @test r.compat
            @test r.cerr == r.nerr
            @test r.msgsame
            r.compat && !r.default && @test ek_defect_possible(q)
        end
    end

    ek_part("generated queries against CSPICE and against the model") do
        badmodel = 0
        for g in generated
            r = ek_compare_query(ctx, g.text)
            if r.signal !== nothing
                # A query that CSPICE cannot encode (seen: NOT over a BETWEEN with times
                # fills its buffer of numbers); both engines raise the same error.
                record("generated, error signalled", g.text, r)
                @test r.compat
                continue
            end
            record("generated, " * g.kind, g.text, r)
            @test r.compat
            @test !r.nerr
            r.nerr && (println("  unexpected error: ", g.text, "\n     ", r.nmsg); continue)
            r.compat && !r.default && @test ek_defect_possible(g.text)
            expect = ekq_expected(g)
            ok = ekq_rowkeys(expect) == ekq_rowkeys(r.values)
            if !ok
                badmodel += 1
                badmodel <= 5 && println("  differs from the model: ", g.text, "\n     model ", length(expect),
                    " rows, NeoSpice ", r.nd, " rows, CSPICE ", r.nc, " rows")
            end
        end
        @test badmodel == 0
    end

    corpus_line(kind, counts) = "  " * rpad(kind, 28) * join(" " * lpad(n, 4) for n in counts)
    total = zeros(Int, 5)
    measured("EK query corpus (queries / equal to CSPICE in compatibility mode / public ekfind equal to " *
             "CSPICE / public ekfind same rows in any order / invalid):")
    for kind in sort(collect(keys(stats)))
        s = stats[kind]
        total .+= s
        measured(corpus_line(kind, s))
    end
    measured(corpus_line("all", total))
    @test total[1] >= 200
    @test total[2] == total[1]
    isempty(examples) || println("  queries where the public ekfind differs from CSPICE (defects of CSPICE):\n",
        join(examples, "\n"))

    ek_part("the defects of CSPICE, against the model") do
        ids = [r[1][1] for r in events.rows]
        # (b) A double value with a fraction against an indexed integer column.
        for (q, pred) in (("SELECT ID FROM EVENTS WHERE ID < 17.4", x -> x < 17.4),
            ("SELECT ID FROM EVENTS WHERE ID <= 17.6", x -> x <= 17.6),
            ("SELECT ID FROM EVENTS WHERE ID > 17.5", x -> x > 17.5),
            ("SELECT ID FROM EVENTS WHERE ID >= 16.5", x -> x >= 16.5),
            ("SELECT ID FROM EVENTS WHERE ID = 17.5", x -> false))
            res, err, _ = ekfind(ctx, q)
            @test !err
            @test res.nmrows == count(pred, ids)
            @test sort(ekcolumn(res, 1)) == sort(filter(pred, ids))
        end
        # (a) A row that two alternatives of an OR give comes once.
        res, _, _ = ekfind(ctx, "SELECT ID FROM EVENTS WHERE ID < 20 OR ID < 10")
        @test sort(ekcolumn(res, 1)) == sort(filter(<(20), ids))
        gids = [r[1][1] for r in groups.rows]
        res, _, _ = ekfind(ctx, "SELECT GID FROM GROUPS WHERE GID > 100 OR GID > 50")
        @test sort(ekcolumn(res, 1)) == sort(filter(>(50), gids))
    end

    ek_part("a number of a query: correctly rounded in ekfind, as nparsd in the compatibility mode") do
        # `nparsd` (zzekscan.c line 900) reads "2.756" as the double after 2.756, so in
        # CSPICE a row with the value 2.756 is not equal to the number 2.756 of a query.
        # The public search compares with the correctly rounded double; the compatibility
        # mode keeps the value of CSPICE. Seed 777 has such a row.
        x = ekq_number("2.756")
        @test x === nextfloat(2.756)
        @test NeoSpice._nparsd("2.756") === x
        @test NeoSpice._ek_qexact(codeunits("2.756"), x) === 2.756
        @test NeoSpice._ek_qexact(codeunits("2.756D0"), x) === 2.756
        @test NeoSpice._ek_qexact(codeunits("17"), 17.0) === 17.0
        ks = NeoSpice.kernelset(ctx)
        vcol = findfirst(c -> c.name == "VAL", events.cols)
        vals = Float64[r[vcol][1] for r in events.rows if r[vcol] !== nothing]
        for (op, f) in (("=", ==), (">=", >=), (">", >))
            q = "SELECT ID FROM EVENTS WHERE VAL $op 2.756"
            eq, err, _ = NeoSpice._ek_encode_query(ks.ektables, ks, NeoSpice.timedefaults(ctx), q)
            @test !err
            con = only(eq.constraints)
            @test con.dval === x && con.dexact === 2.756
            @test ekfind(ctx, q)[1].nmrows == count(v -> f(v, 2.756), vals)
            @test NeoSpice._ek_search(ks, eq; compat=true).nmrows == count(v -> f(v, x), vals) == c_ekfind(q)[1]
        end
    end

    ek_part("fetch calls: indices and data types") do
        q = "SELECT ID, NAME, VAL, EPOCH, IVEC, TAGS FROM EVENTS WHERE ID < 30 ORDER BY ID"
        res, err, msg = ekfind(ctx, q)
        n, cerr, _ = c_ekfind(q)
        @test !err && isempty(msg) && res.nmrows == n == length(res) && n > 0
        @test res isa EKResult
        directi(sel, row, elt) = () -> ccall((:ekgi_c, EKLIB), Cvoid, (Cint, Cint, Cint, Ref{Cint}, Ref{Cint}, Ref{Cint}),
            sel - 1, row - 1, elt - 1, Ref{Cint}(0), Ref{Cint}(0), Ref{Cint}(0))
        directd(sel, row, elt) = () -> ccall((:ekgd_c, EKLIB), Cvoid, (Cint, Cint, Cint, Ref{Float64}, Ref{Cint}, Ref{Cint}),
            sel - 1, row - 1, elt - 1, Ref{Float64}(0), Ref{Cint}(0), Ref{Cint}(0))
        directc(sel, row, elt) = () -> ccall((:ekgc_c, EKLIB), Cvoid, (Cint, Cint, Cint, Cint, Ptr{UInt8}, Ref{Cint},
                Ref{Cint}), sel - 1, row - 1, elt - 1, 100, zeros(UInt8, 100), Ref{Cint}(0), Ref{Cint}(0))
        for (f, cf) in ((() -> ekgi(res, 2, 1), directi(2, 1, 1)),
            (() -> ekgi(res, 3, 1), directi(3, 1, 1)),
            (() -> ekgd(res, 1, 1), directd(1, 1, 1)),
            (() -> ekgd(res, 2, 1), directd(2, 1, 1)),
            (() -> ekgc(res, 1, 1), directc(1, 1, 1)),
            (() -> ekgc(res, 4, 1), directc(4, 1, 1)),
            (() -> ekgi(res, 1, 0), directi(1, 0, 1)),
            (() -> ekgi(res, 1, n + 1), directi(1, n + 1, 1)),
            (() -> ekgi(res, 0, 1), directi(0, 1, 1)),
            (() -> ekgi(res, 7, 1), directi(7, 1, 1)),
            (() -> ekgi(res, 1, 1, 0), directi(1, 1, 0)),
            (() -> ekgd(res, 3, 1, 0), directd(3, 1, 0)),
            (() -> ekgc(res, 2, 1, 0), directc(2, 1, 0)),
            (() -> eknelt(res, 7, 1), () -> ccall((:eknelt_c, EKLIB), Cint, (Cint, Cint), 6, 0)),
            (() -> eknelt(res, 1, n + 1), () -> ccall((:eknelt_c, EKLIB), Cint, (Cint, Cint), 0, n)))
            cf()
            kind = c_pending_kind()
            @test kind !== nothing
            @test neo_kind(f) == kind
        end
        # An element after the last one is "not found" in both engines (for a scalar column
        # CSPICE does not look at the element number).
        @test ekgi(res, 1, 1, 2) == c_ekgi(1, 1, 2)
        for row in 1:n
            k = eknelt(res, 5, row)
            @test ekgi(res, 5, row, k + 1)[2:3] == c_ekgi(5, row, k + 1)[2:3]
        end
        # A time column is read with ekgd.
        @test ekgd(res, 4, 1) == c_ekgd(4, 1)
        # The whole column in one call.
        names = ekcolumn(res, 2)
        @test length(names) == n
        @test all(i -> (v = ekgc(res, 2, i); v[2] ? ismissing(names[i]) : names[i] == v[1]), 1:n)
        ivec = ekcolumn(res, 5)
        @test all(i -> ismissing(ivec[i]) ? ekgi(res, 5, i)[2] : length(ivec[i]) == eknelt(res, 5, i), 1:n)
        @test neo_kind(() -> ekcolumn(res, 7)) == :invalidindex
        # A query with an error gives an empty result and does not throw. (CSPICE can stop
        # the program for a column that is not scalar in a WHERE or ORDER BY clause, when
        # it builds the message; these queries are not sent to it.)
        for q2 in ("SELECT ID FROM EVENTS WHERE VEC < 3.0", "SELECT ID FROM EVENTS ORDER BY VEC",
            "SELECT ID FROM EVENTS WHERE ID < 3 AND TAGS = 'a'", "")
            r2, err2, msg2 = ekfind(ctx, q2)
            @test err2 && !isempty(msg2) && r2.nmrows == 0
        end
        bad, err, msg = ekfind(ctx, "SELECT FOO FROM EVENTS")
        @test err && !isempty(msg) && bad.nmrows == 0
        @test neo_kind(() -> ekgi(bad, 1, 1)) == :invalidindex
    end

    ek_part("a result stays valid after its files are unloaded; unload and reload") do
        q = "SELECT E.ID, G.GNAME, L.TAG FROM EVENTS E, LINKS L, GROUPS G WHERE E.ID = L.EID AND L.GID = G.GID ORDER BY E.ID"
        res, _, _ = ekfind(ctx, q)
        before = neo_ek_result(res)
        unload(ctx, paths[2])
        SPICE.unload(paths[2])
        @test neo_ek_result(res) == before
        @test ekntab(ctx) == c_ekntab() == 2
        @test [ektnam(ctx, n) for n in 1:2] == [c_ektnam(n) for n in 1:2]
        r = ek_compare_query(ctx, q)
        @test r.cerr && r.nerr && r.msgsame
        for q2 in ("SELECT ID, NAME FROM EVENTS WHERE ID < 50 ORDER BY NAME",
            "SELECT E.ID, G.GID FROM EVENTS E, GROUPS G WHERE E.GRP = G.GID")
            @test ek_compare_query(ctx, q2).compat
        end
        # Load again: the file is now the last one, so its segments are searched last.
        furnsh(ctx, paths[2])
        SPICE.furnsh(paths[2])
        @test ekntab(ctx) == c_ekntab() == 3
        for q2 in (q, "SELECT ID, NAME FROM EVENTS WHERE ID < 50", "SELECT EID FROM LINKS")
            @test ek_compare_query(ctx, q2).compat
        end
        # The same file twice in the load list holds its tables once.
        furnsh(ctx, paths[3])
        @test ekntab(ctx) == 3
        res2, _, _ = ekfind(ctx, "SELECT GID FROM GROUPS")
        @test res2.nmrows == length(groups.rows)
        unload(ctx, paths[3])
        @test ekfind(ctx, "SELECT GID FROM GROUPS")[1].nmrows == length(groups.rows)
        kclear(ctx)
        @test ekntab(ctx) == 0
        @test ekfind(ctx, "SELECT GID FROM GROUPS")[2]
    end

    ek_part("load errors: tables that do not fit") do
        d2 = mktempdir()
        rng2 = MersenneTwister(SEED + 3)
        base = c_ek_write(joinpath(d2, "base.bes"), [ek_random_segment(rng2, "LINKS", EK_LINK_COLS, 5)])
        fewer = c_ek_write(joinpath(d2, "fewer.bes"), [ek_random_segment(rng2, "LINKS", EK_LINK_COLS[1:4], 5)])
        other = copy(EK_LINK_COLS)
        other[2] = EKColSpec("OTHER", :INT)
        renamed = c_ek_write(joinpath(d2, "renamed.bes"), [ek_random_segment(rng2, "LINKS", other, 5)])
        other = copy(EK_LINK_COLS)
        other[2] = EKColSpec("GID", :DP)
        retyped = c_ek_write(joinpath(d2, "retyped.bes"), [ek_random_segment(rng2, "LINKS", other, 5)])
        for bad in (fewer, renamed, retyped)
            SPICE.kclear()
            SPICE.furnsh(base)
            ccall((:furnsh_c, EKLIB), Cvoid, (Cstring,), bad)
            kind = c_pending_kind()
            c = neo_context([base])
            @test kind !== nothing
            @test neo_kind(() -> furnsh(c, bad)) == kind
            @test ekntab(c) == 1 && ktotal(c, "EK") == 1
        end
        SPICE.kclear()
    end
end

Base.Experimental.@optlevel -1
