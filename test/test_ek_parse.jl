# Group: the EK query parser (scanner, parser, normalisation of the WHERE clause, name
# resolution, time resolution, semantic check) and the SELECT clause analysis of `ekpsel`.
#
# Feature file: features/ek_query_language.feature
#
# CSPICE is the oracle, through direct calls: `zzekencd_` gives the encoded query (the
# integer cell, the character buffer and the double buffer), the accessors `zzekqtab_`,
# `zzekqcnj_`, `zzekqcon_`, `zzekqord_` and `zzekqsel_` decode it, and `ekfind_c` and
# `ekpsel_c` are the public calls. The group writes its own EK file with the CSPICE
# writers (three small tables with scalar and array columns of every data type), loads
# it in CSPICE with `eklef_c`, and builds the table catalogue of NeoSpice from the
# catalogue calls of CSPICE, so it needs no EK reader of NeoSpice.
#
# For each query of the corpus the group compares the error flag and the error message
# and, for a valid query, the three buffers of the encoded query element by element,
# then the decoded query against `_EKQuery`, and `ekpsel_c` against `_ek_psel`. The
# counts are printed. Queries for which this port differs from CSPICE on purpose (see
# the docstring of `_ek_encode_query`) are in a list of their own with the reason.

isdefined(@__MODULE__, :DSKLIB) || include("dskfixtures.jl")
isdefined(@__MODULE__, :CKFix) || include("ckfixtures.jl")
isdefined(NeoSpice, :_ek_encode_query) ||
    Base.include(NeoSpice, joinpath(@__DIR__, "..", "src", "ek", "query_parse.jl"))

# --- the EK file of the group ----------------------------------------------------------

"Tables of the group: name => columns as (name, declaration, type, size, string length)."
const EKQ_TABLES = [
    "ALPHA" => [
        ("ICOL", "DATATYPE = INTEGER, INDEXED = TRUE, NULLS_OK = TRUE", :INT, 1, 0),
        ("DCOL", "DATATYPE = DOUBLE PRECISION, NULLS_OK = TRUE", :DP, 1, 0),
        ("CCOL", "DATATYPE = CHARACTER*(10), INDEXED = TRUE", :CHR, 1, 10),
        ("VCOL", "DATATYPE = CHARACTER*(*)", :CHR, 1, 30),
        ("TCOL", "DATATYPE = TIME, INDEXED = TRUE", :TIME, 1, 0),
        ("IARR", "DATATYPE = INTEGER, SIZE = 3", :INT, 3, 0),
        ("DARR", "DATATYPE = DOUBLE PRECISION, SIZE = VARIABLE", :DP, -1, 0),
        ("CARR", "DATATYPE = CHARACTER*(8), SIZE = 2", :CHR, 2, 8),
        ("NAME", "DATATYPE = CHARACTER*(20)", :CHR, 1, 20)],
    "BETA" => [
        ("ICOL", "DATATYPE = INTEGER", :INT, 1, 0),
        ("XCOL", "DATATYPE = DOUBLE PRECISION, INDEXED = TRUE", :DP, 1, 0),
        ("TCOL", "DATATYPE = TIME", :TIME, 1, 0),
        ("NAME", "DATATYPE = CHARACTER*(20)", :CHR, 1, 20),
        ("TARR", "DATATYPE = TIME, SIZE = 2", :TIME, 2, 0),
        ("BCHR", "DATATYPE = CHARACTER*(5), NULLS_OK = TRUE", :CHR, 1, 5)],
    "GAMMA_LONG_NAME" => [
        ("KEY", "DATATYPE = INTEGER, INDEXED = TRUE", :INT, 1, 0),
        ("VAL", "DATATYPE = DOUBLE PRECISION", :DP, 1, 0),
        ("LABEL", "DATATYPE = CHARACTER*(16)", :CHR, 1, 16),
        ("EPOCH", "DATATYPE = TIME", :TIME, 1, 0),
        ("ICOL", "DATATYPE = INTEGER", :INT, 1, 0)]]

"`strs` as the fixed-length, null-terminated string array that the CSPICE wrappers take."
function ekq_cstrings(strs, len::Int)
    buf = zeros(UInt8, len * length(strs))
    for (k, s) in enumerate(strs)
        @assert sizeof(s) < len
        copyto!(buf, (k - 1) * len + 1, codeunits(s), 1, sizeof(s))
    end
    return buf
end

"Write the EK file of the group at `path` with the CSPICE writers: three rows in each table."
function ekq_write_ek(path)
    h = Ref{Cint}(0)
    ccall((:ekopn_c, DSKLIB), Cvoid, (Cstring, Cstring, Cint, Ref{Cint}), path, "EKQ", 0, h)
    c_check()
    for (table, cols) in EKQ_TABLES
        segno = Ref{Cint}(0)
        ccall((:ekbseg_c, DSKLIB), Cvoid, (Cint, Cstring, Cint, Cint, Ptr{UInt8}, Cint,
                Ptr{UInt8}, Ref{Cint}), h[], table, length(cols), 33,
            ekq_cstrings(first.(cols), 33), 100, ekq_cstrings(getindex.(cols, 2), 100), segno)
        c_check()
        for r in 1:3
            recno = Ref{Cint}(0)
            ccall((:ekappr_c, DSKLIB), Cvoid, (Cint, Cint, Ref{Cint}), h[], segno[], recno)
            c_check()
            for (name, _, type, size, strlen) in cols
                n = size == -1 ? r : size
                if type === :CHR
                    vals = ekq_cstrings(["S$(r)$(k)" for k in 1:n], strlen + 1)
                    ccall((:ekacec_c, DSKLIB), Cvoid, (Cint, Cint, Cint, Cstring, Cint, Cint,
                            Ptr{UInt8}, Cint), h[], segno[], recno[], name, n, strlen + 1, vals, 0)
                elseif type === :INT
                    ccall((:ekacei_c, DSKLIB), Cvoid, (Cint, Cint, Cint, Cstring, Cint,
                            Ptr{Cint}, Cint), h[], segno[], recno[], name, n,
                        Cint[10r + k for k in 1:n], 0)
                else
                    ccall((:ekaced_c, DSKLIB), Cvoid, (Cint, Cint, Cint, Cstring, Cint,
                            Ptr{Float64}, Cint), h[], segno[], recno[], name, n,
                        Float64[1000.0r + k / 2 for k in 1:n], 0)
                end
                c_check()
            end
        end
    end
    ccall((:ekcls_c, DSKLIB), Cvoid, (Cint,), h[])
    c_check()
    return path
end

"The table catalogue of NeoSpice, built from the catalogue calls of CSPICE."
function ekq_catalogue()
    n = Ref{Cint}(0)
    ccall((:ekntab_c, DSKLIB), Cvoid, (Ref{Cint},), n)
    c_check()
    tables = NeoSpice.EKTable[]
    for t in 0:n[]-1
        buf = zeros(UInt8, 65)
        ccall((:ektnam_c, DSKLIB), Cvoid, (Cint, Cint, Ptr{UInt8}), t, 65, buf)
        c_check()
        table = unsafe_string(pointer(buf))
        cc = Ref{Cint}(0)
        ccall((:ekccnt_c, DSKLIB), Cvoid, (Cstring, Ref{Cint}), table, cc)
        c_check()
        cnames = String[]
        attrs = NeoSpice.EKColumnAttr[]
        for c in 0:cc[]-1
            cbuf = zeros(UInt8, 33)
            att = zeros(Cint, 6)        # SpiceEKAttDsc: class, type, strlen, size, indexd, nullok
            ccall((:ekcii_c, DSKLIB), Cvoid, (Cstring, Cint, Cint, Ptr{UInt8}, Ptr{Cint}),
                table, c, 33, cbuf, att)
            c_check()
            push!(cnames, unsafe_string(pointer(cbuf)))
            push!(attrs, NeoSpice.EKColumnAttr(att[1], (:CHR, :DP, :INT, :TIME)[att[2]+1],
                att[3], att[4], att[5] != 0, att[6] != 0))
        end
        push!(tables, NeoSpice.EKTable(table, cnames, attrs, NTuple{2,Int32}[]))
    end
    return tables
end

# --- the oracle ------------------------------------------------------------------------

const EKQ_MSGLEN = 4096
ekq_fstr(buf) = String(rstrip(String(copy(buf)), ' '))

"""
`zzekencd_` of CSPICE for `query`: `(error, errmsg, kind, i, c, d)` with the elements of
the integer cell, the used part of the character buffer and of the double buffer, and
the kind of the SPICE error that CSPICE signalled (or `nothing`).
"""
function ekq_c_encd(query::String)
    eqryi = zeros(Cint, 27869 + 6)
    eqryc = fill(UInt8(' '), 2000)
    eqryd = zeros(Float64, 100)
    error = Ref{Cint}(0)
    errmsg = fill(UInt8(' '), EKQ_MSGLEN)
    errptr = Ref{Cint}(0)
    ccall((:zzekencd_, DSKLIB), Cint, (Ptr{UInt8}, Ptr{Cint}, Ptr{UInt8}, Ptr{Float64},
            Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Cint, Cint, Cint), query, eqryi, eqryc, eqryd,
        error, errmsg, errptr, sizeof(query), 2000, EKQ_MSGLEN)
    kind = c_pending_kind()
    card = clamp(Int(eqryi[6]), 0, 27869)
    return (error=error[] != 0, errmsg=ekq_fstr(errmsg), kind=kind, i=eqryi[7:6+card],
        c=eqryc, d=eqryd, raw=eqryi)
end

"`ekfind_c` of CSPICE for `query`: `(error, errmsg, kind)`."
function ekq_c_find(query::String)
    nmrows = Ref{Cint}(0)
    error = Ref{Cint}(0)
    errmsg = zeros(UInt8, EKQ_MSGLEN)
    ccall((:ekfind_c, DSKLIB), Cvoid, (Cstring, Cint, Ref{Cint}, Ref{Cint}, Ptr{UInt8}),
        query, EKQ_MSGLEN, nmrows, error, errmsg)
    kind = c_pending_kind()
    return (error=error[] != 0, errmsg=unsafe_string(pointer(errmsg)), kind=kind)
end

"`ekpsel_c` of CSPICE for `query`, with 1-based positions and symbols as `_ek_psel`."
function ekq_c_psel(query::String)
    n = Ref{Cint}(0)
    xbegs = zeros(Cint, 50)
    xends = zeros(Cint, 50)
    xtypes = zeros(Cint, 50)
    xclass = zeros(Cint, 50)
    tabs = zeros(UInt8, 50 * 65)
    cols = zeros(UInt8, 50 * 33)
    error = Ref{Cint}(0)
    errmsg = zeros(UInt8, EKQ_MSGLEN)
    ccall((:ekpsel_c, DSKLIB), Cvoid, (Cstring, Cint, Cint, Cint, Ref{Cint}, Ptr{Cint},
            Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{UInt8}, Ptr{UInt8}, Ref{Cint}, Ptr{UInt8}),
        query, EKQ_MSGLEN, 65, 33, n, xbegs, xends, xtypes, xclass, tabs, cols, error, errmsg)
    kind = c_pending_kind()
    (error[] != 0 || kind !== nothing) && return (n=0, xbegs=Int[], xends=Int[],
        xtypes=Symbol[], xclass=Symbol[], tabs=String[], cols=String[], error=error[] != 0,
        errmsg=unsafe_string(pointer(errmsg)), kind=kind)
    m = Int(n[])
    return (n=m, xbegs=Int.(xbegs[1:m]) .+ 1, xends=Int.(xends[1:m]) .+ 1,
        xtypes=[(:CHR, :DP, :INT, :TIME)[xtypes[k]+1] for k in 1:m],
        xclass=[(:COL, :FUNC, :EXPR)[xclass[k]+1] for k in 1:m],
        tabs=[unsafe_string(pointer(tabs, (k - 1) * 65 + 1)) for k in 1:m],
        cols=[unsafe_string(pointer(cols, (k - 1) * 33 + 1)) for k in 1:m], error=false,
        errmsg="", kind=nothing)
end

"Decode the encoded query `e` of `ekq_c_encd` with the accessors of CSPICE, as an `_EKQuery` without catalogue indices."
function ekq_c_decode(e)
    Q = NeoSpice
    ntab, ncns, ncnj, nord, nsel = (Int(e.i[k]) for k in (8, 9, 10, 11, 12))
    fs() = fill(UInt8(' '), 64)
    tables = String[]
    aliases = String[]
    for n in 1:ntab
        t, a = fs(), fs()
        ccall((:zzekqtab_, DSKLIB), Cint, (Ptr{Cint}, Ptr{UInt8}, Ref{Cint}, Ptr{UInt8},
                Ptr{UInt8}, Cint, Cint, Cint), e.raw, e.c, n, t, a, 2000, 64, 64)
        c_check()
        push!(tables, ekq_fstr(t))
        push!(aliases, ekq_fstr(a))
    end
    conjsizes = Int[]
    for n in 1:ncnj
        size = Ref{Cint}(0)
        ccall((:zzekqcnj_, DSKLIB), Cint, (Ptr{Cint}, Ref{Cint}, Ref{Cint}), e.raw, n, size)
        c_check()
        push!(conjsizes, size[])
    end
    constraints = Q._EKQConstraint[]
    for n in 1:ncns
        lt, lc, rt, rc = fs(), fs(), fs(), fs()
        v = [Ref{Cint}(0) for _ in 1:10]    # cnstyp ltidx lcidx opcode rtidx rcidx dtype cbeg cend ival
        dval = Ref{Float64}(0.0)
        ccall((:zzekqcon_, DSKLIB), Cint, (Ptr{Cint}, Ptr{UInt8}, Ptr{Float64}, Ref{Cint},
                Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Ref{Cint}, Ptr{UInt8},
                Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Float64},
                Ref{Cint}, Cint, Cint, Cint, Cint, Cint), e.raw, e.c, e.d, n, v[1], lt, v[2], lc,
            v[3], v[4], rt, v[5], rc, v[6], v[7], v[8], v[9], dval, v[10], 2000, 64, 64, 64, 64)
        c_check()
        cnstyp, opcode = Int(v[1][]), Int(v[4][])
        iscol = cnstyp == 1
        isval = !iscol && opcode != 9 && opcode != 10
        dtype = isval ? Int(v[7][]) : 0
        push!(constraints, Q._EKQConstraint(cnstyp,
            Q._EKQCol(ekq_fstr(lt), v[2][], ekq_fstr(lc), v[3][], 0, 0), opcode,
            iscol ? Q._EKQCol(ekq_fstr(rt), v[5][], ekq_fstr(rc), v[6][], 0, 0) : Q._EKQ_NOCOL,
            dtype, dtype == 1 ? String(e.c[v[8][]:v[9][]]) : "", isval ? dval[] : 0.0,
            isval ? Int(v[10][]) : 0))
    end
    orderby = Tuple{Q._EKQCol,Int}[]
    for n in 1:nord
        t, c = fs(), fs()
        tabidx, colidx, sense = Ref{Cint}(0), Ref{Cint}(0), Ref{Cint}(0)
        ccall((:zzekqord_, DSKLIB), Cint, (Ptr{Cint}, Ptr{UInt8}, Ref{Cint}, Ptr{UInt8},
                Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Ref{Cint}, Cint, Cint, Cint), e.raw, e.c, n, t,
            tabidx, c, colidx, sense, 2000, 64, 64)
        c_check()
        push!(orderby, (Q._EKQCol(ekq_fstr(t), tabidx[], ekq_fstr(c), colidx[], 0, 0), sense[]))
    end
    select = Q._EKQCol[]
    for n in 1:nsel
        t, c = fs(), fs()
        tabidx, colidx, lxb, lxe = Ref{Cint}(0), Ref{Cint}(0), Ref{Cint}(0), Ref{Cint}(0)
        ccall((:zzekqsel_, DSKLIB), Cint, (Ptr{Cint}, Ptr{UInt8}, Ref{Cint}, Ref{Cint},
                Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Cint, Cint, Cint), e.raw,
            e.c, n, lxb, lxe, t, tabidx, c, colidx, 2000, 64, 64)
        c_check()
        push!(select, Q._EKQCol(ekq_fstr(t), tabidx[], ekq_fstr(c), colidx[], lxb[], lxe[]))
    end
    return Q._EKQuery(tables, aliases, Int[], conjsizes, constraints, orderby, select)
end

"True if the column references are equal; the lexeme bounds count only if `lex`."
ekq_samecol(a, b, lex) = a.table == b.table && a.tabidx == b.tabidx && a.name == b.name &&
                         a.colidx == b.colidx && (!lex || (a.lxbeg == b.lxbeg && a.lxend == b.lxend))

"True if the query `q` of NeoSpice equals the decoded query `c` of CSPICE. `zzekqcon_` gives the double of types 2 and 4 and the integer of type 3."
function ekq_samequery(q, c)
    (q.tables == c.tables && q.aliases == c.aliases && q.conjsizes == c.conjsizes &&
     length(q.constraints) == length(c.constraints) && length(q.orderby) == length(c.orderby) &&
     length(q.select) == length(c.select)) || return false
    for (a, b) in zip(q.constraints, c.constraints)
        (a.cnstyp == b.cnstyp && a.opcode == b.opcode && a.dtype == b.dtype &&
         ekq_samecol(a.lhs, b.lhs, false) && ekq_samecol(a.rhs, b.rhs, false) &&
         a.cval == b.cval && a.ival == b.ival) || return false
        (a.dtype == 2 || a.dtype == 4) && (a.dval === b.dval || return false)
    end
    all(ekq_samecol(a[1], b[1], false) && a[2] == b[2] for (a, b) in zip(q.orderby, c.orderby)) ||
        return false
    return all(ekq_samecol(a, b, true) for (a, b) in zip(q.select, c.select))
end

"`_ek_qencode` of NeoSpice with a thrown `SpiceError` as its kind."
function ekq_n_encd(fx, query)
    try
        q, tabcat, error, errmsg = NeoSpice._ek_qencode(fx.tables, fx.ks, fx.td, query)
        return (q=q, tabcat=tabcat, error=error, errmsg=errmsg, kind=nothing)
    catch e
        e isa NeoSpice.SpiceError || rethrow()
        return (q=nothing, tabcat=Int[], error=false, errmsg="", kind=e.kind)
    end
end

"""
Message of NeoSpice for a non-scalar column where CSPICE has not set the end of the text
that it prints (the local variable `LXE(2)` of `zzeksemc_`): CSPICE then reads memory
at an undefined offset and can end the process, so the oracle is not called.
"""
const EKQ_UNSAFE_NONSCALAR = r"^Non-scalar column <[A-Za-z][A-Za-z0-9_$]*> having size"

"A query that ends with NOT: if the NOT is in the WHERE clause, CSPICE reads the token after the last one, which is what an earlier call left there. The random token strings are not sent to CSPICE in this case."
const EKQ_UNSAFE_NOT = r"(^|[^A-Za-z0-9_$])[Nn][Oo][Tt]\s*$"

"""
Compare the two engines on `query`. Returns `(status, detail)`: `:valid` or `:invalid`
if all results are equal (the query is valid or not), `:signal` if both signal the same
SPICE error, and a status that starts with `diff` otherwise.
"""
function ekq_compare(fx, query::String; guard_not::Bool=false)
    n = ekq_n_encd(fx, query)
    n.kind === :toomanyconjunctions && return :unsafe, "CSPICE stops the program"
    sizeof(rstrip(query, ' ')) > 2000 &&
        return :unsafe, "CSPICE reads past its buffer of 2000 characters"
    (guard_not && occursin(EKQ_UNSAFE_NOT, query)) &&
        return :unsafe, "CSPICE reads the token after the last one"
    (n.error && occursin(EKQ_UNSAFE_NONSCALAR, n.errmsg)) &&
        return :unsafe, "CSPICE reads an undefined position: " * n.errmsg
    c = ekq_c_encd(query)
    if n.kind !== nothing || c.kind !== nothing
        n.kind === c.kind && return :signal, string(n.kind)
        return :diff_signal, "NeoSpice $(n.kind), CSPICE $(c.kind) [$(c.errmsg)]"
    end
    n.error == c.error || return :diff_flag, "NeoSpice $(n.error) [$(n.errmsg)], CSPICE " *
                                             "$(c.error) [$(c.errmsg)]"
    # The semantic check stops at IS NULL. The search of ekfind_c then meets constraints
    # that were not checked (it signals SPICE(BUG) for a non-scalar column), so ekfind_c
    # is called only for a query that is not valid or that was checked to the end.
    f = (n.error || n.q.i[7] == 1) ? ekq_c_find(query) : (error=false, errmsg="", kind=nothing)
    (f.error == c.error && f.errmsg == c.errmsg && f.kind === nothing) ||
        return :diff_ekfind, "ekfind_c: $(f.error) $(f.kind) [$(f.errmsg)]; zzekencd_: [$(c.errmsg)]"
    p = ekq_c_psel(query)
    s = NeoSpice._ek_psel(fx.tables, fx.ks, fx.td, query)
    if n.error
        n.errmsg == c.errmsg || return :diff_msg, "NeoSpice [$(n.errmsg)], CSPICE [$(c.errmsg)]"
        (s.error && p.error && s.errmsg == p.errmsg && s.n == 0) ||
            return :diff_psel, "error case: NeoSpice [$(s.errmsg)], CSPICE [$(p.errmsg)]"
        return :invalid, n.errmsg
    end
    n.q.i == c.i || return :diff_cell, "integer cell differs at " *
                                       "$(findfirst(k -> k > min(length(n.q.i), length(c.i)) || n.q.i[k] != c.i[k], 1:max(length(n.q.i), length(c.i))))"
    nc, nd = Int(c.i[16]) - 1, Int(c.i[14]) - 1
    n.q.c[1:nc] == c.c[1:nc] || return :diff_chars, "character buffer differs"
    all(n.q.d[k] === c.d[k] for k in 1:nd) || return :diff_doubles, "double buffer differs: " *
                                                                    "$(n.q.d[1:nd]) $(c.d[1:nd])"
    q = NeoSpice._ek_qdecode(n.q, n.tabcat)
    ekq_samequery(q, ekq_c_decode(c)) || return :diff_decode, "decoded query differs"
    all(fx.tables[q.tabcat[k]].name == q.tables[k] for k in eachindex(q.tables)) ||
        return :diff_decode, "catalogue indices"
    sum(q.conjsizes; init=0) == length(q.constraints) || return :diff_decode, "invariant"
    (!s.error && !p.error && s.n == p.n && s.xbegs == p.xbegs && s.xends == p.xends &&
     s.xtypes == p.xtypes && s.xclass == p.xclass && s.tabs == p.tabs && s.cols == p.cols) ||
        return :diff_psel, "NeoSpice $(s), CSPICE $(p)"
    return :valid, ""
end

# --- the corpus ------------------------------------------------------------------------

"Valid queries, written by hand: one or more for each rule of the language."
function ekq_valid_corpus()
    q = String[
        # SELECT and FROM
        "SELECT ICOL FROM ALPHA",
        "select icol from alpha",
        "Select Icol , Dcol From Alpha",
        "   SELECT   ICOL   FROM   ALPHA   ",
        "SELECT ICOL,DCOL,CCOL,VCOL,TCOL,IARR,DARR,CARR,NAME FROM ALPHA",
        "SELECT ALPHA.ICOL FROM ALPHA",
        "SELECT ALPHA . ICOL FROM ALPHA",
        "SELECT A.ICOL, A.DCOL FROM ALPHA A",
        "SELECT a.icol FROM alpha a",
        "SELECT ICOL FROM ALPHA ALPHA",
        "SELECT ALPHA.ICOL FROM ALPHA ALPHA",
        "SELECT IARR, DARR, CARR FROM ALPHA",
        "SELECT TARR FROM BETA",
        "SELECT KEY, VAL, LABEL, EPOCH FROM GAMMA_LONG_NAME",
        "SELECT GAMMA_LONG_NAME.KEY FROM GAMMA_LONG_NAME",
        "SELECT G\$_1.KEY FROM GAMMA_LONG_NAME G\$_1",
        "SELECT DCOL FROM ALPHA, BETA B WHERE DCOL < XCOL ORDER BY XCOL SELECT XCOL",
        "SELECT ALPHA.ICOL FROM ALPHA A",
        # the order of the clauses is free
        "FROM ALPHA SELECT ICOL",
        "FROM ALPHA WHERE ICOL > 1 SELECT ICOL",
        "FROM ALPHA WHERE ICOL > 1 ORDER BY ICOL SELECT ICOL",
        "FROM ALPHA ORDER BY ICOL SELECT ICOL",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL WHERE ICOL > 1",
        "WHERE ICOL > 1 FROM ALPHA SELECT ICOL",
        "WHERE ICOL > 1 SELECT ICOL FROM ALPHA",
        "ORDER BY ICOL DESC SELECT ICOL FROM ALPHA WHERE ICOL > 1",
        "ORDER BY ICOL SELECT ICOL FROM ALPHA",
        "SELECT ICOL WHERE DCOL < 2 FROM ALPHA",
        # comparison operators
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL != 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL <> 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL ^= 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL < 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL <= 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL > 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL >= 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL EQ 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL NE 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL LT 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL LE 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL GT 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL GE 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL ge 1 and icol lt 30",
        "SELECT ICOL FROM ALPHA WHERE ICOL<1",
        "SELECT ICOL FROM ALPHA WHERE ICOL>=1AND DCOL<=2.5",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS NOT 1",
        # numbers
        "SELECT ICOL FROM ALPHA WHERE ICOL = -1",
        "SELECT ICOL FROM ALPHA WHERE ICOL = +1",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 0005",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 12345678901",
        "SELECT ICOL FROM ALPHA WHERE ICOL = -99999999999999999999",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 2147483647",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 2147483648",
        "SELECT ICOL FROM ALPHA WHERE ICOL = -2147483648",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1.5",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 2.5",
        "SELECT ICOL FROM ALPHA WHERE ICOL = -2.5",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1.",
        "SELECT ICOL FROM ALPHA WHERE DCOL = .5",
        "SELECT ICOL FROM ALPHA WHERE DCOL = -.5",
        "SELECT ICOL FROM ALPHA WHERE DCOL = +.5",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 1e3",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 1E3",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 1d3",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 1D-3",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 1.e2",
        "SELECT ICOL FROM ALPHA WHERE DCOL = -1.25E+2",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 3.14159265358979323846264338327950288",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 123456789012345678901234567890",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 0.000000000000000000000000000001234",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 1e-400",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 1.7976931348623157e308",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 7",
        "SELECT ICOL FROM ALPHA WHERE DCOL > ICOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL > DCOL",
        # strings
        "SELECT ICOL FROM ALPHA WHERE CCOL = 'ABC'",
        "SELECT ICOL FROM ALPHA WHERE CCOL = \"ABC\"",
        "SELECT ICOL FROM ALPHA WHERE CCOL = 'abc def'",
        "SELECT ICOL FROM ALPHA WHERE CCOL = '  lead and trail  '",
        "SELECT ICOL FROM ALPHA WHERE CCOL = 'it''s'",
        "SELECT ICOL FROM ALPHA WHERE CCOL = ''''",
        "SELECT ICOL FROM ALPHA WHERE CCOL = ''''''",
        "SELECT ICOL FROM ALPHA WHERE CCOL = \"say \"\"hi\"\"\"",
        "SELECT ICOL FROM ALPHA WHERE CCOL = \"it's\"",
        "SELECT ICOL FROM ALPHA WHERE CCOL = 'a \"quoted\" word'",
        "SELECT ICOL FROM ALPHA WHERE CCOL = ' '",
        "SELECT ICOL FROM ALPHA WHERE CCOL = 'SELECT FROM WHERE'",
        "SELECT ICOL FROM ALPHA WHERE CCOL = 'a#b' AND VCOL = '#'",
        "SELECT ICOL FROM ALPHA WHERE CCOL < 'M' AND VCOL >= 'B'",
        "SELECT ICOL FROM ALPHA WHERE CCOL = NAME",
        "SELECT ICOL FROM ALPHA WHERE CCOL = 'ABC'AND VCOL='x'",
        # LIKE
        "SELECT ICOL FROM ALPHA WHERE CCOL LIKE 'A*'",
        "SELECT ICOL FROM ALPHA WHERE CCOL NOT LIKE 'A%'",
        "SELECT ICOL FROM ALPHA WHERE CCOL like '*x*' and vcol not like \"%\"",
        "SELECT ICOL FROM ALPHA WHERE NOT CCOL LIKE 'A*'",
        "SELECT ICOL FROM ALPHA WHERE NOT CCOL NOT LIKE 'A*'",
        "SELECT ICOL FROM ALPHA WHERE CCOL LIKE NAME",
        # NULL
        "SELECT ICOL FROM ALPHA WHERE ICOL IS NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS NOT NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL = NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL != NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL EQ NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL NE NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL <> (NULL)",
        "SELECT ICOL FROM ALPHA WHERE NOT ICOL IS NULL",
        "SELECT ICOL FROM ALPHA WHERE NOT (ICOL IS NOT NULL)",
        "SELECT ICOL FROM ALPHA WHERE DCOL > 1 AND ICOL IS NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS NULL AND CCOL = 5",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS NULL ORDER BY IARR",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS NULL OR IARR = 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS NULL AND IARR = 1 OR ICOL = 'x'",
        "SELECT ICOL FROM BETA WHERE BCHR IS NULL OR BCHR LIKE 'S*'",
        # BETWEEN
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1 AND 10",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 10 AND 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL NOT BETWEEN 1 AND 10",
        "SELECT ICOL FROM ALPHA WHERE DCOL BETWEEN -1.5 AND 1e3",
        "SELECT ICOL FROM ALPHA WHERE CCOL BETWEEN 'A' AND 'M'",
        "SELECT ICOL FROM ALPHA WHERE CCOL NOT BETWEEN 'A' AND 'M'",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN (1) AND (10)",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN DCOL AND 10",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1 AND DCOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN ALPHA.DCOL AND ALPHA.ICOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN DCOL AND ALPHA.ICOL",
        "SELECT ICOL FROM ALPHA WHERE NOT ICOL BETWEEN 1 AND 10",
        "SELECT ICOL FROM ALPHA WHERE NOT ICOL NOT BETWEEN 1 AND 10",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1 AND 10 AND DCOL BETWEEN 2 AND 3",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1 AND 10 OR DCOL NOT BETWEEN 2 AND 3",
        "SELECT A.ICOL FROM ALPHA A, BETA B WHERE A.ICOL BETWEEN B.ICOL AND B.XCOL",
        # logic
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 OR DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 OR DCOL = 2 AND CCOL = 'X'",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND DCOL = 2 OR CCOL = 'X'",
        "SELECT ICOL FROM ALPHA WHERE (ICOL = 1 OR DCOL = 2) AND CCOL = 'X'",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND (DCOL = 2 OR CCOL = 'X')",
        "SELECT ICOL FROM ALPHA WHERE (ICOL = 1 OR DCOL = 2) AND (CCOL = 'X' OR VCOL = 'Y')",
        "SELECT ICOL FROM ALPHA WHERE NOT ICOL = 1",
        "SELECT ICOL FROM ALPHA WHERE NOT NOT ICOL = 1",
        "SELECT ICOL FROM ALPHA WHERE NOT NOT NOT ICOL < 1",
        "SELECT ICOL FROM ALPHA WHERE NOT (ICOL = 1 AND DCOL = 2)",
        "SELECT ICOL FROM ALPHA WHERE NOT (ICOL = 1 OR DCOL = 2)",
        "SELECT ICOL FROM ALPHA WHERE NOT (ICOL = 1 OR DCOL = 2) AND NOT CCOL = 'X'",
        "SELECT ICOL FROM ALPHA WHERE NOT ((ICOL <= 1 OR DCOL > 2) AND (CCOL LIKE 'X' OR VCOL >= 'Y'))",
        "SELECT ICOL FROM ALPHA WHERE NOT ICOL = 1 AND DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE NOT ICOL = 1 OR NOT DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE ((ICOL = 1))",
        "SELECT ICOL FROM ALPHA WHERE (((ICOL = 1) AND ((DCOL = 2))))",
        "SELECT ICOL FROM ALPHA WHERE ICOL = ((1))",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 OR ICOL = 2 OR ICOL = 3 OR ICOL = 4",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND ICOL = 2 AND ICOL = 3 AND ICOL = 4",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 OR ICOL = 2 AND ICOL = 3 OR ICOL = 4 AND ICOL = 5",
        "SELECT ICOL FROM ALPHA WHERE (ICOL=1 OR ICOL=2) AND (ICOL=3 OR ICOL=4) AND (ICOL=5 OR ICOL=6)",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 OR NOT (DCOL = 2 AND NOT (CCOL = 'X' OR VCOL = 'Y'))",
        # joins, aliases, qualified names
        "SELECT A.ICOL, B.ICOL FROM ALPHA A, BETA B WHERE A.ICOL = B.ICOL",
        "SELECT ALPHA.ICOL, BETA.ICOL FROM ALPHA, BETA WHERE ALPHA.ICOL = BETA.ICOL",
        "SELECT DCOL, XCOL FROM ALPHA, BETA WHERE DCOL < XCOL",
        "SELECT A1.ICOL, A2.ICOL FROM ALPHA A1, ALPHA A2 WHERE A1.ICOL = A2.ICOL",
        "SELECT A.ICOL FROM ALPHA A, BETA B, GAMMA_LONG_NAME G WHERE A.ICOL = B.ICOL AND B.ICOL = G.ICOL",
        "SELECT KEY, XCOL, DCOL FROM GAMMA_LONG_NAME, BETA, ALPHA WHERE KEY = DCOL AND XCOL != VAL",
        "SELECT A.NAME FROM ALPHA A, BETA B WHERE A.NAME = B.NAME AND A.TCOL < B.TCOL",
        "SELECT ALPHA.ICOL FROM ALPHA, ALPHA A",
        "SELECT X.DCOL FROM BETA, ALPHA X, GAMMA_LONG_NAME WHERE X.ICOL > BETA.XCOL",
        "SELECT ICOL FROM ALPHA WHERE ALPHA.ICOL = ICOL",
        # ORDER BY
        "SELECT ICOL FROM ALPHA ORDER BY ICOL",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL ASC",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL DESC",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL asc, DCOL desc, CCOL",
        "SELECT ICOL FROM ALPHA ORDER BY ALPHA.ICOL",
        "SELECT ICOL FROM ALPHA ORDER BY ALPHA.ICOL DESC, ALPHA.TCOL ASC",
        "SELECT ICOL FROM ALPHA A ORDER BY A.ICOL, DCOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL > 1 ORDER BY DCOL DESC",
        "SELECT A.ICOL FROM ALPHA A, BETA B WHERE A.ICOL = B.ICOL ORDER BY B.XCOL DESC, A.CCOL",
        "SELECT ICOL FROM ALPHA ORDER BY TCOL, VCOL, NAME, CCOL, DCOL, ICOL, TCOL, VCOL, NAME, CCOL",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL WHERE ICOL > 1",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL DESC WHERE ICOL > 1",
        # times
        "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1 12:00:00'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1 12:00:00 TDB'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '1999-12-31T23:59:59.5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'JD 2451545.0'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > \"jan 5, 1995  11:00 pm\"",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '   1996   mar   3   '",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '2016 DEC 31 23:59:60'",
        "SELECT ICOL FROM ALPHA WHERE TCOL = '1995-123 // 12:00'",
        "SELECT ICOL FROM ALPHA WHERE TCOL BETWEEN '1990 JAN 1' AND '2010 JAN 1'",
        "SELECT ICOL FROM ALPHA WHERE TCOL NOT BETWEEN '1990 JAN 1' AND '2010 JAN 1'",
        "SELECT ICOL FROM ALPHA WHERE TCOL <= '2001 FEB 3 04:05:06.789 UTC' OR TCOL IS NULL",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1 12:00:00.000 PST'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1 12:00:00.000          this is cut off'",
        "SELECT A.ICOL FROM ALPHA A, BETA B WHERE A.TCOL >= B.TCOL",
        "SELECT ICOL FROM GAMMA_LONG_NAME WHERE EPOCH LT '1 JAN 2020' AND EPOCH GE '1 JAN 1980'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 1/1000:5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'ekq_craft  sclk  5000:0'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 2/600000.0'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFTSCLK 77'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_TDT SCLK 1/5:10:10'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_TDT SCLK 2:59.49'",
    ]
    # the limits: 10 tables, 50 SELECT columns, 10 order-by columns, long names
    push!(q, "SELECT T1.ICOL FROM " * join(("ALPHA T$k" for k in 1:10), ", "))
    push!(q, "SELECT " * join(("ICOL" for _ in 1:50), ", ") * " FROM ALPHA")
    push!(q, "SELECT ICOL FROM ALPHA ORDER BY " * join(("ALPHA.ICOL" for _ in 1:10), ", "))
    push!(q, "SELECT " * "A"^64 * ".ICOL FROM ALPHA " * "A"^64)
    push!(q, "SELECT ICOL FROM ALPHA WHERE CCOL = '" * "x"^900 * "' OR VCOL = '" * "y"^900 * "'")
    push!(q, "SELECT ICOL FROM ALPHA WHERE " * join(("ICOL = $k" for k in 1:100), " OR "))
    push!(q, "SELECT ICOL FROM ALPHA WHERE " * join(("ICOL != $k" for k in 1:100), " AND "))
    push!(q, "SELECT ICOL FROM ALPHA WHERE " *
             join(("(ICOL = $k OR DCOL = $k OR ICOL > 9$k)" for k in 1:4), " AND "))
    push!(q, "SELECT ICOL FROM ALPHA WHERE NOT (" *
             join(("(ICOL = $k AND DCOL = $k)" for k in 1:6), " OR ") * ")")
    return q
end

"Queries that are not valid, written by hand: one or more for each error message."
function ekq_invalid_corpus()
    q = String[
        # scanner
        " ",
        "SELECT ICOL FROM ALPHA WHERE ICOL = @",
        "SELECT ICOL FROM ALPHA WHERE ICOL ! 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL ^ 1",
        "SELECT ICOL FROM ALPHA; ",
        "SELECT ICOL\tFROM ALPHA",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 \x01",
        "SELECT ICOL FROM ALPHA WHERE ICOL = \xff",
        "SELECT #ICOL FROM ALPHA",
        "SELECT _ICOL FROM ALPHA",
        "SELECT \$ICOL FROM ALPHA",
        "SELECT ICOL FROM ALPHA WHERE CCOL = 'abc",
        "SELECT ICOL FROM ALPHA WHERE CCOL = \"abc'",
        "SELECT ICOL FROM ALPHA WHERE CCOL = '",
        "SELECT ICOL FROM ALPHA WHERE CCOL = ''",
        "SELECT ICOL FROM ALPHA WHERE CCOL = \"\"",
        "SELECT ICOL FROM ALPHA WHERE CCOL = 'abc''",
        "SELECT ICOL FROM ALPHA WHERE CCOL = '''",
        "SELECT ICOL FROM ALPHA WHERE ICOL = -",
        "SELECT ICOL FROM ALPHA WHERE ICOL = + 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL = -ICOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL = -.",
        "SELECT ICOL FROM ALPHA WHERE ICOL = +.E1",
        "SELECT ICOL FROM ALPHA WHERE ICOL - 1 = 2",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 1e999",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 1e309",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 123456789012345678901234567890e300",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 1.8e308",
        "SELECT ICOL FROM ALPHA WHERE DCOL = -0.0000001E3000",
        # parser
        "SELECT ICOL",
        "ICOL ALPHA",
        "FROM ALPHA",
        "FROM ALPHA WHERE ICOL = 1",
        "FROM ALPHA ORDER BY ICOL",
        "SELECT ICOL FROM",
        "SELECT ICOL FROM 5",
        "SELECT ICOL FROM 'ALPHA'",
        "SELECT ICOL FROM WHERE ICOL = 1",
        "SELECT ICOL FROM ALPHA,",
        "SELECT ICOL FROM ALPHA, , BETA",
        "SELECT ICOL FROM ALPHA A B",
        "SELECT ICOL FROM ALPHA A 5",
        "SELECT ICOL FROM ALPHA.BETA",
        "SELECT ICOL FROM ALPHA 7",
        "SELECT ICOL FROM ALPHA (",
        "SELECT ICOL FROM ALPHA AND BETA",
        "SELECT ICOL FROM ALPHA A NOT",
        "SELECT ICOL FROM ALPHA A FROM BETA",
        "SELECT ICOL FROM ALPHA GROUP BY ICOL",
        "SELECT FROM ALPHA",
        "SELECT * FROM ALPHA",
        "SELECT 5 FROM ALPHA",
        "SELECT 'ICOL' FROM ALPHA",
        "SELECT COUNT(*) FROM ALPHA",
        "SELECT MAX FROM ALPHA",
        "SELECT ICOL, FROM ALPHA",
        "FROM ALPHA SELECT ICOL,",
        "FROM ALPHA SELECT",
        "FROM ALPHA SELECT ICOL.",
        "FROM ALPHA SELECT ICOL. FROM",
        "FROM ALPHA SELECT ALPHA.5",
        "FROM ALPHA SELECT ICOL DCOL",
        "FROM ALPHA SELECT ICOL (",
        "FROM ALPHA SELECT ICOL AND",
        "FROM ALPHA SELECT ALPHA.ICOL DCOL",
        "FROM ALPHA SELECT ALPHA.ICOL AND",
        "FROM ALPHA SELECT ALPHA.ICOL.DCOL",
        "SELECT DISTINCT ICOL FROM ALPHA",
        "SELECT ALL ICOL FROM ALPHA",
        "SELECT ICOL FROM ALPHA ORDER ICOL",
        "SELECT ICOL FROM ALPHA ORDER",
        "SELECT ICOL FROM ALPHA ORDER BY",
        "SELECT ICOL FROM ALPHA ORDER BY 5",
        "SELECT ICOL FROM ALPHA ORDER BY DESC",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL,",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL DCOL",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL 5",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL AND",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL.",
        "SELECT ICOL FROM ALPHA ORDER BY ALPHA.'x'",
        "SELECT ICOL FROM ALPHA ORDER BY ALPHA.ICOL DCOL",
        "SELECT ICOL FROM ALPHA ORDER BY ALPHA.ICOL NOT",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL ASC DESC",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL ASC DCOL",
        "SELECT ICOL FROM ALPHA ORDER BY ICOL DESC BY",
        "ORDER BY ICOL FROM ALPHA",
        # WHERE clause
        "SELECT ICOL FROM ALPHA WHERE",
        "SELECT ICOL FROM ALPHA WHERE ORDER BY ICOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL =",
        "SELECT ICOL FROM ALPHA WHERE = 1 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 GROUP BY ICOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL = MAX",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 , DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL = * AND DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE (ICOL = 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND (",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND ( )",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1)",
        "SELECT ICOL FROM ALPHA WHERE (ICOL = 1))",
        "SELECT ICOL FROM ALPHA WHERE ( ) ICOL = 1",
        "SELECT ICOL FROM ALPHA WHERE (ICOL) = 1",
        "SELECT ICOL FROM ALPHA WHERE 1 = ICOL",
        "SELECT ICOL FROM ALPHA WHERE 1 = 1",
        "SELECT ICOL FROM ALPHA WHERE 1 AND 2",
        "SELECT ICOL FROM ALPHA WHERE 1 OR 2",
        "SELECT ICOL FROM ALPHA WHERE 1 2 3",
        "SELECT ICOL FROM ALPHA WHERE 'a' ICOL DCOL",
        "SELECT ICOL FROM ALPHA WHERE NULL AND 2",
        "SELECT ICOL FROM ALPHA WHERE 1 AND NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN NULL AND 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1 AND NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1 OR 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN AND 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1 AND",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1 AND 2 AND",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1 AND 2 3",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN 1 AND = 3",
        "SELECT ICOL FROM ALPHA WHERE ICOL BETWEEN (1 AND 2)",
        "SELECT ICOL FROM ALPHA WHERE ICOL NOT 1 AND 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL DCOL 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL AND DCOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL = = 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 OR",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND AND DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND OR DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 = 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND 5",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 OR 'x'",
        "SELECT ICOL FROM ALPHA WHERE AND ICOL = 1",
        "SELECT ICOL FROM ALPHA WHERE NOT 5 = ICOL",
        "SELECT ICOL FROM ALPHA WHERE NOT AND ICOL = 1",
        "SELECT ICOL FROM ALPHA WHERE NOT ) ICOL = 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND NOT 5",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND NOT ORDER BY ICOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 OR NOT NOT ORDER BY ICOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 NOT DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE . ICOL = 1",
        "SELECT ICOL FROM ALPHA WHERE ALPHA. = 1",
        "SELECT ICOL FROM ALPHA WHERE ALPHA.ICOL.X = 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL < NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL LIKE NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL NOT LIKE NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS NOT NOT NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS NULL NULL",
        "SELECT ICOL FROM ALPHA WHERE NULL IS NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS",
        "SELECT ICOL FROM ALPHA WHERE ICOL IS NOT",
        "SELECT ICOL FROM ALPHA WHERE ICOL LIKE",
        "SELECT ICOL FROM ALPHA WHERE (ICOL = 1) (DCOL = 2)",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 ( AND DCOL = 2",
        "SELECT ICOL FROM ALPHA WHERE ICOL IN (1, 2)",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 SELECT DCOL",
        # name resolution
        "SELECT A.ICOL FROM ALPHA A, BETA A",
        "SELECT A.ICOL FROM ALPHA A, BETA B, GAMMA_LONG_NAME A",
        "SELECT ALPHA.ICOL FROM ALPHA, BETA ALPHA",
        "SELECT B.ICOL FROM ALPHA BETA, BETA ALPHA",
        "SELECT ICOL FROM NOSUCH",
        "SELECT ICOL FROM ALPHA, NOSUCH N",
        "SELECT ICOL FROM alpha_",
        "SELECT ICOL FROM " * "ALPHA" * "X"^70,
        "SELECT NOSUCH FROM ALPHA",
        "SELECT ICOL FROM ALPHA WHERE NOSUCH = 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL = NOSUCH",
        "SELECT ICOL FROM ALPHA ORDER BY NOSUCH",
        "SELECT XCOL FROM ALPHA",
        "SELECT ICOL FROM ALPHA, BETA",
        "SELECT DCOL FROM ALPHA, BETA WHERE ICOL = 1",
        "SELECT DCOL FROM ALPHA, BETA ORDER BY NAME",
        "SELECT ICOL FROM ALPHA A1, ALPHA A2",
        "SELECT BETA.ICOL FROM ALPHA",
        "SELECT B.ICOL FROM ALPHA BETA, BETA B",
        "SELECT ICOL FROM ALPHA WHERE B.ICOL = 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL = B.ICOL",
        "SELECT ICOL FROM ALPHA ORDER BY B.ICOL",
        "SELECT ALPHA.XCOL FROM ALPHA",
        "SELECT A.XCOL FROM ALPHA A, BETA B",
        "SELECT ICOL FROM ALPHA WHERE ALPHA.NOSUCH = 1",
        "SELECT ICOL FROM ALPHA ORDER BY ALPHA.NOSUCH DESC",
        "SELECT ICOL FROM ALPHA WHERE " * "C"^40 * " = 1",
        "SELECT " * "A"^70 * ".ICOL FROM ALPHA " * "A"^70,
        # time resolution
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'not a time'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1 25:00:00 XYZ'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > ' '",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1 12:00:00.000 this is cut off'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 1",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 1.5e3",
        "SELECT ICOL FROM ALPHA WHERE TCOL BETWEEN '2000 JAN 1' AND 5",
        "SELECT ICOL FROM ALPHA WHERE TCOL LIKE '2000 JAN 1'",
        "SELECT ICOL FROM ALPHA WHERE TCOL LIKE 'garbage time'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'JAN 1'",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND TCOL > '2000 JAN'  AND TCOL < 'tomorrow'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'SCLK 1/1000:5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'NOSUCH_CRAFT SCLK 1/1000:5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '  NO  SUCH   SCLK 1/1000:5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK /5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 1/'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK x/5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 7/5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 2/5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 1:2:3'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 1/abc:5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 99999999999999'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_TDT SCLK 0:10:10'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_TDT SCLK 1/5:10:10:3'",
        # bytes outside ASCII in a time value
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 1/é'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK é/5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 1/5:é'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK 1/5:\xff.\x01'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'EKQ_CRAFT SCLK é:é:é:é:é:é:é:é:é:é:é:é'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > 'é SCLK 5'",
        "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1 é'",
        "SELECT é FROM ALPHA",
        # semantic check
        "SELECT ICOL FROM ALPHA WHERE ICOL LIKE 'A*'",
        "SELECT ICOL FROM ALPHA WHERE DCOL NOT LIKE 'A*'",
        "SELECT ICOL FROM ALPHA WHERE ICOL LIKE CCOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND DCOL LIKE 'x'",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 'A'",
        "SELECT ICOL FROM ALPHA WHERE DCOL = 'A'",
        "SELECT ICOL FROM ALPHA WHERE CCOL = 1",
        "SELECT ICOL FROM ALPHA WHERE CCOL = 1.5",
        "SELECT ICOL FROM ALPHA WHERE CCOL LIKE 5",
        "SELECT ICOL FROM ALPHA WHERE CCOL BETWEEN 'A' AND 5",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 'a#b'",
        "SELECT ICOL FROM ALPHA WHERE ICOL = CCOL",
        "SELECT ICOL FROM ALPHA WHERE CCOL = DCOL",
        "SELECT ICOL FROM ALPHA WHERE TCOL = ICOL",
        "SELECT ICOL FROM ALPHA WHERE DCOL = TCOL",
        "SELECT ICOL FROM ALPHA WHERE CCOL = TCOL",
        "SELECT ICOL FROM ALPHA WHERE TCOL = NAME",
        "SELECT A.ICOL FROM ALPHA A, BETA B WHERE A.ICOL = 1 AND A.NAME = B.XCOL",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND IARR = 1",
        "SELECT ICOL FROM ALPHA WHERE ICOL = DCOL AND DARR > 1.5",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND ICOL = IARR",
        "SELECT ICOL FROM ALPHA WHERE ICOL = IARR",
        "SELECT ICOL FROM ALPHA WHERE CCOL = CARR",
        "SELECT ICOL FROM ALPHA WHERE DCOL < ALPHA.DARR",
        "SELECT B.ICOL FROM ALPHA A, BETA B WHERE B.TCOL = B.TARR",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND CARR LIKE 'x'",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 AND CARR IS NULL",
        "SELECT ICOL FROM ALPHA WHERE ICOL = 1 ORDER BY IARR",
        "SELECT ICOL FROM ALPHA WHERE ICOL = DCOL ORDER BY ICOL, ALPHA.DARR DESC",
    ]
    # the limits
    push!(q, "SELECT T1.ICOL FROM " * join(("ALPHA T$k" for k in 1:11), ", "))
    push!(q, "SELECT " * join(("ICOL" for _ in 1:51), ", ") * " FROM ALPHA")
    push!(q, "SELECT ICOL FROM ALPHA ORDER BY " * join(("ICOL" for _ in 1:11), ", "))
    push!(q, "SELECT ICOL FROM ALPHA WHERE " * join(("ICOL = 1" for _ in 1:125), " AND "))
    push!(q, "SELECT ICOL FROM ALPHA WHERE " * join(("ICOL = 1" for _ in 1:101), " OR "))
    push!(q, "SELECT ICOL FROM ALPHA WHERE " * join(("ICOL = X" for _ in 1:101), " OR ") * " 5 6")
    push!(q, "SELECT ICOL FROM ALPHA WHERE " * "("^250 * "ICOL = X" * ")"^250)
    # distribution of AND over OR: 4^6 conjunctions do not fit
    push!(q, "SELECT ICOL FROM ALPHA WHERE " *
             join(("(ICOL = $k OR DCOL = $k OR ICOL > 9$k OR DCOL < $k)" for k in 1:6), " AND "))
    # negation of 13 conjunctions of 2 relations: 2^13 conjunctions do not fit
    push!(q, "SELECT ICOL FROM ALPHA WHERE NOT (" *
             join(("(ICOL = $k AND DCOL = $k)" for k in 1:13), " OR ") * ")")
    push!(q, "SELECT ICOL FROM ALPHA WHERE NOT (" *
             join(("ICOL BETWEEN $k AND 9$k" for k in 1:12), " OR ") * ")")
    return q
end

"Scalar columns of the tables by type, for the random queries."
const EKQ_SCALAR = Dict(
    "ALPHA" => [("ICOL", :INT), ("DCOL", :DP), ("CCOL", :CHR), ("VCOL", :CHR), ("TCOL", :TIME),
        ("NAME", :CHR)],
    "BETA" => [("ICOL", :INT), ("XCOL", :DP), ("TCOL", :TIME), ("NAME", :CHR), ("BCHR", :CHR)],
    "GAMMA_LONG_NAME" => [("KEY", :INT), ("VAL", :DP), ("LABEL", :CHR), ("EPOCH", :TIME),
        ("ICOL", :INT)])
const EKQ_ARRAY = Dict("ALPHA" => ["IARR", "DARR", "CARR"], "BETA" => ["TARR"],
    "GAMMA_LONG_NAME" => String[])
const EKQ_FROM = [["ALPHA" => ""], ["ALPHA" => "A"], ["BETA" => ""], ["ALPHA" => "A", "BETA" => "B"],
    ["ALPHA" => "", "BETA" => ""], ["ALPHA" => "X", "BETA" => "", "GAMMA_LONG_NAME" => "G"],
    ["ALPHA" => "A1", "ALPHA" => "A2"], ["GAMMA_LONG_NAME" => ""]]
const EKQ_TIMES = ["'2000 JAN 1'", "'1995 MAR 3 12:30:00'", "\"2010-06-15T01:02:03\"",
    "'JD 2451545.5'", "'1 FEB 2003 10:00 PM'", "'2020 DEC 31 23:59:59.999 TDB'", "'tomorrow'"]
const EKQ_STRINGS = ["'A'", "'S11'", "\"x y\"", "'it''s'", "'*S*'", "'%1'", "' '", "'m#n'"]

"A random query with a WHERE clause of nested AND, OR, NOT, BETWEEN and parentheses; most are valid, some have a type error, a wrong name, an array column or a syntax slip."
function ekq_random_query(rng)
    from = rand(rng, EKQ_FROM)
    slip = rand(rng) < 0.2              # allow errors of names and types
    function column(type)
        while true
            (table, alias) = rand(rng, from)
            cands = [c for c in EKQ_SCALAR[table] if type === :ANY || c[2] === type ||
                     (type === :NUM && c[2] in (:INT, :DP))]
            isempty(cands) && continue
            name, t = rand(rng, cands)
            (slip && rand(rng) < 0.15 && !isempty(EKQ_ARRAY[table])) &&
                (name = rand(rng, EKQ_ARRAY[table]))
            unique = count(p -> any(c -> c[1] == name, EKQ_SCALAR[p[1]]), from) == 1
            qual = alias != "" ? alias : table
            ((unique && rand(rng) < 0.6) || (slip && rand(rng) < 0.1)) && return name, t
            return (rand(rng) < 0.5 ? lowercase(qual) : qual) * "." * name, t
        end
    end
    number() = rand(rng, ("1", "-7", "+42", "3.5", "-.25", "1e2", "2.5D-1", "100000", "0"))
    function value(t)
        slip && rand(rng) < 0.2 && (t = rand(rng, (:INT, :CHR, :TIME)))
        t === :CHR && return rand(rng, EKQ_STRINGS)
        t === :TIME && return rand(rng, slip ? EKQ_TIMES : EKQ_TIMES[1:end-1])
        return number()
    end
    function operand(t)
        rand(rng) < 0.25 && return column(t === :INT || t === :DP ? :NUM : t)[1]
        return value(t)
    end
    function atom()
        name, t = column(:ANY)
        r = rand(rng)
        if r < 0.5
            ops = t === :CHR ? ("=", "!=", "<", ">=", "LIKE", "NOT LIKE", "EQ", "LT") :
                  ("=", "!=", "<>", "^=", "<", "<=", ">", ">=", "EQ", "NE", "LT", "LE", "GT", "GE")
            return "$name $(rand(rng, ops)) $(operand(t))"
        elseif r < 0.8
            return "$name $(rand(rng, ("BETWEEN", "NOT BETWEEN"))) $(operand(t)) AND $(operand(t))"
        else
            return "$name $(rand(rng, ("IS NULL", "IS NOT NULL", "= NULL", "NE NULL")))"
        end
    end
    natom = Ref(0)
    function expr(depth)
        r = rand(rng)
        if depth == 0 || natom[] >= 5 || r < 0.35
            natom[] += 1
            return atom()
        end
        wrap(s) = rand(rng) < 0.5 ? "($s)" : s
        r < 0.6 && return wrap(expr(depth - 1)) * " AND " * wrap(expr(depth - 1))
        r < 0.85 && return wrap(expr(depth - 1)) * " OR " * wrap(expr(depth - 1))
        r < 0.95 && return "NOT " * (rand(rng) < 0.7 ? "(" * expr(depth - 1) * ")" : expr(depth - 1))
        return "(" * expr(depth - 1) * ")"
    end
    where = expr(rand(rng, 1:3))
    slip && rand(rng) < 0.2 && (where = replace(where, rand(rng, (" AND " => " ", ")" => "", "(" => "",
        " = " => " ", "NULL" => "NOT")); count=1))
    tabs = join((a == "" ? t : "$t $a" for (t, a) in from), ", ")
    sel = join((column(:ANY)[1] for _ in 1:rand(rng, 1:3)), ", ")
    order = rand(rng) < 0.3 ? " ORDER BY " * join((column(:ANY)[1] * rand(rng, ("", " ASC", " DESC"))
                                                   for _ in 1:rand(rng, 1:2)), ", ") : ""
    return "SELECT $sel FROM $tabs WHERE $where$order"
end

const EKQ_VOCABULARY = vcat(collect(NeoSpice._EKQ_KEYWORDS), ["ICOL", "DCOL", "CCOL", "TCOL",
        "IARR", "NAME", "XCOL", "ALPHA", "BETA", "A", "B", "NOSUCH", "alpha", "icol"],
    ["1", "2.5", "-3", "1e2", ".5", "'x'", "\"y z\"", "'2000 JAN 1'", "''", "'"],
    ["(", ")", ",", ".", "*", "=", "!=", "<>", "<", "<=", ">", ">=", "^=", "(", ")", ",", ".",
        "-", "@", ";"])

"A random string of tokens of the query language: a skeleton of clauses filled with random tokens, or a valid query with one token removed, doubled, swapped or replaced."
function ekq_fuzz_query(rng, valid::Vector{String})
    words(n) = join((rand(rng, EKQ_VOCABULARY) for _ in 1:n), " ")
    names() = join((rand(rng, ("ICOL", "DCOL", "A.ICOL", "ALPHA.CCOL", "B.XCOL", "TCOL", "NAME",
            "IARR", "NOSUCH", words(1))) for _ in 1:rand(rng, 1:3)), rand(rng, (", ", ", ", " ")))
    r = rand(rng)
    if r < 0.15
        return words(rand(rng, 1:12))
    elseif r < 0.55
        from = rand(rng, ("ALPHA", "ALPHA A", "ALPHA A, BETA B", "ALPHA, BETA", words(rand(rng, 1:3))))
        where = rand(rng) < 0.8 ? " WHERE " * words(rand(rng, 0:9)) : ""
        order = rand(rng) < 0.4 ? " ORDER BY " * (rand(rng) < 0.5 ? names() : words(rand(rng, 0:4))) : ""
        parts = ["SELECT " * (rand(rng) < 0.7 ? names() : words(rand(rng, 0:3))), " FROM " * from,
            where, order]
        rand(rng) < 0.2 && shuffle!(rng, parts)
        return strip(join(parts))
    end
    toks = [m.match for m in eachmatch(r"'[^']*'|\"[^\"]*\"|[A-Za-z][A-Za-z0-9_$]*|[0-9.eE+-]+|\S",
        rand(rng, valid))]
    length(toks) > 60 && return words(5)
    k = rand(rng, eachindex(toks))
    r = rand(rng)
    if r < 0.3
        deleteat!(toks, k)
    elseif r < 0.5
        insert!(toks, k, toks[k])
    elseif r < 0.7
        j = rand(rng, eachindex(toks))
        toks[k], toks[j] = toks[j], toks[k]
    else
        toks[k] = rand(rng, EKQ_VOCABULARY)
    end
    return join(toks, " ")
end

# --- the group -------------------------------------------------------------------------

group("ek_parse: EK query parser against CSPICE", Kernels.LSK) do
    dir = mktempdir()
    ekpath = ekq_write_ek(joinpath(dir, "ekq.bes"))
    names = joinpath(dir, "ekq_names.tk")
    write(names, """
    KPL/FK
    \\begindata
    NAIF_BODY_NAME += ( 'EKQ_CRAFT', 'EKQ_TDT' )
    NAIF_BODY_CODE += ( $(CKFix.SC_TDB), $(CKFix.SC_TDT) )
    \\begintext
    """)
    text = [Kernels.LSK, CKFix.sclk_kernel(), names]
    oracle_load(text)
    handle = Ref{Cint}(0)
    ccall((:eklef_c, DSKLIB), Cvoid, (Cstring, Ref{Cint}), ekpath, handle)
    c_check()
    try
        ctx = neo_context(text)
        fx = (tables=ekq_catalogue(), ks=NeoSpice.kernelset(ctx), td=NeoSpice.TimeDefaults())
        report(name, list) = for (query, status, detail) in first(list, 12)
            println("  ", name, " ", status, ": ", first(query, 150), "\n      ", detail)
        end
        function run(queries; guard_not=false)
            res = [(query, ekq_compare(fx, query; guard_not)...) for query in queries]
            return res, [r for r in res if !(r[2] in (:valid, :invalid, :signal, :unsafe))]
        end

        @testset "catalogue" begin
            @test [t.name for t in fx.tables] == first.(EKQ_TABLES)
            for (t, (_, cols)) in zip(fx.tables, EKQ_TABLES)
                @test t.cnames == first.(cols)
                @test [a.dtype for a in t.attrs] == getindex.(cols, 3)
                @test [a.size for a in t.attrs] == getindex.(cols, 4)
            end
        end

        @testset "valid queries" begin
            valid = ekq_valid_corpus()
            res, bad = run(valid)
            report("valid", bad)
            notvalid = [r for r in res if r[2] !== :valid]
            report("not valid", setdiff(notvalid, bad))
            measured(@sprintf("EK parser: %d valid queries by hand, %d equal to CSPICE in the encoded query and in ekpsel",
                length(valid), count(r -> r[2] === :valid, res)))
            @test length(valid) >= 150
            @test isempty(notvalid)
        end

        @testset "invalid queries" begin
            invalid = ekq_invalid_corpus()
            res, bad = run(invalid)
            report("invalid", bad)
            notinvalid = [r for r in res if r[2] !== :invalid]
            report("not invalid", setdiff(notinvalid, bad))
            messages = unique(replace(r[3], r"\d+" => "0") for r in res if r[2] === :invalid)
            measured(@sprintf("EK parser: %d invalid queries by hand, %d with flag and message equal to CSPICE, %d distinct messages",
                length(invalid), count(r -> r[2] === :invalid, res), length(messages)))
            @test length(invalid) >= 100
            @test isempty(notinvalid)
            @test length(messages) >= 60
        end

        @testset "random queries" begin
            rng = MersenneTwister(SEED)
            queries = unique(ekq_random_query(rng) for _ in 1:2500)
            res, bad = run(queries)
            report("random", bad)
            nvalid = count(r -> r[2] === :valid, res)
            ninvalid = count(r -> r[2] === :invalid, res)
            measured(@sprintf("EK parser: %d random queries, %d valid and %d invalid equal to CSPICE, %d different, %d not sent to CSPICE (undefined read)",
                length(queries), nvalid, ninvalid, length(bad), count(r -> r[2] === :unsafe, res)))
            @test isempty(bad)
            @test nvalid >= 1000
            @test ninvalid >= 100
        end

        @testset "random token strings" begin
            rng = MersenneTwister(SEED + 1)
            valid = ekq_valid_corpus()
            queries = unique(String(ekq_fuzz_query(rng, valid)) for _ in 1:4000)
            res, bad = run(queries; guard_not=true)
            report("fuzz", bad)
            nvalid = count(r -> r[2] === :valid, res)
            ninvalid = count(r -> r[2] === :invalid, res)
            messages = unique(replace(r[3], r"\d+" => "0", r"<[^>]*>" => "<>", r": .*$" => "") for r in res if r[2] === :invalid)
            measured(@sprintf("EK parser: %d random token strings, %d valid and %d invalid equal to CSPICE, %d different, %d signal, %d not sent to CSPICE, %d message forms",
                length(queries), nvalid, ninvalid, length(bad), count(r -> r[2] === :signal, res),
                count(r -> r[2] === :unsafe, res), length(messages)))
            get(ENV, "EKQ_MESSAGES", "") == "" || foreach(println, sort(messages))
            @test isempty(bad)
            @test ninvalid >= 2000
        end

        @testset "threads" begin
            # The parser has no shared state: the same queries on all threads give the
            # results of one thread.
            queries = vcat(ekq_valid_corpus(), ekq_invalid_corpus())
            enc(query) = NeoSpice._ek_encode_query(fx.tables, fx.ks, fx.td, query)
            same(a, b) = a[2] == b[2] && a[3] == b[3] && (a[1] === nothing) == (b[1] === nothing) &&
                         (a[1] === nothing || (a[1].conjsizes == b[1].conjsizes &&
                           [c.dval for c in a[1].constraints] == [c.dval for c in b[1].constraints] &&
                           a[1].select == b[1].select))
            serial = map(enc, queries)
            parallel = Vector{Any}(undef, 4 * length(queries))
            Threads.@threads for k in eachindex(parallel)
                parallel[k] = enc(queries[mod1(k, length(queries))])
            end
            @test all(same(parallel[k], serial[mod1(k, length(queries))]) for k in eachindex(parallel))
        end

        @testset "SPICE errors" begin
            # 8 x 8 x 8 conjunctions of 3 relations: the integer cell is too small.
            big = "SELECT ICOL FROM ALPHA WHERE " * join(("(" * join(("ICOL = $(10k + j)"
                                                                      for j in 1:8), " OR ") * ")" for k in 1:3), " AND ")
            @test ekq_compare(fx, big) == (:signal, "celltoosmall")
            # 100 numbers and one time: the double buffer is too small.
            nums = "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1' AND (" *
                   join(("ICOL = $k" for k in 1:100), " OR ") * ")"
            @test ekq_compare(fx, nums) == (:signal, "buffertoosmall")
            # A qualified order-by column before FROM, WHERE or SELECT: CSPICE writes no
            # sense, the SELECT columns are shifted, and name resolution signals a bug.
            for query in ("SELECT ICOL FROM ALPHA ORDER BY ALPHA.ICOL WHERE ICOL > 1",
                "FROM ALPHA ORDER BY ALPHA.ICOL SELECT ICOL")
                @test ekq_compare(fx, query) == (:signal, "bug")
            end
            # More than 1000 conjunctions: CSPICE stops the program (array bound check of
            # f2c), so only NeoSpice is called.
            huge = "SELECT ICOL FROM ALPHA WHERE (" * join(("ICOL = $k" for k in 1:32), " OR ") *
                   ") AND (" * join(("DCOL = $k" for k in 1:32), " OR ") * ")"
            @test ekq_compare(fx, huge)[1] === :unsafe
        end

        @testset "documented differences" begin
            Q = NeoSpice
            enc(query) = Q._ek_encode_query(fx.tables, fx.ks, fx.td, query)
            # More than 2000 characters: CSPICE reads past its buffer.
            long = "SELECT ICOL FROM ALPHA WHERE CCOL = '" * "x"^2000 * "'"
            @test enc(long) == (nothing, true,
                "Query has $(sizeof(long)) characters; the maximum allowed is 2000.")
            @test enc("SELECT ICOL FROM ALPHA" * " "^3000)[2] == false
            # A string of length zero: ekfind_c signals SPICE(EMPTYSTRING).
            @test enc("") == (nothing, true,
                "Every query must contain a FROM clause. The FROM keyword was not found.")
            # A non-scalar column on the left of the first constraint, or in ORDER BY
            # without a constraint before it: CSPICE prints the query up to an undefined
            # position (a segmentation fault was seen); NeoSpice prints the column name.
            # After IS NULL the semantic check of both engines has stopped already.
            for (query, msg) in (
                "SELECT ICOL FROM ALPHA WHERE IARR = 1" =>
                    "Non-scalar column <IARR> having size 3 found in query constraint.",
                "SELECT ICOL FROM ALPHA WHERE ALPHA.DARR > 1.5" =>
                    "Non-scalar column <DARR> having size -1 found in query constraint.",
                "SELECT ICOL FROM ALPHA ORDER BY CARR" =>
                    "Non-scalar column <CARR> having size 2 found in order-by column.",
                "SELECT ICOL FROM ALPHA WHERE ICOL IS NULL ORDER BY CARR" => "")
                r = enc(query)
                @test msg == "" ? !r[2] : r == (nothing, true, msg)
                @test ekq_compare(fx, query)[1] === (msg == "" ? :valid : :unsafe)
            end
            # A time that tpartv accepts and str2et rejects, and a missing kernel: CSPICE
            # sets the flag and the message and leaves a SPICE error signalled as well.
            for query in ("SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1 12:00:00 TDB PST'",
                "SELECT ICOL FROM ALPHA WHERE TCOL > '1997 JAN 32 12:29:29'",
                "SELECT ICOL FROM ALPHA WHERE TCOL > 'MOON SCLK 1/1000:5'")
                n = enc(query)
                c = ekq_c_encd(query)
                @test n[2] && c.error && c.kind !== nothing
                @test n[3] == c.errmsg
                @test startswith(n[3], "Unexpected SPICELIB error encountered")
            end
            none = Q._ek_encode_query(fx.tables, NeoSpice.empty_kernelset(), fx.td,
                "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1'")
            @test none == (nothing, true, "Unexpected SPICELIB error encountered while " *
                                          "attempting to parse the string <2000 JAN 1>")
        end

        @testset "decoded query" begin
            Q = NeoSpice
            query = "select a.icol, vcol, b.tarr from ALPHA a, BETA b where a.icol >= 2 and " *
                    "(ccol like 'S*' or b.xcol not between 1 and a.dcol) and " *
                    "a.tcol > '2000 JAN 1 12:00:00 TDB' order by b.icol desc, ccol"
            q, error, errmsg = Q._ek_encode_query(fx.tables, fx.ks, fx.td, query)
            @test !error && errmsg == ""
            @test q.tables == ["ALPHA", "BETA"] && q.aliases == ["A", "B"] && q.tabcat == [1, 2]
            @test q.conjsizes == [3, 4, 4]
            @test [(c.lhs.name, c.opcode) for c in q.constraints] == [
                ("ICOL", Q._EKQ_GE), ("CCOL", Q._EKQ_LIKE), ("TCOL", Q._EKQ_GT),
                ("ICOL", Q._EKQ_GE), ("XCOL", Q._EKQ_LT), ("XCOL", Q._EKQ_LT), ("TCOL", Q._EKQ_GT),
                ("ICOL", Q._EKQ_GE), ("XCOL", Q._EKQ_GT), ("XCOL", Q._EKQ_GT), ("TCOL", Q._EKQ_GT)]
            c = q.constraints
            @test c[1].cnstyp == Q._EKQ_EQVAL && c[1].dtype == 3 && c[1].ival == 2 && c[1].dval == 2.0
            at = first(findfirst("a.icol >=", query))
            @test c[1].lhs == Q._EKQCol("A", 1, "ICOL", 1, at, at + 5)
            @test c[2].dtype == 1 && c[2].cval == "S*" && c[2].lhs.table == ""
            @test c[3].dtype == 4 && c[3].dval == 0.0          # J2000 TDB
            @test c[5].dtype == 3 && c[5].ival == 1
            @test c[6].cnstyp == Q._EKQ_EQCOL && c[6].rhs.name == "DCOL" && c[6].rhs.tabidx == 1 &&
                  c[6].rhs.colidx == 2 && c[6].dtype == 0
            @test [(o[1].name, o[1].tabidx, o[2]) for o in q.orderby] ==
                  [("ICOL", 2, Q._EKQ_DESCEND), ("CCOL", 1, Q._EKQ_ASCEND)]
            @test [(s.table, s.tabidx, s.name, s.colidx, s.lxbeg, s.lxend) for s in q.select] ==
                  [("A", 1, "ICOL", 1, 8, 13), ("", 1, "VCOL", 4, 16, 19), ("B", 2, "TARR", 5, 22, 27)]
            s = Q._ek_psel(fx.tables, fx.ks, fx.td, "SELECT B.TCOL, DCOL, CARR, B.ICOL FROM ALPHA, BETA B")
            @test s.n == 4 && s.xtypes == [:TIME, :DP, :CHR, :INT] && s.xclass == fill(:COL, 4)
            @test s.tabs == ["BETA", "ALPHA", "ALPHA", "BETA"] && s.xbegs == [8, 16, 22, 28]
            # No WHERE clause: no conjunction, all rows are selected.
            q0 = Q._ek_encode_query(fx.tables, fx.ks, fx.td, "SELECT ICOL FROM ALPHA")[1]
            @test isempty(q0.conjsizes) && isempty(q0.constraints) && isempty(q0.orderby)
            # The integer of a value that does not fit in 32 bits, as i_dnnt gives it.
            big = Q._ek_encode_query(fx.tables, fx.ks, fx.td, "SELECT ICOL FROM ALPHA WHERE ICOL = 12345678901")[1]
            @test big.constraints[1].ival == typemin(Int32) && big.constraints[1].dval == 12345678901.0
            # The time defaults of the caller are used.
            tdb = Q.TimeDefaults("GREGORIAN", "TDB", "")
            qt = Q._ek_encode_query(fx.tables, fx.ks, tdb, "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1 12:00:00'")[1]
            @test qt.constraints[1].dval == 0.0
        end

        @testset "column names (zzekpcol)" begin
            Q = NeoSpice
            query = "SELECT A.ICOL FROM ALPHA A, BETA, GAMMA_LONG_NAME G WHERE A.ICOL > 0"
            q = Q._ek_encode_query(fx.tables, fx.ks, fx.td, query)[1]
            e = ekq_c_encd(query)
            for qcol in ("DCOL", "dcol", " xcol ", "A.ICOL", "a . icol", "ALPHA.NAME", "BETA.ICOL",
                "G.KEY", "g.epoch", "GAMMA_LONG_NAME.VAL", "ICOL", "NAME", "NOSUCH", "A.NOSUCH",
                "B.ICOL", "Z.ICOL", "A.", ".ICOL", "A.ICOL.X", "A ICOL", "5", "A.5", "'x'", "A,ICOL",
                "A.'x'", "1.ICOL", "@", "A.ICOL B", "KEY", "G . VAL", "X"^70 * ".ICOL")
                table, alias, column = (fill(UInt8(' '), 64) for _ in 1:3)
                tabidx, colidx, error = Ref{Cint}(0), Ref{Cint}(0), Ref{Cint}(0)
                errmsg = fill(UInt8(' '), EKQ_MSGLEN)
                ccall((:zzekpcol_, DSKLIB), Cint, (Ptr{UInt8}, Ptr{Cint}, Ptr{UInt8}, Ptr{UInt8},
                        Ptr{UInt8}, Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Ref{Cint}, Ptr{UInt8}, Cint,
                        Cint, Cint, Cint, Cint, Cint), qcol, e.raw, e.c, table, alias, tabidx,
                    column, colidx, error, errmsg, sizeof(qcol), 2000, 64, 64, 64, EKQ_MSGLEN)
                c_check()
                r = Q._ek_qpcol(fx.tables, q, qcol)
                @test r[6] == (error[] != 0)
                @test r[7] == ekq_fstr(errmsg)
                if error[] == 0 && !r[6]
                    @test r[1:5] == (tabidx[], ekq_fstr(table), ekq_fstr(alias), colidx[],
                        ekq_fstr(column))
                end
            end
        end
    finally
        ccall((:ekuef_c, DSKLIB), Cvoid, (Cint,), handle[])
        c_pending_kind()
    end
end
