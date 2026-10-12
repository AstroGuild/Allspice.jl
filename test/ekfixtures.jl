# Fixtures of the EK groups: direct calls into CSPICE and EK files that the tests generate
# (no EK file is needed on the machine).
#
# A logical model (`EKSegSpec`: one segment with its column declarations and its rows)
# is written by CSPICE (`c_ek_write`) or by NeoSpice (`neo_ek_write`) with the same
# sequence of calls. CSPICE is called from the main task only.
#
# The C interface of CSPICE numbers segments, rows and elements from 0; the wrappers
# here take and return numbers from 1, as NeoSpice and SPICE.jl do.

isdefined(@__MODULE__, :DSKLIB) || include("dskfixtures.jl")

import Libdl
const EKLIB = SPICE.libcspice
const EKDL = Libdl.dlopen(EKLIB)

"One column of a generated segment. `strlen`: 0 not a string, -1 variable. `size`: -1 variable."
struct EKColSpec
    name::String
    dtype::Symbol          # :CHR, :DP, :INT, :TIME
    strlen::Int
    size::Int
    indexed::Bool
    nullok::Bool
    fixed::Bool            # FIXED_COUNT = TRUE (column class 7 to 9)
end
EKColSpec(name, dtype; strlen=0, size=1, indexed=false, nullok=false, fixed=false) =
    EKColSpec(name, dtype, strlen, size, indexed, nullok, fixed)

"The declaration string of `c` for `ekbseg` and `ekifld`."
function ekdecl(c::EKColSpec)
    t = c.dtype === :CHR ? "CHARACTER*(" * (c.strlen == -1 ? "*" : string(c.strlen)) * ")" :
        c.dtype === :DP ? "DOUBLE PRECISION" : c.dtype === :INT ? "INTEGER" : "TIME"
    d = "DATATYPE = $t"
    c.size == 1 || (d *= ", SIZE = " * (c.size == -1 ? "VARIABLE" : string(c.size)))
    c.indexed && (d *= ", INDEXED = TRUE")
    c.nullok && (d *= ", NULLS_OK = TRUE")
    c.fixed && (d *= ", FIXED_COUNT = TRUE")
    return d
end

"""
One generated segment: table name, columns, and the rows (`rows[r][c]` is `nothing` for a
null entry, else the vector of the values of the entry). `fast` selects the fast writers.
"""
struct EKSegSpec
    table::String
    cols::Vector{EKColSpec}
    rows::Vector{Vector{Any}}
    fast::Bool
end

"A C array of `length(strs)` strings of `len` bytes each (NUL-terminated)."
function c_strarray(strs, len::Int)
    buf = zeros(UInt8, max(1, length(strs)) * len)
    for (i, s) in enumerate(strs)
        cu = codeunits(String(s))
        n = min(length(cu), len - 1)
        copyto!(buf, (i - 1) * len + 1, cu, 1, n)
    end
    return buf
end

function c_ekopn(path, ifname="EK", ncomch=0)
    h = Ref{Cint}(0)
    ccall((:ekopn_c, EKLIB), Cvoid, (Cstring, Cstring, Cint, Ref{Cint}), path, ifname, ncomch, h)
    c_check()
    return h[]
end
function c_ekopr(path)
    h = Ref{Cint}(0)
    ccall((:ekopr_c, EKLIB), Cvoid, (Cstring, Ref{Cint}), path, h)
    c_check()
    return h[]
end
function c_ekopw(path)
    h = Ref{Cint}(0)
    ccall((:ekopw_c, EKLIB), Cvoid, (Cstring, Ref{Cint}), path, h)
    c_check()
    return h[]
end
c_ekcls(h) = (ccall((:ekcls_c, EKLIB), Cvoid, (Cint,), h); c_check())
c_eknseg(h) = (n = ccall((:eknseg_c, EKLIB), Cint, (Cint,), h); c_check(); Int(n))

function c_ekbseg(h, table, cnames, decls)
    segno = Ref{Cint}(0)
    ccall((:ekbseg_c, EKLIB), Cvoid, (Cint, Cstring, Cint, Cint, Ptr{UInt8}, Cint, Ptr{UInt8},
            Ref{Cint}), h, table, length(cnames), 33, c_strarray(cnames, 33), 201,
        c_strarray(decls, 201), segno)
    c_check()
    return segno[] + 1
end
function c_ekappr(h, segno)
    recno = Ref{Cint}(0)
    ccall((:ekappr_c, EKLIB), Cvoid, (Cint, Cint, Ref{Cint}), h, segno - 1, recno)
    c_check()
    return recno[] + 1
end
c_ekinsr(h, segno, recno) = (ccall((:ekinsr_c, EKLIB), Cvoid, (Cint, Cint, Cint), h, segno - 1,
    recno - 1); c_check())
c_ekdelr(h, segno, recno) = (ccall((:ekdelr_c, EKLIB), Cvoid, (Cint, Cint, Cint), h, segno - 1,
    recno - 1); c_check())

"Largest string length of `vals` plus one, at least 2 (the `vallen` of the C calls)."
c_vallen(vals) = max(2, maximum(sizeof, vals; init=0) + 1)

function c_ekacec(h, segno, recno, column, vals, isnull; fn=:ekacec_c)
    len = c_vallen(vals)
    ccall(Libdl.dlsym(EKDL, fn), Cvoid, (Cint, Cint, Cint, Cstring, Cint, Cint, Ptr{UInt8}, Cint), h,
        segno - 1, recno - 1, column, length(vals), len, c_strarray(vals, len), isnull)
    c_check()
end
function c_ekaced(h, segno, recno, column, vals, isnull; fn=:ekaced_c)
    ccall(Libdl.dlsym(EKDL, fn), Cvoid, (Cint, Cint, Cint, Cstring, Cint, Ptr{Float64}, Cint), h,
        segno - 1, recno - 1, column, length(vals), Float64[v for v in vals], isnull)
    c_check()
end
function c_ekacei(h, segno, recno, column, vals, isnull; fn=:ekacei_c)
    ccall(Libdl.dlsym(EKDL, fn), Cvoid, (Cint, Cint, Cint, Cstring, Cint, Ptr{Cint}, Cint), h,
        segno - 1, recno - 1, column, length(vals), Cint[v for v in vals], isnull)
    c_check()
end
c_ekucec(h, segno, recno, column, vals, isnull) = c_ekacec(h, segno, recno, column, vals, isnull; fn=:ekucec_c)
c_ekuced(h, segno, recno, column, vals, isnull) = c_ekaced(h, segno, recno, column, vals, isnull; fn=:ekuced_c)
c_ekucei(h, segno, recno, column, vals, isnull) = c_ekacei(h, segno, recno, column, vals, isnull; fn=:ekucei_c)

function c_ekifld(h, table, nrows, cnames, decls)
    segno = Ref{Cint}(0)
    rcptrs = zeros(Cint, nrows)
    ccall((:ekifld_c, EKLIB), Cvoid, (Cint, Cstring, Cint, Cint, Cint, Ptr{UInt8}, Cint,
            Ptr{UInt8}, Ref{Cint}, Ptr{Cint}), h, table, length(cnames), nrows, 33,
        c_strarray(cnames, 33), 201, c_strarray(decls, 201), segno, rcptrs)
    c_check()
    return segno[] + 1, rcptrs
end
function c_ekaclc(h, segno, column, vals, entszs, nlflgs, rcptrs)
    # CSPICE signals INVALIDSIZE if the strings are shorter than the string length of the column.
    len = max(c_vallen(vals), 101)
    wk = zeros(Cint, length(rcptrs))
    ccall((:ekaclc_c, EKLIB), Cvoid, (Cint, Cint, Cstring, Cint, Ptr{UInt8}, Ptr{Cint},
            Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), h, segno - 1, column, len, c_strarray(vals, len),
        Cint[e for e in entszs], Cint[f for f in nlflgs], rcptrs, wk)
    c_check()
end
function c_ekacld(h, segno, column, vals, entszs, nlflgs, rcptrs)
    wk = zeros(Cint, length(rcptrs))
    ccall((:ekacld_c, EKLIB), Cvoid, (Cint, Cint, Cstring, Ptr{Float64}, Ptr{Cint}, Ptr{Cint},
            Ptr{Cint}, Ptr{Cint}), h, segno - 1, column, Float64[v for v in vals],
        Cint[e for e in entszs], Cint[f for f in nlflgs], rcptrs, wk)
    c_check()
end
function c_ekacli(h, segno, column, vals, entszs, nlflgs, rcptrs)
    wk = zeros(Cint, length(rcptrs))
    ccall((:ekacli_c, EKLIB), Cvoid, (Cint, Cint, Cstring, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
            Ptr{Cint}, Ptr{Cint}), h, segno - 1, column, Cint[v for v in vals],
        Cint[e for e in entszs], Cint[f for f in nlflgs], rcptrs, wk)
    c_check()
end
c_ekffld(h, segno, rcptrs) = (ccall((:ekffld_c, EKLIB), Cvoid, (Cint, Cint, Ptr{Cint}), h,
    segno - 1, rcptrs); c_check())

"`ekrcec_c`: (values without trailing blanks, isnull). `maxn` and `maxlen` size the buffer."
function c_ekrcec(h, segno, recno, column; maxn=1100, maxlen=2100)
    n = Ref{Cint}(0)
    isnull = Ref{Cint}(0)
    buf = zeros(UInt8, maxn * maxlen)
    ccall((:ekrcec_c, EKLIB), Cvoid, (Cint, Cint, Cint, Cstring, Cint, Ref{Cint}, Ptr{UInt8},
            Ref{Cint}), h, segno - 1, recno - 1, column, maxlen, n, buf, isnull)
    c_check()
    isnull[] != 0 && return String[], true
    out = Vector{String}(undef, n[])
    for i in 1:n[]
        # The C interface returns each string padded with blanks to `maxlen - 1` characters.
        first = (i - 1) * maxlen + 1
        last = something(findnext(==(0x00), buf, first), first + maxlen) - 1
        out[i] = String(rstrip(String(buf[first:last]), ' '))
    end
    return out, false
end
function c_ekrced(h, segno, recno, column; maxn=1100)
    n = Ref{Cint}(0)
    isnull = Ref{Cint}(0)
    buf = zeros(Float64, maxn)
    ccall((:ekrced_c, EKLIB), Cvoid, (Cint, Cint, Cint, Cstring, Ref{Cint}, Ptr{Float64},
            Ref{Cint}), h, segno - 1, recno - 1, column, n, buf, isnull)
    c_check()
    isnull[] != 0 && return Float64[], true
    return buf[1:n[]], false
end
function c_ekrcei(h, segno, recno, column; maxn=1100)
    n = Ref{Cint}(0)
    isnull = Ref{Cint}(0)
    buf = zeros(Cint, maxn)
    ccall((:ekrcei_c, EKLIB), Cvoid, (Cint, Cint, Cint, Cstring, Ref{Cint}, Ptr{Cint},
            Ref{Cint}), h, segno - 1, recno - 1, column, n, buf, isnull)
    c_check()
    isnull[] != 0 && return Int[], true
    return Int[buf[i] for i in 1:n[]], false
end

"`ekssum_` (f2c): (tabnam, nrows, ncols, cnames, dtypes, sizes, strlns, indexd, nullok)."
function c_ekssum(h, segno)
    tabnam = fill(UInt8(' '), 64)
    cnames = fill(UInt8(' '), 32 * 100)
    dtypes = fill(UInt8(' '), 4 * 100)
    nrows = Ref{Cint}(0)
    ncols = Ref{Cint}(0)
    sizes = zeros(Cint, 100)
    strlns = zeros(Cint, 100)
    indexd = zeros(Cint, 100)
    nullok = zeros(Cint, 100)
    ccall((:ekssum_, EKLIB), Cint, (Ref{Cint}, Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Ref{Cint},
            Ptr{UInt8}, Ptr{UInt8}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Cint, Cint, Cint),
        h, segno, tabnam, nrows, ncols, cnames, dtypes, sizes, strlns, indexd, nullok, 64, 32, 4)
    c_check()
    nc = Int(ncols[])
    fstr(buf, i, len) = String(rstrip(String(buf[(i-1)*len+1:i*len])))
    return (tabnam=String(rstrip(String(tabnam))), nrows=Int(nrows[]), ncols=nc,
        cnames=[fstr(cnames, i, 32) for i in 1:nc], dtypes=[Symbol(fstr(dtypes, i, 4)) for i in 1:nc],
        sizes=Int.(sizes[1:nc]), strlns=Int.(strlns[1:nc]), indexd=indexd[1:nc] .!= 0,
        nullok=nullok[1:nc] .!= 0)
end

# ---------------------------------------------------------------------------------------
# Writing a model with either engine
# ---------------------------------------------------------------------------------------

"Flat values, entry sizes and null flags of column `c` of `seg`, as the fast writers take them."
function ek_flat_column(seg::EKSegSpec, c::Int)
    col = seg.cols[c]
    vals = Any[]
    entszs = Int[]
    nlflgs = Bool[]
    filler = col.dtype === :CHR ? " " : col.dtype === :INT ? 0 : 0.0
    for row in seg.rows
        e = row[c]
        if e === nothing
            # A null entry still takes the room of one entry in the value array.
            n = col.size == -1 ? 1 : col.size
            append!(vals, fill(filler, n))
            push!(entszs, n)
            push!(nlflgs, true)
        else
            append!(vals, e)
            push!(entszs, length(e))
            push!(nlflgs, false)
        end
    end
    return vals, entszs, nlflgs
end

"""
Write `segs` to the new EK file `path` with the writers `w` (a NamedTuple of functions with
the signatures of the `c_ek*` wrappers above). The call sequence is the same for both engines.
"""
function ek_write_model(w, path, segs::Vector{EKSegSpec})
    h = w.ekopn(path, "EK_TEST", 0)
    for seg in segs
        cnames = [c.name for c in seg.cols]
        decls = [ekdecl(c) for c in seg.cols]
        if seg.fast
            segno, rcptrs = w.ekifld(h, seg.table, length(seg.rows), cnames, decls)
            for (c, col) in enumerate(seg.cols)
                vals, entszs, nlflgs = ek_flat_column(seg, c)
                f = col.dtype === :CHR ? w.ekaclc : col.dtype === :INT ? w.ekacli : w.ekacld
                f(h, segno, col.name, vals, entszs, nlflgs, rcptrs)
            end
            w.ekffld(h, segno, rcptrs)
        else
            segno = w.ekbseg(h, seg.table, cnames, decls)
            for row in seg.rows
                recno = w.ekappr(h, segno)
                for (c, col) in enumerate(seg.cols)
                    e = row[c]
                    f = col.dtype === :CHR ? w.ekacec : col.dtype === :INT ? w.ekacei : w.ekaced
                    if e === nothing
                        n = col.size == -1 ? 1 : col.size
                        f(h, segno, recno, col.name, fill(col.dtype === :CHR ? " " :
                            col.dtype === :INT ? 0 : 0.0, n), true)
                    else
                        f(h, segno, recno, col.name, e, false)
                    end
                end
            end
        end
    end
    w.ekcls(h)
    return path
end

const C_EK_WRITERS = (ekopn=c_ekopn, ekcls=c_ekcls, ekbseg=c_ekbseg, ekappr=c_ekappr,
    ekacec=c_ekacec, ekaced=c_ekaced, ekacei=c_ekacei, ekifld=c_ekifld, ekaclc=c_ekaclc,
    ekacld=c_ekacld, ekacli=c_ekacli, ekffld=c_ekffld)

c_ek_write(path, segs) = ek_write_model(C_EK_WRITERS, path, segs)

# ---------------------------------------------------------------------------------------
# Random models
# ---------------------------------------------------------------------------------------

const EK_WORDS = ["alpha", "Alpha", "ALPHA", "beta", "Beta gamma", "gamma", "delta", "DELTA",
    "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron",
    "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega", " lead", "a%b", "a*b",
    "x_y", "it's", "MiXeD cAsE", "0123", "z", "Zeta", "a", "A", "ab", "aB", "abc", "", "  two"]

"A random string for a column of string length `strlen` (-1: variable, sometimes long)."
function ek_random_string(rng, strlen::Int)
    if strlen == -1
        r = rand(rng)
        r < 0.02 && return join(rand(rng, 'a':'z', rand(rng, 1000:2500)))   # crosses pages
        r < 0.75 && return rand(rng, EK_WORDS)
        s = join(rand(rng, EK_WORDS, rand(rng, 1:6)), ' ')
        return isempty(strip(s)) ? "w" : s
    end
    s = rand(rng) < 0.7 ? rand(rng, EK_WORDS) : join(rand(rng, ['a', 'b', 'c', 'A', 'B', ' ', 'x', '%'],
        rand(rng, 1:strlen)))
    return s[1:min(end, strlen)]
end

"A random entry (vector of values, or `nothing` for null) of column `col`; `row` is the row number."
function ek_random_entry(rng, col::EKColSpec, row::Int; keyrange=1:200)
    col.nullok && rand(rng) < 0.12 && return nothing
    n = col.size == -1 ? (rand(rng) < 0.03 ? rand(rng, 250:700) : rand(rng, 1:6)) : col.size
    if col.dtype === :CHR
        return Any[ek_random_string(rng, col.strlen) for _ in 1:n]
    elseif col.dtype === :INT
        return Any[rand(rng) < 0.5 ? rand(rng, keyrange) : rand(rng, -1000:1000) for _ in 1:n]
    elseif col.dtype === :DP
        return Any[rand(rng) < 0.3 ? Float64(rand(rng, -20:20)) : round(1e3 * randn(rng); digits=3) for _ in 1:n]
    else
        # Times near the years 1990 to 2010, as seconds past J2000; many on whole days.
        return Any[rand(rng) < 0.5 ? 86400.0 * rand(rng, -3650:3650) : 3.15e8 * (2 * rand(rng) - 1) for _ in 1:n]
    end
end

"A segment of `nrows` random rows for the table `table` with the columns `cols`."
function ek_random_segment(rng, table, cols, nrows; fast=false, keyrange=1:200)
    rows = Vector{Any}[Any[ek_random_entry(rng, c, r; keyrange=keyrange) for c in cols] for r in 1:nrows]
    return EKSegSpec(table, cols, rows, fast)
end

"Columns of the table EVENTS: every data type, scalar and array entries, classes 1 to 6."
const EK_EVENT_COLS = [
    EKColSpec("ID", :INT; indexed=true),
    EKColSpec("NAME", :CHR; strlen=16, indexed=true, nullok=true),
    EKColSpec("DESCR", :CHR; strlen=-1, nullok=true),
    EKColSpec("VAL", :DP; indexed=true, nullok=true),
    EKColSpec("EPOCH", :TIME; indexed=true),
    EKColSpec("GRP", :INT; nullok=true),
    EKColSpec("LABEL", :CHR; strlen=-1, indexed=true),
    EKColSpec("T2", :TIME; nullok=true),
    EKColSpec("X", :DP),
    EKColSpec("VEC", :DP; size=3),
    EKColSpec("IVEC", :INT; size=-1, nullok=true),
    EKColSpec("TAGS", :CHR; strlen=8, size=-1),
    EKColSpec("WORDS", :CHR; strlen=12, size=4, nullok=true),
    EKColSpec("TVEC", :TIME; size=-1),
    EKColSpec("DVEC", :DP; size=-1, nullok=true),
    EKColSpec("IPAIR", :INT; size=2)]

"Columns of the table GROUPS: fixed-count columns (classes 7 to 9) of a type 2 segment."
const EK_GROUP_COLS = [
    EKColSpec("GID", :INT; indexed=true, fixed=true),
    EKColSpec("GNAME", :CHR; strlen=10, indexed=true, nullok=true, fixed=true),
    EKColSpec("WEIGHT", :DP; indexed=true, nullok=true, fixed=true),
    EKColSpec("T0", :TIME; fixed=true),
    EKColSpec("NOTE", :CHR; strlen=6, fixed=true),
    EKColSpec("K", :INT; nullok=true, fixed=true)]

"Columns of the small table LINKS, for joins of three tables."
const EK_LINK_COLS = [
    EKColSpec("EID", :INT; indexed=true),
    EKColSpec("GID", :INT),
    EKColSpec("TAG", :CHR; strlen=-1, indexed=true, nullok=true),
    EKColSpec("SCORE", :DP; nullok=true),
    EKColSpec("WHEN", :TIME; nullok=true)]

"""
The standard model of the EK groups, as (file name => segments) pairs: EVENTS in three
segments over two files (one made with the fast writers), GROUPS in two type 2 segments,
LINKS in one segment. `scale` multiplies the row counts.
"""
function ek_standard_model(seed::Int; scale::Float64=1.0)
    rng = MersenneTwister(seed)
    n(k) = max(3, round(Int, k * scale))
    ev1 = ek_random_segment(rng, "EVENTS", EK_EVENT_COLS, n(1500))
    ev2 = ek_random_segment(rng, "EVENTS", EK_EVENT_COLS, n(600))
    ev3 = ek_random_segment(rng, "EVENTS", EK_EVENT_COLS, n(700); fast=true)
    gr1 = ek_random_segment(rng, "GROUPS", EK_GROUP_COLS, n(1300); fast=true)
    gr2 = ek_random_segment(rng, "GROUPS", EK_GROUP_COLS, n(250); fast=true)
    lk = ek_random_segment(rng, "LINKS", EK_LINK_COLS, n(120))
    return ["events_a.bes" => [ev1, gr2], "events_b.bes" => [ev2, ev3, lk], "groups.bes" => [gr1]]
end

"Compare the entry `e` of the model (`nothing` or values) with what an engine read: `(vals, isnull)`."
function ek_entry_equal(col::EKColSpec, e, got)
    vals, isnull = got
    e === nothing && return isnull
    isnull && return false
    length(vals) == length(e) || return false
    if col.dtype === :CHR
        return all(String(rstrip(String(a), ' ')) == b for (a, b) in zip(e, vals))
    end
    return all(a == b for (a, b) in zip(e, vals))
end

# ---------------------------------------------------------------------------------------
# Loaded EK data and queries in CSPICE (direct calls; rows, columns, elements from 1)
# ---------------------------------------------------------------------------------------

"The C string in `buf` (up to the first NUL)."
c_string(buf::AbstractVector{UInt8}) = String(buf[1:something(findfirst(iszero, buf), length(buf) + 1)-1])

c_ekntab() = (n = Ref{Cint}(0); ccall((:ekntab_c, EKLIB), Cvoid, (Ref{Cint},), n); c_check(); Int(n[]))
function c_ektnam(n)
    buf = zeros(UInt8, 65)
    ccall((:ektnam_c, EKLIB), Cvoid, (Cint, Cint, Ptr{UInt8}), n - 1, 65, buf)
    c_check()
    return c_string(buf)
end
c_ekccnt(table) = (n = Ref{Cint}(0); ccall((:ekccnt_c, EKLIB), Cvoid, (Cstring, Ref{Cint}), table, n);
    c_check(); Int(n[]))

"`ekcii_c`: (column name, (cclass, dtype symbol, strlen, size, indexd, nullok))."
function c_ekcii(table, cindex)
    buf = zeros(UInt8, 33)
    dsc = zeros(Cint, 6)
    ccall((:ekcii_c, EKLIB), Cvoid, (Cstring, Cint, Cint, Ptr{UInt8}, Ptr{Cint}), table, cindex - 1, 33, buf, dsc)
    c_check()
    return c_string(buf), (Int(dsc[1]), (:CHR, :DP, :INT, :TIME)[dsc[2]+1], Int(dsc[3]), Int(dsc[4]),
        dsc[5] != 0, dsc[6] != 0)
end

"`ekfind_c`: (nmrows, error, errmsg). The result stays in CSPICE for the fetch calls."
function c_ekfind(query)
    n = Ref{Cint}(0)
    err = Ref{Cint}(0)
    buf = zeros(UInt8, 1841)
    ccall((:ekfind_c, EKLIB), Cvoid, (Cstring, Cint, Ref{Cint}, Ref{Cint}, Ptr{UInt8}), query, 1841, n, err, buf)
    c_check()
    return Int(n[]), err[] != 0, c_string(buf)
end

"`ekpsel_c` with positions from 1: the same named tuple as `ekpsel` of NeoSpice."
function c_ekpsel(query)
    n = Ref{Cint}(0)
    err = Ref{Cint}(0)
    xb = zeros(Cint, 50); xe = zeros(Cint, 50); xt = zeros(Cint, 50); xc = zeros(Cint, 50)
    tabs = zeros(UInt8, 50 * 65); cols = zeros(UInt8, 50 * 33)
    buf = zeros(UInt8, 1841)
    ccall((:ekpsel_c, EKLIB), Cvoid, (Cstring, Cint, Cint, Cint, Ref{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
            Ptr{Cint}, Ptr{UInt8}, Ptr{UInt8}, Ref{Cint}, Ptr{UInt8}), query, 1841, 65, 33, n, xb, xe, xt, xc,
        tabs, cols, err, buf)
    c_check()
    k = err[] != 0 ? 0 : Int(n[])
    return (n=k, xbegs=Int.(xb[1:k]) .+ 1, xends=Int.(xe[1:k]) .+ 1,
        xtypes=Symbol[(:CHR, :DP, :INT, :TIME)[xt[i]+1] for i in 1:k],
        xclass=Symbol[(:COL, :FUNC, :EXPR)[xc[i]+1] for i in 1:k],
        tabs=String[c_string(view(tabs, (i-1)*65+1:i*65)) for i in 1:k],
        cols=String[c_string(view(cols, (i-1)*33+1:i*33)) for i in 1:k], error=err[] != 0, errmsg=c_string(buf))
end

c_eknelt(selidx, row) = (n = ccall((:eknelt_c, EKLIB), Cint, (Cint, Cint), selidx - 1, row - 1); c_check(); Int(n))

const C_EKGC_BUF = zeros(UInt8, 2600)
function c_ekgc(selidx, row, elment=1)
    null = Ref{Cint}(0); found = Ref{Cint}(0)
    ccall((:ekgc_c, EKLIB), Cvoid, (Cint, Cint, Cint, Cint, Ptr{UInt8}, Ref{Cint}, Ref{Cint}), selidx - 1,
        row - 1, elment - 1, length(C_EKGC_BUF), C_EKGC_BUF, null, found)
    c_check()
    isnull = found[] != 0 && null[] != 0
    return (found[] != 0 && !isnull) ? String(rstrip(c_string(C_EKGC_BUF), ' ')) : "", isnull, found[] != 0
end
function c_ekgd(selidx, row, elment=1)
    v = Ref{Float64}(0.0); null = Ref{Cint}(0); found = Ref{Cint}(0)
    ccall((:ekgd_c, EKLIB), Cvoid, (Cint, Cint, Cint, Ref{Float64}, Ref{Cint}, Ref{Cint}), selidx - 1, row - 1,
        elment - 1, v, null, found)
    c_check()
    isnull = found[] != 0 && null[] != 0
    return (found[] != 0 && !isnull) ? v[] : 0.0, isnull, found[] != 0
end
function c_ekgi(selidx, row, elment=1)
    v = Ref{Cint}(0); null = Ref{Cint}(0); found = Ref{Cint}(0)
    ccall((:ekgi_c, EKLIB), Cvoid, (Cint, Cint, Cint, Ref{Cint}, Ref{Cint}, Ref{Cint}), selidx - 1, row - 1,
        elment - 1, v, null, found)
    c_check()
    isnull = found[] != 0 && null[] != 0
    return (found[] != 0 && !isnull) ? Int(v[]) : 0, isnull, found[] != 0
end

"""
All values of the query result that CSPICE holds after `c_ekfind`: `out[row][sel]` is
`nothing` for a null entry, else the vector of the elements. `types` are the data types
of the SELECT columns.
"""
function c_ek_result(nmrows::Int, types::Vector{Symbol})
    out = Vector{Vector{Any}}(undef, nmrows)
    for row in 1:nmrows
        vals = Vector{Any}(undef, length(types))
        for (sel, t) in enumerate(types)
            get = t === :CHR ? c_ekgc : t === :INT ? c_ekgi : c_ekgd
            n = c_eknelt(sel, row)
            v, isnull, found = get(sel, row, 1)
            vals[sel] = isnull ? nothing : Any[k == 1 ? v : get(sel, row, k)[1] for k in 1:n]
        end
        out[row] = vals
    end
    return out
end

"The same table of values for an `EKResult` of NeoSpice, read with `eknelt`, `ekgc`, `ekgd`, `ekgi`."
function neo_ek_result(res)
    types = res.seltype
    out = Vector{Vector{Any}}(undef, res.nmrows)
    for row in 1:res.nmrows
        vals = Vector{Any}(undef, length(types))
        for (sel, t) in enumerate(types)
            get = t === :CHR ? ekgc : t === :INT ? ekgi : ekgd
            n = eknelt(res, sel, row)
            v, isnull, found = get(res, sel, row, 1)
            vals[sel] = isnull ? nothing : Any[k == 1 ? v : get(res, sel, row, k)[1] for k in 1:n]
        end
        out[row] = vals
    end
    return out
end
