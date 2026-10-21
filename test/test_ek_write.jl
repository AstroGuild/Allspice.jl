# Group EK 3: writing EK files.
#
# CSPICE is the oracle. A seeded random model is written by CSPICE and by Allspice with
# the same sequence of calls; after the close the three DAS address spaces of the two
# files (characters, doubles, integers) must be equal word for word. The same holds
# for random edit scripts (insert, delete, update, new segments) that either engine
# applies to a file that either engine wrote. CSPICE reads and queries the files that
# Allspice wrote, and every input check gives the error of CSPICE.
#
# The calls are qualified (`Allspice.ekopn`), although the writers are exported.

# The code of the EK groups is long and most of it runs once, so Julia compiles the
# module of the tests without optimization while this file runs (5 s less for the five
# EK groups); the last line of the file gives the default level back. The level of
# Allspice does not change.
Base.Experimental.@optlevel 0

isdefined(@__MODULE__, :EKLIB) || include("ekfixtures.jl")

# ---------------------------------------------------------------------------------------
# The two engines behind one set of signatures
# ---------------------------------------------------------------------------------------

"""
`ekaclc_c` with a string length of at least 129. `c_ekaclc` of the fixtures passes the
longest string of the column plus one, and CSPICE then refuses a fixed-count column
whose string length is greater (`zzekac09`), and can lose a character of an array
column (`zzekac06`).
"""
function ekw_c_ekaclc(h, segno, column, vals, entszs, nlflgs, rcptrs; minlen=129)
    len = max(c_vallen(vals), minlen)
    wk = zeros(Cint, length(rcptrs))
    ccall((:ekaclc_c, EKLIB), Cvoid, (Cint, Cint, Cstring, Cint, Ptr{UInt8}, Ptr{Cint},
            Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), h, segno - 1, column, len, c_strarray(vals, len),
        Cint[e for e in entszs], Cint[f for f in nlflgs], rcptrs, wk)
    c_check()
end

"`ekbseg_c` with room for names of 80 characters (`c_ekbseg` cuts a name after 32)."
function ekw_c_ekbseg(h, table, cnames, decls)
    segno = Ref{Cint}(0)
    ccall((:ekbseg_c, EKLIB), Cvoid, (Cint, Cstring, Cint, Cint, Ptr{UInt8}, Cint, Ptr{UInt8},
            Ref{Cint}), h, table, length(cnames), 81, c_strarray(cnames, 81), 201,
        c_strarray(decls, 201), segno)
    c_check()
    return segno[] + 1
end

const EKW_C = (ekopn=c_ekopn, ekopw=c_ekopw, ekopr=c_ekopr, ekcls=c_ekcls, ekbseg=ekw_c_ekbseg,
    ekappr=c_ekappr, ekinsr=c_ekinsr, ekdelr=c_ekdelr, ekacec=c_ekacec, ekaced=c_ekaced,
    ekacei=c_ekacei, ekucec=c_ekucec, ekuced=c_ekuced, ekucei=c_ekucei, ekifld=c_ekifld,
    ekaclc=ekw_c_ekaclc, ekacld=c_ekacld, ekacli=c_ekacli, ekffld=c_ekffld)

const EKW_N = (ekopn=Allspice.ekopn, ekopw=Allspice.ekopw, ekopr=ekopr, ekcls=ekcls,
    ekbseg=Allspice.ekbseg, ekappr=Allspice.ekappr, ekinsr=Allspice.ekinsr,
    ekdelr=Allspice.ekdelr, ekacec=Allspice.ekacec, ekaced=Allspice.ekaced,
    ekacei=Allspice.ekacei, ekucec=Allspice.ekucec, ekuced=Allspice.ekuced,
    ekucei=Allspice.ekucei, ekifld=Allspice.ekifld, ekaclc=Allspice.ekaclc,
    ekacld=Allspice.ekacld, ekacli=Allspice.ekacli, ekffld=Allspice.ekffld)

"""
One call `f(args...)` of an engine as a value. `c_kind` and `neo_kind` are compiled for
each type of their argument; with this type they are compiled once for all input checks.
"""
struct EKWCall
    f::Any
    args::Any
end
(c::EKWCall)() = c.f(c.args...)

"Factor of the row counts of the standard model: the whole model in the exhaustive run only."
const EKW_SCALE = FULL ? 1.0 : 0.25

"The three DAS address spaces of the file `path`."
function ekw_spaces(path)
    h = dasopr(path)
    s = (copy(h.chars), copy(h.dps), copy(h.ints))
    dascls(h)
    return s
end

"Number of words of the files `a` and `b` that differ (a length difference counts), and the word count of `a`."
function ekw_diff(a, b)
    sa = ekw_spaces(a)
    sb = ekw_spaces(b)
    nd = 0
    for k in 1:3
        nd += ekw_ndiff(sa[k], sb[k])       # a call for each element type, so the loop is typed
    end
    return nd, sum(length, sa)
end

"Number of elements of `x` and `y` that differ; a length difference counts."
function ekw_ndiff(x::AbstractVector, y::AbstractVector)
    n = min(length(x), length(y))
    nd = abs(length(x) - length(y))
    for i in 1:n
        nd += !isequal(x[i], y[i])
    end
    @assert 0 <= nd <= max(length(x), length(y))
    return nd
end

"Depths of the trees of the EK `path`: (segment tree, record trees, index trees, data page trees)."
function ekw_depths(path)
    h = ekopr(path)
    depth(tree) = tree > 0 ? Int(h.ints[Allspice._ek_ibase(tree)+Allspice._EK_TRDPTH]) : 0
    rec = 0
    idx = 0
    dat = [0, 0, 0]
    for segno in 1:eknseg(h)
        seg = Allspice._ek_segment(h, segno)
        rec = max(rec, depth(seg.segdsc[Allspice._EK_RTIDX]))
        for t in 1:3
            dat[t] = max(dat[t], depth(seg.segdsc[Allspice._EK_CPTIDX+t-1]))
        end
        for col in seg.columns
            col.ixtype == 1 && (idx = max(idx, depth(col.ixptr)))
        end
    end
    sgt = depth(Allspice._ek_segtree(h))
    ekcls(h)
    return (seg=sgt, rec=rec, idx=idx, chr=dat[1], dp=dat[2], int=dat[3])
end

ekw_depth_line(name, d) = "ek_write: tree depths of $name: segment $(d.seg), record $(d.rec), " *
                          "index $(d.idx), data pages chr $(d.chr), dp $(d.dp), int $(d.int)"

"""
Read every entry of the EK `path` with CSPICE and with Allspice and compare both with
the model `segs`; returns (entries compared, entries that differ).
"""
function ekw_read_check(path, segs::Vector{EKSegSpec})
    hc = c_ekopr(path)
    hn = ekopr(path)
    count = 0
    bad = 0
    try
        @test eknseg(hn) == c_eknseg(hc) == length(segs)
        for (segno, seg) in enumerate(segs)
            sc = c_ekssum(hc, segno)
            sn = ekssum(hn, segno)
            @test sn.tabnam == sc.tabnam == seg.table
            @test sn.nrows == sc.nrows == length(seg.rows)
            @test sn.cnames == sc.cnames == [c.name for c in seg.cols]
            @test sc.dtypes == [c.dtype for c in seg.cols]
            @test sc.sizes == [c.size for c in seg.cols]
            @test sc.indexd == [c.indexed for c in seg.cols]
            @test sc.nullok == [c.nullok for c in seg.cols]
            for (c, col) in enumerate(seg.cols), (r, row) in enumerate(seg.rows)
                got, ref = if col.dtype === :CHR
                    ekrcec(hn, segno, r, col.name), c_ekrcec(hc, segno, r, col.name;
                        maxn=col.size == -1 ? 720 : col.size,
                        maxlen=col.strlen == -1 ? 2600 : col.strlen + 1)
                elseif col.dtype === :INT
                    ekrcei(hn, segno, r, col.name), c_ekrcei(hc, segno, r, col.name)
                else
                    ekrced(hn, segno, r, col.name), c_ekrced(hc, segno, r, col.name)
                end
                count += 1
                (got == ref && ek_entry_equal(col, row[c], ref)) || (bad += 1)
            end
        end
    finally
        c_ekcls(hc)
        ekcls(hn)
    end
    return count, bad
end

# ---------------------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------------------

"Columns whose entries fill many data pages of each type, so the data page trees get deep."
const EKW_DEEP_COLS = [
    EKColSpec("K", :INT; indexed=true),
    EKColSpec("S", :CHR; strlen=-1, indexed=true, nullok=true),
    EKColSpec("V", :DP; indexed=true, nullok=true),
    EKColSpec("TXT", :CHR; strlen=-1),
    EKColSpec("IV", :INT; size=230),
    EKColSpec("DV", :DP; size=115)]

"A segment of `EKW_DEEP_COLS` with `nrows` rows; the strings of TXT have 880 to 900 characters."
function ekw_deep_segment(rng, table, nrows; fast=false)
    rows = Vector{Any}[]
    for r in 1:nrows
        push!(rows, Any[Any[rand(rng, 1:nrows)],
            rand(rng) < 0.1 ? nothing : Any[rand(rng, EK_WORDS)],
            rand(rng) < 0.1 ? nothing : Any[Float64(rand(rng, -500:500))],
            Any[join(rand(rng, 'a':'z', rand(rng, 880:900)))],
            Any[rand(rng, -9:9) for _ in 1:230],
            Any[Float64(rand(rng, -9:9)) for _ in 1:115]])
    end
    return EKSegSpec(table, EKW_DEEP_COLS, rows, fast)
end

"Columns of a narrow table with four indexes, for many rows and for queries."
const EKW_KEY_COLS = [
    EKColSpec("K", :INT; indexed=true),
    EKColSpec("S", :CHR; strlen=12, indexed=true, nullok=true),
    EKColSpec("V", :DP; indexed=true, nullok=true),
    EKColSpec("T", :TIME; indexed=true),
    EKColSpec("U", :CHR; strlen=-1, indexed=true)]

"The fixed-count form of `EKW_KEY_COLS` (a type 2 segment) without the column U."
const EKW_FIXED_COLS = [EKColSpec(c.name, c.dtype, c.strlen, c.size, c.indexed, c.nullok, true)
                        for c in EKW_KEY_COLS[1:4]]

# ---------------------------------------------------------------------------------------
# Edit scripts
# ---------------------------------------------------------------------------------------

"The entry `e` (`nothing` for null) as the (values, isnull) arguments of an `ekace*` call."
function ekw_args(col::EKColSpec, e)
    e === nothing || return e, false
    n = col.size == -1 ? 1 : col.size
    return fill(col.dtype === :CHR ? " " : col.dtype === :INT ? 0 : 0.0, n), true
end

"""
A random edit script for the segments `segs` and the model after the script. An
operation is `(:append, segno, row)`, `(:insert, segno, recno, row)`, `(:delete,
segno, recno)`, `(:update, segno, recno, c, entry)`, or `(:segment, spec)` for a new
segment. `weights` are the relative counts of append, insert, delete and update.

A segment that the fast writers made gets new rows only, unless `fastedits` is set:
CSPICE cannot delete its rows and cannot update all of its entries (see the group on
the edits of fast-loaded segments).
"""
function ekw_script(rng, segs::Vector{EKSegSpec}, nops::Int; weights=(1, 2, 2, 5),
        newsegs=EKSegSpec[], fastedits=false)
    model = [EKSegSpec(s.table, s.cols, copy(s.rows), s.fast) for s in segs]
    editable = [i for (i, s) in enumerate(model) if !any(c -> c.fixed, s.cols)]
    ops = Any[]
    total = sum(weights)
    for k in 1:nops
        if !isempty(newsegs) && k % max(1, nops ÷ (length(newsegs) + 1)) == 0 &&
           length(model) < length(segs) + length(newsegs)
            spec = newsegs[length(model)-length(segs)+1]
            push!(ops, (:segment, spec))
            push!(model, EKSegSpec(spec.table, spec.cols, copy(spec.rows), spec.fast))
            any(c -> c.fixed, spec.cols) || push!(editable, length(model))
            continue
        end
        segno = rand(rng, editable)
        seg = model[segno]
        n = length(seg.rows)
        x = rand(rng, 1:total)
        (seg.fast && !fastedits) && (x = min(x, weights[1] + weights[2]))
        newrow() = Any[ek_random_entry(rng, c, n + 1) for c in seg.cols]
        if x <= weights[1] || n == 0
            row = newrow()
            push!(ops, (:append, segno, row))
            push!(seg.rows, row)
        elseif x <= weights[1] + weights[2]
            recno = rand(rng, 1:n+1)
            row = newrow()
            push!(ops, (:insert, segno, recno, row))
            insert!(seg.rows, recno, row)
        elseif x <= weights[1] + weights[2] + weights[3]
            recno = rand(rng, 1:n)
            push!(ops, (:delete, segno, recno))
            deleteat!(seg.rows, recno)
        else
            recno = rand(rng, 1:n)
            c = rand(rng, 1:length(seg.cols))
            e = ek_random_entry(rng, seg.cols[c], recno)
            push!(ops, (:update, segno, recno, c, e))
            seg.rows[recno] = copy(seg.rows[recno])
            seg.rows[recno][c] = e
        end
    end
    return ops, model
end

"Write the segment `seg` to the open file `h` with the writers `w`; returns its number."
function ekw_write_segment(w, h, seg::EKSegSpec)
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
        return segno
    end
    segno = w.ekbseg(h, seg.table, cnames, decls)
    for row in seg.rows
        ekw_fill(w, h, segno, w.ekappr(h, segno), seg.cols, row)
    end
    return segno
end

"Write `segs` to the new file `path` with the writers `w` (the fast `ekaclc` gets wide strings)."
function ekw_write(w, path, segs::Vector{EKSegSpec})
    h = w.ekopn(path, "EK_TEST", 0)
    foreach(seg -> ekw_write_segment(w, h, seg), segs)
    w.ekcls(h)
    return path
end

"Add all entries of `row` to the new record `recno`."
function ekw_fill(w, h, segno, recno, cols, row)
    for (c, col) in enumerate(cols)
        f = col.dtype === :CHR ? w.ekacec : col.dtype === :INT ? w.ekacei : w.ekaced
        vals, isnull = ekw_args(col, row[c])
        f(h, segno, recno, col.name, vals, isnull)
    end
end

"Apply the script `ops` to the file `path` with the writers `w`; `cols[segno]` are the columns."
function ekw_apply(w, path, segs::Vector{EKSegSpec}, ops)
    cols = [s.cols for s in segs]
    h = w.ekopw(path)
    for op in ops
        if op[1] === :append
            ekw_fill(w, h, op[2], w.ekappr(h, op[2]), cols[op[2]], op[3])
        elseif op[1] === :insert
            w.ekinsr(h, op[2], op[3])
            ekw_fill(w, h, op[2], op[3], cols[op[2]], op[4])
        elseif op[1] === :delete
            w.ekdelr(h, op[2], op[3])
        elseif op[1] === :update
            col = cols[op[2]][op[4]]
            f = col.dtype === :CHR ? w.ekucec : col.dtype === :INT ? w.ekucei : w.ekuced
            vals, isnull = ekw_args(col, op[5])
            f(h, op[2], op[3], col.name, vals, isnull)
        else
            ekw_write_segment(w, h, op[2])
            push!(cols, op[2].cols)
        end
    end
    w.ekcls(h)
    return path
end

"""
Write `segs` with both engines, then let both engines apply `ops` to both files. All
four results must be equal word for word and must read back as `model`. Returns the
path of one result, the number of words compared and the number of entries read.
"""
function ekw_cross(dir, tag, segs, ops, model)
    base = Dict(e => ekw_write(w, joinpath(dir, "$(tag)_base_$e.bes"), segs)
                for (e, w) in (("c", EKW_C), ("n", EKW_N)))
    nd, words = ekw_diff(base["c"], base["n"])
    @test nd == 0
    outs = String[]
    for b in ("c", "n"), (e, w) in (("c", EKW_C), ("n", EKW_N))
        out = joinpath(dir, "$(tag)_$(b)_$e.bes")
        cp(base[b], out)
        push!(outs, ekw_apply(w, out, segs, ops))
    end
    for o in outs[2:end]
        nd, nw = ekw_diff(outs[1], o)
        nd == 0 || println("ek_write: $tag: $(basename(o)) differs from $(basename(outs[1])) in $nd words")
        @test nd == 0
        words += nw
    end
    count, bad = ekw_read_check(outs[4], model)
    @test bad == 0
    return outs[4], words, count
end

# ---------------------------------------------------------------------------------------
# (a) The same model with both engines, word for word
# ---------------------------------------------------------------------------------------

group("EK write: files of Allspice are equal to files of CSPICE word for word") do
    dir = mktempdir()
    words = 0
    files = 0
    entries = 0
    for (name, segs) in ek_standard_model(SEED; scale=EKW_SCALE)
        pc = ekw_write(EKW_C, joinpath(dir, "c_" * name), segs)
        pn = ekw_write(EKW_N, joinpath(dir, "n_" * name), segs)
        nd, nw = ekw_diff(pc, pn)
        @test nd == 0
        words += nw
        files += 1
        count, bad = ekw_read_check(pn, segs)
        @test bad == 0
        entries += count
        measured(ekw_depth_line(name, ekw_depths(pn)))
    end
    # An empty file and a file with an empty segment.
    for (k, segs) in enumerate((EKSegSpec[], [EKSegSpec("EMPTY", EK_EVENT_COLS, Vector{Any}[], false)]))
        pc = ekw_write(EKW_C, joinpath(dir, "c_empty$k.bes"), segs)
        pn = ekw_write(EKW_N, joinpath(dir, "n_empty$k.bes"), segs)
        nd, nw = ekw_diff(pc, pn)
        @test nd == 0
        words += nw
        files += 1
    end
    # A file with comment area and many segments, so that the segment tree is used.
    rng = MersenneTwister(SEED + 2)
    many = [ek_random_segment(rng, "T$(k % 7)", k % 3 == 0 ? EK_GROUP_COLS : EK_LINK_COLS,
        rand(rng, 1:12); fast=k % 3 == 0) for k in 1:120]
    paths = String[]
    for (e, w) in (("c", EKW_C), ("n", EKW_N))
        path = joinpath(dir, "$(e)_many.bes")
        h = w.ekopn(path, "MANY SEGMENTS", 3000)
        foreach(seg -> ekw_write_segment(w, h, seg), many)
        w.ekcls(h)
        push!(paths, path)
    end
    nd, nw = ekw_diff(paths...)
    @test nd == 0
    # Also the file record, the comment area and the directory are equal. CSPICE leaves
    # records that it does not use any more at the end of the file; Allspice does not.
    bc = read(paths[1])
    bn = read(paths[2])
    @test length(bn) <= length(bc)
    @test view(bc, 1:length(bn)) == bn
    words += nw
    files += 1
    count, bad = ekw_read_check(paths[2], many)
    @test bad == 0
    measured(ekw_depth_line("many.bes", ekw_depths(paths[2])))
    measured("ek_write: $files files written by both engines, $words words equal, " *
             "$(entries + count) entries read back by CSPICE and Allspice")
end

# ---------------------------------------------------------------------------------------
# Deep trees
# ---------------------------------------------------------------------------------------

group("EK write: trees of depth 3 (records, indexes, data pages)") do
    dir = mktempdir()
    rng = MersenneTwister(SEED + 3)
    nrows = 5600
    # CSPICE needs 0.6 ms for each row that `ekappr` and `ekace*` add to KEYS (4 s for the
    # segment) and 9 ms for each row of DEEP (51 s). The default run therefore compares
    # word for word only the segments of the fast writers, which have trees of depth 3
    # of every kind. The segments that grow row by row are compared word for word in the
    # exhaustive run. In the default run Allspice alone writes the KEYS segment row by
    # row, and both engines then edit and read it (below); the DEEP segment that grows
    # row by row is not written.
    deepslow = FULL ? [ekw_deep_segment(rng, "DEEP", nrows)] : EKSegSpec[]
    deep = [deepslow..., ekw_deep_segment(rng, "DEEP", nrows; fast=true)]
    keyslow = ek_random_segment(rng, "KEYS", EKW_KEY_COLS, 6000; keyrange=1:3000)
    keyfast = ek_random_segment(rng, "KEYS", EKW_KEY_COLS, 6000; fast=true, keyrange=1:3000)
    keys = [(FULL ? [keyslow] : EKSegSpec[])..., keyfast,
        ek_random_segment(rng, "FIXED", EKW_FIXED_COLS, 6000; fast=true, keyrange=1:3000)]
    words = 0
    for (name, segs) in (("deep.bes", deep), ("keys.bes", keys))
        pc = ekw_write(EKW_C, joinpath(dir, "c_" * name), segs)
        pn = ekw_write(EKW_N, joinpath(dir, "n_" * name), segs)
        nd, nw = ekw_diff(pc, pn)
        @test nd == 0
        words += nw
        d = ekw_depths(pn)
        measured(ekw_depth_line(name, d))
        @test d.rec == 3
        @test d.idx == 3
        name == "deep.bes" && @test (d.chr, d.dp, d.int) == (3, 3, 3)
        count, bad = ekw_read_check(pn, segs)
        @test bad == 0
    end
    # Deletes and inserts in the deep trees: both engines edit the file of Allspice.
    segs = [keyslow, keyfast]
    ops, model = ekw_script(rng, segs, nsample(500, 5000); weights=(0, 1, 6, 1))
    base = ekw_write(EKW_N, joinpath(dir, "keys12.bes"), segs)
    d = ekw_depths(base)
    @test d.rec == 3 && d.idx == 3
    outs = String[]
    for (e, w) in (("c", EKW_C), ("n", EKW_N))
        out = joinpath(dir, "keys12_$e.bes")
        cp(base, out)
        push!(outs, ekw_apply(w, out, segs, ops))
    end
    nd, nw = ekw_diff(outs...)
    @test nd == 0
    words += nw
    count, bad = ekw_read_check(outs[2], model)
    @test bad == 0
    measured(ekw_depth_line("keys after $(length(ops)) edits ($(length(model[1].rows)) and " *
                            "$(length(model[2].rows)) rows left)", ekw_depths(outs[2])))
    measured("ek_write: deep files, $words words equal")
end

# ---------------------------------------------------------------------------------------
# (b) CSPICE loads and queries the files of Allspice
# ---------------------------------------------------------------------------------------

group("EK write: CSPICE loads and queries files that Allspice wrote") do
    dir = mktempdir()
    rng = MersenneTwister(SEED + 4)
    k1 = ek_random_segment(rng, "KEYS", EKW_KEY_COLS, 900; keyrange=1:300)
    k2 = ek_random_segment(rng, "KEYS", EKW_KEY_COLS, 700; fast=true, keyrange=1:300)
    fx = ek_random_segment(rng, "FIXED", EKW_FIXED_COLS, 800; fast=true, keyrange=1:300)
    path = ekw_write(EKW_N, joinpath(dir, "query.bes"), [k1, k2, fx])
    nent, bad = ekw_read_check(path, [k1, k2, fx])
    @test bad == 0
    # The index searches of the write layer (`zzekillt`, `zzekille`) count as the model does.
    h = ekopr(path)
    searches = 0
    for (segno, seg) in enumerate((k1, k2, fx)), (c, key) in ((1, 150), (1, 0), (1, 301), (3, -3.0),
        (3, 2.5e3), (2, "beta"), (2, "a"), (2, "~"), (4, 0.0))
        _, mbase = Allspice._ek_mloc(h, segno)
        segdsc = h.ints[mbase+1:mbase+24]
        coldsc = h.ints[mbase+24+(c-1)*11+1:mbase+24+c*11]
        vals = [r[c][1] for r in seg.rows if r[c] !== nothing]
        nnull = length(seg.rows) - length(vals)
        key isa String && (vals = [String(rstrip(v)) for v in vals])
        @test Allspice._ekw_illt(h, segdsc, coldsc, length(seg.rows), key) == nnull + count(<(key), vals)
        @test Allspice._ekw_ille(h, segdsc, coldsc, length(seg.rows), key) == nnull + count(<=(key), vals)
        searches += 2
    end
    @test neo_kind(() -> Allspice._ekw_illt(h, h.ints[1:24], h.ints[1:11], 0, 1)) === :invalidsize
    ekcls(h)
    SPICE.kclear()
    SPICE.furnsh(path)
    queries = 0
    try
        @test c_ekntab() == 2
        for (table, rows) in (("KEYS", vcat(k1.rows, k2.rows)), ("FIXED", fx.rows))
            val(row, c) = row[c] === nothing ? nothing : row[c][1]
            # Equality on the integer index.
            for key in (17, 1000)                           # in the column, and not in it
                n, err, msg = c_ekfind("SELECT K, V FROM $table WHERE K = $key")
                @test !err
                got = c_ek_result(n, [:INT, :DP])
                want = [r for r in rows if val(r, 1) == key]
                @test n == length(want)
                @test sort([something(g[2], [1.0e99])[1] for g in got]) ==
                      sort([something(val(r, 3), 1.0e99) for r in want])
                queries += 1
            end
            # A range on the double index, in the order of the index.
            for (lo, hi) in ((-5.0, 5.0), (0.0, 1.0e9))     # a narrow range, and one without end
                n, err, msg = c_ekfind("SELECT V FROM $table WHERE V GT $lo AND V LE $hi ORDER BY V")
                @test !err
                got = [g[1][1] for g in c_ek_result(n, [:DP])]
                want = sort([val(r, 3) for r in rows if val(r, 3) !== nothing && lo < val(r, 3) <= hi])
                @test got == want
                queries += 1
            end
            # Equality and a range on the character index.
            for word in ("Beta gamma", "nothing")           # in the column, and not in it
                n, err, msg = c_ekfind("SELECT S, K FROM $table WHERE S EQ \"$word\"")
                @test !err
                @test n == count(r -> val(r, 2) !== nothing && rstrip(val(r, 2)[1:min(end, 12)]) == word, rows)
                queries += 1
            end
            n, err, msg = c_ekfind("SELECT S FROM $table WHERE S GE \"m\" ORDER BY S DESC")
            @test !err
            got = [g[1][1] for g in c_ek_result(n, [:CHR])]
            want = sort([String(rstrip(val(r, 2))) for r in rows if val(r, 2) !== nothing &&
                         String(rstrip(val(r, 2))) >= "m"]; rev=true)
            @test got == want
            # ORDER BY over all rows on the integer and on the time index; nulls come first.
            n, err, msg = c_ekfind("SELECT K FROM $table ORDER BY K")
            @test !err
            @test [g[1][1] for g in c_ek_result(n, [:INT])] == sort([val(r, 1) for r in rows])
            n, err, msg = c_ekfind("SELECT T FROM $table WHERE K LT 40 ORDER BY T DESC")
            @test !err
            @test [g[1][1] for g in c_ek_result(n, [:DP])] ==
                  sort([val(r, 4) for r in rows if val(r, 1) < 40]; rev=true)
            n, err, msg = c_ekfind("SELECT V FROM $table WHERE V IS NULL")
            @test !err
            @test n == count(r -> val(r, 3) === nothing, rows)
            queries += 4
        end
        # The variable-length string index of the table KEYS.
        rows = vcat(k1.rows, k2.rows)
        n, err, msg = c_ekfind("SELECT U FROM KEYS WHERE U BETWEEN \"b\" AND \"p\" ORDER BY U")
        @test !err
        strs = [String(rstrip(r[5][1])) for r in rows]
        @test [g[1][1] for g in c_ek_result(n, [:CHR])] == sort([s for s in strs if "b" <= s <= "p"])
        queries += 1
    finally
        SPICE.kclear()
    end
    measured("ek_write: CSPICE answered $queries queries on a file of Allspice, $nent entries " *
             "read, $searches index searches of the write layer")
end

# ---------------------------------------------------------------------------------------
# (c), (e) Edit scripts on files of both engines by both engines
# ---------------------------------------------------------------------------------------

group("EK write: edit scripts by both engines on files of both engines") do
    dir = mktempdir()
    rng = MersenneTwister(SEED + 5)
    words = 0
    entries = 0
    # Mixed edits of a segment with every column class and of a fast-loaded segment.
    segs = [ek_random_segment(rng, "EVENTS", EK_EVENT_COLS, nsample(120, 250)),
        ek_random_segment(rng, "EVENTS", EK_EVENT_COLS, nsample(100, 200); fast=true),
        ek_random_segment(rng, "GROUPS", EK_GROUP_COLS, 300; fast=true)]
    ops, model = ekw_script(rng, segs, nsample(400, 1500))
    out, nw, ne = ekw_cross(dir, "mixed", segs, ops, model)
    words += nw
    entries += ne
    measured(ekw_depth_line("mixed edits", ekw_depths(out)))
    # Reopen and append: new segments (both writers, both segment types) and new rows.
    segs = [ek_random_segment(rng, "LINKS", EK_LINK_COLS, 150)]
    newsegs = [ek_random_segment(rng, "EVENTS", EK_EVENT_COLS, 60),
        ek_random_segment(rng, "LINKS", EK_LINK_COLS, 90; fast=true),
        ek_random_segment(rng, "GROUPS", EK_GROUP_COLS, 70; fast=true),
        EKSegSpec("EVENTS", EK_EVENT_COLS, Vector{Any}[], false)]
    ops, model = ekw_script(rng, segs, nsample(150, 400); weights=(6, 1, 0, 1), newsegs=newsegs)
    @test length(model) == 5
    out, nw, ne = ekw_cross(dir, "append", segs, ops, model)
    words += nw
    entries += ne
    # Delete all rows, then add rows again; the pages come back from the free lists.
    segs = [ek_random_segment(rng, "EVENTS", EK_EVENT_COLS, nsample(90, 180)),
        ek_random_segment(rng, "LINKS", EK_LINK_COLS, nsample(110, 220))]
    ops = Any[]
    model = [EKSegSpec(s.table, s.cols, Vector{Any}[], s.fast) for s in segs]
    for (segno, seg) in enumerate(segs)
        n = length(seg.rows)
        while n > 0
            push!(ops, (:delete, segno, rand(rng, 1:n)))
            n -= 1
        end
    end
    path = Dict(e => ekw_apply(w, cp(ekw_write(EKW_N, joinpath(dir, "all_base_$e.bes"), segs),
            joinpath(dir, "all_$e.bes")), segs, ops) for (e, w) in (("c", EKW_C), ("n", EKW_N)))
    nd, nw = ekw_diff(path["c"], path["n"])
    @test nd == 0
    nent, bad = ekw_read_check(path["n"], model)
    @test bad == 0 && nent == 0
    readd, model2 = ekw_script(rng, model, nsample(200, 500); weights=(3, 3, 1, 3))
    out, nw, ne = ekw_cross(dir, "readd", segs, vcat(ops, readd), model2)
    words += nw
    entries += ne
    measured("ek_write: edit scripts, $words words equal over 4 engine combinations, " *
             "$entries entries read back")
end

# ---------------------------------------------------------------------------------------
# Edits that CSPICE cannot do
# ---------------------------------------------------------------------------------------

group("EK write: edits of fast-loaded segments and of cut strings, which fail in CSPICE") do
    dir = mktempdir()
    decl = ["DATATYPE = CHARACTER*(4), INDEXED = TRUE"]
    scol = [EKColSpec("S", :CHR; strlen=4, indexed=true)]
    # 1. The fast writer orders equal strings of an indexed column with a sort that is
    #    not stable; CSPICE then looks for a row at the wrong place of the index.
    function fastfile(w, name)
        h = w.ekopn(joinpath(dir, name), "X", 0)
        segno, rcptrs = w.ekifld(h, "T", 4, ["S"], decl)
        w.ekaclc(h, segno, "S", ["b", "c", "a", "b"], fill(1, 4), fill(false, 4), rcptrs)
        w.ekffld(h, segno, rcptrs)
        return h
    end
    hc = fastfile(EKW_C, "c_dup.bes")
    @test c_kind(() -> c_ekdelr(hc, 1, 1)) === :itemnotfound          # the defect of CSPICE
    @test c_kind(() -> c_ekucec(hc, 1, 4, "S", ["z"], false)) === :itemnotfound
    c_ekcls(hc)
    for del in 1:4
        hn = fastfile(EKW_N, "n_dup$del.bes")
        Allspice.ekucec(hn, 1, del, "S", ["bb"], false)
        Allspice.ekdelr(hn, 1, del)
        ekcls(hn)
        rows = Vector{Any}[Any[Any[s]] for s in deleteat!(["b", "c", "a", "b"], del)]
        nent, bad = ekw_read_check(joinpath(dir, "n_dup$del.bes"), [EKSegSpec("T", scol, rows, true)])
        @test bad == 0 && nent == 3
    end
    # 2. A string that is longer than the string length of the column goes into the index
    #    at the place of the whole string, but the column holds the cut string.
    function longfile(w, name)
        h = w.ekopn(joinpath(dir, name), "X", 0)
        segno = w.ekbseg(h, "T", ["S"], decl)
        w.ekacec(h, segno, w.ekappr(h, segno), "S", ["abcd"], false)
        w.ekacec(h, segno, w.ekappr(h, segno), "S", ["abcd"], false)
        w.ekucec(h, segno, 1, "S", ["abcdZ"], false)
        w.ekcls(h)
        return joinpath(dir, name)
    end
    pc = longfile(EKW_C, "c_long.bes")
    pn = longfile(EKW_N, "n_long.bes")
    @test ekw_diff(pc, pn)[1] == 0                  # the same index up to here
    hc = c_ekopw(pc)
    hn = Allspice.ekopw(pn)
    @test c_kind(() -> c_ekdelr(hc, 1, 2)) === :itemnotfound          # the defect of CSPICE
    Allspice.ekdelr(hn, 1, 2)
    c_ekcls(hc)
    ekcls(hn)
    nent, bad = ekw_read_check(pn, [EKSegSpec("T", scol, Vector{Any}[Any[Any["abcd"]]], false)])
    @test bad == 0 && nent == 1
    # 3. The record pages of a fast load have the link count 0. CSPICE frees such a page
    #    with the first row that is deleted, and new rows then take the place of the
    #    other rows of the page.
    kcol = [EKColSpec("K", :INT)]
    function recfile(w, name)
        path = ekw_write(w, joinpath(dir, name), [EKSegSpec("T", kcol,
            Vector{Any}[Any[Any[10]], Any[Any[20]], Any[Any[30]]], true)])
        h = w.ekopw(path)
        w.ekdelr(h, 1, 1)
        for k in (40, 50)
            w.ekacei(h, 1, w.ekappr(h, 1), "K", [k], false)
        end
        w.ekcls(h)
        hr = c_ekopr(path)
        got = try
            [c_ekrcei(hr, 1, r, "K")[1][1] for r in 1:4]
        catch
            c_pending_kind()
            Int[]
        end
        c_ekcls(hr)
        return got
    end
    @test recfile(EKW_C, "c_rec.bes") != [20, 30, 40, 50]             # the defect of CSPICE
    @test recfile(EKW_N, "n_rec.bes") == [20, 30, 40, 50]
    # 4. The fast writers of the array columns give a data page no link for an entry
    #    that continues on it. CSPICE frees such a page when the entries that start on
    #    it are replaced, and the next entry is then written over the end of another.
    acol = [EKColSpec("IV", :INT; size=200)]
    function arrfile(w, name)
        path = ekw_write(w, joinpath(dir, name), [EKSegSpec("T", acol,
            Vector{Any}[Any[Any[k for _ in 1:200]] for k in 1:3], true)])
        h = w.ekopw(path)
        w.ekucei(h, 1, 3, "IV", fill(33, 200), false)
        w.ekucei(h, 1, 1, "IV", fill(11, 200), false)
        w.ekcls(h)
        hr = c_ekopr(path)
        got = c_ekrcei(hr, 1, 2, "IV")[1]
        c_ekcls(hr)
        return got
    end
    @test arrfile(EKW_C, "c_arr.bes") != fill(2, 200)                 # the defect of CSPICE
    @test arrfile(EKW_N, "n_arr.bes") == fill(2, 200)
    # 5. Random edits by Allspice of fast-loaded segments: equal strings in indexed
    #    columns, array entries over several pages, all rows deleted and rows added
    #    again. CSPICE reads the result and uses the indexes.
    rng = MersenneTwister(SEED + 8)
    segs = [ek_random_segment(rng, "KEYS", EKW_KEY_COLS, 6000; fast=true, keyrange=1:60),
        ek_random_segment(rng, "EVENTS", EK_EVENT_COLS, 500; fast=true),
        ek_random_segment(rng, "LINKS", EK_LINK_COLS, 300; fast=true)]
    ops, model = ekw_script(rng, segs, nsample(1500, 5000); weights=(1, 2, 4, 6), fastedits=true)
    for recno in length(model[3].rows):-1:1
        push!(ops, (:delete, 3, rand(rng, 1:recno)))
    end
    empty!(model[3].rows)
    more, model = ekw_script(rng, model, 600; weights=(2, 2, 1, 3), fastedits=true)
    path = ekw_apply(EKW_N, ekw_write(EKW_N, joinpath(dir, "fastedit.bes"), segs), segs, vcat(ops, more))
    nent, bad = ekw_read_check(path, model)
    @test bad == 0
    measured(ekw_depth_line("fast-loaded segments after $(length(ops) + length(more)) edits",
        ekw_depths(path)))
    SPICE.kclear()
    SPICE.furnsh(path)
    try
        for (table, c, name) in (("KEYS", 2, "S"), ("KEYS", 5, "U"), ("LINKS", 3, "TAG"), ("EVENTS", 7, "LABEL"))
            rows = model[findfirst(s -> s.table == table, model)].rows
            n, err, msg = c_ekfind("SELECT $name FROM $table WHERE $name IS NOT NULL ORDER BY $name")
            @test !err
            @test [g[1][1] for g in c_ek_result(n, [:CHR])] ==
                  sort([String(rstrip(r[c][1])) for r in rows if r[c] !== nothing])
            n, err, msg = c_ekfind("SELECT $name FROM $table WHERE $name EQ \"alpha\"")
            @test !err
            @test n == count(r -> r[c] !== nothing && rstrip(r[c][1]) == "alpha", rows)
        end
        n, err, msg = c_ekfind("SELECT K FROM KEYS ORDER BY K")
        @test !err
        @test [g[1][1] for g in c_ek_result(n, [:INT])] == sort([r[1][1] for r in model[1].rows])
    finally
        SPICE.kclear()
    end
    measured("ek_write: Allspice applied $(length(ops) + length(more)) edits to fast-loaded " *
             "segments ($nent entries read back by CSPICE and Allspice)")
end

# ---------------------------------------------------------------------------------------
# (d) The fast writers, column classes 1 to 9
# ---------------------------------------------------------------------------------------

group("EK write: fast writers for the column classes 1 to 9") do
    dir = mktempdir()
    rng = MersenneTwister(SEED + 6)
    words = 0
    # Row counts around the page sizes, columns stored in another order than declared.
    for (k, nrows) in enumerate((1, 2, 14, 126, 127, 253, 254, 255, 1014, 1015, 2100))
        t1 = ek_random_segment(rng, "EVENTS", EK_EVENT_COLS, min(nrows, 300); fast=true)
        t2 = ek_random_segment(rng, "GROUPS", EK_GROUP_COLS, nrows; fast=true)
        t3 = ek_random_segment(rng, "KEYS", EKW_KEY_COLS, nrows; fast=true, keyrange=1:50)
        paths = String[]
        for (e, w) in (("c", EKW_C), ("n", EKW_N))
            path = joinpath(dir, "$(e)_fast$k.bes")
            h = w.ekopn(path, "FAST", 0)
            for seg in (t1, t2, t3)
                segno, rcptrs = w.ekifld(h, seg.table, length(seg.rows), [c.name for c in seg.cols],
                    [ekdecl(c) for c in seg.cols])
                for c in reverse(1:length(seg.cols))
                    col = seg.cols[c]
                    vals, entszs, nlflgs = ek_flat_column(seg, c)
                    f = col.dtype === :CHR ? w.ekaclc : col.dtype === :INT ? w.ekacli : w.ekacld
                    f(h, segno, col.name, vals, entszs, nlflgs, rcptrs)
                end
                w.ekffld(h, segno, rcptrs)
            end
            w.ekcls(h)
            push!(paths, path)
        end
        nd, nw = ekw_diff(paths...)
        nd == 0 || println("ek_write: fast load with $nrows rows differs in $nd words")
        @test nd == 0
        words += nw
        count, bad = ekw_read_check(paths[2], [t1, t2, t3])
        @test bad == 0
    end
    # All rows null, and no row null, in columns that allow nulls.
    for allnull in (true, false)
        cols = [EKColSpec(c.name, c.dtype, c.strlen, c.size, c.indexed, true, c.fixed) for c in EK_EVENT_COLS]
        fcols = [EKColSpec(c.name, c.dtype, c.strlen, c.size, c.indexed, true, c.fixed) for c in EK_GROUP_COLS]
        mk(cs, n) = EKSegSpec("N", cs, Vector{Any}[Any[allnull ? nothing :
                                                       something(ek_random_entry(rng, c, r), ek_random_entry(rng,
                    EKColSpec(c.name, c.dtype, c.strlen, c.size, false, false, c.fixed), r)) for c in cs] for r in 1:n], true)
        segs = [mk(cols, 40), mk(fcols, 300)]
        pc = ekw_write(EKW_C, joinpath(dir, "c_null$allnull.bes"), segs)
        pn = ekw_write(EKW_N, joinpath(dir, "n_null$allnull.bes"), segs)
        nd, nw = ekw_diff(pc, pn)
        @test nd == 0
        words += nw
        count, bad = ekw_read_check(pn, segs)
        @test bad == 0
    end
    # Two fast loads at the same time on two handles: every load has its own state.
    a = ek_random_segment(rng, "EVENTS", EK_EVENT_COLS, 90; fast=true)
    b = ek_random_segment(rng, "LINKS", EK_LINK_COLS, 130; fast=true)
    pa = joinpath(dir, "two_a.bes")
    pb = joinpath(dir, "two_b.bes")
    ha = Allspice.ekopn(pa, "A")
    hb = Allspice.ekopn(pb, "B")
    sa, ra = Allspice.ekifld(ha, a.table, length(a.rows), [c.name for c in a.cols], [ekdecl(c) for c in a.cols])
    sb, rb = Allspice.ekifld(hb, b.table, length(b.rows), [c.name for c in b.cols], [ekdecl(c) for c in b.cols])
    @test ra isa Allspice.EKFastLoad && ra isa AbstractVector{Int32} && length(ra) == 90
    for c in 1:max(length(a.cols), length(b.cols)), (h, s, r, seg) in ((ha, sa, ra, a), (hb, sb, rb, b))
        c <= length(seg.cols) || continue
        col = seg.cols[c]
        vals, entszs, nlflgs = ek_flat_column(seg, c)
        f = col.dtype === :CHR ? Allspice.ekaclc : col.dtype === :INT ? Allspice.ekacli : Allspice.ekacld
        f(h, s, col.name, vals, entszs, nlflgs, r)
    end
    @test neo_kind(() -> Allspice.ekffld(ha, sa, rb)) === :wrongsegment
    Allspice.ekffld(hb, sb, rb)
    Allspice.ekffld(ha, sa, ra)
    ekcls(ha)
    ekcls(hb)
    for (p, seg) in ((pa, a), (pb, b))
        nd, nw = ekw_diff(ekw_write(EKW_C, p * ".c", [seg]), p)
        @test nd == 0
        words += nw
    end
    measured("ek_write: fast writers, $words words equal")
end

# ---------------------------------------------------------------------------------------
# (f) Input checks
# ---------------------------------------------------------------------------------------

group("EK write: input checks give the errors of CSPICE") do
    dir = mktempdir()
    rng = MersenneTwister(SEED + 7)
    cols = [EKColSpec("I", :INT; indexed=true), EKColSpec("D", :DP; nullok=true),
        EKColSpec("C", :CHR; strlen=8, indexed=true), EKColSpec("IA", :INT; size=3),
        EKColSpec("DA", :DP; size=-1, nullok=true), EKColSpec("CA", :CHR; strlen=6, size=2),
        EKColSpec("T", :TIME)]
    seg = ek_random_segment(rng, "CHECKS", cols, 12)
    fseg = ek_random_segment(rng, "FIXED", EKW_FIXED_COLS, 9; fast=true)
    pc = ekw_write(EKW_C, joinpath(dir, "c_checks.bes"), [seg, fseg])
    pn = ekw_write(EKW_N, joinpath(dir, "n_checks.bes"), [seg, fseg])
    hc = c_ekopw(pc)
    hn = Allspice.ekopw(pn)
    checked = 0
    mismatch = String[]
    # Both engines must fail with the same kind. The call is `w.fname(h, args...)`, or
    # `w.fname(args...)` with `handle=false`. It is given as data and not as a closure:
    # the 100 closures of the checks took 2.3 s to compile.
    # `oracle=false` omits the call of CSPICE and takes `expect` as its answer: `ekacei_`
    # and `ekucei_` read the column descriptor through a segment descriptor that
    # `zzeksdsc_` did not set when the segment number is bad, and only then look at the
    # error state. The read uses the garbage of the stack as a DAS address, and its
    # search took 525 s in one run on Julia 1.10.
    function same(name, fname::Symbol, @nospecialize(args...); expect=nothing, handle::Bool=true,
            oracle::Bool=true)
        @assert oracle || expect !== nothing "a check without CSPICE needs the expected kind"
        kc = oracle ? c_kind(EKWCall(EKW_C[fname], handle ? (hc, args...) : args)) : expect
        kn = neo_kind(EKWCall(EKW_N[fname], handle ? (hn, args...) : args))
        (kc == kn && kc !== nothing) || push!(mismatch, "$name: CSPICE $kc, Allspice $kn")
        @test kc == kn
        @test kc !== nothing
        expect === nothing || @test kn === expect
        checked += 1
    end
    one = ["DATATYPE = INTEGER"]
    # ekbseg: names.
    same("table name too long", :ekbseg, "T"^65, ["A"], one; expect=:idstringtoolong)
    same("table name not printable", :ekbseg, "T\tX", ["A"], one; expect=:nonprintablechars)
    same("table name starts with digit", :ekbseg, "1T", ["A"], one; expect=:invalidname)
    same("table name with blank", :ekbseg, "T X", ["A"], one; expect=:invalidname)
    same("table name with leading blank", :ekbseg, " TX", ["A"], one; expect=:invalidname)
    same("table name with minus", :ekbseg, "T-X", ["A"], one; expect=:invalidname)
    same("table name empty", :ekbseg, "", ["A"], one; expect=:emptystring)
    same("column name too long", :ekbseg, "T", ["A"^33], one; expect=:idstringtoolong)
    same("column name invalid", :ekbseg, "T", ["A.B"], one; expect=:invalidname)
    same("column name blank", :ekbseg, "T", [" "], one; expect=:invalidname)
    same("101 columns", :ekbseg, "T", ["C$k" for k in 1:101], fill(one[1], 101);
        expect=:invalidcount)
    # ekbseg: declarations.
    for (name, decl, kind) in (
        ("no DATATYPE", "SIZE = 3", :badcolumdecl),
        ("blank declaration", " ", :badcolumdecl),
        ("DATATYPE without value", "INDEXED = TRUE, DATATYPE", :badcolumndecl),
        ("unknown type", "DATATYPE = REAL", :badcolumndecl),
        ("INTEGER with suffix", "DATATYPE = INTEGER*4", :badcolumndecl),
        ("CHARACTER without asterisk", "DATATYPE = CHARACTER(10)", :badcolumndecl),
        ("CHARACTER alone", "DATATYPE = CHARACTER", :badcolumndecl),
        ("CHARACTER length not a number", "DATATYPE = CHARACTER*(TEN)", :badcolumndecl),
        ("CHARACTER length missing", "DATATYPE = CHARACTER*()", :badcolumndecl),
        ("SIZE without value", "DATATYPE = INTEGER, SIZE", :badcolumndecl),
        ("SIZE not a number", "DATATYPE = INTEGER, SIZE = MANY", :badcolumndecl),
        ("VARIABLE size of CHARACTER*(*)", "DATATYPE = CHARACTER*(*), SIZE = VARIABLE", :badcolumndecl),
        ("NULLS_OK without value", "DATATYPE = TIME, NULLS_OK", :badcolumndecl),
        ("NULLS_OK not a boolean", "DATATYPE = TIME, NULLS_OK = YES", :badcolumndecl),
        ("INDEXED without value", "DATATYPE = TIME, INDEXED", :badcolumndecl),
        ("INDEXED not a boolean", "DATATYPE = TIME, INDEXED = 1", :badcolumndecl),
        ("index on array column", "DATATYPE = INTEGER, SIZE = 3, INDEXED = TRUE", :badcolumndecl),
        ("FIXED_COUNT without value", "DATATYPE = TIME, FIXED_COUNT", :badcolumndecl),
        ("FIXED_COUNT not a boolean", "DATATYPE = TIME, FIXED_COUNT = NO", :badcolumndecl),
        ("FIXED_COUNT on array column", "DATATYPE = DOUBLE PRECISION, SIZE = 2, FIXED_COUNT = TRUE",
            :badcolumndecl))
        same("declaration: $name", :ekbseg, "T", ["A"], [decl]; expect=kind)
    end
    same("fixed-count and other columns", :ekbseg, "T", ["A", "B"],
            ["DATATYPE = INTEGER", "DATATYPE = INTEGER, FIXED_COUNT = TRUE"]; expect=:badattributes)
    # Segment and row numbers.
    # The record calls check the segment number in `_ekw_segment` and the six entry
    # calls check segment, column and row in `_ekw_locate`: both bounds are checked
    # through one call of each kind, and one bad number through each of the others.
    for segno in (0, 3)
        same("ekappr segment $segno", :ekappr, segno; expect=:invalidindex)
        same("ekacei segment $segno", :ekacei, segno, 1, "I", [1], false; expect=:invalidindex, oracle=false)
    end
    same("ekinsr segment 3", :ekinsr, 3, 1; expect=:invalidindex)
    same("ekdelr segment 3", :ekdelr, 3, 1; expect=:invalidindex)
    same("ekucei segment 3", :ekucei, 3, 1, "I", [1], false; expect=:invalidindex, oracle=false)
    for recno in (0, 14)
        same("ekinsr row $recno", :ekinsr, 1, recno; expect=:invalidindex)
    end
    for recno in (0, 13)
        same("ekdelr row $recno", :ekdelr, 1, recno; expect=:invalidindex)
        same("ekacei row $recno", :ekacei, 1, recno, "IA", [1, 2, 3], false;
            expect=:indexoutofrange)
    end
    same("ekuced row 13", :ekuced, 1, 13, "D", [1.0], false; expect=:indexoutofrange)
    # Columns and data types: each call has its own data types.
    for (fname, bad) in ((:ekacec, "I"), (:ekaced, "C"), (:ekacei, "T"), (:ekucec, "D"),
        (:ekuced, "IA"), (:ekucei, "DA"))
        vals = fname in (:ekacec, :ekucec) ? ["x", "y"] : fname in (:ekaced, :ekuced) ? [1.0] : [1, 2, 3]
        same("$fname wrong type", fname, 1, 1, bad, vals, false; expect=:wrongdatatype)
    end
    for (fname, vals) in ((:ekacec, ["x", "y"]), (:ekucei, [1, 2, 3]))
        same("$fname unknown column", fname, 1, 1, "NOPE", vals, false; expect=:bug)
        same("$fname empty column name", fname, 1, 1, "", vals, false; expect=:emptystring)
    end
    # Entries.
    same("value in a full scalar entry (integer)", :ekacei, 1, 1, "I", [5], false;
        expect=:nonemptyentry)
    same("value in a full scalar entry (double)", :ekaced, 1, 1, "T", [5.0], false;
        expect=:nonemptyentry)
    both(f) = (f(EKW_C, hc); f(EKW_N, hn))
    both((w, h) -> w.ekappr(h, 1))
    same("null in a column without nulls (integer)", :ekacei, 1, 13, "I", [0], true;
        expect=:badattribute)
    same("null in a column without nulls (string)", :ekacec, 1, 13, "C", [" "], true;
        expect=:badattribute)
    same("null in a column without nulls (array)", :ekacei, 1, 13, "IA", [0, 0, 0], true;
        expect=:badattribute)
    same("null in a column without nulls (strings)", :ekacec, 1, 13, "CA", [" ", " "], true;
        expect=:badattribute)
    same("wrong entry size (integer)", :ekacei, 1, 13, "IA", [1, 2], false;
        expect=:invalidcount)
    same("wrong entry size (strings)", :ekacec, 1, 13, "CA", ["a", "b", "c"], false;
        expect=:invalidcount)
    # (fault fixed: an update that throws deleted the old entry first, as CSPICE does;
    # the entry is now as it was. A valid update then brings the two files together
    # again: CSPICE has deleted the old entry already, Allspice deletes it now.)
    old_ia = Allspice.ekrcei(hn, 1, 2, "IA")
    old_c = Allspice.ekrcec(hn, 1, 3, "C")
    old_ca = Allspice.ekrcec(hn, 1, 3, "CA")
    same("wrong entry size in update", :ekucei, 1, 2, "IA", [1, 2, 3, 4], false;
        expect=:invalidcount)
    same("update to null without nulls", :ekucec, 1, 3, "C", [" "], true;
        expect=:badattribute)
    same("wrong entry size in update (strings)", :ekucec, 1, 3, "CA", ["a"], false;
        expect=:invalidcount)
    @test Allspice.ekrcei(hn, 1, 2, "IA") == old_ia
    @test Allspice.ekrcec(hn, 1, 3, "C") == old_c
    @test Allspice.ekrcec(hn, 1, 3, "CA") == old_ca
    both((w, h) -> w.ekucei(h, 1, 2, "IA", [4, 5, 6], false))
    both((w, h) -> w.ekucec(h, 1, 3, "C", ["new"], false))
    both((w, h) -> w.ekucec(h, 1, 3, "CA", ["p", "q"], false))
    # (fault fixed: a second `ekacec`, or a second `ekacei`/`ekaced` of an array column,
    # wrote the entry again and left the pages of the first one in use, as CSPICE does;
    # the three calls now throw for every column class. Allspice only: CSPICE has no check.)
    for (fname, col, vals) in ((:ekacec, "C", ["x"]), (:ekacec, "CA", ["x", "y"]),
        (:ekacei, "IA", [1, 2, 3]), (:ekaced, "DA", [1.0]))
        @test neo_kind(EKWCall(EKW_N[fname], (hn, 1, 1, col, vals, false))) === :nonemptyentry
    end
    # Fast writers.
    same("ekifld without rows", :ekifld, "T", 0, ["A"], one; expect=:invalidcount)
    same("ekifld with bad name", :ekifld, "T T", 3, ["A"], one; expect=:invalidname)
    # After the errors both files hold the same words, and both handles still work.
    both((w, h) -> ekw_fill(w, h, 1, 13, cols, [ek_random_entry(MersenneTwister(5), c, 1) for c in cols]))
    both((w, h) -> w.ekdelr(h, 1, 4))
    c_ekcls(hc)
    ekcls(hn)
    nd, nw = ekw_diff(pc, pn)
    @test nd == 0
    # Read handles and closed handles.
    hc = c_ekopr(pc)
    hn = ekopr(pn)
    calls = Any[
        (:ekbseg, ("T", ["A"], one)),
        (:ekappr, (1,)),
        (:ekinsr, (1, 1)),
        (:ekdelr, (1, 1)),
        (:ekifld, ("T", 2, ["A"], one))]
    for (fname, args) in calls
        same("$fname with a read handle", fname, args...; expect=:dasinvalidaccess)
    end
    c_ekcls(hc)
    ekcls(hn)
    for (fname, args) in calls
        same("$fname with a closed handle", fname, args...; expect=:dasnosuchhandle)
    end
    same("ekopn with a negative comment count", :ekopn, joinpath(dir, "neg.bes"), "X", -1;
        expect=:invalidcount, handle=false)
    @test !ispath(joinpath(dir, "neg.bes"))
    same("ekopw of a missing file", :ekopw, joinpath(dir, "missing.bes"); handle=false)
    isempty(mismatch) || foreach(println, mismatch)

    # Checks of Allspice whose CSPICE counterpart differs or does not exist.
    path = joinpath(dir, "own.bes")
    h = Allspice.ekopn(path, "OWN")
    s1, r1 = Allspice.ekifld(h, "A", 3, ["K", "D"], ["DATATYPE = INTEGER", "DATATYPE = DOUBLE PRECISION"])
    @test neo_kind(() -> Allspice.ekacli(h, s1, "NOPE", [1, 2, 3], [1, 1, 1], falses(3), r1)) === :bug
    @test neo_kind(() -> Allspice.ekacli(h, s1, "D", [1, 2, 3], [1, 1, 1], falses(3), r1)) === :wrongdatatype
    @test neo_kind(() -> Allspice.ekacld(h, s1, "K", [1.0, 2, 3], [1, 1, 1], falses(3), r1)) === :wrongdatatype
    @test neo_kind(() -> Allspice.ekaclc(h, s1, "K", ["a", "b", "c"], [1, 1, 1], falses(3), r1)) === :wrongdatatype
    @test neo_kind(() -> Allspice.ekacli(h, s1, "K", [1, 2], [1, 1, 1], falses(3), r1)) === :invalidcount
    @test neo_kind(() -> Allspice.ekacli(h, s1 + 1, "K", [1, 2, 3], [1, 1, 1], falses(3), r1)) === :invalidindex
    Allspice.ekacli(h, s1, "K", [1, 2, 3], [1, 1, 1], falses(3), r1)
    @test neo_kind(() -> Allspice.ekffld(h, s1, r1)) === :invalidaddress     # a column is missing
    # (fault fixed: a column stored twice was not detected when every column was stored
    # before; `ekffld` of CSPICE then reads the pointers of one column for another)
    @test neo_kind(() -> Allspice.ekacli(h, s1, "K", [1, 2, 3], [1, 1, 1], falses(3), r1)) === :columnstoredtwice
    Allspice.ekacld(h, s1, "D", [1.0, 2.0, 3.0], [1, 1, 1], falses(3), r1)
    Allspice.ekffld(h, s1, r1)
    s2, r2 = Allspice.ekifld(h, "B", 2, ["K"], ["DATATYPE = INTEGER, FIXED_COUNT = TRUE"])
    @test neo_kind(() -> Allspice.ekacli(h, s2, "K", [1, 2], [1, 1], falses(2), r1)) === :wrongsegment
    @test neo_kind(() -> Allspice.ekffld(h, s2, r1)) === :wrongsegment
    @test neo_kind(() -> Allspice.ekffld(h, s2, r2)) === :invalidaddress     # fixed count: a column is missing
    Allspice.ekacli(h, s2, "K", [7, 8], Int[], Bool[], r2)
    @test neo_kind(() -> Allspice.ekacli(h, s2, "K", [7, 8], Int[], Bool[], r2)) === :columnstoredtwice
    Allspice.ekffld(h, s2, r2)
    # Records of a segment with fixed-count columns cannot be changed (CSPICE damages the file).
    @test neo_kind(() -> Allspice.ekappr(h, s2)) === :invalidtype
    @test neo_kind(() -> Allspice.ekinsr(h, s2, 1)) === :invalidtype
    @test neo_kind(() -> Allspice.ekdelr(h, s2, 1)) === :invalidtype
    @test neo_kind(() -> Allspice.ekacei(h, s2, 1, "K", [1], false)) === :indexoutofrange
    @test neo_kind(() -> Allspice.ekbseg(h, "T", ["A", "B"], one)) === :invalidcount
    ekcls(h)
    @test neo_kind(() -> Allspice.ekopn(path, "OWN")) === :fileopenfailed
    hr = c_ekopr(path)
    @test c_eknseg(hr) == 2
    @test c_ekrcei(hr, 1, 3, "K") == ([3], false)
    @test c_ekrced(hr, 1, 2, "D") == ([2.0], false)
    @test c_ekrcei(hr, 2, 2, "K") == ([8], false)
    # The refused record changes left the segment with fixed-count columns as it was.
    @test c_ekssum(hr, 2).nrows == 2 && c_ekrcei(hr, 2, 1, "K") == ([7], false)
    c_ekcls(hr)
    # (fault fixed: an entry size below 1 for a non-null row of a numeric SIZE = VARIABLE
    # column, below 0 for a character column, or below 0 for a null row threw a
    # BoundsError. CSPICE stops the process with BADSUBSCRIPT in `zzekac04` and
    # `zzekac05`, and signals OUTOFRANGE in `zzekac06`, so the kinds are not compared live.)
    hv = Allspice.ekopn(joinpath(dir, "sizes.bes"), "OWN")
    sv, rv = Allspice.ekifld(hv, "V", 3, ["I", "D", "C"],
        ["DATATYPE = INTEGER, SIZE = VARIABLE, NULLS_OK = TRUE",
            "DATATYPE = DOUBLE PRECISION, SIZE = VARIABLE", "DATATYPE = CHARACTER*(8), SIZE = VARIABLE"])
    words = (copy(hv.chars), copy(hv.dps), copy(hv.ints))
    @test neo_kind(() -> Allspice.ekacli(hv, sv, "I", [1, 2, 3, 4], [2, 0, 2], falses(3), rv)) === :badsubscript
    @test neo_kind(() -> Allspice.ekacld(hv, sv, "D", [1.0, 2, 3, 4], [2, -1, 3], falses(3), rv)) === :badsubscript
    @test neo_kind(() -> Allspice.ekaclc(hv, sv, "C", ["a", "b", "c", "d"], [2, -1, 3], falses(3), rv)) === :outofrange
    @test neo_kind(() -> Allspice.ekacli(hv, sv, "I", [1, 2, 3, 4], [-1, 2, 2], [true, false, false], rv)) === :invalidcount
    @test (hv.chars, hv.dps, hv.ints) == words
    # A null row of size 0 and an empty character entry are valid, as in CSPICE.
    Allspice.ekacli(hv, sv, "I", [1, 2, 3, 4], [2, 0, 2], [false, true, false], rv)
    Allspice.ekacld(hv, sv, "D", [1.0, 2, 3, 4], [2, 1, 1], falses(3), rv)
    Allspice.ekaclc(hv, sv, "C", ["a", "b", "c", "d"], [2, 0, 2], falses(3), rv)
    Allspice.ekffld(hv, sv, rv)
    @test Allspice.ekrcec(hv, sv, 2, "C") == (String[], false)
    ekcls(hv)
    measured("ek_write: $checked input checks with the error kind of CSPICE, " *
             "$(length(mismatch)) mismatches; $nw words equal after the failed calls")
end

# ---------------------------------------------------------------------------------------
# Entries that are not written yet
# ---------------------------------------------------------------------------------------

group("EK write: a new record has no entries until they are written") do
    dir = mktempdir()
    cols = [EKColSpec("I", :INT; indexed=true), EKColSpec("C", :CHR; strlen=-1, indexed=true),
        EKColSpec("DA", :DP; size=-1), EKColSpec("CA", :CHR; strlen=6, size=2)]
    rows = Vector{Any}[Any[Any[7], Any["seven"], Any[1.5, 2.5], Any["a", "b"]],
        Any[Any[9], Any["nine"], Any[3.5], Any["c", "d"]]]
    kinds = Dict{String,Any}()
    paths = Dict{String,String}()
    for (e, w, rdi, rdc) in (("c", EKW_C, c_ekrcei, c_ekrcec), ("n", EKW_N, ekrcei, ekrcec))
        kind = e == "c" ? c_kind : neo_kind
        paths[e] = joinpath(dir, "$(e)_new.bes")
        h = w.ekopn(paths[e], "NEW", 0)
        segno = w.ekbseg(h, "T", [c.name for c in cols], [ekdecl(c) for c in cols])
        for _ in 1:3
            w.ekappr(h, segno)
        end
        # A read of an entry that was never written is an error.
        kinds[e] = (kind(() -> rdi(h, segno, 1, "I")), kind(() -> rdc(h, segno, 1, "C")))
        # A row with such entries can be deleted, and the update calls write the entry.
        w.ekdelr(h, segno, 2)
        w.ekucei(h, segno, 1, "I", rows[1][1], false)
        w.ekucec(h, segno, 1, "C", rows[1][2], false)
        w.ekuced(h, segno, 1, "DA", rows[1][3], false)
        w.ekucec(h, segno, 1, "CA", rows[1][4], false)
        ekw_fill(w, h, segno, 2, cols, rows[2])
        w.ekcls(h)
    end
    @test kinds["n"] == kinds["c"] == (:uninitializedvalue, :uninitialized)
    # (regression: an integer outside `Int32` threw `InexactError`; it is a `SpiceError`
    # now, and the entry stays unwritten)
    h = ekopw(cp(paths["n"], joinpath(dir, "n_range.bes")))
    r = ekappr(h, 1)
    @test neo_kind(() -> ekacei(h, 1, r, "I", [2^31], false)) === :valueoutofrange
    @test neo_kind(() -> ekrcei(h, 1, r, "I")) === :uninitializedvalue
    dascls(h)
    @test ekw_diff(paths["c"], paths["n"])[1] == 0
    nent, bad = ekw_read_check(paths["n"], [EKSegSpec("T", cols, rows, false)])
    @test bad == 0 && nent == 8
end

# ---------------------------------------------------------------------------------------
# Declarations that are valid but unusual
# ---------------------------------------------------------------------------------------

group("EK write: column declarations are read as CSPICE reads them") do
    dir = mktempdir()
    decls = ["datatype = integer", "DATATYPE=DOUBLE PRECISION,NULLS_OK=TRUE",
        "DATATYPE = double   precision , INDEXED = true", "DATATYPE = DOUBLEPRECISION",
        "DATATYPE = T I M E, INDEXED = FALSE, NULLS_OK = FALSE",
        "DATATYPE = CHARACTER * ( 10 ), INDEXED = TRUE", "DATATYPE = character*(*)",
        "DATATYPE = CHARACTER*20, SIZE = VARIABLE", "DATATYPE = CHARACTER*(5), SIZE = 4",
        "SIZE = 3, DATATYPE = INTEGER", "NULLS_OK = TRUE, SIZE = variable, DATATYPE = TIME",
        "DATATYPE = INTEGER, FIXED_COUNT = FALSE", "DATATYPE = INTEGER, SIZE = 1, INDEXED = TRUE",
        "DATATYPE = INTEGER, SIZE = 2.9", "DATATYPE = CHARACTER*(1.2E1)",
        "DATATYPE = INTEGER, SIZE = 1, SIZE = 2", "DATATYPE = INTEGER,, SIZE = 2",
        "DATATYPE = TIME, UNKNOWN = 3", "DATATYPE = DOUBLE PRECISION, SIZE = -1"]
    fixed = ["DATATYPE = INTEGER, FIXED_COUNT = TRUE, INDEXED = TRUE",
        "FIXED_COUNT = TRUE, DATATYPE = CHARACTER*(7), NULLS_OK = TRUE",
        "DATATYPE = TIME, FIXED_COUNT = true"]
    paths = String[]
    for (e, w) in (("c", EKW_C), ("n", EKW_N))
        path = joinpath(dir, "$(e)_decl.bes")
        h = w.ekopn(path, "DECL", 0)
        w.ekbseg(h, "mixed_Case\$1", ["Col_$k" for k in 1:length(decls)], decls)
        w.ekbseg(h, "F", ["x$k" for k in 1:length(fixed)], fixed)
        w.ekbseg(h, "WIDE", ["C$k" for k in 1:100], [decls[1+k%length(decls)] for k in 1:100])
        w.ekcls(h)
        push!(paths, path)
    end
    nd, nw = ekw_diff(paths...)
    @test nd == 0
    h = ekopr(paths[2])
    s = ekssum(h, 1)
    @test s.tabnam == "MIXED_CASE\$1"
    @test s.cnames[3] == "COL_3"
    @test [a.cclass for a in s.cdescrs[1:9]] == [1, 2, 2, 2, 2, 3, 3, 6, 6]
    @test s.cdescrs[14].size == 2 && s.cdescrs[15].strlen == 12
    @test [a.cclass for a in ekssum(h, 2).cdescrs] == [7, 9, 8]
    @test ekssum(h, 3).ncols == 100
    ekcls(h)
    measured("ek_write: $(length(decls) + length(fixed)) unusual declarations, $nw words equal")
end

group("EK write: a scratch EK (ekops) holds what the file of ekopn holds") do
    dir = mktempdir()
    before = readdir(dir)
    words = 0
    for (name, segs) in ek_standard_model(SEED; scale=EKW_SCALE)
        h = ekops()
        foreach(seg -> ekw_write_segment(EKW_N, h, seg), segs)
        @test eknseg(h) == length(segs)
        f = ekopr(ekw_write(EKW_N, joinpath(dir, "f_" * name), segs))
        @test h.ints == f.ints && h.dps == f.dps && h.chars == f.chars
        words += length(h.ints) + length(h.dps)
        # The read calls and the record writers take the scratch handle.
        @test ekssum(h, 1).tabnam == ekssum(f, 1).tabnam
        ekcls(f)
        ekcls(h)
        @test throws_kind(() -> eknseg(h), :dasnosuchhandle)
        ekcls(h)                                    # a closed handle: nothing to do
    end
    @test Allspice.dashfn(ekops()) == ""
    @test all(startswith("f_"), setdiff(readdir(dir), before))      # no file of the scratch EK
    measured("ek_write: scratch EK, $words words equal to the files")
end

Base.Experimental.@optlevel -1
