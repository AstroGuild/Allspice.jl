# Group EK 1: read access to an EK file through an owned handle.
#
# Feature file: features/ek_files.feature
#
# CSPICE is the oracle. The EK files are generated with the CSPICE writers (`ekbseg`,
# `ekappr`, `ekacec`/`ekaced`/`ekacei`, and the fast writers) from a seeded random model
# with every column class and data type, null entries, and entries that continue over
# several data pages. Then `eknseg`, `ekssum` and every column entry are read with both
# engines and compared with each other and with the model.

isdefined(@__MODULE__, :EKLIB) || include("ekfixtures.jl")

"Read every entry of the EK `path` with both engines; returns the number of entries compared."
function ek_compare_file(path, segs::Vector{EKSegSpec})
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
            @test sn.ncols == sc.ncols == length(seg.cols)
            @test sn.cnames == sc.cnames
            @test [a.dtype for a in sn.cdescrs] == sc.dtypes == [c.dtype for c in seg.cols]
            @test [a.size for a in sn.cdescrs] == sc.sizes
            @test [a.strlen for a in sn.cdescrs] == sc.strlns
            @test [a.indexd for a in sn.cdescrs] == sc.indexd == [c.indexed for c in seg.cols]
            @test [a.nullok for a in sn.cdescrs] == sc.nullok == [c.nullok for c in seg.cols]
            for (c, col) in enumerate(seg.cols), (r, row) in enumerate(seg.rows)
                got, ref = if col.dtype === :CHR
                    ekrcec(hn, segno, r, col.name), c_ekrcec(hc, segno, r, col.name;
                        maxn=col.size == -1 ? 720 : col.size, maxlen=col.strlen == -1 ? 2600 : col.strlen + 1)
                elseif col.dtype === :INT
                    ekrcei(hn, segno, r, col.name), c_ekrcei(hc, segno, r, col.name)
                else
                    ekrced(hn, segno, r, col.name), c_ekrced(hc, segno, r, col.name)
                end
                count += 1
                (got == ref && ek_entry_equal(col, row[c], got)) || (bad += 1)
            end
        end
    finally
        c_ekcls(hc)
        ekcls(hn)
    end
    @test bad == 0
    return count
end

group("EK files: segments, summaries and column entries through a handle") do
    dir = mktempdir()
    model = ek_standard_model(SEED)
    total = 0
    for (name, segs) in model
        path = c_ek_write(joinpath(dir, name), segs)
        total += ek_compare_file(path, segs)
    end
    measured("  ek_file: $total column entries read with both engines")
    @test total > 30_000

    @testset "input checks" begin
        name, segs = model[1]
        path = joinpath(dir, name)
        hn = ekopr(path)
        hc = c_ekopr(path)
        ck(f) = c_kind(f)
        nk(f) = neo_kind(f)
        @test nk(() -> ekssum(hn, 0)) == ck(() -> c_ekssum(hc, 0)) == :indexoutofrange
        @test nk(() -> ekssum(hn, 99)) == ck(() -> c_ekssum(hc, 99))
        @test nk(() -> ekrcei(hn, 9, 1, "ID")) == ck(() -> c_ekrcei(hc, 9, 1, "ID"))
        @test nk(() -> ekrcei(hn, 1, 1, "NAME")) == ck(() -> c_ekrcei(hc, 1, 1, "NAME")) == :wrongdatatype
        @test nk(() -> ekrced(hn, 1, 1, "ID")) == ck(() -> c_ekrced(hc, 1, 1, "ID")) == :wrongdatatype
        @test nk(() -> ekrcec(hn, 1, 1, "VAL")) == ck(() -> c_ekrcec(hc, 1, 1, "VAL"; maxn=4, maxlen=40)) == :wrongdatatype
        @test nk(() -> ekrcei(hn, 1, 1, "NOSUCH")) == ck(() -> c_ekrcei(hc, 1, 1, "NOSUCH")) == :bug
        @test nk(() -> ekrcei(hn, 1, 0, "ID")) == ck(() -> c_ekrcei(hc, 1, 0, "ID")) == :indexoutofrange
        @test nk(() -> ekrcei(hn, 1, 10^6, "ID")) == ck(() -> c_ekrcei(hc, 1, 10^6, "ID"))
        # Column names are compared without regard to case and blanks.
        @test ekrcei(hn, 1, 1, " i d ") == ekrcei(hn, 1, 1, "ID")
        @test ekrced(hn, 1, 2, "epoch") == c_ekrced(hc, 1, 2, "epoch")
        # A row number outside a type 2 segment: CSPICE reads whatever is there.
        @test nk(() -> ekrcei(hn, 2, 10^6, "GID")) == :indexoutofrange
        c_ekcls(hc)
        ekcls(hn)
        ekcls(hn)
        @test nk(() -> eknseg(hn)) == :dasnosuchhandle
        @test nk(() -> ekrcei(hn, 1, 1, "ID")) == :dasnosuchhandle
        @test nk(() -> ekopr(joinpath(dir, "missing.bes"))) == :filenotfound
        # A DAS file that is not an EK.
        other = joinpath(dir, "other.das")
        h = dasonw(other, "TEST", "x", 0)
        dasadi(h, Int32[1, 2, 3])
        dascls(h)
        @test nk(() -> ekopr(other)) == ck(() -> c_ekopr(other)) == :wrongarchitecture
    end
end
