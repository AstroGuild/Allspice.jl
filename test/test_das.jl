# Group 11: DAS files and the DLA segment list.
#
# Feature file: features/das_files.feature
#
# CSPICE is the oracle. No DAS file is needed on the machine: the group makes its files
# with the CSPICE writers and with the Allspice writers in a temporary directory, and
# each engine reads the files of the other.
#
#   (a) a DSK file that CSPICE wrote: file record, summary, comments, and every integer,
#       double and character, read with both engines;
#   (b) a generic DAS file that Allspice wrote, read with CSPICE;
#   (c) updates and comment changes through `dasopw`;
#   (d) a file that is not segregated (CSPICE closes it with `daswbr_`, `dasllc_`);
#   (e) files of the other byte order, made in the test by a byte swap of each record;
#   (f) the DLA segment list in both directions;
#   (g) every input check, with the error kind of CSPICE where CSPICE has one.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

isdefined(@__MODULE__, :DSKLIB) || include("dskfixtures.jl")

function c_dasonw(path, ftype, ifname, ncomr)
    h = Ref{Cint}(0)
    ccall((:dasonw_, DSKLIB), Cint, (Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}, Ref{Cint}, Ref{Cint},
            Cint, Cint, Cint), path, ftype, ifname, ncomr, h, sizeof(path), sizeof(ftype),
        sizeof(ifname))
    c_check()
    return h[]
end
c_dasadi(h, data::Vector{Cint}) = (ccall((:dasadi_, DSKLIB), Cint, (Ref{Cint}, Ref{Cint},
        Ptr{Cint}), h, length(data), data); c_check())
c_dasadd(h, data::Vector{Float64}) = (ccall((:dasadd_, DSKLIB), Cint, (Ref{Cint}, Ref{Cint},
        Ptr{Float64}), h, length(data), data); c_check())
c_dasadc(h, data::String) = (ccall((:dasadc_, DSKLIB), Cint, (Ref{Cint}, Ref{Cint}, Ref{Cint},
        Ref{Cint}, Ptr{UInt8}, Cint), h, sizeof(data), 1, sizeof(data), data, sizeof(data)); c_check())

"Close the CSPICE write handle `h` without segregation of the file."
function c_close_unsegregated(h)
    ccall((:daswbr_, DSKLIB), Cint, (Ref{Cint},), h)
    c_check()
    ccall((:dasllc_, DSKLIB), Cint, (Ref{Cint},), h)
    c_check()
end

"""
Copy of the native DAS file `src` in the other byte order at `dst`: the integers of the
file record, every directory record and every integer and double data record are
swapped, and the format string names the other order. The directory chain is followed
in the test, without Allspice.
"""
function byteswapped_das(src, dst)
    raw = read(src)
    out = copy(raw)
    i32(off) = reinterpret(Int32, raw[off+1:off+4])[1]
    swap!(off, n, w) = for k in 0:n-1
        out[off+k*w+1:off+k*w+w] = reverse(raw[off+k*w+1:off+k*w+w])
    end
    swap!(68, 4, 4)
    native = String(raw[85:92])
    out[85:92] = codeunits(native == "LTL-IEEE" ? "BIG-IEEE" : "LTL-IEEE")
    rec = i32(68) + i32(76) + 2
    next = (2, 3, 1)
    prev = (3, 1, 2)
    while rec > 0                        # invariant: the records before `rec` are converted
        base = (rec - 1) * 1024
        swap!(base, 256, 4)
        typ = i32(base + 32)
        datarec = rec + 1
        prvtyp = typ in 1:3 ? prev[typ] : 0
        pos = 10
        while typ in 1:3 && pos <= 256 && i32(base + 4 * (pos - 1)) != 0
            count = i32(base + 4 * (pos - 1))
            cur = count > 0 ? next[prvtyp] : prev[prvtyp]
            for r in datarec:datarec+abs(count)-1
                cur == 2 && swap!((r - 1) * 1024, 128, 8)
                cur == 3 && swap!((r - 1) * 1024, 256, 4)
            end
            datarec += abs(count)
            prvtyp = cur
            pos += 1
        end
        rec = i32(base + 4)
    end
    write(dst, out)
    return dst
end

group("11. DAS files and DLA segment lists") do
    rng = MersenneTwister(SEED + 11)
    dir = mktempdir()
    full = (-pi, pi, -pi / 2, pi / 2)
    ev, ep = ellipsoid_shape(3.0, 2.5, 2.0, 24, 12)
    cuv, cup = cube_shape(1.0)
    specs = [DskSpec(499, 1, 2, "IAU_MARS", 1, Float64[], full, -1e9, 1e9, 5.0, 4, ev, ep),
        DskSpec(499, 2, 2, "IAU_MARS", 3, Float64[], (-1.0, 1.0, -1.0, 1.0), 0.0, 1e3, 1.0, 1, cuv, cup),
        DskSpec(401, 7, 1, "IAU_PHOBOS", 1, Float64[], full, -5.0, 5.0, 4.0, 3, ev, ep)]
    comments = ["First comment line", "", "  indented line with trailing blanks   ",
        "x"^300, "last line"]
    cfile = c_write_dsk(joinpath(dir, "c.bds"), specs; ifname="internal name", ncomch=5000,
        comments)

    @apart @testset "(a) a file that CSPICE wrote" begin
        h = dasopr(cfile)
        ch = c_dasopr(cfile)
        @test isopen(h)
        @test dashfn(h) == cfile
        @test dasrfr(h) == SPICE.dasrfr(ch)
        @test dasrfr(h)[1] == "DAS/DSK" && dasrfr(h)[2] == "internal name"
        @test daslla(h) == c_daslla(ch)
        @test dashfs(h) == c_dashfs(ch)
        @test dasec(h) == SPICE.dasec(ch)
        @test dasec(h) == rstrip.(comments)
        lastc, lastd, lasti = daslla(h)
        @test lastc == 0 && lastd > 0 && lasti > 0
        @test dasrdi(h, 1, lasti) == c_dasrdi(ch, 1, lasti)
        @test bits_equal(dasrdd(h, 1, lastd), c_dasrdd(ch, 1, lastd))
        for _ in 1:200
            a = rand(rng, 1:lasti); b = rand(rng, a:min(lasti, a + 700))
            @test dasrdi(h, a, b) == c_dasrdi(ch, a, b)
            a = rand(rng, 1:lastd); b = rand(rng, a:min(lastd, a + 700))
            @test bits_equal(dasrdd(h, a, b), c_dasrdd(ch, a, b))
        end
        @test dasrdi(h, 5, 4) == Int[] && dasrdd(h, 5, 4) == Float64[] && dasrdc(h, 5, 4) == ""
        c_dascls(ch)
        dascls(h)
        @test !isopen(h)
        dascls(h)                        # a second close does nothing
        @test getfat(cfile) == ("DAS", "DSK") == SPICE.getfat(cfile)
        @test Allspice.kernel_arch(cfile) === :DSK
    end

    @apart @testset "(b) a generic DAS file that Allspice wrote" begin
        path = joinpath(dir, "generic.das")
        h = dasonw(path, "TEST", "generic file", 2)
        @test isfile(path)
        ints = Cint[]
        dps = Float64[]
        chars = ""
        for _ in 1:40                    # appends of random sizes, the three types mixed
            a = rand(rng, Cint, rand(rng, 0:700)); append!(ints, a); dasadi(h, a)
            b = randn(rng, rand(rng, 0:400)); append!(dps, b); dasadd(h, b)
            c = randstring(rng, rand(rng, 0:900)); chars *= c; dasadc(h, c)
        end
        push!(dps, NaN, -0.0, Inf, floatmin(Float64), 5e-324)
        dasadd(h, dps[end-4:end])
        dasac(h, comments)
        @test daslla(h) == (sizeof(chars), length(dps), length(ints))
        @test dasrdi(h, 1, length(ints)) == ints
        @test dasrdc(h, 1, sizeof(chars)) == chars
        summary = dashfs(h)
        dascls(h)
        ch = c_dasopr(path)
        @test SPICE.dasrfr(ch) == ("DAS/TEST", "generic file", 0, 0, 2, sum(length, rstrip.(comments)) + 5)
        @test c_dashfs(ch) == summary
        @test c_dasrdi(ch, 1, length(ints)) == ints
        @test bits_equal(c_dasrdd(ch, 1, length(dps)), dps)
        @test c_dasrdc(ch, 1, sizeof(chars)) == chars
        @test SPICE.dasec(ch) == rstrip.(comments)
        c_dascls(ch)
        @test SPICE.getfat(path) == ("DAS", "TEST") == getfat(path)
        @test neo_das_data(path) == (ints, dps, chars) || bits_equal(neo_das_data(path)[2], dps)

        # An empty file, and a file with one data type only.
        e = joinpath(dir, "empty.das")
        dascls(dasonw(e, "TEST", "", 0))
        ch = c_dasopr(e)
        @test c_daslla(ch) == (0, 0, 0)
        @test c_dashfs(ch) == (h2 = dasopr(e); s = dashfs(h2); dascls(h2); s)
        c_dascls(ch)
        o = joinpath(dir, "onlydp.das")
        h = dasonw(o, "TEST", "doubles", 0)
        dasadd(h, dps[1:129])
        s = dashfs(h)
        dascls(h)
        ch = c_dasopr(o)
        @test c_dashfs(ch) == s
        @test bits_equal(c_dasrdd(ch, 1, 129), dps[1:129])
        c_dascls(ch)
    end

    @apart @testset "(c) updates and comments through dasopw" begin
        path = joinpath(dir, "update.das")
        cp(joinpath(dir, "generic.das"), path)
        ints0, dps0, chars0 = c_das_data(path)
        h = dasopw(path)
        dasudi(h, 10, 14, Cint[1, 2, 3, 4, 5])
        dasudd(h, 3, 4, [1.5, -2.5])
        dasudi(h, 7, 6, Cint[])           # an empty range changes nothing
        dasadi(h, Cint[77, 78])
        dasac(h, ["one more line"])
        dascls(h)
        ints1, dps1, chars1 = c_das_data(path)
        @test ints1 == vcat(ints0[1:9], Cint[1, 2, 3, 4, 5], ints0[15:end], Cint[77, 78])
        @test bits_equal(dps1, vcat(dps0[1:2], [1.5, -2.5], dps0[5:end]))
        @test chars1 == chars0
        ch = c_dasopr(path)
        @test SPICE.dasec(ch) == vcat(rstrip.(comments), "one more line")
        c_dascls(ch)
        # CSPICE changes the file; Allspice reads the change. Then the comments are deleted.
        ch = c_dasopw(path)
        SPICE.dasac(ch, ["from CSPICE"])
        c_dascls(ch)
        h = dasopw(path)
        @test dasec(h)[end] == "from CSPICE"
        dasdc(h)
        @test dasec(h) == String[]
        dascls(h)
        ch = c_dasopr(path)
        @test SPICE.dasrfr(ch)[5:6] == (0, 0)
        @test c_dasrdi(ch, 1, length(ints1)) == ints1
        c_dascls(ch)
        # Comments that need new comment records.
        many = ["line $i" * " y"^(i % 70) for i in 1:400]
        h = dasopw(path)
        dasac(h, many)
        dascls(h)
        ch = c_dasopr(path)
        @test SPICE.dasec(ch; bufsiz=500) == many
        @test c_das_data(path)[1] == ints1
        c_dascls(ch)
        h = dasopr(path)
        @test dasec(h) == many
        dascls(h)
    end

    unseg = joinpath(dir, "unseg.das")
    uints = Cint[]
    udps = Float64[]
    uchars = ""
    @apart @testset "(d) a file that is not segregated" begin
        ch = c_dasonw(unseg, "TEST", "not segregated", 1)
        for _ in 1:120                   # 360 clusters: more than one directory record holds
            a = rand(rng, Cint, rand(rng, 1:600)); append!(uints, a); c_dasadi(ch, a)
            b = randn(rng, rand(rng, 1:300)); append!(udps, b); c_dasadd(ch, b)
            c = randstring(rng, rand(rng, 1:2000)); uchars *= c; c_dasadc(ch, c)
        end
        c_close_unsegregated(ch)
        @test c_das_data(unseg) == (uints, udps, uchars)      # CSPICE reads its own file
        h = dasopr(unseg)
        @test daslla(h) == (sizeof(uchars), length(udps), length(uints))
        @test dasrdi(h, 1, length(uints)) == uints
        @test bits_equal(dasrdd(h, 1, length(udps)), udps)
        @test dasrdc(h, 1, sizeof(uchars)) == uchars
        ch = c_dasopr(unseg)
        @test dashfs(h) == c_dashfs(ch)
        c_dascls(ch)
        dascls(h)
        # Opened for write and closed, the file is segregated and holds the same data.
        seg = joinpath(dir, "reseg.das")
        cp(unseg, seg)
        dascls(dasopw(seg))
        @test c_das_data(seg) == (uints, udps, uchars)
        @test filesize(seg) < filesize(unseg)
    end

    @apart @testset "(e) the other byte order" begin
        for (name, src) in (("segregated DSK", cfile), ("not segregated", unseg),
                ("generic", joinpath(dir, "generic.das")))
            swapped = byteswapped_das(src, joinpath(dir, "swapped_" * basename(src)))
            @test read(swapped) != read(src)
            ci, cd, cc = c_das_data(src)
            @test c_das_data(swapped)[1] == ci           # CSPICE accepts the converted file
            ni, nd, nc = neo_das_data(swapped)
            @test ni == ci
            @test bits_equal(nd, cd)
            @test nc == cc
            h = dasopr(swapped)
            ch = c_dasopr(swapped)
            @test dasrfr(h) == SPICE.dasrfr(ch)
            @test dashfs(h) == c_dashfs(ch)
            @test dasec(h) == SPICE.dasec(ch; bufsiz=500)
            c_dascls(ch)
            dascls(h)
        end
        # A DSK of the other byte order can be loaded and searched.
        swapped = joinpath(dir, "swapped_c.bds")
        @test dskobj(swapped) == [401, 499]
        ctx = SpiceContext()
        furnsh(ctx, swapped)
        ref = SpiceContext()
        furnsh(ref, cfile)
        for _ in 1:500
            v = V3(randn(rng, 3)) * 8
            d = V3(randn(rng, 3)) * 2 - v
            @test dskxsi(ctx, false, 499, (), 0.0, "IAU_MARS", v, d) ==
                  dskxsi(ref, false, 499, (), 0.0, "IAU_MARS", v, d)
        end
    end

    @apart @testset "(f) the DLA segment list" begin
        h = dasopr(cfile)
        ch = c_dasopr(cfile)
        d = dlabfs(h)
        cd = SPICE.dlabfs(ch)
        n = 0
        fwd = DLADescr[]
        while d !== nothing
            n += 1
            push!(fwd, d)
            @test dla_tuple(d) == dla_tuple(cd)
            d = dlafns(h, d)
            cd = SPICE.dlafns(ch, cd)
        end
        @test cd === nothing && n == length(specs)
        d = dlabbs(h)
        cd = SPICE.dlabbs(ch)
        bwd = DLADescr[]
        while d !== nothing
            push!(bwd, d)
            @test dla_tuple(d) == dla_tuple(cd)
            d = dlafps(h, d)
            cd = SPICE.dlafps(ch, cd)
        end
        @test cd === nothing && bwd == reverse(fwd)
        @test fwd[1].bwdptr == -1 && fwd[end].fwdptr == -1
        c_dascls(ch)
        dascls(h)
        # A DLA file without segments.
        e = joinpath(dir, "nosegment.bds")
        dskcls(dskopn(e, "no segment", 0))
        h = dasopr(e)
        @test dlabfs(h) === nothing && dlabbs(h) === nothing
        @test dasrdi(h, 1, 3) == [1000000, -1, -1]
        dascls(h)
        ch = c_dasopr(e)
        @test SPICE.dlabfs(ch) === nothing
        c_dascls(ch)
        @test dskobj(e) == Int[]
    end

    @apart @testset "(g) input checks" begin
        missing_file = joinpath(dir, "does_not_exist.das")
        c_open_kind(p) = c_kind(() -> ccall((:dasopr_c, DSKLIB), Cvoid, (Cstring, Ref{Cint}), p, Ref{Cint}(0)))
        @test neo_kind(() -> dasopr(missing_file)) === :filenotfound === c_open_kind(missing_file)
        @test neo_kind(() -> dasopr(" ")) === :blankfilename
        text = joinpath(dir, "text.tk")
        write(text, "KPL/FK\n\\begindata\nA = 1\n" * " "^1100)
        @test neo_kind(() -> dasopr(text)) in (:filarchmismatch, :unknownfilarc)
        short = joinpath(dir, "short.das")
        write(short, read(cfile)[1:500])
        @test neo_kind(() -> dasopr(short)) === :filereadfailed
        # A damaged FTP validation string, as a text-mode transfer leaves it.
        damaged = joinpath(dir, "damaged.bds")
        bytes = read(cfile)
        bytes[710] = 0x20
        write(damaged, bytes)
        @test neo_kind(() -> dasopr(damaged)) === :ftpxfererror === c_open_kind(damaged)
        # A directory that points outside the file.
        cut = joinpath(dir, "cut.bds")
        write(cut, read(cfile)[1:end-4096])
        @test neo_kind(() -> dasopr(cut)) === :baddasdirectory
        bytes = read(cfile)
        bytes[85:92] = codeunits("VAX-GFLT")
        vax = joinpath(dir, "vax.bds")
        write(vax, bytes)
        @test neo_kind(() -> dasopr(vax)) === :unknownbff

        h = dasopr(cfile)
        lastc, lastd, lasti = daslla(h)
        @test neo_kind(() -> dasrdi(h, 0, 3)) === :dasnosuchaddress
        @test neo_kind(() -> dasrdi(h, 1, lasti + 1)) === :dasnosuchaddress
        @test neo_kind(() -> dasrdd(h, lastd, lastd + 1)) === :dasnosuchaddress
        @test neo_kind(() -> dasrdc(h, 1, 1)) === :dasnosuchaddress
        @test neo_kind(() -> dasadi(h, Cint[1])) === :dasinvalidaccess
        @test neo_kind(() -> dasadd(h, [1.0])) === :dasinvalidaccess
        @test neo_kind(() -> dasadc(h, "a")) === :dasinvalidaccess
        @test neo_kind(() -> dasudi(h, 1, 1, Cint[1])) === :dasinvalidaccess
        @test neo_kind(() -> dasac(h, ["a"])) === :dasinvalidaccess
        @test neo_kind(() -> dasdc(h)) === :dasinvalidaccess
        @test neo_kind(() -> dlabns(h)) === :dasinvalidaccess
        dascls(h)
        for f in (dashfn, dasrfr, daslla, dashfs, dasec, dlabfs, dlabbs)
            @test neo_kind(() -> f(h)) === :dasnosuchhandle
        end
        @test neo_kind(() -> dasrdi(h, 1, 1)) === :dasnosuchhandle
        @test occursin("closed", sprint(show, h))

        @test neo_kind(() -> dasonw(cfile, "TEST", "x", 0)) === :fileopenfailed
        @test isfile(cfile) && getfat(cfile) == ("DAS", "DSK")       # the file is not damaged
        @test neo_kind(() -> dasonw(" ", "TEST", "x", 0)) === :blankfilename
        @test neo_kind(() -> dasonw(joinpath(dir, "n1.das"), " ", "x", 0)) === :blankfiletype
        @test neo_kind(() -> dasonw(joinpath(dir, "n2.das"), "T\x01", "x", 0)) === :illegalcharacter
        @test neo_kind(() -> dasonw(joinpath(dir, "n3.das"), "TEST", "x", -1)) === :invalidcount
        @test neo_kind(() -> dlaopn(joinpath(dir, "n4.das"), "DSK", "x", -1)) === :badrecordcount
        @test !any(isfile, joinpath.(dir, ["n1.das", "n2.das", "n3.das", "n4.das"]))
        w = dasonw(joinpath(dir, "w.das"), "TEST", "x", 0)
        @test neo_kind(() -> dasac(w, String[])) === :invalidargument
        @test neo_kind(() -> dasac(w, ["tab\there"])) === :illegalcharacter
        @test neo_kind(() -> dasudd(w, 1, 2, [1.0, 2.0])) === :dasnosuchaddress
        @test_throws InexactError dasadi(w, [2^40])
        @test occursin("write", sprint(show, w))
        dascls(w)
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
