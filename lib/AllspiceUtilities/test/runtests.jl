# Tests of AllspiceUtilities. The oracle is the NAIF program of the same name from
# CSPICE_jll; its text report is parsed and compared with the structured result.

using Test
using Allspice: Allspice, SpiceError, Window, spkezr, furnsh, SpiceContext
using AllspiceUtilities
import CSPICE_jll

include(joinpath(@__DIR__, "..", "..", "..", "test", "kernels.jl"))
Kernels.fetch_kernels()
const K = Kernels
const NEED = (K.LSK, K.DE440S, K.DE405, K.MOON_BPC, K.MOON_TF)

"Run a NAIF program in `dir` and return its standard output."
naif(prog, args...; dir=pwd()) = cd(() -> read(ignorestatus(`$(prog()) $args`), String), dir)

"The lines of the comment area as `commnt -r` prints them, without trailing blanks."
naifcomments(f) = rstrip.(split(strip(naif(CSPICE_jll.commnt, "-r", f)), '\n'))

"The floating-point numbers of a line of NAIF output."
numbers(line) = [parse(Float64, m.match) for m in
    eachmatch(r"[-+]?\d+\.\d+E[-+]\d+|[-+]?\d+\.\d+(?![\d.E])", line)]

"Rows of numbers of a NAIF dump report (the lines that start with a number)."
dumprows(text) = [numbers(l) for l in split(text, '\n') if occursin(r"^[-+]?\d", l)]

"The numbers of the first line of `text` that contains `key`."
after(text, key) = numbers(split(text, '\n')[findfirst(l -> occursin(key, l), split(text, '\n'))])

if !isempty(K.missing_kernels(NEED...))
    K.REQUIRED && error("kernels missing: $(K.missing_kernels(NEED...))")
    @warn "AllspiceUtilities tests skipped: kernels missing" K.missing_kernels(NEED...)
else

@testset "AllspiceUtilities" begin
    tmp = mktempdir()

    @testset "brief" begin
        cov = brief(K.DE440S, K.MOON_BPC, K.LSK)
        # The NAIF program takes SPK and PCK files in separate runs only.
        spkout = naif(CSPICE_jll.brief, K.DE440S)
        out = spkout * naif(CSPICE_jll.brief, K.MOON_BPC)
        ids = [parse(Int, m[1]) for m in eachmatch(r"\((\d+)\)", spkout)]
        @test [c.id for c in cov if c.kind === :spk] == sort(ids)
        @test [c.id for c in cov if c.kind === :pck] == [31008]
        printed = sprint(printbrief, cov)
        naiflines = [strip(l) for l in split(out, '\n') if occursin(r"^\s+\d{4} [A-Z]{3}", l)]
        @test [strip(l) for l in split(printed, '\n') if occursin(r"^\s+\d{4} [A-Z]{3}", l)] ==
            naiflines
        @test_throws ArgumentError brief(joinpath(@__DIR__, "runtests.jl"))
    end

    @testset "commnt" begin
        f = joinpath(tmp, "c.bsp")
        AllspiceUtilities.spkmerge(f, [MergeSource(K.DE440S; bodies=[301], stop=0.0)])
        @test isempty(comments(f))
        addcomments(f, ["first line", "", "  third line"])
        @test naifcomments(f) == ["first line", "", "  third line"]
        txt = joinpath(tmp, "c.txt")
        write(txt, "added by NAIF\nsecond\n")
        naif(CSPICE_jll.commnt, "-a", f, txt)
        @test comments(f)[end-1:end] == ["added by NAIF", "second"]
        deletecomments(f)
        @test occursin("no comments", naif(CSPICE_jll.commnt, "-r", f))
        # A DAS file (DSK) goes through the DAS calls.
        g = joinpath(tmp, "c.bds")
        Allspice.dascls(Allspice.dasonw(g, "DSK", "c.bds", 0))
        addcomments(g, ["dsk comment"])
        @test naifcomments(g) == ["dsk comment"]
        deletecomments(g)
        @test isempty(comments(g))
        @test_throws ArgumentError comments(K.LSK)
    end

    @testset "spkmerge against the NAIF program" begin
        # Earlier sources have priority; the later one fills the gaps, per body.
        y(yr) = (yr - 2000) * 365.25 * 86400.0
        spec = [MergeSource(K.DE440S; bodies=[301, 3], start=y(2000), stop=y(2015)),
            MergeSource(K.DE405; bodies=[301, 399], start=y(2010), stop=y(2020))]
        ours = joinpath(tmp, "ours.bsp")
        segs = spkmerge(ours, spec; comments=["merged by AllspiceUtilities"])
        setup = joinpath(tmp, "setup.txt")
        write(setup, """
            LEAPSECONDS_KERNEL = $(K.LSK)
            SPK_KERNEL = naif.bsp
            """ * join(["""
              SOURCE_SPK_KERNEL = $(s.path)
                BODIES = $(join(s.bodies, ", "))
                BEGIN_TIME = $(Allspice.etcal(s.start)) TDB
                END_TIME = $(Allspice.etcal(s.stop)) TDB
            """ for s in spec]))
        naif(CSPICE_jll.spkmerge, setup; dir=tmp)
        h = Allspice.dafopr(joinpath(tmp, "naif.bsp"))
        theirs = []
        Allspice.dafbfs(h)
        while Allspice.daffna(h)
            body, center, frame, type, a, b = Allspice.spkuds(Allspice.dafgs(h))
            push!(theirs, (; body, center, frame, type, start=a, stop=b,
                ident=Allspice.dafgn(h)))
        end
        Allspice.dafcls(h)
        @test [Base.structdiff(s, NamedTuple{(:source,)}) for s in segs] == theirs
        @test comments(ours)[1] == "merged by AllspiceUtilities"
        c1, c2 = SpiceContext(), SpiceContext()
        furnsh(c1, ours)
        furnsh(c2, joinpath(tmp, "naif.bsp"))
        for t in range(y(2000), y(2020); length=41)
            @test spkezr(c1, 301, t, "J2000", "NONE", 3) == spkezr(c2, 301, t, "J2000", "NONE", 3)
        end
        @test_throws SpiceError spkmerge(joinpath(tmp, "none.bsp"),
            [MergeSource(K.DE440S; bodies=[-999])])
        @test !isfile(joinpath(tmp, "none.bsp"))
    end

    @testset "spkdiff against the NAIF program" begin
        args = ("-k", K.LSK, "-b1", "301", "-c1", "3", "-r1", "J2000", "-n", "7",
            "-b", "1990-01-01", "-e", "2030-01-01")
        d = spkdiff(K.DE440S, K.DE405; body="MOON", center=3, kernels=K.LSK, n=7,
            start="1990-01-01", stop="2030-01-01")
        dump = dumprows(naif(CSPICE_jll.spkdiff, args..., "-t", "dump", K.DE440S, K.DE405))
        @test [r[1] for r in dump] ≈ d.times rtol = 1e-13
        @test all(isapprox(r[2:7], a - b; rtol=1e-12) for (r, a, b) in zip(dump, d.state1, d.state2))
        basic = naif(CSPICE_jll.spkdiff, args..., "-t", "basic", K.DE440S, K.DE405)
        s = basicstats(d)
        @test [s.maxrelpos, s.avgrelpos, s.maxrelvel, s.avgrelvel] ≈
            [after(basic, "Position:"); after(basic, "Velocity:")] rtol = 1e-12
        @test [s.maxpos, s.avgpos, s.maxvel, s.avgvel] ≈
            [after(basic, "Position (km)"); after(basic, "Velocity (km/s)")] rtol = 1e-12
        stats = naif(CSPICE_jll.spkdiff, args..., "-t", "stats", K.DE440S, K.DE405)
        v = viewstats(d)
        for (key, tup) in (("1a)", v.mean), ("2a)", v.meanabs), ("3a)", v.rms),
                ("4a)", v.atmaxrel), ("5a)", v.atmaxabs))
            n = key[1]
            ref = [after(stats, "$(n)$(c))")[1] for c in "abcd"]
            @test collect(tup) ≈ ref rtol = 1e-9
        end
        @test v.epochmaxabs ≈ after(stats, "5e)")[1] rtol = 1e-13
        @test_throws ArgumentError spkdiff(K.DE440S, K.DE405; body=301, center=3,
            kernels=K.LSK, start=1e12, stop=2e12)
    end

    @testset "frmdiff against the NAIF program" begin
        ks = [K.LSK, K.MOON_BPC, K.MOON_TF]
        args = ("-k", join(ks, " "), "-f1", "MOON_PA_DE440", "-t1", "J2000",
            "-f2", "MOON_ME_DE440_ME421", "-t2", "J2000", "-s", "5000000",
            "-b", "2000-01-01", "-e", "2001-01-01", "-a", "yes")
        for m in (:from, :to)
            d = frmdiff(; kernels=ks, from="MOON_PA_DE440", to="J2000",
                from2="MOON_ME_DE440_ME421", start="2000-01-01", stop="2001-01-01",
                step=5e6, angvel=true, avframe=m)
            rows = dumprows(naif(CSPICE_jll.frmdiff, args..., "-m", string(m), "-t", "dumpaa"))
            # The NAIF program leaves out the last epoch before the stop.
            @test length(rows) == length(d.times) - 1
            idx = Int[]
            for r in rows
                k = findfirst(t -> isapprox(t, r[1]; rtol=1e-13), d.times)
                @test k !== nothing
                push!(idx, k)
                @test r[2] ≈ d.angle[k] rtol = 1e-11
                @test r[3:5] ≈ d.axis[k] rtol = 1e-9
                @test r[6:8] ≈ d.avdiff[k] rtol = 1e-7 atol = 1e-20
            end
            basic = naif(CSPICE_jll.frmdiff, args..., "-m", string(m), "-t", "basic")
            s = basicstats(d)
            @test s.maxangle ≈ after(basic, "Rotation (rad)")[1] rtol = 1e-11
            @test maximum(k -> AllspiceUtilities.norm(d.avdiff[k]), idx) ≈
                after(basic, "Angular Velocity")[1] rtol = 1e-9
        end
        # Without start and stop the window is the coverage of the frame in the PCK.
        d = frmdiff(K.MOON_BPC; kernels=[K.LSK, K.MOON_TF], from="MOON_PA_DE440",
            to="J2000", from2="MOON_ME_DE440_ME421", n=3)
        @test (first(d.times), last(d.times)) == Allspice.wnfetd(Allspice.pckcov(K.MOON_BPC, 31008), 1)
        @test_throws ArgumentError frmdiff(; kernels=ks, from="MOON_PA_DE440", to="J2000")
    end

    @testset "epochs of a window with a gap" begin
        w = Window([0.0, 10.0, 20.0, 30.0])
        @test AllspiceUtilities._sample(w, 5, nothing) == [0.0, 5.0, 10.0, 25.0, 30.0]
        @test AllspiceUtilities._sample(w, 0, 4.0) == [0.0, 4.0, 8.0, 10.0, 20.0, 24.0, 28.0, 30.0]
    end
end

end
