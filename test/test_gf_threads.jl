# Geometry finder, concurrent use. A search holds no global state: the step, the
# tolerance and the callbacks are arguments, and the search reads one immutable snapshot
# of the loaded kernels. This group checks three properties, each bit for bit:
#
# 1. the same search made on many tasks at the same time gives the window of the serial
#    search;
# 2. a search that uses tasks inside (`threads=true`) gives the window of the serial
#    search, for every relation;
# 3. searches on a context give the serial window while another task loads and unloads
#    a kernel that the search does not use.
#
# CSPICE is not used here, because it keeps one process-wide state.
#
# The default run has 35 searches (each of the five calls with each of the seven
# relations once), each made twice on many tasks and once on each of three tasks during
# the loads and unloads; the exhaustive run (`NEOSPICE_TEST_FULL=1`) has 70 searches,
# four times on many tasks and three times on each task but one of the process.
#
# Behaviour specification: test/features/gf_threads.feature

isdefined(@__MODULE__, :GFO) || include("gforacle.jl")

"One geometry finder search of the thread tests: a label and a function of a context."
struct GFThreadJob{F}
    label::String
    run::F
end

"""
    gf_thread_jobs(rng, n) -> Vector{GFThreadJob}

`n` searches with random confinement windows of 2000-2001, one for each entry of a
cycle through the stage 1 calls and the relations. Postcondition: the length is `n`.
"""
function gf_thread_jobs(rng::AbstractRNG, n::Integer)
    day = 86400.0
    rels = ("<", "=", ">", "LOCMIN", "LOCMAX", "ABSMIN", "ABSMAX")
    jobs = GFThreadJob[]
    # Invariant: `jobs` holds one job for each index before `i`.
    for i in 1:n
        cn = gf_random_window(rng, rand(rng, (1, 3, 8)), 0.0, 400 * day)
        rel = rels[mod1(i, length(rels))]
        adjust = rel in ("ABSMIN", "ABSMAX") && isodd(i ÷ 7)
        kind = mod1(i, 5)
        opt = (; threads=isodd(i ÷ 5))
        job = if kind == 1
            GFThreadJob("gfdist $rel", c -> gfdist(c, "MOON", "LT+S", "EARTH", rel, 3.8e5,
                adjust ? 2000.0 : 0.0, day, cn; opt...))
        elseif kind == 2
            GFThreadJob("gfsep $rel", c -> gfsep(c, "MOON", "SPHERE", "IAU_MOON", "SUN",
                "SPHERE", "IAU_SUN", "LT+S", "EARTH", rel, 1.2, adjust ? 0.05 : 0.0, day, cn;
                opt...))
        elseif kind == 3
            GFThreadJob("gfrr $rel", c -> gfrr(c, "MOON", "NONE", "EARTH", rel, 0.01,
                adjust ? 0.005 : 0.0, day, cn; opt...))
        elseif kind == 4
            GFThreadJob("gfpa $rel", c -> gfpa(c, "MOON", "SUN", "LT+S", "EARTH", rel, 1.0,
                adjust ? 0.05 : 0.0, day, cn; opt...))
        else
            GFThreadJob("gfposc $rel", c -> gfposc(c, "MOON", "IAU_EARTH", "LT+S", "EARTH",
                "LATITUDINAL", "LONGITUDE", rel, 3.0, adjust ? 0.02 : 0.0, 3600.0, cn;
                opt...))
        end
        push!(jobs, job)
    end
    @assert length(jobs) == n
    return jobs
end

"User function of the searches with `gfuds` and `gfudb`, and its derivative."
gft_f(t) = sin(t) + 0.3 * sin(2.7 * t + 0.4)
gft_d(t) = cos(t) + 0.81 * cos(2.7 * t + 0.4)

@gf_light group("gf threads: parallel searches, threads keyword, load and unload") do
    rng = MersenneTwister(SEED + 6400)
    ctx = neo_context(Kernels.full_set())
    day = 86400.0
    jobs = gf_thread_jobs(rng, nsample(35, 70))
    serial = [collect(j.run(ctx)) for j in jobs]

    @testset "the searches are not trivial" begin
        @test count(!isempty, serial) > length(serial) ÷ 2
        @test all(issorted, serial)
    end

    @testset "threads keyword: bit for bit, every relation" begin
        full = Window([0.0, 3 * 365 * day])
        many = gf_random_window(rng, 24, 0.0, 3 * 365 * day)
        for cn in (full, many), rel in ("<", "=", ">", "LOCMIN", "LOCMAX", "ABSMIN", "ABSMAX"),
            adjust in (0.0, 1500.0)
            adjust != 0 && !(rel in ("ABSMIN", "ABSMAX")) && continue
            a = gfdist(ctx, "MOON", "LT+S", "EARTH", rel, 3.8e5, adjust, day, cn)
            b = gfdist(ctx, "MOON", "LT+S", "EARTH", rel, 3.8e5, adjust, day, cn; threads=true)
            @test collect(a) == collect(b)
            a = gfposc(ctx, "MOON", "IAU_EARTH", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION",
                rel, 1.0, adjust * 1.0e-5, 3600.0, wnintd(cn, Window([0.0, 60 * day])))
            b = gfposc(ctx, "MOON", "IAU_EARTH", "NONE", "EARTH", "RA/DEC", "RIGHT ASCENSION",
                rel, 1.0, adjust * 1.0e-5, 3600.0, wnintd(cn, Window([0.0, 60 * day]));
                threads=true)
            @test collect(a) == collect(b)
        end
        f = GFScalar(gft_f)
        d = gf_decreasing(gft_d)
        above = GFState(t -> gft_f(t) > 0.2)
        cn = gf_random_window(rng, 16, 0.0, 200.0)
        for rel in ("<", "=", ">", "LOCMIN", "LOCMAX", "ABSMIN", "ABSMAX")
            @test collect(gfuds(f, d, rel, 0.2, 0.0, 0.1, cn)) ==
                  collect(gfuds(f, d, rel, 0.2, 0.0, 0.1, cn; threads=true))
        end
        @test collect(gfudb(above, 0.1, cn)) == collect(gfudb(above, 0.1, cn; threads=true))
    end

    if Threads.nthreads() < 2
        skip_notice("gf threads", "Threads.nthreads() == 1")
    else
        @testset "the same searches on many tasks" begin
            out = Vector{Vector{Float64}}(undef, nsample(2, 4) * length(jobs))
            Threads.@threads for k in eachindex(out)
                out[k] = collect(jobs[mod1(k, length(jobs))].run(ctx))
            end
            @test all(k -> out[k] == serial[mod1(k, length(jobs))], eachindex(out))
        end

        @testset "the tolerance of one context does not reach another" begin
            other = neo_context(Kernels.full_set())
            gfstol(other, 1.0e-2)
            cn = Window([0.0, 200 * day])
            ref = collect(gfdist(ctx, "MOON", "NONE", "EARTH", "=", 3.8e5, 0.0, day, cn))
            coarse = collect(gfdist(other, "MOON", "NONE", "EARTH", "=", 3.8e5, 0.0, day, cn))
            res = Vector{Bool}(undef, 64)
            Threads.@threads for k in eachindex(res)
                c, want = isodd(k) ? (ctx, ref) : (other, coarse)
                res[k] = collect(gfdist(c, "MOON", "NONE", "EARTH", "=", 3.8e5, 0.0, day, cn)) == want
            end
            @test all(res)
            @test ref != coarse && length(ref) == length(coarse)
            @test maximum(abs.(ref .- coarse)) <= 2.0e-2
            @test gfstol(ctx) == 1.0e-6 && gfstol(other) == 1.0e-2
        end

        @testset "searches while another task loads and unloads a kernel" begin
            extra = Kernels.GM_TPC
            shared = neo_context(Kernels.full_set())
            stop = Threads.Atomic{Bool}(false)
            bad = Threads.Atomic{Int}(0)
            total = Threads.Atomic{Int}(0)
            writer = Threads.@spawn begin
                cycles = 0
                while !stop[]
                    furnsh(shared, extra)
                    yield()
                    unload(shared, extra)
                    cycles += 1
                end
                cycles
            end
            nrep = nsample(1, 3)
            nreaders = FULL ? max(2, Threads.nthreads() - 1) : clamp(Threads.nthreads() - 1, 2, 3)
            readers = [Threads.@spawn begin
                for rep in 1:nrep, k in eachindex(jobs)
                    ok = try
                        collect(jobs[k].run(shared)) == serial[k]
                    catch
                        false
                    end
                    ok || Threads.atomic_add!(bad, 1)
                    Threads.atomic_add!(total, 1)
                end
            end for _ in 1:nreaders]
            foreach(wait, readers)
            stop[] = true
            cycles = fetch(writer)
            measured("  gf threads: $(total[]) searches on $(length(readers)) tasks during " *
                     "$cycles load/unload cycles; $(bad[]) differ from the serial search")
            @test bad[] == 0
            @test cycles > 0
            @test total[] == nrep * length(jobs) * length(readers)
            @test [collect(j.run(shared)) for j in jobs] == serial
        end
    end
end
