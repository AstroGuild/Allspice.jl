# Group EK 3: EK queries on many threads while another task loads and unloads an EK file.
#
# CSPICE is called only to write the EK files, before any task starts; it has one global
# query result and cannot run this test. The reference answers come from Allspice,
# computed serially (group `ek_query` compares the same kinds of queries with CSPICE).
#
# The writer task cycles one context through two load lists: with and without the file
# `groups.bes`, which holds one of the two segments of the table GROUPS. A query that
# reads GROUPS has exactly two valid answers; a query that mixes two snapshots (the
# table catalogue of one, the files of the other) would give neither or would fail.

# The code of the EK groups is long and most of it runs once, so Julia compiles the
# module of the tests without optimization while this file runs (5 s less for the five
# EK groups); the last line of the file gives the default level back. The level of
# Allspice does not change.
Base.Experimental.@optlevel 0

isdefined(@__MODULE__, :EKLIB) || include("ekfixtures.jl")

group("EK queries: thread safety", Kernels.LSK) do
    if Threads.nthreads() < 2
        skip_notice("EK queries: thread safety", "Threads.nthreads() == 1; run with " *
            "`julia --threads=8 --project=test test/runtests.jl`")
        return
    end
    _, paths = ek_standard_files(SEED + 1; scale=0.15)      # the files of the group `ek_query`
    base = [Kernels.LSK, paths[1], paths[2]]
    extra = paths[3]
    queries = ["SELECT GID, GNAME, WEIGHT FROM GROUPS WHERE GID < 120 ORDER BY GNAME, GID",
        "SELECT GID FROM GROUPS WHERE GNAME LIKE '*a*' OR WEIGHT > 0.5 OR K IS NULL",
        "SELECT E.ID, G.GID, G.NOTE FROM EVENTS E, GROUPS G WHERE E.GRP = G.GID ORDER BY E.ID DESC, G.GID",
        "SELECT E.ID, L.TAG, G.T0 FROM EVENTS E, LINKS L, GROUPS G WHERE E.ID = L.EID AND L.GID <= G.GID AND G.GID < 40",
        "SELECT A.GID, B.GID FROM GROUPS A, GROUPS B WHERE A.GID = B.K AND A.GID BETWEEN 10 AND 90",
        "SELECT GID, T0 FROM GROUPS WHERE T0 > '2000 JAN 1' ORDER BY T0 DESC",
        "SELECT ID, NAME, IVEC, TAGS FROM EVENTS WHERE ID BETWEEN 20 AND 60 ORDER BY NAME DESC, ID",
        "SELECT ID, GNAME FROM EVENTS, GROUPS WHERE VAL < WEIGHT AND ID < 12 AND GID < 30"]
    answer(ctx, q) = begin
        res, err, msg = ekfind(ctx, q)
        err && error("query failed: $msg")
        (ntab=ekntab(ctx), rows=neo_ek_result(res))
    end
    ctxA = neo_context(base)
    ctxB = neo_context([base..., extra])
    ansA = [answer(ctxA, q) for q in queries]
    ansB = [answer(ctxB, q) for q in queries]
    differ = [ansA[i].rows != ansB[i].rows for i in eachindex(queries)]
    @test count(differ) >= 6                  # the queries on GROUPS see the extra segment
    @test !differ[7]                          # the query on EVENTS alone does not
    nq = length(queries)

    ek_part("parallel queries on one snapshot equal the serial results") do
        reps = 12
        par = Vector{Any}(undef, nq * reps)
        Threads.@threads for k in 1:nq*reps
            par[k] = answer(ctxB, queries[(k-1)%nq+1])
        end
        @test all(par[k] == ansB[(k-1)%nq+1] for k in 1:nq*reps)
        # One result object read by many tasks at the same time.
        res, _, _ = ekfind(ctxB, queries[3])
        serial = neo_ek_result(res)
        same = Vector{Bool}(undef, 64)
        Threads.@threads for k in 1:64
            same[k] = neo_ek_result(res) == serial
        end
        @test all(same)
    end

    ek_part("queries during furnsh and unload of an EK see one valid snapshot") do
        ctx = neo_context(base)
        stop = Threads.Atomic{Bool}(false)
        cycles = Threads.Atomic{Int}(0)
        writer = Threads.@spawn begin
            while !stop[]
                furnsh(ctx, extra)
                yield()
                unload(ctx, extra)
                Threads.atomic_add!(cycles, 1)
                yield()
            end
        end
        nreaders = 8
        target_cycles = 40
        deadline = time() + 120
        readers = map(1:nreaders) do r
            Threads.@spawn begin
                bad = 0; sawA = 0; sawB = 0; n = 0; round = 0; held = 0
                while (cycles[] < target_cycles && time() < deadline) || round < 6
                    round += 1
                    for i in 1:nq
                        j = (i + 3 * r + round - 1) % nq + 1
                        res, err, _ = ekfind(ctx, queries[j])
                        # The result is read after the writer had time to unload its files.
                        yield()
                        rows = err ? nothing : neo_ek_result(res)
                        n += 1
                        if rows == ansA[j].rows
                            sawA += !differ[j] ? 0 : 1
                        elseif rows == ansB[j].rows
                            sawB += 1
                        else
                            bad += 1
                        end
                        held += res.nmrows
                    end
                end
                (; bad, sawA, sawB, n, held)
            end
        end
        results = fetch.(readers)
        stop[] = true
        fetch(writer)
        nbad = sum(r.bad for r in results)
        nA = sum(r.sawA for r in results)
        nB = sum(r.sawB for r in results)
        measured("  EK: $(nreaders) readers on $(Threads.nthreads()) threads, $(cycles[]) load/unload cycles, " *
                 "$(sum(r.n for r in results)) queries, $(sum(r.held for r in results)) rows read: " *
                 "without/with the extra segment = $nA/$nB, invalid = $nbad")
        @test nbad == 0
        @test cycles[] >= target_cycles
        @test nA > 0 && nB > 0
        @test ktotal(ctx, "EK") == 2
        @test [answer(ctx, q) for q in queries] == ansA
    end
end

Base.Experimental.@optlevel -1
