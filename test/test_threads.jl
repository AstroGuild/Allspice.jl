# Group 7: thread safety by immutable snapshots.
#
# Feature file: features/concurrent_use.feature
#
# CSPICE is not called in this group: the reference answers come from Allspice itself,
# computed serially before any task starts.
#
# The writer task cycles one context through four load lists:
#
#   base                         states A, rotations A
#   base + de405                 states B, rotations A
#   base + de405 + pck00010      states B, rotations B
#   base + de405                 states B, rotations A
#
# `de405.bsp` is loaded after `de440s.bsp` and takes precedence, so a state query has
# exactly two valid answers; `pck00010.tpc` replaces the Neptune orientation constants of
# `pck00011.tpc`, so an `IAU_NEPTUNE` rotation has exactly two valid answers. A query that
# mixes two snapshots (segment list of one, data of the other) would give neither.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

group("7. thread safety", Kernels.FULL_SET_FILES..., Kernels.DE405, Kernels.PCK10) do
    if Threads.nthreads() < 2
        skip_notice("7. thread safety", "Threads.nthreads() == 1; run with " *
            "`julia --threads=8 --project=test test/runtests.jl`")
        return
    end
    rng = MersenneTwister(SEED + 7)
    base = Kernels.full_set()
    epochs = random_epochs(rng, 64, DE440S_COVERAGE...)
    pairs = [(301, 399), (399, 10), (10, 301), (5, 399), (299, 199), (4, 0)]
    nq = length(epochs) * length(pairs)
    qstate(ctx, i) = begin
        et = epochs[(i-1)%length(epochs)+1]
        targ, obs = pairs[(i-1)÷length(epochs)+1]
        spkezr(ctx, targ, et, "J2000", "NONE", obs)
    end
    qrot(ctx, i) = pxform(ctx, "J2000", "IAU_NEPTUNE", epochs[(i-1)%length(epochs)+1])
    qder(ctx, i) = begin
        et = epochs[(i-1)%length(epochs)+1]
        targ, obs = pairs[(i-1)÷length(epochs)+1]
        spkder(ctx, targ, et, "IAU_EARTH", obs, Val(2))
    end

    # Reference answers of the two snapshots, computed serially.
    ctxA = neo_context(base)
    ctxB = neo_context([base..., Kernels.DE405, Kernels.PCK10])
    stateA = [qstate(ctxA, i) for i in 1:nq]
    stateB = [qstate(ctxB, i) for i in 1:nq]
    rotA = [qrot(ctxA, i) for i in 1:nq]
    rotB = [qrot(ctxB, i) for i in 1:nq]
    @test all(stateA[i] != stateB[i] for i in 1:nq)
    @test all(rotA[i] != rotB[i] for i in 1:nq)

    @apart @testset "pure-query parallel loop equals the serial results bit for bit" begin
        serial_state = stateA
        serial_rot = rotA
        serial_der = [qder(ctxA, i) for i in 1:nq]
        par_state = similar(serial_state); par_rot = similar(serial_rot)
        par_der = similar(serial_der)
        reps = 20
        ok = Threads.Atomic{Bool}(true)
        for _ in 1:reps
            Threads.@threads for i in 1:nq
                par_state[i] = qstate(ctxA, i)
                par_rot[i] = qrot(ctxA, i)
                par_der[i] = qder(ctxA, i)
            end
            (par_state == serial_state && par_rot == serial_rot && par_der == serial_der) ||
                (ok[] = false)
        end
        @test ok[]
        println("  ", reps * nq * 3, " parallel queries on ", Threads.nthreads(), " threads")
    end

    @apart @testset "queries during furnsh and unload see one valid snapshot" begin
        ctx = neo_context(base)
        stop = Threads.Atomic{Bool}(false)
        cycles = Threads.Atomic{Int}(0)
        writer = Threads.@spawn begin
            while !stop[]
                furnsh(ctx, Kernels.DE405)
                furnsh(ctx, Kernels.PCK10)
                unload(ctx, Kernels.PCK10)
                unload(ctx, Kernels.DE405)
                Threads.atomic_add!(cycles, 1)
                yield()
            end
        end
        nreaders = max(2, Threads.nthreads() - 1)
        target_cycles = 30                    # writer cycles that the readers must span
        deadline = time() + 120
        readers = map(1:nreaders) do r
            Threads.@spawn begin
                bad = 0; sawA = 0; sawB = 0; rA = 0; rB = 0; round = 0
                while (cycles[] < target_cycles && time() < deadline) || round < 20
                  round += 1
                  for i in 1:nq
                    # Spread the readers over the queries.
                    j = (i + 37 * r + round - 1) % nq + 1
                    s = qstate(ctx, j)
                    if s == stateA[j]
                        sawA += 1
                    elseif s == stateB[j]
                        sawB += 1
                    else
                        bad += 1
                    end
                    R = qrot(ctx, j)
                    if R == rotA[j]
                        rA += 1
                    elseif R == rotB[j]
                        rB += 1
                    else
                        bad += 1
                    end
                  end
                end
                (; bad, sawA, sawB, rA, rB, n=round * nq)
            end
        end
        results = fetch.(readers)                 # rethrows an exception of a reader
        stop[] = true
        fetch(writer)
        total = (; bad=sum(r.bad for r in results), sawA=sum(r.sawA for r in results),
            sawB=sum(r.sawB for r in results), rA=sum(r.rA for r in results),
            rB=sum(r.rB for r in results))
        measured("  $(nreaders) readers, $(cycles[]) writer cycles, " *
            "$(2 * sum(r.n for r in results)) queries: states A/B = " *
            "$(total.sawA)/$(total.sawB), rotations A/B = $(total.rA)/$(total.rB), " *
            "invalid = $(total.bad)")
        @test total.bad == 0
        @test total.sawA + total.sawB == sum(r.n for r in results)
        @test total.rA + total.rB == sum(r.n for r in results)
        @test cycles[] >= target_cycles
        # Both snapshots were observed; if not, the run did not exercise the race.
        @test total.sawA > 0 && total.sawB > 0
        # The writer leaves the context as it was.
        @test ktotal(ctx) == length(base)
        @test [qstate(ctx, i) for i in 1:nq] == stateA
        @test [qrot(ctx, i) for i in 1:nq] == rotA
    end

    @apart @testset "independent contexts in parallel" begin
        answers = Vector{Any}(undef, 16)
        Threads.@threads for k in 1:16
            c = SpiceContext()
            furnsh(c, isodd(k) ? Kernels.DE440S : Kernels.DE405)
            answers[k] = [spkezr(c, 301, et, "J2000", "NONE", 10) for et in epochs]
            kclear(c)
        end
        only440 = neo_context([Kernels.DE440S]); only405 = neo_context([Kernels.DE405])
        a = [spkezr(only440, 301, et, "J2000", "NONE", 10) for et in epochs]
        b = [spkezr(only405, 301, et, "J2000", "NONE", 10) for et in epochs]
        @test all(answers[k] == (isodd(k) ? a : b) for k in 1:16)
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
