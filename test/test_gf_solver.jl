# Geometry finder, solver core: `gfuds`, `gfudb`, `gfrefn`, `uddf`, `uddc`, `gfstol` and
# the callbacks. No kernel is needed. The user functions are the same Julia
# functions for Allspice and for CSPICE, so the two searches evaluate the same numbers
# and every end point must be equal bit for bit: this group tests the step-and-bracket
# logic, the bisection, the confinement end points, the adjust value and the absolute
# extrema without any noise from the ephemeris.
#
# Behaviour specification: test/features/gf_solver.feature

isdefined(@__MODULE__, :GFO) || include("gforacle.jl")

# Scalar functions with their analytical derivatives.
gfv_f1(t) = sin(t)
gfv_d1(t) = cos(t)
gfv_f2(t) = sin(t) + 0.3 * sin(2.7 * t + 0.4)
gfv_d2(t) = cos(t) + 0.81 * cos(2.7 * t + 0.4)
gfv_f3(t) = 1.0e-3 * (t - 3.0) * (t - 11.0) * (t - 17.5) + 0.05 * t
gfv_d3(t) = 1.0e-3 * ((t - 11.0) * (t - 17.5) + (t - 3.0) * (t - 17.5) + (t - 3.0) * (t - 11.0)) + 0.05

@gf_light group("gf solver: user-defined functions, tolerance, callbacks") do
    rng = MersenneTwister(SEED + 6100)
    relations = ("<", "=", ">", "LOCMIN", "LOCMAX", "ABSMIN", "ABSMAX")

    # The scalar functions behind one type (`GFScalar`), with their analytical derivatives.
    f1, f2, f3 = GFScalar(gfv_f1), GFScalar(gfv_f2), GFScalar(gfv_f3)
    d1, d2, d3 = gfv_d1, gfv_d2, gfv_d3
    funcs = ((f1, d1, "sin"), (f2, d2, "two sines"), (f3, d3, "cubic"))

    @testset "gfrefn, uddf, uddc" begin
        for _ in 1:200
            a, b = 100 * randn(rng), 100 * randn(rng)
            @test gfrefn(a, b, true, false) === SPICE.gfrefn(a, b, true, false)
            lo, hi = minmax(a, b)
            @test lo <= gfrefn(a, b, false, true) <= hi
            x, dx = 10 * randn(rng), 10.0^rand(rng, -6:1)
            @test uddf(f2, x, dx) === GFO.uddf(f2, x, dx)
            @test uddc(f2, x, dx) === (GFO.uddf(f2, x, dx) < 0)
        end
        @test throws_kind(() -> uddf(f1, 1.0, 0.0), :dividebyzero)
        @test GFO.kind(() -> GFO.uddf(f1, 1.0, 0.0)) === :dividebyzero
    end

    @testset "gfuds: every relation, bit for bit" begin
        n = 0
        for (f, d, name) in funcs, rel in relations, trial in 1:6
            nint = rand(rng, (1, 1, 2, 4))
            cn = gf_random_window(rng, nint, 0.0, 25.0)
            step = rand(rng, (0.05, 0.2, 0.37))
            vals = [f(t) for t in range(0.0, 25.0; length=400)]
            refval = rel in ("<", "=", ">") ? rand(rng) * (maximum(vals) - minimum(vals)) * 0.8 +
                                              minimum(vals) : 0.0
            adjust = (rel in ("ABSMIN", "ABSMAX") && isodd(trial)) ?
                     0.2 * rand(rng) * (maximum(vals) - minimum(vals)) : 0.0
            qdec = gf_decreasing(d)
            w = gfuds(f, qdec, rel, refval, adjust, step, cn)
            ref = GFO.gfuds(f, qdec, rel, refval, adjust, step, collect(cn))
            @test collect(w) == ref
            gf_compare("uds " * gf_family(rel, adjust), w, ref; atol=0.0)
            @test all(t -> wnelmd(cn, t), w)
            n += 1
        end
        @test n == 3 * 7 * 6
    end

    @testset "gfuds: numerical test \"is decreasing\" (uddc)" begin
        for (f, _, name) in funcs, rel in relations
            cn = Window([1.0, 22.0])
            qdec = GFState(t -> uddc(f, t, 1.0e-3))
            refval = rel in ("<", "=", ">") ? 0.2 : 0.0
            w = gfuds(f, qdec, rel, refval, 0.0, 0.25, cn)
            @test collect(w) == GFO.gfuds(f, qdec, rel, refval, 0.0, 0.25, collect(cn))
        end
    end

    @testset "gfuds: events on the end points of the confinement window" begin
        # sin has roots at 0 and at pi (the double nearest to pi gives a small positive
        # sine) and extrema at pi/2 and 3 pi/2.
        qdec = gf_decreasing(d1)
        for cn in (Window([0.0, pi]), Window([0.0, 2pi]), Window([pi / 2, 3pi / 2]),
                Window([pi / 2, pi / 2]), Window([0.0, pi / 2, pi, 3pi / 2]),
                Window([pi, 2pi, 2pi, 3pi])), rel in relations
            refs = rel in ("<", "=", ">") ? (0.0, 1.0, -1.0, 0.5) : (0.0,)
            for refval in refs
                w = gfuds(f1, qdec, rel, refval, 0.0, 0.3, cn)
                @test collect(w) == GFO.gfuds(f1, qdec, rel, refval, 0.0, 0.3, collect(cn))
            end
        end
        # The interior root is always found. A root at an end of the window is found
        # only if the bisection, which works on the window expanded by 0.5 s, puts it
        # inside the window; CSPICE behaves in the same way (compared above).
        w = gfuds(f1, qdec, "=", 0.0, 0.0, 0.3, Window([0.0, 2pi]))
        @test any(t -> abs(t - pi) < 1e-6, w)
        @test all(t -> min(abs(t), abs(t - pi), abs(t - 2pi)) < 1e-6, w)
    end

    @testset "gfuds: searches with no result" begin
        qdec = gf_decreasing(d1)
        cn = Window([0.0, 20.0])
        for (rel, refval) in (("=", 2.0), ("<", -2.0), (">", 2.0), ("=", -1.5))
            w = gfuds(f1, qdec, rel, refval, 0.0, 0.3, cn)
            @test isempty(w)
            @test isempty(GFO.gfuds(f1, qdec, rel, refval, 0.0, 0.3, collect(cn)))
        end
        # A monotone function has no local extremum inside the window.
        for rel in ("LOCMIN", "LOCMAX")
            w = gfuds(GFScalar(t -> 2.0 * t), GFState(t -> false), rel, 0.0, 0.0, 0.5, cn)
            @test isempty(w)
            @test isempty(GFO.gfuds(t -> 2.0 * t, t -> false, rel, 0.0, 0.0, 0.5, collect(cn)))
        end
        @test isempty(gfuds(f1, qdec, "=", 0.0, 0.0, 0.3, Window()))
        @test wnsumd(gfuds(f1, qdec, "<", 2.0, 0.0, 0.3, cn))[1] == 20.0
    end

    @testset "gfudb: bit for bit" begin
        for trial in 1:40
            cn = gf_random_window(rng, rand(rng, (1, 2, 5)), -10.0, 30.0)
            step = rand(rng, (0.05, 0.3, 0.9))
            c = 0.6 * (2 * rand(rng) - 1)
            b = GFState(t -> f2(t) > c)
            w = gfudb(b, step, cn)
            ref = GFO.gfudb(b, step, collect(cn))
            @test collect(w) == ref
            gf_compare("udb state", w, ref; atol=0.0)
        end
        @test collect(gfudb(GFState(t -> true), 1.0, Window([0.0, 3.0, 5.0, 5.0]))) == [0.0, 3.0, 5.0, 5.0]
        @test isempty(gfudb(GFState(t -> false), 1.0, Window([0.0, 3.0])))
        @test isempty(gfudb(GFState(t -> true), 1.0, Window()))
    end

    @testset "gfstol: default of a context, keyword, CSPICE" begin
        ctx = SpiceContext()
        other = SpiceContext()
        @test gfstol(ctx) == 1.0e-6
        qdec = gf_decreasing(d2)
        cn = Window([0.0, 25.0])
        base = gfuds(ctx, f2, qdec, "=", 0.1, 0.0, 0.2, cn)
        try
            for tol in (1.0e-2, 1.0e-4, 1.0e-9)
                GFO.gfstol(tol)
                @test gfstol(ctx, tol) == tol
                @test gfstol(ctx) == tol
                @test gfstol(other) == 1.0e-6          # the default belongs to one context
                ref = GFO.gfuds(f2, qdec, "=", 0.1, 0.0, 0.2, collect(cn))
                @test collect(gfuds(ctx, f2, qdec, "=", 0.1, 0.0, 0.2, cn)) == ref
                # the keyword overrides the default, on every form
                @test collect(gfuds(other, f2, qdec, "=", 0.1, 0.0, 0.2, cn; tol=tol)) == ref
                @test collect(gfuds(Allspice.kernelset(other), f2, qdec, "=", 0.1, 0.0, 0.2, cn;
                    tol=tol)) == ref
                @test gfuds(other, f2, qdec, "=", 0.1, 0.0, 0.2, cn) == base
                refb = GFO.gfudb(t -> f2(t) > 0.1, 0.2, collect(cn))
                @test collect(gfudb(ctx, GFState(t -> f2(t) > 0.1), 0.2, cn)) == refb
            end
            kclear(ctx)
            @test gfstol(ctx) == 1.0e-9                # `kclear` keeps it, as in CSPICE
            @test throws_kind(() -> gfstol(ctx, 0.0), :invalidtolerance)
            @test throws_kind(() -> gfstol(ctx, -1.0), :invalidtolerance)
            @test GFO.kind(() -> GFO.gfstol(0.0)) === :invalidtolerance
            @test gfstol(ctx) == 1.0e-9
            # a coarse tolerance moves the root by less than the tolerance
            coarse = gfuds(ctx, f2, qdec, "=", 0.1, 0.0, 0.2, cn; tol=1.0e-2)
            @test length(coarse) == length(base)
            @test maximum(abs.(collect(coarse) .- collect(base))) <= 1.0e-2
        finally
            GFO.gfstol(1.0e-6)
        end
        # the default context
        old = gfstol()
        try
            @test gfstol(1.0e-3) == 1.0e-3
            GFO.gfstol(1.0e-3)
            @test collect(gfuds(f2, qdec, "=", 0.1, 0.0, 0.2, cn)) ==
                  GFO.gfuds(f2, qdec, "=", 0.1, 0.0, 0.2, collect(cn))
        finally
            gfstol(old)
            GFO.gfstol(1.0e-6)
        end
    end

    @testset "errors, with the kinds of CSPICE" begin
        qdec = gf_decreasing(d1)
        cn = Window([0.0, 10.0])
        cases = (
            (:invalidstep, ("=", 0.0, 0.0, 0.0)),
            (:invalidstep, ("=", 0.0, 0.0, -1.0)),
            (:notrecognized, ("==", 0.0, 0.0, 0.3)),
            (:notrecognized, ("MAX", 0.0, 0.0, 0.3)),
            (:valueoutofrange, ("ABSMAX", 0.0, -0.1, 0.3)),
            (:invalidvalue, ("<", 0.0, 0.1, 0.3)),
            (:invalidvalue, ("LOCMAX", 0.0, 0.1, 0.3)),
        )
        for (kind, (rel, refval, adjust, step)) in cases
            @test throws_kind(() -> gfuds(f1, qdec, rel, refval, adjust, step, cn), kind)
            @test GFO.kind(() -> GFO.gfuds(f1, qdec, rel, refval, adjust, step, collect(cn))) === kind
        end
        @test throws_kind(() -> gfudb(GFState(t -> true), 0.0, cn), :invalidstep)
        @test GFO.kind(() -> GFO.gfudb(t -> true, 0.0, collect(cn))) === :invalidstep
        @test throws_kind(() -> gfuds(f1, qdec, "=", 0.0, 0.0, 0.3, cn; tol=0.0), :invalidtolerance)
        @test throws_kind(() -> gfudb(GFState(t -> true), 0.3, cn; tol=-1.0), :invalidtolerance)
        # relations are read without regard to case and outer blanks
        @test gfuds(f1, qdec, " locmax ", 0.0, 0.0, 0.3, cn) == gfuds(f1, qdec, "LOCMAX", 0.0, 0.0, 0.3, cn)
        @test collect(gfuds(f1, qdec, " = ", 0.0, 0.0, 0.3, cn)) ==
              GFO.gfuds(f1, qdec, " = ", 0.0, 0.0, 0.3, collect(cn))
        # the confinement window is not modified
        before = copy(cn)
        gfuds(f1, qdec, "=", 0.0, 0.0, 0.3, cn)
        @test cn == before
    end

    @testset "progress report and interrupt callbacks" begin
        qdec = gf_decreasing(d1)
        cn = Window([0.0, 10.0, 12.0, 20.0])
        events = Any[]
        rep = (ev, args...) -> push!(events, (ev, args...))
        w = gfuds(f1, qdec, "<", 0.3, 0.0, 0.3, cn; report=rep)
        @test w == gfuds(f1, qdec, "<", 0.3, 0.0, 0.3, cn)
        kinds = first.(events)
        @test count(==(:init), kinds) == 2 && count(==(:final), kinds) == 2    # two passes
        @test kinds[1] === :init && kinds[end] === :final
        @test all(e -> e[1] !== :update || e[2] <= e[4] <= e[3], events)
        @test events[1][2] == cn                       # the window of the pass
        empty!(events)
        gfuds(f1, qdec, "LOCMAX", 0.0, 0.0, 0.3, cn; report=rep)
        @test count(==(:init), first.(events)) == 1    # one pass
        empty!(events)
        gfudb(GFState(t -> f1(t) > 0), 0.3, cn; report=rep)
        @test events[1][1:3] == (:init, cn, "User defined boolean event search ")
        # interrupt: the search stops at once
        calls = Ref(0)
        @test_throws InterruptException gfuds(f1, qdec, "<", 0.3, 0.0, 0.3, cn;
            bail=() -> (calls[] += 1; calls[] > 5))
        @test calls[] == 6
        @test gfuds(f1, qdec, "<", 0.3, 0.0, 0.3, cn; bail=() -> false) == w
    end

    @testset "threads keyword: the same window" begin
        qdec = gf_decreasing(d2)
        for rel in relations, nint in (1, 3, 9)
            cn = gf_random_window(rng, nint, 0.0, 60.0)
            refval = rel in ("<", "=", ">") ? 0.2 : 0.0
            adjust = rel in ("ABSMIN", "ABSMAX") ? 0.1 : 0.0
            @test gfuds(f2, qdec, rel, refval, adjust, 0.2, cn; threads=true) ==
                  gfuds(f2, qdec, rel, refval, adjust, 0.2, cn)
        end
        cn = gf_random_window(rng, 12, 0.0, 60.0)
        @test gfudb(GFState(t -> f2(t) > 0.1), 0.2, cn; threads=true) == gfudb(GFState(t -> f2(t) > 0.1), 0.2, cn)
    end

    @testset "window helpers of the solver" begin
        for _ in 1:100
            w = gf_random_window(rng, rand(rng, 1:6), 0.0, 10.0)
            amount = rand(rng, (0.0, 0.5, 2.0))
            c = GFO.cell(collect(w))
            SPICE.wnexpd!(c, amount, amount)
            @test collect(Allspice._gf_wnexpd(w, amount, amount)) == collect(Float64, c)
        end
        @test isempty(Allspice._gf_wnexpd(Window(), 0.5, 0.5))
        inner = Window([1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 7.0, 7.0])
        @test collect(Allspice._gf_wsts(inner, Window([1.0, 3.0, 5.0, 8.0]))) == [2.0, 2.0, 7.0, 7.0]
    end

    gf_report("ud")
end
