# Group 15: ellipsoid, plane and ellipse geometry against CSPICE.
#
# Feature file: features/math_ellipsoid.feature
#
# The ellipsoids have a random size from 1e-10 to 1e10 and axis ratios up to 10; some
# are spheres or spheroids and a few have a semi-axis that is zero or negative. The
# points are outside, inside, on the surface, on the axes, in the principal planes and
# at the center. The planes and the ellipses include the degenerate ones: a plane
# through the origin, an ellipse that is a segment or a point, a ray parallel to a plane.

import ForwardDiff
isdefined(@__MODULE__, :CS) || include("mathoracle.jl")

"A NeoSpice result with a missing value as the symbol `:none`."
nf(x) = x === nothing ? :none : x
"A CSPICE result with a false found flag (the last output) as the symbol `:none`."
cf(t::Tuple) = t[end] ? (length(t) == 2 ? t[1] : t[1:end-1]) : :none
"An intersection result without the points that CSPICE does not set."
pts(r::Tuple, least::Int) = r[1] >= least ? r : r[1]

group("15. math: ellipsoids, planes and ellipses", Kernels.FULL_SET_FILES...) do
    rng = MersenneTwister(SEED + 15)
    n = 400
    t = Tally()
    r3() = SVector(randn(rng), randn(rng), randn(rng))
    "Semi-axes: triaxial, spheroid or sphere, of any size; rarely not valid."
    function raxes()
        s = rand(rng) < 0.5 ? 1.0 : rmag(rng)
        u = rand(rng)
        ax = u < 0.15 ? (1.0, 1.0, 1.0) : u < 0.3 ? (1.0, 1.0, 0.1 + 2 * rand(rng)) :
             u < 0.4 ? (1.0, 0.1 + 2 * rand(rng), 1.0) : (10.0^(2 * rand(rng) - 1), 10.0^(2 * rand(rng) - 1), 10.0^(2 * rand(rng) - 1))
        v = rand(rng)
        v < 0.02 && (ax = (0.0, ax[2], ax[3]))
        v > 0.98 && (ax = (ax[1], ax[2], -ax[3]))
        return s * ax[1], s * ax[2], s * ax[3]
    end
    "Point relative to an ellipsoid: far, near, inside, on an axis, in a plane, center."
    function rpoint(a, b, c)
        m = max(abs(a), abs(b), abs(c))
        m == 0.0 && (m = 1.0)
        u = rand(rng)
        p = u < 0.35 ? m * (1.0 + 3 * rand(rng)) * normalize(r3()) :
            u < 0.55 ? m * rmag(rng)^0.3 * r3() :
            u < 0.8 ? m * rand(rng) * 0.5 * r3() : m * 2 * r3()
        w = rand(rng)
        w < 0.05 && (p = SVector(0.0, 0.0, 0.0))
        0.05 <= w < 0.10 && (p = SVector(p[1], 0.0, 0.0))
        0.10 <= w < 0.15 && (p = SVector(0.0, 0.0, p[3]))
        0.15 <= w < 0.20 && (p = SVector(0.0, p[2], 0.0))
        0.20 <= w < 0.25 && (p = SVector(p[1], p[2], 0.0))
        0.25 <= w < 0.30 && (p = SVector(0.0, p[2], p[3]))
        0.30 <= w < 0.35 && (p = SVector(p[1], 0.0, p[3]))
        return p
    end
    "Plane (validated constructor): random, through the origin, or normal to an axis."
    function rplane(scale=1.0)
        u = rand(rng)
        nrm = u < 0.15 ? rand(rng, (SVector(1.0, 0, 0), SVector(0.0, 1, 0), SVector(0.0, 0, -1))) : rmag(rng) * r3()
        k = rand(rng) < 0.15 ? 0.0 : scale * randn(rng) * vnorm(nrm)
        return nvc2pl(nrm, k)
    end
    "Ellipse from generating vectors: general, a circle, a segment or a point."
    function rellipse(scale=1.0)
        u = rand(rng)
        c = rand(rng) < 0.1 ? SVector(0.0, 0.0, 0.0) : scale * r3()
        v1 = scale * r3()
        v2 = u < 0.1 ? 2.0 * v1 : u < 0.15 ? SVector(0.0, 0.0, 0.0) : scale * r3()
        u > 0.95 && (v1 = SVector(0.0, 0.0, 0.0); v2 = v1)
        return cgv2el(c, v1, v2)
    end

    @testset "surface normal, ray intercept, nearest point" begin
        for _ in 1:2n
            a, b, c = raxes()
            p = rpoint(a, b, c)
            u = rand(rng) < 0.6 ? -p + 0.7 * max(abs(a), abs(b), abs(c)) * r3() : rvec(rng)
            versus!(t, :surfnm, () -> surfnm(a, b, c, p), () -> CS.surfnm(a, b, c, p))
            versus!(t, :surfpt, () -> nf(surfpt(p, u, a, b, c)), () -> cf(CS.surfpt(p, u, a, b, c)))
            versus!(t, :edpnt, () -> edpnt(p, a, b, c), () -> CS.edpnt(p, a, b, c))
            versus!(t, :ednmpt, () -> ednmpt(a, b, c, u), () -> CS.ednmpt(a, b, c, u))
            d = versus!(t, :nearpt, () -> nearpt(p, a, b, c), () -> CS.nearpt(p, a, b, c))
            st = vcat(p, rand(rng) < 0.1 ? SVector(0.0, 0.0, 0.0) : r3())
            versus!(t, :dnearp, () -> nf(dnearp(st, a, b, c)), () -> cf(CS.dnearp(st, a, b, c)))
            sd = vcat(u, rand(rng) < 0.1 ? SVector(0.0, 0.0, 0.0) : vnorm(u) * 1e-3 * r3())
            versus!(t, :surfpv, () -> nf(surfpv(st, sd, a, b, c)), () -> cf(CS.surfpv(st, sd, a, b, c)))
            versus!(t, :npedln, () -> npedln(a, b, c, p, u), () -> CS.npedln(a, b, c, p, u))
            if a > 0 && b > 0 && c > 0 && isfinite(d)
                m = max(a, b, c)
                np, alt = nearpt(p, a, b, c)
                level = (np[1] / a)^2 + (np[2] / b)^2 + (np[3] / c)^2
                @test abs(level - 1.0) <= 1e-12                 # on the ellipsoid
                inside = (p[1] / a)^2 + (p[2] / b)^2 + (p[3] / c)^2 < 1.0
                @test inside ? alt <= 0.0 : alt >= 0.0
                @test abs(abs(alt) - vdist(p, np)) <= 1e-13 * max(m, vnorm(p))
                x = surfpt(p, u, a, b, c)
                if x !== nothing
                    @test abs((x[1] / a)^2 + (x[2] / b)^2 + (x[3] / c)^2 - 1.0) <= 1e-9
                end
            end
        end
        # the foot of the normal: outside, the offset is along the outward normal
        p = SVector(3.0, -4.0, 5.0)
        np, alt = nearpt(p, 3.0, 2.0, 1.0)
        @test vsep(p - np, surfnm(3.0, 2.0, 1.0, np)) <= 1e-12
        @test ednmpt(3.0, 2.0, 1.0, surfnm(3.0, 2.0, 1.0, np)) ≈ np rtol = 1e-13
        @test nearpt([0.0, 0.0, 0.0], 3.0, 2.0, 1.0) == ([0.0, 0.0, 1.0], -1.0)
        @test nearpt([0.0, 0.0, 7.0], 3.0, 2.0, 1.0) == ([0.0, 0.0, 1.0], 6.0)
        @test surfpt([10.0, 0, 0], [-1.0, 0, 0], 3.0, 2.0, 1.0) == [3.0, 0.0, 0.0]
        @test surfpt([10.0, 0, 0], [1.0, 0, 0], 3.0, 2.0, 1.0) === nothing
        @test surfpt([10.0, 0, 5.0], [-1.0, 0, 0], 3.0, 2.0, 1.0) === nothing
        @test edpnt([6.0, 0.0, 0.0], 3.0, 2.0, 1.0) == [3.0, 0.0, 0.0]
        @test throws_kind(() -> nearpt(p, 3.0, 0.0, 1.0), :badaxislength)
        @test throws_kind(() -> surfnm(3.0, -2.0, 1.0, p), :badaxislength)
        @test throws_kind(() -> surfpt(p, [0.0, 0, 0], 3.0, 2.0, 1.0), :zerovector)
        @test throws_kind(() -> edpnt([0.0, 0, 0], 3.0, 2.0, 1.0), :zerovector)
        @test throws_kind(() -> ednmpt(3.0, 2.0, 1.0, [0.0, 0, 0]), :zerovector)
        @test throws_kind(() -> npedln(3.0, 2.0, 1.0, p, [0.0, 0, 0]), :zerovector)
        @test throws_kind(() -> surfpv([3.0, 0, 0, 0, 0, 0], [-1.0, 0, 0, 0, 0, 0], 3.0, 2.0, 1.0), :invalidvertex)
    end

    @testset "derivatives" begin
        a, b, c = 6378.0, 6300.0, 6356.0
        for _ in 1:40
            p = (0.2 + 3 * rand(rng)) * 6378.0 * normalize(r3())
            abs((p[1] / a)^2 + (p[2] / b)^2 + (p[3] / c)^2 - 1.0) < 0.05 && continue
            v = r3()
            # the velocity of the near point is the Jacobian of `nearpt` times the velocity
            jac = ForwardDiff.jacobian(x -> nearpt(x, a, b, c)[1], p)
            dalt = ForwardDiff.gradient(x -> nearpt(x, a, b, c)[2], p)
            dn = dnearp(vcat(p, v), a, b, c)
            @test dn !== nothing
            @test dn[1][SVector(4, 5, 6)] ≈ jac * v rtol = 1e-9
            @test dn[2][2] ≈ dalt ⋅ v rtol = 1e-9 atol = 1e-12
            # a central difference of the plain iteration
            h = 1e-3
            fd = (nearpt(p + h * v, a, b, c)[1] - nearpt(p - h * v, a, b, c)[1]) / 2h
            @test jac * v ≈ fd rtol = 1e-6 atol = 1e-8
            # the derivative with respect to the semi-axes
            da = ForwardDiff.derivative(s -> nearpt(p, a + s, b, c)[1], 0.0)
            fda = (nearpt(p, a + h, b, c)[1] - nearpt(p, a - h, b, c)[1]) / 2h
            @test da ≈ fda rtol = 1e-6 atol = 1e-8
            # the velocity of a ray intercept is the time derivative of `surfpt`
            po = 3 * 6378.0 * normalize(r3())
            u = -po + 3000.0 * r3()
            du = 10.0 * r3()
            sv = surfpv(vcat(po, v), vcat(u, du), a, b, c)
            sv === nothing && continue
            d = ForwardDiff.derivative(τ -> surfpt(po + τ * v, u + τ * du, a, b, c), 0.0)
            @test sv[SVector(4, 5, 6)] ≈ d rtol = 1e-9
            g = ForwardDiff.jacobian(x -> surfnm(a, b, c, x), p)
            @test norm(g * p) <= 1e-12            # the normal does not change along the ray from the center
            tt = NeoSpice.taylor_var(0.0, Val(2))
            x = surfpt(po .+ tt .* v, u .+ tt .* du, a, b, c)
            @test map(q -> NeoSpice.derivatives(q)[2], x) ≈ d rtol = 1e-9
            y = nearpt(p .+ tt .* v, a, b, c)[1]
            @test map(q -> NeoSpice.derivatives(q)[2], y) ≈ jac * v rtol = 1e-9
        end
    end

    @testset "planes" begin
        for _ in 1:n
            nrm, pnt = rvec(rng), rvec(rng; special=0.2)
            k = rand(rng) < 0.1 ? 0.0 : rmag(rng) * randn(rng)
            s1, s2 = rpair(rng)
            versus!(t, :nvc2pl, () -> nvc2pl(nrm, k), () -> CS.nvc2pl(nrm, k))
            versus!(t, :nvp2pl, () -> nvp2pl(nrm, pnt), () -> CS.nvp2pl(nrm, pnt))
            versus!(t, :psv2pl, () -> psv2pl(pnt, s1, s2), () -> CS.psv2pl(pnt, s1, s2))
            pl = rplane(rmag(rng))
            q = rplane()
            raw = flat(pl)
            versus!(t, :pl2nvc, () -> pl2nvc(pl), () -> CS.pl2nvc(raw))
            versus!(t, :pl2nvp, () -> pl2nvp(pl), () -> CS.pl2nvp(raw))
            versus!(t, :pl2psv, () -> pl2psv(pl), () -> CS.pl2psv(raw))
            versus!(t, :vprjp, () -> vprjp(pnt, pl), () -> CS.vprjp(pnt, raw))
            versus!(t, :vprjpi, () -> nf(vprjpi(pnt, pl, q)), () -> cf(CS.vprjpi(pnt, raw, flat(q))))
            dir = rand(rng) < 0.15 ? vcrss(pl.normal, r3()) : rvec(rng)
            vtx = rand(rng) < 0.15 ? pl2nvp(pl)[2] : pl.constant * r3() + r3()
            versus!(t, :inrypl, () -> pts(inrypl(vtx, dir, pl), 1), () -> pts(CS.inrypl(vtx, dir, raw), 1))
            versus!(t, :nplnpt, () -> nplnpt(pnt, dir, vtx), () -> CS.nplnpt(pnt, dir, vtx))
            # invariants: unit normal, constant not negative, the spanning vectors
            @test abs(vnorm(pl.normal) - 1.0) <= 4eps() && pl.constant >= 0.0
            pt, u1, u2 = pl2psv(pl)
            @test abs(u1 ⋅ u2) <= 4eps() && abs(u1 ⋅ pl.normal) <= 4eps() && abs(u2 ⋅ pl.normal) <= 4eps()
            @test abs(pt ⋅ pl.normal - pl.constant) <= 4eps() * max(1.0, pl.constant)
            pr = vprjp(pnt, pl)
            isfinite(vnorm(pr)) && @test abs(pr ⋅ pl.normal - pl.constant) <= 1e-13 * max(vnorm(pnt), pl.constant, floatmin())
        end
        @test flat(nvc2pl([0.0, 0, 2.0], -6.0)) == [0.0, 0.0, -1.0, 3.0]
        @test throws_kind(() -> nvc2pl([0.0, 0, 0], 1.0), :zerovector)
        @test throws_kind(() -> nvp2pl([0.0, 0, 0], [1.0, 2, 3]), :zerovector)
        @test throws_kind(() -> psv2pl([1.0, 2, 3], [1.0, 0, 0], [2.0, 0, 0]), :degeneratecase)
        @test throws_kind(() -> nplnpt([1.0, 2, 3], [0.0, 0, 0], [2.0, 0, 0]), :zerovector)
        @test inrypl([0.0, 0, 5.0], [0.0, 0, -1.0], nvc2pl([0.0, 0, 1.0], 2.0)) == (1, [0.0, 0.0, 2.0])
        @test inrypl([0.0, 0, 5.0], [0.0, 0, 1.0], nvc2pl([0.0, 0, 1.0], 2.0))[1] == 0
        @test inrypl([0.0, 0, 2.0], [1.0, 0, 0.0], nvc2pl([0.0, 0, 1.0], 2.0))[1] == -1
        @test vprjpi([1.0, 2, 0], nvc2pl([0.0, 0, 1.0], 0.0), nvc2pl([1.0, 0, 0.0], 3.0)) === nothing
        @test nplnpt([0.0, 0, 0], [1.0, 0, 0], [5.0, 3.0, 4.0]) == ([5.0, 0.0, 0.0], 5.0)
    end

    @testset "ellipses" begin
        for _ in 1:n
            ctr = rvec(rng)
            v1, v2 = rpair(rng)
            versus!(t, :saelgv, () -> saelgv(v1, v2), () -> CS.saelgv(v1, v2))
            versus!(t, :cgv2el, () -> cgv2el(ctr, v1, v2), () -> CS.cgv2el(ctr, v1, v2))
            el = rellipse()
            pl = rplane()
            versus!(t, :el2cgv, () -> el2cgv(el), () -> CS.el2cgv(flat(el)))
            versus!(t, :pjelpl, () -> pjelpl(el, pl), () -> CS.pjelpl(flat(el), flat(pl)))
            versus!(t, :inelpl, () -> pts(inelpl(el, pl), 1), () -> pts(CS.inelpl(flat(el), flat(pl)), 1))
            p = rand(rng) < 0.1 ? el.center : 2.0 * r3()
            versus!(t, :npelpt, () -> npelpt(p, el), () -> CS.npelpt(p, flat(el)))
            a, b, c = raxes()
            m = max(abs(a), abs(b), abs(c))
            pe = rand(rng) < 0.2 ? nvc2pl(rand(rng, (SVector(1.0, 0, 0), SVector(0.0, 1, 0), SVector(0.0, 0, 1))), m * rand(rng)) :
                 nvc2pl(r3(), m * 1.5 * rand(rng))
            versus!(t, :inedpl, () -> nf(inedpl(a, b, c, pe)), () -> cf(CS.inedpl(a, b, c, flat(pe))))
            vp = rpoint(a, b, c)
            versus!(t, :edlimb, () -> edlimb(a, b, c, vp), () -> CS.edlimb(a, b, c, vp))
            # invariants of an ellipse: orthogonal semi-axes, the major one is not shorter
            sma, smi = saelgv(v1, v2)
            if isfinite(vnorm(sma)) && vnorm(sma) < 1e150
                @test abs(sma ⋅ smi) <= 1e-13 * max(vnorm(sma)^2, floatmin())
                @test vnorm(sma) >= vnorm(smi)
            end
            e = (a > 0 && b > 0 && c > 0) ? inedpl(a, b, c, pe) : nothing
            if e !== nothing
                for th in (0.3, 2.0, 4.4)
                    x = e.center + cos(th) * e.semimajor + sin(th) * e.semiminor
                    @test abs((x[1] / a)^2 + (x[2] / b)^2 + (x[3] / c)^2 - 1.0) <= 1e-10
                    @test abs(x ⋅ pe.normal - pe.constant) <= 1e-12 * m
                end
            end
            if a > 0 && b > 0 && c > 0 && (vp[1] / a)^2 + (vp[2] / b)^2 + (vp[3] / c)^2 > 1.0 + 1e-6
                lb = edlimb(a, b, c, vp)
                for th in (0.3, 2.0, 4.4)
                    x = lb.center + cos(th) * lb.semimajor + sin(th) * lb.semiminor
                    @test abs((x[1] / a)^2 + (x[2] / b)^2 + (x[3] / c)^2 - 1.0) <= 1e-9
                    # the line of sight is tangent to the surface at a limb point
                    @test abs(x[1] * vp[1] / a^2 + x[2] * vp[2] / b^2 + x[3] * vp[3] / c^2 - 1.0) <= 1e-9
                end
            end
        end
        # a plane that is not valid, an ellipse with semi-axes that are not orthogonal
        el = cgv2el([0.0, 0, 0], [2.0, 0, 0], [0.0, 1, 0])
        bad = SpicePlane(SVector(0.0, 0.0, 1.0), -1.0)
        skew = SpiceEllipse(SVector(0.0, 0, 0), SVector(2.0, 0, 0), SVector(1.0, 1, 0))
        versus!(t, :inelpl, () -> pts(inelpl(el, bad), 1), () -> pts(CS.inelpl(flat(el), flat(bad)), 1))
        versus!(t, :inelpl, () -> pts(inelpl(skew, nvc2pl([1.0, 0, 0], 0.5)), 1), () -> pts(CS.inelpl(flat(skew), [1.0, 0, 0, 0.5]), 1))
        @test throws_kind(() -> inelpl(skew, nvc2pl([1.0, 0, 0], 0.5)), :invalidellipse)
        @test inelpl(el, nvc2pl([1.0, 0, 0], 1.0))[1] == 2
        @test inelpl(el, nvc2pl([1.0, 0, 0], 2.0))[1] == 1
        @test inelpl(el, nvc2pl([1.0, 0, 0], 3.0))[1] == 0
        @test inelpl(el, nvc2pl([0.0, 0, 1.0], 0.0))[1] == -1
        @test npelpt([5.0, 0, 3.0], el)[1] == [2.0, 0.0, 0.0]
        @test npelpt([5.0, 0, 3.0], el)[2] ≈ sqrt(18.0) rtol = 4eps()
        @test inedpl(3.0, 2.0, 1.0, nvc2pl([0.0, 0, 1.0], 2.0)) === nothing
        e = inedpl(3.0, 2.0, 1.0, nvc2pl([0.0, 0, 1.0], 0.0))
        @test vnorm(e.semimajor) ≈ 3.0 && vnorm(e.semiminor) ≈ 2.0 && vnorm(e.center) == 0.0
        @test throws_kind(() -> edlimb(3.0, 2.0, 1.0, [0.1, 0, 0]), :degeneratecase)
        @test throws_kind(() -> edlimb(3.0, 0.0, 1.0, [9.0, 0, 0]), :degeneratecase)
    end

    @testset "terminator" begin
        oracle_load(Kernels.FULL_SET_FILES)
        ctx = neo_context(Kernels.FULL_SET_FILES)
        cases = (("EARTH", "IAU_EARTH"), ("MOON", "IAU_MOON"), ("MOON", "MOON_PA_DE440"), ("399", "ITRF93"),
            ("SUN", "IAU_SUN"), ("EARTH", "IAU_MOON"), ("NOBODY", "IAU_EARTH"), ("EARTH", "NOFRAME"))
        for _ in 1:150
            target, fixref = rand(rng, cases)
            trmtyp = rand(rng, ("UMBRAL", "PENUMBRAL", " umbral", "Penumbral ", "SHADOW"))
            source = rand(rng, ("SUN", "10", "SUN", "MOON", "NOBODY"))
            obs = rand(rng, ("EARTH", "MOON", "SUN", "EARTH BARYCENTER", "NOBODY", "3"))
            abcorr = rand(rng, ("NONE", "LT", "LT+S", "CN", "CN+S", " lt+s", "XLT", "BAD"))
            et = 6.0e8 * (2 * rand(rng) - 1)
            npts = rand(rng, (1, 2, 7, 40, 0, -3))
            versus!(t, :edterm, () -> edterm(ctx, trmtyp, source, target, et, fixref, abcorr, obs, npts),
                () -> CS.edterm(trmtyp, source, target, et, fixref, abcorr, obs, npts))
        end
        trgepc, obspos, trm = edterm(ctx, "UMBRAL", "SUN", "MOON", 1.0e8, "IAU_MOON", "LT+S", "EARTH", 12)
        @test length(trm) == 12 && trm isa Vector{SVector{3,Float64}}
        rad = bodvrd(ctx, "MOON", "RADII")
        sun, _ = spkpos(ctx, "SUN", trgepc, "IAU_MOON", "LT+S", "MOON")
        rsun = maximum(bodvrd(ctx, "SUN", "RADII"))
        for x in trm
            @test abs((x[1] / rad[1])^2 + (x[2] / rad[2])^2 + (x[3] / rad[3])^2 - 1.0) <= 1e-12
            # at a point of the umbral terminator the upper limb of the Sun is on the
            # horizon: the center of the Sun is below it by the angular radius
            nrm = surfnm(rad[1], rad[2], rad[3], x)
            @test asin(nrm ⋅ normalize(sun - x)) ≈ -asin(rsun / vnorm(sun - x)) rtol = 1e-9
        end
        @test trgepc < 1.0e8
        @test throws_kind(() -> edterm(ctx, "SHADOW", "SUN", "MOON", 0.0, "IAU_MOON", "NONE", "EARTH", 3), :notsupported)
        @test throws_kind(() -> edterm(ctx, "UMBRAL", "SUN", "MOON", 0.0, "IAU_EARTH", "NONE", "EARTH", 3), :invalidfixref)
        @test throws_kind(() -> edterm(ctx, "UMBRAL", "SUN", "MOON", 0.0, "IAU_MOON", "NONE", "EARTH", 0), :invalidsize)
        SPICE.kclear()
    end

    @testset "differences against CSPICE" begin
        # The ports follow the operation order of CSPICE, so most results are equal bit
        # for bit; the others have a sine, a cosine or an arc tangent of the C library.
        report!(t, "ellipsoid"; tol=0.0, tols=Dict(:ednmpt => 4eps(), :surfpv => 4eps(),
            :inelpl => 8eps(), :edterm => 4eps()))
    end

    @testset "inferred types and no allocation for static arguments" begin
        p = SVector(3.0, -4.0, 5.0)
        u = SVector(-3.0, 4.2, -5.1)
        st = SVector(3.0, -4.0, 5.0, 0.1, 0.2, 0.3)
        sd = SVector(-3.0, 4.2, -5.1, 0.01, 0.0, 0.02)
        a, b, c = 3.0, 2.0, 1.0
        pl = nvc2pl(SVector(1.0, 2.0, 3.0), 0.7)
        q = nvc2pl(SVector(1.0, -2.0, 2.0), 0.2)
        el = cgv2el(SVector(0.1, 0.2, 0.3), SVector(1.0, 1.0, 0.0), SVector(0.0, 0.5, 2.0))
        V3, V6, P, E = SVector{3,Float64}, SVector{6,Float64}, SpicePlane{Float64}, SpiceEllipse{Float64}
        calls = Any[(surfnm, (a, b, c, p), V3), (surfpt, (p, u, a, b, c), V3), (edpnt, (p, a, b, c), V3),
            (ednmpt, (a, b, c, u), V3), (nearpt, (p, a, b, c), Tuple{V3,Float64}),
            (dnearp, (st, a, b, c), Tuple{V6,SVector{2,Float64}}), (surfpv, (st, sd, a, b, c), V6),
            (npedln, (a, b, c, p, u), Tuple{V3,Float64}), (nvc2pl, (u, 2.0), P), (nvp2pl, (u, p), P),
            (psv2pl, (p, u, SVector(1.0, 0.0, 0.0)), P), (pl2nvc, (pl,), Tuple{V3,Float64}),
            (pl2nvp, (pl,), Tuple{V3,V3}), (pl2psv, (pl,), Tuple{V3,V3,V3}), (vprjp, (p, pl), V3),
            (vprjpi, (p, pl, q), V3), (inrypl, (p, u, pl), Tuple{Int,V3}), (nplnpt, (p, u, st[SVector(4, 5, 6)]), Tuple{V3,Float64}),
            (saelgv, (p, u), Tuple{V3,V3}), (cgv2el, (p, u, SVector(1.0, 0.0, 0.0)), E), (el2cgv, (el,), Tuple{V3,V3,V3}),
            (pjelpl, (el, pl), E), (inelpl, (el, pl), Tuple{Int,V3,V3}), (npelpt, (p, el), Tuple{V3,Float64}),
            (inedpl, (a, b, c, pl), E), (edlimb, (a, b, c, p), E)]
        for (fn, args, T) in calls
            @test fn(args...) isa T
            @test Base.promote_op(fn, map(typeof, args)...) <: Union{T,Nothing}
            @test alloc_of(fn, args...) == 0
        end
    end
end
