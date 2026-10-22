# Conic types: type 5 (discrete states, two-body propagation), type 15 (precessing
# conic) and type 17 (equinoctial elements).

@noinline _conic_error(msg::String, kind::Symbol=:badorbit) = throw(SpiceError(kind, msg))

"Reciprocals `1/(i (i+1))`, `i = 1:20`, of the series of `stmp03`."
const _STUMPFF_PAIRS = ntuple(i -> 1.0 / (Float64(i) * Float64(i + 1)), 20)

"Lowest argument of `_stumpff`: below it `cosh` overflows."
const _STUMPFF_LBOUND = -(log(2.0) + log(floatmax(Float64)))^2

"""
    _stumpff(x::T, kind=:badorbit) -> (c0, c1, c2, c3)

Stumpff functions of `x`, as `stmp03`: the closed forms with `cosh`, `sinh` for
`x < -1` and with `cos`, `sin` for `x > 1`, and the nested series of 20 terms
between. The branch follows `primal(x)`. Throws `SpiceError(kind)` if
`x <= _STUMPFF_LBOUND`.
"""
@inline function _stumpff(x::T, kind::Symbol=:badorbit) where {T<:Real}
    x0 = primal(x)
    x0 <= _STUMPFF_LBOUND &&
        _conic_error("the argument of the Stumpff functions is too low", kind)
    if x0 < -1.0
        z = sqrt(-x)
        c0 = cosh(z)
        c1 = sinh(z) / z
        return c0, c1, (1 - c0) / x, (1 - c1) / x
    elseif x0 > 1.0
        z = sqrt(x)
        s, c0 = sincos(z)
        c1 = s / z
        return c0, c1, (1 - c0) / x, (1 - c1) / x
    end
    p = _STUMPFF_PAIRS
    c3 = one(T)
    for i in 20:-2:4          # invariant: c3 is the nested sum of the terms above i
        c3 = 1.0 - x * p[i] * c3
    end
    c3 = p[2] * c3
    c2 = one(T)
    for i in 19:-2:3
        c2 = 1.0 - x * p[i] * c2
    end
    c2 = p[1] * c2
    return 1.0 - x * c2, 1.0 - x * c3, c2, c3
end

"""
    _prop2b(gm, pos, vel, dt::T, kind=:badorbit) -> (pos::Vec3{T}, vel::Vec3{T})

Two-body state at the time `dt` after the state `(pos, vel)` about a mass with
gravitational parameter `gm`, by the universal-variable formulation of `prop2b`.

Contract: `gm > 0`, `pos` and `vel` are not zero and not parallel, else
`SpiceError(kind)`, which is also the error for a time interval that is too long. The
universal variable `x` is the plain `Float64` root of the Kepler equation
`kfun(x) = primal(dt)`, found by the bracketing and bisection of CSPICE, so a plain `dt`
gives the CSPICE state. For another number type `x`
gets its derivative parts from Newton steps with `d kfun / dx = b r` (see
`_with_primal`), and the closing formulas are evaluated in that type.
"""
function _prop2b(gm::Float64, pos::Vec3{Float64}, vel::Vec3{Float64}, dt::T,
        kind::Symbol=:badorbit) where {T<:Real}
    gm > 0.0 || _conic_error("the gravitational parameter is not positive", kind)
    (_maxabs(pos) == 0.0 || _maxabs(vel) == 0.0) &&
        _conic_error("the initial position or velocity is zero", kind)
    r0 = _vnorm(pos)
    rv = _dot(pos, vel)
    hvec = _cross(pos, vel)
    h2 = _dot(hvec, hvec)
    h2 == 0.0 && _conic_error("the initial position and velocity are parallel", kind)
    eqvec = (1.0 / gm) * _cross(vel, hvec) + (-1.0 / r0) * pos
    e = _vnorm(eqvec)
    q = h2 / (gm * (e + 1.0))
    f = 1.0 - e
    b = sqrt(q / gm)
    br0 = b * r0
    b2rv = b * b * rv
    bq = b * q
    qovr0 = q / r0
    # Bound of |x| that keeps the Stumpff functions and the Kepler function finite.
    maxc = max(1.0, abs(br0), abs(b2rv), abs(bq), abs(qovr0 / bq))
    bound = if f < 0.0
        fixed = log(floatmax(Float64) / 2.0) - log(maxc)
        rootf = sqrt(-f)
        min(fixed / rootf, (fixed + log(-f) * 1.5) / rootf)
    else
        exp((log(1.5) + log(floatmax(Float64)) - log(maxc)) / 3.0)
    end

    t = primal(dt)
    x = 0.0
    if t != 0.0
        x = clamp(t / bq, -bound, bound)
        _, c1, c2, c3 = _stumpff(f * x * x, kind)
        kfun = x * (br0 * c1 + x * (b2rv * c2 + x * (bq * c3)))
        lower = upper = 0.0
        if t < 0.0
            lower = x
            # invariant: kfun(upper) > t; `lower` doubles until kfun(lower) <= t
            while kfun > t
                upper = lower
                lower *= 2.0
                oldx = x
                x = clamp(lower, -bound, bound)
                x == oldx && _conic_error("the time interval is too long for two-body " *
                    "propagation", kind)
                _, c1, c2, c3 = _stumpff(f * x * x, kind)
                kfun = x * (br0 * c1 + x * (b2rv * c2 + x * (bq * c3)))
            end
        else
            upper = x
            # invariant: kfun(lower) < t; `upper` doubles until kfun(upper) >= t
            while kfun < t
                lower = upper
                upper *= 2.0
                oldx = x
                x = clamp(upper, -bound, bound)
                x == oldx && _conic_error("the time interval is too long for two-body " *
                    "propagation", kind)
                _, c1, c2, c3 = _stumpff(f * x * x, kind)
                kfun = x * (br0 * c1 + x * (b2rv * c2 + x * bq * c3))
            end
        end
        # Bisection. The count starts again at 64 when both ends are away from zero.
        x = min(upper, max(lower, (lower + upper) / 2.0))
        _, c1, c2, c3 = _stumpff(f * x * x, kind)
        lcount = 0
        mostc = 1000
        # invariant: kfun(lower) <= t <= kfun(upper)
        while x > lower && x < upper && lcount < mostc
            kfun = x * (br0 * c1 + x * (b2rv * c2 + x * bq * c3))
            if kfun > t
                upper = x
            elseif kfun < t
                lower = x
            else
                upper = x
                lower = x
            end
            if mostc > 64 && upper != 0.0 && lower != 0.0
                mostc = 64
                lcount = 0
            end
            x = min(upper, max(lower, (lower + upper) / 2.0))
            _, c1, c2, c3 = _stumpff(f * x * x, kind)
            lcount += 1
        end
    elseif T === Float64
        return pos, vel
    end

    # Derivative parts of x: Newton steps on kfun(x) - dt with its plain value removed.
    xt = convert(T, x)
    for _ in 1:_ad_order(T)       # invariant: one more derivative order of xt is exact
        d0, d1, d2, d3 = _stumpff(f * xt * xt, kind)
        kf = xt * (br0 * d1 + xt * (b2rv * d2 + xt * bq * d3))
        brt = br0 * d0 + xt * (b2rv * d1 + xt * (bq * d2))
        xt -= _with_primal(kf - dt, 0.0) / brt
    end
    c0, c1, c2, c3 = _stumpff(f * xt * xt, kind)
    x2 = xt * xt
    x3 = x2 * xt
    br = br0 * c0 + xt * (b2rv * c1 + xt * (bq * c2))
    pc = 1.0 - qovr0 * x2 * c2
    vc = dt - bq * x3 * c3
    pcdot = -(qovr0 / br) * xt * c1
    vcdot = 1.0 - bq / br * x2 * c2
    return pc * pos + vc * vel, pcdot * pos + vcdot * vel
end

"""
    _spk5(data, seg, et::T) -> (pos, vel)

State of a type 5 segment, as `spkr05` and `spke05`. The segment holds `n` states,
their epochs, an epoch directory, the gravitational parameter and `n`. Between two
epochs the state is the weighted mean of the two-body propagations of the two
neighbor states, with the weight `(1 + cos(pi (t - t1)/(t2 - t1)))/2` for the
earlier one, and the velocity has the term of the weight rate. Before the first
epoch and after the last one the state is the propagation of the end state alone.
"""
function _spk5(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    len >= 9 || _badsegment(seg, "the segment is too short")
    @inbounds gm = data[e-1]
    n = _seg_count(seg, @inbounds(data[e]), 1, _MAXCOUNT, "state count")
    7 * n + n ÷ 100 + 2 == len || _badsegment(seg, "the state count does not agree with the length")
    ebase = b - 1 + 6 * n                       # epoch i is data[ebase+i]
    i = _lstltd(data, ebase, n, primal(et))
    if i == 0 || i == n
        k = max(i, 1)
        a = b + 6 * (k - 1)
        return _prop2b(gm, _vec3(data, a), _vec3(data, a + 3), et - @inbounds(data[ebase+k]))
    end
    a = b + 6 * (i - 1)
    @inbounds t1, t2 = data[ebase+i], data[ebase+i+1]
    p1, v1 = _prop2b(gm, _vec3(data, a), _vec3(data, a + 3), et - t1)
    p2, v2 = _prop2b(gm, _vec3(data, a + 6), _vec3(data, a + 9), et - t2)
    denom = t2 - t1
    arg = (et - t1) * _PI / denom
    dargdt = _PI / denom
    s, c = sincos(arg)
    w = c * 0.5 + 0.5
    dwdt = s * -0.5 * dargdt
    pos = w * p1 + (1.0 - w) * p2
    vel = (w * v1 + (1.0 - w) * v2) + (dwdt * p1 + (-dwdt) * p2)
    return pos, vel
end

"""
    _spk15(data, seg, et::T) -> (pos, vel)

State of a type 15 segment (precessing conic), as `spke15`. The segment is one
record: epoch of periapsis, unit trajectory pole TP, unit periapsis vector PA,
semi-latus rectum P, eccentricity, J2 flag, central body pole PV, GM, J2 and the
equatorial radius of the central body.

The state is the two-body state from periapsis. For an ellipse outside the body
with J2 not zero it is then turned about TP by the drift of the periapsis (unless
the flag is 1) and about PV by the drift of the node (unless the flag is 2); a flag
of 3 gives no turn. Both drifts are in proportion to the true anomaly counted from
the periapsis epoch. CSPICE turns the two-body velocity only, so the velocity is
not the derivative of the position; the derivative of the true anomaly comes from
its `atan` form in the orbit plane.
"""
function _spk15(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    b = seg.baddr
    seg.eaddr - b + 1 == 16 || _badsegment(seg, "the segment does not have 16 doubles")
    @inbounds begin
        epoch = data[b]
        tp, pa = _vec3(data, b + 1), _vec3(data, b + 4)
        p, ecc = data[b+7], data[b+8]
        j2flg = _trunc_index(data[b+9])
        pv = _vec3(data, b + 10)
        gm, oj2, rpl = data[b+13], data[b+14], data[b+15]
    end
    (p > 0.0 && ecc >= 0.0 && gm > 0.0 && rpl >= 0.0) || _badsegment(seg,
        "P = $p, ECC = $ecc, GM = $gm, radius = $rpl are not all in their ranges")
    (_maxabs(tp) == 0.0 || _maxabs(pa) == 0.0 || _maxabs(pv) == 0.0) &&
        _badsegment(seg, "the trajectory pole, the periapsis vector or the body pole is zero")
    pa, tp, pv = _vhat(pa), _vhat(tp), _vhat(pv)
    abs(_dot(pa, tp)) > 1.0e-5 &&
        _badsegment(seg, "the periapsis vector and the trajectory pole are not orthogonal")
    near = p / (ecc + 1.0)
    speed = sqrt(gm / p) * (ecc + 1.0)
    q = _cross(tp, pa)
    dt = et - epoch
    pos, vel = _prop2b(gm, near * pa, speed * q, dt)
    (j2flg != 3 && oj2 != 0.0 && ecc < 1.0 && near > rpl) || return pos, vel

    # Mean anomaly since periapsis, split in a whole number of revolutions and the
    # rest in [-pi, pi].
    oneme2 = 1.0 - ecc * ecc
    dmdt = oneme2 / p * sqrt(gm * oneme2 / p)
    twopi = 2.0 * _PI
    manom = dmdt * primal(dt)
    theta = rem(manom, twopi)
    if abs(theta) > _PI
        theta -= _dsign(twopi, theta)
    end
    k2pi = manom - theta
    pos0 = map(primal, pos)
    ta0 = _dsign(_vsep(pa, pos0), theta) + k2pi
    # The same angle with derivative parts: components of the position along PA and
    # along the in-plane normal to PA.
    ta = _with_primal(atan(_dot(pos, _vhat(q)), _dot(pos, pa)), ta0)
    cosinc = _dot(pv, tp)
    d = rpl / p
    z = ta * 1.5 * oj2 * (d * d)
    dnode = -z * cosinc
    dperi = z * (cosinc * cosinc * 2.5 - 0.5)
    if j2flg != 1
        pos, vel = _vrotv(pos, tp, dperi), _vrotv(vel, tp, dperi)
    end
    if j2flg != 2
        pos, vel = _vrotv(pos, pv, dnode), _vrotv(vel, pv, dnode)
    end
    return pos, vel
end

"""
    _kepler_equinoctial(h::Float64, k::Float64) -> Float64

Root `x` of `x = h cos(x) + k sin(x)` for `h² + k² < 1`, by the bisection and the
five Newton steps of `kpsolv`. Throws `SpiceError(:badorbit)` if `h² + k² >= 1`.
"""
function _kepler_equinoctial(h::Float64, k::Float64)
    ecc2 = h * h + k * k
    ecc2 >= 1.0 && _conic_error("the equinoctial Kepler equation has an eccentricity >= 1")
    y0 = -h
    ecc = sqrt(ecc2)
    xl, xu = if y0 > 0.0
        -ecc, 0.0
    elseif y0 < 0.0
        0.0, ecc
    else
        return 0.0
    end
    maxit = min(32, max(1, _round_index(1.0 / (1.0 - ecc))))
    xm = 0.0
    for _ in 1:maxit          # invariant: the root is in [xl, xu]
        xm = max(xl, min(xu, (xl + xu) * 0.5))
        yxm = xm - h * cos(xm) - k * sin(xm)
        if yxm > 0.0
            xu = xm
        else
            xl = xm
        end
    end
    x = xm
    for _ in 1:5
        cosx = cos(x)
        sinx = sin(x)
        yx = x - h * cosx - k * sinx
        ypx = h * sinx + 1.0 - k * cosx
        x -= yx / ypx
    end
    return x
end

"`x - y * trunc(x / y)`, the remainder of `d_mod` of the f2c library."
@inline function _f2c_dmod(x::Float64, y::Float64)
    q = x / y
    return x - y * (q >= 0.0 ? floor(q) : -floor(-q))
end

"""
    _eqn_state(et::T, epoch, a, h0, k0, l, p0, q0, dlpdt, mldt, nodedt, rapol, decpol)
        -> (pos::Vec3{T}, vel::Vec3{T})

State at `et` from equinoctial elements, as `eqncpv`, without its checks of the
semi-major axis and of the eccentricity. The elements are moved to `et` with their
rates, the Kepler equation is solved in its equinoctial form (plain root as `kpsolv`,
derivative parts by Newton steps, see `_with_primal`), and the state in the equatorial
frame of the pole `(rapol, decpol)` is turned to the frame of the pole's coordinates.
Throws `SpiceError(:eccoutofbounds)` if `h² + k² >= 0.81` at `et`, as `kepleq`.
"""
function _eqn_state(et::T, epoch::Float64, a::Float64, h0::Float64, k0::Float64, l::Float64,
        p0::Float64, q0::Float64, dlpdt::Float64, mldt::Float64, nodedt::Float64,
        rapol::Float64, decpol::Float64) where {T<:Real}
    sa, ca = sincos(rapol)
    sd, cd = sincos(decpol)
    dt = et - epoch
    # Elements at `et`: H, K turn with the periapsis, P, Q with the node.
    san, can = sincos(dt * dlpdt)
    h = h0 * can + k0 * san
    k = k0 * can - h0 * san
    sn, cn = sincos(dt * nodedt)
    p = p0 * cn + q0 * sn
    q = q0 * cn - p0 * sn
    prate = dlpdt - nodedt
    b = sqrt(1.0 - h * h - k * k)
    b = 1.0 / (b + 1.0)
    di = 1.0 / (p * p + 1.0 + q * q)
    vf = Vec3((1.0 - p * p + q * q) * di, p * 2.0 * q * di, p * -2.0 * di)
    vg = Vec3(p * 2.0 * q * di, (p * p + 1.0 - q * q) * di, q * 2.0 * di)
    # Mean longitude, with the whole revolutions of the rate term removed.
    mlrate = mldt * dt
    ml = l + _with_primal(mlrate, _f2c_dmod(primal(mlrate), 2.0 * _PI))
    primal(h)^2 + primal(k)^2 >= 0.81 && throw(SpiceError(:eccoutofbounds,
        "the eccentricity of the equinoctial elements at the epoch is not below 0.9"))
    # Eccentric longitude: ml plus the root of the Kepler equation in the elements
    # turned by ml (`kepleq`).
    sml, cml = sincos(ml)
    ev1 = -h * cml + k * sml
    ev2 = h * sml + k * cml
    x = convert(T, _kepler_equinoctial(primal(ev1), primal(ev2)))
    for _ in 1:_ad_order(T)       # invariant: one more derivative order of x is exact
        sx, cx = sincos(x)
        yx = x - ev1 * cx - ev2 * sx
        ypx = ev1 * sx + 1.0 - ev2 * cx
        x -= _with_primal(yx, 0.0) / ypx
    end
    sf, cf = sincos(ml + x)
    x1 = a * ((1.0 - b * (h * h)) * cf + (h * k * b * sf - k))
    y1 = a * ((1.0 - b * (k * k)) * sf + (h * k * b * cf - h))
    rb = h * sf + k * cf
    r = a * (1.0 - rb)
    ra = mldt * a * a / r
    dx1 = ra * (-sf + h * b * rb)
    dy1 = ra * (cf - k * b * rb)
    nfac = 1.0 - dlpdt / mldt
    dx = nfac * dx1 - prate * y1
    dy = nfac * dy1 + prate * x1
    xh = x1 * vf + y1 * vg
    temp = Vec3(-nodedt * xh[2], nodedt * xh[1], zero(xh[1]))
    vh = temp + dx * vf + dy * vg
    # Rows of the matrix from the frame of the pole to the frame of the segment.
    r1 = Vec3(-sa, -ca * sd, ca * cd)
    r2 = Vec3(ca, -sa * sd, sa * cd)
    r3 = Vec3(0.0, cd, sd)
    return Vec3(_dot(r1, xh), _dot(r2, xh), _dot(r3, xh)),
        Vec3(_dot(r1, vh), _dot(r2, vh), _dot(r3, vh))
end

"""
    _spk17(data, seg, et::T) -> (pos, vel)

State of a type 17 segment (equinoctial elements), as `spke17` with `eqncpv`
(`_eqn_state`). The segment is one record: epoch, semi-major axis, H, K, mean
longitude, P, Q, rate of the longitude of periapsis, mean longitude rate, rate of the
node, and the right ascension and declination of the pole of the reference plane.
Throws `SpiceError(:baddaf)` unless the semi-major axis is positive and the
eccentricity is not above 0.9.
"""
function _spk17(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    bs = seg.baddr
    seg.eaddr - bs + 1 == 12 || _badsegment(seg, "the segment does not have 12 doubles")
    @inbounds begin
        epoch = data[bs]
        a, h0, k0, l = data[bs+1], data[bs+2], data[bs+3], data[bs+4]
        p0, q0 = data[bs+5], data[bs+6]
        dlpdt, mldt, nodedt = data[bs+7], data[bs+8], data[bs+9]
        rapol, decpol = data[bs+10], data[bs+11]
    end
    (a > 0.0 && sqrt(h0 * h0 + k0 * k0) <= 0.9) || _badsegment(seg,
        "the semi-major axis is not positive or the eccentricity is above 0.9")
    return _eqn_state(et, epoch, a, h0, k0, l, p0, q0, dlpdt, mldt, nodedt, rapol, decpol)
end
