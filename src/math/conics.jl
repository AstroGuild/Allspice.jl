# Two-body (conic) calls of CSPICE: Stumpff functions, propagation of a state, and the
# conversions between a state and osculating elements.
#
# Elements are in the CSPICE order: perifocal distance `rp`, eccentricity `ecc`,
# inclination `inc`, longitude of the ascending node `lnode`, argument of periapsis
# `argp`, mean anomaly at epoch `m0`, epoch `t0` and gravitational parameter `mu`.
# The propagator is the private `_prop2b` of `spk_types/conic.jl`, which SPK types 5
# and 15 use; nothing of it is repeated here.

"""
    stmp03(x) -> (c0, c1, c2, c3)

Return the Stumpff functions `c0` to `c3` of `x`, as `stmp03_`: closed forms with
hyperbolic functions for `x < -1` and with circular functions for `x > 1`, and a
series of 20 terms between. Requires `x` above `-(log(2) + log(floatmax))^2`
(`SpiceError(:valueoutofrange)`), where `c0` would overflow. Postcondition:
`c0 = 1 - x c2` and `c1 = 1 - x c3` to rounding.
"""
function stmp03(x::Real)
    xx = _m_f(x)
    primal(xx) <= _STUMPFF_LBOUND && throw(SpiceError(:valueoutofrange,
        "the argument $(primal(xx)) of the Stumpff functions is below the lower bound $_STUMPFF_LBOUND"))
    return _stumpff(xx)
end

"""
    prop2b(gm, pvinit, dt) -> SVector{6}

Return the state at the time `dt` after the state `pvinit` of a body that moves about a
mass with the gravitational parameter `gm` with no other force, as `prop2b_c`
(universal variables; the Kepler equation is solved by the bracketing and bisection of
CSPICE). `gm` and `pvinit` are plain numbers; `dt` may be a `ForwardDiff` or `Taylor`
number, and the result then carries the derivatives with respect to `dt`. Requires
`gm > 0` (`SpiceError(:nonpositivemass)`), a position and a velocity that are not zero
(`SpiceError(:zeroposition)`, `SpiceError(:zerovelocity)`) and not parallel
(`SpiceError(:nonconicmotion)`); `SpiceError(:dtoutofrange)` if `dt` is too long for
the arithmetic. Postcondition: the specific energy and the angular momentum are kept
to rounding.
"""
function prop2b(gm::Real, pvinit::AbstractVector, dt::Real)
    s = _m_v6(pvinit)
    mu = Float64(gm)
    pos = SVector{3,Float64}(s[1], s[2], s[3])
    vel = SVector{3,Float64}(s[4], s[5], s[6])
    mu > 0.0 || throw(SpiceError(:nonpositivemass, "the gravitational parameter is $mu; it must be positive"))
    _m_iszero(pos) && throw(SpiceError(:zeroposition, "the initial position is the zero vector"))
    _m_iszero(vel) && throw(SpiceError(:zerovelocity, "the initial velocity is the zero vector"))
    h = _m_cross(pos, vel)
    _m_dot(h, h) == 0.0 && throw(SpiceError(:nonconicmotion,
        "the initial position and velocity are parallel: the motion is not conic"))
    p, v = _m_prop2b(mu, pos, vel, _m_f(dt))
    return _m_state(p, v)
end

"`_prop2b` with its remaining failure reported as `SpiceError(:dtoutofrange)`."
@inline function _m_prop2b(mu::Float64, pos::SVector{3,Float64}, vel::SVector{3,Float64}, dt::Real)
    try
        return _prop2b(mu, pos, vel, dt)
    catch e
        (e isa SpiceError && e.kind === :badorbit) || rethrow()
        throw(SpiceError(:dtoutofrange, "the time interval is too long for two-body propagation: " * e.msg))
    end
end

"""
    conics(elts, et) -> SVector{6}

Return the state at the epoch `et` of a body on the conic orbit with the elements
`elts = (rp, ecc, inc, lnode, argp, m0, t0, mu)`, as `conics_c`. `et` may be a
`ForwardDiff` or `Taylor` number. Requires `ecc >= 0`
(`SpiceError(:badeccentricity)`), `rp > 0` (`SpiceError(:badperiapsevalue)`) and
`mu > 0` (`SpiceError(:badgm)`). Postcondition: `conics(oscelt(state, et, mu), et)` is
`state` to rounding.
"""
function conics(elts::AbstractVector, et::Real)
    length(elts) == 8 || _m_baddim("a set of conic elements", 8)
    rp, ecc, inc, lnode, argp, m0, t0, mu = ntuple(i -> Float64(elts[begin+i-1]), Val(8))
    ecc >= 0.0 || throw(SpiceError(:badeccentricity, "the eccentricity is $ecc; it must not be negative"))
    rp > 0.0 || throw(SpiceError(:badperiapsevalue, "the perifocal distance is $rp; it must be positive"))
    mu > 0.0 || throw(SpiceError(:badgm, "the gravitational parameter is $mu; it must be positive"))
    sini, cosi = sincos(inc)
    sinn, cosn = sincos(lnode)
    sinw, cosw = sincos(argp)
    snci = sinn * cosi
    cnci = cosn * cosi
    basisp = SVector(cosn * cosw - snci * sinw, sinn * cosw + cnci * sinw, sini * sinw)
    basisq = SVector(-cosn * sinw - snci * cosw, -sinn * sinw + cnci * cosw, sini * cosw)
    v = sqrt(mu * (ecc + 1.0) / rp)
    t = _m_f(et)
    if ecc < 1.0
        ainvrs = (1.0 - ecc) / rp
        n = sqrt(mu * ainvrs) * ainvrs
        period = _M_TWOPI / n
        w = t - t0 + m0 / n
        # the remainder of `d_mod`: the subtracted number of periods is a constant
        q = primal(w) / period
        dt = w - period * (q >= 0.0 ? floor(q) : -floor(-q))
    elseif ecc > 1.0
        ainvrs = (ecc - 1.0) / rp
        n = sqrt(mu * ainvrs) * ainvrs
        dt = t - t0 + m0 / n
    else
        n = sqrt(mu / (rp * 2.0)) / rp
        dt = t - t0 + m0 / n
    end
    p, vv = _m_prop2b(mu, rp * basisp, v * basisq, dt)
    return _m_state(p, vv)
end

"`|a|` with the sign of `b`, as `d_sign` of the f2c library (a zero `b` counts as positive)."
@inline _m_dsign(a, b) = (x = _m_abs(a); primal(b) >= 0.0 ? x : -x)

@inline function _m_oscelt(state::SVector{6}, et, mu)
    primal(mu) > 0.0 || throw(SpiceError(:nonpositivemass,
        "the gravitational parameter is $(primal(mu)); it must be positive"))
    r, v = _m_pv(state)
    _m_iszero(r) && throw(SpiceError(:degeneratecase, "the position is the zero vector"))
    _m_iszero(v) && throw(SpiceError(:degeneratecase, "the velocity is the zero vector"))
    rmag = _m_norm(r)
    vmag = _m_norm(v)
    h = _m_cross(r, v)
    _m_iszero(h) && throw(SpiceError(:degeneratecase,
        "the position and the velocity are parallel: the angular momentum is zero"))
    z0 = zero(h[1])
    n = SVector(-h[2], h[1], z0)
    e = (1.0 / mu) * ((vmag * vmag - mu / rmag) * r + (-_m_dot(r, v)) * v)
    ecc = _m_norm(e)
    abs(primal(ecc) - 1.0) <= 1.0e-10 && (ecc = one(ecc))
    p = _m_dot(h, h) / mu
    rp = p / (ecc + 1.0)
    inc = vsep(h, SVector(0.0, 0.0, 1.0))
    if abs(primal(inc)) < 1.0e-10
        inc = zero(inc)
        n = SVector(one(z0), z0, z0)
    elseif abs(primal(inc) - Float64(pi)) < 1.0e-10
        inc = Float64(pi) + zero(inc)
        n = SVector(one(z0), z0, z0)
    end
    lnode = atan(n[2], n[1])
    primal(lnode) < 0.0 && (lnode += _M_TWOPI)
    argp = zero(inc)
    if primal(ecc) != 0.0
        argp = vsep(n, e) + zero(inc)
        if primal(argp) != 0.0
            if primal(inc) == 0.0 || primal(inc) == Float64(pi)
                primal(_m_dot(e, _m_ucrss(h, n))) < 0.0 && (argp = _M_TWOPI - argp)
            elseif primal(e[3]) < 0.0
                argp = _M_TWOPI - argp
            end
        end
    end
    perix = primal(ecc) == 0.0 ? _m_hat(n) : _m_hat(e)
    periy = _m_ucrss(h, perix)
    nu = atan(_m_dot(r, periy), _m_dot(r, perix))
    snu, cnu = sincos(nu)
    if primal(ecc) < 1.0
        cosea = (ecc + cnu) / (ecc * cnu + 1.0)
        sinea = rmag / rp * sqrt((1.0 - ecc) / (ecc + 1.0)) * snu
        ea = atan(sinea, cosea)
        m0 = _m_dsign(ea - ecc * sin(ea), nu)
        primal(m0) < 0.0 && (m0 += _M_TWOPI)
    elseif primal(ecc) > 1.0
        coshf = (ecc + cnu) / (ecc * cnu + 1.0)
        x = primal(coshf) > 1.0 ? coshf : one(coshf)
        ea = log(x + x * sqrt(1.0 - 1.0 / x / x))
        m0 = _m_dsign(ecc * sinh(ea) - ea, nu)
    else
        ea = tan(nu / 2.0)
        m0 = _m_dsign(ea + ea * (ea * ea) / 3.0, nu)
    end
    return SVector(promote(rp, ecc, inc, lnode, argp, m0, _m_f(et), mu)...)
end

"""
    oscelt(state, et, mu) -> SVector{8}

Return the osculating conic elements `(rp, ecc, inc, lnode, argp, m0, t0, mu)` of the
state `state` at the epoch `et` about a mass with the gravitational parameter `mu`, as
`oscelt_c`. The angles are in radians: `inc` in [0, pi], `lnode`, `argp` and (for an
ellipse) `m0` in [0, 2 pi). An eccentricity within 1e-10 of 1 is set to 1, and an
inclination within 1e-10 of 0 or pi to that value, as in CSPICE. Requires `mu > 0`
(`SpiceError(:nonpositivemass)`) and a position and a velocity that are not zero and
not parallel (`SpiceError(:degeneratecase)`).
"""
oscelt(state::AbstractVector, et::Real, mu::Real) = _m_oscelt(_m_v6(state), et, _m_f(mu))

"""
    oscltx(state, et, mu) -> SVector{20}

Return the elements of `oscelt` in places 1 to 8, the true anomaly in [0, 2 pi) in
place 9, the semi-major axis in place 10 (negative for a hyperbola, 0 where it is not
defined) and the period in place 11 (0 where it is not defined or would overflow), as
`oscltx_c`; places 12 to 20 are zero. Errors as `oscelt`.
"""
function oscltx(state::AbstractVector, et::Real, mu::Real)
    s = _m_v6(state)
    mu = _m_f(mu)
    el = _m_oscelt(s, et, mu)
    T = eltype(el)
    limit = floatmax(Float64) / 200.0
    pos, vel = _m_pv(s)
    rhat, rmag = _m_unorm(pos)
    vmag = _m_norm(vel)
    ecc0 = abs(primal(el[2])) <= 1.0e-10 ? 0.0 : primal(el[2])
    nu = zero(T)
    if ecc0 != 0.0
        c1 = 1.0 / mu * (rmag * (vmag * vmag) - mu)
        c2 = -(1.0 / mu) * _m_dot(pos, vel)
        eccvec = c1 * rhat + c2 * vel
        pole = _m_ucrss(pos, vel)
        # `twovec_(eccvec, 1, pole, 3)`: x along the eccentricity vector, y along
        # pole × x, z along x × y
        ax2 = _m_ucrss(pole, eccvec)
        _m_iszero(ax2) && throw(SpiceError(:dependentvectors,
            "the eccentricity vector and the orbit pole are linearly dependent"))
        xmat = _m_rows(_m_hat(eccvec), ax2, _m_ucrss(eccvec, ax2))
        rperi = _m_mxv(xmat, pos)
        _, lon, _ = _m_reclat(rperi, false)
        nu = (primal(lon) < 0.0 ? lon + _M_TWOPI : lon) + zero(T)
    else
        nu = el[6]
    end
    a = zero(T)
    tau = zero(T)
    if primal(el[2]) != 1.0 && primal(rmag) > primal(mu) / limit
        e = vmag * vmag * 0.5 - mu / rmag
        if abs(primal(e)) >= abs(primal(mu)) / limit
            a = -mu / (e * 2.0) + zero(T)
            if primal(el[2]) < 1.0
                b = (limit / _M_TWOPI)^0.66666666666666663
                mucubr = mu^0.33333333333333331
                cmptau = primal(mu) >= 1.0 ? primal(a / mucubr) < b : primal(a) < b * primal(mucubr)
                cmptau && (tau = _M_TWOPI * (a / mucubr)^1.5 + zero(T))
            end
        end
    end
    z = zero(T)
    return SVector{20,T}(el[1], el[2], el[3], el[4], el[5], el[6], el[7], el[8], nu, a, tau,
        z, z, z, z, z, z, z, z, z)
end
