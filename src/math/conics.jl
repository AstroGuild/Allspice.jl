# Two-body (conic) calls of CSPICE: Stumpff functions, propagation of a state, the
# conversions between a state and osculating elements, and the propagation of
# equinoctial elements with secular rates.
#
# Elements are in the CSPICE order: perifocal distance `rp`, eccentricity `ecc`,
# inclination `inc`, longitude of the ascending node `lnode`, argument of periapsis
# `argp`, mean anomaly at epoch `m0`, epoch `t0` and gravitational parameter `mu`.
# The propagator is the private `_prop2b` of `spk_types/conic.jl`, which SPK types 5
# and 15 use, and the state from equinoctial elements of `eqncpv` is the private
# `_eqn_state` of SPK type 17; nothing of them is repeated here.

"""
    stmp03(x) -> (c0, c1, c2, c3)

Return the Stumpff functions `c0`, `c1`, `c2`, `c3` of `x`, as `stmp03` of
SPICELIB (CSPICE has no `_c` wrapper).

`c0 = 1 - x * c2` and `c1 = 1 - x * c3` hold to rounding. For `x > 0`, `c0` is
`cos(sqrt(x))` and `c1` is `sin(sqrt(x)) / sqrt(x)`.

Throws `SpiceError(:valueoutofrange)` if `x <= -(log(2) + log(floatmax(Float64)))^2`,
about -5.05e5, where `c0` would overflow.

# Examples
```jldoctest
julia> stmp03(0)
(1.0, 1.0, 0.5, 0.16666666666666666)
```

See also [`prop2b`](@ref).
"""
function stmp03(x::Real)
    xx = _m_f(x)
    primal(xx) <= _STUMPFF_LBOUND && throw(SpiceError(:valueoutofrange,
        "the argument $(primal(xx)) of the Stumpff functions is below the lower bound $_STUMPFF_LBOUND"))
    return _stumpff(xx)
end

"""
    prop2b(gm, pvinit, dt) -> SVector{6}

Return the state of a body at the time `dt` after the state `pvinit`, for two-body
motion about a mass with the gravitational parameter `gm`, as `prop2b_c` of CSPICE.

The units must be consistent, for example km, km/s, seconds and km³/s². `gm` and
`pvinit` are converted to `Float64`. `dt` can be a `ForwardDiff.Dual` or
`Allspice.Taylor` number; the result then holds the derivatives with respect to `dt`.

Throws:
- `SpiceError(:nonpositivemass)` if `gm <= 0`.
- `SpiceError(:zeroposition)` or `SpiceError(:zerovelocity)` if that part of `pvinit`
  is zero.
- `SpiceError(:nonconicmotion)` if the position and the velocity are parallel.
- `SpiceError(:dtoutofrange)` if `dt` is so large that the solution would overflow.
  The bound depends on the orbit, as in CSPICE.
- `SpiceError(:badarg)` if `pvinit` does not have 6 elements.

# Examples
A quarter of a circular orbit with unit radius and `gm = 1`:
```jldoctest
julia> prop2b(1, [1, 0, 0, 0, 1, 0], pi / 2) ≈ [0, 1, 0, -1, 0, 0]
true
```

See also [`conics`](@ref), [`oscelt`](@ref).
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
`elts`, as `conics_c` of CSPICE.

| Index | Name    | Meaning                                    |
|:------|:--------|:-------------------------------------------|
| 1     | `rp`    | perifocal distance                         |
| 2     | `ecc`   | eccentricity                               |
| 3     | `inc`   | inclination, rad                           |
| 4     | `lnode` | longitude of the ascending node, rad       |
| 5     | `argp`  | argument of periapsis, rad                 |
| 6     | `m0`    | mean anomaly at the epoch `t0`, rad        |
| 7     | `t0`    | epoch of the elements                      |
| 8     | `mu`    | gravitational parameter of the central mass |

The units must be consistent, for example km, seconds and km³/s². The elements are
converted to `Float64`. `et` can be a `ForwardDiff.Dual` or `Allspice.Taylor` number.

Throws:
- `SpiceError(:badeccentricity)` if `ecc < 0`.
- `SpiceError(:badperiapsevalue)` if `rp <= 0`.
- `SpiceError(:badgm)` if `mu <= 0`.
- `SpiceError(:dtoutofrange)` if `et - t0` is so large that the solution would
  overflow. The bound depends on the orbit, as in CSPICE.
- `SpiceError(:badarg)` if `elts` does not have 8 elements.

# Examples
A circular orbit with unit radius and `mu = 1`, a quarter period after the epoch:
```jldoctest
julia> conics([1, 0, 0, 0, 0, 0, 0, 1], pi / 2) ≈ [0, 1, 0, -1, 0, 0]
true
```

See also [`oscelt`](@ref), [`prop2b`](@ref).
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
        inc = Float64(pi) - zero(inc)
        n = SVector(one(z0), z0, z0)
    end
    lnode = atan(n[2], n[1])
    primal(lnode) < 0.0 && (lnode += _M_TWOPI)
    argp = zero(inc)
    if primal(ecc) != 0.0
        argp = vsep(n, e) - zero(inc)
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

Return the osculating conic elements `(rp, ecc, inc, lnode, argp, m0, t0, mu)` of a
state, as `oscelt_c` of CSPICE.

`state` is the position and the velocity at the epoch `et`, relative to a mass with the
gravitational parameter `mu`. The table in [`conics`](@ref) defines the elements; `t0`
is `et`. The angles are in radians: `inc` in [0, pi], and `lnode` and `argp` in
[0, 2 pi). `m0` is in [0, 2 pi) for an ellipse.

As in CSPICE, an eccentricity within 1e-10 of 1 becomes 1. An inclination within 1e-10
of 0 or pi becomes that value, and `lnode` is then 0.

Throws `SpiceError(:nonpositivemass)` if `mu <= 0`. Throws
`SpiceError(:degeneratecase)` if the position or the velocity is zero or if they are
parallel. Throws `SpiceError(:badarg)` if `state` does not have 6 elements.

# Examples
The apoapsis of an orbit with `mu = 1`, semi-major axis 4/3 and eccentricity 0.5:
```jldoctest
julia> round.(Tuple(oscelt([0, 2, 0, -0.5, 0, 0], 0, 1)); digits=6)
(0.666667, 0.5, 0.0, 0.0, 4.712389, 3.141593, 0.0, 1.0)
```

See also [`oscltx`](@ref), [`conics`](@ref).
"""
oscelt(state::AbstractVector, et::Real, mu::Real) = _m_oscelt(_m_v6(state), et, _m_f(mu))

"""
    oscltx(state, et, mu) -> SVector{20}

Return the osculating elements of [`oscelt`](@ref) with the true anomaly, the
semi-major axis and the period, as `oscltx_c` of CSPICE.

| Index | Content                                                                    |
|:------|:---------------------------------------------------------------------------|
| 1:8   | the elements of [`oscelt`](@ref)                                           |
| 9     | true anomaly, rad, in [0, 2 pi)                                            |
| 10    | semi-major axis; negative for a hyperbola, 0 where it is not defined       |
| 11    | period; 0 for an orbit that is not an ellipse, or if it would overflow     |
| 12:20 | zero                                                                       |

Throws as [`oscelt`](@ref).

# Examples
```jldoctest
julia> elts = oscltx([0, 2, 0, -0.5, 0, 0], 0, 1);

julia> round.(Tuple(elts[9:11]); digits=6)
(3.141593, 1.333333, 9.673597)
```

See also [`oscelt`](@ref), [`conics`](@ref).
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
        nu = (primal(lon) < 0.0 ? lon + _M_TWOPI : lon) - zero(T)
    else
        nu = el[6]
    end
    a = zero(T)
    tau = zero(T)
    if primal(el[2]) != 1.0 && primal(rmag) > primal(mu) / limit
        e = vmag * vmag * 0.5 - mu / rmag
        if abs(primal(e)) >= abs(primal(mu)) / limit
            a = -mu / (e * 2.0) - zero(T)
            if primal(el[2]) < 1.0
                b = (limit / _M_TWOPI)^0.66666666666666663
                mucubr = mu^0.33333333333333331
                cmptau = primal(mu) >= 1.0 ? primal(a / mucubr) < b : primal(a) < b * primal(mucubr)
                cmptau && (tau = _M_TWOPI * (a / mucubr)^1.5 - zero(T))
            end
        end
    end
    z = zero(T)
    return SVector{20,T}(el[1], el[2], el[3], el[4], el[5], el[6], el[7], el[8], nu, a, tau,
        z, z, z, z, z, z, z, z, z)
end

"""
    eqncpv(et, epoch, eqel, rapol, decpol) -> SVector{6}

Return the state at `et` of a body whose orbit is given by equinoctial elements with
the rates of the node and of the periapsis, as `eqncpv_c` of CSPICE.

`eqel` holds the elements at `epoch` (seconds past J2000 TDB):

| Index | Meaning                                                  |
|:------|:---------------------------------------------------------|
| 1     | semi-major axis `a`, km                                  |
| 2, 3  | `h = e sin(argp + node)`, `k = e cos(argp + node)`       |
| 4     | mean longitude `argp + node + M0`, rad                   |
| 5, 6  | `p = tan(inc/2) sin(node)`, `q = tan(inc/2) cos(node)`   |
| 7     | rate of the longitude of periapsis, rad/s                |
| 8     | mean longitude rate, rad/s                               |
| 9     | rate of the longitude of the ascending node, rad/s       |

The elements are relative to the equator of a pole with the right ascension `rapol`
and the declination `decpol` (rad), given in the frame of the result. This is the
model of SPK type 17. The result is in km and km/s. The arguments other than `et` are
converted to `Float64`; `et` can be a `ForwardDiff.Dual` or `Allspice.Taylor` number.
A zero mean longitude rate gives a velocity of NaN, as in CSPICE.

Throws `SpiceError`:
- `:badsemiaxis` if `a <= 0`;
- `:eccoutofrange` if `sqrt(h² + k²) > 0.9`;
- `:eccoutofbounds` if `h² + k² >= 0.81` at `et` (an eccentricity of 0.9 that the
  first check passes);
- `:badarg` if `eqel` does not have 9 elements.

# Examples
A circular equatorial orbit (pole on the z axis) a quarter period after the epoch:
```jldoctest
julia> n = 1e-3; s = eqncpv(pi / 2 / n, 0, [7000, 0, 0, 0, 0, 0, 0, n, 0], -pi / 2, pi / 2);

julia> s ≈ [0, 7000, 0, -7, 0, 0]
true
```

See also [`conics`](@ref).
"""
function eqncpv(et::Real, epoch::Real, eqel::AbstractVector, rapol::Real, decpol::Real)
    length(eqel) == 9 || _m_baddim("a set of equinoctial elements", 9)
    a, h, k, l, p, q, dlpdt, mldt, nodedt = ntuple(i -> Float64(eqel[begin+i-1]), Val(9))
    a > 0.0 || throw(SpiceError(:badsemiaxis, "the semi-major axis is $a; it must be positive"))
    ecc = sqrt(h * h + k * k)
    ecc > 0.9 && throw(SpiceError(:eccoutofrange,
        "the eccentricity of the elements is $ecc; it must not be above 0.9"))
    pos, vel = _eqn_state(_m_f(et), Float64(epoch), a, h, k, l, p, q, dlpdt, mldt, nodedt,
        Float64(rapol), Float64(decpol))
    return _m_state(pos, vel)
end
