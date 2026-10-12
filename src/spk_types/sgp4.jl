# SGP4 and SDP4 propagation of two-line elements, as the routines `zzsgp4.f`
# (`xxsgp4i`, `xxsgp4e`), `zzinil.f`, `zzdscm.f`, `zzdsin.f`, `zzdspc.f` and `zzdspr.f` of
# CSPICE N0067, which follow the code of Vallado, Crawford, Hujsak and Kelso,
# "Revisiting Spacetrack Report #3".
#
# CSPICE keeps the initialized propagator in the saved variables of `zzsgp4`. Here
# there is no saved state: `_sgp4_init` returns the constants of one element set as a
# plain-bits record (all `Float64`), and `_sgp4_eval` propagates that record to a time
# that can have any number type, so a `Taylor` time gives the time derivatives of the
# state. The operations have the order of the CSPICE source, so a plain time gives the
# CSPICE state to rounding level. Every discrete decision uses `primal`.

@noinline _sgp4_error(kind::Symbol, msg::String) = throw(SpiceError(kind, msg))

"The value of `twopi_()` of CSPICE."
const _SGP4_TWOPI = 2.0 * Float64(pi)
"The constant 2/3 as the source writes it."
const _SGP4_X2O3 = 0.66666666666666663
"Earth rotation rate in radians per minute (`RPTIM`)."
const _SGP4_RPTIM = 0.00437526908801129966
"Operation mode of `zzinil`: sidereal time of Spacetrack Report #3 (the mode of `spke10` and `evsgp4`)."
const _SGP4_AFSPC = 1
"Operation mode of `zzinil`: sidereal time of the IAU 1982 model."
const _SGP4_IMPRVD = 2

"""
    _sgp4_pow(x, p::Float64)

`x` to the power `p` (`pow_dd` of f2c). For a `Taylor` value the coefficients come
from the recurrence `r_k = (1/(k x_0)) Σ_{j=1}^{k} (p j - (k - j)) x_j r_{k-j}`.
Requires `primal(x) > 0` (not checked).
"""
@inline _sgp4_pow(x::Real, p::Float64) = x^p
@generated function _sgp4_pow(x::Taylor{N,T}, p::Float64) where {N,T<:Real}
    ex = _unpack(:a, :x, N)
    push!(ex, :(r0 = _sgp4_pow(a0, p)))
    for k in 1:N
        s = _sum([:((p * $j - $(k - j)) * $(_sym(:a, j)) * $(_sym(:r, k - j))) for j in 1:k])
        push!(ex, :($(_sym(:r, k)) = $s / ($k * a0)))
    end
    return Expr(:block, Expr(:meta, :inline), ex..., _pack(:r, N))
end

"""
    _sgp4_mod(x::T, m::Float64) -> T

Remainder of `x` divided by `m` with the sign of `x` (`d_mod` of f2c). The plain value
is `rem(primal(x), m)`; the derivative parts of `x` do not change.
"""
@inline _sgp4_mod(x::Real, m::Float64) = _with_primal(x, rem(primal(x), m))

"""
    _sgp4_lsterms(zcosg, zsing, zcosi, zsini, zcosh, zsinh, cc, c) -> NamedTuple

One pass of the loop of `zzdscm` over the solar and the lunar terms: the values
`s1 ... s7`, `z1 ... z3`, `z11 ... z33` for the perturbing body with the given
orientation constants and the coefficient `cc`. `c` holds the values of the orbit
(`cosim`, `sinim`, `cosomm`, `sinomm`, `emsq`, `betasq`, `rtemsq`, `eccm`, `xnoi`).
"""
@inline function _sgp4_lsterms(zcosg::Float64, zsing::Float64, zcosi::Float64,
        zsini::Float64, zcosh::Float64, zsinh::Float64, cc::Float64, c)
    cosim, sinim, cosomm, sinomm = c.cosim, c.sinim, c.cosomm, c.sinomm
    emsq, betasq, rtemsq = c.emsq, c.betasq, c.rtemsq
    a1 = zcosg * zcosh + zsing * zcosi * zsinh
    a3 = -zsing * zcosh + zcosg * zcosi * zsinh
    a7 = -zcosg * zsinh + zsing * zcosi * zcosh
    a8 = zsing * zsini
    a9 = zsing * zsinh + zcosg * zcosi * zcosh
    a10 = zcosg * zsini
    a2 = cosim * a7 + sinim * a8
    a4 = cosim * a9 + sinim * a10
    a5 = -sinim * a7 + cosim * a8
    a6 = -sinim * a9 + cosim * a10
    x1 = a1 * cosomm + a2 * sinomm
    x2 = a3 * cosomm + a4 * sinomm
    x3 = -a1 * sinomm + a2 * cosomm
    x4 = -a3 * sinomm + a4 * cosomm
    x5 = a5 * sinomm
    x6 = a6 * sinomm
    x7 = a5 * cosomm
    x8 = a6 * cosomm
    z31 = x1 * 12.0 * x1 - x3 * 3.0 * x3
    z32 = x1 * 24.0 * x2 - x3 * 6.0 * x4
    z33 = x2 * 12.0 * x2 - x4 * 3.0 * x4
    z1 = (a1 * a1 + a2 * a2) * 3.0 + z31 * emsq
    z2 = (a1 * a3 + a2 * a4) * 6.0 + z32 * emsq
    z3 = (a3 * a3 + a4 * a4) * 3.0 + z33 * emsq
    z11 = a1 * -6.0 * a5 + emsq * (x1 * -24.0 * x7 - x3 * 6.0 * x5)
    z12 = (a1 * a6 + a3 * a5) * -6.0 +
          emsq * ((x2 * x7 + x1 * x8) * -24.0 - (x3 * x6 + x4 * x5) * 6.0)
    z13 = a3 * -6.0 * a6 + emsq * (x2 * -24.0 * x8 - x4 * 6.0 * x6)
    z21 = a2 * 6.0 * a5 + emsq * (x1 * 24.0 * x5 - x3 * 6.0 * x7)
    z22 = (a4 * a5 + a2 * a6) * 6.0 +
          emsq * ((x2 * x5 + x1 * x6) * 24.0 - (x4 * x7 + x3 * x8) * 6.0)
    z23 = a4 * 6.0 * a6 + emsq * (x2 * 24.0 * x6 - x4 * 6.0 * x8)
    z1 = z1 + z1 + betasq * z31
    z2 = z2 + z2 + betasq * z32
    z3 = z3 + z3 + betasq * z33
    s3 = cc * c.xnoi
    s2 = s3 * -0.5 / rtemsq
    s4 = s3 * rtemsq
    s1 = c.eccm * -15.0 * s4
    s5 = x1 * x3 + x2 * x4
    s6 = x2 * x3 + x1 * x4
    s7 = x2 * x4 - x1 * x3
    return (; s1, s2, s3, s4, s5, s6, s7, z1, z2, z3, z11, z12, z13, z21, z22, z23,
        z31, z32, z33)
end

"""
    _sgp4_dscom(epoch, ecco, argpo, inclo, nodeo, no) -> NamedTuple

Constants of the lunar and solar perturbations of a deep-space orbit, as `zzdscm` at
the element epoch (`TC = 0`): the solar terms `sun`, the lunar terms `moon` (each from
`_sgp4_lsterms`), the coefficients of the long-period periodics (`se2` ... `xh3`), the
mean anomalies `zmos` and `zmol` of the Sun and the Moon, and `sinim`, `cosim`,
`emsq` of the orbit.

Contract: `epoch` is in days from 1949 December 31 00:00 UTC, the angles are in
radians, `no` is the mean motion in radians per minute after the removal of the
Kozai correction, with `no != 0` and `0 <= ecco < 1`.
"""
function _sgp4_dscom(epoch::Float64, ecco::Float64, argpo::Float64, inclo::Float64,
        nodeo::Float64, no::Float64)
    zes = 0.01675
    zel = 0.0549
    c1ss = 2.9864797e-6
    c1l = 4.7968065e-7
    zsinis = 0.39785416
    zcosis = 0.91744867
    zcosgs = 0.1945905
    zsings = -0.98088458
    eccm = ecco
    snodm, cnodm = sincos(nodeo)
    sinomm, cosomm = sincos(argpo)
    sinim, cosim = sincos(inclo)
    emsq = eccm * eccm
    betasq = 1.0 - emsq
    rtemsq = sqrt(betasq)
    # EPOCH + 18261.5 is the count of days from 1899 December 31 12:00.
    day = epoch + 18261.5 + 0.0 / 1440.0
    xnodce = rem(4.523602 - day * 9.2422029e-4, _SGP4_TWOPI)
    stem, ctem = sincos(xnodce)
    zcosil = 0.91375164 - ctem * 0.03568096
    zsinil = sqrt(1.0 - zcosil * zcosil)
    zsinhl = stem * 0.089683511 / zsinil
    zcoshl = sqrt(1.0 - zsinhl * zsinhl)
    gam = day * 0.001944368 + 5.8351514
    zx = stem * 0.39785416 / zsinil
    zy = zcoshl * ctem + zsinhl * 0.91744867 * stem
    zx = atan(zx, zy)
    zx = gam + zx - xnodce
    zsingl, zcosgl = sincos(zx)
    c = (; cosim, sinim, cosomm, sinomm, emsq, betasq, rtemsq, eccm, xnoi=1.0 / no)
    sun = _sgp4_lsterms(zcosgs, zsings, zcosis, zsinis, cnodm, snodm, c1ss, c)
    moon = _sgp4_lsterms(zcosgl, zsingl, zcosil, zsinil, zcoshl * cnodm + zsinhl * snodm,
        snodm * zcoshl - cnodm * zsinhl, c1l, c)
    zmol = rem(day * 0.2299715 + 4.7199672 - gam, _SGP4_TWOPI)
    zmos = rem(day * 0.017201977 + 6.2565837, _SGP4_TWOPI)
    # Solar terms.
    se2 = sun.s1 * 2.0 * sun.s6
    se3 = sun.s1 * 2.0 * sun.s7
    si2 = sun.s2 * 2.0 * sun.z12
    si3 = sun.s2 * 2.0 * (sun.z13 - sun.z11)
    sl2 = sun.s3 * -2.0 * sun.z2
    sl3 = sun.s3 * -2.0 * (sun.z3 - sun.z1)
    sl4 = sun.s3 * -2.0 * (-21.0 - emsq * 9.0) * zes
    sgh2 = sun.s4 * 2.0 * sun.z32
    sgh3 = sun.s4 * 2.0 * (sun.z33 - sun.z31)
    sgh4 = sun.s4 * -18.0 * zes
    sh2 = sun.s2 * -2.0 * sun.z22
    sh3 = sun.s2 * -2.0 * (sun.z23 - sun.z21)
    # Lunar terms.
    ee2 = moon.s1 * 2.0 * moon.s6
    e3 = moon.s1 * 2.0 * moon.s7
    xi2 = moon.s2 * 2.0 * moon.z12
    xi3 = moon.s2 * 2.0 * (moon.z13 - moon.z11)
    xl2 = moon.s3 * -2.0 * moon.z2
    xl3 = moon.s3 * -2.0 * (moon.z3 - moon.z1)
    xl4 = moon.s3 * -2.0 * (-21.0 - emsq * 9.0) * zel
    xgh2 = moon.s4 * 2.0 * moon.z32
    xgh3 = moon.s4 * 2.0 * (moon.z33 - moon.z31)
    xgh4 = moon.s4 * -18.0 * zel
    xh2 = moon.s2 * -2.0 * moon.z22
    xh3 = moon.s2 * -2.0 * (moon.z23 - moon.z21)
    return (; sun, moon, sinim, cosim, emsq, zmol, zmos, se2, se3, si2, si3, sl2, sl3,
        sl4, sgh2, sgh3, sgh4, sh2, sh3, ee2, e3, xi2, xi3, xl2, xl3, xl4, xgh2, xgh3,
        xgh4, xh2, xh3)
end

"""
    _sgp4_dsinit(xke, c, argpo, gsto, mo, mdot, no, nodeo, nodedot, xpidot, ecco, eccsq, inclo) -> NamedTuple

Secular rates of the lunar and solar perturbations (`dedt`, `didt`, `dmdt`, `dnodt`,
`domdt`) and the constants of the resonance of a deep-space orbit, as `zzdsin` at the
element epoch (`T = TC = 0`). `irez` is 0 for no resonance, 1 for a one-day orbit
(with `del1`, `del2`, `del3`) and 2 for a half-day orbit of eccentricity 0.5 or more
(with `d2201` ... `d5433`); `xlamo` and `xfact` start the resonance integration.

Contract: `c` is the result of `_sgp4_dscom` for the same elements; `xke` is the
fourth geophysical constant.
"""
function _sgp4_dsinit(xke::Float64, c, argpo::Float64, gsto::Float64, mo::Float64,
        mdot::Float64, no::Float64, nodeo::Float64, nodedot::Float64, xpidot::Float64,
        ecco::Float64, eccsq::Float64, inclo::Float64)
    q22 = 1.7891679e-6
    q31 = 2.1460748e-6
    q33 = 2.2123015e-7
    root22 = 1.7891679e-6
    root44 = 7.3636953e-9
    root54 = 2.1765803e-9
    rptim = _SGP4_RPTIM
    root32 = 3.7393792e-7
    root52 = 1.1428639e-7
    znl = 1.5835218e-4
    zns = 1.19459e-5
    sun, moon = c.sun, c.moon
    cosim, sinim = c.cosim, c.sinim
    emsq = c.emsq
    xn = no
    eccm = ecco
    inclm = inclo
    irez = 0
    (xn < 0.0052359877 && xn > 0.0034906585) && (irez = 1)
    (xn >= 0.00826 && xn <= 0.00924 && eccm >= 0.5) && (irez = 2)
    # Solar terms.
    ses = sun.s1 * zns * sun.s5
    sis = sun.s2 * zns * (sun.z11 + sun.z13)
    sls = -zns * sun.s3 * (sun.z1 + sun.z3 - 14.0 - emsq * 6.0)
    sghs = sun.s4 * zns * (sun.z31 + sun.z33 - 6.0)
    shs = -zns * sun.s2 * (sun.z21 + sun.z23)
    polar = inclm < 0.052359877 || inclm > Float64(pi) - 0.052359877
    polar && (shs = 0.0)
    sinim != 0.0 && (shs /= sinim)
    sgs = sghs - cosim * shs
    # Lunar terms.
    dedt = ses + moon.s1 * znl * moon.s5
    didt = sis + moon.s2 * znl * (moon.z11 + moon.z13)
    dmdt = sls - znl * moon.s3 * (moon.z1 + moon.z3 - 14.0 - emsq * 6.0)
    sghl = moon.s4 * znl * (moon.z31 + moon.z33 - 6.0)
    shl = -znl * moon.s2 * (moon.z21 + moon.z23)
    polar && (shl = 0.0)
    domdt = sgs + sghl
    dnodt = shs
    if sinim != 0.0
        domdt -= cosim / sinim * shl
        dnodt += shl / sinim
    end
    # Resonance. With T = 0 the mean elements of the epoch do not change.
    theta = rem(gsto + 0.0 * rptim, _SGP4_TWOPI)
    eccm += dedt * 0.0
    emsq = eccm * eccm
    d2201 = d2211 = d3210 = d3222 = d4410 = d4422 = 0.0
    d5220 = d5232 = d5421 = d5433 = 0.0
    del1 = del2 = del3 = xfact = xlamo = 0.0
    if irez != 0
        aonv = (xn / xke)^_SGP4_X2O3
        if irez == 2
            # Geopotential resonance of the half-day orbits.
            cosisq = cosim * cosim
            eccm = ecco
            emsq = eccsq
            eoc = eccm * emsq
            g201 = -0.306 - (eccm - 0.64) * 0.44
            if eccm <= 0.65
                g211 = 3.616 - eccm * 13.247 + emsq * 16.29
                g310 = eccm * 117.39 - 19.302 - emsq * 228.419 + eoc * 156.591
                g322 = eccm * 109.7927 - 18.9068 - emsq * 214.6334 + eoc * 146.5816
                g410 = eccm * 242.694 - 41.122 - emsq * 471.094 + eoc * 313.953
                g422 = eccm * 841.88 - 146.407 - emsq * 1629.014 + eoc * 1083.435
                g520 = eccm * 3017.977 - 532.114 - emsq * 5740.032 + eoc * 3708.276
            else
                g211 = eccm * 331.819 - 72.099 - emsq * 508.738 + eoc * 266.724
                g310 = eccm * 1582.851 - 346.844 - emsq * 2415.925 + eoc * 1246.113
                g322 = eccm * 1554.908 - 342.585 - emsq * 2366.899 + eoc * 1215.972
                g410 = eccm * 4758.686 - 1052.797 - emsq * 7193.992 + eoc * 3651.957
                g422 = eccm * 16178.11 - 3581.69 - emsq * 24462.77 + eoc * 12422.52
                if eccm > 0.715
                    g520 = eccm * 29936.92 - 5149.66 - emsq * 54087.36 + eoc * 31324.56
                else
                    g520 = 1464.74 - eccm * 4664.75 + emsq * 3763.64
                end
            end
            if eccm < 0.7
                g533 = eccm * 4988.61 - 919.2277 - emsq * 9064.77 + eoc * 5542.21
                g521 = eccm * 4568.6173 - 822.71072 - emsq * 8491.4146 + eoc * 5337.524
                g532 = eccm * 4690.25 - 853.666 - emsq * 8624.77 + eoc * 5341.4
            else
                g533 = eccm * 161616.52 - 37995.78 - emsq * 229838.2 + eoc * 109377.94
                g521 = eccm * 218913.95 - 51752.104 - emsq * 309468.16 + eoc * 146349.42
                g532 = eccm * 170470.89 - 40023.88 - emsq * 242699.48 + eoc * 115605.82
            end
            sini2 = sinim * sinim
            f220 = (cosim * 2.0 + 1.0 + cosisq) * 0.75
            f221 = sini2 * 1.5
            f321 = sinim * 1.875 * (1.0 - cosim * 2.0 - cosisq * 3.0)
            f322 = sinim * -1.875 * (cosim * 2.0 + 1.0 - cosisq * 3.0)
            f441 = sini2 * 35.0 * f220
            f442 = sini2 * 39.375 * sini2
            f522 = sinim * 9.84375 * (sini2 * (1.0 - cosim * 2.0 - cosisq * 5.0) +
                                      (cosim * 4.0 - 2.0 + cosisq * 6.0) * 0.33333333)
            f523 = sinim * (sini2 * 4.92187512 * (-2.0 - cosim * 4.0 + cosisq * 10.0) +
                            (cosim * 2.0 + 1.0 - cosisq * 3.0) * 6.56250012)
            f542 = sinim * 29.53125 *
                   (2.0 - cosim * 8.0 + cosisq * (cosim * 8.0 - 12.0 + cosisq * 10.0))
            f543 = sinim * 29.53125 *
                   (-2.0 - cosim * 8.0 + cosisq * (cosim * 8.0 + 12.0 - cosisq * 10.0))
            xno2 = xn * xn
            ainv2 = aonv * aonv
            temp1 = xno2 * 3.0 * ainv2
            temp = temp1 * root22
            d2201 = temp * f220 * g201
            d2211 = temp * f221 * g211
            temp1 *= aonv
            temp = temp1 * root32
            d3210 = temp * f321 * g310
            d3222 = temp * f322 * g322
            temp1 *= aonv
            temp = temp1 * 2.0 * root44
            d4410 = temp * f441 * g410
            d4422 = temp * f442 * g422
            temp1 *= aonv
            temp = temp1 * root52
            d5220 = temp * f522 * g520
            d5232 = temp * f523 * g532
            temp = temp1 * 2.0 * root54
            d5421 = temp * f542 * g521
            d5433 = temp * f543 * g533
            xlamo = rem(mo + nodeo + nodeo - theta - theta, _SGP4_TWOPI)
            xfact = mdot + dmdt + (nodedot + dnodt - rptim) * 2.0 - no
        else
            # Resonance of the one-day orbits.
            g200 = emsq * (emsq * 0.8125 - 2.5) + 1.0
            g310 = emsq * 2.0 + 1.0
            g300 = emsq * (emsq * 6.60937 - 6.0) + 1.0
            f220 = (cosim + 1.0) * 0.75 * (cosim + 1.0)
            f311 = sinim * 0.9375 * sinim * (cosim * 3.0 + 1.0) - (cosim + 1.0) * 0.75
            f330 = cosim + 1.0
            f330 = f330 * 1.875 * f330 * f330
            del1 = xn * 3.0 * xn * aonv * aonv
            del2 = del1 * 2.0 * f220 * g200 * q22
            del3 = del1 * 3.0 * f330 * g300 * q33 * aonv
            del1 = del1 * f311 * g310 * q31 * aonv
            xlamo = rem(mo + nodeo + argpo - theta, _SGP4_TWOPI)
            xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no
        end
    end
    return (; irez, dedt, didt, dmdt, dnodt, domdt, d2201, d2211, d3210, d3222, d4410,
        d4422, d5220, d5232, d5421, d5433, del1, del2, del3, xfact, xlamo)
end

"""
    _sgp4_gsto(epoch::Float64, opmode::Int) -> Float64

Greenwich sidereal angle in `[0, 2π)` at `epoch` (days from 1949 December 31 00:00
UTC), as the end of `zzinil`: the formula of Spacetrack Report #3 for
`opmode == _SGP4_AFSPC` and the IAU 1982 formula for `_SGP4_IMPRVD`. Another mode
throws `SpiceError(:unknownmode)`.
"""
function _sgp4_gsto(epoch::Float64, opmode::Int)
    gsto = if opmode == _SGP4_IMPRVD
        radday = _SGP4_TWOPI * 1.002737909350795
        temp = epoch + 2433281.5
        tut1 = (trunc(temp - 0.5) + 0.5 - 2451545.0) / 36525.0
        tut1 * 628.331970688841 + 1.75336855923327 + tut1 * 6.77071394490334e-6 * tut1 -
        tut1 * 4.50876723431868e-10 * tut1 * tut1 + radday * (temp - 0.5 - trunc(temp - 0.5))
    elseif opmode == _SGP4_AFSPC
        ts70 = epoch - 7305.0
        ids70 = Float64(unsafe_trunc(Int, clamp(ts70 + 1e-8, -1.0e15, 1.0e15)))
        tfrac = ts70 - ids70
        c1 = 0.0172027916940703639
        thgr70 = 1.7321343856509374
        fk5r = 5.07551419432269442e-15
        c1p2p = c1 + _SGP4_TWOPI
        thgr70 + c1 * ids70 + c1p2p * tfrac + ts70 * ts70 * fk5r
    else
        _sgp4_error(:unknownmode, "the SGP4 operation mode $opmode is not 1 or 2")
    end
    gsto = rem(gsto, _SGP4_TWOPI)
    return gsto < 0.0 ? gsto + _SGP4_TWOPI : gsto
end

"""
    _sgp4_init(geophs::NTuple{8,Float64}, elems::NTuple{10,Float64}, epoch::Float64, opmode::Int) -> NamedTuple

Constants of the propagation of one element set, as `xxsgp4i`. `geophs` holds J2, J3,
J4, KE, QO, SO, ER, AE and `elems` holds NDT20, NDD60, BSTAR, INCL, NODE0, ECC,
OMEGA, M0, N0 and the epoch, in the units of `getelm`. `epoch` is the element epoch
in days from 1949 December 31 00:00 UTC (Julian date UTC minus 2433281.5), which
CSPICE computes with `ttrans` from `elems[10]`.

Contract: the result is plain bits (`isbits`), with all the values that `_sgp4_eval`
needs: `dodeep` tells whether the period is 225 minutes or more (SDP4), and `dosimp`
whether the drag terms of higher order are left out. Throws
`SpiceError(:suborbital)` if the perigee radius is below one Earth radius, and
`SpiceError(:invalidelements)` for an eccentricity above 1 with a negative mean
motion, where CSPICE leaves the constants undefined. As in CSPICE, QO, SO and AE
are not read: the atmosphere bounds 120 km and 78 km are fixed in the code.
"""
function _sgp4_init(geophs::NTuple{8,Float64}, elems::NTuple{10,Float64}, epoch::Float64,
        opmode::Int)
    bstar, inclo, nodeo, ecco = elems[3], elems[4], elems[5], elems[6]
    argpo, mo, no = elems[7], elems[8], elems[9]
    j2, j3, j4, xke, er = geophs[1], geophs[2], geophs[3], geophs[4], geophs[7]
    j3oj2 = j3 / j2
    x2o3 = _SGP4_X2O3
    ss = 78.0 / er + 1.0
    d = 42.0 / er
    d *= d
    qzms2t = d * d
    temp4 = 1.5e-12
    # ZZINIL: auxiliary quantities of the epoch and the mean motion without the Kozai
    # correction.
    eccsq = ecco * ecco
    omeosq = 1.0 - eccsq
    rteosq = sqrt(max(omeosq, 0.0))
    omeosq < 0.0 && (rteosq = NaN)
    cosio = cos(inclo)
    cosio2 = cosio * cosio
    ak = _sgp4_realpow(xke / no, x2o3)
    d1 = j2 * 0.75 * (cosio2 * 3.0 - 1.0) / (rteosq * omeosq)
    del = d1 / (ak * ak)
    adel = ak * (1.0 - del * del - del * (del * 134.0 * del / 81.0 + 0.33333333333333331))
    del = d1 / (adel * adel)
    no = no / (del + 1.0)
    ao = _sgp4_realpow(xke / no, x2o3)
    sinio = sin(inclo)
    po = ao * omeosq
    con42 = 1.0 - cosio2 * 5.0
    con41 = -con42 - cosio2 - cosio2
    posq = po * po
    rp = ao * (1.0 - ecco)
    gsto = _sgp4_gsto(epoch, opmode)
    rp < 1.0 && _sgp4_error(:suborbital, "TLE elements suborbital.")
    (omeosq >= 0.0 || no >= 0.0) || _sgp4_error(:invalidelements,
        "the eccentricity is above 1 and the mean motion is negative")

    dosimp = rp < 220.0 / er + 1.0
    sfour = ss
    qzms24 = qzms2t
    perige = (rp - 1.0) * er
    # For perigees below 156 km the constants S and QOMS2T change.
    if perige < 156.0
        sfour = perige - 78.0
        perige <= 98.0 && (sfour = 20.0)
        d = (120.0 - sfour) / er
        d *= d
        qzms24 = d * d
        sfour = sfour / er + 1.0
    end
    pinvsq = 1.0 / posq
    tsi = 1.0 / (ao - sfour)
    eta = ao * ecco * tsi
    etasq = eta * eta
    eeta = ecco * eta
    psisq = abs(1.0 - etasq)
    d = tsi
    d *= d
    coef = qzms24 * (d * d)
    coef1 = coef / _sgp4_realpow(psisq, 3.5)
    cc2 = coef1 * no * (ao * (etasq * 1.5 + 1.0 + eeta * (etasq + 4.0)) +
                        j2 * 0.375 * tsi / psisq * con41 * (etasq * 3.0 * (etasq + 8.0) + 8.0))
    cc1 = bstar * cc2
    cc3 = 0.0
    ecco > 1e-4 && (cc3 = coef * -2.0 * tsi * j3oj2 * no * sinio / ecco)
    x1mth2 = 1.0 - cosio2
    cc4 = no * 2.0 * coef1 * ao * omeosq *
          (eta * (etasq * 0.5 + 2.0) + ecco * (etasq * 2.0 + 0.5) -
           j2 * tsi / (ao * psisq) *
           (con41 * -3.0 * (1.0 - eeta * 2.0 + etasq * (1.5 - eeta * 0.5)) +
            x1mth2 * 0.75 * (etasq * 2.0 - eeta * (etasq + 1.0)) * cos(argpo * 2.0)))
    cc5 = coef1 * 2.0 * ao * omeosq * ((etasq + eeta) * 2.75 + 1.0 + eeta * etasq)
    cosio4 = cosio2 * cosio2
    temp1 = j2 * 1.5 * pinvsq * no
    temp2 = temp1 * 0.5 * j2 * pinvsq
    temp3 = j4 * -0.46875 * pinvsq * pinvsq * no
    mdot = no + temp1 * 0.5 * rteosq * con41 +
           temp2 * 0.0625 * rteosq * (13.0 - cosio2 * 78.0 + cosio4 * 137.0)
    argpdot = temp1 * -0.5 * con42 + temp2 * 0.0625 * (7.0 - cosio2 * 114.0 + cosio4 * 395.0) +
              temp3 * (3.0 - cosio2 * 36.0 + cosio4 * 49.0)
    xhdot1 = -temp1 * cosio
    nodedot = xhdot1 +
              (temp2 * 0.5 * (4.0 - cosio2 * 19.0) + temp3 * 2.0 * (3.0 - cosio2 * 7.0)) * cosio
    xpidot = argpdot + nodedot
    omgcof = bstar * cc3 * cos(argpo)
    xmcof = 0.0
    ecco > 1e-4 && (xmcof = -x2o3 * coef * bstar / eeta)
    xnodcf = omeosq * 3.5 * xhdot1 * cc1
    t2cof = cc1 * 1.5
    # Guard against a division by zero at an inclination of 180 degrees.
    xlcof = if abs(cosio + 1.0) > 1.5e-12
        j3oj2 * -0.25 * sinio * (cosio * 5.0 + 3.0) / (cosio + 1.0)
    else
        j3oj2 * -0.25 * sinio * (cosio * 5.0 + 3.0) / temp4
    end
    aycof = j3oj2 * -0.5 * sinio
    d = eta * cos(mo) + 1.0
    delmo = d * (d * d)
    sinmao = sin(mo)
    x7thm1 = cosio2 * 7.0 - 1.0

    # Deep space: period of 225 minutes or more.
    dodeep = _SGP4_TWOPI / no >= 225.0
    z = 0.0
    ds = (; irez=0, dedt=z, didt=z, dmdt=z, dnodt=z, domdt=z, d2201=z, d2211=z, d3210=z,
        d3222=z, d4410=z, d4422=z, d5220=z, d5232=z, d5421=z, d5433=z, del1=z, del2=z,
        del3=z, xfact=z, xlamo=z)
    pr = (; zmol=z, zmos=z, se2=z, se3=z, si2=z, si3=z, sl2=z, sl3=z, sl4=z, sgh2=z,
        sgh3=z, sgh4=z, sh2=z, sh3=z, ee2=z, e3=z, xi2=z, xi3=z, xl2=z, xl3=z, xl4=z,
        xgh2=z, xgh3=z, xgh4=z, xh2=z, xh3=z)
    if dodeep
        dosimp = true
        c = _sgp4_dscom(epoch, ecco, argpo, inclo, nodeo, no)
        # The call of ZZDSPR at this place in CSPICE changes nothing: the offsets
        # PEO, PINCO, PLO, PGHO and PHO of the periodics stay at zero.
        ds = _sgp4_dsinit(xke, c, argpo, gsto, mo, mdot, no, nodeo, nodedot, xpidot, ecco,
            eccsq, inclo)
        pr = (; c.zmol, c.zmos, c.se2, c.se3, c.si2, c.si3, c.sl2, c.sl3, c.sl4, c.sgh2,
            c.sgh3, c.sgh4, c.sh2, c.sh3, c.ee2, c.e3, c.xi2, c.xi3, c.xl2, c.xl3, c.xl4,
            c.xgh2, c.xgh3, c.xgh4, c.xh2, c.xh3)
    end
    d2 = d3 = d4 = t3cof = t4cof = t5cof = 0.0
    if !dosimp
        cc1sq = cc1 * cc1
        d2 = ao * 4.0 * tsi * cc1sq
        temp = d2 * tsi * cc1 / 3.0
        d3 = (ao * 17.0 + sfour) * temp
        d4 = temp * 0.5 * ao * tsi * (ao * 221.0 + sfour * 31.0) * cc1
        t3cof = d2 + cc1sq * 2.0
        t4cof = (d3 * 3.0 + cc1 * (d2 * 12.0 + cc1sq * 10.0)) * 0.25
        t5cof = (d4 * 3.0 + cc1 * 12.0 * d3 + d2 * 6.0 * d2 + cc1sq * 15.0 * (d2 * 2.0 + cc1sq)) * 0.2
    end
    return (; opmode, dodeep, dosimp, j2, j3oj2, xke, er, bstar, inclo, nodeo, ecco, argpo,
        mo, no, gsto, con41, x1mth2, x7thm1, cc1, cc4, cc5, mdot, argpdot, nodedot, omgcof,
        xmcof, xnodcf, t2cof, xlcof, aycof, delmo, sinmao, eta, d2, d3, d4, t3cof, t4cof,
        t5cof, ds, pr)
end

"`x` to the power `p` for plain numbers; a negative `x` gives NaN, as `pow` of C, and no error."
@inline _sgp4_realpow(x::Float64, p::Float64) = x < 0.0 ? NaN : x^p

"""
    _sgp4_dspace(r, t::T, eccm, argpm, inclm, mm, nodem) -> (eccm, argpm, inclm, mm, nodem, xn)

Secular effects of the Sun and the Moon and the resonance effects on the mean
elements of a deep-space orbit at `t` minutes from the epoch, as `zzdspc`. `xn` is the
mean motion with the resonance term.

For a resonant orbit the mean longitude and the mean motion are integrated from the
epoch in steps of 720 minutes with plain numbers (the count of steps follows
`primal(t)`) and the last part of a step is a polynomial in the type of `t`. CSPICE
can continue the integration from the time of the previous call; the steps are the
same from the epoch, so the result is the same. Loop invariant: `atime` is a whole
number of steps from the epoch, and `xli`, `xni` are the integrated values at `atime`.
"""
@inline function _sgp4_dspace(r, t::T, eccm, argpm, inclm, mm, nodem) where {T<:Real}
    ds = r.ds
    fasx2 = 0.13130908
    fasx4 = 2.8843198
    fasx6 = 0.37448087
    g22 = 5.7686396
    g32 = 0.95240898
    g44 = 1.8014998
    g52 = 1.050833
    g54 = 4.4108898
    stepp = 720.0
    step2 = 259200.0
    theta = _sgp4_mod(r.gsto + t * _SGP4_RPTIM, _SGP4_TWOPI)
    eccm = eccm + ds.dedt * t
    inclm = inclm + ds.didt * t
    argpm = argpm + ds.domdt * t
    nodem = nodem + ds.dnodt * t
    mm = mm + ds.dmdt * t
    xn = convert(T, r.no)
    ds.irez == 0 && return eccm, argpm, inclm, mm, nodem, xn
    t0 = primal(t)
    atime = 0.0
    xni = r.no
    xli = ds.xlamo
    delt = t0 > 0.0 ? stepp : -stepp
    xndt = xldot = xnddt = 0.0
    while true
        if ds.irez != 2
            # Terms of the near-synchronous resonance.
            xndt = ds.del1 * sin(xli - fasx2) + ds.del2 * sin((xli - fasx4) * 2.0) +
                   ds.del3 * sin((xli - fasx6) * 3.0)
            xldot = xni + ds.xfact
            xnddt = ds.del1 * cos(xli - fasx2) + ds.del2 * 2.0 * cos((xli - fasx4) * 2.0) +
                    ds.del3 * 3.0 * cos((xli - fasx6) * 3.0)
            xnddt *= xldot
        else
            # Terms of the near-half-day resonance.
            xomi = r.argpo + r.argpdot * atime
            x2omi = xomi + xomi
            x2li = xli + xli
            xndt = ds.d2201 * sin(x2omi + xli - g22) + ds.d2211 * sin(xli - g22) +
                   ds.d3210 * sin(xomi + xli - g32) + ds.d3222 * sin(-xomi + xli - g32) +
                   ds.d4410 * sin(x2omi + x2li - g44) + ds.d4422 * sin(x2li - g44) +
                   ds.d5220 * sin(xomi + xli - g52) + ds.d5232 * sin(-xomi + xli - g52) +
                   ds.d5421 * sin(xomi + x2li - g54) + ds.d5433 * sin(-xomi + x2li - g54)
            xldot = xni + ds.xfact
            xnddt = ds.d2201 * cos(x2omi + xli - g22) + ds.d2211 * cos(xli - g22) +
                    ds.d3210 * cos(xomi + xli - g32) + ds.d3222 * cos(-xomi + xli - g32) +
                    ds.d5220 * cos(xomi + xli - g52) + ds.d5232 * cos(-xomi + xli - g52) +
                    (ds.d4410 * cos(x2omi + x2li - g44) + ds.d4422 * cos(x2li - g44) +
                     ds.d5421 * cos(xomi + x2li - g54) + ds.d5433 * cos(-xomi + x2li - g54)) * 2.0
            xnddt *= xldot
        end
        abs(t0 - atime) >= stepp || break
        xli = xli + xldot * delt + xndt * step2
        xni = xni + xndt * delt + xnddt * step2
        atime += delt
    end
    ft = t - atime
    xn = xni + xndt * ft + xnddt * ft * ft * 0.5
    xl = xli + xldot * ft + xndt * ft * ft * 0.5
    if ds.irez != 1
        mm = xl - nodem * 2.0 + theta * 2.0
    else
        mm = xl - nodem - argpm + theta
    end
    dndt = xn - r.no
    xn = r.no + dndt
    return eccm, argpm, inclm, mm, nodem, xn
end

"""
    _sgp4_dpper(r, t::T, eccp, inclp, nodep, argpp, mp) -> (eccp, inclp, nodep, argpp, mp)

Long-period periodic effects of the Sun and the Moon on the elements of a deep-space
orbit at `t` minutes from the epoch, as `zzdspr` outside the initialization. For an
inclination below 0.2 radians the node and the argument of perigee change by the
Lyddane form; the choice follows the plain value of the inclination.
"""
@inline function _sgp4_dpper(r, t::T, eccp, inclp, nodep, argpp, mp) where {T<:Real}
    p = r.pr
    zes = 0.01675
    zel = 0.0549
    zns = 1.19459e-5
    znl = 1.5835218e-4
    zm = p.zmos + zns * t
    zf = zm + zes * 2.0 * sin(zm)
    sinzf, coszf = sincos(zf)
    f2 = sinzf * 0.5 * sinzf - 0.25
    f3 = sinzf * -0.5 * coszf
    ses = p.se2 * f2 + p.se3 * f3
    sis = p.si2 * f2 + p.si3 * f3
    sls = p.sl2 * f2 + p.sl3 * f3 + p.sl4 * sinzf
    sghs = p.sgh2 * f2 + p.sgh3 * f3 + p.sgh4 * sinzf
    shs = p.sh2 * f2 + p.sh3 * f3
    zm = p.zmol + znl * t
    zf = zm + zel * 2.0 * sin(zm)
    sinzf, coszf = sincos(zf)
    f2 = sinzf * 0.5 * sinzf - 0.25
    f3 = sinzf * -0.5 * coszf
    sel = p.ee2 * f2 + p.e3 * f3
    sil = p.xi2 * f2 + p.xi3 * f3
    sll = p.xl2 * f2 + p.xl3 * f3 + p.xl4 * sinzf
    sghl = p.xgh2 * f2 + p.xgh3 * f3 + p.xgh4 * sinzf
    shl = p.xh2 * f2 + p.xh3 * f3
    pe = ses + sel
    pinc = sis + sil
    pl = sls + sll
    pgh = sghs + sghl
    ph = shs + shl
    inclp = inclp + pinc
    eccp = eccp + pe
    sinip, cosip = sincos(inclp)
    if primal(inclp) >= 0.2
        ph = ph / sinip
        pgh = pgh - cosip * ph
        argpp = argpp + pgh
        nodep = nodep + ph
        mp = mp + pl
    else
        # Periodics with the Lyddane modification.
        sinop, cosop = sincos(nodep)
        alfdp = sinip * sinop
        betdp = sinip * cosop
        dalf = ph * cosop + pinc * cosip * sinop
        dbet = -ph * sinop + pinc * cosip * cosop
        alfdp = alfdp + dalf
        betdp = betdp + dbet
        nodep = _sgp4_mod(nodep, _SGP4_TWOPI)
        (primal(nodep) < 0.0 && r.opmode == _SGP4_AFSPC) && (nodep = nodep + _SGP4_TWOPI)
        xls = mp + argpp + cosip * nodep
        dls = pl + pgh - pinc * nodep * sinip
        xls = xls + dls
        xnoh = primal(nodep)
        nodep = atan(alfdp, betdp)
        (primal(nodep) < 0.0 && r.opmode == _SGP4_AFSPC) && (nodep = nodep + _SGP4_TWOPI)
        if abs(xnoh - primal(nodep)) > Float64(pi)
            if primal(nodep) < xnoh
                nodep = nodep + _SGP4_TWOPI
            else
                nodep = nodep - _SGP4_TWOPI
            end
        end
        mp = mp + pl
        argpp = xls - mp - cosip * nodep
    end
    return eccp, inclp, nodep, argpp, mp
end

"""
    _sgp4_kepler(u::Float64, axnl::Float64, aynl::Float64) -> Float64

Argument `e` at which CSPICE takes the sine and the cosine of the eccentric longitude
in the Kepler iteration of `xxsgp4e`: Newton steps on
`u - aynl cos(e) + axnl sin(e) - e` from `e = u`, each limited to 1 in size, until a
step is below 1e-12 or 10 steps are made. CSPICE uses the sine and the cosine that it
computed before the last step, so the result is the iterate before that step.
Loop invariant: `e` is the iterate that precedes the step `tem5`.
"""
@inline function _sgp4_kepler(u::Float64, axnl::Float64, aynl::Float64)
    eo1 = u
    e = u
    temp = 9999.9
    iter = 0
    while temp >= 1e-12 && iter < 10
        iter += 1
        e = eo1
        sineo1, coseo1 = sincos(eo1)
        tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl
        tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5
        temp = abs(tem5)
        temp > 1.0 && (tem5 /= temp)
        eo1 += tem5
    end
    return e
end

"""
    _sgp4_eval(r, t::T) -> (pos::Vec3{T}, vel::Vec3{T})

Position (km) and velocity (km/s) in the TEME frame at `t` minutes from the element
epoch, for the record `r` of `_sgp4_init`, as `xxsgp4e`.

Contract: with a plain `t` the result is the CSPICE state to rounding level. With a
series-valued `t` both results are the series of the CSPICE position and velocity
functions: the Kepler equation is solved with plain numbers as CSPICE does, and the
eccentric longitude then gets its derivative parts from Newton steps on the residual
with its plain value removed (`_with_primal`). Errors, with the kinds of CSPICE:
`:badmeanmotion` (mean motion not positive), `:badmeccentricity` (mean eccentricity
outside `[-0.001, 1)`), `:badmsemimajor` (mean semi-major axis below 0.95 Earth
radii), `:badpeccentricity` (perturbed eccentricity outside `[0, 1]`),
`:badsemilatus` (negative semi-latus rectum) and `:orbitdecay` (radius below one
Earth radius). No allocation.
"""
@inline function _sgp4_eval(r, t::T) where {T<:Real}
    x2o3 = _SGP4_X2O3
    temp4 = 1.5e-12
    kps = r.er * r.xke / 60.0
    # Secular gravity and atmospheric drag.
    xmdf = r.mo + r.mdot * t
    omgadf = r.argpo + r.argpdot * t
    xnoddf = r.nodeo + r.nodedot * t
    argpm = omgadf
    mm = xmdf
    t2 = t * t
    nodem = xnoddf + r.xnodcf * t2
    tempa = 1.0 - r.cc1 * t
    tempe = r.bstar * r.cc4 * t
    templ = r.t2cof * t2
    if !r.dosimp
        delomg = r.omgcof * t
        d = r.eta * cos(xmdf) + 1.0
        delm = r.xmcof * (d * (d * d) - r.delmo)
        temp = delomg + delm
        mm = xmdf + temp
        argpm = omgadf - temp
        t3 = t2 * t
        t4 = t3 * t
        tempa = tempa - r.d2 * t2 - r.d3 * t3 - r.d4 * t4
        tempe = tempe + r.bstar * r.cc5 * (sin(mm) - r.sinmao)
        templ = templ + r.t3cof * t3 + t4 * (r.t4cof + t * r.t5cof)
    end
    xn = convert(T, r.no)
    eccm = convert(T, r.ecco)
    inclm = convert(T, r.inclo)
    if r.dodeep
        eccm, argpm, inclm, mm, nodem, xn = _sgp4_dspace(r, t, eccm, argpm, inclm, mm, nodem)
    end
    primal(xn) <= 0.0 && _sgp4_error(:badmeanmotion,
        "Mean motion less-than zero. This error may indicate a bad TLE set.")
    am = _sgp4_pow(r.xke / xn, x2o3) * (tempa * tempa)
    xn = r.xke / _sgp4_pow(am, 1.5)
    eccm = eccm - tempe
    (primal(eccm) >= 1.0 || primal(eccm) < -0.001) && _sgp4_error(:badmeccentricity,
        "Mean eccentricity value beyond allowed bounds [-0.001,1.0). This error may " *
        "indicate a bad TLE set.")
    primal(am) < 0.95 && _sgp4_error(:badmsemimajor,
        "Mean semi-major axis value below allowed minimum of 0.95. This error may " *
        "indicate a bad TLE set or a decayed orbit.")
    primal(eccm) < 1e-6 && (eccm = convert(T, 1e-6))
    mm = mm + r.no * templ
    xlm = mm + argpm + nodem
    nodem = _sgp4_mod(nodem, _SGP4_TWOPI)
    argpm = _sgp4_mod(argpm, _SGP4_TWOPI)
    xlm = _sgp4_mod(xlm, _SGP4_TWOPI)
    mm = _sgp4_mod(xlm - argpm - nodem, _SGP4_TWOPI)
    # Lunar and solar periodics.
    eccp = eccm
    xincp = inclm
    argpp = argpm
    nodep = nodem
    mp = mm
    sinip, cosip = sincos(inclm)
    aycof = convert(T, r.aycof)
    xlcof = convert(T, r.xlcof)
    con41 = convert(T, r.con41)
    x1mth2 = convert(T, r.x1mth2)
    x7thm1 = convert(T, r.x7thm1)
    if r.dodeep
        eccp, xincp, nodep, argpp, mp = _sgp4_dpper(r, t, eccp, xincp, nodep, argpp, mp)
        if primal(xincp) < 0.0
            xincp = -xincp
            nodep = nodep + Float64(pi)
            argpp = argpp - Float64(pi)
        end
        (primal(eccp) < 0.0 || primal(eccp) > 1.0) && _sgp4_error(:badpeccentricity,
            "Perturbed eccentricity value beyond allowed bounds [0,1]. This error may " *
            "indicate a bad TLE set.")
        # Long-period periodics of a deep-space orbit.
        sinip, cosip = sincos(xincp)
        aycof = r.j3oj2 * -0.5 * sinip
        xlcof = if abs(primal(cosip) + 1.0) > 1.5e-12
            r.j3oj2 * -0.25 * sinip * (cosip * 5.0 + 3.0) / (cosip + 1.0)
        else
            r.j3oj2 * -0.25 * sinip * (cosip * 5.0 + 3.0) / temp4
        end
    end
    sargpp, cargpp = sincos(argpp)
    axnl = eccp * cargpp
    temp = 1.0 / (am * (1.0 - eccp * eccp))
    aynl = eccp * sargpp + temp * aycof
    xl = mp + argpp + nodep + temp * xlcof * axnl
    # Kepler equation: plain solution as CSPICE, then the derivative parts.
    u = _sgp4_mod(xl - nodep, _SGP4_TWOPI)
    eo1 = convert(T, _sgp4_kepler(primal(u), primal(axnl), primal(aynl)))
    for _ in 1:_ad_order(T)       # invariant: one more derivative order of eo1 is exact
        s, c = sincos(eo1)
        res = u - aynl * c + axnl * s - eo1
        eo1 = eo1 + _with_primal(res, 0.0) / (1.0 - c * axnl - s * aynl)
    end
    sineo1, coseo1 = sincos(eo1)
    # Short-period preliminary quantities.
    ecose = axnl * coseo1 + aynl * sineo1
    esine = axnl * sineo1 - aynl * coseo1
    el2 = axnl * axnl + aynl * aynl
    pl = am * (1.0 - el2)
    primal(pl) < 0.0 && _sgp4_error(:badsemilatus, "Semi-latus rectum less-than zero.")
    rl = am * (1.0 - ecose)
    rdotl = sqrt(am) * esine / rl
    rvdotl = sqrt(pl) / rl
    betal = sqrt(1.0 - el2)
    temp = esine / (betal + 1.0)
    sinu = am / rl * (sineo1 - aynl - axnl * temp)
    cosu = am / rl * (coseo1 - axnl + aynl * temp)
    su = atan(sinu, cosu)
    sin2u = (cosu + cosu) * sinu
    cos2u = 1.0 - sinu * 2.0 * sinu
    temp = 1.0 / pl
    temp1 = r.j2 * 0.5 * temp
    temp2 = temp1 * temp
    # Short-period periodics of a deep-space orbit.
    if r.dodeep
        cosisq = cosip * cosip
        con41 = cosisq * 3.0 - 1.0
        x1mth2 = 1.0 - cosisq
        x7thm1 = cosisq * 7.0 - 1.0
    end
    mr = rl * (1.0 - temp2 * 1.5 * betal * con41) + temp1 * 0.5 * x1mth2 * cos2u
    su = su - temp2 * 0.25 * x7thm1 * sin2u
    xnode = nodep + temp2 * 1.5 * cosip * sin2u
    xinc = xincp + temp2 * 1.5 * cosip * sinip * cos2u
    mv = rdotl - xn * temp1 * x1mth2 * sin2u / r.xke
    rvdot = rvdotl + xn * temp1 * (x1mth2 * cos2u + con41 * 1.5) / r.xke
    # Orientation vectors.
    sinsu, cossu = sincos(su)
    snod, cnod = sincos(xnode)
    sini, cosi = sincos(xinc)
    xmx = -snod * cosi
    xmy = cnod * cosi
    ux = xmx * sinsu + cnod * cossu
    uy = xmy * sinsu + snod * cossu
    uz = sini * sinsu
    vx = xmx * cossu - cnod * sinsu
    vy = xmy * cossu - snod * sinsu
    vz = sini * cossu
    pos = Vec3{T}(mr * ux * r.er, mr * uy * r.er, mr * uz * r.er)
    vel = Vec3{T}((mv * ux + rvdot * vx) * kps, (mv * uy + rvdot * vy) * kps,
        (mv * uz + rvdot * vz) * kps)
    primal(mr) < 1.0 && _sgp4_error(:orbitdecay, "Satellite has decayed.")
    return pos, vel
end
