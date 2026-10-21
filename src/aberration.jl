# Aberration corrections: one-way light time and stellar aberration, for reception and
# transmission, as CSPICE N0067 applies them (`spkez`, `spkezp`, `spkacs`, `spkaps`,
# `spkltc`, `spkapo`, `zzstelab`, `stelab`), and the time derivatives of any order of
# the corrected position.
#
# Every function is a pure function of the `KernelSet` snapshot and of its arguments,
# works on `StaticArrays` values and plain bits, and does not allocate. The numeric
# functions are generic over the epoch type (`Float64`, a dual number); all discrete
# decisions use `primal`.

# --- the correction as a plain-bits value -------------------------------------------

const _AB_LT = 0x01     # light time correction
const _AB_CN = 0x02     # converged Newtonian light time
const _AB_S = 0x04      # stellar aberration
const _AB_X = 0x08      # transmission (the default is reception)

"""
    Aberration(abcorr::AbstractString)
    Aberration(; lt=false, converged=false, stellar=false, transmit=false)

An aberration correction as an `isbits` value. Every call that takes a correction
string also takes an `Aberration`; the string is then not parsed at each call.

| String             | Correction                                                          |
|:-------------------|:--------------------------------------------------------------------|
| `"NONE"`           | none: the geometric state                                           |
| `"LT"`             | reception, one-way light time, one iteration                        |
| `"LT+S"`           | `"LT"` and stellar aberration                                       |
| `"CN"`             | reception, light time iterated to convergence                       |
| `"CN+S"`           | `"CN"` and stellar aberration                                       |
| `"XLT"`, `"XLT+S"` | the `"LT"` forms for transmission                                   |
| `"XCN"`, `"XCN+S"` | the `"CN"` forms for transmission                                   |

Reception gives the target as the observer sees it at `et`: the target at
`et - lt`. Transmission gives the target where a signal sent at `et` reaches it: the
target at `et + lt`. Letter case and blanks in the string are ignored.

The keyword form builds the value from flags; `converged`, `stellar` and `transmit`
each imply `lt`. A value `c` has the properties `c.uselt` (the flag `lt`),
`c.converged`, `c.stellar`, `c.transmit` and `c.geometric` (true for `"NONE"` only).
`string(c)` gives the CSPICE spelling.

Throws `SpiceError(:invalidoption)` if the string is not one of the nine above. This includes
the values that CSPICE recognizes and rejects: `RL`, `RL+S`, `XRL`, `XRL+S`, `S`, `XS`.

# Examples
```jldoctest
julia> c = Aberration(" lt + s ")
Aberration("LT+S")

julia> c.stellar, c.transmit, c.geometric
(true, false, false)

julia> Aberration(converged=true, transmit=true)
Aberration("XCN")
```

See also [`spkezr`](@ref), [`spkpos`](@ref).
"""
struct Aberration
    flags::UInt8
    global _aberration(flags::UInt8) = new(flags)
end

Aberration(c::Aberration) = c

function Aberration(; lt::Bool=false, converged::Bool=false, stellar::Bool=false,
        transmit::Bool=false)
    uselt = lt | converged | stellar | transmit
    return _aberration((uselt ? _AB_LT : 0x00) | (converged ? _AB_CN : 0x00) |
                       (stellar ? _AB_S : 0x00) | (transmit ? _AB_X : 0x00))
end

@inline function Base.getproperty(c::Aberration, s::Symbol)
    f = getfield(c, :flags)
    s === :uselt && return f & _AB_LT != 0
    s === :converged && return f & _AB_CN != 0
    s === :stellar && return f & _AB_S != 0
    s === :transmit && return f & _AB_X != 0
    s === :geometric && return f == 0x00
    return getfield(c, s)
end
Base.propertynames(::Aberration) = (:uselt, :converged, :stellar, :transmit, :geometric, :flags)

"The CSPICE spelling of the correction `c`."
function Base.string(c::Aberration)
    c.geometric && return "NONE"
    return (c.transmit ? "X" : "") * (c.converged ? "CN" : "LT") * (c.stellar ? "+S" : "")
end
Base.show(io::IO, c::Aberration) = print(io, "Aberration(\"", string(c), "\")")

"`+1` for transmission, `-1` for reception: the sign of the light time in the epoch."
@inline _ltsign(c::Aberration) = c.transmit ? 1 : -1

"The characters of `s` (at most 8, ASCII) packed in one integer, most significant first."
function _abcode(s::AbstractString)
    code = UInt64(0)
    for ch in s     # invariant: `code` holds the characters before `ch`
        code = (code << 8) | UInt64(UInt8(ch))
    end
    return code
end

@noinline _throw_abcorr(s, why) = throw(SpiceError(:invalidoption,
    "aberration correction \"$s\" $why; the valid values are NONE, LT, LT+S, CN, CN+S, " *
    "XLT, XLT+S, XCN, XCN+S"))

"""
    _spk_abcorr(abcorr) -> Aberration

The correction for `spkpos` and `spkapo`: CSPICE signals SPICE(SPKINVALIDOPTION) there
for a string that is not valid, where `spkezr` signals SPICE(INVALIDOPTION).
"""
@inline _spk_abcorr(c::Aberration) = c
@inline function _spk_abcorr(s::AbstractString)
    try
        return Aberration(s)
    catch e
        (e isa SpiceError && e.kind === :invalidoption) || rethrow()
        throw(SpiceError(:spkinvalidoption, e.msg))
    end
end

"Parse a correction string that is not one of the nine common spellings."
@noinline function _parse_abcorr(s::AbstractString)
    code = UInt64(0)
    n = 0
    # invariant: `code` holds the `n` non-blank characters read so far, in upper case
    for ch in s
        ch == ' ' && continue
        (isascii(ch) && n < 5) || _throw_abcorr(s, "is not recognized")
        code = (code << 8) | UInt64(UInt8(uppercase(ch)))
        n += 1
    end
    code == _abcode("NONE") && return _aberration(0x00)
    code == _abcode("LT") && return _aberration(_AB_LT)
    code == _abcode("LT+S") && return _aberration(_AB_LT | _AB_S)
    code == _abcode("CN") && return _aberration(_AB_LT | _AB_CN)
    code == _abcode("CN+S") && return _aberration(_AB_LT | _AB_CN | _AB_S)
    code == _abcode("XLT") && return _aberration(_AB_X | _AB_LT)
    code == _abcode("XLT+S") && return _aberration(_AB_X | _AB_LT | _AB_S)
    code == _abcode("XCN") && return _aberration(_AB_X | _AB_LT | _AB_CN)
    code == _abcode("XCN+S") && return _aberration(_AB_X | _AB_LT | _AB_CN | _AB_S)
    (code == _abcode("RL") || code == _abcode("RL+S") || code == _abcode("XRL") ||
     code == _abcode("XRL+S")) &&
        _throw_abcorr(s, "calls for relativistic corrections, which are not supported")
    (code == _abcode("S") || code == _abcode("XS")) && _throw_abcorr(s,
        "calls for stellar aberration without light time, which is not supported")
    _throw_abcorr(s, "is not recognized")
end

@inline function Aberration(s::AbstractString)
    # Fast path: the common spellings, compared as they are.
    s == "NONE" && return _aberration(0x00)
    s == "LT" && return _aberration(_AB_LT)
    s == "LT+S" && return _aberration(_AB_LT | _AB_S)
    s == "CN" && return _aberration(_AB_LT | _AB_CN)
    s == "CN+S" && return _aberration(_AB_LT | _AB_CN | _AB_S)
    s == "XLT" && return _aberration(_AB_X | _AB_LT)
    s == "XLT+S" && return _aberration(_AB_X | _AB_LT | _AB_S)
    s == "XCN" && return _aberration(_AB_X | _AB_LT | _AB_CN)
    s == "XCN+S" && return _aberration(_AB_X | _AB_LT | _AB_CN | _AB_S)
    return _parse_abcorr(s)
end

const AbcorrArg = Union{AbstractString,Aberration}

# --- small vector helpers -----------------------------------------------------------
#
# Written out so that they hold for every element type (`Taylor` has no comparisons).

@inline _dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _norm3(a) = sqrt(_dot3(a, a))

"Unit vector along `a`, or the zero vector if `a` is zero, as `vhat`."
@inline function _vhat(a)
    n = _norm3(a)
    return primal(n) == 0 ? zero(a) : a / n
end

"Component of `a` perpendicular to `b`; `a` itself if `b` is zero, as `vperp`."
@inline function _vperp(a, b)
    bb = _dot3(b, b)
    return primal(bb) == 0 ? a - zero(b) : a - (_dot3(a, b) / bb) * b
end

"Unit vector along `r` and its time derivative for the velocity `v`, as `dvhat`."
@inline function _dvhat(r, v)
    n = _norm3(r)
    primal(n) == 0 && return zero(r) - zero(v), v - zero(r)
    u = r / n
    return u - zero(v), _vperp(v, u) / n
end

@inline _ssb(ks::KernelSet, body::Int32, et, iref::Int32) =
    spk_jet(ks, body, et, iref, Int32(0), Val(1))

@inline _ssbpos(ks::KernelSet, body::Int32, et, iref::Int32) =
    spk_jet(ks, body, et, iref, Int32(0), Val(0))[1]

"True if frame `ref` is inertial (class 1): the test `type == 1` of `spkez`."
@inline _is_inertial(ks::KernelSet, ref::Int32) =
    ref == Int32(1) || frame_def(ks, ref).class == Int32(1)

@noinline _throw_not_inertial(ks, ref) = throw(SpiceError(:badarg,
    "the frame with ID $ref is not an inertial frame; this call needs one"))

@inline function _need_inertial(ks::KernelSet, ref::Int32)
    _is_inertial(ks, ref) || _throw_not_inertial(ks, ref)
    return ref
end

@noinline _throw_range_rate() = throw(SpiceError(:badarg,
    "the range rate of the target is near the speed of light: the light time " *
    "derivative cannot be computed"))

# --- light time ---------------------------------------------------------------------

"""
    _ltc(ks, targ, et, iref, corr, robs, vobs) -> (; r, v, lt, dlt, ltprev, k)

The computation of `spkltc`: state `(r, v)` of `targ` relative to an observer whose
barycentric state at `et` is `(robs, vobs)`, in the inertial frame `iref`, corrected
for light time as `corr` says; the light time `lt` and its rate `dlt`.

Contract: `iref` is inertial. `NONE` gives the geometric state. `LT` forms make one
iteration, `CN` forms at most five and stop when the light time changes by less than
`1e-17 * max(1, |epoch|)`. Postconditions: `r` is the target position at
`et + ltsign * ltprev` minus `robs`, `lt = |r| / c`, and `k = ltsign * (r̂ · v_targ) / c`
is the derivative of the light-time map with respect to the light time (`ltprev` and
`k` are zero when no iteration is made). Stellar aberration is not applied here.
"""
function _ltc(ks::KernelSet, targ::Int32, et, iref::Int32, corr::Aberration, robs, vobs)
    rt, vt = _ssb(ks, targ, et, iref)
    r = rt - robs
    v = vt - vobs
    dist = _norm3(r)
    lt = dist / CLIGHT
    z = zero(lt)
    primal(lt) == 0 && return (; r, v, lt, dlt=z, ltprev=z, k=z)
    corr.uselt || return (; r, v, lt, dlt=_dot3(r, v) / (dist * CLIGHT), ltprev=z, k=z)
    s = _ltsign(corr)
    numitr = corr.converged ? 5 : 1
    ltprev = lt
    lterr = 1.0
    i = 0
    # invariant: `lt` is the light time after `i` iterations, and for i >= 1 `r`, `vt`
    # are at the epoch `et + s * ltprev`
    while i < numitr && lterr > 1e-17
        epoch = et + s * lt
        rt, vt = _ssb(ks, targ, epoch, iref)
        r = rt - robs
        v = vt - vobs
        ltprev = lt
        lt = _norm3(r) / CLIGHT
        lterr = abs(primal(lt) - primal(ltprev)) / max(1.0, abs(primal(epoch)))
        i += 1
    end
    a = 1 / (CLIGHT * _norm3(r))
    b = _dot3(r, v)
    k = s * _dot3(r, vt) * a
    primal(k) > 1 - 1e-10 && _throw_range_rate()
    dlt = a * b / (1 - k)
    v = (1 + s * dlt) * vt - vobs
    return (; r, v, lt, dlt, ltprev, k)
end

"""
    _apo(ks, targ, et, iref, corr, robs, vobs) -> (r, lt)

The computation of `spkapo`: apparent position of `targ` relative to an observer with
the barycentric state `(robs, vobs)` at `et`, in the inertial frame `iref`.

Contract: `iref` is inertial. `LT` forms make one iteration and `CN` forms exactly
three (no convergence test); the stellar aberration is the rotation of `stelab`.
Postcondition: `lt` is the norm of the position before the stellar aberration,
divided by `c`.
"""
function _apo(ks::KernelSet, targ::Int32, et, iref::Int32, corr::Aberration, robs, vobs)
    r = _ssbpos(ks, targ, et, iref) - robs
    lt = _norm3(r) / CLIGHT
    s = _ltsign(corr)
    maxitr = !corr.uselt ? 0 : corr.converged ? 3 : 1
    # invariant: `lt` is the light time after i - 1 iterations
    for _ in 1:maxitr
        r = _ssbpos(ks, targ, et + s * lt, iref) - robs
        lt = _norm3(r) / CLIGHT
    end
    corr.stellar || return r, lt
    return _stelab(r, corr.transmit ? -vobs : vobs), lt
end

# --- stellar aberration -------------------------------------------------------------

@noinline _throw_vobs() = throw(SpiceError(:badarg,
    "the speed of the observer is not below the speed of light"))

"""
    _stelab(pobj, vobs) -> Vec3

The computation of `stelab`: `pobj` rotated by the angle `asin(|u × vobs/c|)` about
`u × vobs/c`, with `u` the unit vector along `pobj`. Requires `|vobs| < c`
(`SpiceError(:badarg)` otherwise). Postcondition: the norm of `pobj` is kept.
"""
@inline function _stelab(pobj, vobs)
    u = _vhat(pobj)
    vbyc = (1 / CLIGHT) * vobs
    primal(_dot3(vbyc, vbyc)) >= 1 && _throw_vobs()
    h = cross(u, vbyc)
    sinphi = _norm3(h)
    primal(sinphi) == 0 && return pobj - zero(h)
    phi = asin(sinphi)
    # `vrotv`: rotate `pobj` about the axis `h` by `phi`.
    x = h / sinphi
    p = _dot3(pobj, x) * x
    v1 = pobj - p
    v2 = cross(x, v1)
    sn, cs = sincos(phi)
    return cs * v1 + sn * v2 + p
end

"The stellar aberration offset of `zzstelab` for the position `r`, the unit vector `rhat` and the observer velocity `vobs`."
@inline function _stel_offset(r, rhat, vobs)
    vp = _vperp(vobs, rhat)
    s = _norm3(vp) / CLIGHT
    c2 = 1 - s * s
    c = primal(c2) > 0 ? sqrt(c2) : zero(c2)
    vphat = _vhat(vp)
    ptgmag = _norm3(r)
    return (ptgmag * s) * vphat + (ptgmag * (c - 1)) * rhat, vp, vphat, s, c, ptgmag
end

@noinline _throw_stelab_cos() = throw(SpiceError(:badarg,
    "the cosine of the stellar aberration angle is zero: the observer velocity is not " *
    "realistic"))

"""
    _zzstelab(transmit, accobs, vobs, r, v) -> (scorr, dscorr)

The computation of `zzstelab`: the offset `scorr` that stellar aberration adds to the
position `r` (velocity `v`) for an observer with barycentric velocity `vobs` and
acceleration `accobs`, and its time derivative `dscorr`. For transmission the velocity
and the acceleration change sign. As in CSPICE the derivative is analytical if the
sine of the aberration angle is at least `1e-6`, and a central difference of the
offset over ±1 s of linear motion otherwise.
"""
function _zzstelab(transmit::Bool, accobs, vobs, r, v)
    lcvobs = transmit ? -vobs : +vobs
    lcacc = transmit ? -accobs : +accobs
    rhat, drhat = _dvhat(r, v)
    dptmag = _dot3(v, rhat)
    scorr, vp, vphat, s, c, ptgmag = _stel_offset(r, rhat, lcvobs)
    primal(c) == 0 && _throw_stelab_cos()
    if primal(s) >= 1e-6
        dvp = lcacc + (-_dot3(lcvobs, drhat) - _dot3(lcacc, rhat)) * rhat +
              (-_dot3(lcvobs, rhat)) * drhat
        _, dvphat = _dvhat(vp, dvp)
        dphi = 1 / (c * CLIGHT) * _dot3(dvp, vphat)
        term1 = s * dvphat + (c * dphi) * vphat
        term2 = (c - 1) * drhat + (-s * dphi) * rhat
        dscorr = ptgmag * (term1 + term2) + (dptmag * s) * vphat + (dptmag * (c - 1)) * rhat
        return scorr, dscorr
    end
    em = r - v
    ep = r + v
    lo = _stel_offset(em, _vhat(em), lcvobs - lcacc)[1]
    hi = _stel_offset(ep, _vhat(ep), lcvobs + lcacc)[1]
    return scorr, 0.5 * hi - 0.5 * lo
end

"""
    _aps(ks, targ, et, iref, corr, robs, vobs, accobs) -> (; r, v, lt, dlt, ltprev, k)

The computation of `spkaps`: `_ltc`, and for the `+S` forms the stellar aberration
offset of `_zzstelab` and its derivative added to the state. `lt`, `dlt`, `ltprev`
and `k` are those of `_ltc`.
"""
@inline function _aps(ks::KernelSet, targ::Int32, et, iref::Int32, corr::Aberration, robs,
        vobs, accobs)
    q = _ltc(ks, targ, et, iref, corr, robs, vobs)
    corr.stellar || return q
    scorr, dscorr = _zzstelab(corr.transmit, accobs, vobs, q.r, q.v)
    return (; r=q.r + scorr, v=q.v + dscorr, q.lt, q.dlt, q.ltprev, q.k)
end

"""
    _obs_state(ks, obs, et, iref, corr) -> (robs, vobs, accobs)

Barycentric state of the observer at `et` and, for the `+S` forms, the acceleration
that `spkacs` uses: the central difference of the velocity over ±1 s. The
acceleration is zero for the other forms.
"""
@inline function _obs_state(ks::KernelSet, obs::Int32, et, iref::Int32, corr::Aberration)
    robs, vobs = _ssb(ks, obs, et, iref)
    corr.stellar || return robs, vobs, zero(vobs)
    vm = _ssb(ks, obs, et - 1.0, iref)[2]
    vp = _ssb(ks, obs, et + 1.0, iref)[2]
    return robs, vobs, 0.5 * vp - 0.5 * vm
end

# --- corrected states in any frame --------------------------------------------------

"""
    abcorr_state(ks, targ::Int32, et, ref::Int32, corr, obs::Int32) -> (SVector{6}, lt)

The computation of `spkez` for a correction other than `NONE`: state of `targ`
relative to `obs` in frame `ref`, and the light time.

Contract: for an inertial `ref` this is `spkacs`. For another frame the state is
computed in J2000 and rotated with the frame evaluated at `et + ltsign * ltc`, where
`ltc` is the light time between the observer and the center of the frame (zero if the
center is the observer, the target light time if it is the target); the derivative
block of the transformation is multiplied by `1 + ltsign * d ltc/dt`. No allocation.
"""
function abcorr_state(ks::KernelSet, targ::Int32, et, ref::Int32, corr::Aberration,
        obs::Int32)
    inertial = _is_inertial(ks, ref)
    iref = inertial ? ref : Int32(1)
    robs, vobs, accobs = _obs_state(ks, obs, et, iref, corr)
    q = _aps(ks, targ, et, iref, corr, robs, vobs, accobs)
    inertial && return vcat(q.r, q.v), q.lt
    center = frame_def(ks, ref).center
    s = _ltsign(corr)
    if center == obs
        ltc = zero(q.lt)
        dltc = zero(q.dlt)
    elseif center == targ
        ltc = q.lt
        dltc = q.dlt
    else
        qc = _ltc(ks, center, et, iref, corr, robs, vobs)
        ltc = qc.lt
        dltc = qc.dlt
    end
    R, dR = frame_jet(ks, Int32(1), ref, et + s * ltc, Val(1))
    return vcat(_m_mxv(R, q.r), (1 + s * dltc) * _m_mxv(dR, q.r) + _m_mxv(R, q.v)), q.lt
end

"""
    abcorr_position(ks, targ::Int32, et, ref::Int32, corr, obs::Int32) -> (SVector{3}, lt)

The computation of `spkezp` for a correction other than `NONE`: position of `targ`
relative to `obs` in frame `ref`, and the light time.

Contract: as `abcorr_state`, with the light-time and stellar-aberration rules of
`spkapo` (`_apo`). The position can therefore differ from the position part of
`abcorr_state` at rounding level, as in CSPICE. No allocation.
"""
function abcorr_position(ks::KernelSet, targ::Int32, et, ref::Int32, corr::Aberration,
        obs::Int32)
    inertial = _is_inertial(ks, ref)
    iref = inertial ? ref : Int32(1)
    robs, vobs = corr.stellar ? _ssb(ks, obs, et, iref) :
                 (r0 = _ssbpos(ks, obs, et, iref); (r0, zero(r0)))
    r, lt = _apo(ks, targ, et, iref, corr, robs, vobs)
    inertial && return r, lt
    center = frame_def(ks, ref).center
    ltc = center == obs ? zero(lt) : center == targ ? lt :
          _apo(ks, center, et, iref, corr, robs, vobs)[2]
    R = frame_rotation(ks, Int32(1), ref, et + _ltsign(corr) * ltc)
    return _m_mxv(R, r), lt
end

# --- derivatives of any order -------------------------------------------------------

"Taylor series of order N of each component from the jet `jet`, starting at element `first`."
@inline function _series(jet::NTuple{M,Vec3{T}}, first::Int, ::Val{N}) where {M,T,N}
    return Vec3{Taylor{N,T}}(ntuple(i -> Taylor{N,T}(ntuple(
        j -> jet[j+first-1][i] / factorial(j - 1), Val(N + 1))), Val(3)))
end

"""
    _lt_series(ks, body, et, t, iref, corr, robs, vobs, ro) -> Vec3{Taylor}

Series of the light-time corrected position of `body` relative to the observer whose
barycentric position series is `ro` (state `(robs, vobs)` at `et`); `t` is the epoch
as a series.

Contract: the result is `p(t) = r_body(t + ltsign * ℓ(t)) - r_obs(t)`, and its
coefficient 0 is the position of `_ltc`. For the `LT` forms `ℓ` is the geometric
light time `|r_body(t) - r_obs(t)| / c`, because the single iteration of CSPICE
evaluates the body at that light time. For the `CN` forms coefficient 0 of `ℓ` is the
light time `ltprev` at which `_ltc` evaluated the body, and the coefficients of order
1 to N satisfy the light-time equation `ℓ = |p| / c` order by order.
"""
function _lt_series(ks::KernelSet, body::Int32, et, t::Taylor{N,T}, iref::Int32,
        corr::Aberration, robs, vobs, ro) where {N,T}
    s = _ltsign(corr)
    geo = _ssbpos(ks, body, t, iref) - ro
    # Coincident bodies: the geometric series, as `spkltc` returns the geometric state.
    primal(_dot3(geo, geo)) == 0 && return geo
    if !corr.converged
        return _ssbpos(ks, body, t + s * (_norm3(geo) / CLIGHT), iref) - ro
    end
    q = _ltc(ks, body, et, iref, corr, robs, vobs)
    g = inv(1 - q.k)
    l0 = convert(T, q.ltprev)
    L = Taylor{N,T}(l0)
    # Newton steps on the series. Coefficient n of the map F(ℓ) = |p(ℓ)| / c is
    # `k * ℓ_n` plus a function of the lower coefficients, so one step with the factor
    # 1 / (1 - k) makes one more coefficient final.
    # invariant: before step j the coefficients 0 .. j-1 of L are final
    for _ in 1:N
        F = _norm3(_ssbpos(ks, body, t + s * L, iref) - ro) / CLIGHT
        L = _taylor(l0, map((l, f) -> l + (f - l) * g, L.d, F.d))
    end
    return _ssbpos(ks, body, t + s * L, iref) - ro
end

"Series of the position `p` with stellar aberration for the observer velocity series `vo`: the formula of `zzstelab`, `p + |p| s v̂p + |p| (c - 1) r̂`."
@inline function _stelab_series(transmit::Bool, vo, p)
    v = transmit ? -vo : +vo
    pm = _norm3(p)
    rhat = p / pm
    vp = v - _dot3(v, rhat) * rhat
    c = sqrt(1 - _dot3(vp, vp) / CLIGHT^2)
    return c * p + (pm / CLIGHT) * vp
end

"""
    abcorr_jet(ks, targ::Int32, et, ref::Int32, corr, obs::Int32, Val(N)) -> NTuple{N+1,Vec3}

Position of `targ` relative to `obs` in frame `ref` with the correction `corr`, and
its time derivatives to order N. See `spkder` for what is differentiated. No
allocation.
"""
function abcorr_jet(ks::KernelSet, targ::Int32, et::T, ref::Int32, corr::Aberration,
        obs::Int32, ::Val{N}) where {T<:Real,N}
    corr.geometric && return spk_jet(ks, targ, et, ref, obs, Val(N))
    # A body relative to itself: the zero state of `spkltc`, in every frame.
    targ == obs && return ntuple(_ -> zero(Vec3{T}), Val(N + 1))
    inertial = _is_inertial(ks, ref)
    iref = inertial ? ref : Int32(1)
    s = _ltsign(corr)
    ojet = spk_jet(ks, obs, et, iref, Int32(0), Val(N + 1))
    robs, vobs = ojet[1], ojet[2]
    ro = _series(ojet, 1, Val(N))
    t = taylor_var(et, Val(N))
    p = _lt_series(ks, targ, et, t, iref, corr, robs, vobs, ro)
    moving = primal(_dot3(p, p)) > 0
    q = corr.stellar && moving ? _stelab_series(corr.transmit, _series(ojet, 2, Val(N)), p) : p
    if !inertial
        center = frame_def(ks, ref).center
        tau = t
        if center == targ
            moving && (tau = t + s * (_norm3(p) / CLIGHT))
        elseif center != obs
            pc = _lt_series(ks, center, et, t, iref, corr, robs, vobs, ro)
            primal(_dot3(pc, pc)) > 0 && (tau = t + s * (_norm3(pc) / CLIGHT))
        end
        q = _m_mxv(frame_rotation(ks, Int32(1), ref, tau), q)
    end
    d = map(derivatives, q)
    return ntuple(k -> map(x -> x[k], d), Val(N + 1))
end

# --- public interface ---------------------------------------------------------------

function spkezr(ks::KernelSet, targ::NameOrID, et::Real, ref::NameOrID, corr::Aberration,
        obs::NameOrID)
    t, o = body_code(ks, targ), body_code(ks, obs)
    # As in CSPICE, a geometric state of a body relative to itself is zero before `ref`
    # and `et` are checked.
    corr.geometric && t == o && return zero(SVector{6,typeof(float(et))}), zero(float(et))
    e, f = check_epoch(et), frame_id(ks, ref)
    if corr.geometric
        r, v = spk_jet(ks, t, e, f, o, Val(1))
        return vcat(r, v), norm(r) / CLIGHT
    end
    return abcorr_state(ks, t, e, f, corr, o)
end
spkezr(ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID, corr::Aberration,
    obs::NameOrID) = spkezr(kernelset(ctx), targ, et, ref, corr, obs)
spkezr(targ::NameOrID, et::Real, ref::NameOrID, corr::Aberration, obs::NameOrID) =
    spkezr(default_context(), targ, et, ref, corr, obs)

function spkpos(ks::KernelSet, targ::NameOrID, et::Real, ref::NameOrID, corr::Aberration,
        obs::NameOrID)
    t, o = body_code(ks, targ), body_code(ks, obs)
    corr.geometric && t == o && return zero(Vec3{typeof(float(et))}), zero(float(et))
    e, f = check_epoch(et), frame_id(ks, ref)
    if corr.geometric
        (r,) = spk_jet(ks, t, e, f, o, Val(0))
        return r, norm(r) / CLIGHT
    end
    return abcorr_position(ks, t, e, f, corr, o)
end
spkpos(ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID, corr::Aberration,
    obs::NameOrID) = spkpos(kernelset(ctx), targ, et, ref, corr, obs)
spkpos(targ::NameOrID, et::Real, ref::NameOrID, corr::Aberration, obs::NameOrID) =
    spkpos(default_context(), targ, et, ref, corr, obs)

"""
    spkez([ctx,] targ::Integer, et, ref, abcorr, obs::Integer) -> (state, lt)

Return the state of body `targ` relative to body `obs` in frame `ref` at `et`, and the
one-way light time, as `spkez_c`.

This is [`spkezr`](@ref), which takes names or IDs, restricted to NAIF ID codes for the
bodies, as in CSPICE; the other arguments, the result and the errors are those of
`spkezr`.

See also [`spkezp`](@ref), [`spkgeo`](@ref).
"""
spkez(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg,
    obs::Integer) = spkezr(ks, targ, et, ref, abcorr, obs)
spkez(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg,
    obs::Integer) = spkezr(kernelset(ctx), targ, et, ref, abcorr, obs)
spkez(targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg, obs::Integer) =
    spkez(default_context(), targ, et, ref, abcorr, obs)

"""
    spkezp([ctx,] targ::Integer, et, ref, abcorr, obs::Integer) -> (pos, lt)

Return the position of body `targ` relative to body `obs` in frame `ref` at `et`, and
the one-way light time, as `spkezp_c`.

This is [`spkpos`](@ref), which takes names or IDs, restricted to NAIF ID codes for the
bodies, as in CSPICE; the other arguments, the result and the errors are those of
`spkpos`.

See also [`spkez`](@ref), [`spkgps`](@ref).
"""
spkezp(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg,
    obs::Integer) = spkpos(ks, targ, et, ref, abcorr, obs)
spkezp(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg,
    obs::Integer) = spkpos(kernelset(ctx), targ, et, ref, abcorr, obs)
spkezp(targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg, obs::Integer) =
    spkezp(default_context(), targ, et, ref, abcorr, obs)

"""
    spkder([ctx,] targ, et, ref, obs, Val(N)) -> NTuple{N+1,SVector{3}}
    spkder([ctx,] targ, et, ref, abcorr, obs, Val(N)) -> NTuple{N+1,SVector{3}}

Return the position of `targ` relative to `obs` in frame `ref` at `et` and its time
derivatives to order `N`. This call has no CSPICE counterpart.

Element `k + 1` of the result is the derivative of order `k`, in km/s^k. `et` is in
seconds past J2000 TDB. `targ`, `obs` and `ref` are names or integer IDs. The form
without `abcorr` is geometric, as is `abcorr = "NONE"`; `abcorr` is a string or an
[`Aberration`](@ref), as for [`spkezr`](@ref). With ForwardDiff loaded, `et` can be a
dual number. The call does not allocate under the condition that [`spkezr`](@ref)
states.

Order 0 is the position of [`spkezr`](@ref) and of [`spkpos`](@ref), to rounding.
Every higher order is the exact derivative of that position as a function of `et`, and
agrees with a finite difference of `spkpos`. The exception is an SPK segment that stores
the velocity apart from the position (for example types 3 and 9): there order 1 is the
stored velocity and the higher orders are its derivatives. A dual-number epoch gives
the same derivatives as this call.

Order 1 is the velocity of `spkezr` for a geometric state. With an `LT` correction it
differs from that velocity, by about 1e-8 relative for the Moon seen from the Earth,
because CSPICE approximates the rate of the light time.
[Exact derivatives and the rates of CSPICE](@ref) gives the function that is
differentiated.

Throws `SpiceError(:badarg)` if `N` is not an integer `>= 0`, `SpiceError(:invalidoption)`
if `abcorr` is not valid,
and the errors of [`spkezr`](@ref).

# Examples
```julia
et = str2et("2025-03-01T12:00:00")
r, v, a = spkder("MOON", et, "J2000", "EARTH", Val(2))          # km, km/s, km/s^2
r, v, a, j = spkder("MOON", et, "J2000", "LT+S", "EARTH", Val(3))
```

See also [`pxder`](@ref), [`spkezr`](@ref).
"""
function spkder(ks::KernelSet, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbcorrArg,
        obs::NameOrID, ::Val{N}) where {N}
    N isa Int && N >= 0 || throw(SpiceError(:badarg, "spkder: the order must be an integer >= 0"))
    return abcorr_jet(ks, body_code(ks, targ), check_epoch(et), frame_id(ks, ref),
        Aberration(abcorr), body_code(ks, obs), Val(N))
end
spkder(ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbcorrArg,
    obs::NameOrID, v::Val) = spkder(kernelset(ctx), targ, et, ref, abcorr, obs, v)
spkder(targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbcorrArg, obs::NameOrID, v::Val) =
    spkder(default_context(), targ, et, ref, abcorr, obs, v)

"""
    spkssb([ctx,] targ::Integer, et, ref) -> SVector{6}

Return the geometric state of body `targ` relative to the solar system barycenter in
frame `ref` at `et`, as `spkssb_c`.

Position in km, then velocity in km/s. `et` is in seconds past J2000 TDB. `targ` is a
NAIF ID code; `ref` is a frame name or a frame ID. The result equals the state of
`spkgeo(targ, et, ref, 0)`.

Throws the errors of [`spkgeo`](@ref).
"""
function spkssb(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID)
    r, v = spk_jet(ks, Int32(targ), check_epoch(et), frame_id(ks, ref), Int32(0), Val(1))
    return vcat(r, v)
end
spkssb(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID) =
    spkssb(kernelset(ctx), targ, et, ref)
spkssb(targ::Integer, et::Real, ref::NameOrID) = spkssb(default_context(), targ, et, ref)

@inline function _split_state(state::AbstractVector)
    length(state) == 6 || throw(SpiceError(:badarg, "a state vector must have 6 elements"))
    @inbounds return Vec3(state[1], state[2], state[3]), Vec3(state[4], state[5], state[6])
end

@inline function _vec3(a::AbstractVector)
    length(a) == 3 || throw(SpiceError(:badarg, "a vector must have 3 elements"))
    @inbounds return Vec3(a[1], a[2], a[3])
end

"""
    spkltc([ctx,] targ::Integer, et, ref, abcorr, stobs) -> (state, lt, dlt)

Return the state of body `targ` relative to an observer with a given barycentric
state, corrected for light time, with the light time and its rate, as `spkltc_c`.

- `et`: epoch of the observer in seconds past J2000 TDB.
- `stobs`: state of the observer relative to the solar system barycenter at `et`, in
  `ref`; 6 elements, km and km/s.
- `ref`: an inertial frame, by name or ID.
- `abcorr`: as for [`spkezr`](@ref). Stellar aberration is not applied, whatever
  `abcorr` says.

`state` is an `SVector{6}` in km and km/s, `lt` is in seconds and `dlt` is the time
derivative of `lt` (no unit). The `LT` forms make one iteration. The `CN` forms make
at most five and stop when the light time changes by less than `1e-17 * max(1, |epoch|)`.

Throws `SpiceError(:badarg)` if `ref` is not inertial, `stobs` does not have 6
elements, `et` is not finite, or the range rate of the target
is above `(1 - 1e-10)` times the speed of light. Throws `SpiceError(:invalidoption)` if
`abcorr` is not valid, and the SPK errors of [`spkgeo`](@ref).

See also [`spkaps`](@ref), [`spkacs`](@ref).
"""
function spkltc(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg,
        stobs::AbstractVector)
    robs, vobs = _split_state(stobs)
    q = _ltc(ks, Int32(targ), check_epoch(et), _need_inertial(ks, frame_id(ks, ref)),
        Aberration(abcorr), robs, vobs)
    return vcat(q.r, q.v), q.lt, q.dlt
end
spkltc(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg,
    stobs::AbstractVector) = spkltc(kernelset(ctx), targ, et, ref, abcorr, stobs)
spkltc(targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg, stobs::AbstractVector) =
    spkltc(default_context(), targ, et, ref, abcorr, stobs)

"""
    spkaps([ctx,] targ::Integer, et, ref, abcorr, stobs, accobs) -> (state, lt, dlt)

Return the state of body `targ` relative to an observer with a given barycentric state
and acceleration, corrected for light time and stellar aberration, as `spkaps_c`.

`stobs` (6 elements, km and km/s) and `accobs` (3 elements, km/s^2) are relative to
the solar system barycenter at `et`, in the inertial frame `ref`. The other arguments
are those of [`spkltc`](@ref). `abcorr` is as for
[`spkezr`](@ref). The acceleration is used for the `+S` forms only: it enters the
velocity part of the stellar aberration.

`state` is an `SVector{6}` in km and km/s, `lt` is in seconds and `dlt` is the time
derivative of `lt`.

Throws `SpiceError(:badarg)` if `accobs` does not have 3 elements, and the errors of
[`spkltc`](@ref).

See also [`spkacs`](@ref).
"""
function spkaps(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg,
        stobs::AbstractVector, accobs::AbstractVector)
    robs, vobs = _split_state(stobs)
    q = _aps(ks, Int32(targ), check_epoch(et), _need_inertial(ks, frame_id(ks, ref)),
        Aberration(abcorr), robs, vobs, _vec3(accobs))
    return vcat(q.r, q.v), q.lt, q.dlt
end
spkaps(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg,
    stobs::AbstractVector, accobs::AbstractVector) =
    spkaps(kernelset(ctx), targ, et, ref, abcorr, stobs, accobs)
spkaps(targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg, stobs::AbstractVector,
    accobs::AbstractVector) = spkaps(default_context(), targ, et, ref, abcorr, stobs, accobs)

"""
    spkacs([ctx,] targ::Integer, et, ref, abcorr, obs::Integer) -> (state, lt, dlt)

Return the state of body `targ` relative to body `obs` in the inertial frame `ref` at
`et` with the correction `abcorr`, the light time and its rate, as `spkacs_c`.

`state` is an `SVector{6}` in km and km/s, `lt` is in seconds and `dlt` is the time
derivative of `lt`. `targ` and `obs` are NAIF ID codes. For the `+S` forms the
acceleration of the observer is the central difference of its velocity over ±1 s, as
in CSPICE. `et` is in seconds past J2000 TDB. The state and `lt` are those of
[`spkezr`](@ref).

Throws the errors of [`spkltc`](@ref).
"""
function spkacs(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg,
        obs::Integer)
    corr = Aberration(abcorr)
    e = check_epoch(et)
    iref = _need_inertial(ks, frame_id(ks, ref))
    robs, vobs, accobs = _obs_state(ks, Int32(obs), e, iref, corr)
    q = _aps(ks, Int32(targ), e, iref, corr, robs, vobs, accobs)
    return vcat(q.r, q.v), q.lt, q.dlt
end
spkacs(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg,
    obs::Integer) = spkacs(kernelset(ctx), targ, et, ref, abcorr, obs)
spkacs(targ::Integer, et::Real, ref::NameOrID, abcorr::AbcorrArg, obs::Integer) =
    spkacs(default_context(), targ, et, ref, abcorr, obs)

"""
    spkapo([ctx,] targ::Integer, et, ref, sobs, abcorr) -> (pos, lt)

Return the position of body `targ` relative to an observer with a given barycentric
state, with the correction `abcorr`, and the light time, as `spkapo_c`.

`sobs` is the state of the observer relative to the solar system barycenter at `et`
(seconds past J2000 TDB) in the inertial frame `ref` (6 elements, km and km/s). `pos`
is an `SVector{3}` in km
and `lt` is in seconds. The `LT` forms make one light-time iteration and the `CN`
forms exactly three, as in CSPICE. Note the argument order: `sobs` is before `abcorr`.

Throws `SpiceError(:badarg)` if `ref` is not inertial, `sobs` does not have 6
elements, `et` is not finite, or (for `+S`) the speed of the
observer is not below the speed of light. Throws `SpiceError(:spkinvalidoption)` if
`abcorr` is not valid, and the SPK errors of [`spkgeo`](@ref).

See also [`spkpos`](@ref), [`spkltc`](@ref).
"""
function spkapo(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID, sobs::AbstractVector,
        abcorr::AbcorrArg)
    robs, vobs = _split_state(sobs)
    return _apo(ks, Int32(targ), check_epoch(et), _need_inertial(ks, frame_id(ks, ref)),
        _spk_abcorr(abcorr), robs, vobs)
end
spkapo(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID, sobs::AbstractVector,
    abcorr::AbcorrArg) = spkapo(kernelset(ctx), targ, et, ref, sobs, abcorr)
spkapo(targ::Integer, et::Real, ref::NameOrID, sobs::AbstractVector, abcorr::AbcorrArg) =
    spkapo(default_context(), targ, et, ref, sobs, abcorr)

"""
    ltime([ctx,] etobs, obs::Integer, dir, targ::Integer) -> (ettarg, elapsd)

Return the epoch at which a signal is at body `targ`, given that it is at body `obs`
at `etobs`, and the light time between the two, as `ltime_c`.

`dir = "->"`: the signal leaves `obs` at `etobs` and reaches `targ` at `ettarg`.
`dir = "<-"`: the signal reaches `obs` at `etobs` and left `targ` at `ettarg`. Epochs
are in seconds past J2000 TDB; `elapsd = |ettarg - etobs|` is in seconds. `obs` and
`targ` are NAIF ID codes. The call makes three iterations after the geometric start,
as CSPICE does.

Throws `SpiceError(:badarg)` if `dir` is not `"->"` or `"<-"` (blanks are not
ignored) or `etobs` is not finite, and the SPK errors of [`spkgeo`](@ref).
"""
function ltime(ks::KernelSet, etobs::Real, obs::Integer, dir::AbstractString, targ::Integer)
    s = dir == "->" ? 1 : dir == "<-" ? -1 :
        throw(SpiceError(:badarg, "ltime: the direction must be \"->\" or \"<-\", not \"$dir\""))
    et = check_epoch(etobs)
    robs = _ssbpos(ks, Int32(obs), et, Int32(1))
    elapsd = _norm3(_ssbpos(ks, Int32(targ), et, Int32(1)) - robs) / CLIGHT
    ettarg = et + s * elapsd
    # invariant: `ettarg = et + s * elapsd`, with `elapsd` from i - 1 iterations
    for _ in 1:3
        elapsd = _norm3(_ssbpos(ks, Int32(targ), ettarg, Int32(1)) - robs) / CLIGHT
        ettarg = et + s * elapsd
    end
    return ettarg, elapsd
end
ltime(ctx::SpiceContext, etobs::Real, obs::Integer, dir::AbstractString, targ::Integer) =
    ltime(kernelset(ctx), etobs, obs, dir, targ)
ltime(etobs::Real, obs::Integer, dir::AbstractString, targ::Integer) =
    ltime(default_context(), etobs, obs, dir, targ)

"""
    stelab(pobj, vobs) -> SVector{3}

Return the position `pobj` corrected for stellar aberration at reception, as
`stelab_c`.

`pobj` is the position of the target relative to the observer and `vobs` is the
velocity of the observer relative to the solar system barycenter, in the same
inertial frame (3 elements each; km and km/s). The result is `pobj` rotated toward
`vobs`; its norm is the norm of `pobj`.

Throws `SpiceError(:badarg)` if a vector does not have 3 elements or the speed
`|vobs|` is not below the speed of light.

# Examples
```jldoctest
julia> p = stelab([1.0e8, 0.0, 0.0], [0.0, 30.0, 0.0]);

julia> round(p[2] / p[1] * 1.0e4; digits=4)    # tangent of the angle is near v/c = 1.0007e-4
1.0007
```

See also [`stlabx`](@ref).
"""
stelab(pobj::AbstractVector, vobs::AbstractVector) = _stelab(_vec3(pobj), _vec3(vobs))

"""
    stlabx(pobj, vobs) -> SVector{3}

Return the position `pobj` corrected for stellar aberration at transmission, as
`stlabx_c`: [`stelab`](@ref) with the velocity negated.

Units and errors are those of [`stelab`](@ref).
"""
stlabx(pobj::AbstractVector, vobs::AbstractVector) = _stelab(_vec3(pobj), -_vec3(vobs))
