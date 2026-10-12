# Local solar time (`et2lst`) and planetocentric longitude of the Sun (`lspcn`).

"Quotient and remainder of `num / den` with a remainder in `[0, den)`, as `rmaind_`."
function _lst_rmaind(num::Float64, den::Float64)
    q = trunc(num / den)
    r = num - q * den
    if r < 0.0
        q -= 1.0
        r += den
    end
    return q, r
end

"""
    _lst_sense(ks, body) -> Int

Contract: the sign that turns a planetographic longitude of `body` into a
planetocentric one, as `pgrrec_`: `BODY<body>_PGR_POSITIVE_LON` (`EAST` gives 1, `WEST`
gives -1) if it exists; otherwise 1 for a retrograde body and -1 for a prograde body
(sign of the prime meridian rate), and 1 for the Earth, the Moon and the Sun. Throws
`SpiceError(:invalidoption)` for another value of the variable and
`SpiceError(:missingdata)` if the rate is needed and `BODY<body>_PM` is not available.
"""
function _lst_sense(ks::KernelSet, body::Integer)
    key = "BODY$(body)_PGR_POSITIVE_LON"
    v = get(ks.pool, key, nothing)
    if v isa Vector{String} && !isempty(v)
        s = uppercase(replace(v[1], ' ' => ""))
        s == "EAST" && return 1
        s == "WEST" && return -1
        throw(SpiceError(:invalidoption,
            "kernel variable $key may have the values EAST or WEST; the value was $(v[1])"))
    end
    body == 399 && return 1
    pm = get(ks.pool, "BODY$(body)_PM", nothing)
    (pm isa Vector{Float64} && length(pm) >= 2) || throw(SpiceError(:missingdata,
        "the prime meridian rate of kernel variable BODY$(body)_PM is required but not " *
        "available for body $body"))
    (body == 301 || body == 10) && return 1
    return pm[2] < 0.0 ? 1 : -1
end

"Two decimal digits of the integer `n` in 0:99."
_lst_two(n::Integer) = string(n ÷ 10, n % 10)

"""
    et2lst(et, body, lon, type) -> (hr, mn, sc, time, ampm)

Local solar time at ET `et` for the longitude `lon` (radians) on the body with ID
`body`, as `et2lst_c`: hour, minute and second on a 24 hour clock, the string
`HH:MM:SS`, and the string `HH:MM:SS A.M.` or `HH:MM:SS P.M.` on a 12 hour clock.
`type` is `"PLANETOCENTRIC"` or `"PLANETOGRAPHIC"` (case ignored) and says how `lon`
is measured. The local solar time is the angle between the meridian of `lon` and the
meridian of the Sun (position with the correction `LT+S` in the body-fixed frame of
`cidfrm`), 12:00:00 at the sub-solar meridian, with 24 "hours" in one rotation; the
sign of `BODY<body>_PM[2]` gives the direction of the rotation. For the Sun the
result is always 12:00:00.

Throws `SpiceError(:unknownsystem)` for another `type`, `SpiceError(:cantfindframe)`
if the body has no body-fixed frame, and `SpiceError(:cantgetrotationtype)` without a
numeric `BODY<body>_PM` of two or more values. Forms: on a `KernelSet`, on a
`SpiceContext`, and on the default context.
"""
function et2lst(ks::KernelSet, et::Real, body::Integer, lon::Real, type::AbstractString)
    system = uppercase(strip(type))
    if system == "PLANETOGRAPHIC"
        g = _lst_sense(ks, body) * Float64(lon)
        x, y = cos(g), sin(g)
        mylong = (x == 0.0 && y == 0.0) ? 0.0 : atan(y, x)
    elseif system == "PLANETOCENTRIC"
        mylong = Float64(lon)
    else
        throw(SpiceError(:unknownsystem, "the coordinate system '$type' is not a " *
            "recognized system of longitude; the recognized systems are " *
            "'PLANETOCENTRIC' and 'PLANETOGRAPHIC'"))
    end
    body == 10 && return 12, 0, 0, "12:00:00", "12:00:00 P.M."
    fix = cidfrm(ks, body)
    fix === nothing && throw(SpiceError(:cantfindframe, "the body-fixed frame " *
        "associated with body $body could not be determined"))
    state, _ = spkezr(ks, 10, et, fix[1], Aberration("LT+S"), body)
    slong = (state[1] == 0.0 && state[2] == 0.0) ? 0.0 : atan(state[2], state[1])
    _, angle = _lst_rmaind(mylong - slong, _TWOPI)
    angle > pi && (angle -= _TWOPI)
    pm = get(ks.pool, "BODY$(body)_PM", nothing)
    (pm isa Vector{Float64} && length(pm) >= 2) || throw(SpiceError(:cantgetrotationtype,
        "the rotation type for the body $body could not be determined because the " *
        "BODY$(body)_PM keyword was not found in the pool or was not of the expected " *
        "type or dimension"))
    pm[2] < 0.0 && (angle = -angle)
    secnds = clamp(angle * 86400.0 / _TWOPI, -43200.0, 43200.0)
    hours, rest = _lst_rmaind(secnds, 3600.0)
    mins, secnds = _lst_rmaind(rest, 60.0)
    hr, mn, sc = trunc(Int, hours) + 12, trunc(Int, mins), trunc(Int, secnds)
    hr == 24 && (hr = 0)
    h12 = hr == 0 ? 12 : hr > 12 ? hr - 12 : hr
    suffix = hr >= 12 ? "P.M." : "A.M."
    tail = string(":", _lst_two(mn), ":", _lst_two(sc))
    return hr, mn, sc, _lst_two(hr) * tail, string(_lst_two(h12), tail, " ", suffix)
end
et2lst(ctx::SpiceContext, et::Real, body::Integer, lon::Real, type::AbstractString) =
    et2lst(kernelset(ctx), et, body, lon, type)
et2lst(et::Real, body::Integer, lon::Real, type::AbstractString) =
    et2lst(default_context(), et, body, lon, type)

"""
    lspcn(body, et, abcorr) -> Float64

Planetocentric longitude of the Sun (L_s, radians in `[0, 2pi)`) for the body `body`
(name or ID) at ET `et`, as `lspcn_c`: the right ascension of the Sun as seen from the
body, with the correction `abcorr`, in the frame whose Z axis is the orbital angular
momentum of the body about the Sun and whose X axis is the cross product of the north
pole of the body and that Z axis (the vernal equinox of the body). Needs PCK data for
the pole and SPK data for the body and the Sun. Forms: on a `KernelSet`, on a
`SpiceContext`, and on the default context.
"""
function lspcn(ks::KernelSet, body::NameOrID, et::Real, abcorr::AbcorrArg)
    id = body_code(ks, body)
    tipm = tipbod(ks, 1, id, et)
    npole = Vec3{Float64}(tipm[3, 1], tipm[3, 2], tipm[3, 3])
    r, v = spk_jet(ks, id, check_epoch(et), _J2000_FRAME, Int32(10), Val(1))
    z = _dyn_ucross(r, v)
    x = _dyn_ucross(npole, z)
    primal(_dot3(x, x)) == 0 && throw(SpiceError(:dependentvectors,
        "lspcn: the pole of body $id is parallel to its orbital angular momentum"))
    y = _dyn_ucross(z, x)
    s, _ = spkezr(ks, 10, et, 1, Aberration(abcorr), id)
    p = Vec3{Float64}(s[1], s[2], s[3])
    px, py = dot(x, p), dot(y, p)
    ra = (px == 0.0 && py == 0.0) ? 0.0 : atan(py, px)
    return ra < 0.0 ? ra + _TWOPI : ra
end
lspcn(ctx::SpiceContext, body::NameOrID, et::Real, abcorr::AbcorrArg) =
    lspcn(kernelset(ctx), body, et, abcorr)
lspcn(body::NameOrID, et::Real, abcorr::AbcorrArg) = lspcn(default_context(), body, et, abcorr)
