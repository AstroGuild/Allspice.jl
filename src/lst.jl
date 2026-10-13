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
    et2lst([ctx,] et, body, lon, type) -> (hr, mn, sc, time, ampm)

Return the local solar time at `et` for the longitude `lon` on the body with ID
`body`, as `et2lst_c`.

- `et`: seconds past J2000 TDB.
- `lon`: longitude in radians.
- `type`: `"PLANETOCENTRIC"` or `"PLANETOGRAPHIC"`, in any letter case; it says how
  `lon` is measured.

`hr`, `mn` and `sc` are the integer hour, minute and second on a 24 hour clock. `time`
is the string `HH:MM:SS`, and `ampm` is the string `HH:MM:SS A.M.` or `HH:MM:SS P.M.`
on a 12 hour clock.

The local solar time is the angle between the meridian of `lon` and the meridian of
the Sun, with 12:00:00 at the sub-solar meridian and 24 "hours" in one rotation; an
"hour" is not 3600 s. The position of the Sun has the correction `"LT+S"` and is in
the body-fixed frame of [`cidfrm`](@ref). The sign of the second value of
`BODY<body>_PM` gives the direction of the rotation. For the Sun (`body = 10`) the
result is always 12:00:00. The call needs SPK data for the body and the Sun and PCK
data for the body.

Throws `SpiceError`:
- `:unknownsystem` if `type` is not one of the two;
- `:cantfindframe` if the body has no body-fixed frame;
- `:cantgetrotationtype` if `BODY<body>_PM` is absent or is not a numeric variable of
  two or more values;
- `:invalidoption`, for `"PLANETOGRAPHIC"`, if `BODY<body>_PGR_POSITIVE_LON` is not
  `EAST` or `WEST`;
- `:missingdata`, for `"PLANETOGRAPHIC"`, if `BODY<body>_PGR_POSITIVE_LON` and
  `BODY<body>_PM` are both absent;
- the kinds of [`spkezr`](@ref).

# Examples
```julia
furnsh("naif0012.tls"); furnsh("pck00011.tpc"); furnsh("de440s.bsp")
et = str2et("2025-03-01T12:00:00")
hr, mn, sc, time, ampm = et2lst(et, 399, deg2rad(-118.0), "PLANETOCENTRIC")
```
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
    lspcn([ctx,] body, et, abcorr) -> Float64

Return the planetocentric longitude of the Sun (`Ls`) for `body` at `et`, in radians
in `[0, 2pi)`, as `lspcn_c`.

`Ls` is the season angle of the body: 0 at its spring equinox of the north, `pi/2` at
the summer solstice of the north. It is the right ascension of the Sun as seen from
the body, in this frame:

- the Z axis is the orbital angular momentum of the body about the Sun;
- the X axis is the cross product of the north pole of the body and the Z axis.

`body` is a name or a NAIF ID, `et` is in seconds past J2000 TDB, and `abcorr` is as
for [`spkezr`](@ref). The call needs PCK data for the pole and SPK data for the body
and the Sun.

Throws `SpiceError(:dependentvectors)` if the pole is parallel to the orbital angular
momentum, and the errors of [`spkezr`](@ref) and [`tipbod`](@ref).

# Examples
```julia
furnsh("naif0012.tls"); furnsh("pck00011.tpc"); furnsh("de440s.bsp")
ls = rad2deg(lspcn("EARTH", str2et("2025-03-01"), "LT+S"))    # degrees
```
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
