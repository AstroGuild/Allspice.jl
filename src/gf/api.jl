# Geometry finder, public interface: `gfdist`, `gfsep`, `gfposc`, `gfrr`, `gfpa`, `gfuds`,
# `gfudb`, `gfevnt` and `gfstol`, as the calls of CSPICE N0067 with these differences
# (the searches that need surface, occultation or field-of-view geometry are in
# `surface.jl`, `illumination.jl`, `occultation.jl` and `fov.jl`, with the same rules):
#
# - The result is a new `Window`; there is no result cell, no workspace and no interval
#   count (`nintvls`), because a `Window` has no size limit.
# - CSPICE keeps the step (`gfsstp`), the tolerance (`gfstol`) and the progress reporter
#   in saved variables. Here the step is an argument, the tolerance is the keyword
#   `tol`, and the callbacks are the keywords `bail` and `report`. `gfstol` sets the
#   default of `tol` for the searches on one `SpiceContext`; it is an atomic field of
#   the context, as the defaults of `timdef` are, and `kclear` does not change it. A
#   search on a `KernelSet` has the default of CSPICE, $(_GF_CNVTOL) s.
# - A search reads the kernel snapshot once, at its start, and is a pure function of
#   that snapshot and of its arguments, so searches can run on many threads at once and
#   `furnsh` on another thread does not disturb them.
#
# Every call has a form on a `KernelSet`, on a `SpiceContext` and on the default context.

"""
    gfstol([ctx,] value) -> Float64
    gfstol([ctx]) -> Float64

Set the default convergence tolerance of the geometry finder searches on the context
`ctx` to `value` seconds and return it, as `gfstol_c`. Without `value`, return the
current default.

The default of a new context is $(_GF_CNVTOL) s. The keyword `tol` of a search overrides
it for that search. A search reads the default once, at its start, so a change does not
reach a search that is running. `kclear`, `furnsh` and `unload` do not change it.

Throws `SpiceError(:invalidtolerance)` if `value <= 0`; the default is then unchanged.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> gfstol(ctx)
1.0e-6

julia> gfstol(ctx, 1e-3)
0.001

julia> gfstol(ctx)
0.001
```

See also [`gfdist`](@ref).
"""
function gfstol(ctx::SpiceContext, value::Real)
    value > 0 || throw(SpiceError(:invalidtolerance,
        "Convergence tolerance must be greater-than zero. Input VALUE = $value."))
    @atomic ctx.gftol = Float64(value)
    return Float64(value)
end
gfstol(value::Real) = gfstol(default_context(), value)
gfstol(ctx::SpiceContext) = @atomic ctx.gftol
gfstol() = gfstol(default_context())

"""
    gfdist([ctx,] target, abcorr, obsrvr, relate, refval, adjust, step, cnfine;
           tol=gfstol(ctx), bail=nothing, report=nothing, threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which the distance in km
between `target` and `obsrvr` satisfies the relation `relate` with the reference value
`refval`, as `gfdist_c`.

The result is a new `Window` that is a subset of the confinement window `cnfine`.
`cnfine` is not modified, and an empty `cnfine` gives an empty window. CSPICE needs a
result cell, a workspace and an interval count; this call needs none.

# Arguments
- `target`, `obsrvr`: names or IDs of two different bodies.
- `abcorr`: any aberration correction of [`spkezr`](@ref), a string or an `Aberration`.
- `relate`: `"<"`, `"="`, `">"`, `"LOCMIN"`, `"LOCMAX"`, `"ABSMIN"` or `"ABSMAX"`, in any
  case. An inequality gives the intervals on which it holds. `"="` and the extrema with
  `adjust = 0` give intervals of zero length, one for each event. A local extremum is
  reported only if it is strictly inside an interval of `cnfine`; an absolute extremum
  can be at an end.
- `refval`: reference value in the unit of the quantity. The extrema do not use it.
- `adjust`: zero, except for `"ABSMIN"` and `"ABSMAX"`. There a positive value gives the
  window on which the quantity is within `adjust` of the absolute extremum.
- `step`: step of the search in seconds, positive. The search samples the direction of
  change of the quantity at this step to find its local extrema, so `step` must be
  shorter than the shortest interval, inside `cnfine`, on which the quantity is
  monotone. A longer step can miss events. A shorter step costs time and gives no
  accuracy: `tol` sets the accuracy.
- `cnfine`: the `Window` to search in. For `"="` the quantity is also evaluated up to
  $(_GF_ADDWIN) s outside each interval, so the kernels must cover that margin.

# Keywords
Every search takes these keywords; [`gfevnt`](@ref), [`gfocce`](@ref) and
[`gffove`](@ref) take `tol` as an argument instead.
- `tol`: convergence tolerance in seconds, positive. Each end point of the result is
  bisected until its bracket is not longer than `tol`. The default is the value of
  [`gfstol`](@ref) for the context.
- `bail`: `nothing`, or a function without arguments that the search calls in its
  loops. When it returns `true` the search stops with an `InterruptException`.
- `report`: `nothing`, or a progress callback (`gfrepi`, `gfrepu`, `gfrepf` of CSPICE).
  The search calls `report(:init, window, prefix, suffix)` at the start of a pass over
  `window`, `report(:update, ivbeg, ivend, et)` while it searches [`ivbeg`, `ivend`],
  and `report(:final)` at the end of the pass. `prefix` and `suffix` are the two strings
  of the progress message of CSPICE, such as `"Distance pass 1 of 2"` and `"done."`.
- `threads`: if `true`, each interval that is searched runs in its own task. The result
  does not change. The callbacks, and the functions of [`gfuds`](@ref) and
  [`gfudb`](@ref), are then called from several threads.

A search uses the kernels that are loaded when it starts. A `furnsh` or an `unload` on
another thread does not change a search that runs, and several searches can run at the
same time.

Throws `SpiceError` with the kind
- `:invalidstep` if `step <= 0`, `:invalidtolerance` if `tol <= 0`;
- `:notrecognized` if `relate` is not a relation;
- `:valueoutofrange` if `adjust < 0`, `:invalidvalue` if `adjust` is not zero with a
  relation other than `"ABSMIN"` and `"ABSMAX"`;
- `:idcodenotfound` if a body name is not known, `:bodiesnotdistinct` if the two bodies
  are the same;
- `:invalidoption` if `abcorr` is not a correction.

The errors of [`spkezr`](@ref), such as missing ephemeris data, come from the epochs
that the search samples.

# Examples
The example needs a leap seconds kernel and an SPK with the Earth and the Moon.
```julia
furnsh("naif0012.tls")
furnsh("de440s.bsp")
cnfine = Window([str2et("2025 JAN 01"), str2et("2026 JAN 01")])

# Perigees of the Moon: local minima of the distance, searched with a step of one day.
perigees = gfdist("MOON", "NONE", "EARTH", "LOCMIN", 0.0, 0.0, 86400.0, cnfine)
for i in 1:wncard(perigees)
    et, _ = wnfetd(perigees, i)          # an event is an interval of zero length
    println(et2utc(et, "C", 0))
end

# Intervals on which the Moon is nearer than 370000 km, and their total length.
near = gfdist("MOON", "NONE", "EARTH", "<", 370000.0, 0.0, 86400.0, cnfine)
wnsumd(near)[1]                          # seconds
```

See also [`gfrr`](@ref), [`gfsep`](@ref), [`gfpa`](@ref), [`gfposc`](@ref),
[`gfuds`](@ref), [`gfstol`](@ref), [`Window`](@ref).
"""
function gfdist(ks::KernelSet, target::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
        relate::AbstractString, refval::Real, adjust::Real, step::Real, cnfine::Window;
        tol::Real=_GF_CNVTOL, bail=nothing, report=nothing, threads::Bool=false)
    run = _gf_run(step, tol; bail, report, threads)
    _gf_relcode(relate)
    q = GFDistance(ks, target, abcorr, obsrvr)
    return _gf_search(q, relate, refval, adjust, cnfine, run)
end

"""
    gfrr([ctx,] target, abcorr, obsrvr, relate, refval, adjust, step, cnfine;
         tol=gfstol(ctx), bail=nothing, report=nothing, threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which the range rate in km/s
of `target` relative to `obsrvr` satisfies the relation `relate` with `refval`, as
`gfrr_c`.

The arguments, the keywords, the result and the errors are those of [`gfdist`](@ref).
The search also reads the state 1 s before and 1 s after each epoch that it samples,
as CSPICE does, so the kernels must cover `cnfine` with that margin.

See also [`gfdist`](@ref).
"""
function gfrr(ks::KernelSet, target::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
        relate::AbstractString, refval::Real, adjust::Real, step::Real, cnfine::Window;
        tol::Real=_GF_CNVTOL, bail=nothing, report=nothing, threads::Bool=false)
    run = _gf_run(step, tol; bail, report, threads)
    _gf_relcode(relate)
    q = GFRangeRate(ks, target, abcorr, obsrvr)
    return _gf_search(q, relate, refval, adjust, cnfine, run)
end

for f in (:gfdist, :gfrr)
    @eval begin
        $f(ctx::SpiceContext, target::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
            relate::AbstractString, refval::Real, adjust::Real, step::Real, cnfine::Window;
            tol::Real=gfstol(ctx), kw...) =
            $f(kernelset(ctx), target, abcorr, obsrvr, relate, refval, adjust, step, cnfine;
                tol, kw...)
        $f(target::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, relate::AbstractString,
            refval::Real, adjust::Real, step::Real, cnfine::Window; kw...) =
            $f(default_context(), target, abcorr, obsrvr, relate, refval, adjust, step,
                cnfine; kw...)
    end
end

"""
    gfsep([ctx,] targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, relate,
          refval, adjust, step, cnfine; tol=gfstol(ctx), bail=nothing, report=nothing,
          threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which the angular separation
in radians of `targ1` and `targ2`, as seen from `obsrvr`, satisfies the relation
`relate` with `refval`, as `gfsep_c`.

Each shape is `"POINT"` or `"SPHERE"`, in any case. For a sphere the separation is
measured from the limb, so it is negative when the disks overlap. The radius of the
sphere is the largest of the radii of the body (`BODY<id>_RADII`, km). `frame1` and
`frame2` are not read, as in CSPICE; pass `" "`. The three bodies must differ. The
other arguments, the keywords and the result are those of [`gfdist`](@ref).

Throws the errors of [`gfdist`](@ref), and `SpiceError` with the kind
- `:notrecognized` if a shape is not one of the two;
- the kinds of [`srfrec`](@ref) if the radii of a sphere are missing or not valid;
- `:insidebody` or `:badgeometry` if the observer is not outside a sphere at an epoch
  of the search.

See also [`gfdist`](@ref), [`gfoclt`](@ref).
"""
function gfsep(ks::KernelSet, targ1::NameOrID, shape1::AbstractString, frame1::NameOrID,
        targ2::NameOrID, shape2::AbstractString, frame2::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, relate::AbstractString, refval::Real, adjust::Real, step::Real,
        cnfine::Window; tol::Real=_GF_CNVTOL, bail=nothing, report=nothing,
        threads::Bool=false)
    run = _gf_run(step, tol; bail, report, threads)
    _gf_relcode(relate)
    q = GFSeparation(ks, targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr)
    return _gf_search(q, relate, refval, adjust, cnfine, run)
end
gfsep(ctx::SpiceContext, targ1::NameOrID, shape1::AbstractString, frame1::NameOrID,
    targ2::NameOrID, shape2::AbstractString, frame2::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, relate::AbstractString, refval::Real, adjust::Real, step::Real,
    cnfine::Window; tol::Real=gfstol(ctx), kw...) =
    gfsep(kernelset(ctx), targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr,
        relate, refval, adjust, step, cnfine; tol, kw...)
gfsep(targ1::NameOrID, shape1::AbstractString, frame1::NameOrID, targ2::NameOrID,
    shape2::AbstractString, frame2::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
    relate::AbstractString, refval::Real, adjust::Real, step::Real, cnfine::Window; kw...) =
    gfsep(default_context(), targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr,
        relate, refval, adjust, step, cnfine; kw...)

"""
    gfpa([ctx,] target, illmn, abcorr, obsrvr, relate, refval, adjust, step, cnfine;
         tol=gfstol(ctx), bail=nothing, report=nothing, threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which the phase angle in
radians satisfies the relation `relate` with `refval`, as `gfpa_c`. The phase angle is
the angle at `target` between the directions to the illuminator `illmn` and to the
observer `obsrvr`.

`abcorr` is `"NONE"` or a reception correction. The three bodies must differ. The other
arguments, the keywords and the result are those of [`gfdist`](@ref).

Throws the errors of [`gfdist`](@ref), and `SpiceError(:invalidoption)` if `abcorr` is a
transmission correction.

See also [`gfdist`](@ref), [`gfilum`](@ref), [`phaseq`](@ref).
"""
function gfpa(ks::KernelSet, target::NameOrID, illmn::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, relate::AbstractString, refval::Real, adjust::Real, step::Real,
        cnfine::Window; tol::Real=_GF_CNVTOL, bail=nothing, report=nothing,
        threads::Bool=false)
    run = _gf_run(step, tol; bail, report, threads)
    _gf_relcode(relate)
    q = GFPhase(ks, target, illmn, abcorr, obsrvr)
    return _gf_search(q, relate, refval, adjust, cnfine, run)
end
gfpa(ctx::SpiceContext, target::NameOrID, illmn::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, relate::AbstractString, refval::Real, adjust::Real, step::Real,
    cnfine::Window; tol::Real=gfstol(ctx), kw...) =
    gfpa(kernelset(ctx), target, illmn, abcorr, obsrvr, relate, refval, adjust, step,
        cnfine; tol, kw...)
gfpa(target::NameOrID, illmn::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
    relate::AbstractString, refval::Real, adjust::Real, step::Real, cnfine::Window; kw...) =
    gfpa(default_context(), target, illmn, abcorr, obsrvr, relate, refval, adjust, step,
        cnfine; kw...)

"""
    _gf_position_coordinate(ks, target, frame, abcorr, obsrvr, crdsys, coord) -> GFCoordinate

The coordinate `coord` of the system `crdsys` of the position of `target` relative to
`obsrvr` in `frame`, with the checks of `zzgfcoin` in their order: bodies
(`:idcodenotfound`, `:bodiesnotdistinct`), correction, system and coordinate
(`:notsupported`), frame (`:unknownframe`; `:noframe` for the systems that need the
radii of its center).
"""
function _gf_position_coordinate(ks::KernelSet, target::NameOrID, frame::NameOrID,
        abcorr::AbcorrArg, obsrvr::NameOrID, crdsys::AbstractString, coord::AbstractString)
    targ = _gf_body(ks, target, "target object")
    obs = _gf_body(ks, obsrvr, "observer")
    targ == obs && throw(SpiceError(:bodiesnotdistinct,
        "The observer and target must be distinct objects, but are not: " *
        "OBSRVR = $obsrvr; TARGET = $target."))
    corr = Aberration(abcorr)
    ref = Int32(0)                       # 0: the frame is not known
    try
        ref = _gf_frame(ks, frame isa AbstractString ? String(strip(frame, ' ')) : frame)
    catch e
        (e isa SpiceError && e.kind === :unknownframe) || rethrow()
    end
    q = GFCoordinate(GFPositionVector(ks, targ, obs, corr, ref), ks, crdsys, coord)
    ref == 0 && throw(SpiceError(:unknownframe, "the frame name '$frame' is not known"))
    return q
end

"""
    gfposc([ctx,] target, frame, abcorr, obsrvr, crdsys, coord, relate, refval, adjust,
           step, cnfine; tol=gfstol(ctx), bail=nothing, report=nothing,
           threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which the coordinate `coord`
of the position of `target` relative to `obsrvr` in the frame `frame` satisfies the
relation `relate` with `refval`, as `gfposc_c`.

The coordinate systems and their coordinates, in any case, are

| `crdsys` | `coord` |
|---|---|
| `"RECTANGULAR"` | `"X"`, `"Y"`, `"Z"` |
| `"LATITUDINAL"` | `"RADIUS"`, `"LONGITUDE"`, `"LATITUDE"` |
| `"RA/DEC"` | `"RANGE"`, `"RIGHT ASCENSION"`, `"DECLINATION"` |
| `"SPHERICAL"` | `"RADIUS"`, `"COLATITUDE"`, `"LONGITUDE"` |
| `"CYLINDRICAL"` | `"RADIUS"`, `"LONGITUDE"`, `"Z"` |
| `"GEODETIC"` | `"LONGITUDE"`, `"LATITUDE"`, `"ALTITUDE"` |
| `"PLANETOGRAPHIC"` | `"LONGITUDE"`, `"LATITUDE"`, `"ALTITUDE"` |

Angles, and `refval` and `adjust` for an angle, are in radians; lengths are in km.
Longitude is in (-pi, pi] in the latitudinal, spherical and geodetic systems. Longitude
in the cylindrical and planetographic systems, and right ascension, are in [0, 2 pi).
The geodetic and planetographic systems use the radii of the center of `frame`, whose
two equatorial radii must be equal.

A search on a longitude or on a right ascension accounts for the jump at the branch
cut: the jump gives no root and no extremum. Such a search converges to
`max(1e-7, tol / 10)` seconds, not to `tol`.

The other arguments, the keywords and the result are those of [`gfdist`](@ref).

Throws the errors of [`gfdist`](@ref), and `SpiceError` with the kind
- `:notsupported` if the system or the coordinate is not in the table, or if the
  equatorial radii differ;
- `:unknownframe` if `frame` is not known (`:noframe` with the geodetic and
  planetographic systems);
- `:badframe` if, for these two systems, the center of `frame` has no radii loaded.

# Examples
The example needs a leap seconds kernel and an SPK with the Earth and the Sun.
```julia
furnsh("naif0012.tls")
furnsh("de440s.bsp")
cnfine = Window([str2et("2025 JAN 01"), str2et("2026 JAN 01")])

# Epochs at which the Sun crosses the equator of J2000, with a step of 10 days.
w = gfposc("SUN", "J2000", "LT+S", "EARTH", "RA/DEC", "DECLINATION", "=", 0.0, 0.0,
    10 * 86400.0, cnfine)
[et2utc(wnfetd(w, i)[1], "C", 0) for i in 1:wncard(w)]

# Intervals on which the right ascension of the Sun is below 6 h (pi/2 rad).
gfposc("SUN", "J2000", "LT+S", "EARTH", "RA/DEC", "RIGHT ASCENSION", "<", pi / 2, 0.0,
    10 * 86400.0, cnfine)
```

See also [`gfdist`](@ref), [`gfsubc`](@ref), [`gfsntc`](@ref).
"""
function gfposc(ks::KernelSet, target::NameOrID, frame::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, crdsys::AbstractString, coord::AbstractString,
        relate::AbstractString, refval::Real, adjust::Real, step::Real, cnfine::Window;
        tol::Real=_GF_CNVTOL, bail=nothing, report=nothing, threads::Bool=false)
    run = _gf_run(step, tol; bail, report, threads)
    _gf_relcode(relate)
    q = _gf_position_coordinate(ks, target, frame, abcorr, obsrvr, crdsys, coord)
    return _gf_search(q, relate, refval, adjust, cnfine, run)
end
gfposc(ctx::SpiceContext, target::NameOrID, frame::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, crdsys::AbstractString, coord::AbstractString,
    relate::AbstractString, refval::Real, adjust::Real, step::Real, cnfine::Window;
    tol::Real=gfstol(ctx), kw...) =
    gfposc(kernelset(ctx), target, frame, abcorr, obsrvr, crdsys, coord, relate, refval,
        adjust, step, cnfine; tol, kw...)
gfposc(target::NameOrID, frame::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
    crdsys::AbstractString, coord::AbstractString, relate::AbstractString, refval::Real,
    adjust::Real, step::Real, cnfine::Window; kw...) =
    gfposc(default_context(), target, frame, abcorr, obsrvr, crdsys, coord, relate, refval,
        adjust, step, cnfine; kw...)

"""
    gfuds([ctx,] udfuns, udqdec, relate, refval, adjust, step, cnfine;
          tol=gfstol(ctx), bail=nothing, report=nothing, threads=false) -> Window

Return the window of `cnfine` on which the scalar function `udfuns(et)` satisfies the
relation `relate` with `refval`, as `gfuds_c`.

`udfuns(et)` returns a real number. `udqdec(et)` returns `true` if `udfuns` decreases at
`et`, and must return a `Bool`. For a function without an analytical derivative, use
`et -> uddc(udfuns, et, dx)`. Each function takes the epoch alone. It must return the
same value each time for the same epoch. CSPICE also passes `udfuns` to `udqdec`; here
a closure does that. `step` must be shorter than the shortest interval, inside
`cnfine`, on which `udfuns` is monotone. The other arguments, the keywords and the
result are those of [`gfdist`](@ref).

The search reads no kernel, and `ctx` only gives the default of `tol`. A kernel call
in `udfuns` reads the kernels that are loaded at the time of that call, not those of
the start of the search.

Throws `SpiceError(:invalidstep)` if `step <= 0`, `SpiceError(:invalidtolerance)` if
`tol <= 0`, `SpiceError(:notrecognized)` if `relate` is not a relation,
`SpiceError(:valueoutofrange)` if `adjust < 0`, and `SpiceError(:invalidvalue)` if
`adjust` is not zero with a relation other than `"ABSMIN"` and `"ABSMAX"`.

# Examples
```jldoctest
julia> cnfine = Window([0.0, 10.0]);

julia> w = gfuds(sin, et -> cos(et) < 0, ">", 0.5, 0.0, 1.0, cnfine);

julia> [round.(wnfetd(w, i); digits=4) for i in 1:wncard(w)]
2-element Vector{Tuple{Float64, Float64}}:
 (0.5236, 2.618)
 (6.8068, 8.9012)

julia> peaks = gfuds(sin, et -> uddc(sin, et, 1e-3), "LOCMAX", 0.0, 0.0, 1.0, cnfine);

julia> [round(wnfetd(peaks, i)[1]; digits=4) for i in 1:wncard(peaks)]
2-element Vector{Float64}:
 1.5708
 7.854
```

See also [`gfudb`](@ref), [`uddc`](@ref), [`uddf`](@ref), [`gfdist`](@ref).
"""
function gfuds(::KernelSet, udfuns::F, udqdec::D, relate::AbstractString, refval::Real,
        adjust::Real, step::Real, cnfine::Window; tol::Real=_GF_CNVTOL, bail=nothing,
        report=nothing, threads::Bool=false) where {F,D}
    run = _gf_run(step, tol; bail, report, threads)
    return _gf_search(GFFunction(udfuns, udqdec), relate, refval, adjust, cnfine, run)
end
gfuds(ctx::SpiceContext, udfuns::F, udqdec::D, relate::AbstractString, refval::Real,
    adjust::Real, step::Real, cnfine::Window; tol::Real=gfstol(ctx), kw...) where {F,D} =
    gfuds(kernelset(ctx), udfuns, udqdec, relate, refval, adjust, step, cnfine; tol, kw...)
gfuds(udfuns::F, udqdec::D, relate::AbstractString, refval::Real, adjust::Real,
    step::Real, cnfine::Window; kw...) where {F,D} =
    gfuds(default_context(), udfuns, udqdec, relate, refval, adjust, step, cnfine; kw...)

"""
    gfudb([ctx,] udfunb, step, cnfine; tol=gfstol(ctx), bail=nothing, report=nothing,
          threads=false) -> Window

Return the window of `cnfine` on which the boolean function `udfunb(et)` is `true`, as
`gfudb_c`.

`udfunb` takes the epoch alone and must return a `Bool`; the scalar function that CSPICE
passes through to it is not needed. The search samples `udfunb` from the start of each
interval of `cnfine` in steps of `step` seconds, and bisects each change of value until
the bracket is not longer than `tol`. `step` must therefore be shorter than the shortest
interval on which `udfunb` keeps one value. The keywords are those of [`gfdist`](@ref),
and `ctx` acts as in [`gfuds`](@ref).

Throws `SpiceError(:invalidstep)` if `step <= 0` and `SpiceError(:invalidtolerance)` if
`tol <= 0`.

# Examples
```jldoctest
julia> w = gfudb(et -> sin(et) > 0.5, 1.0, Window([0.0, 10.0]));

julia> [round.(wnfetd(w, i); digits=4) for i in 1:wncard(w)]
2-element Vector{Tuple{Float64, Float64}}:
 (0.5236, 2.618)
 (6.8068, 8.9012)
```

See also [`gfuds`](@ref), [`gfoclt`](@ref), [`gftfov`](@ref).
"""
function gfudb(::KernelSet, udfunb::F, step::Real, cnfine::Window; tol::Real=_GF_CNVTOL,
        bail=nothing, report=nothing, threads::Bool=false) where {F}
    run = _gf_run(step, tol; bail, report, threads)
    _gf_repi(run, cnfine, "User defined boolean event search ", "done.")
    result = _gf_scan!(Window(), et -> udfunb(et)::Bool, cnfine, run, false)
    _gf_repf(run)
    return result
end
gfudb(ctx::SpiceContext, udfunb::F, step::Real, cnfine::Window; tol::Real=gfstol(ctx),
    kw...) where {F} = gfudb(kernelset(ctx), udfunb, step, cnfine; tol, kw...)
gfudb(udfunb::F, step::Real, cnfine::Window; kw...) where {F} =
    gfudb(default_context(), udfunb, step, cnfine; kw...)

# --- gfevnt ------------------------------------------------------------------------------

function _gf_qdpars(qdpars::AbstractVector{<:Real}, name::String)
    length(qdpars) >= 3 || throw(SpiceError(:invalidcount, "The parameter $name needs " *
        "3 numeric values in QDPARS; $(length(qdpars)) were supplied."))
    return Vec3{Float64}(qdpars[begin], qdpars[begin+1], qdpars[begin+2])
end

"Names of the quantities of `gfevnt` and, for each, the parameters that must be supplied."
const _GF_QPARS = (
    "ANGULAR SEPARATION" => ("TARGET1", "FRAME1", "SHAPE1", "TARGET2", "FRAME2", "SHAPE2",
        "OBSERVER", "ABCORR"),
    "DISTANCE" => ("TARGET", "OBSERVER", "ABCORR"),
    "COORDINATE" => ("TARGET", "OBSERVER", "ABCORR", "COORDINATE SYSTEM", "COORDINATE",
        "REFERENCE FRAME", "VECTOR DEFINITION", "METHOD", "DVEC", "DREF"),
    "RANGE RATE" => ("TARGET", "OBSERVER", "ABCORR"),
    "PHASE ANGLE" => ("TARGET", "OBSERVER", "ILLUM", "ABCORR"),
    "ILLUMINATION ANGLE" => ("TARGET", "ILLUM", "OBSERVER", "ABCORR", "REFERENCE FRAME",
        "ANGTYP", "METHOD", "SPOINT"),
)

"""
    gfevnt([ctx,] udstep, udrefn, gquant, qpnams, qcpars, qdpars, qipars, qlpars, op,
           refval, tol, adjust, cnfine; bail=nothing, report=nothing, threads=false) -> Window

Return the window of `cnfine` on which the quantity `gquant` satisfies the relation `op`
with the reference value `refval`, as `gfevnt_c`. This is the general entry of the
geometry finder, with a step and a refinement function that the caller supplies.

# Arguments
- `udstep`: the step in seconds, a positive number or a function `et -> step`. A step
  function that returns a value that is not positive throws `SpiceError(:invalidstep)`
  during the search.
- `udrefn`: the refinement function `(t1, t2, s1, s2) -> t`. It returns an epoch in the
  bracket [`t1`, `t2`]. `s1` and `s2` are the values (`Bool`) of the condition of the
  search at `t1` and `t2`. [`gfrefn`](@ref) is the bisection of CSPICE.
- `gquant`: `"DISTANCE"`, `"ANGULAR SEPARATION"`, `"COORDINATE"`, `"RANGE RATE"`,
  `"PHASE ANGLE"` or `"ILLUMINATION ANGLE"`, in any case.
- `qpnams`, `qcpars`: the names of the parameters of the quantity and their string
  values, in the same order. The parameters that each quantity needs are in the table
  below. A value is the string that the search of the quantity takes for that argument,
  for example `ILLUM` is the `illmn` of [`gfpa`](@ref).
- `qdpars`: the three numbers of the vector `DVEC` of a surface intercept, or of
  `SPOINT` of an illumination angle. The name `DVEC` or `SPOINT` is in `qpnams` with a
  blank string value. Another quantity reads no number: pass `Float64[]`.
- `qipars`, `qlpars`: not read; they keep the argument order of CSPICE. Pass `Int[]`
  and `Bool[]`.
- `op`, `refval`, `adjust`, `cnfine`: as `relate`, `refval`, `adjust` and `cnfine` of
  [`gfdist`](@ref).
- `tol`: convergence tolerance in seconds, positive. It is an argument, as in CSPICE;
  the default of [`gfstol`](@ref) is not used.

| `gquant` | Parameters in `qpnams` |
|---|---|
| `"DISTANCE"`, `"RANGE RATE"` | `TARGET`, `OBSERVER`, `ABCORR` |
| `"PHASE ANGLE"` | `TARGET`, `OBSERVER`, `ILLUM`, `ABCORR` |
| `"ANGULAR SEPARATION"` | `TARGET1`, `FRAME1`, `SHAPE1`, `TARGET2`, `FRAME2`, `SHAPE2`, `OBSERVER`, `ABCORR` |
| `"COORDINATE"` | `TARGET`, `OBSERVER`, `ABCORR`, `COORDINATE SYSTEM`, `COORDINATE`, `REFERENCE FRAME`, `VECTOR DEFINITION`, `METHOD`, `DVEC`, `DREF` |
| `"ILLUMINATION ANGLE"` | `TARGET`, `ILLUM`, `OBSERVER`, `ABCORR`, `REFERENCE FRAME`, `ANGTYP`, `METHOD`, `SPOINT` |

`VECTOR DEFINITION` is `POSITION`, `SUB-OBSERVER POINT` or `SURFACE INTERCEPT POINT`.
Every parameter in the row of the quantity must be in `qpnams`, also one that the
quantity does not use (a blank value is then enough). The keywords `bail`, `report` and
`threads` are those of [`gfdist`](@ref).

Throws `SpiceError` with the kind
- `:notrecognized` if `gquant` or `op` is not known;
- `:invalidcount` if there are more than 10 parameter names, fewer values in `qcpars`
  than names, or fewer than 3 values in `qdpars` where a vector is needed;
- `:missingvalue` if a parameter of the quantity is not in `qpnams`;
- `:notsupported` if the vector definition is not known.

Also throws the errors of the search for the quantity: [`gfdist`](@ref),
[`gfrr`](@ref), [`gfpa`](@ref), [`gfsep`](@ref), [`gfposc`](@ref), [`gfsubc`](@ref),
[`gfsntc`](@ref), [`gfilum`](@ref).

# Examples
The example needs a leap seconds kernel and an SPK with the Earth and the Moon, loaded
with [`furnsh`](@ref).
```julia
cnfine = Window([str2et("2025 JAN 01"), str2et("2025 MAR 01")])

# The search of `gfdist("MOON", "LT+S", "EARTH", "<", 370000.0, 0.0, 86400.0, cnfine)`,
# with a tolerance of 1e-3 s.
gfevnt(86400.0, gfrefn, "DISTANCE", ["TARGET", "OBSERVER", "ABCORR"],
    ["MOON", "EARTH", "LT+S"], Float64[], Int[], Bool[], "<", 370000.0, 1e-3, 0.0, cnfine)
```

See also [`gfocce`](@ref), [`gffove`](@ref).
"""
function gfevnt(ks::KernelSet, udstep, udrefn, gquant::AbstractString,
        qpnams::AbstractVector{<:AbstractString}, qcpars::AbstractVector{<:AbstractString},
        qdpars::AbstractVector{<:Real}, qipars::AbstractVector{<:Integer},
        qlpars::AbstractVector{Bool}, op::AbstractString, refval::Real, tol::Real,
        adjust::Real, cnfine::Window; bail=nothing, report=nothing, threads::Bool=false)
    quant = uppercase(strip(gquant, ' '))
    k = findfirst(p -> p.first == quant, _GF_QPARS)
    k === nothing && throw(SpiceError(:notrecognized, "The geometric quantity, $gquant " *
        "is not recognized. Supported quantities are: DISTANCE, PHASE ANGLE, " *
        "COORDINATE, RANGE RATE, ANGULAR SEPARATION,ILLUMINATION ANGLE."))
    n = length(qpnams)
    (n <= 10 && length(qcpars) >= n) || throw(SpiceError(:invalidcount,
        "Number of quantity definition parameters = $n;  must be in range 0:10."))
    pnames = String[uppercase(strip(s, ' ')) for s in qpnams]
    for p in _GF_QPARS[k].second         # invariant: the parameters before p are supplied
        p in pnames || throw(SpiceError(:missingvalue, "The parameter $p is required in " *
            "order to compute events pertaining to the quantity $quant; this parameter " *
            "was not supplied."))
    end
    par(name) = String(strip(qcpars[findfirst(==(name), pnames)], ' '))
    _gf_relcode(op)
    run = _gf_run(udstep, tol; refine=udrefn, bail, report, threads)
    if quant == "DISTANCE"
        q = GFDistance(ks, par("TARGET"), par("ABCORR"), par("OBSERVER"))
        return _gf_search(q, op, refval, adjust, cnfine, run)
    elseif quant == "RANGE RATE"
        q = GFRangeRate(ks, par("TARGET"), par("ABCORR"), par("OBSERVER"))
        return _gf_search(q, op, refval, adjust, cnfine, run)
    elseif quant == "PHASE ANGLE"
        q = GFPhase(ks, par("TARGET"), par("ILLUM"), par("ABCORR"), par("OBSERVER"))
        return _gf_search(q, op, refval, adjust, cnfine, run)
    elseif quant == "ANGULAR SEPARATION"
        q = GFSeparation(ks, par("TARGET1"), par("SHAPE1"), par("FRAME1"), par("TARGET2"),
            par("SHAPE2"), par("FRAME2"), par("ABCORR"), par("OBSERVER"))
        return _gf_search(q, op, refval, adjust, cnfine, run)
    elseif quant == "COORDINATE"
        rawdef = qcpars[findfirst(==("VECTOR DEFINITION"), pnames)]
        vecdef = _gf_compress(rawdef)
        if vecdef == "POSITION"
            q = _gf_position_coordinate(ks, par("TARGET"), par("REFERENCE FRAME"),
                par("ABCORR"), par("OBSERVER"), par("COORDINATE SYSTEM"), par("COORDINATE"))
            return _gf_search(q, op, refval, adjust, cnfine, run)
        end
        intercept = vecdef == "SURFACE INTERCEPT POINT"
        (intercept || vecdef == "SUB-OBSERVER POINT") || throw(SpiceError(:notsupported,
            "The vector definition $vecdef is not supported."))
        # CSPICE compares the definition as the caller wrote it with the upper-case name
        # before it checks the center of the frame, so another spelling skips that check
        # (the geometry call then makes it at the first evaluation); this is kept.
        chkctr = rstrip(rawdef, ' ') == vecdef
        dvec = intercept ? _gf_qdpars(qdpars, "DVEC") : _GF_ZERO3
        q = _gf_surface_coordinate(ks, intercept, chkctr, par("TARGET"),
            par("REFERENCE FRAME"), par("METHOD"), par("ABCORR"), par("OBSERVER"),
            par("DREF"), dvec, par("COORDINATE SYSTEM"), par("COORDINATE"))
        return _gf_search(q, op, refval, adjust, cnfine, run)
    end
    q = GFIllumination(ks, par("METHOD"), par("ANGTYP"), par("TARGET"), par("ILLUM"),
        par("REFERENCE FRAME"), par("ABCORR"), par("OBSERVER"), _gf_qdpars(qdpars, "SPOINT"))
    return _gf_search(q, op, refval, adjust, cnfine, run)
end
gfevnt(ctx::SpiceContext, udstep, udrefn, gquant::AbstractString,
    qpnams::AbstractVector{<:AbstractString}, qcpars::AbstractVector{<:AbstractString},
    qdpars::AbstractVector{<:Real}, qipars::AbstractVector{<:Integer},
    qlpars::AbstractVector{Bool}, op::AbstractString, refval::Real, tol::Real,
    adjust::Real, cnfine::Window; kw...) =
    gfevnt(kernelset(ctx), udstep, udrefn, gquant, qpnams, qcpars, qdpars, qipars, qlpars,
        op, refval, tol, adjust, cnfine; kw...)
gfevnt(udstep, udrefn, gquant::AbstractString, qpnams::AbstractVector{<:AbstractString},
    qcpars::AbstractVector{<:AbstractString}, qdpars::AbstractVector{<:Real},
    qipars::AbstractVector{<:Integer}, qlpars::AbstractVector{Bool}, op::AbstractString,
    refval::Real, tol::Real, adjust::Real, cnfine::Window; kw...) =
    gfevnt(default_context(), udstep, udrefn, gquant, qpnams, qcpars, qdpars, qipars,
        qlpars, op, refval, tol, adjust, cnfine; kw...)
