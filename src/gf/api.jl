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
    gfstol([ctx], value) -> Float64
    gfstol([ctx]) -> Float64

Set the default convergence tolerance, in seconds, of the geometry-finder searches on
the context `ctx` (the default context if it is omitted), as `gfstol_c`; the form
without `value` returns the current default. The default of a new context is
$(_GF_CNVTOL) s.

Contract: precondition `value > 0`, else `SpiceError(:invalidtolerance)` and the default
is unchanged. The default is an atomic field of the context: a search reads it once at
its start, so a change never reaches a search that is running, and `kclear`, `furnsh`
and `unload` do not change it. The keyword `tol` of a search overrides it for that
search. Postcondition: `gfstol(ctx) == value`.
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
    gfdist(target, abcorr, obsrvr, relate, refval, adjust, step, cnfine; keywords) -> Window

Window of the confinement window `cnfine` on which the distance (km) between `target`
and `obsrvr` satisfies the relation `relate` with the reference value `refval`, as
`gfdist_c`. The bodies are names or IDs; `abcorr` is any correction of `spkezr`, a
string or an `Aberration`.

`relate` is `"<"`, `"="`, `">"`, `"LOCMIN"`, `"LOCMAX"`, `"ABSMIN"` or `"ABSMAX"`, in any
case. `adjust` is zero but for an absolute extremum, where a positive value gives the
window on which the quantity is within `adjust` of the extremum. `step` (seconds) is
the step of the search for the local extrema of the quantity: it must be shorter than
the shortest interval, inside `cnfine`, on which the quantity is monotone. `cnfine` is
not modified.

Keywords, the same for every search of the geometry finder:

- `tol`: convergence tolerance in seconds. The default is the value of `gfstol` for the
  context ($(_GF_CNVTOL) s unless changed, and for the form on a `KernelSet`).
- `bail`: `nothing` or a function without arguments that the search calls in its loops;
  when it returns `true` the search stops with an `InterruptException`.
- `report`: `nothing` or a progress callback, called as `report(:init, window, prefix,
  suffix)` at the start of a pass over `window`, `report(:update, ivbeg, ivend, et)`
  while the interval [`ivbeg`, `ivend`] is searched, and `report(:final)` at the end of
  the pass (`gfrepi`, `gfrepu`, `gfrepf`).
- `threads`: if `true`, the intervals of a window that is searched (the intervals of
  `cnfine` in the first pass, the intervals on which the quantity is monotone in the
  second) each run in a task. The result does not change. The callbacks, and the
  functions of `gfuds` and `gfudb`, are then called from several threads.

Contract: the search is the search of CSPICE (`zzgfrelx`, `zzgfsolvx`): the same
samples, the same brackets and the same bisection, so an end point differs from the one
of CSPICE by rounding of the states only, and by no more than about `tol`. The result
is a new valid window that is a subset of `cnfine`; an empty `cnfine` gives an empty
window. Errors: `:invalidstep` (`step <= 0`), `:invalidtolerance`, `:notrecognized`
(relation), `:idcodenotfound`, `:bodiesnotdistinct`, `:badarg` (correction),
`:valueoutofrange` (`adjust < 0`), `:invalidvalue` (`adjust` with another relation),
and the errors of `spkezr`.
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
    gfrr(target, abcorr, obsrvr, relate, refval, adjust, step, cnfine; keywords) -> Window

Window of `cnfine` on which the range rate (km/s) of `target` relative to `obsrvr`
satisfies the relation, as `gfrr_c`. The arguments and the keywords are those of
`gfdist`.

Contract: as `gfdist`. The range rate is the one of the state of `spkezr` in J2000; the
test "is decreasing" takes the acceleration from the states one second before and after
the epoch, as CSPICE does.
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
    gfsep(targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, relate, refval,
          adjust, step, cnfine; keywords) -> Window

Window of `cnfine` on which the angular separation (radians) of `targ1` and `targ2`
seen from `obsrvr` satisfies the relation, as `gfsep_c`. Each shape is `"POINT"` or
`"SPHERE"`; for a sphere the separation is measured from the limb, with the largest
radius of the body (`BODY<id>_RADII`). The frames are not used by these shapes, as in
CSPICE. The other arguments and the keywords are those of `gfdist`.

Contract: as `gfdist`, and `:notrecognized` for another shape, `:insidebody` or
`:badgeometry` if the observer is inside a sphere, and the errors of `bodvcd` for a
sphere without radii.
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
    gfpa(target, illmn, abcorr, obsrvr, relate, refval, adjust, step, cnfine; keywords) -> Window

Window of `cnfine` on which the phase angle (radians) at `target`, between the
illuminator `illmn` and the observer `obsrvr`, satisfies the relation, as `gfpa_c`. The
other arguments and the keywords are those of `gfdist`.

Contract: as `gfdist`; `abcorr` is a reception correction, and a transmission
correction throws `SpiceError(:invalidoption)`.
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
    gfposc(target, frame, abcorr, obsrvr, crdsys, coord, relate, refval, adjust, step,
           cnfine; keywords) -> Window

Window of `cnfine` on which the coordinate `coord` of the position of `target` relative
to `obsrvr` in `frame` satisfies the relation, as `gfposc_c`. The systems and their
coordinates are

| `crdsys` | `coord` |
|---|---|
| `RECTANGULAR` | `X`, `Y`, `Z` |
| `LATITUDINAL` | `RADIUS`, `LONGITUDE`, `LATITUDE` |
| `RA/DEC` | `RANGE`, `RIGHT ASCENSION`, `DECLINATION` |
| `SPHERICAL` | `RADIUS`, `COLATITUDE`, `LONGITUDE` |
| `CYLINDRICAL` | `RADIUS`, `LONGITUDE`, `Z` |
| `GEODETIC` | `LONGITUDE`, `LATITUDE`, `ALTITUDE` |
| `PLANETOGRAPHIC` | `LONGITUDE`, `LATITUDE`, `ALTITUDE` |

in any case; angles are in radians and lengths in km. Longitude is in (-pi, pi] in the
latitudinal, spherical and geodetic systems and in [0, 2 pi) in the cylindrical and
planetographic systems, as right ascension is. The geodetic and planetographic systems
use the radii of the center of `frame`. The other arguments and the keywords are those
of `gfdist`.

Contract: as `gfdist`. A longitude or a right ascension is searched with the proxy
functions of `zzgflong`, so that the jump at the branch cut gives no false root and no
false extremum; these searches use the tolerance `max(1e-7, tol / 10)`. Errors as
`gfdist`, and `:notsupported` (system, coordinate, or unequal equatorial radii),
`:unknownframe`, `:noframe`, `:badframe` (no radii for the center of the frame).
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
    gfuds(udfuns, udqdec, relate, refval, adjust, step, cnfine; keywords) -> Window

Window of `cnfine` on which the scalar function `udfuns(et)` satisfies the relation, as
`gfuds_c`. `udqdec(et)` returns `true` if the function decreases at `et`; for a
function without an analytical derivative it can be `et -> uddc(udfuns, et, dx)`. The
functions are Julia functions or closures of the epoch alone (CSPICE passes `udfuns` to
`udqdec`; a closure makes this needless). The other arguments and the keywords are
those of `gfdist`.

Forms: `gfuds(udfuns, ...)` and `gfuds(ctx, udfuns, ...)` take the default of `tol` from
the default context and from `ctx`; `gfuds(ks::KernelSet, udfuns, ...)` has the default
$(_GF_CNVTOL) s. No kernel is read.

Contract: as `gfdist`, without the errors of the bodies; `udfuns` returns a real number
and `udqdec` a `Bool`, and each is a function of the epoch alone (the search may skip an
evaluation that CSPICE repeats at the same epoch).
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
    gfudb(udfunb, step, cnfine; keywords) -> Window

Window of `cnfine` on which the boolean function `udfunb(et)` is true, as `gfudb_c`.
`step` (seconds) must be shorter than the shortest interval on which the function keeps
one value. `udfunb` is a Julia function or closure of the epoch alone (the scalar
function `udfuns` that CSPICE passes through to `udfunb` is not needed). The keywords
are those of `gfdist`, and the forms those of `gfuds`.

Contract: the search is `zzgfsolvx` on each interval of `cnfine`: the state is sampled
from the start of the interval in steps of `step`, and a change of state is bisected
until the bracket is not longer than `tol`. The result is a new valid window that is a
subset of `cnfine`. Errors: `:invalidstep`, `:invalidtolerance`.
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

"The first three numeric parameters of `gfevnt`, the vector `name` (`SpiceError(:invalidcount)` if there are fewer)."
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
    gfevnt(udstep, udrefn, gquant, qpnams, qcpars, qdpars, qipars, qlpars, op, refval, tol,
           adjust, cnfine; bail=nothing, report=nothing, threads=false) -> Window

General entry of the geometry finder, as `gfevnt_c`: window of `cnfine` on which the
quantity `gquant` satisfies the relation `op` with the reference value `refval`.

- `udstep`: the step in seconds, a positive number or a function `et -> step`.
- `udrefn`: the refinement function `(t1, t2, s1, s2) -> t`, which gives an epoch in the
  bracket [`t1`, `t2`] whose ends have the states `s1` and `s2`; `gfrefn` is the
  bisection of CSPICE.
- `gquant`: `"DISTANCE"`, `"ANGULAR SEPARATION"`, `"COORDINATE"`, `"RANGE RATE"`,
  `"PHASE ANGLE"` or `"ILLUMINATION ANGLE"`, in any case.
- `qpnams`, `qcpars`: names of the parameters of the quantity and their character
  values, in the same order. The names are those of CSPICE: `TARGET`, `OBSERVER`,
  `ABCORR` (distance, range rate); with `ILLUM` (phase angle); `TARGET1`, `FRAME1`,
  `SHAPE1`, `TARGET2`, `FRAME2`, `SHAPE2`, `OBSERVER`, `ABCORR` (angular separation);
  `TARGET`, `OBSERVER`, `ABCORR`, `COORDINATE SYSTEM`, `COORDINATE`, `REFERENCE FRAME`,
  `VECTOR DEFINITION`, `METHOD`, `DVEC`, `DREF` (coordinate; the vector definition is
  `POSITION`, `SUB-OBSERVER POINT` or `SURFACE INTERCEPT POINT`); `TARGET`, `ILLUM`,
  `OBSERVER`, `ABCORR`, `REFERENCE FRAME`, `ANGTYP`, `METHOD`, `SPOINT` (illumination
  angle).
- `qdpars`, `qipars`, `qlpars`: numeric, integer and logical parameter values. Only the
  first three values of `qdpars` are used, for the vector `DVEC` of a surface intercept
  or `SPOINT` of an illumination angle; the others are in the list for the order of
  CSPICE.
- `tol`: convergence tolerance in seconds, an argument here as in CSPICE.

The keywords are those of `gfdist`.

Contract: as `gfdist`, with the refinement of the caller in the place of the bisection.
Errors: `:notrecognized` (quantity or relation), `:invalidcount` (more than 10
parameters, or fewer than 3 numeric values where a vector is needed), `:missingvalue`
(a parameter of the quantity is not supplied), `:notsupported` (vector definition),
and the errors of the search of the quantity ([`gfsubc`](@ref), [`gfsntc`](@ref),
[`gfilum`](@ref) for the surface quantities).
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
