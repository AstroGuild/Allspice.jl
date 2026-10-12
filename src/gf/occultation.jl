# Geometry finder, occultations and transits: `gfoclt` and `gfocce` of CSPICE N0067,
# with the initializer `zzgfocin` and the state `zzgfocst` of `zzgfocu`.
#
# An occultation search has no scalar quantity: it is the search for the intervals on
# which a binary state is true, so it is one call of the step-and-bracket search
# (`zzgfsolv`) for each interval of the confinement window, as `gfudb` is. The state at
# an epoch is the classification of `occult` for the roles "front" and "back"
# (`_geo_occ_eval` of `geometry/occultation.jl`).

"Names of the occultation types, and the code of `_geo_occ_eval` that each one accepts (0: any code but 0)."
const _GF_OCCTYPES = ("ANNULAR" => 2, "ANY" => 0, "PARTIAL" => 1, "FULL" => 3)

"""
    _GFOccultation

The checked inputs of an occultation search, the saved values of `zzgfocin`: the two
targets with their shapes and frames, the surface list for the bounding sphere of the
back target, the observer, the correction without stellar aberration, the code `want`
of the occultation type and the DSK tolerances.
"""
struct _GFOccultation
    ks::KernelSet
    tols::DSKTolerances
    front::_GeoOccBody
    fframe::Int32
    back::_GeoOccBody
    bframe::Int32
    bsrf::Vector{Int}
    obs::Int32
    corr::Aberration
    want::Int
end

"""
    _gf_occultation(ks, tols, occtyp, front, fshape, fframe, back, bshape, bframe,
                    abcorr, obsrvr) -> _GFOccultation

The checks of `gfocce` and `zzgfocin`, in their order: the two shapes are not both
`POINT` (`:invalidshapecombo`); the three bodies are known (`:idcodenotfound`) and
distinct (`:bodiesnotdistinct`); each shape is `POINT`, `ELLIPSOID` or a DSK method
string (`:invalidshape`, and the errors of the method string); a DSK shape needs a
point as the other target (`:invalidshapecombo`); `occtyp` is `ANY`, `FULL`, `ANNULAR`
or `PARTIAL` (`:invalidocctype`), and `ANY` if a target is a point
(`:badtypeshapecombo`); `abcorr` is a correction of `spkezr` (`:invalidoption`); an
extended target has radii or DSK data, and a body-fixed frame that is known and
centered on it (`:invalidframe`).
"""
function _gf_occultation(ks::KernelSet, tols::DSKTolerances, occtyp::AbstractString,
        front::NameOrID, fshape::AbstractString, fframe::NameOrID, back::NameOrID,
        bshape::AbstractString, bframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID)
    (_geo_occ_upper(fshape) == "POINT" && _geo_occ_upper(bshape) == "POINT") &&
        throw(SpiceError(:invalidshapecombo, "The front and back target shape " *
            "specifications are both POINT; at least one of these targets must be an " *
            "extended object."))
    idf = _gf_body(ks, front, "front target object")
    idb = _gf_body(ks, back, "back target object")
    obs = _gf_body(ks, obsrvr, "observer")
    (idf == idb || idf == obs || idb == obs) && throw(SpiceError(:bodiesnotdistinct,
        "The observer and both targets must be distinct objects, but are not: " *
        "OBSRVR = $obsrvr; FRONT = $front; BACK = $back."))
    pf = _geo_occ_parse(ks, idf, fshape, "front target")
    pf.kind == 0x00 && _geo_occ_bad_shape(fshape)
    pb = _geo_occ_parse(ks, idb, bshape, "back target")
    pb.kind == 0x00 && _geo_occ_bad_shape(bshape)
    fb = _GeoOccBody(idf, pf.kind, pf.m, Int32(0))
    bb = _GeoOccBody(idb, pb.kind, pb.m, Int32(0))
    # CSPICE reads the surface list of the back target with the ID of the front target.
    bsrf = (bb.kind == _GEO_OCC_DSK && !isempty(bb.m.srflst)) ?
           _geo_method(ks, idf, _geo_occ_upper(bshape)).srflst : bb.m.srflst
    ((fb.kind == _GEO_OCC_DSK && bb.kind != _GEO_OCC_POINT) ||
     (bb.kind == _GEO_OCC_DSK && fb.kind != _GEO_OCC_POINT)) && throw(SpiceError(
        :invalidshapecombo, "Front target shape from FSHAPE string was <$fshape>; back " *
        "target shape from BSHAPE was <$bshape>. When one shape is DSK, the other must " *
        "be POINT."))
    typ = uppercase(lstrip(rstrip(occtyp, ' '), ' '))
    k = findfirst(p -> p.first == typ, _GF_OCCTYPES)
    k === nothing && throw(SpiceError(:invalidocctype, "The occultation type $occtyp is " *
        "not recognized.  Supported types are: ANNULAR, ANY, PARTIAL,  FULL."))
    ((fb.kind == _GEO_OCC_POINT || bb.kind == _GEO_OCC_POINT) && typ != "ANY") &&
        throw(SpiceError(:badtypeshapecombo, "Occultation type $occtyp is not allowed " *
            "when either target body is modeled as a point. Set OCCTYP to ANY for use " *
            "with point targets."))
    full = _geo_abcorr(abcorr)
    corr = full.uselt ? Aberration(; lt=true, converged=full.converged, transmit=full.transmit) :
           full
    # Radii, frames and bounding spheres, the front target first; the frames are kept
    # as IDs for the evaluations.
    ffrm = _geo_occ_extent(ks, fb, fframe, fb.m.srflst)[4]
    bfrm = _geo_occ_extent(ks, bb, bframe, bsrf)[4]
    return _GFOccultation(ks, tols, fb, ffrm, bb, bfrm, bsrf, obs, corr,
        _GF_OCCTYPES[k].second)
end

"""
The state of `zzgfocst`: true if at `et` the front target occults the back target in
the way that the occultation type asks for.
"""
@inline function _gf_occulted(o::_GFOccultation, et::Float64)
    code = _geo_occ_eval(o.ks, o.tols, o.front, o.fframe, o.back, o.bframe, o.bsrf, o.obs,
        o.corr, et)
    return o.want == 0 ? code != 0 : code == o.want
end

"""
    _gf_oclt(ks, tols, occtyp, front, fshape, fframe, back, bshape, bframe, abcorr,
             obsrvr, cnfine, run) -> Window

`gfocce` after the checks of the step and of the tolerance: the checks of
`_gf_occultation`, then the search of each interval of `cnfine` for the intervals on
which the state `_gf_occulted` is true. Postcondition: a new window, a subset of
`cnfine`.
"""
function _gf_oclt(ks::KernelSet, tols::DSKTolerances, occtyp::AbstractString,
        front::NameOrID, fshape::AbstractString, fframe::NameOrID, back::NameOrID,
        bshape::AbstractString, bframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
        cnfine::Window, run::_GFRun)
    o = _gf_occultation(ks, tols, occtyp, front, fshape, fframe, back, bshape, bframe,
        abcorr, obsrvr)
    _gf_repi(run, cnfine, "Occultation/transit search ", "done.")
    result = _gf_scan!(Window(), et -> _gf_occulted(o, et), cnfine, run, false)
    _gf_repf(run)
    return result
end

"""
    gfoclt(occtyp, front, fshape, fframe, back, bshape, bframe, abcorr, obsrvr, step,
           cnfine; keywords) -> Window

Window of the confinement window `cnfine` on which the body `front` occults the body
`back` as seen from `obsrvr` (or, which is the same event, `front` is in transit across
`back`), as `gfoclt_c`.

- `occtyp`: `"FULL"` (the back target is hidden), `"ANNULAR"` (the front target lies
  wholly within the disc of the back target), `"PARTIAL"` (the discs overlap in part) or
  `"ANY"` (one of the three), in any case. With a point target only `"ANY"` is valid.
- `fshape`, `bshape`: `"ELLIPSOID"`, `"POINT"` or `"DSK/UNPRIORITIZED[/SURFACES =
  <list>]"`. One target at least is extended, and a target with DSK data needs a point
  as the other target.
- `fframe`, `bframe`: the body-fixed frames of the extended targets (a blank string for
  a point).
- `abcorr`: a correction of `spkezr`; stellar aberration is not used (`"LT+S"` acts as
  `"LT"`), as in CSPICE.
- `step`: the step of the search in seconds. It must be shorter than the shortest
  occultation of the kind asked for, and than the shortest gap between two of them,
  inside `cnfine`.

The keywords are those of [`gfdist`](@ref); the forms are on a `KernelSet`, on a
`SpiceContext` (which gives the default of `tol` and the DSK tolerances) and on the
default context.

Contract: the search is the search of CSPICE (`zzgfsolv` on each interval of `cnfine`
with the state of `zzgfocst`): the same samples and the same bisection, so an end point
differs from the one of CSPICE by rounding of the states only, and by no more than
about `tol`. The result is a new valid window that is a subset of `cnfine`; `cnfine` is
not modified. Errors: `:invalidstep`, `:invalidtolerance`, `:invalidshapecombo`,
`:idcodenotfound`, `:bodiesnotdistinct`, `:invalidshape`, `:invalidocctype`,
`:badtypeshapecombo`, `:invalidoption` (correction), `:invalidframe`, the errors of the
radii lookup and of the method string, `:notdisjoint` (the observer is inside a target,
or the targets overlap), and the errors of `spkpos` and of [`sincpt`](@ref).
"""
function gfoclt(ks::KernelSet, occtyp::AbstractString, front::NameOrID,
        fshape::AbstractString, fframe::NameOrID, back::NameOrID, bshape::AbstractString,
        bframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window;
        tol::Real=_GF_CNVTOL, bail=nothing, report=nothing, threads::Bool=false)
    return _gf_oclt(ks, _DSK_DEFAULT_TOLS, occtyp, front, fshape, fframe, back, bshape,
        bframe, abcorr, obsrvr, cnfine, _gf_run(step, tol; bail, report, threads))
end
gfoclt(ctx::SpiceContext, occtyp::AbstractString, front::NameOrID, fshape::AbstractString,
    fframe::NameOrID, back::NameOrID, bshape::AbstractString, bframe::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window; tol::Real=gfstol(ctx),
    bail=nothing, report=nothing, threads::Bool=false) =
    _gf_oclt(kernelset(ctx), dsktolerances(ctx), occtyp, front, fshape, fframe, back,
        bshape, bframe, abcorr, obsrvr, cnfine, _gf_run(step, tol; bail, report, threads))
gfoclt(occtyp::AbstractString, front::NameOrID, fshape::AbstractString, fframe::NameOrID,
    back::NameOrID, bshape::AbstractString, bframe::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, step::Real, cnfine::Window; kw...) =
    gfoclt(default_context(), occtyp, front, fshape, fframe, back, bshape, bframe, abcorr,
        obsrvr, step, cnfine; kw...)

"""
    gfocce(occtyp, front, fshape, fframe, back, bshape, bframe, abcorr, obsrvr, tol,
           udstep, udrefn, cnfine; bail=nothing, report=nothing, threads=false) -> Window

The occultation search of [`gfoclt`](@ref) with the step and the refinement of the
caller, as `gfocce_c`: `udstep` is the step in seconds, a positive number or a function
`et -> step`; `udrefn` is the refinement function `(t1, t2, s1, s2) -> t`, which gives
an epoch in the bracket [`t1`, `t2`] whose ends have the states `s1` and `s2`
(`gfrefn` is the bisection of CSPICE); `tol` is the convergence tolerance in seconds,
an argument here as in CSPICE.

Contract: as `gfoclt`; `gfocce(..., tol, step, gfrefn, cnfine)` is `gfoclt(..., step,
cnfine; tol)`.
"""
function gfocce(ks::KernelSet, occtyp::AbstractString, front::NameOrID,
        fshape::AbstractString, fframe::NameOrID, back::NameOrID, bshape::AbstractString,
        bframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, tol::Real, udstep, udrefn,
        cnfine::Window; bail=nothing, report=nothing, threads::Bool=false)
    return _gf_oclt(ks, _DSK_DEFAULT_TOLS, occtyp, front, fshape, fframe, back, bshape,
        bframe, abcorr, obsrvr, cnfine,
        _gf_run(udstep, tol; refine=udrefn, bail, report, threads))
end
gfocce(ctx::SpiceContext, occtyp::AbstractString, front::NameOrID, fshape::AbstractString,
    fframe::NameOrID, back::NameOrID, bshape::AbstractString, bframe::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, tol::Real, udstep, udrefn, cnfine::Window;
    bail=nothing, report=nothing, threads::Bool=false) =
    _gf_oclt(kernelset(ctx), dsktolerances(ctx), occtyp, front, fshape, fframe, back,
        bshape, bframe, abcorr, obsrvr, cnfine,
        _gf_run(udstep, tol; refine=udrefn, bail, report, threads))
gfocce(occtyp::AbstractString, front::NameOrID, fshape::AbstractString, fframe::NameOrID,
    back::NameOrID, bshape::AbstractString, bframe::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, tol::Real, udstep, udrefn, cnfine::Window; kw...) =
    gfocce(default_context(), occtyp, front, fshape, fframe, back, bshape, bframe, abcorr,
        obsrvr, tol, udstep, udrefn, cnfine; kw...)
