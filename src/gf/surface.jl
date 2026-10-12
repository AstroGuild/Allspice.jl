# Geometry finder, coordinates of a surface point: the vector sources `SUB-OBSERVER
# POINT` and `SURFACE INTERCEPT POINT` of the coordinate quantity, and the calls `gfsubc`
# and `gfsntc` (CSPICE N0067: `gfsubc`, `gfsntc`, `zzgfcoin`, `zzgfcoq`, `zzgfcost`,
# `zzgfssob`, `zzgfssin`, and the existence window of `zzgfcslv`).
#
# The value of a coordinate comes from the position that `subpnt` or `sincpt` gives. The
# test "is decreasing" comes from the state of the surface point, and that state is the
# one of CSPICE (`zzgfssob`, `zzgfssin`), not the exact derivative of the position:
# CSPICE holds the target epoch fixed, scales the velocity of the target and the
# derivative block of the frame transformation by the rate of the light time, and
# repeats this two or three times. The same formula is used here, so that the brackets
# of the extrema are the brackets of CSPICE.
#
# CSPICE accepts one shape model only in these searches, the reference ellipsoid, and
# only the old spellings of the method: `Near point: ellipsoid` and
# `Intercept: ellipsoid` for the sub-observer point, `Ellipsoid` for the intercept. Any
# other method throws `SpiceError(:notsupported)` at the first evaluation of the state,
# as in CSPICE.

"Apply the state transformation with the rotation `R` and its derivative `dR` to the state `(r, v)`."
@inline _gf_xf(R, dR, r, v) = (R * r, dR * r + R * v)

"Apply the inverse of the state transformation with the blocks `R` and `dR` to the state `(r, v)` (`invstm`)."
@inline _gf_xfi(R, dR, r, v) = (transpose(R) * r, transpose(dR) * r + transpose(R) * v)

"The 3-vector zero and the state zero of the surface sources."
const _GF_ZERO3 = zero(Vec3{Float64})
const _GF_ZERO6 = zero(SVector{6,Float64})

@noinline _gf_bad_method(what::String, method) = throw(SpiceError(:notsupported,
    "$what computation method $method is not supported by this routine."))

"""
    _gf_fixed_frame(ks, ref, targ)

The frame check of `zzgfssob` and `zzgfssin` for a geometric state: the frame `ref` is
known (`SpiceError(:noframe)`) and centered on `targ` (`:invalidframe`).
"""
function _gf_fixed_frame(ks::KernelSet, ref::Int32, targ::Int32)
    center, _, found = _geo_frinfo(ks, ref)
    found || throw(SpiceError(:noframe, "Input reference frame $ref was not recognized."))
    center == targ || throw(SpiceError(:invalidframe, "Input reference frame $ref is " *
        "centered on body $center instead of body $targ."))
    return nothing
end

"""
    _gf_ltstate(ks, targ, obs, ref, corr, et, spoint, trgepc, nitr, point) -> (state, found)

The light-time part of `zzgfssob` and `zzgfssin`: the state of a surface point of `targ`
in the body-fixed frame `ref`, for the observer `obs` at `et`, from the point `spoint`
at the target epoch `trgepc` (the outputs of `subpnt` or `sincpt`).

`point(obspos, obsvel, R, dR)` gives the state `(position, velocity)` of the surface
point for the body-fixed state of the observer, or `nothing`; `R` and `dR` are the
blocks of the corrected state transformation from J2000 to `ref`. The iteration of
CSPICE is kept: the target epoch is held at `trgepc`; in each of the `nitr` passes the
velocity of the target and `dR` are multiplied by `1 - dlt` (`1 + dlt` for
transmission), where `dlt` is the rate of the light time of the pass before (zero in
the first pass), and for the `+S` forms the stellar aberration offset of the surface
point and its rate (`zzstelab`) are added to the observer-target state.

Contract: `corr.uselt`. `found` is false if `point` returns `nothing`. No allocation
for a `point` that does not allocate.
"""
function _gf_ltstate(ks::KernelSet, targ::Int32, obs::Int32, ref::Int32, corr::Aberration,
        et::Float64, spoint::Vec3{Float64}, trgepc::Float64, nitr::Int, point::F) where {F}
    xmit = corr.transmit
    robs, vobs = _ssb(ks, obs, et, _GF_J2000)
    rt0, vt0 = _ssb(ks, targ, trgepc, _GF_J2000)
    R, dR = frame_jet(ks, _GF_J2000, ref, trgepc, Val(1))
    acc = zero(vobs)
    if corr.stellar
        # `qderiv` of the velocity of the observer over one second on each side
        acc = 0.5 * _ssb(ks, obs, et + 1.0, _GF_J2000)[2] +
              (-0.5) * _ssb(ks, obs, et - 1.0, _GF_J2000)[2]
    end
    pr = spoint
    pv = _GF_ZERO3
    dlt = 0.0
    # invariant: `(pr, pv)` is the state of the surface point after the passes made, and
    # `dlt` the rate of the light time that this state gives
    for _ in 1:nitr
        scale = xmit ? dlt + 1.0 : 1.0 - dlt
        otr = rt0 - robs
        otv = scale * vt0 - vobs
        dRs = scale * dR
        jr, jv = _gf_xfi(R, dRs, pr, pv)
        opr = otr + jr
        opv = otv + jv
        sr, sv = otr, otv
        if corr.stellar
            sa, savel = _zzstelab(xmit, acc, vobs, opr, opv)
            sr = otr + sa
            sv = otv + savel
        end
        tr, tv = _gf_xf(R, dRs, sr, sv)
        st = point(-tr, -tv, R, dRs)
        st === nothing && return _GF_ZERO6, false
        pr, pv = st
        jr, jv = _gf_xfi(R, dRs, pr, pv)
        opr = otr + jr
        opv = otv + jv
        dlt = vdot(opv, vhat(opr)) / clight()
    end
    return vcat(pr, pv), true
end

# --- sub-observer point ------------------------------------------------------------------

"""
    GFSubObserverVector

The vector `SUB-OBSERVER POINT`: the sub-observer point on the reference ellipsoid of
`targ` in its body-fixed frame `ref` (`zzgfcoq` with `subpnt`, `zzgfssob`). `method` is
the method string in upper case; `near` is 1 for the nearest point, 0 for the intercept
of the line to the center, and -1 for a method that CSPICE does not accept here.
`radii` are the radii of the center of `ref`.
"""
struct GFSubObserverVector <: GFVector
    ks::KernelSet
    targ::Int32
    obs::Int32
    corr::Aberration
    ref::Int32
    method::String
    near::Int8
    radii::NTuple{3,Float64}
end

_gf_ks(v::GFSubObserverVector) = v.ks
_gf_vecframe(v::GFSubObserverVector) = v.ref
_gf_vecexist(::GFSubObserverVector, cnfine::Window, run::_GFRun) = copy(cnfine)

@inline function _gf_vecpos(v::GFSubObserverVector, et::Float64)
    spoint, _, _ = _geo_subpnt(v.ks, _DSK_DEFAULT_TOLS, v.method, v.targ, et, v.ref, v.corr,
        v.obs)
    return spoint, true
end

"""
    _gf_subpoint_state(near, or, ov, radii) -> (position, velocity)

State of the sub-point for the body-fixed state `(or, ov)` of the observer: `dnearp`
for the nearest point, `surfpv` with the ray toward the center otherwise. Throws
`SpiceError(:degeneratecase)` where the velocity is not defined.
"""
@inline function _gf_subpoint_state(near::Bool, or::Vec3{Float64}, ov::Vec3{Float64},
        radii::NTuple{3,Float64})
    a, b, c = radii
    fxo = vcat(or, ov)
    if near
        res = dnearp(fxo, a, b, c)
        res === nothing && throw(SpiceError(:degeneratecase, "The sub-observer state " *
            "could could not be computed because the velocity was not well defined. " *
            "DNEARP returned \"not found.\""))
        return _gf_pv(res[1])
    end
    st = surfpv(fxo, -fxo, a, b, c)
    st === nothing && throw(SpiceError(:degeneratecase, "The sub-observer state could " *
        "could not be computed because the velocity was not well defined. SURFPV " *
        "returned \"not found.\""))
    return _gf_pv(st)
end

"""
State of the sub-observer point, as `zzgfssob`. Without correction it is the state of
the sub-point for the geometric state of the observer in the body-fixed frame; with
light time it is the result of two passes of `_gf_ltstate` from the point of `subpnt`.
Throws `SpiceError(:notsupported)` for a method that is not `Near point: ellipsoid` or
`Intercept: ellipsoid` (case and blanks ignored).
"""
function _gf_vecstate(v::GFSubObserverVector, et::Float64)
    v.near < 0 && _gf_bad_method("Sub-observer point", v.method)
    near = v.near == 1
    radii = v.radii
    if v.corr.geometric
        _gf_fixed_frame(v.ks, v.ref, v.targ)
        tr, tv = _gf_pv(_gf_spkez(v.ks, v.targ, et, v.ref, v.corr, v.obs)[1])
        pr, pv = _gf_subpoint_state(near, -tr, -tv, radii)
        return vcat(pr, pv), true
    end
    spoint, trgepc, _ = _geo_subpnt(v.ks, _DSK_DEFAULT_TOLS, v.method, v.targ, et, v.ref,
        v.corr, v.obs)
    return _gf_ltstate(v.ks, v.targ, v.obs, v.ref, v.corr, et, spoint, trgepc, 2,
        (or, ov, R, dR) -> _gf_subpoint_state(near, or, ov, radii))
end

# --- surface intercept point -------------------------------------------------------------

"""
    GFInterceptVector

The vector `SURFACE INTERCEPT POINT`: the point where the ray from `obs` along `dvec`
(frame `dref`, whose center is `dctr`) meets the reference ellipsoid of `targ`, in the
body-fixed frame `ref` (`zzgfcoq` with `sincpt`, `zzgfssin`). `ok` is false for a
method that CSPICE does not accept here. The vector does not exist where the ray misses
the target.
"""
struct GFInterceptVector <: GFVector
    ks::KernelSet
    targ::Int32
    obs::Int32
    corr::Aberration
    ref::Int32
    method::String
    ok::Bool
    radii::NTuple{3,Float64}
    dref::Int32
    dctr::Int32
    dvec::Vec3{Float64}
end

_gf_ks(v::GFInterceptVector) = v.ks
_gf_vecframe(v::GFInterceptVector) = v.ref

@inline function _gf_vecpos(v::GFInterceptVector, et::Float64)
    res = _geo_sincpt(v.ks, _DSK_DEFAULT_TOLS, v.method, v.targ, et, v.ref, v.corr, v.obs,
        v.dref, v.dvec)
    res === nothing && return _GF_ZERO3, false
    return res[1], true
end

"""
State of the intercept point, as `zzgfssin`; `found` is false where the ray misses the
target. The state of the ray direction in J2000 comes from the state transformation of
`dref` at `et`, or, for a frame that is not centered on the observer and a correction
with light time, at the epoch corrected for the light time to the center of the frame
(`spkacs`), with the derivative block multiplied by the rate of that epoch. With light
time the result is three passes of `_gf_ltstate` from the point of `sincpt`. Throws
`SpiceError(:notsupported)` for a method that is not `Ellipsoid` (case and blanks
ignored).
"""
function _gf_vecstate(v::GFInterceptVector, et::Float64)
    v.ok || _gf_bad_method("Surface intercept point", v.method)
    ks = v.ks
    a, b, c = v.radii
    corr = v.corr
    if corr.geometric
        Rd, dRd = frame_jet(ks, v.dref, _GF_J2000, et, Val(1))
        j2r, j2v = _gf_xf(Rd, dRd, v.dvec, _GF_ZERO3)
        _gf_fixed_frame(ks, v.ref, v.targ)
        tr, tv = _gf_pv(_gf_spkez(ks, v.targ, et, v.ref, corr, v.obs)[1])
        R, dR = frame_jet(ks, _GF_J2000, v.ref, et, Val(1))
        dr, dv = _gf_xf(R, dR, j2r, j2v)
        st = surfpv(vcat(-tr, -tv), vcat(dr, dv), a, b, c)
        st === nothing && return _GF_ZERO6, false
        return SVector{6,Float64}(st), true
    end
    if v.obs == v.dctr
        Rd, dRd = frame_jet(ks, v.dref, _GF_J2000, et, Val(1))
    else
        s = corr.transmit ? 1.0 : -1.0
        _, ltctr, dltctr = spkacs(ks, v.dctr, et, _GF_J2000, corr, v.obs)
        Rd, dRd = frame_jet(ks, v.dref, _GF_J2000, et + s * ltctr, Val(1))
        dRd = (s * dltctr + 1.0) * dRd
    end
    j2r, j2v = _gf_xf(Rd, dRd, v.dvec, _GF_ZERO3)
    res = _geo_sincpt(ks, _DSK_DEFAULT_TOLS, v.method, v.targ, et, v.ref, corr, v.obs,
        v.dref, v.dvec)
    res === nothing && return _GF_ZERO6, false
    spoint, trgepc, _ = res
    return _gf_ltstate(ks, v.targ, v.obs, v.ref, corr, et, spoint, trgepc, 3,
        function (or, ov, R, dR)
            dr, dv = _gf_xf(R, dR, j2r, j2v)
            st = surfpv(vcat(or, ov), vcat(dr, dv), a, b, c)
            st === nothing && return nothing
            return _gf_pv(SVector{6,Float64}(st))
        end)
end

"""
    _gf_wncond(w, amount) -> Window

Copy of `w` with each interval contracted by `amount` at both ends, as `wncond_c`: an
interval that is shorter than `2 amount` is removed.
"""
function _gf_wncond(w::Window, amount::Float64)
    d = Float64[]
    for i in 1:wncard(w)                 # invariant: d holds the contracted intervals before i
        a = w.data[2i-1] + amount
        b = w.data[2i] - amount
        a > b || push!(d, a, b)
    end
    return _window_of(d)
end

"""
The existence window of `zzgfcslv`: the window of `cnfine` on which the state of the
intercept exists, found with the step of the search, with each interval contracted by
`tol + 1` seconds at both ends, so that the coordinate search stays away from the
epochs at which the ray is tangent to the target.
"""
function _gf_vecexist(v::GFInterceptVector, cnfine::Window, run::_GFRun)
    _gf_repi(run, cnfine, "Intercept existence pass 1 of 1", "done.")
    w = _gf_scan!(Window(), et -> _gf_vecstate(v, et)[2], cnfine, run, false)
    _gf_repf(run)
    return _gf_wncond(w, run.tol + 1.0)
end

# --- the coordinate of a surface point ---------------------------------------------------

"""
    _gf_surface_coordinate(ks, intercept, chkctr, target, frame, method, abcorr, obsrvr,
                           dref, dvec, crdsys, coord) -> GFCoordinate

The coordinate `coord` of the system `crdsys` of the sub-observer point (`intercept`
false) or of the surface intercept of the ray `dvec` of frame `dref` (`intercept`
true), with the checks of `zzgfcoin` in their order: bodies (`:idcodenotfound`,
`:bodiesnotdistinct`), correction (`:invalidoption`), system and coordinate
(`:notsupported`), frame (`:noframe`; `:invalidframe` if `chkctr` and the frame is not
centered on the target), radii of the center of the frame (`:badframe`, `:invalidcount`,
`:badaxislength`; `:notsupported` for unequal equatorial radii with the geodetic and
planetographic systems), and for the intercept the ray (`:zerovector`, `:noframe`).

`chkctr` is the comparison that CSPICE makes with the vector definition as the caller
wrote it: it is true for `gfsubc` and `gfsntc`.
"""
function _gf_surface_coordinate(ks::KernelSet, intercept::Bool, chkctr::Bool,
        target::NameOrID, frame::NameOrID, method::AbstractString, abcorr::AbcorrArg,
        obsrvr::NameOrID, dref::NameOrID, dvec::AbstractVector, crdsys::AbstractString,
        coord::AbstractString)
    targ = _gf_body(ks, target, "target object")
    obs = _gf_body(ks, obsrvr, "observer")
    targ == obs && throw(SpiceError(:bodiesnotdistinct,
        "The observer and target must be distinct objects, but are not: " *
        "OBSRVR = $obsrvr; TARGET = $target."))
    corr = _geo_abcorr(abcorr)
    sys, idx = _gf_sysidx(crdsys, coord)
    ref = _geo_namfrm(ks, frame isa AbstractString ? strip(frame, ' ') : frame)
    center, _, found = _geo_frinfo(ks, ref)
    found || throw(SpiceError(:noframe, "Frame system did not recognize frame $frame."))
    (chkctr && center != targ) && throw(SpiceError(:invalidframe, "Vector definition " *
        "method is $(intercept ? "SURFACE INTERCEPT POINT" : "SUB-OBSERVER POINT"), but " *
        "input reference frame $frame has center $center. For this vector definition, " *
        "the frame must be centered on the target body $target."))
    bodfnd(ks, center, "RADII") || throw(SpiceError(:badframe, "No RADII data in kernel " *
        "pool for frame '$frame' center body $center. Confirm the proper input frame. " *
        "Bodies {0,..,9} represent barycenters and so lack physical properties."))
    radii = _gf_radii(ks, center)
    re = 0.0
    f = 0.0
    sense = 0
    if sys == _GF_GEODETIC || sys == _GF_PLANETOGRAPHIC
        radii[1] == radii[2] || throw(SpiceError(:notsupported, "Central body $center " *
            "of reference frame $frame has radii $(radii[1]) $(radii[2]) $(radii[3]). " *
            "Unequal equatorial ellipsoid radii are not supported for " *
            "$(_GF_SYSNMS[sys]) coordinates."))
        re = radii[1]
        f = (radii[1] - radii[3]) / radii[1]
        if sys == _GF_PLANETOGRAPHIC
            lon = recpgr(ks, center, SVector(0.0, 1.0, 0.0), re, f)[1]
            sense = lon > Float64(pi) ? -1 : 1
        end
    end
    meth = uppercase(lstrip(method, ' '))
    rad = (radii[1], radii[2], radii[3])
    if !intercept
        near = _dsk_eqstr(meth, "Near point: ellipsoid") ? Int8(1) :
               _dsk_eqstr(meth, "Intercept: ellipsoid") ? Int8(0) : Int8(-1)
        vec = GFSubObserverVector(ks, targ, obs, corr, ref, meth, near, rad)
        return GFCoordinate(vec, sys, idx, re, f, sense, center)
    end
    dv = _geo_pos(dvec)
    (dv[1] == 0.0 && dv[2] == 0.0 && dv[3] == 0.0) && throw(SpiceError(:zerovector,
        "Ray's direction vector is the zero vector. This variable might be uninitialized."))
    dfr = _geo_namfrm(ks, dref isa AbstractString ? strip(dref, ' ') : dref)
    dctr, _, dfound = _geo_frinfo(ks, dfr)
    dfound || throw(SpiceError(:noframe, "Frame system did not recognize frame $dref."))
    vec = GFInterceptVector(ks, targ, obs, corr, ref, meth, _dsk_eqstr(meth, "Ellipsoid"),
        rad, dfr, dctr, dv)
    return GFCoordinate(vec, sys, idx, re, f, sense, center)
end

# --- public interface --------------------------------------------------------------------

"""
    gfsubc(target, fixref, method, abcorr, obsrvr, crdsys, coord, relate, refval, adjust,
           step, cnfine; keywords) -> Window

Window of the confinement window `cnfine` on which the coordinate `coord` of the
sub-observer point on `target`, as seen from `obsrvr`, satisfies the relation `relate`
with the reference value `refval`, as `gfsubc_c`. The point is given in the body-fixed
frame `fixref`, which must be centered on `target`. `method` is `"Near point:
ellipsoid"` (the point of the reference ellipsoid nearest to the observer) or
`"Intercept: ellipsoid"` (the point where the line from the observer to the center of
the target meets the ellipsoid), in any case. The coordinate systems and coordinates
are those of [`gfposc`](@ref); the other arguments and the keywords are those of
[`gfdist`](@ref).

Contract: as `gfdist`: the search is the search of CSPICE (`zzgfcslv`, `zzgflong`,
`zzgfrelx`), the value is the coordinate of the point of [`subpnt`](@ref), and the test
"is decreasing" uses the state of the point as `zzgfssob` computes it, so an end point
differs from the one of CSPICE by rounding only, and by no more than about `tol`. The
result is a new valid window that is a subset of `cnfine`. Errors: `:invalidstep`,
`:invalidtolerance`, `:notrecognized` (relation), `:idcodenotfound`,
`:bodiesnotdistinct`, `:invalidoption` (correction), `:notsupported` (system,
coordinate, unequal equatorial radii, or a method other than the two above; the method
is read at the first evaluation, so an empty `cnfine` gives no error, as in CSPICE),
`:noframe`, `:invalidframe` (`fixref` is not centered on `target`), `:badframe` (no
radii), `:valueoutofrange`, `:invalidvalue`, `:degeneratecase` (the velocity of the
point is not defined), and the errors of `subpnt`.
"""
function gfsubc(ks::KernelSet, target::NameOrID, fixref::NameOrID, method::AbstractString,
        abcorr::AbcorrArg, obsrvr::NameOrID, crdsys::AbstractString, coord::AbstractString,
        relate::AbstractString, refval::Real, adjust::Real, step::Real, cnfine::Window;
        tol::Real=_GF_CNVTOL, bail=nothing, report=nothing, threads::Bool=false)
    run = _gf_run(step, tol; bail, report, threads)
    _gf_relcode(relate)
    q = _gf_surface_coordinate(ks, false, true, target, fixref, method, abcorr, obsrvr, 0,
        _GF_ZERO3, crdsys, coord)
    return _gf_search(q, relate, refval, adjust, cnfine, run)
end
gfsubc(ctx::SpiceContext, target::NameOrID, fixref::NameOrID, method::AbstractString,
    abcorr::AbcorrArg, obsrvr::NameOrID, crdsys::AbstractString, coord::AbstractString,
    relate::AbstractString, refval::Real, adjust::Real, step::Real, cnfine::Window;
    tol::Real=gfstol(ctx), kw...) =
    gfsubc(kernelset(ctx), target, fixref, method, abcorr, obsrvr, crdsys, coord, relate,
        refval, adjust, step, cnfine; tol, kw...)
gfsubc(target::NameOrID, fixref::NameOrID, method::AbstractString, abcorr::AbcorrArg,
    obsrvr::NameOrID, crdsys::AbstractString, coord::AbstractString, relate::AbstractString,
    refval::Real, adjust::Real, step::Real, cnfine::Window; kw...) =
    gfsubc(default_context(), target, fixref, method, abcorr, obsrvr, crdsys, coord, relate,
        refval, adjust, step, cnfine; kw...)

"""
    gfsntc(target, fixref, method, abcorr, obsrvr, dref, dvec, crdsys, coord, relate,
           refval, adjust, step, cnfine; keywords) -> Window

Window of `cnfine` on which the coordinate `coord` of the surface intercept on `target`
of the ray from `obsrvr` along `dvec` (a vector of 3 elements in the frame `dref`)
satisfies the relation, as `gfsntc_c`. The point is given in the body-fixed frame
`fixref`, which must be centered on `target`. `method` is `"Ellipsoid"`, in any case.
The coordinate systems and coordinates are those of [`gfposc`](@ref); the other
arguments and the keywords are those of [`gfdist`](@ref).

The search is confined to the epochs at which the ray meets the target: a first pass
with the step `step` finds the window on which the intercept exists, and each interval
of that window is made shorter by `tol + 1` seconds at both ends (`zzgfcslv`). `step`
must therefore also be shorter than the shortest interval on which the intercept
exists, and than the shortest gap between two such intervals.

Contract: as [`gfsubc`](@ref), with the point of [`sincpt`](@ref) and the state of
`zzgfssin`. The result is a subset of the existence window. Errors as `gfsubc`, and
`:zerovector` (`dvec`), `:noframe` (`dref`), `:badarg` (`dvec` has not 3 elements).
"""
function gfsntc(ks::KernelSet, target::NameOrID, fixref::NameOrID, method::AbstractString,
        abcorr::AbcorrArg, obsrvr::NameOrID, dref::NameOrID, dvec::AbstractVector,
        crdsys::AbstractString, coord::AbstractString, relate::AbstractString,
        refval::Real, adjust::Real, step::Real, cnfine::Window; tol::Real=_GF_CNVTOL,
        bail=nothing, report=nothing, threads::Bool=false)
    run = _gf_run(step, tol; bail, report, threads)
    _gf_relcode(relate)
    q = _gf_surface_coordinate(ks, true, true, target, fixref, method, abcorr, obsrvr,
        dref, dvec, crdsys, coord)
    return _gf_search(q, relate, refval, adjust, cnfine, run)
end
gfsntc(ctx::SpiceContext, target::NameOrID, fixref::NameOrID, method::AbstractString,
    abcorr::AbcorrArg, obsrvr::NameOrID, dref::NameOrID, dvec::AbstractVector,
    crdsys::AbstractString, coord::AbstractString, relate::AbstractString, refval::Real,
    adjust::Real, step::Real, cnfine::Window; tol::Real=gfstol(ctx), kw...) =
    gfsntc(kernelset(ctx), target, fixref, method, abcorr, obsrvr, dref, dvec, crdsys,
        coord, relate, refval, adjust, step, cnfine; tol, kw...)
gfsntc(target::NameOrID, fixref::NameOrID, method::AbstractString, abcorr::AbcorrArg,
    obsrvr::NameOrID, dref::NameOrID, dvec::AbstractVector, crdsys::AbstractString,
    coord::AbstractString, relate::AbstractString, refval::Real, adjust::Real, step::Real,
    cnfine::Window; kw...) =
    gfsntc(default_context(), target, fixref, method, abcorr, obsrvr, dref, dvec, crdsys,
        coord, relate, refval, adjust, step, cnfine; kw...)
