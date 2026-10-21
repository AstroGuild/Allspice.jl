# Geometry finder, visibility in a field of view: `gftfov`, `gfrfov` and `gffove` of
# CSPICE N0067, with the initializer `zzgffvin` and the state `zzgffvst` of `zzgffvu`.
#
# A visibility search has no scalar quantity: it is the search for the intervals on
# which a binary state is true, as `gfoclt` is. `fovtrg` and `fovray`
# (`geometry/fov.jl`) make the checks and the test in one call for one epoch; a search
# makes the checks once and the test at many epochs, so the two parts of `zzgffvu` are
# separate here: `_gf_fov` reads the field of view and makes every check, and
# `_gf_visible` is the test at an epoch. The test is the one of `_geo_fovvis`, step for
# step.

"""
    _GFFov

The checked inputs of a visibility search, the saved values of `zzgffvin`.

- `obs`, `targ`: observer and target (0 for a ray); `tshape`: `_GEO_FOV_POINT`,
  `_GEO_FOV_ELLIPSOID` or `_GEO_FOV_RAY`; `tfid`: frame of the target or of the ray.
- `corr`: correction of a target; `ustl`, `xmit`: stellar aberration of a ray and its
  direction; `trad`: radii of an ellipsoid; `rdir`: unit direction of a ray.
- `ifid`, `ifrname`: frame of the instrument; `poly`, `circle`: shape of the field of
  view; `faxi`: its unit axis; `arad`: its angular radius; `cang`: the half angle of a
  circular field of view; `bounds`: the boundary vectors; `plane`: the plane normal to
  the axis at distance 1; `fvct`: the center of the field of view in that plane.
- `semi1`, `semi2`, `xmag1`, `xmag2`: semi-axes of an elliptical field of view in the
  plane and their lengths; `edct`, `fsmx`: center and semi-axes of the small ellipsoid
  that stands for an elliptical field of view in the test against an ellipsoid.
- `fx`, `fy`, `fpol`: the two axes of the plane and the vertices of a polygonal field
  of view in it, for the test of a point or a ray.
"""
struct _GFFov
    ks::KernelSet
    obs::Int32
    targ::Int32
    tshape::Int
    tfid::Int32
    corr::Aberration
    ustl::Bool
    xmit::Bool
    trad::NTuple{3,Float64}
    rdir::Vec3{Float64}
    ifid::Int32
    ifrname::String
    poly::Bool
    circle::Bool
    faxi::Vec3{Float64}
    arad::Float64
    cang::Float64
    bounds::Vector{Vec3{Float64}}
    plane::SpicePlane{Float64}
    fvct::Vec3{Float64}
    semi1::Vec3{Float64}
    semi2::Vec3{Float64}
    xmag1::Float64
    xmag2::Float64
    edct::Vec3{Float64}
    fsmx::Mat3{Float64}
    fx::Vec3{Float64}
    fy::Vec3{Float64}
    fpol::Vector{SVector{2,Float64}}
end

"""
    _gf_fov(ks, inst, tshape, raydir, target, tframe, abcorr, obsrvr) -> _GFFov

`zzgffvin`: the checks of a visibility search, in the order of CSPICE, and the values
that do not depend on the epoch.

Contract: the observer is known (`:idcodenotfound`); `tshape` is `POINT`, `ELLIPSOID`
or `RAY` (`:invalidshape`); a target is known and differs from the observer
(`:idcodenotfound`, `:bodiesnotdistinct`); the frame of an ellipsoid or of a ray is
not blank, is known (`:invalidframe`, `:frameinfonotfound`) and, for an ellipsoid, is
centered on the target (`:invalidframe`); the correction is one of `spkezr` for a
target and `NONE`, `S` or `XS` for a ray (`:invalidoption`); an ellipsoid has three
positive radii; a ray is not zero (`:zerovector`); the instrument is known
(`:idcodenotfound`) and its field of view is in the kernel pool (the errors of
`getfov`), with an angular radius below 90 degrees minus 1e-6 radians (`:fovtoowide`)
and the errors of the axis of a polygon (`:degeneratecase`, `:facenotfound`,
`:notsupported`). The frame of the instrument is looked up at the first evaluation.
"""
function _gf_fov(ks::KernelSet, inst::NameOrID, tshape::AbstractString,
        raydir::Vec3{Float64}, target::NameOrID, tframe::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID)
    obs = _geo_body(ks, obsrvr, "observer")
    shp = _geo_fov_tshape(tshape)
    uray = shp == _GEO_FOV_RAY
    etrg = shp == _GEO_FOV_ELLIPSOID
    targ = Int32(0)
    if !uray
        targ = _geo_body(ks, target, "target")
        targ == obs && _geo_same_bodies("field of view test", obsrvr)
    end
    tfid = Int32(0)
    if uray || etrg
        (tframe isa AbstractString && all(==(' '), tframe)) && _geo_fov_error(:invalidframe,
            "the name of the reference frame is blank")
        tfid = _geo_namfrm(ks, tframe)
        tfid == Int32(0) && _geo_fov_error(:invalidframe,
            "the reference frame $tframe is not recognized")
        center, _, found = _geo_frinfo(ks, tfid)
        found || _geo_fov_error(:frameinfonotfound,
            "no data are available for the reference frame $tframe")
        (etrg && center != targ) && _geo_bad_center(tframe, target, center)
    end
    ustl = false
    xmit = false
    corr = _aberration(0x00)
    if uray
        ustl, xmit = _geo_fov_raycorr(abcorr)
    else
        corr = _geo_abcorr(abcorr)
    end
    trad = etrg ? _geo_radii(ks, targ) : (0.0, 0.0, 0.0)
    rdir = _GF_ZERO3
    if uray
        _m_iszero(raydir) && _geo_fov_error(:zerovector,
            "the direction of the ray is the zero vector")
        rdir = _m_hat(raydir)
    end
    instid = _geo_body(ks, inst, "instrument")
    shape, ifrname, bsite, bounds = _geo_getfov(ks, instid, _GEO_MAXVRT)
    ifid = _geo_namfrm(ks, ifrname)
    poly = shape == "POLYGON" || shape == "RECTANGLE"
    circle = shape == "CIRCLE"
    faxi = poly ? _geo_fovaxi(bounds) : _m_hat(bsite)
    arad = 0.0
    for b in bounds                      # invariant: arad is the largest angle so far
        arad = max(arad, vsep(b, faxi))
    end
    arad > _M_HALFPI - 1e-6 && _geo_fov_error(:fovtoowide,
        "the field of view of the instrument $inst has the angular radius " *
        "$(rad2deg(arad)) degrees; the limit is $(rad2deg(_M_HALFPI - 1e-6)) degrees")
    plane = nvc2pl(faxi, 1.0)
    fvct = _m_hat(faxi)
    origin = _GF_ZERO3
    semi1 = _GF_ZERO3
    semi2 = _GF_ZERO3
    xmag1 = 0.0
    xmag2 = 0.0
    if shape == "ELLIPSE" || (circle && etrg)
        b2 = circle ? _m_vrotv(bounds[1], faxi, _M_HALFPI) : bounds[2]
        nx, p = inrypl(origin, bounds[1], plane)
        nx == 1 || _geo_fov_error(:degeneratecase,
            "the first boundary vector does not meet the plane of the field of view")
        semi1 = p - fvct
        xmag1 = _m_norm(semi1)
        xmag1 == 0.0 && _geo_fov_error(:degeneratecase, "the field of view has a zero semi-axis")
        nx, p = inrypl(origin, b2, plane)
        nx == 1 || _geo_fov_error(:degeneratecase,
            "the second boundary vector does not meet the plane of the field of view")
        semi2 = p - fvct
        xmag2 = _m_norm(semi2)
        xmag2 == 0.0 && _geo_fov_error(:degeneratecase, "the field of view has a zero semi-axis")
    end
    edct = _GF_ZERO3
    fsmx = zero(Mat3{Float64})
    if etrg && !poly
        # The field of view as a small ellipsoid, for the occultation test.
        edct = (xmag1 * xmag1 + 1.0) * fvct
        r3 = xmag1 * sqrt(xmag1^2 + 1.0)
        c3 = r3 * fvct
        c1 = r3 * _m_hat(semi1)
        r2 = xmag2 / xmag1 * r3
        c2 = r2 * _m_ucrss(c3, c1)
        escale = 1e-5 / (_m_norm(edct) + max(r3, r2))
        fsmx = hcat(escale * c1, escale * c2, escale * c3)
        edct = escale * edct
    end
    cang = circle ? vsep(faxi, bounds[1]) : 0.0
    fx = _GF_ZERO3
    fy = _GF_ZERO3
    fpol = SVector{2,Float64}[]
    if poly && !etrg
        _, fx, fy = _geo_frame(faxi)
        for i in eachindex(bounds)       # invariant: fpol holds the vertices before i
            nx, p = inrypl(origin, bounds[i], plane)
            nx == 1 || _geo_fov_error(:degeneratecase,
                "the boundary vector $i does not meet the plane of the field of view")
            q = p - fvct
            push!(fpol, SVector(_m_dot(fx, q), _m_dot(fy, q)))
        end
    end
    return _GFFov(ks, obs, targ, shp, tfid, corr, ustl, xmit, trad, rdir, ifid,
        String(ifrname), poly, circle, faxi, arad, Float64(cang), bounds,
        SpicePlane{Float64}(plane.normal, plane.constant), fvct, semi1, semi2, xmag1, xmag2,
        edct, fsmx, fx, fy, fpol)
end

"""
    _gf_visible(v, et) -> Bool

`zzgffvst`: true if at `et` the target of `v` (an ellipsoid, a point or a ray) is in
the field of view. The frame of the instrument is evaluated at `et`, and the frame of
an ellipsoid at the epoch corrected for the light time.

Contract: the frame of the instrument is known (`SpiceError(:unknownframe)`). The test
of an ellipsoid against a polygonal field of view allocates (`_geo_elvupy`); the other
tests do not.
"""
function _gf_visible(v::_GFFov, et::Float64)
    ks = v.ks
    ifid = v.ifid
    ifid == Int32(0) && _geo_unknown_frame(v.ifrname)
    origin = _GF_ZERO3
    if v.tshape == _GEO_FOV_ELLIPSOID
        trad = v.trad
        pos, lt = _geo_spkezp(ks, v.targ, et, v.tfid, v.corr, v.obs)
        ettarg = et + _geo_sign(v.corr) * lt
        insmat = _m_mxm(frame_rotation(ks, _GF_J2000, ifid, et),
            frame_rotation(ks, v.tfid, _GF_J2000, ettarg))
        if v.poly
            limb = edlimb(trad[1], trad[2], trad[3], -pos)
            lc, s1, s2 = el2cgv(limb)
            fvlimb = _m_cgv2el(_m_mxv(insmat, pos + lc), _m_mxv(insmat, s1), _m_mxv(insmat, s2))
            return _geo_elvupy(fvlimb, v.faxi, v.bounds)
        end
        trgsmx = hcat(trad[1] * insmat[:, 1], trad[2] * insmat[:, 2], trad[3] * insmat[:, 3])
        return _geo_occed(origin, v.edct, v.fsmx, _m_mxv(insmat, pos), trgsmx) != 0
    end
    if v.tshape == _GEO_FOV_RAY
        dir = _m_mxv(frame_rotation(ks, v.tfid, ifid, et), v.rdir)
        if v.ustl
            _, vobs = _ssb(ks, v.obs, et, _GF_J2000)
            m1 = frame_rotation(ks, ifid, _GF_J2000, et)
            dir = _m_mtxv(m1, _stelab(_m_mxv(m1, dir), v.xmit ? -vobs : vobs))
        end
    else
        p, _ = _geo_spkezp(ks, v.targ, et, _GF_J2000, v.corr, v.obs)
        dir = _m_mxv(frame_rotation(ks, _GF_J2000, ifid, et), p)
    end
    sep = vsep(dir, v.faxi)
    v.circle && return sep <= v.cang
    sep > v.arad && return false
    nx, xpt = inrypl(origin, dir, v.plane)
    nx == 0 && return false
    nx == 1 || _geo_fov_error(:bug, "the ray is in the plane of the field of view")
    fovpt = xpt - v.fvct
    if !v.poly
        c1 = _m_dot(fovpt, v.semi1) / v.xmag1
        c2 = _m_dot(fovpt, v.semi2) / v.xmag2
        return (c1 / v.xmag1)^2 + (c2 / v.xmag2)^2 <= 1.0
    end
    return _geo_wind2d(v.fpol, SVector(_m_dot(v.fx, fovpt), _m_dot(v.fy, fovpt))) != 0
end

"""
    _gf_fove(ks, inst, tshape, raydir, target, tframe, abcorr, obsrvr, cnfine, run) -> Window

`gffove` after the checks of the step and of the tolerance: the checks of `_gf_fov`,
then the search of each interval of `cnfine` for the intervals on which the state
`_gf_visible` is true. Postcondition: a new window, a subset of `cnfine`.
"""
function _gf_fove(ks::KernelSet, inst::NameOrID, tshape::AbstractString,
        raydir::Vec3{Float64}, target::NameOrID, tframe::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, cnfine::Window, run::_GFRun)
    v = _gf_fov(ks, inst, tshape, raydir, target, tframe, abcorr, obsrvr)
    _gf_repi(run, cnfine, "Target visibility search ", "done.")
    result = _gf_scan!(Window(), et -> _gf_visible(v, et), cnfine, run, false)
    _gf_repf(run)
    return result
end

@noinline _gf_ray_in_gftfov() = throw(SpiceError(:invalidoption, "The target shape RAY " *
    "is not supported by this routine. Use the routine GFRFOV instead."))

"""
    gftfov([ctx,] inst, target, tshape, tframe, abcorr, obsrvr, step, cnfine;
           tol=gfstol(ctx), bail=nothing, report=nothing, threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which the body `target` is in
the field of view of the instrument `inst`, as seen from `obsrvr`, as `gftfov_c`.

The result is a new `Window` that is a subset of the confinement window `cnfine`;
`cnfine` is not modified. The test at each epoch is the one of [`fovtrg`](@ref).

# Arguments
- `inst`: name or ID of an instrument whose field of view is in the kernel pool (see
  [`getfov`](@ref)). The instrument frame is evaluated at the epoch under test, as in
  [`fovtrg`](@ref).
- `tshape`: `"ELLIPSOID"` or `"POINT"`, in any case. For an ellipsoid, any part of the
  reference ellipsoid in the field of view counts.
- `tframe`: body-fixed frame of an ellipsoid, centered on `target`. It is not read for a
  point.
- `abcorr`: any aberration correction of [`spkezr`](@ref).
- `step`: step of the search in seconds, positive. The search samples the visibility
  from the start of each interval of `cnfine` at this step, so `step` must be shorter
  than the shortest interval on which the target is visible, and than the shortest gap
  between two such intervals. Each change is then bisected until the bracket is not
  longer than `tol`.

The keywords are those of [`gfdist`](@ref).

Throws `SpiceError(:invalidstep)` if `step <= 0`, `SpiceError(:invalidtolerance)` if
`tol <= 0`, and
- `SpiceError(:invalidoption)` if `tshape` is `"RAY"` (use [`gfrfov`](@ref)) or `abcorr`
  is not a correction;
- `SpiceError(:invalidshape)` if `tshape` is not known;
- `SpiceError(:idcodenotfound)` if the observer, the target or the instrument is not
  known, and `SpiceError(:bodiesnotdistinct)` if the target is the observer;
- `SpiceError(:invalidframe)` if `tframe` is blank, not known or not centered on the
  target, and `SpiceError(:frameinfonotfound)` if the frame has no data (ellipsoid only);
- `SpiceError(:fovtoowide)` if the angular radius of the field of view about its axis
  is not below 90 degrees minus 1e-6 rad;
- `SpiceError(:unknownframe)` if the frame of the instrument is not known. This check is
  made at the first epoch that is searched.

Also throws the kinds of [`srfrec`](@ref) if the radii of an ellipsoid are missing or
not valid, the errors of [`getfov`](@ref), and the kinds of [`fovtrg`](@ref) for the
shape of the field of view and for an observer inside the target. The errors of
[`spkpos`](@ref) come from the epochs that the search samples.

See also [`gfrfov`](@ref), [`gffove`](@ref), [`fovtrg`](@ref).
"""
function gftfov(ks::KernelSet, inst::NameOrID, target::NameOrID, tshape::AbstractString,
        tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window;
        tol::Real=_GF_CNVTOL, bail=nothing, report=nothing, threads::Bool=false)
    _dsk_eqstr(tshape, "RAY") && _gf_ray_in_gftfov()
    run = _gf_run(step, tol; bail, report, threads)
    return _gf_fove(ks, inst, tshape, _GF_ZERO3, target, tframe, abcorr, obsrvr, cnfine, run)
end
gftfov(ctx::SpiceContext, inst::NameOrID, target::NameOrID, tshape::AbstractString,
    tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window;
    tol::Real=gfstol(ctx), kw...) =
    gftfov(kernelset(ctx), inst, target, tshape, tframe, abcorr, obsrvr, step, cnfine;
        tol, kw...)
gftfov(inst::NameOrID, target::NameOrID, tshape::AbstractString, tframe::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window; kw...) =
    gftfov(default_context(), inst, target, tshape, tframe, abcorr, obsrvr, step, cnfine;
        kw...)

"""
    gfrfov([ctx,] inst, raydir, rframe, abcorr, obsrvr, step, cnfine;
           tol=gfstol(ctx), bail=nothing, report=nothing, threads=false) -> Window

Return the window of epochs, in seconds past J2000 TDB, on which the ray from `obsrvr`
along `raydir` is in the field of view of the instrument `inst`, as `gfrfov_c`.

`raydir` is a vector of 3 elements in the frame `rframe`. `abcorr` is `"NONE"`, `"S"`
(the direction is corrected for the stellar aberration of the observer, as for light
that is received) or `"XS"` (as for light that is sent). The test at each epoch is the
one of [`fovray`](@ref). `inst`, `step`, `cnfine`, the keywords and the result are those
of [`gftfov`](@ref).

Throws the errors of [`gftfov`](@ref) for `step`, `tol`, the observer, the instrument
and its field of view, and `SpiceError` with the kind
- `:badarg` if `raydir` does not have 3 elements, `:zerovector` if it is zero;
- `:invalidframe` if `rframe` is blank or not known;
- `:invalidoption` if `abcorr` is not one of the three, such as a correction with light
  time.

See also [`gftfov`](@ref), [`gffove`](@ref), [`fovray`](@ref).
"""
function gfrfov(ks::KernelSet, inst::NameOrID, raydir::AbstractVector, rframe::NameOrID,
        abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window;
        tol::Real=_GF_CNVTOL, bail=nothing, report=nothing, threads::Bool=false)
    run = _gf_run(step, tol; bail, report, threads)
    return _gf_fove(ks, inst, "RAY", _geo_pos(raydir), 0, rframe, abcorr, obsrvr, cnfine, run)
end
gfrfov(ctx::SpiceContext, inst::NameOrID, raydir::AbstractVector, rframe::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, step::Real, cnfine::Window; tol::Real=gfstol(ctx),
    kw...) =
    gfrfov(kernelset(ctx), inst, raydir, rframe, abcorr, obsrvr, step, cnfine; tol, kw...)
gfrfov(inst::NameOrID, raydir::AbstractVector, rframe::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, step::Real, cnfine::Window; kw...) =
    gfrfov(default_context(), inst, raydir, rframe, abcorr, obsrvr, step, cnfine; kw...)

"""
    gffove([ctx,] inst, tshape, raydir, target, tframe, abcorr, obsrvr, tol, udstep,
           udrefn, cnfine; bail=nothing, report=nothing, threads=false) -> Window

Return the visibility window of [`gftfov`](@ref) or [`gfrfov`](@ref), found with a
step and a refinement function that the caller supplies, as `gffove_c`.

`tshape` is `"ELLIPSOID"`, `"POINT"` or `"RAY"`. For a ray, `raydir` is the direction,
`tframe` its frame, and `target` is not read. For a body, `raydir` is not read, but it
must have 3 elements. `tol`, `udstep` and `udrefn` are those of [`gfevnt`](@ref).

The other arguments, the keywords `bail`, `report` and `threads`, the result and the
errors are those of [`gftfov`](@ref) and [`gfrfov`](@ref).
`gffove(inst, tshape, raydir, target, tframe, abcorr, obsrvr, tol, step, gfrefn, cnfine)`
gives the window of `gftfov` or `gfrfov` with the step `step` and the tolerance `tol`.

See also [`gftfov`](@ref), [`gfrfov`](@ref), [`gfevnt`](@ref).
"""
function gffove(ks::KernelSet, inst::NameOrID, tshape::AbstractString,
        raydir::AbstractVector, target::NameOrID, tframe::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, tol::Real, udstep, udrefn, cnfine::Window; bail=nothing,
        report=nothing, threads::Bool=false)
    run = _gf_run(udstep, tol; refine=udrefn, bail, report, threads)
    return _gf_fove(ks, inst, tshape, _geo_pos(raydir), target, tframe, abcorr, obsrvr,
        cnfine, run)
end
gffove(ctx::SpiceContext, inst::NameOrID, tshape::AbstractString, raydir::AbstractVector,
    target::NameOrID, tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, tol::Real,
    udstep, udrefn, cnfine::Window; kw...) =
    gffove(kernelset(ctx), inst, tshape, raydir, target, tframe, abcorr, obsrvr, tol,
        udstep, udrefn, cnfine; kw...)
gffove(inst::NameOrID, tshape::AbstractString, raydir::AbstractVector, target::NameOrID,
    tframe::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, tol::Real, udstep, udrefn,
    cnfine::Window; kw...) =
    gffove(default_context(), inst, tshape, raydir, target, tframe, abcorr, obsrvr, tol,
        udstep, udrefn, cnfine; kw...)
