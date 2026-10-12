# Illumination angles at a surface point and the phase angle between bodies: `ilumin`,
# `illumg`, `illumf` and `phaseq` of CSPICE N0067 (with `zzgfpaq`).
#
# The functions are generic over the epoch number type (a dual number passes through;
# the angles need the arcsine, which the `Taylor` type does not have).

"""
    _geo_illum_setup(ks, method, target, fixref, abcorr, obsrvr)
        -> (corr, trgcde, obscde, fixfid, m)

The checks that `illumg` and `illumf` share, in the order of CSPICE.
"""
@inline function _geo_illum_setup(ks::KernelSet, method::AbstractString, target::NameOrID,
        fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID)
    corr = _geo_abcorr(abcorr)
    trgcde = _geo_body(ks, target, "target")
    m = _geo_method(ks, trgcde, method)
    m.subtyp == 0x00 || _geo_sub_in_method(method)
    obscde = _geo_body(ks, obsrvr, "observer")
    obscde == trgcde && _geo_same_bodies("illumination angles", obsrvr)
    fixfid = _geo_fixfid(ks, fixref, trgcde, target)
    return corr, trgcde, obscde, fixfid, m
end

"""
    _geo_illum_states(ks, corr, trgcde, obscde, srccde, fixfid, spoint, et)
        -> (trgepc, srfvec, tipos)

The two states of `illumg`: the surface point as seen from the observer (`spkcpt` with
the locus `TARGET`), which gives `srfvec` and the epoch `trgepc` at the point, and the
illumination source as seen from the surface point at `trgepc` (`spkcpo` with the
locus `OBSERVER`), which gives `tipos`.
"""
@inline function _geo_illum_states(ks::KernelSet, corr::Aberration, trgcde::Int32,
        obscde::Int32, srccde::Int32, fixfid::Int32, spoint::Vec3{Float64}, et)
    f = _GeoFixed(spoint, zero(Vec3{Float64}), 0.0, trgcde, fixfid)
    opstat, lt = _geo_cvt(ks, f, et, fixfid, _GEO_AT_TARGET, corr, obscde)
    trgepc = et + _geo_sign(corr) * lt
    tistat, _ = _geo_cvo(ks, srccde, trgepc, fixfid, _GEO_AT_OBSERVER, corr, f)
    return trgepc, SVector(opstat[1], opstat[2], opstat[3]), SVector(tistat[1], tistat[2], tistat[3])
end

function _geo_illumg(ks::KernelSet, tols::DSKTolerances, method::AbstractString,
        target::NameOrID, ilusrc::NameOrID, et::Real, fixref::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, spoint::AbstractVector)
    corr, trgcde, obscde, fixfid, m = _geo_illum_setup(ks, method, target, fixref, abcorr, obsrvr)
    sp = _geo_pos(spoint)
    e = check_epoch(et)
    srccde = _geo_body(ks, ilusrc, "illumination source")
    trgepc, srfvec, tipos = _geo_illum_states(ks, corr, trgcde, obscde, srccde, fixfid, sp, e)
    normal = if m.dsk
        _dsk_bun(ks, tols, trgcde, m.srflst, _geo_plain(e), fixfid, sp)
    else
        a, b, c = _geo_radii(ks, trgcde)
        _m_surfnm(a, b, c, sp)
    end
    obspos = -srfvec
    return trgepc, srfvec, vsep(obspos, tipos), vsep(normal, tipos), vsep(normal, obspos)
end

function _geo_illumf(ks::KernelSet, tols::DSKTolerances, method::AbstractString,
        target::NameOrID, ilusrc::NameOrID, et::Real, fixref::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID, spoint::AbstractVector)
    corr, trgcde, obscde, fixfid, m = _geo_illum_setup(ks, method, target, fixref, abcorr, obsrvr)
    sp = _geo_pos(spoint)
    e = check_epoch(et)
    srccde = _geo_body(ks, ilusrc, "illumination source")
    trgepc, srfvec, tipos = _geo_illum_states(ks, corr, trgcde, obscde, srccde, fixfid, sp, e)
    obspos = -srfvec
    if !m.dsk
        a, b, c = _geo_radii(ks, trgcde)
        normal = _m_surfnm(a, b, c, sp)
        phase = vsep(obspos, tipos)
        incdnc = vsep(normal, tipos)
        emissn = vsep(normal, obspos)
        return trgepc, srfvec, phase, incdnc, emissn, primal(emissn) <= _M_HALFPI,
            primal(incdnc) <= _M_HALFPI
    end
    # As in CSPICE, the DSK branch of `illumf` uses every surface of the target (the
    # surface list of the method string is not used) and the epoch at the point.
    sh = _geo_dskshape(ks, tols, trgcde, fixfid, _GEO_NOSURF)
    tp = _geo_plain(trgepc)
    normal = _dsk_vhat(_dsk_bun(ks, tols, trgcde, _GEO_NOSURF, tp, fixfid, sp))
    phase = vsep(obspos, tipos)
    incdnc = vsep(normal, tipos)
    emissn = vsep(normal, obspos)
    visibl = emissn <= _M_HALFPI
    lit = incdnc <= _M_HALFPI
    if lit || visibl
        # A point is hidden or in shadow if the ray from just above it meets the surface.
        vertex = sp + (sh.maxrad * 1e-10) * normal
        visibl && (visibl = !_geo_rayx(ks, sh, vertex, obspos, tp)[1])
        lit && (lit = !_geo_rayx(ks, sh, vertex, tipos, tp)[1])
    end
    return trgepc, srfvec, phase, incdnc, emissn, visibl, lit
end

"""
    illumg(method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
        -> (trgepc, srfvec, phase, incdnc, emissn)

Return the illumination angles at the surface point `spoint` of `target` (km, body-fixed
frame `fixref`) for the illumination source `ilusrc` and the observer `obsrvr`, as
`illumg_c`: the epoch `trgepc` at the point, the vector `srfvec` from the observer to
the point, the phase angle, the incidence angle of the source and the emission angle
(radians). Forms: on a `KernelSet`, on a `SpiceContext` (which gives the DSK
tolerances), and on the default context.

`method` is `"ELLIPSOID"` or `"DSK/UNPRIORITIZED[/SURFACES = <list>]"`; it selects the
surface normal. The arguments are in the order of CSPICE.

Contract: `spoint` has 3 elements and is taken as a surface point without a check;
`fixref` is centered on `target`. For an ellipsoid the epoch can be a dual number, and
the call does not allocate for integer IDs. Throws `SpiceError`: `:invalidoption`,
`:idcodenotfound`, the kinds of the method string, `:invalidmethod`,
`:bodiesnotdistinct`, `:noframe`, `:invalidframe`, and the errors of `spkcpt`, of the
radii lookup and of [`srfnrm`](@ref) for DSK data.
"""
function illumg end

"""
    illumf(method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
        -> (trgepc, srfvec, phase, incdnc, emissn, visibl, lit)

As [`illumg`](@ref), with the two flags of `illumf_c`: `visibl` tells whether the point
is seen from the observer and `lit` whether the source lights it. For an ellipsoid the
flags compare the emission and the incidence angle with 90 degrees; for DSK data a
point also fails if the surface blocks the ray toward the observer or the source.

Contract: as [`illumg`](@ref). As in CSPICE, the DSK branch uses all loaded surfaces of
the target and the epoch `trgepc`, whatever the surface list of `method` says.
"""
function illumf end

for (f, g) in ((:illumg, :_geo_illumg), (:illumf, :_geo_illumf))
    @eval begin
        $f(ks::KernelSet, method::AbstractString, target::NameOrID, ilusrc::NameOrID, et::Real,
            fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, spoint::AbstractVector) =
            $g(ks, _DSK_DEFAULT_TOLS, method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
        $f(ctx::SpiceContext, method::AbstractString, target::NameOrID, ilusrc::NameOrID,
            et::Real, fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
            spoint::AbstractVector) = $g(kernelset(ctx), dsktolerances(ctx), method, target,
            ilusrc, et, fixref, abcorr, obsrvr, spoint)
        $f(method::AbstractString, target::NameOrID, ilusrc::NameOrID, et::Real,
            fixref::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID, spoint::AbstractVector) =
            $f(default_context(), method, target, ilusrc, et, fixref, abcorr, obsrvr, spoint)
    end
end

"""
    ilumin(method, target, et, fixref, abcorr, obsrvr, spoint)
        -> (trgepc, srfvec, phase, incdnc, emissn)

Return the illumination angles at the surface point `spoint` of `target` with the Sun
as the source, as `ilumin_c`: [`illumg`](@ref) with `ilusrc = 10`. The arguments are in
the order of CSPICE.

Contract: as [`illumg`](@ref). For an ellipsoid the epoch can be a dual number, and the
call does not allocate for integer IDs.
"""
ilumin(ks::KernelSet, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, spoint::AbstractVector) =
    _geo_illumg(ks, _DSK_DEFAULT_TOLS, method, target, 10, et, fixref, abcorr, obsrvr, spoint)
ilumin(ctx::SpiceContext, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    abcorr::AbcorrArg, obsrvr::NameOrID, spoint::AbstractVector) =
    _geo_illumg(kernelset(ctx), dsktolerances(ctx), method, target, 10, et, fixref, abcorr,
        obsrvr, spoint)
ilumin(method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID, abcorr::AbcorrArg,
    obsrvr::NameOrID, spoint::AbstractVector) =
    ilumin(default_context(), method, target, et, fixref, abcorr, obsrvr, spoint)

"""
    phaseq(et, target, illmn, obsrvr, abcorr) -> phase

Return the phase angle (radians) at `target` between the directions to the
illumination source `illmn` and to the observer `obsrvr`, at ET `et`, as `phaseq_c`.
With a correction the source is taken at the epoch `et - lt`, with `lt` the light time
between the observer and the target. Forms: on a `KernelSet`, on a `SpiceContext`, and
on the default context.

Contract: the three bodies differ (`SpiceError(:bodiesnotdistinct)`); `abcorr` is a
reception correction (`SpiceError(:invalidoption)` for a transmission correction or an
unknown value). Postcondition: the result is in [0, pi]. The epoch can be a dual
number. No allocation for integer IDs. Throws also `SpiceError(:idcodenotfound)` and
the errors of `spkpos`.
"""
function phaseq(ks::KernelSet, et::Real, target::NameOrID, illmn::NameOrID, obsrvr::NameOrID,
        abcorr::AbcorrArg)
    targ = _geo_body(ks, target, "target")
    illum = _geo_body(ks, illmn, "illumination source")
    obs = _geo_body(ks, obsrvr, "observer")
    corr = _geo_abcorr(abcorr)
    corr.transmit && throw(SpiceError(:invalidoption,
        "the aberration correction $(string(corr)) calls for transmission; phaseq " *
        "accepts reception corrections only"))
    (targ == obs || targ == illum || obs == illum) && throw(SpiceError(:bodiesnotdistinct,
        "the target $target, the illumination source $illmn and the observer $obsrvr " *
        "must be three different bodies"))
    e = check_epoch(et)
    pv1, lt = _geo_spkezp(ks, targ, e, Int32(1), corr, obs)
    pv2, _ = _geo_spkezp(ks, illum, corr.geometric ? e : e - lt, Int32(1), corr, targ)
    return pi - vsep(pv1, pv2)
end
phaseq(ctx::SpiceContext, et::Real, target::NameOrID, illmn::NameOrID, obsrvr::NameOrID,
    abcorr::AbcorrArg) = phaseq(kernelset(ctx), et, target, illmn, obsrvr, abcorr)
phaseq(et::Real, target::NameOrID, illmn::NameOrID, obsrvr::NameOrID, abcorr::AbcorrArg) =
    phaseq(default_context(), et, target, illmn, obsrvr, abcorr)
