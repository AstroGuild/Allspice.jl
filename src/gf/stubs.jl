# Geometry finder, stage 2: the searches that need the high-level geometry calls. Those
# calls (`subpnt`, `sincpt`, `ilumin`, the occultation classifier and the field-of-view
# test) are written in `src/geometry/`; until they are in the package each search here
# throws `SpiceError(:unsupported)`. The solver core, the coordinate solver and the
# longitude solver are complete, so each stub is replaced by a quantity (or a binary
# state) and one call of the core; the comment at each stub says how.

"Throw the error of a search that waits for the geometry call `needs`."
@noinline _gf_stub(what::AbstractString, needs::AbstractString) = throw(SpiceError(:unsupported,
    "$what is not available yet: it needs the geometry call `$needs`"))

"""
    gfsubc(target, fixref, method, abcorr, obsrvr, crdsys, coord, relate, refval, adjust,
           step, cnfine; keywords) -> Window

Window on which a coordinate of the sub-observer point satisfies a relation, as
`gfsubc_c`. Contract: not available yet; throws `SpiceError(:unsupported)`.
"""
# Needs `subpnt` (value, `zzgfcoq`) and the state of the sub-observer point (`zzgfssob`,
# for the test "is decreasing"). Wiring: a `GFVector` subtype `GFSubObserverVector` with
# `_gf_vecpos` = `subpnt(...)[1]`, `_gf_vecstate` = `zzgfssob`, `_gf_vecframe` = `fixref`
# and `_gf_vecexist` = `copy(cnfine)`; then `GFCoordinate(vec, ks, crdsys, coord)` and
# `_gf_search`, as `gfposc` does. `zzgfcoin` also requires that `fixref` is centered on
# the target (`:invalidframe`).
gfsubc(args...; kwargs...) = _gf_stub("gfsubc", "subpnt")

"""
    gfsntc(target, fixref, method, abcorr, obsrvr, dref, dvec, crdsys, coord, relate,
           refval, adjust, step, cnfine; keywords) -> Window

Window on which a coordinate of the surface intercept of a ray satisfies a relation, as
`gfsntc_c`. Contract: not available yet; throws `SpiceError(:unsupported)`.
"""
# Needs `sincpt` (value) and the state of the intercept point (`zzgfssin`). Wiring: a
# `GFVector` subtype `GFInterceptVector` whose `_gf_vecpos` and `_gf_vecstate` return
# `found = false` where the ray misses the target, and whose `_gf_vecexist` is the
# existence window of `zzgfcslv`: `_gf_scan!` on the binary state "the state exists"
# (`zzgfcoex`) over `cnfine` with the step of the search, then each interval contracted
# by `tol + 1` seconds at both ends (`wncond`). `zzgfcoin` also requires a non-zero
# `dvec` (`:zerovector`) and that `fixref` is centered on the target (`:invalidframe`).
gfsntc(args...; kwargs...) = _gf_stub("gfsntc", "sincpt")

"""
    gfilum(method, angtyp, target, illmn, fixref, abcorr, obsrvr, spoint, relate, refval,
           adjust, step, cnfine; keywords) -> Window

Window on which an illumination angle (phase, incidence or emission) at a surface point
satisfies a relation, as `gfilum_c`. Contract: not available yet; throws
`SpiceError(:unsupported)`.
"""
# Needs the illumination angles and their rates at a surface point (`zzilusta`, which is
# `ilumin`/`illumg` with derivatives). Wiring: a `GFQuantity` subtype `GFIllumination`
# with `_gf_value` = the angle `angtyp` and `_gf_decreasing` = "its rate is negative"
# (`zzgfilu`: `zzgfilgq`, `zzgfildc`), then `_gf_search`; the branch `ILLUMINATION ANGLE`
# of `gfevnt` builds the same quantity.
gfilum(args...; kwargs...) = _gf_stub("gfilum", "ilumin")

"""
    gfoclt(occtyp, front, fshape, fframe, back, bshape, bframe, abcorr, obsrvr, step,
           cnfine; keywords) -> Window

Window on which one body occults another, or is in transit across it, as `gfoclt_c`.
Contract: not available yet; throws `SpiceError(:unsupported)`.
"""
# Needs the occultation classifier (`zzocced`, the core of `occult`). Wiring: a binary
# state `et -> the occultation of type occtyp is in progress` (`zzgfocst` of `zzgfocu`)
# and `_gf_scan!(Window(), state, cnfine, run, false)`, as `gfudb` does (`gfocce` calls
# `zzgfsolv` on each interval of `cnfine`).
gfoclt(args...; kwargs...) = _gf_stub("gfoclt", "occult")

"""
    gftfov(inst, target, tshape, tframe, abcorr, obsrvr, step, cnfine; keywords) -> Window

Window on which an ephemeris object is in the field of view of an instrument, as
`gftfov_c`. Contract: not available yet; throws `SpiceError(:unsupported)`.
"""
# Needs the field-of-view test for a target body (`fovtrg`, `zzgffvu`). Wiring: a binary
# state `et -> the target is visible` (`zzgffvst`) and `_gf_scan!`, as `gfudb` does
# (`gffove` calls `zzgfsolv` on each interval of `cnfine`).
gftfov(args...; kwargs...) = _gf_stub("gftfov", "fovtrg")

"""
    gfrfov(inst, raydir, rframe, abcorr, obsrvr, step, cnfine; keywords) -> Window

Window on which a ray is in the field of view of an instrument, as `gfrfov_c`.
Contract: not available yet; throws `SpiceError(:unsupported)`.
"""
# Needs the field-of-view test for a ray (`fovray`, `zzgffvu`). Wiring: as `gftfov`, with
# the target shape `RAY`.
gfrfov(args...; kwargs...) = _gf_stub("gfrfov", "fovray")
