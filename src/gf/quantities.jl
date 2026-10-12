# Geometry finder, quantities: for each scalar quantity of `gfevnt` the value and the
# test "is the quantity decreasing", as the utility packages of CSPICE N0067 compute
# them (`zzgfdiu`/`zzgfdiq` distance, `zzgfspu`/`zzsepq` angular separation,
# `zzgfrru`/`zzgfrrq` range rate, `zzgfpau`/`zzgfpaq` phase angle, `zzgfcou`/`zzgfcoq`/
# `zzgfcost`/`zzgfcprx` coordinates of a vector).
#
# A quantity is a small immutable object that holds the `KernelSet` snapshot and the
# checked inputs; the initializer of CSPICE (`zzgf**in`) is its constructor. Two methods
# make it usable by the solver:
#
#     _gf_value(q, et::Float64) -> Float64
#     _gf_decreasing(q, et::Float64) -> Bool
#
# and `_gf_label(q)` gives the progress message. A new quantity is a new type with these
# methods; `solver.jl` does not change. The value uses the position of `spkezp` and the
# test uses the state of `spkez`, as in CSPICE, so that the brackets are the same.

"A scalar quantity of the geometry finder; see the head of `gf/quantities.jl`."
abstract type GFQuantity end

"J2000, the frame of every quantity that does not depend on a frame."
const _GF_J2000 = Int32(1)

"""
    _gf_body(ks, name, what) -> Int32

ID of the body `name` (a name, a string that holds an integer, or an integer), as
`bods2c`. Throws `SpiceError(:idcodenotfound)` for an unknown name; `what` names the
role of the body in the message.
"""
function _gf_body(ks::KernelSet, name::AbstractString, what::AbstractString)
    id = bods2c(ks, name)
    id === nothing && throw(SpiceError(:idcodenotfound, "The $what, '$name', is not a " *
        "recognized name for an ephemeris object. The cause of this problem may be that " *
        "you need an updated version of the SPICE Toolkit."))
    return Int32(id)
end
_gf_body(ks::KernelSet, id::Integer, what::AbstractString) = Int32(id)

"State (km, km/s) of `targ` relative to `obs` in frame `ref`, and the light time (`spkez`)."
@inline _gf_spkez(ks::KernelSet, targ::Int32, et::Float64, ref::Int32, corr::Aberration,
    obs::Int32) = spkezr(ks, targ, et, ref, corr, obs)

"Position (km) of `targ` relative to `obs` in frame `ref`, and the light time (`spkezp`)."
@inline _gf_spkezp(ks::KernelSet, targ::Int32, et::Float64, ref::Int32, corr::Aberration,
    obs::Int32) = spkpos(ks, targ, et, ref, corr, obs)

"Position and velocity parts of a state."
@inline _gf_pv(s::SVector{6,Float64}) =
    (SVector(s[1], s[2], s[3]), SVector(s[4], s[5], s[6]))

"""
    _gf_radii(ks, body) -> SVector{3,Float64}

The three radii of `body` from the kernel pool, as `zzgftreb`. Throws
`SpiceError(:invalidcount)` if the variable does not hold three values and
`SpiceError(:badaxislength)` if one is not positive, and the errors of `bodvcd`.
"""
function _gf_radii(ks::KernelSet, body::Integer)
    r = bodvcd(ks, body, "RADII")
    length(r) == 3 || throw(SpiceError(:invalidcount,
        "Only $(length(r)) axes were found  for ID $body. Three axes expected."))
    for i in 1:3
        r[i] > 0.0 || throw(SpiceError(:badaxislength, "Degenerate case. The $i axis of " *
            "body $body is negative or zero.  Please check the text PCK file."))
    end
    return SVector(r[1], r[2], r[3])
end

# --- distance --------------------------------------------------------------------------

"Distance between observer and target (`zzgfdiu`)."
struct GFDistance <: GFQuantity
    ks::KernelSet
    targ::Int32
    obs::Int32
    corr::Aberration
end

"""
    GFDistance(ks, target, abcorr, obsrvr)

Contract (`zzgfdiin`): the bodies are known (`:idcodenotfound`) and distinct
(`:bodiesnotdistinct`), and `abcorr` is a correction of `spkezr`.
"""
function GFDistance(ks::KernelSet, target::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID)
    targ = _gf_body(ks, target, "target object")
    obs = _gf_body(ks, obsrvr, "observer object")
    targ == obs && throw(SpiceError(:bodiesnotdistinct,
        "The observer and target must be distinct objects, but are not: " *
        "OBSRVR = $obsrvr; TARGET = $target."))
    return GFDistance(ks, targ, obs, Aberration(abcorr))
end

_gf_label(::GFDistance) = "Distance pass # of # "

@inline _gf_value(q::GFDistance, et::Float64) =
    vnorm(_gf_spkezp(q.ks, q.targ, et, _GF_J2000, q.corr, q.obs)[1])

"The distance decreases if and only if the position and the velocity have a negative dot product."
@inline function _gf_decreasing(q::GFDistance, et::Float64)
    r, v = _gf_pv(_gf_spkez(q.ks, q.targ, et, _GF_J2000, q.corr, q.obs)[1])
    return vdot(r, v) < 0.0
end

# --- angular separation ----------------------------------------------------------------

"Angular separation of two targets, each a point or a sphere (`zzgfspu`)."
struct GFSeparation <: GFQuantity
    ks::KernelSet
    bod1::Int32
    bod2::Int32
    obs::Int32
    corr::Aberration
    rad1::Float64
    rad2::Float64
end

"Radius of the shape `shape` of `body`: 0 for `POINT`, the largest radius for `SPHERE`."
function _gf_shape_radius(ks::KernelSet, body::Int32, shape::AbstractString)
    s = uppercase(strip(shape, ' '))
    s == "POINT" && return 0.0
    s == "SPHERE" && return maximum(_gf_radii(ks, body))
    throw(SpiceError(:notrecognized,
        "The body shape, $shape is not recognized.  Supported quantities are: POINT, SPHERE."))
end

"""
    GFSeparation(ks, targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr)

Contract (`zzgfspin`): the three bodies are known (`:idcodenotfound`) and distinct
(`:bodiesnotdistinct`), `abcorr` is a correction of `spkezr`, and each shape is `POINT`
or `SPHERE` (`:notrecognized`); a sphere takes the largest of the radii of its body.
The frames are not used by these shapes, as in CSPICE.
"""
function GFSeparation(ks::KernelSet, targ1::NameOrID, shape1::AbstractString, frame1,
        targ2::NameOrID, shape2::AbstractString, frame2, abcorr::AbcorrArg, obsrvr::NameOrID)
    b1 = _gf_body(ks, targ1, "object name for target 1")
    b2 = _gf_body(ks, targ2, "object name for target 2")
    obs = _gf_body(ks, obsrvr, "object name for the observer")
    (obs == b1 || obs == b2 || b1 == b2) && throw(SpiceError(:bodiesnotdistinct,
        "All three objects associated with an ANGULAR SEPARATION search must be " *
        "distinct. The objects whose angular separation is of interest were $b1 and " *
        "$b2. The observer was $obs."))
    corr = Aberration(abcorr)
    r1 = _gf_shape_radius(ks, b1, shape1)
    r2 = _gf_shape_radius(ks, b2, shape2)
    return GFSeparation(ks, b1, b2, obs, corr, r1, r2)
end

_gf_label(::GFSeparation) = "Angular separation pass # of #"

"Angular radius of a sphere of radius `r` at the range `range` (`zzsepq`)."
@inline function _gf_angrad(q::GFSeparation, body::Int32, r::Float64, range::Float64)
    range > r || throw(SpiceError(:insidebody, "Observer object $(q.obs) located within " *
        "surface of target object $body. Range to target $range, radius of target $r"))
    return asin(clamp(r / range, -1.0, 1.0))
end

"Separation of the limbs: the angle between the positions minus the two angular radii (`zzsepq`)."
@inline function _gf_value(q::GFSeparation, et::Float64)
    p1 = _gf_spkezp(q.ks, q.bod1, et, _GF_J2000, q.corr, q.obs)[1]
    p2 = _gf_spkezp(q.ks, q.bod2, et, _GF_J2000, q.corr, q.obs)[1]
    a1 = _gf_angrad(q, q.bod1, q.rad1, vnorm(p1))
    a2 = _gf_angrad(q, q.bod2, q.rad2, vnorm(p2))
    return vsep(p1, p2) - a1 - a2
end

"""
    _gf_dhfa(state, bodyr) -> Float64

Time derivative of the half angle that a sphere of radius `bodyr` subtends at the
observer, for the state `state` of its center, as `dhfa`. Throws
`SpiceError(:degeneratecase)` for a zero position and `SpiceError(:badgeometry)` if the
observer is not outside the sphere.
"""
@inline function _gf_dhfa(state::SVector{6,Float64}, bodyr::Float64)
    r, v = _gf_pv(state)
    rmag = vnorm(r)
    rmag == 0.0 && throw(SpiceError(:degeneratecase,
        "The position component of the input state vector equals the zero vector."))
    p = vhat(r)
    rngrat = vdot(p, v)
    base = rmag * rmag - bodyr * bodyr
    base <= 0.0 && throw(SpiceError(:badgeometry, "Invalid case. The body radius, " *
        "$bodyr, equals or exceeds the range to the target, $rmag."))
    return -(rngrat * bodyr) / (sqrt(base) * rmag)
end

@inline function _gf_decreasing(q::GFSeparation, et::Float64)
    s1 = _gf_spkez(q.ks, q.bod1, et, _GF_J2000, q.corr, q.obs)[1]
    s2 = _gf_spkez(q.ks, q.bod2, et, _GF_J2000, q.corr, q.obs)[1]
    dtheta = dvsep(s1, s2)
    dtheta = dtheta - _gf_dhfa(s1, q.rad1) - _gf_dhfa(s2, q.rad2)
    return dtheta < 0.0
end

# --- range rate ------------------------------------------------------------------------

"Range rate of the target relative to the observer (`zzgfrru`); `dt` is the half interval of the numerical derivative."
struct GFRangeRate <: GFQuantity
    ks::KernelSet
    targ::Int32
    obs::Int32
    corr::Aberration
    dt::Float64
end

"""
    GFRangeRate(ks, target, abcorr, obsrvr, dt=1.0)

Contract (`zzgfrrin`): the bodies are known (`:idcodenotfound`) and distinct
(`:bodiesnotdistinct`), `abcorr` is a correction of `spkezr`, and `dt != 0`
(`:invalidvalue`); `gfevnt` of CSPICE uses one second.
"""
function GFRangeRate(ks::KernelSet, target::NameOrID, abcorr::AbcorrArg, obsrvr::NameOrID,
        dt::Real=1.0)
    targ = _gf_body(ks, target, "target object")
    obs = _gf_body(ks, obsrvr, "observer object")
    targ == obs && throw(SpiceError(:bodiesnotdistinct,
        "The observer and target must be distinct objects, but are not: " *
        "OBSRVR = $obsrvr; TARGET = $target."))
    corr = Aberration(abcorr)
    dt == 0 && throw(SpiceError(:invalidvalue,
        "The timestep value for the numerical derivative is zero."))
    return GFRangeRate(ks, targ, obs, corr, Float64(dt))
end

_gf_label(::GFRangeRate) = "Range Rate pass # of #"

@inline _gf_value(q::GFRangeRate, et::Float64) =
    dvnorm(_gf_spkez(q.ks, q.targ, et, _GF_J2000, q.corr, q.obs)[1])

"""
The range rate decreases if its derivative is negative. CSPICE has no acceleration, so
the derivative of the state is the difference quotient of the states at `et - dt` and
`et + dt` (`qderiv`), and the derivative of the range rate is
`acc . rhat + vel . d(rhat)/dt`.
"""
@inline function _gf_decreasing(q::GFRangeRate, et::Float64)
    s0 = _gf_spkez(q.ks, q.targ, et - q.dt, _GF_J2000, q.corr, q.obs)[1]
    s2 = _gf_spkez(q.ks, q.targ, et + q.dt, _GF_J2000, q.corr, q.obs)[1]
    dfdt = (0.5 / q.dt) * s2 + (-0.5 / q.dt) * s0
    s = _gf_spkez(q.ks, q.targ, et, _GF_J2000, q.corr, q.obs)[1]
    srhat = SVector{6,Float64}(dvhat(s))
    drvel = vdot(_gf_pv(dfdt)[2], _gf_pv(srhat)[1]) + vdot(_gf_pv(s)[2], _gf_pv(srhat)[2])
    return drvel < 0.0
end

# --- phase angle -----------------------------------------------------------------------

"Phase angle at the target between the illuminator and the observer (`zzgfpau`)."
struct GFPhase <: GFQuantity
    ks::KernelSet
    targ::Int32
    illum::Int32
    obs::Int32
    corr::Aberration
end

"""
    GFPhase(ks, target, illmn, abcorr, obsrvr)

Contract (`zzgfpain`): the three bodies are known (`:idcodenotfound`) and distinct
(`:bodiesnotdistinct`), and `abcorr` is a reception correction of `spkezr`; a
transmission correction throws `SpiceError(:invalidoption)`.
"""
function GFPhase(ks::KernelSet, target::NameOrID, illmn::NameOrID, abcorr::AbcorrArg,
        obsrvr::NameOrID)
    targ = _gf_body(ks, target, "target object")
    illum = _gf_body(ks, illmn, "illuminator object")
    obs = _gf_body(ks, obsrvr, "observer object")
    (targ == obs || targ == illum || obs == illum) && throw(SpiceError(:bodiesnotdistinct,
        "The observer, illuminator, and target must be distinct objects, but are not: " *
        "OBSRVR = $obsrvr, TARGET = $target, are not: ILLUM= $illmn."))
    corr = Aberration(abcorr)
    corr.transmit && throw(SpiceError(:invalidoption, "Invalid aberration correction " *
        "'$(string(corr))'. Phase angle geometry calculations currently restricted to " *
        "reception cases."))
    return GFPhase(ks, targ, illum, obs, corr)
end

_gf_label(::GFPhase) = "Phase angle search pass # of #"

"""
The phase angle as `zzgfpaq` computes it: pi minus the angle between the position of the
target seen from the observer and the position of the illuminator seen from the target
at the epoch `et - lt` (at `et` without correction).
"""
@inline function _gf_value(q::GFPhase, et::Float64)
    p1, lt = _gf_spkezp(q.ks, q.targ, et, _GF_J2000, q.corr, q.obs)
    ett = q.corr.geometric ? et : et - lt
    p2 = _gf_spkezp(q.ks, q.illum, ett, _GF_J2000, q.corr, q.targ)[1]
    return Float64(pi) - vsep(p1, p2)
end

"The phase angle decreases if the angle of `_gf_value` increases (`dvsep`); the velocity of the illuminator is scaled by the rate of the epoch `et - lt`."
@inline function _gf_decreasing(q::GFPhase, et::Float64)
    s1, lt = _gf_spkez(q.ks, q.targ, et, _GF_J2000, q.corr, q.obs)
    if q.corr.geometric
        s2 = _gf_spkez(q.ks, q.illum, et, _GF_J2000, q.corr, q.targ)[1]
    else
        s = _gf_spkez(q.ks, q.illum, et - lt, _GF_J2000, q.corr, q.targ)[1]
        r1, v1 = _gf_pv(s1)
        dlt = vdot(vhat(r1), v1) / clight()
        r2, v2 = _gf_pv(s)
        s2 = vcat(r2, (1.0 - dlt) * v2)
    end
    return -dvsep(s1, s2) < 0.0
end

# --- user-defined scalar function ------------------------------------------------------

"A quantity given by two functions of the epoch: the value and the test \"is decreasing\" (`gfuds`)."
struct GFFunction{F,D} <: GFQuantity
    udfuns::F
    udqdec::D
end

_gf_label(::GFFunction) = ""
@inline _gf_value(q::GFFunction, et::Float64) = Float64(q.udfuns(et))::Float64
@inline _gf_decreasing(q::GFFunction, et::Float64) = q.udqdec(et)::Bool

# --- coordinates of a vector -----------------------------------------------------------

# Coordinate systems, in the order of `zzgfcoq`.
const _GF_RECTANGULAR = 1
const _GF_LATITUDINAL = 2
const _GF_RADEC = 3
const _GF_SPHERICAL = 4
const _GF_CYLINDRICAL = 5
const _GF_GEODETIC = 6
const _GF_PLANETOGRAPHIC = 7

"Names of the coordinate systems, blanks removed, in the order of the codes."
const _GF_SYSNMS = ("RECTANGULAR", "LATITUDINAL", "RA/DEC", "SPHERICAL", "CYLINDRICAL",
    "GEODETIC", "PLANETOGRAPHIC")

"Names of the three coordinates of each system, in the order of the conversion call."
const _GF_CRDNMS = (("X", "Y", "Z"), ("RADIUS", "LONGITUDE", "LATITUDE"),
    ("RANGE", "RIGHT ASCENSION", "DECLINATION"), ("RADIUS", "COLATITUDE", "LONGITUDE"),
    ("RADIUS", "LONGITUDE", "Z"), ("LONGITUDE", "LATITUDE", "ALTITUDE"),
    ("LONGITUDE", "LATITUDE", "ALTITUDE"))

"""
    GFVector

Source of the vector whose coordinate is searched (the `VECTOR DEFINITION` of `gfevnt`).
A source has the methods

    _gf_vecpos(v, et) -> (position::SVector{3,Float64}, found::Bool)
    _gf_vecstate(v, et) -> (state::SVector{6,Float64}, found::Bool)
    _gf_vecframe(v) -> Int32          # frame of the vector
    _gf_vecexist(v, cnfine, run) -> Window

`_gf_vecexist` is the window of `cnfine` on which the vector exists (`cnfine` itself
for a vector that always exists). `GFPositionVector` is the source `POSITION`; the
sources `SUB-OBSERVER POINT` and `SURFACE INTERCEPT POINT` are added with `gfsubc` and
`gfsntc` (`gf/stubs.jl`).
"""
abstract type GFVector end

"The vector `POSITION`: the position of the target relative to the observer (`zzgfcost`, `zzgfcoq`)."
struct GFPositionVector <: GFVector
    ks::KernelSet
    targ::Int32
    obs::Int32
    corr::Aberration
    ref::Int32
end

@inline _gf_vecpos(v::GFPositionVector, et::Float64) =
    (_gf_spkezp(v.ks, v.targ, et, v.ref, v.corr, v.obs)[1], true)
@inline _gf_vecstate(v::GFPositionVector, et::Float64) =
    (_gf_spkez(v.ks, v.targ, et, v.ref, v.corr, v.obs)[1], true)
_gf_vecframe(v::GFPositionVector) = v.ref
_gf_vecexist(::GFPositionVector, cnfine::Window, run::_GFRun) = copy(cnfine)

"Left-justified upper-case copy of `s` with every blank removed (`cmprss` with 0)."
_gf_squeeze(s::AbstractString) = uppercase(replace(s, ' ' => ""))
"Left-justified upper-case copy of `s` with blank runs compressed to one (`cmprss` with 1)."
_gf_compress(s::AbstractString) = uppercase(replace(strip(s, ' '), r" +" => " "))

"""
    GFCoordinate

One coordinate of a vector in a coordinate system (`zzgfcou`). `sys` is the system code,
`idx` the place of the coordinate among the three of the system, `re` and `f` the
equatorial radius and the flattening of the center of the frame (geodetic and
planetographic systems; else 0), `sense` +1 or -1 for planetographic longitude positive
east or west (else 0), and `center` the center of the frame of the vector.
"""
struct GFCoordinate{V<:GFVector} <: GFQuantity
    vec::V
    sys::Int
    idx::Int
    re::Float64
    f::Float64
    sense::Int
    center::Int32
end

"""
    _gf_frame(ks, ref) -> Int32

Frame ID of `ref` for a coordinate search. Throws `SpiceError(:unknownframe)` for an
unknown frame, the error of `spkez` in CSPICE.
"""
_gf_frame(ks::KernelSet, ref::NameOrID) = frame_id(ks, ref)

"""
    GFCoordinate(vec, ks, crdsys, crdnam)

Contract (the part of `zzgfcoin` after the checks of the vector): `crdsys` is one of the
seven systems and `crdnam` one of its three coordinates, in any case
(`SpiceError(:notsupported)` otherwise). For the geodetic and planetographic systems the
center of the frame of the vector has three positive radii in the kernel pool
(`:noframe`, `:badframe`, `:invalidcount`, `:badaxislength`) with equal equatorial radii
(`:notsupported`).
"""
function GFCoordinate(vec::GFVector, ks::KernelSet, crdsys::AbstractString,
        crdnam::AbstractString)
    csys = _gf_squeeze(crdsys)
    sys = something(findfirst(==(csys), _GF_SYSNMS), 0)
    sys == 0 && throw(SpiceError(:notsupported,
        "The coordinate system $crdsys is not supported."))
    crd = _gf_compress(crdnam)
    idx = something(findfirst(==(crd), _GF_CRDNMS[sys]), 0)
    idx == 0 && throw(SpiceError(:notsupported, "The coordinate name $crdnam belonging " *
        "to the coordinate system $crdsys is not recognized."))
    re = 0.0
    f = 0.0
    sense = 0
    center = Int32(0)
    if sys == _GF_GEODETIC || sys == _GF_PLANETOGRAPHIC
        ref = _gf_vecframe(vec)
        info = frinfo(ks, ref)
        info === nothing && throw(SpiceError(:noframe,
            "Frame system did not recognize frame $ref."))
        center = Int32(info[1])
        bodfnd(ks, center, "RADII") || throw(SpiceError(:badframe, "No RADII data in " *
            "kernel pool for frame '$ref' center body, ID $center. Geodetic and " *
            "planetographic coordinates require a reference frame centered on a finite " *
            "body. Confirm the proper input frame. Bodies $center may lack radii data in " *
            "the kernel pool."))
        radii = _gf_radii(ks, center)
        radii[1] == radii[2] || throw(SpiceError(:notsupported, "Central body $center " *
            "of reference frame $ref has radii $(radii[1]) $(radii[2]) $(radii[3]). " *
            "Unequal equatorial ellipsoid radii are not supported for $csys coordinates."))
        re = radii[1]
        f = (radii[1] - radii[3]) / radii[1]
        if sys == _GF_PLANETOGRAPHIC
            lon = recpgr(ks, center, SVector(0.0, 1.0, 0.0), re, f)[1]
            sense = lon > Float64(pi) ? -1 : 1
        end
    end
    return GFCoordinate(vec, sys, idx, re, f, sense, center)
end

"Kernel snapshot of a vector source."
_gf_ks(v::GFPositionVector) = v.ks

_gf_label(::GFCoordinate) = "Coordinate pass # of #"

"The three coordinates of the position `pos` in the system of `q`."
@inline function _gf_coords(q::GFCoordinate, pos::SVector{3,Float64})
    sys = q.sys
    sys == _GF_RECTANGULAR && return (pos[1], pos[2], pos[3])
    sys == _GF_LATITUDINAL && return reclat(pos)
    sys == _GF_RADEC && return recrad(pos)
    sys == _GF_SPHERICAL && return recsph(pos)
    sys == _GF_CYLINDRICAL && return reccyl(pos)
    sys == _GF_GEODETIC && return recgeo(pos, q.re, q.f)
    return recpgr(_gf_ks(q.vec), q.center, pos, q.re, q.f)
end

@noinline _gf_notcomputable(et::Float64) = throw(SpiceError(:notcomputable,
    "Coordinate could not be computed at $(etcal(et)) TDB"))

"Value of the coordinate (`zzgfcog`); throws `SpiceError(:notcomputable)` if the vector does not exist."
@inline function _gf_value(q::GFCoordinate, et::Float64)
    pos, found = _gf_vecpos(q.vec, et)
    found || _gf_notcomputable(et)
    return _gf_coords(q, pos)[q.idx]::Float64
end

"`zzrtnmat`: the matrix from the base frame to the frame of the radial, east and north directions of `v`."
@inline function _gf_rtnmat(v::SVector{3,Float64})
    (v[1] == 0.0 && v[2] == 0.0) && throw(SpiceError(:degeneratecase,
        "Input vector ($(v[1]) $(v[2]) $(v[3])) lies on Z-axis; tangential and normal " *
        "directions are undefined."))
    lon = atan(v[2], v[1])
    vlon = SVector(cos(lon), sin(lon), 0.0)
    east = SVector{3,Float64}(ucrss(SVector(0.0, 0.0, 1.0), vlon))
    north = SVector{3,Float64}(ucrss(v, east))
    rad = SVector{3,Float64}(vhat(v))
    return SMatrix{3,3,Float64,9}(rad[1], east[1], north[1], rad[2], east[2], north[2],
        rad[3], east[3], north[3])
end

"Sign of `x` as an integer: -1, 0 or 1."
@inline _gf_isign(x::Float64) = x == 0.0 ? 0 : (x < 0.0 ? -1 : 1)

"""
    _gf_cprx(state, sys, re, f, sense) -> NTuple{3,Int}

Signs (-1, 0, 1) of the time derivatives of the three coordinates of the system `sys`
for the Cartesian state `state`, as `zzgfcprx`: the velocity is expressed in the frame
of the local radial, east and north directions, whose components have the signs of the
derivatives of radius (or altitude), longitude and latitude. On the Z axis the signs of
the longitude and of the latitude are zero.
"""
function _gf_cprx(state::SVector{6,Float64}, sys::Int, re::Float64, f::Float64, sense::Int)
    if sys == _GF_PLANETOGRAPHIC && sense != 1 && sense != -1
        throw(SpiceError(:valueoutofrange, "Longitude sense $sense should be 1 or -1."))
    end
    pos, v = _gf_pv(state)
    (v[1] == 0.0 && v[2] == 0.0 && v[3] == 0.0) && return (0, 0, 0)
    vel = SVector{3,Float64}(vhat(v))
    sys == _GF_RECTANGULAR &&
        return (_gf_isign(vel[1]), _gf_isign(vel[2]), _gf_isign(vel[3]))
    if pos[1] == 0.0 && pos[2] == 0.0
        dpsign = pos[3] != 0.0 ? _gf_isign(vdot(pos, vel)) : 0
        (sys == _GF_LATITUDINAL || sys == _GF_SPHERICAL || sys == _GF_RADEC) &&
            return (dpsign, 0, 0)
        (sys == _GF_GEODETIC || sys == _GF_PLANETOGRAPHIC) && return (0, 0, dpsign)
        return (0, 0, _gf_isign(vel[3]))
    end
    if sys == _GF_GEODETIC || sys == _GF_PLANETOGRAPHIC
        lon, lat, _ = recgeo(pos, re, f)
        normal = SVector{3,Float64}(latrec(1.0, lon, lat))
    elseif sys == _GF_CYLINDRICAL
        normal = SVector{3,Float64}(vhat(SVector(pos[1], pos[2], 0.0)))
    else
        normal = SVector{3,Float64}(vhat(pos))
    end
    rtnvel = _gf_rtnmat(normal) * vel
    r = (_gf_isign(rtnvel[1]), _gf_isign(rtnvel[2]), _gf_isign(rtnvel[3]))
    sys == _GF_SPHERICAL && return (r[1], -r[3], r[2])
    sys == _GF_GEODETIC && return (r[2], r[3], r[1])
    sys == _GF_PLANETOGRAPHIC && return (r[2] * sense, r[3], r[1])
    return r
end

"State of the vector of `q` at `et` (`zzgfcost`); throws `SpiceError(:notcomputable)` if it does not exist."
@inline function _gf_cstate(q::GFCoordinate, et::Float64)
    state, found = _gf_vecstate(q.vec, et)
    found || _gf_notcomputable(et)
    return state
end

"The coordinate decreases if the sign of its derivative is -1 (`zzgfcodc`)."
@inline function _gf_decreasing(q::GFCoordinate, et::Float64)
    state = _gf_cstate(q, et)
    return _gf_cprx(state, q.sys, q.re, q.f, q.sense)[q.idx] == -1
end

# The cosine and the sine of a coordinate, the proxy functions of the longitude searches
# (`zzgfcocg`, `zzgfcosg`, `zzgfcocd`, `zzgfcosd`). The tests take the coordinate from the
# position part of the state, as CSPICE does.
@inline _gf_cosvalue(q::GFCoordinate, et::Float64) = cos(_gf_value(q, et))
@inline _gf_sinvalue(q::GFCoordinate, et::Float64) = sin(_gf_value(q, et))

@inline function _gf_cosdecreasing(q::GFCoordinate, et::Float64)
    state = _gf_cstate(q, et)
    value = _gf_coords(q, _gf_pv(state)[1])[q.idx]::Float64
    return -sin(value) * _gf_cprx(state, q.sys, q.re, q.f, q.sense)[q.idx] < 0.0
end

@inline function _gf_sindecreasing(q::GFCoordinate, et::Float64)
    state = _gf_cstate(q, et)
    value = _gf_coords(q, _gf_pv(state)[1])[q.idx]::Float64
    return cos(value) * _gf_cprx(state, q.sys, q.re, q.f, q.sense)[q.idx] < 0.0
end

# --- the search of a quantity ------------------------------------------------------------

"The two progress messages of a search of `q`: `#` is the pass, then the number of passes."
function _gf_prefix(q::GFQuantity, rel::Int, adjust::Real)
    one = rel == _GF_LOCMIN || rel == _GF_LOCMAX ||
          (adjust == 0 && (rel == _GF_ABSMIN || rel == _GF_ABSMAX))
    n = one ? 1 : 2
    t = _gf_label(q)
    p(i) = replace(replace(t, "#" => string(i); count=1), "#" => string(n); count=1)
    return (p(1), n == 2 ? p(2) : "")
end

"""
    _gf_search(q, relate, refval, adjust, cnfine, run) -> Window

Window of `cnfine` on which the quantity `q` satisfies the relation: the call of
`zzgfrelx` that `gfevnt` makes for `q`. The method for a `GFCoordinate` is in
`gf/longitude.jl`.
"""
function _gf_search(q::GFQuantity, relate::AbstractString, refval::Real, adjust::Real,
        cnfine::Window, run::_GFRun)
    rel = _gf_relcode(relate)
    return _gf_relx(et -> _gf_value(q, et), et -> _gf_decreasing(q, et), relate, refval,
        adjust, cnfine, run, _gf_prefix(q, rel, adjust), "done.")
end
