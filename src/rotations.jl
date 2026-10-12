# Frame orientation: text-PCK and binary-PCK body frames (class 2), fixed-offset
# frames (class 4), built-in inertial frames (class 1), and the chain logic that
# joins them. Everything that depends on the epoch is generic over the epoch type, so
# an epoch of type `Taylor{N,Float64}` gives the time derivatives to order N.

const _SEC_PER_DAY = 86400.0                 # seconds per day, `spd_`
const _SEC_PER_CENTURY = _SEC_PER_DAY * 36525.0          # seconds per Julian century
const _J2000_JED = 2451545.0         # `j2000_`
const _RPD = pi / 180                # radians per degree, `rpd_`
const _TWOPI = 2 * Float64(pi)       # `twopi_`
const _HALFPI = Float64(pi) / 2      # `halfpi_`
const _J2000_FRAME = Int32(1)

"Maximum number of frames between a frame and J2000 that `frame_rotation` follows."
const MAX_FRAME_CHAIN = 10

# --- small rotation helpers -----------------------------------------------------------

"`_axis_indices(axis)`: the cyclic index triple `(i1, i2, i3)` of `rotate_` for `axis` in 1:3."
@inline _axis_indices(axis::Integer) =
    axis == 1 ? (1, 2, 3) : axis == 2 ? (2, 3, 1) : (3, 1, 2)

"""
    _rotate(angle, axis) -> Mat3

Coordinate-system rotation by `angle` radians about axis `axis` (1, 2 or 3), as
`rotate_` of CSPICE: the matrix that takes components in the original frame to
components in the rotated frame.
"""
function _rotate(angle::Float64, axis::Integer)
    s, c = sincos(angle)
    i1, i2, i3 = _axis_indices(axis)
    m = zeros(MMatrix{3,3,Float64,9})
    m[i1, i1] = 1.0
    m[i2, i2] = c
    m[i3, i2] = -s
    m[i2, i3] = s
    m[i3, i3] = c
    return Mat3{Float64}(m)
end

"""
    _eul2m(a3, a2, a1, axis3, axis2, axis1) -> Mat3{Float64}

`[a3]_axis3 * [a2]_axis2 * [a1]_axis1`, as `eul2m_` of CSPICE. Requires every axis in 1:3.
"""
function _eul2m(a3::Float64, a2::Float64, a1::Float64, axis3::Integer, axis2::Integer,
    axis1::Integer)
    all(a -> 1 <= a <= 3, (axis3, axis2, axis1)) ||
        throw(SpiceError(:badaxis, "Euler axis numbers are $axis3, $axis2, $axis1; each must be 1, 2 or 3"))
    return _rotate(a3, axis3) * _rotate(a2, axis2) * _rotate(a1, axis1)
end

"""
    _eul313(w, delta, phi) -> Mat3

`[w]_3 * [delta]_1 * [phi]_3`, the same matrix as `eul2m_(w, delta, phi, 3, 1, 3)` with
the arithmetic of `rotate_` and `rotmat_` (the products with the constant 0 and 1
elements are left out, which does not change a result). Generic over the angle type.
"""
@inline function _eul313(w::T, delta::T, phi::T) where {T<:Real}
    sp, cp = sincos(phi)
    sd, cd = sincos(delta)
    sw, cw = sincos(w)
    # rows 2 and 3 after the rotation about axis 1
    m21 = -(cd * sp)
    m22 = cd * cp
    m31 = sd * sp
    m32 = -(sd * cp)
    return Mat3{T}(
        cw * cp + sw * m21, cw * m21 - sw * cp, m31,
        cw * sp + sw * m22, cw * m22 - sw * sp, m32,
        sw * sd, cw * sd, cd)
end

"""
    _reduce_twopi(w)

`w` reduced by a whole number of turns with the arithmetic of `d_mod(w, twopi_())` in
`tisbod_` (the f2c library computes `w - twopi * trunc(w / twopi)` in floating point,
which is not the exact remainder; the same operations are used here so that the
order-0 value equals the one of CSPICE). The subtracted quantity is a constant, so
derivatives are not changed. Generic over the type of `w`.
"""
@inline _reduce_twopi(w::Real) = w - _TWOPI * trunc(primal(w) / _TWOPI)

# --- text PCK constants ---------------------------------------------------------------

"""
    _bodbry(body) -> Int

Barycenter code of the planetary system of `body` (planets and satellites `PNN` and
`PXNNN` give `P`), or `body` itself for any other body, as `zzbodbry_` of CSPICE.
"""
_bodbry(body::Integer) =
    100 <= body <= 999 ? Int(body) ÷ 100 : 10000 <= body <= 99999 ? Int(body) ÷ 10000 : Int(body)

"Numeric pool variable `name`, or `nothing` if it is absent or not numeric."
function _pool_numbers(pool::AbstractDict, name::AbstractString)
    v = get(pool, name, nothing)
    return v isa Vector{Float64} ? v : nothing
end

"Numeric pool variable `name`, or an empty vector."
_pool_numbers_or_empty(pool, name) = something(_pool_numbers(pool, name), Float64[])

"First value of the first of the numeric pool variables `names` that exists, else `default`."
function _pool_first(pool, names, default::Float64)
    for name in names
        v = _pool_numbers(pool, name)
        v !== nothing && !isempty(v) && return v[1]
    end
    return default
end

"""
    build_bodyrot(pool) -> Dict{Int32,BodyRot}

Orientation constants of every body that has numeric `BODY<id>_PM`, `_POLE_RA` and
`_POLE_DEC` in `pool`, read as `tisbod_` reads them: the nutation-precession angles,
the reference frame and the reference epoch come from the barycenter of the body
system (`_bodbry`), the amplitudes from the body itself. `nut_prec_angles` is the
stored vector; its layout is `phase_degree + 1` coefficients per angle, with
`phase_degree` from `BODY<sys>_MAX_PHASE_DEGREE` (default 1).
"""
function build_bodyrot(pool::AbstractDict)
    out = Dict{Int32,BodyRot}()
    for (name, value) in pool
        value isa Vector{Float64} || continue
        m = match(r"^BODY(-?\d+)_PM$", name)
        m === nothing && continue
        id = tryparse(Int32, m.captures[1])
        id === nothing && continue
        ra = _pool_numbers(pool, "BODY$(id)_POLE_RA")
        dec = _pool_numbers(pool, "BODY$(id)_POLE_DEC")
        (ra === nothing || dec === nothing) && continue
        sys = _bodbry(id)
        epoch = _pool_first(pool,
            ("BODY$(sys)_CONSTANTS_JED_EPOCH", "BODY$(sys)_CONSTS_JED_EPOCH"), _J2000_JED)
        ref = _pool_first(pool,
            ("BODY$(sys)_CONSTANTS_REF_FRAME", "BODY$(sys)_CONSTS_REF_FRAME"), 1.0)
        out[id] = BodyRot(ra, dec, value,
            _pool_numbers_or_empty(pool, "BODY$(id)_NUT_PREC_RA"),
            _pool_numbers_or_empty(pool, "BODY$(id)_NUT_PREC_DEC"),
            _pool_numbers_or_empty(pool, "BODY$(id)_NUT_PREC_PM"),
            _pool_numbers_or_empty(pool, "BODY$(sys)_NUT_PREC_ANGLES"),
            round(Int32, ref), epoch, _phase_degree(pool, sys))
    end
    return out
end

"""
    _phase_degree(pool, sys) -> Int32

Value of `BODY<sys>_MAX_PHASE_DEGREE` rounded as CSPICE rounds it, or 1 if the variable
is absent. The value is not checked here; `_phase_ncoef` checks it at evaluation.
"""
function _phase_degree(pool::AbstractDict, sys::Integer)
    v = _pool_numbers(pool, "BODY$(sys)_MAX_PHASE_DEGREE")
    (v === nothing || isempty(v)) && return Int32(1)
    return round(Int32, clamp(v[1], -1e9, 1e9), RoundNearestTiesAway)
end

"""
    _phase_ncoef(br, body) -> Int

Number of polynomial coefficients per nutation-precession angle: `br.phase_degree + 1`.
Throws `SpiceError(:degreeoutofrange)` unless the degree is in 1:3.
"""
@inline function _phase_ncoef(br::BodyRot, body::Integer)
    deg = Int(br.phase_degree)
    1 <= deg <= 3 || throw(SpiceError(:degreeoutofrange, "maximum phase angle degree " *
        "for body $(_bodbry(body)) must be in the range 1:3 but was $deg"))
    return deg + 1
end

"Value at `x` of the polynomial with coefficients `c` (lowest order first) by Horner's rule."
@inline function _horner(c::Vector{Float64}, x::T) where {T<:Real}
    n = length(c)
    n == 0 && return zero(T)
    acc = zero(T) + @inbounds c[n]
    for k in (n - 1):-1:1
        acc = @inbounds(c[k]) + x * acc
    end
    return acc
end

"""
    _text_pck_rotation(br, body, et::T) -> Mat3{T}

Matrix from J2000 to the body-fixed frame of `body` from the text-PCK constants `br`,
as `tisbod_` of CSPICE computes it (same order of operations): RA, DEC and W
polynomials, trigonometric terms, `eul2m(W, pi/2 - DEC, pi/2 + RA, 3, 1, 3)`, then the
rotation from J2000 to the reference frame of the constants. Throws
`SpiceError(:insufficientangles)` if an amplitude vector is longer than the angle list.
"""
function _text_pck_rotation(br::BodyRot, body::Integer, et::T) where {T<:Real}
    epoch = et - _SEC_PER_DAY * (br.jed_epoch - _J2000_JED)
    td = epoch / _SEC_PER_DAY
    tc = epoch / _SEC_PER_CENTURY
    ra = _horner(br.pole_ra, tc)
    dec = _horner(br.pole_dec, tc)
    w = _horner(br.pm, td)

    ac, dc, wc, th = br.nut_prec_ra, br.nut_prec_dec, br.nut_prec_pm, br.nut_prec_angles
    na, nd, nw = length(ac), length(dc), length(wc)
    nterm = max(na, nd, nw)
    if nterm > 0
        nco = _phase_ncoef(br, body)
        nphase = length(th) ÷ nco
        nterm <= nphase || throw(SpiceError(:insufficientangles,
            "insufficient number of nutation/precession angles for body $body: " *
            "number of angles is $nphase, number required is $nterm"))
        sa = zero(T)
        sd = zero(T)
        sw = zero(T)
        @inbounds for i in 1:nterm    # invariant: sa, sd, sw hold the sums of terms 1:i-1
            k = nco * (i - 1)
            if nco == 2
                theta = (th[k+1] + tc * th[k+2]) * _RPD
            else
                theta = zero(T) + th[k+1]
                p = tc
                for j in 2:nco
                    theta += p * th[k+j]
                    p = tc * p
                end
                theta *= _RPD
            end
            s, c = sincos(theta)
            i <= na && (sa += ac[i] * s)
            i <= nd && (sd += dc[i] * c)
            i <= nw && (sw += wc[i] * s)
        end
        ra += sa
        dec += sd
        w += sw
    end

    ra *= _RPD
    dec *= _RPD
    w = _reduce_twopi(w * _RPD)
    tipm = _eul313(w, _HALFPI - dec, ra + _HALFPI)
    br.ref == _J2000_FRAME && return tipm
    return tipm * inertial_rotation(br.ref)
end

"""
    _binary_pck_rotation(data, seg, et::T) -> Mat3{T}

Matrix from J2000 to the body frame of the binary PCK segment `seg` at `et`, as
`pcke02_`, `pckmat_` and `tisbod_` give it: 3-1-3 Euler angles of the segment (W
reduced by whole turns), then the rotation
from J2000 to the reference frame of the segment. Requires `seg.t0 <= primal(et) <= seg.t1`.
"""
@inline function _binary_pck_rotation(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    ang = pck_eval(data, seg, et, Val(0))[1]
    # `pcke02_` reduces W with `d_mod`; the same arithmetic keeps the value of CSPICE.
    tipm = _eul313(_reduce_twopi(ang[3]), ang[2], ang[1])
    seg.frame == _J2000_FRAME && return tipm
    return tipm * inertial_rotation(seg.frame)
end

"""
    _body_rotation(ks, classid, et::T) -> Mat3{T}

Matrix from J2000 to the class 2 frame with class ID `classid`: from a binary PCK
segment if one covers `primal(et)`, otherwise from the text-PCK constants, as
`tisbod_` does. Throws `SpiceError(:framedatanotfound)` if neither exists.
"""
function _body_rotation(ks::KernelSet, classid::Int32, et::T) where {T<:Real}
    segs = get(ks.pck, classid, nothing)
    if segs !== nothing
        i = find_segment(segs, primal(et))
        if i > 0
            seg = @inbounds segs[i]
            return _binary_pck_rotation(ks.daf[seg.file], seg, et)::Mat3{T}
        end
    end
    br = get(ks.bodyrot, classid, nothing)
    br === nothing && throw(SpiceError(:framedatanotfound,
        "PCK data required to compute the orientation of the body-fixed frame with " *
        "class ID $classid at ET $(primal(et)) were not found (no binary PCK segment " *
        "covers the epoch and no text PCK constants are loaded)"))
    return _text_pck_rotation(br, classid, et)::Mat3{T}
end

# --- fixed-offset (TK) frames ---------------------------------------------------------

"""
    _sharpr(m) -> Mat3{Float64}

`m` with its columns made orthonormal as `sharpr_` of CSPICE does: column 1 scaled to
unit length, column 3 the unit cross product of columns 1 and 2, column 2 the unit
cross product of columns 3 and 1.
"""
function _sharpr(m::Mat3{Float64})
    unit(v) = (n = norm(v); n > 0 ? v / n : v)
    c1 = unit(m[:, 1])
    c3 = unit(cross(c1, m[:, 2]))
    c2 = unit(cross(c3, c1))
    return hcat(c1, c2, c3)
end

"""
    _q2m(q) -> Mat3{Float64}

Rotation matrix of the SPICE-style quaternion `q` (scalar first), as `q2m_` of CSPICE
after the normalization that `tkfram_` applies (`vhatg_`). Requires `q` of length 4.
"""
function _q2m(q::AbstractVector{Float64})
    n = sqrt(q[1]^2 + q[2]^2 + q[3]^2 + q[4]^2)
    q0, q1, q2, q3 = n > 0 ? (q[1] / n, q[2] / n, q[3] / n, q[4] / n) : (q[1], q[2], q[3], q[4])
    q01, q02, q03 = q0 * q1, q0 * q2, q0 * q3
    q12, q13, q23 = q1 * q2, q1 * q3, q2 * q3
    q1s, q2s, q3s = q1 * q1, q2 * q2, q3 * q3
    l2 = q0 * q0 + q1s + q2s + q3s
    if l2 != 1.0 && l2 != 0.0
        f = 1.0 / l2
        q01 *= f; q02 *= f; q03 *= f
        q12 *= f; q13 *= f; q23 *= f
        q1s *= f; q2s *= f; q3s *= f
    end
    return Mat3{Float64}(
        1.0 - 2.0 * (q2s + q3s), 2.0 * (q12 + q03), 2.0 * (q13 - q02),
        2.0 * (q12 - q03), 1.0 - 2.0 * (q1s + q3s), 2.0 * (q23 + q01),
        2.0 * (q13 + q02), 2.0 * (q23 - q01), 1.0 - 2.0 * (q1s + q2s))
end

"Factor that converts an angle in `units` (a `convrt_` angle unit name) to radians."
function _angle_factor(units::AbstractString)
    u = uppercase(strip(units))
    table = ("RADIANS" => 1.0, "DEGREES" => pi / 180, "ARCMINUTES" => pi / 10800,
        "ARCSECONDS" => pi / 648000, "HOURANGLE" => pi / 12, "MINUTEANGLE" => pi / 720,
        "SECONDANGLE" => pi / 43200)
    for (name, f) in table
        u == name && return f
    end
    throw(SpiceError(:badunits, "unit '$units' of a TK frame angle is not an angle unit"))
end

"""
    _tk_item(pool, idkey, namekey, suffix)

Value of the TK frame variable `suffix`, from the form keyed by frame name if it
exists and otherwise from the form keyed by ID; `nothing` if neither exists. Throws
`SpiceError(:competingframespec)` if both exist, as `tkfram_` does.
"""
function _tk_item(pool::AbstractDict, idkey::AbstractString, namekey::AbstractString,
    suffix::AbstractString)
    a = get(pool, idkey * suffix, nothing)
    b = namekey == idkey ? nothing : get(pool, namekey * suffix, nothing)
    a !== nothing && b !== nothing && throw(SpiceError(:competingframespec,
        "frame definition keywords '$idkey$suffix' and '$namekey$suffix' are both present in the pool"))
    return b !== nothing ? b : a
end

"""
    _tk_relative_id(pool, frames, name) -> Union{Int32,Nothing}

ID of the frame called `name`: a frame of `frames` (defined in the pool), a
`FRAME_<name>` pool variable, a built-in frame of `frames.jl`, or `nothing`.
"""
function _tk_relative_id(pool::AbstractDict, frames::AbstractDict, name::AbstractString)
    key = normalize_name(name)
    for fd in values(frames)
        normalize_name(fd.name) == key && return fd.id
    end
    v = _pool_numbers(pool, "FRAME_" * key)
    v !== nothing && !isempty(v) && return round(Int32, v[1])
    return get(BUILTIN_FRAMENAME2ID, key, nothing)
end

"""
    _tk_matrix(pool, idkey, namekey) -> Mat3{Float64}

Constant matrix of one TK frame (from the frame to its relative frame) for the
specification `MATRIX`, `ANGLES` or `QUATERNION`, as `tkfram_` of CSPICE builds it.
Throws `SpiceError` on a missing or malformed variable.
"""
function _tk_matrix(pool::AbstractDict, idkey::AbstractString, namekey::AbstractString)
    bad(what) = throw(SpiceError(:badtkframe, "TK frame $idkey: $what"))
    spec = _tk_item(pool, idkey, namekey, "_SPEC")
    (spec isa Vector{String} && length(spec) == 1) || bad("_SPEC must be one string")
    kind = uppercase(strip(spec[1]))
    if kind == "MATRIX"
        v = _tk_item(pool, idkey, namekey, "_MATRIX")
        (v isa Vector{Float64} && length(v) == 9) || bad("_MATRIX must have 9 numbers")
        m = Mat3{Float64}(v)                    # stored column by column
        r = _sharpr(m)
        c2 = dot(r[:, 2], m[:, 2]) < 0 ? -r[:, 2] : r[:, 2]
        c3 = dot(r[:, 3], m[:, 3]) < 0 ? -r[:, 3] : r[:, 3]
        return hcat(r[:, 1], c2, c3)
    elseif kind == "ANGLES"
        ang = _tk_item(pool, idkey, namekey, "_ANGLES")
        axes = _tk_item(pool, idkey, namekey, "_AXES")
        units = _tk_item(pool, idkey, namekey, "_UNITS")
        (ang isa Vector{Float64} && length(ang) == 3) || bad("_ANGLES must have 3 numbers")
        (axes isa Vector{Float64} && length(axes) == 3) || bad("_AXES must have 3 numbers")
        f = units isa Vector{String} && !isempty(units) ? _angle_factor(units[1]) : 1.0
        ax = round.(Int, axes)
        return _eul2m(ang[1] * f, ang[2] * f, ang[3] * f, ax[1], ax[2], ax[3])
    elseif kind == "QUATERNION"
        q = _tk_item(pool, idkey, namekey, "_Q")
        (q isa Vector{Float64} && length(q) == 4) || bad("_Q must have 4 numbers")
        return _q2m(q)
    end
    bad("specification '$(spec[1])' is not MATRIX, ANGLES or QUATERNION")
end

"""
    build_tk(pool, frames) -> Dict{Int32,TKFrame}

Fixed-offset data of every class 4 frame of `frames` and of the built-in table, keyed
by frame ID: the relative
frame (`TKFRAME_<classid>_RELATIVE`, or the same variable keyed by frame name) and the
constant matrix from the frame to its relative frame, as `tkfram_` of CSPICE. A frame
whose definition is missing or malformed is left out, so that loading a kernel never
fails; `frame_rotation` then reports that the frame has no data.
"""
function build_tk(pool::AbstractDict, frames::AbstractDict)
    out = Dict{Int32,TKFrame}()
    builtin = (fd for fd in values(BUILTIN_FRAMES) if !haskey(frames, fd.id))
    for fd in Iterators.flatten((values(frames), builtin))
        fd.class == 4 || continue
        idkey = "TKFRAME_$(fd.classid)"
        namekey = "TKFRAME_" * fd.name
        try
            rel = _tk_item(pool, idkey, namekey, "_RELATIVE")
            (rel isa Vector{String} && length(rel) == 1) || continue
            relid = _tk_relative_id(pool, frames, rel[1])
            (relid === nothing || relid == fd.id) && continue
            out[fd.id] = TKFrame(relid, _tk_matrix(pool, idkey, namekey))
        catch err
            err isa SpiceError || rethrow()
        end
    end
    return out
end

# --- frame chains -----------------------------------------------------------------------

"""
    _frame_step(ks, id, et::T) -> (parent::Int32, rot::Mat3{T})

One link of the frame chain: the frame that `id` is defined relative to and the matrix
from frame `id` to that parent at `et`. Requires `id != 1`. Classes 3 and 5 throw
`SpiceError(:unsupported)`; a class 4 frame without data throws
`SpiceError(:framedatanotfound)`.
"""
function _frame_step(ks::KernelSet, id::Int32, et::T) where {T<:Real}
    fd = frame_def(ks, id)
    class = fd.class
    if class == 1
        return _J2000_FRAME, Mat3{T}(transpose(inertial_rotation(id)))
    elseif class == 2
        return _J2000_FRAME, Mat3{T}(transpose(_body_rotation(ks, fd.classid, et)))
    elseif class == 4
        tk = get(ks.tk, id, nothing)
        tk === nothing && throw(SpiceError(:framedatanotfound,
            "no valid TKFRAME_* definition is loaded for the fixed-offset frame $id"))
        return tk.relative, Mat3{T}(tk.rot)
    end
    throw(SpiceError(:unsupported,
        "frame $id has class $class; only classes 1, 2 and 4 are supported"))
end

"""
    _frame_parent(ks, id) -> Int32

The frame that `id != 1` is defined relative to, without any evaluation: the relative
frame of a class 4 frame, J2000 for classes 1 and 2. Throws as `_frame_step` does.
"""
function _frame_parent(ks::KernelSet, id::Int32)
    class = frame_def(ks, id).class
    (class == 1 || class == 2) && return _J2000_FRAME
    if class == 4
        tk = get(ks.tk, id, nothing)
        tk === nothing && throw(SpiceError(:framedatanotfound,
            "no valid TKFRAME_* definition is loaded for the fixed-offset frame $id"))
        return tk.relative
    end
    throw(SpiceError(:unsupported,
        "frame $id has class $class; only classes 1, 2 and 4 are supported"))
end

@noinline _chain_too_long(id) = throw(SpiceError(:framechain,
    "the chain of frame $id does not reach J2000 in $MAX_FRAME_CHAIN steps"))

"""
    _rotation_up(ks, id, stop, et::T) -> Mat3{T}

Matrix from frame `id` to its ancestor `stop` at `et`. Requires `stop` to be `id`
itself or a frame on the chain from `id` to J2000.
"""
function _rotation_up(ks::KernelSet, id::Int32, stop::Int32, et::T) where {T<:Real}
    rot = one(Mat3{T})
    cur = id
    n = 0
    while cur != stop      # invariant: rot takes vectors from frame id to frame cur
        n += 1
        n > MAX_FRAME_CHAIN && _chain_too_long(id)
        cur, step = _frame_step(ks, cur, et)
        rot = n == 1 ? step : step * rot
    end
    return rot
end

"""
    frame_rotation(ks, from, to, et::T) -> Mat3{T}

Matrix that rotates vectors from frame `from` to frame `to` at `et` (ET seconds),
generic over the epoch type. The two frame chains (class 4 to its relative frame,
classes 1 and 2 to J2000) are followed to their first common frame, as `refchg_` of
CSPICE does, and at most `MAX_FRAME_CHAIN` links each. `from == to` gives the identity.
Classes 3 and 5 throw `SpiceError(:unsupported)`.
"""
function frame_rotation(ks::KernelSet, from::Int32, to::Int32, et::T) where {T<:Real}
    from == to && return one(Mat3{T})
    # The chain of `to`, without evaluation: to, parent(to), ..., J2000.
    chain = ntuple(_ -> Int32(0), Val(MAX_FRAME_CHAIN + 1))
    cur = to
    len = 1
    chain = Base.setindex(chain, cur, 1)
    while cur != _J2000_FRAME    # invariant: chain[1:len] is the chain of `to` down to cur
        len > MAX_FRAME_CHAIN && _chain_too_long(to)
        cur = _frame_parent(ks, cur)
        len += 1
        chain = Base.setindex(chain, cur, len)
    end
    # The first frame on the chain of `from` that is also on the chain of `to`.
    node = from
    n = 0
    while !_in_chain(chain, len, node)
        n += 1
        n > MAX_FRAME_CHAIN && _chain_too_long(from)
        node = _frame_parent(ks, node)
    end
    up = _rotation_up(ks, from, node, et)       # from -> node
    node == to && return up
    down = _rotation_up(ks, to, node, et)       # to -> node
    node == from && return Mat3{T}(transpose(down))
    return Mat3{T}(transpose(down) * up)
end

frame_rotation(ks::KernelSet, from::Integer, to::Integer, et::Real) =
    frame_rotation(ks, Int32(from), Int32(to), et)

"True if `id` is one of the first `len` entries of `chain`."
@inline function _in_chain(chain::NTuple{M,Int32}, len::Int, id::Int32) where {M}
    for k in 1:len
        @inbounds chain[k] == id && return true
    end
    return false
end

"""
    rotation_jet(ks, from, to, et::Real, Val(N)) -> NTuple{N+1,Mat3}

The matrix from frame `from` to frame `to` at `et` and its ordinary time derivatives
of orders 1 to N (element `k+1` is the derivative of order `k`, per second^k).
"""
@inline rotation_jet(ks::KernelSet, from::Integer, to::Integer, et::Real, ::Val{0}) =
    (frame_rotation(ks, Int32(from), Int32(to), et),)

function rotation_jet(ks::KernelSet, from::Integer, to::Integer, et::Real, v::Val{N}) where {N}
    N isa Integer && N >= 0 || throw(ArgumentError("derivative order must be a non-negative integer"))
    r = frame_rotation(ks, Int32(from), Int32(to), taylor_var(et, v))
    d = map(derivatives, r)
    return ntuple(k -> map(x -> x[k], d), Val(N + 1))
end
