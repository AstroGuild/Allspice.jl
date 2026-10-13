# Rotation calls of CSPICE: Euler angles, axis and angle, quaternions, angular velocity
# and state transformations.
#
# A matrix argument or result has the mathematical index order (`m[i, j]` is row `i`,
# column `j`), which is the order of the CSPICE C interface. The f2c sources index a
# matrix by one number in column-major order; a comment gives that index where the
# branch logic of CSPICE reads single elements.

"`acos` that returns NaN outside [-1, 1], as the C library does."
@inline _m_acos(x::Real) = abs(primal(x)) > 1.0 ? x * NaN : acos(x)

"`asin` that returns NaN outside [-1, 1], as the C library does."
@inline _m_asin(x::Real) = abs(primal(x)) > 1.0 ? x * NaN : asin(x)

@noinline _m_badaxes(a, b, c) =
    throw(SpiceError(:badaxisnumbers, "axis numbers are $a, $b, $c; each must be 1, 2 or 3"))

@inline _m_axisok(a::Integer) = 1 <= a <= 3

"Successor of the axis `a` in the cycle 1, 2, 3."
@inline _m_next(a::Integer) = Int(a) % 3 + 1

@inline function _m_eul2m(angle3, angle2, angle1, axis3::Integer, axis2::Integer, axis1::Integer)
    (_m_axisok(axis3) && _m_axisok(axis2) && _m_axisok(axis1)) || _m_badaxes(axis3, axis2, axis1)
    return _m_rotmat(_m_rotmat(_m_rotate(angle1, axis1), angle2, axis2), angle3, axis3)
end

"""
    eul2m(angle3, angle2, angle1, axis3, axis2, axis1) -> SMatrix{3,3}

Return the rotation matrix `[angle3]_axis3 * [angle2]_axis2 * [angle1]_axis1` of three
Euler angles in radians, as `eul2m_c`. Requires each axis in 1:3
(`SpiceError(:badaxisnumbers)` otherwise). Postcondition: the result is orthogonal with
determinant 1.
"""
function eul2m(angle3::Real, angle2::Real, angle1::Real, axis3::Integer, axis2::Integer,
        axis1::Integer)
    a3, a2, a1 = promote(_m_f(angle3), _m_f(angle2), _m_f(angle1))
    return _m_eul2m(a3, a2, a1, axis3, axis2, axis1)
end

"Axis change of `m2eul_`: column `j` has the value `sj` in row `rj` and zeros elsewhere."
@inline function _m_change(r1::Int, s1::Float64, r2::Int, s2::Float64, r3::Int, s3::Float64)
    return SMatrix{3,3,Float64,9}(ntuple(Val(9)) do k
        i, j = _m_rowcol(k)
        j == 1 ? (i == r1 ? s1 : 0.0) : j == 2 ? (i == r2 ? s2 : 0.0) : (i == r3 ? s3 : 0.0)
    end)
end

@inline _m_p0(x) = primal(x) == 0.0

@inline function _m_m2eul(r::SMatrix{3,3}, axis3::Integer, axis2::Integer, axis1::Integer)
    (_m_axisok(axis3) && _m_axisok(axis2) && _m_axisok(axis1)) || _m_badaxes(axis3, axis2, axis1)
    (axis3 == axis2 || axis1 == axis2) && throw(SpiceError(:badaxisnumbers,
        "the middle axis matches a neighbour: $axis3 $axis2 $axis1"))
    _m_isrot(r, 0.1, 0.1) || throw(SpiceError(:notarotation, "the input matrix is not a rotation"))
    t0 = _m_cols(_m_hat(_m_col(r, 1)), _m_hat(_m_col(r, 2)), _m_hat(_m_col(r, 3)))
    sign = axis2 == _m_next(axis3) ? 1.0 : -1.0
    if axis3 == axis1
        # a-b-a sequence: the factorization of a 3-1-3 rotation after the axis change
        c = 6 - Int(axis3) - Int(axis2)
        change = _m_change(Int(axis2), 1.0, c, sign, Int(axis3), 1.0)
        t = _m_mtxm(change, _m_mxm(t0, change))
        # f2c elements 6, 7 are (1,3), (2,3); 2, 5 are (3,1), (3,2); 8 is (3,3)
        degen = (_m_p0(t[7]) && _m_p0(t[8])) || (_m_p0(t[3]) && _m_p0(t[6])) ||
            abs(primal(t[9])) == 1.0
        angle2 = _m_acos(t[9])
        if degen
            return zero(angle2), angle2, atan(t[4], t[1]) - zero(angle2)
        end
        return atan(t[7], t[8]) - zero(angle2), angle2, atan(t[3], -t[6]) - zero(angle2)
    end
    # a-b-c sequence: the factorization of a 1-2-3 rotation after the axis change
    change = _m_change(Int(axis3), 1.0, Int(axis2), 1.0, Int(axis1), sign)
    t = _m_mtxm(change, _m_mxm(t0, change))
    # f2c elements 0, 3 are (1,1), (1,2); 7, 8 are (2,3), (3,3); 6 is (1,3)
    degen = (_m_p0(t[1]) && _m_p0(t[4])) || (_m_p0(t[8]) && _m_p0(t[9])) ||
        abs(primal(t[7])) == 1.0
    angle2 = _m_asin(-t[7])
    if degen
        return zero(angle2), angle2, sign * atan(-t[2], t[5]) - zero(angle2)
    end
    return atan(t[8], t[9]) - zero(angle2), angle2, sign * atan(t[4], t[1]) - zero(angle2)
end

"""
    m2eul(r, axis3, axis2, axis1) -> (angle3, angle2, angle1)

Return the Euler angles in radians that factor the rotation matrix `r` as
`[angle3]_axis3 * [angle2]_axis2 * [angle1]_axis1`, as `m2eul_c`. `angle3` and `angle1`
are in (-pi, pi]; `angle2` is in [0, pi] if `axis3 == axis1` and in [-pi/2, pi/2]
otherwise. In the degenerate case (`angle2` at an end of its range) `angle3` is 0, with
the test of CSPICE. Requires axes in 1:3 with the middle one different from its
neighbours (`SpiceError(:badaxisnumbers)`) and a rotation matrix within the tolerance
0.1 of `isrot` (`SpiceError(:notarotation)`).
"""
m2eul(r::AbstractMatrix, axis3::Integer, axis2::Integer, axis1::Integer) =
    _m_m2eul(_m_m3(r), axis3, axis2, axis1)

"The factor `delta(a, b)` of `xf2eul_`: 1 if `b` follows `a` in the cycle, else -1."
@inline _m_delta(a::Integer, b::Integer) = b == _m_next(a) ? 1.0 : -1.0

"""
    eul2xf(eulang, axisa, axisb, axisc) -> SMatrix{6,6}

Return the state transformation of the Euler angles `eulang[1:3]` (for the axes
`axisa`, `axisb`, `axisc`, in the order of the matrix product) and their rates
`eulang[4:6]`, as `eul2xf_c`. Requires 6 elements and axes in 1:3
(`SpiceError(:badaxisnumbers)`). Equal neighbour axes are accepted, as in CSPICE: the
two angles add. Postcondition: the upper left and lower right blocks are
`eul2m(eulang[1], eulang[2], eulang[3], axisa, axisb, axisc)`.
"""
function eul2xf(eulang::AbstractVector, axisa::Integer, axisb::Integer, axisc::Integer)
    e = _m_v6(eulang)
    (_m_axisok(axisa) && _m_axisok(axisb) && _m_axisok(axisc)) || _m_badaxes(axisa, axisb, axisc)
    a1, a2, a3, d1, d2, d3 = e
    locaxb = Int(axisb)
    if axisb == axisa || axisb == axisc
        if axisb == axisa
            a1 += a2
            d1 += d2
        else
            a3 += a2
            d3 += d2
        end
        a2 = zero(a2)
        d2 = zero(d2)
        locaxb = axisc == _m_next(axisa) ? _m_next(axisc) : _m_next(axisa)
    end
    r = _m_eul2m(a1, a2, a3, axisa, locaxb, axisc)
    a = Int(axisa)
    b = locaxb
    l = 6 - a - b
    d = _m_delta(a, b)
    sa, ca = sincos(a1)
    s2, c2 = sincos(a2)
    u, v = axisa == axisc ? (c2, d * s2) : (-d * s2, c2)
    # domega = solutn * rates, with the matrix `solutn` of `eul2xf_`
    w1 = -d * d1 + (-d * u) * d3
    w2 = (-d * ca) * d2 + (-sa * v) * d3
    w3 = sa * d2 + (-d * ca * v) * d3
    z = zero(w1)
    drdtrt = SMatrix{3,3}(ntuple(Val(9)) do k
        i, j = _m_rowcol(k)
        (i == l && j == b) ? w1 : (i == b && j == l) ? -w1 :
        (i == a && j == l) ? w2 : (i == l && j == a) ? -w2 :
        (i == b && j == a) ? w3 : (i == a && j == b) ? -w3 : z
    end)
    return _m_xform(r, _m_mxm(drdtrt, r))
end

"""
    xf2eul(xform, axisa, axisb, axisc) -> (eulang::SVector{6}, unique::Bool)

Return the Euler angles and their rates of the state transformation `xform`, as
`xf2eul_c`; `unique` is false in the degenerate case, where the first rate is set to 0.
Requirements and errors as `m2eul` for the rotation block.
"""
function xf2eul(xform::AbstractMatrix, axisa::Integer, axisb::Integer, axisc::Integer)
    r, drdt = _m_blocks(_m_m6(xform))
    e1, e2, e3 = _m_m2eul(r, axisa, axisb, axisc)
    a = Int(axisa)
    b = Int(axisb)
    l = 6 - a - b
    d = _m_delta(a, b)
    drdtrt = _m_mxmt(drdt, r)
    o1 = d * drdtrt[l, b]
    o2 = d * drdtrt[a, l]
    o3 = d * drdtrt[b, a]
    sa, ca = sincos(e1)
    s2, c2 = sincos(e2)
    u, v = axisa == axisc ? (c2, d * s2) : (-d * s2, c2)
    abs(primal(u)) == 1.0 && (v = zero(v))
    if primal(v) == 0.0
        r5 = abs(primal(ca)) > abs(primal(sa)) ? -o2 / ca : d * o3 / sa
        return SVector(e1 - zero(o1), e2 - zero(o1), e3 - zero(o1), zero(o1), r5, -u * o1), false
    end
    r4 = -o1 + (u * d * sa / v) * o2 + (u * ca / v) * o3
    r5 = -ca * o2 + (d * sa) * o3
    r6 = (-d * sa / v) * o2 + (-ca / v) * o3
    return SVector(e1 - zero(r4), e2 - zero(r4), e3 - zero(r4), r4, r5, r6), true
end

"""
    axisar(axis, angle) -> SMatrix{3,3}

Return the matrix that rotates vectors by `angle` radians about the vector `axis`
(right-hand rule), as `axisar_c`: column `i` is `vrotv` of unit vector `i`.
Postcondition: a zero `axis` gives the identity.
"""
function axisar(axis::AbstractVector, angle::Real)
    ax = _m_v3(axis)
    th = _m_f(angle)
    return _m_cols(_m_vrotv(SVector(1.0, 0.0, 0.0), ax, th), _m_vrotv(SVector(0.0, 1.0, 0.0), ax, th),
        _m_vrotv(SVector(0.0, 0.0, 1.0), ax, th))
end

@inline function _m_m2q(r::SMatrix{3,3})
    _m_isrot(r, 0.1, 0.1) || throw(SpiceError(:notarotation, "the input matrix is not a rotation"))
    trace = r[1] + r[5] + r[9]
    mtrace = 1.0 - trace
    cc4 = 1.0 + trace
    s114 = mtrace + 2.0 * r[1]
    s224 = mtrace + 2.0 * r[5]
    s334 = mtrace + 2.0 * r[9]
    # The first of the four squares that is at least 1/4 gives the divisor; the order
    # of the tests is the one of `m2q_`.
    if 1.0 <= primal(cc4)
        c = sqrt(cc4 * 0.25)
        factor = 1.0 / (c * 4.0)
        s1 = (r[6] - r[8]) * factor
        s2 = (r[7] - r[3]) * factor
        s3 = (r[2] - r[4]) * factor
    elseif 1.0 <= primal(s114)
        s1 = sqrt(s114 * 0.25)
        factor = 1.0 / (s1 * 4.0)
        c = (r[6] - r[8]) * factor
        s2 = (r[4] + r[2]) * factor
        s3 = (r[7] + r[3]) * factor
    elseif 1.0 <= primal(s224)
        s2 = sqrt(s224 * 0.25)
        factor = 1.0 / (s2 * 4.0)
        c = (r[7] - r[3]) * factor
        s1 = (r[4] + r[2]) * factor
        s3 = (r[8] + r[6]) * factor
    else
        s3 = sqrt(s334 * 0.25)
        factor = 1.0 / (s3 * 4.0)
        c = (r[2] - r[4]) * factor
        s1 = (r[7] + r[3]) * factor
        s2 = (r[8] + r[6]) * factor
    end
    l2 = c * c + s1 * s1 + s2 * s2 + s3 * s3
    if primal(l2) != 1.0
        polish = 1.0 / sqrt(l2)
        c *= polish
        s1 *= polish
        s2 *= polish
        s3 *= polish
    end
    return primal(c) > 0.0 ? SVector(c, s1, s2, s3) : SVector(-c, -s1, -s2, -s3)
end

"""
    m2q(r) -> SVector{4}

Return the unit quaternion in the SPICE convention (scalar part first) of the rotation
matrix `r`, as `m2q_c`. The branch on the trace and on the diagonal is the one of
CSPICE. Requires a rotation matrix within the tolerance 0.1 of `isrot`
(`SpiceError(:notarotation)`). Postcondition: unit norm and a scalar part that is not
negative; `q2m(m2q(r))` is `r` to rounding.
"""
m2q(r::AbstractMatrix) = _m_m2q(_m_m3(r))

"""
    q2m(q) -> SMatrix{3,3}

Return the rotation matrix of the SPICE-convention quaternion `q` (scalar part first),
as `q2m_c`. A quaternion that is not of unit length and not zero is scaled, as in
CSPICE. Contract: `q` has 4 elements.
"""
function q2m(q::AbstractVector)
    v = _m_v4(q)
    q01 = v[1] * v[2]
    q02 = v[1] * v[3]
    q03 = v[1] * v[4]
    q12 = v[2] * v[3]
    q13 = v[2] * v[4]
    q23 = v[3] * v[4]
    q1s = v[2] * v[2]
    q2s = v[3] * v[3]
    q3s = v[4] * v[4]
    l2 = v[1] * v[1] + q1s + q2s + q3s
    if primal(l2) != 1.0 && primal(l2) != 0.0
        sharpn = 1.0 / l2
        q01 *= sharpn
        q02 *= sharpn
        q03 *= sharpn
        q12 *= sharpn
        q13 *= sharpn
        q23 *= sharpn
        q1s *= sharpn
        q2s *= sharpn
        q3s *= sharpn
    end
    return SMatrix{3,3}(
        1.0 - (q2s + q3s) * 2.0, (q12 + q03) * 2.0, (q13 - q02) * 2.0,
        (q12 - q03) * 2.0, 1.0 - (q1s + q3s) * 2.0, (q23 + q01) * 2.0,
        (q13 + q02) * 2.0, (q23 - q01) * 2.0, 1.0 - (q1s + q2s) * 2.0)
end

@inline function _m_qxq(q1::SVector{4}, q2::SVector{4})
    a = SVector(q1[2], q1[3], q1[4])
    b = SVector(q2[2], q2[3], q2[4])
    v = q1[1] * b + q2[1] * a + 1.0 * _m_cross(a, b)
    return SVector(q1[1] * q2[1] - _m_dot(a, b), v[1], v[2], v[3])
end

"""
    qxq(q1, q2) -> SVector{4}

Return the product of two SPICE-convention quaternions, as `qxq_c`. Postcondition:
`q2m(qxq(q1, q2))` is `q2m(q1) * q2m(q2)` to rounding for unit quaternions.
"""
qxq(q1::AbstractVector, q2::AbstractVector) = _m_qxq(_m_v4(q1), _m_v4(q2))

"""
    qdq2av(q, dq) -> SVector{3}

Return the angular velocity of a frame from its SPICE-convention quaternion `q` and the
time derivative `dq`, as `qdq2av_c`: `-2` times the vector part of `conj(q/|q|) * dq`.
The unit of the result is radians per time unit of `dq`.
"""
function qdq2av(q::AbstractVector, dq::AbstractVector)
    qq = _m_v4(q)
    qhat = unormg(qq)[1]
    qstar = SVector(qhat[1], -qhat[2], -qhat[3], -qhat[4])
    t = _m_qxq(qstar, _m_v4(dq))
    return SVector(-2.0 * t[2], -2.0 * t[3], -2.0 * t[4])
end

"""
    rav2xf(rot, av) -> SMatrix{6,6}

Return the state transformation of the rotation matrix `rot` and the angular velocity
`av` of the frame, as `rav2xf_c`. Postcondition: the upper left and lower right blocks
are `rot`, the upper right block is zero.
"""
function rav2xf(rot::AbstractMatrix, av::AbstractVector)
    r = _m_m3(rot)
    w = _m_v3(av)
    z = zero(w[1])
    omegat = SMatrix{3,3}(z, -w[3], w[2], w[3], z, -w[1], -w[2], w[1], z)
    return _m_xform(r, _m_mxm(r, omegat))
end

"""
    xf2rav(xform) -> (rot::SMatrix{3,3}, av::SVector{3})

Return the rotation matrix and the angular velocity of the state transformation
`xform`, as `xf2rav_c`. Postcondition: `rav2xf(xf2rav(xform)...)` is `xform` to
rounding if the upper left block is a rotation.
"""
function xf2rav(xform::AbstractMatrix)
    r, drdt = _m_blocks(_m_m6(xform))
    omega = _m_mtxm(drdt, r)
    return r, SVector(omega[3, 2], omega[1, 3], omega[2, 1])
end

"""
    raxisa(matrix) -> (axis::SVector{3}, angle)

Return the unit axis and the angle in [0, pi] of the rotation `matrix`, as `raxisa_c`.
Postcondition: the identity gives the axis (0, 0, 1) and the angle 0. Requires a
rotation matrix (`SpiceError(:notarotation)`).
"""
function raxisa(matrix::AbstractMatrix)
    q = _m_m2q(_m_m3(matrix))
    v = SVector(q[2], q[3], q[4])
    if _m_iszero(v)
        return SVector(zero(q[1]), zero(q[1]), one(q[1])), zero(atan(q[1]))
    elseif primal(q[1]) == 0.0
        return v, Float64(pi) - zero(atan(q[1]))
    end
    return _m_hat(v), 2.0 * atan(_m_norm(v), q[1])
end

"""
    sharpr(rot) -> SMatrix{3,3}

Return `rot` with orthonormal columns, as `sharpr_` of SPICELIB: column 1 scaled to
unit length, column 3 the unit cross product of columns 1 and 2, column 2 the unit
cross product of columns 3 and 1.
"""
function sharpr(rot::AbstractMatrix)
    m = _m_m3(rot)
    c1 = _m_hat(_m_col(m, 1))
    c3 = _m_ucrss(c1, _m_col(m, 2))
    c2 = _m_ucrss(c3, c1)
    return _m_cols(c1, c2, c3)
end
