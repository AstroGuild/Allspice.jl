# Frame interface: rotation matrices between reference frames and their derivatives.

"""
    pxform(from, to, et) -> SMatrix{3,3}

Return the matrix that rotates position vectors from frame `from` to frame `to` at ET
`et`, as `pxform_c`. Frames are names or integer IDs. Forms: on a `KernelSet`, on a
`SpiceContext`, and on the default context.
"""
pxform(ks::KernelSet, from::NameOrID, to::NameOrID, et::Real) =
    frame_rotation(ks, frame_id(ks, from), frame_id(ks, to), check_epoch(et))
pxform(ctx::SpiceContext, from::NameOrID, to::NameOrID, et::Real) =
    pxform(kernelset(ctx), from, to, et)
pxform(from::NameOrID, to::NameOrID, et::Real) = pxform(default_context(), from, to, et)

"""
    pxder(from, to, et, Val(N)) -> NTuple{N+1,SMatrix{3,3}}

Return the rotation matrix from frame `from` to frame `to` at ET `et` and its time
derivatives to order N.

For a chain with a CK frame (class 3) the function that is differentiated is the
rotation that `pxform` returns: the interpolation of the CK segment (constant for type
1 and at single pointing instances, a rotation at the constant stored rate for type 2,
a rotation about a fixed axis that is linear in the clock for type 3, normalized
Chebyshev, Hermite or Lagrange polynomials of the quaternion components for types 4, 5
and 6), composed with the piecewise linear clock conversion `sce2c`. The stored angular
velocity is not used, so order 1 can differ from the derivative block of `sxform`,
which is what CSPICE returns.

For a chain with a dynamic frame (class 5) the function that is differentiated is also
the rotation that `pxform` returns, the one of `zzdynrot_`:

- of-date frames: the precession, nutation and obliquity series as functions of `et`
  (all terms, also the slow drift of the nutation coefficients that the rates of
  `zzwahr_` omit);
- two-vector frames: the positions of `spkpos`, the unit velocities of `spkezr`, the
  near points and the rotated constant vectors, each with its aberration correction as
  a function of `et` (light-time corrections `LT` with the one iteration of CSPICE, `CN`
  with the converged light-time equation, stellar aberration with the velocity of the
  observer as a function of `et`);
- Euler frames: the polynomials of the angles;
- product frames: the product of the factors.

A frozen frame (`FREEZE_EPOCH`) is constant, so its derivatives relative to its base
are zero. `ROTATION_STATE = 'INERTIAL'` does not change the rotation, so it does not
change `pxder`: the derivatives are those of the rotating frame, while `sxform` gives a
zero derivative block for such a frame. A velocity that is given in a non-inertial
frame takes the frame derivative of `sxform`, as `spkezr` does. For a switch frame
(class 6) the base that is selected at `et` is differentiated; at an epoch at which
the selection changes the result is the one of the base that `pxform` selects there.
"""
function pxder(ks::KernelSet, from::NameOrID, to::NameOrID, et::Real, ::Val{N}) where {N}
    N isa Int && N >= 0 || throw(SpiceError(:badarg, "pxder: the order must be an integer >= 0"))
    return _exact_frame_jet(ks, frame_id(ks, from), frame_id(ks, to), check_epoch(et), Val(N))
end
pxder(ctx::SpiceContext, from::NameOrID, to::NameOrID, et::Real, v::Val) =
    pxder(kernelset(ctx), from, to, et, v)
pxder(from::NameOrID, to::NameOrID, et::Real, v::Val) =
    pxder(default_context(), from, to, et, v)

"""
    sxform(from, to, et) -> SMatrix{6,6}

Return the state transformation matrix `[R 0; dR R]` from frame `from` to frame `to`
at ET `et`, as `sxform_c`; `dR` is the time derivative of the rotation matrix `R`.

For a chain with a CK frame (class 3) `dR` is built from the angular velocity of the CK
segment, as CSPICE does, and only segments with angular velocity are used; see `pxder`
for the derivative of the interpolated rotation.

For a chain with a dynamic frame (class 5) or a switch frame (class 6) `dR` is the
derivative block of `zzdynfrm_`, with the choices of CSPICE: velocities of `spkezr`
for the position vectors, a central difference of the unit velocity over
`±max(1, et * 2^-27)` seconds for the velocity vectors, the rates of `zzwahr_` for the
nutation, a zero derivative for constant vectors and for their stellar aberration, a
zero block for a frozen frame and for `ROTATION_STATE = 'INERTIAL'`. It can therefore
differ from order 1 of `pxder`, which is the exact derivative of `R`.
"""
function sxform(ks::KernelSet, from::NameOrID, to::NameOrID, et::Real)
    R, dR = frame_jet(ks, frame_id(ks, from), frame_id(ks, to), check_epoch(et), Val(1))
    return vcat(hcat(R, zero(R)), hcat(dR, R))
end
sxform(ctx::SpiceContext, from::NameOrID, to::NameOrID, et::Real) =
    sxform(kernelset(ctx), from, to, et)
sxform(from::NameOrID, to::NameOrID, et::Real) = sxform(default_context(), from, to, et)
