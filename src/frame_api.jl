# Frame interface: rotation matrices between reference frames and their derivatives.

"""
    pxform([ctx,] from, to, et) -> SMatrix{3,3}

Return the matrix that rotates a position vector from frame `from` to frame `to` at
`et`, as `pxform_c`.

For the matrix `R`, `R * v` gives in `to` the components of a vector `v` given in
`from`. `from` and `to` are frame names or frame IDs. `et` is in seconds past J2000
TDB. The frame classes of CSPICE are supported: inertial, PCK, CK, fixed offset,
dynamic and switch frames. Two inertial frames need no kernel. With ForwardDiff
loaded, `et` can be a dual number.

Throws `SpiceError(:unknownframe)` if a frame name is not known,
`SpiceError(:badarg)` if `et` is not finite, and `SpiceError(:framedatanotfound)` or
`SpiceError(:noframeconnect)` if the loaded kernels do not hold the orientation data
that connect the two frames at `et`.

# Examples
```jldoctest
julia> R = pxform("J2000", "ECLIPJ2000", 0.0);

julia> round.(R; digits=6)
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1.0   0.0       0.0
 0.0   0.917482  0.397777
 0.0  -0.397777  0.917482
```

With a planetary constants kernel loaded:

```julia
furnsh("pck00011.tpc")
R = pxform("J2000", "IAU_EARTH", et)
r_fixed = R * r_inertial
```

See also [`sxform`](@ref), [`pxder`](@ref).
"""
pxform(ks::KernelSet, from::NameOrID, to::NameOrID, et::Real) =
    frame_rotation(ks, frame_id(ks, from), frame_id(ks, to), check_epoch(et))
pxform(ctx::SpiceContext, from::NameOrID, to::NameOrID, et::Real) =
    pxform(kernelset(ctx), from, to, et)
pxform(from::NameOrID, to::NameOrID, et::Real) = pxform(default_context(), from, to, et)

"""
    pxder([ctx,] from, to, et, Val(N)) -> NTuple{N+1,SMatrix{3,3}}

Return the rotation matrix from frame `from` to frame `to` at `et` and its time
derivatives to order `N`. This call has no CSPICE counterpart.

Element `k + 1` of the result is the derivative of order `k`, in 1/s^k. Order 0 is
the matrix of [`pxform`](@ref), and every higher order is the exact derivative of that
matrix as a function of `et`. Order 1 can therefore differ from the derivative block
of [`sxform`](@ref), which follows the choices of CSPICE:

- CK frames (class 3): the interpolation of the CK segment is differentiated,
  composed with the piecewise linear clock conversion of [`sce2c`](@ref). The
  interpolation is constant for type 1 and at single pointing instances, a rotation
  at the constant stored rate for type 2, a rotation about a fixed axis that is
  linear in the clock for type 3, and normalized Chebyshev, Hermite or Lagrange
  polynomials of the quaternion components for types 4, 5 and 6. The stored angular
  velocity is not used.
- Dynamic frames (class 5), of-date: the precession, nutation and obliquity series as
  functions of `et`, with all terms, also the slow drift of the nutation coefficients
  that the rates of CSPICE omit.
- Dynamic frames, two-vector: the positions of [`spkpos`](@ref), the unit velocities
  of [`spkezr`](@ref), the near points and the rotated constant vectors, each with
  its aberration correction as a function of `et`, as in [`spkder`](@ref).
- Dynamic frames, Euler: the polynomials of the angles. Product frames: the product
  of the factors.
- A frozen frame (`FREEZE_EPOCH`) is constant, so its derivatives relative to its
  base are zero.
- `ROTATION_STATE = 'INERTIAL'` does not change the rotation, so it does not change
  `pxder`: the derivatives are those of the rotating frame, while `sxform` gives a
  zero derivative block for such a frame.
- Switch frames (class 6): the base that is selected at `et` is differentiated. At an
  epoch at which the selection changes, the result is the one of the base that
  `pxform` selects there.

A velocity that a frame definition gives in a frame that is not inertial takes the
frame derivative of `sxform`, as `spkezr` does.

Throws `SpiceError(:badarg)` if `N` is not an integer `>= 0`, and the errors of
[`pxform`](@ref).

# Examples
```jldoctest
julia> R, dR = pxder("J2000", "ECLIPJ2000", 0.0, Val(1));

julia> iszero(dR)    # two inertial frames
true
```

See also [`spkder`](@ref).
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
    sxform([ctx,] from, to, et) -> SMatrix{6,6}

Return the matrix that transforms a state (position, velocity) from frame `from` to
frame `to` at `et`, as `sxform_c`.

The matrix is `[R 0; dR R]`: `R` is the rotation of [`pxform`](@ref) and `dR` is its
time derivative in 1/s. `from` and `to` are frame names or frame IDs; `et` is in
seconds past J2000 TDB.

`dR` follows the choices of CSPICE, so it can differ from order 1 of
[`pxder`](@ref), which is the exact derivative of `R`:

- CK frames (class 3): `dR` is built from the angular velocity of the CK segment, and
  only segments with angular velocity are used.
- Dynamic frames (class 5) and switch frames (class 6): velocities of
  [`spkezr`](@ref) for the position vectors, a central difference of the unit
  velocity over `±max(1, et * 2^-27)` seconds for the velocity vectors, the nutation
  rates of CSPICE, a zero derivative for constant vectors and for their stellar
  aberration, and a zero block for a frozen frame and for
  `ROTATION_STATE = 'INERTIAL'`.

Throws the errors of [`pxform`](@ref).

# Examples
```julia
furnsh("naif0012.tls")
furnsh("pck00011.tpc")
et = str2et("2025-03-01T12:00:00")
X = sxform("J2000", "IAU_MARS", et)
state_fixed = X * state_inertial       # km, km/s
```

See also [`pxform`](@ref), [`pxder`](@ref).
"""
function sxform(ks::KernelSet, from::NameOrID, to::NameOrID, et::Real)
    R, dR = frame_jet(ks, frame_id(ks, from), frame_id(ks, to), check_epoch(et), Val(1))
    return vcat(hcat(R, zero(R)), hcat(dR, R))
end
sxform(ctx::SpiceContext, from::NameOrID, to::NameOrID, et::Real) =
    sxform(kernelset(ctx), from, to, et)
sxform(from::NameOrID, to::NameOrID, et::Real) = sxform(default_context(), from, to, et)
