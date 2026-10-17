# Frame interface: rotation matrices between reference frames, at one epoch or at two
# epochs, and their derivatives.

"""
    pxform([ctx,] from, to, et) -> SMatrix{3,3}

Return the matrix that rotates a position vector from frame `from` to frame `to` at
`et`, as `pxform_c`.

For the matrix `R`, `R * v` gives in `to` the components of a vector `v` given in
`from`. `from` and `to` are frame names or frame IDs. `et` is in seconds past J2000
TDB. The frame classes of CSPICE are supported: inertial, PCK, CK, fixed offset,
dynamic and switch frames. Two inertial frames need no kernel. With ForwardDiff
loaded, `et` can be a dual number.

Throws `SpiceError`:
- `:unknownframe` if a frame name is not known;
- `:badarg` if `et` is not finite;
- `:framedatanotfound` if the definition or the constants of a frame are not loaded
  (the `BODY<id>_PM` variables of a PCK frame, the matrix of a fixed-offset frame);
- `:noframeconnect` if a CK frame has no pointing or no clock kernel at `et`, or a
  switch frame has no base there.

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
    pxfrm2([ctx,] from, to, etfrom, etto) -> SMatrix{3,3}

Return the matrix that rotates a position vector from frame `from` at the epoch
`etfrom` to frame `to` at the epoch `etto`, as `pxfrm2_c` of CSPICE.

The matrix is the product of the rotation from `from` to J2000 at `etfrom` and the
rotation from J2000 to `to` at `etto`. It is the rotation that the light-time
corrections need, where the frames of the observer and of the target are evaluated
at different epochs. The arguments and the errors are those of [`pxform`](@ref); with
`etfrom == etto` the result is the matrix of `pxform` to rounding.

# Examples
```jldoctest
julia> R = pxfrm2(SpiceContext(), "J2000", "ECLIPJ2000", 0.0, 1.0e9);

julia> R == pxform(SpiceContext(), "J2000", "ECLIPJ2000", 0.0)    # inertial frames
true
```

See also [`pxform`](@ref), [`spkcpt`](@ref).
"""
function pxfrm2(ks::KernelSet, from::NameOrID, to::NameOrID, etfrom::Real, etto::Real)
    fcode = frame_id(ks, from)
    tcode = frame_id(ks, to)
    jf = frame_rotation(ks, fcode, Int32(1), check_epoch(etfrom))
    tj = frame_rotation(ks, Int32(1), tcode, check_epoch(etto))
    return _m_mxm(tj, jf)
end
pxfrm2(ctx::SpiceContext, from::NameOrID, to::NameOrID, etfrom::Real, etto::Real) =
    pxfrm2(kernelset(ctx), from, to, etfrom, etto)
pxfrm2(from::NameOrID, to::NameOrID, etfrom::Real, etto::Real) =
    pxfrm2(default_context(), from, to, etfrom, etto)

"""
    pxder([ctx,] from, to, et, Val(N)) -> NTuple{N+1,SMatrix{3,3}}

Return the rotation matrix from frame `from` to frame `to` at `et` and its time
derivatives to order `N`. This call has no CSPICE counterpart.

Element `k + 1` of the result is the derivative of order `k`, in 1/s^k. The
arguments `from`, `to` and `et` are those of [`pxform`](@ref). Order 0 is the matrix of
`pxform`, and every higher order is the exact derivative of that matrix as a function
of `et`.

Order 1 is the derivative block of [`sxform`](@ref) for inertial, PCK and fixed-offset
frames. It can differ for these frames, where `sxform` follows the rates of CSPICE:

- CK frames (class 3);
- dynamic frames (class 5), also a frozen frame and one with
  `ROTATION_STATE = 'INERTIAL'`;
- switch frames (class 6) whose base is one of the above.

[Exact derivatives and the rates of CSPICE](@ref) says what is differentiated for each
class.

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

`dR` is the rate that CSPICE gives. For a CK frame it comes from the angular velocity
of the CK segment, and only segments with angular velocity are used. For the frame
classes that [`pxder`](@ref) lists, `dR` is not the exact derivative of `R`; use
order 1 of `pxder` for that.

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
