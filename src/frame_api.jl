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
"""
function pxder(ks::KernelSet, from::NameOrID, to::NameOrID, et::Real, ::Val{N}) where {N}
    N isa Int && N >= 0 || throw(SpiceError(:badarg, "pxder: the order must be an integer >= 0"))
    return frame_jet(ks, frame_id(ks, from), frame_id(ks, to), check_epoch(et), Val(N))
end
pxder(ctx::SpiceContext, from::NameOrID, to::NameOrID, et::Real, v::Val) =
    pxder(kernelset(ctx), from, to, et, v)
pxder(from::NameOrID, to::NameOrID, et::Real, v::Val) =
    pxder(default_context(), from, to, et, v)

"""
    sxform(from, to, et) -> SMatrix{6,6}

Return the state transformation matrix `[R 0; dR R]` from frame `from` to frame `to`
at ET `et`, as `sxform_c`; `dR` is the time derivative of the rotation matrix `R`.
"""
function sxform(ks::KernelSet, from::NameOrID, to::NameOrID, et::Real)
    R, dR = frame_jet(ks, frame_id(ks, from), frame_id(ks, to), check_epoch(et), Val(1))
    return vcat(hcat(R, zero(R)), hcat(dR, R))
end
sxform(ctx::SpiceContext, from::NameOrID, to::NameOrID, et::Real) =
    sxform(kernelset(ctx), from, to, et)
sxform(from::NameOrID, to::NameOrID, et::Real) = sxform(default_context(), from, to, et)
