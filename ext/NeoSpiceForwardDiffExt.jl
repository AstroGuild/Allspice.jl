module NeoSpiceForwardDiffExt

using NeoSpice
using NeoSpice: KernelSet, spk_jet, frame_jet, frame_rotation, rotation_jet, primal
using ForwardDiff
using ForwardDiff: Dual, value, partials
import StaticArrays

# A query with a dual-number epoch does not run the evaluators on dual numbers. It takes
# the analytical time derivatives at the plain epoch (the code that the package image
# already holds) and assembles the dual result from them by the chain rule. The first
# call with a new `Dual` type then compiles only the small assembly functions below, and
# a query costs about as much as `spkder` of the same order.
#
# The partial of a result is therefore the next order of `spkder`/`pxder`. For an SPK
# type with a separate velocity series (type 3) the partial of the position is the
# stored velocity, as in `spkder`, and not the derivative of the position series.
#
# The generic evaluators stay reachable with `invoke` on the `et::Real` signature; the
# tests use this as an independent check of the analytical derivatives.

"`primal(x::Dual)`: the plain `Float64` value of a dual number, with all partials removed."
NeoSpice.primal(x::Dual) = NeoSpice.primal(value(x))

"Nesting depth of a dual-number type: 0 for a plain real, 1 + depth of the value type."
depth(::Type{<:Real}) = 0
depth(::Type{Dual{T,V,P}}) where {T,V,P} = 1 + depth(V)

"""
    lift(k, jet, x)

Return the derivative of order `k` of a function of time at the epoch `x`, where
`jet[m+1]` is the derivative of order `m` at `primal(x)`. Requires
`k + depth(typeof(x)) < length(jet)`. For a dual `x` the value is order `k` at
`value(x)` and each partial is order `k+1` at `value(x)` times the partial of `x`.
"""
@inline lift(k::Int, jet::Tuple, ::Real) = jet[k+1]
@inline function lift(k::Int, jet::Tuple, x::Dual{T}) where {T}
    v = lift(k, jet, value(x))
    d = lift(k + 1, jet, value(x))
    p = partials(x)
    return map((a, b) -> Dual{T}(a, b * p), v, d)
end

"Orders 0 to N at the dual epoch `x` from the jet of order `N + depth` at `primal(x)`."
@inline lift_jet(jet::Tuple, x::Dual, ::Val{N}) where {N} =
    ntuple(i -> lift(i - 1, jet, x), Val(N + 1))

function NeoSpice.spk_jet(ks::KernelSet, targ::Int32, et::D, ref::Int32, obs::Int32,
        v::Val{N}) where {D<:Dual,N}
    jet = spk_jet(ks, targ, primal(et), ref, obs, Val(N + depth(D)))
    return lift_jet(jet, et, v)
end

function NeoSpice.frame_rotation(ks::KernelSet, from::Int32, to::Int32, et::D) where {D<:Dual}
    jet = rotation_jet(ks, from, to, primal(et), Val(depth(D)))
    return lift(0, jet, et)
end

NeoSpice.frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Dual, ::Val{0}) =
    (frame_rotation(ks, from, to, et),)
function NeoSpice.frame_jet(ks::KernelSet, from::Int32, to::Int32, et::D,
        v::Val{N}) where {D<:Dual,N}
    jet = rotation_jet(ks, from, to, primal(et), Val(N + depth(D)))
    return lift_jet(jet, et, v)
end

"`_ad_order(Dual)`: a dual number carries one derivative order more than its value type."
NeoSpice._ad_order(::Type{ForwardDiff.Dual{Tg,V,N}}) where {Tg,V,N} = 1 + NeoSpice._ad_order(V)

"`_with_primal(x::Dual, p)`: `x` with the plain value `p` and the same partials."
NeoSpice._with_primal(x::ForwardDiff.Dual{Tg}, p::Float64) where {Tg} =
    ForwardDiff.Dual{Tg}(NeoSpice._with_primal(ForwardDiff.value(x), p), ForwardDiff.partials(x))

# --- aberration corrections ---------------------------------------------------------
#
# A corrected position with a dual-number epoch is assembled in the same way, from the
# jet of `spkder` with the correction at the plain epoch. The value is the position of
# the plain call, and the partials are the exact derivatives of the position function
# (see the docstring of `spkder`). For `spkezr` the value of the velocity is the
# velocity that CSPICE defines and its partials are the next orders of the same jet;
# for the `LT` forms that velocity is not the derivative of the position, so the
# partial of the position and the value of the velocity differ there, as they do
# between `spkder` and `spkezr`. The light time is `|r| / c`, because the stellar
# aberration and the frame change are rotations.

"Light time `|r| / c` of a dual position; zero for the zero vector (a body from itself)."
@inline function _dual_lt(r::NeoSpice.Vec3{D}) where {D<:Dual}
    iszero(primal(r[1])) && iszero(primal(r[2])) && iszero(primal(r[3])) && return zero(D)
    return sqrt(r[1] * r[1] + r[2] * r[2] + r[3] * r[3]) / NeoSpice.CLIGHT
end

function NeoSpice.abcorr_jet(ks::KernelSet, targ::Int32, et::D, ref::Int32,
        corr::NeoSpice.Aberration, obs::Int32, v::Val{N}) where {D<:Dual,N}
    jet = NeoSpice.abcorr_jet(ks, targ, primal(et), ref, corr, obs, Val(N + depth(D)))
    return lift_jet(jet, et, v)
end

function NeoSpice.abcorr_position(ks::KernelSet, targ::Int32, et::D, ref::Int32,
        corr::NeoSpice.Aberration, obs::Int32) where {D<:Dual}
    t = primal(et)
    jet = NeoSpice.abcorr_jet(ks, targ, t, ref, corr, obs, Val(depth(D)))
    r0, _ = NeoSpice.abcorr_position(ks, targ, t, ref, corr, obs)
    r = lift(0, Base.setindex(jet, r0, 1), et)
    return r, _dual_lt(r)
end

function NeoSpice.abcorr_state(ks::KernelSet, targ::Int32, et::D, ref::Int32,
        corr::NeoSpice.Aberration, obs::Int32) where {D<:Dual}
    t = primal(et)
    jet = NeoSpice.abcorr_jet(ks, targ, t, ref, corr, obs, Val(depth(D) + 1))
    s0, _ = NeoSpice.abcorr_state(ks, targ, t, ref, corr, obs)
    r = lift(0, Base.setindex(jet, s0[StaticArrays.SVector(1, 2, 3)], 1), et)
    v = lift(1, Base.setindex(jet, s0[StaticArrays.SVector(4, 5, 6)], 2), et)
    return vcat(r, v), _dual_lt(r)
end

end # module
