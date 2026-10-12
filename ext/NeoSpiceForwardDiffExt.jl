module NeoSpiceForwardDiffExt

using NeoSpice
using ForwardDiff

"`primal(x::Dual)`: the plain `Float64` value of a dual number, with all partials removed."
NeoSpice.primal(x::ForwardDiff.Dual) = NeoSpice.primal(ForwardDiff.value(x))

"`_ad_order(Dual)`: a dual number carries one derivative order more than its value type."
NeoSpice._ad_order(::Type{ForwardDiff.Dual{Tg,V,N}}) where {Tg,V,N} = 1 + NeoSpice._ad_order(V)

"`_with_primal(x::Dual, p)`: `x` with the plain value `p` and the same partials."
NeoSpice._with_primal(x::ForwardDiff.Dual{Tg}, p::Float64) where {Tg} =
    ForwardDiff.Dual{Tg}(NeoSpice._with_primal(ForwardDiff.value(x), p), ForwardDiff.partials(x))

end # module
