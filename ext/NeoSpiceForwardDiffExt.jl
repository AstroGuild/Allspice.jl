module NeoSpiceForwardDiffExt

using NeoSpice
using ForwardDiff

"`primal(x::Dual)`: the plain `Float64` value of a dual number, with all partials removed."
NeoSpice.primal(x::ForwardDiff.Dual) = NeoSpice.primal(ForwardDiff.value(x))

end # module
