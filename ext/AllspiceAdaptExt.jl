module AllspiceAdaptExt

using Allspice: EphemerisTable
import Adapt

# Move the coefficient array of an ephemeris table to another array type and keep the
# scalar fields. With this method `adapt(CuArray, tab)` gives a table on the GPU, and
# CUDA.jl converts that table to an `isbits` value when it is an argument of a kernel
# or is captured by a broadcast. `cu(tab)` also converts the coefficients to `Float32`.
Adapt.adapt_structure(to, tab::EphemerisTable) = EphemerisTable(Adapt.adapt(to, tab.coef),
    tab.t0, tab.dt, tab.nint, tab.ncoef, tab.target, tab.observer, tab.frame, tab.maxerr)

end # module
