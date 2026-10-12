module NeoSpice

using StaticArrays
using Mmap
using LinearAlgebra

export SpiceContext, SpiceError, default_context
export furnsh, unload, kclear, ktotal
export spkezr, spkpos, spkgeo, spkder
export Aberration, spkssb, spkltc, spkaps, spkacs, spkapo, ltime, stelab, stlabx
export pxform, sxform, pxder
export gdpool, gcpool, gipool, bodvrd, bodn2c, bodc2n, bods2c, namfrm, frmnam
export deltet, utc2et, str2et, et2utc

include("types.jl")
include("taylor.jl")        # Taylor{N,T} truncated series, `primal`
include("daf.jl")           # open_daf, segment descriptors
include("spk_eval.jl")      # spk_eval, pck_eval, find_segment
include("textkernel.jl")    # parse_text_kernel, kernel_arch, metakernel_files
include("names.jl")         # body names and codes, built-in and pool-defined
include("frames.jl")        # frame definitions, built-in and pool-defined
include("time.jl")          # leap seconds, deltet, utc2et, et2utc
include("rotations.jl")     # build_tk, build_bodyrot, frame_rotation
include("context.jl")       # build_kernelset, SpiceContext, furnsh, unload, kclear
include("spk_api.jl")       # spkgeo, spkezr, spkpos, spkder
include("aberration.jl")    # Aberration, light time, stellar aberration, spkltc, ...
include("frame_api.jl")     # pxform, sxform, pxder
include("pool_api.jl")      # gdpool, gcpool, gipool, bodvrd, bodn2c, ...
include("precompile.jl")    # precompile directives for the common signatures

end # module
