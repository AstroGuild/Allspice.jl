module NeoSpice

using StaticArrays
using Mmap
using LinearAlgebra

export SpiceContext, SpiceError, default_context
export furnsh, unload, kclear, ktotal
export spkezr, spkpos, spkgeo, spkder
export pxform, sxform, pxder
export gdpool, gcpool, gipool, bodvrd, bodn2c, bodc2n, bods2c, namfrm, frmnam
export deltet, utc2et, str2et, et2utc
export DAFHandle, dafopr, dafopw, dafonw, dafcls, dafhsf, dafrfr, getfat
export dafbfs, dafbbs, daffna, daffpa, dafgs, dafgn, dafps, dafus, dafgda, dafgsr
export dafbna, dafada, dafena, dafws, dafrs, dafrn, dafac, dafec, dafdc
export spkopn, spkopa, spkcls, spkpds, spkw01, spkw02, spkw03, spkw05, spkw08, spkw09
export spkw12, spkw13, spkw15, spkw17, spkw18, spkw19, spkw20, spkw21
export ckopn, ckcls, ckw01, ckw02, ckw03, ckw05, ckw06
export spkobj, spkobj!, spkcov, spkcov!, kinfo, kdata, tkvrsn
export Window, card, wncard, wninsd!, wnvald!, wnunid, wnintd, wndifd, wncomd
export wnelmd, wnincd, wnfetd, wnsumd

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
include("frame_api.jl")     # pxform, sxform, pxder
include("pool_api.jl")      # gdpool, gcpool, gipool, bodvrd, bodn2c, ...
include("windows.jl")       # Window, wninsd!, wnunid, wnintd, wndifd, wncomd, ...
include("daf_access.jl")    # DAFHandle, dafopr, dafbfs, daffna, dafgs, dafgda, dafec, getfat
include("daf_write.jl")     # dafonw, dafopw, dafbna, dafada, dafena, dafac, dafdc
include("spk_write.jl")     # spkopn, spkopa, spkcls, spkw01 ... spkw21
include("ck_write.jl")      # ckopn, ckcls, ckw01, ckw02, ckw03, ckw05, ckw06
include("kernel_info.jl")   # spkobj, spkcov, ktotal(kind), kdata, kinfo, tkvrsn
include("precompile.jl")    # precompile directives for the common signatures

end # module
