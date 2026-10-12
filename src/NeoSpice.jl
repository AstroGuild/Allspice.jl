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
export deltet, utc2et, str2et, et2utc, etcal, timout, tpictr, tparse, timdef, unitim
export jul2gr, gr2jul
export j1900, j1950, j2000, j2100, jyear, tyear, spd, b1900, b1950, clight
export halfpi, twopi, rpd, dpr
export sct2e, sce2c, sce2t, sce2s, scs2e, scencd, scdecd, sctiks, scfmt, scpart, sctype
export ckgp, ckgpav, ckcov, ckobj, ckmeta
export DAFHandle, dafopr, dafopw, dafonw, dafcls, dafhsf, dafrfr, getfat
export dafbfs, dafbbs, daffna, daffpa, dafgs, dafgn, dafps, dafus, dafgda, dafgsr
export dafbna, dafada, dafena, dafws, dafrs, dafrn, dafac, dafec, dafdc
export spkopn, spkopa, spkcls, spkpds, spkw01, spkw02, spkw03, spkw05, spkw08, spkw09
export spkw12, spkw13, spkw15, spkw17, spkw18, spkw19, spkw20, spkw21
export ckopn, ckcls, ckw01, ckw02, ckw03, ckw05, ckw06
export spkobj, spkobj!, spkcov, spkcov!, kinfo, kdata, tkvrsn
export Window, card, wncard, wninsd!, wnvald!, wnunid, wnintd, wndifd, wncomd
export wnelmd, wnincd, wnfetd, wnsumd
export vadd, vsub, vscl, vminus, vdot, vcrss, vnorm, vhat, vdist, vsep, vrel, vproj, vperp
export vlcom, vlcom3, vrotv, vzero, vequ, vpack, vupack, ucrss, unorm, vtmv
export mxv, mtxv, mxm, mtxm, mxmt, xpose, xpose6, invert, invort, trace, isrot, ident
export rotate, rotmat, rotvec, twovec, twovxf, dvhat, dvdot, dvcrss, ducrss, dvnorm, dvsep
export vaddg, vsubg, vsclg, vminug, vequg, vdotg, vnormg, vhatg, unormg, vdistg, vsepg, vrelg
export vprojg, vlcomg, vzerog, vtmvg, mxvg, mtxvg, mxmg, mtxmg, mxmtg, xposeg
export eul2m, m2eul, eul2xf, xf2eul, axisar, raxisa, m2q, q2m, qxq, qdq2av, rav2xf, xf2rav
export sharpr

include("types.jl")
include("taylor.jl")        # Taylor{N,T} truncated series, `primal`
include("daf.jl")           # open_daf, segment descriptors
include("spk_eval.jl")      # spk_eval, pck_eval, find_segment
include("spk_more.jl")      # SPK types 1, 5, 8, 9, 12 to 15, 17 to 21
include("textkernel.jl")    # parse_text_kernel, kernel_arch, metakernel_files
include("names.jl")         # body names and codes, built-in and pool-defined
include("frames.jl")        # frame definitions, built-in and pool-defined
include("time.jl")          # leap seconds, deltet, unitim, ttrans, calendars, constants
include("timeparse.jl")     # tpartv, str2et, utc2et, tparse, tpictr, timdef rules
include("timeout.jl")       # number formats, timout, et2utc, etcal
include("rotations.jl")     # build_tk, build_bodyrot, frame_rotation
include("sclk.jl")          # type 1 spacecraft clocks: build_sclk, sct2e, sce2c, ...
include("ck.jl")            # CK segment evaluation, types 1 to 6; class 3 frames
include("context.jl")       # build_kernelset, SpiceContext, furnsh, unload, kclear
include("spk_api.jl")       # spkgeo, spkezr, spkpos, spkder
include("aberration.jl")    # Aberration, light time, stellar aberration, spkltc, ...
include("frame_api.jl")     # pxform, sxform, pxder
include("pool_api.jl")      # gdpool, gcpool, gipool, bodvrd, bodn2c, ...
include("time_api.jl")      # context forms of the time calls, timdef
include("ck_api.jl")        # ckgp, ckgpav, ckcov, ckobj, ckmeta; CSPICE state rotation
include("windows.jl")       # Window, wninsd!, wnunid, wnintd, wndifd, wncomd, ...
include("daf_access.jl")    # DAFHandle, dafopr, dafbfs, daffna, dafgs, dafgda, dafec, getfat
include("daf_write.jl")     # dafonw, dafopw, dafbna, dafada, dafena, dafac, dafdc
include("spk_write.jl")     # spkopn, spkopa, spkcls, spkw01 ... spkw21
include("ck_write.jl")      # ckopn, ckcls, ckw01, ckw02, ckw03, ckw05, ckw06
include("kernel_info.jl")   # spkobj, spkcov, ktotal(kind), kdata, kinfo, tkvrsn
include("math/vectors.jl")      # vadd ... vsep, mxv, mxm, invert, rotate, twovec, dvhat, g forms
include("math/rotations.jl")    # eul2m, m2eul, m2q, q2m, qxq, axisar, raxisa, rav2xf, xf2eul
include("precompile.jl")    # precompile directives for the common signatures

end # module
