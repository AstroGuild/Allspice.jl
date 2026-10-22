"""
    AllspiceUtilities

Julia versions of the NAIF toolkit programs `brief`, `commnt`, `spkmerge`, `spkdiff` and
`frmdiff`, built on Allspice.jl. Each program is a set of functions that return
structured results; a printer gives a text report in the layout of the NAIF program.
"""
module AllspiceUtilities

using Allspice: Allspice, SpiceContext, Window, furnsh, getfat, wncard, wnfetd,
    wnintd, wnunid, wndifd, spkobj, spkcov, pckfrm, pckcov, ckcov, bodc2n, bodn2c, namfrm,
    frinfo, str2et, etcal, spkezr, pxform, sxform, xf2rav, raxisa, dafopr, dafopw, dafcls,
    dafbbs, daffpa, dafgs, dafgn, dafec, dafac, dafdc, dasopr, dasopw, dascls, dasec,
    dasac, dasdc, spkuds, spkopn, spksub, spkcls
using LinearAlgebra: norm, cross, dot
using Printf: @sprintf
using StaticArrays: SVector, SMatrix

export Coverage, brief, printbrief
export comments, addcomments, deletecomments
export MergeSource, spkmerge
export SpkDiff, spkdiff, basicstats, viewstats
export FrmDiff, frmdiff

include("common.jl")
include("brief.jl")
include("commnt.jl")
include("spkmerge.jl")
include("spkdiff.jl")
include("frmdiff.jl")

end
