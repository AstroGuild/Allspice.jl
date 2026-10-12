# Kernel paths for the tests and the benchmarks, in one place.
#
# Environment overrides:
#   NEOSPICE_KERNEL_DIR  directory of the text kernels and the binary PCK
#   NEOSPICE_DE440S      path of de440s.bsp
#   NEOSPICE_DE405       path of de405.bsp
#
# A test group names the kernels that it needs; `missing_kernels` gives the ones that
# are not on disk, and the group is skipped with a visible message if there is one.

module Kernels

const DIR = get(ENV, "NEOSPICE_KERNEL_DIR", "/home/alex/.tudat/resource/spice_kernels")

const DE440S = get(ENV, "NEOSPICE_DE440S", "/home/alex/data/ephemerides/de440s.bsp")
const DE405 = get(ENV, "NEOSPICE_DE405", "/home/alex/malto/data/de405.bsp")
const LSK = joinpath(DIR, "naif0012.tls")
const PCK11 = joinpath(DIR, "pck00011.tpc")
const PCK10 = joinpath(DIR, "pck00010.tpc")
const MOON_BPC = joinpath(DIR, "moon_pa_de440_200625.bpc")
const MOON_TF = joinpath(DIR, "moon_de440_200625.tf")
const MOON_TF_OLD = joinpath(DIR, "moon_080317.tf")
"Optional kernel with `BODY<id>_GM`; the GM checks are skipped if it is absent."
const GM_TPC = joinpath(DIR, "gm_de431.tpc")

"Paths of `paths` that are not regular files."
missing_kernels(paths...) = String[p for p in paths if !isfile(p)]

const ALIAS_PATH = Ref("")

"""
    moon_alias_kernel() -> String

Write (once per process) a frame kernel that defines `MOON_PA` and `MOON_ME` as
fixed-offset aliases of `MOON_PA_DE440` and `MOON_ME_DE440_ME421`, and return its path.

`moon_080317.tf` defines these two names as aliases of the DE421 frames, and the binary
PCK for DE421 is not on this machine, so CSPICE itself cannot evaluate them from the
kernels at hand. The alias kernel keeps the generic names, uses the DE440 orientation
that is on disk, and exercises a chain of two fixed-offset frames (`MOON_ME`).
"""
function moon_alias_kernel()
    isfile(ALIAS_PATH[]) && return ALIAS_PATH[]
    path = joinpath(mktempdir(), "moon_alias_de440.tf")
    write(path, """
    KPL/FK

    Aliases MOON_PA and MOON_ME for the DE440 lunar frames (written by the tests).

    \\begindata

       FRAME_MOON_PA                 = 31000
       FRAME_31000_NAME              = 'MOON_PA'
       FRAME_31000_CLASS             = 4
       FRAME_31000_CLASS_ID          = 31000
       FRAME_31000_CENTER            = 301

       TKFRAME_31000_SPEC            = 'MATRIX'
       TKFRAME_31000_RELATIVE        = 'MOON_PA_DE440'
       TKFRAME_31000_MATRIX          = ( 1 0 0
                                         0 1 0
                                         0 0 1 )

       FRAME_MOON_ME                 = 31001
       FRAME_31001_NAME              = 'MOON_ME'
       FRAME_31001_CLASS             = 4
       FRAME_31001_CLASS_ID          = 31001
       FRAME_31001_CENTER            = 301

       TKFRAME_31001_SPEC            = 'MATRIX'
       TKFRAME_31001_RELATIVE        = 'MOON_ME_DE440_ME421'
       TKFRAME_31001_MATRIX          = ( 1 0 0
                                         0 1 0
                                         0 0 1 )

    \\begintext
    """)
    ALIAS_PATH[] = path
    @assert isfile(path)
    return path
end

"Text kernels of the standard set, in load order (the alias kernel is written on demand)."
text_set() = String[LSK, PCK11, MOON_TF, moon_alias_kernel()]

"Standard set: text kernels, the lunar binary PCK and `de440s.bsp`, in load order."
full_set() = String[text_set()..., MOON_BPC, DE440S]

"Kernels that must exist on disk for `full_set()`."
const FULL_SET_FILES = (LSK, PCK11, MOON_TF, MOON_BPC, DE440S)
const TEXT_SET_FILES = (LSK, PCK11, MOON_TF)

end # module
