# Kernel paths for the tests and the benchmarks, in one place.
#
# The kernels are NAIF generic kernels. `fetch_kernels` downloads
# the missing ones into `test/data/` (not under version control) and checks their
# SHA-256. Environment variables:
#   ALLSPICE_KERNEL_DIR        directory of the text kernels and the binary PCK
#   ALLSPICE_DE440S            path of de440s.bsp
#   ALLSPICE_DE430             path of de430.bsp
#   ALLSPICE_DOWNLOAD_KERNELS  0: do not download
#   ALLSPICE_REQUIRE_KERNELS   1: a group whose kernels are missing fails instead of
#                              being skipped
#
# A test group names the kernels that it needs; `missing_kernels` gives the ones that
# are not on disk, and the group is skipped with a visible message if there is one.

module Kernels

import Downloads, SHA

const DIR = get(ENV, "ALLSPICE_KERNEL_DIR", joinpath(@__DIR__, "data"))

const DE440S = get(ENV, "ALLSPICE_DE440S", joinpath(DIR, "de440s.bsp"))
const DE430 = get(ENV, "ALLSPICE_DE430", joinpath(DIR, "de430.bsp"))
const LSK = joinpath(DIR, "naif0012.tls")
const PCK11 = joinpath(DIR, "pck00011.tpc")
const PCK10 = joinpath(DIR, "pck00010.tpc")
const MOON_BPC = joinpath(DIR, "moon_pa_de440_200625.bpc")
const MOON_TF = joinpath(DIR, "moon_de440_200625.tf")
const MOON_TF_OLD = joinpath(DIR, "moon_080317.tf")
"Optional kernel with `BODY<id>_GM`; the GM checks are skipped if it is absent."
const GM_TPC = joinpath(DIR, "gm_de431.tpc")

"True if a group whose kernels are missing must fail (ALLSPICE_REQUIRE_KERNELS=1)."
const REQUIRED = get(ENV, "ALLSPICE_REQUIRE_KERNELS", "0") == "1"

const NAIF = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/"

"""
URLs and SHA-256 of each kernel, by file name. The URLs of a kernel hold the same file
and are tried in order.
"""
const SOURCES = Dict(
    "naif0012.tls" => ([NAIF * "lsk/naif0012.tls"],
        "678e32bdb5a744117a467cd9601cd6b373f0e9bc9bbde1371d5eee39600a039b"),
    "pck00011.tpc" => ([NAIF * "pck/pck00011.tpc"],
        "3dff7b1dbeceaa01f25467767d3fa25816051c85d162d1edf04acb310ee28bb1"),
    "pck00010.tpc" => ([NAIF * "pck/pck00010.tpc"],
        "59468328349aa730d18bf1f8d7e86efe6e40b75dfb921908f99321b3a7a701d2"),
    "moon_pa_de440_200625.bpc" => ([NAIF * "pck/moon_pa_de440_200625.bpc"],
        "60cd55aa401ea2ea97360636f567554bfe4e37bb829f901b4460a455dfaf783f"),
    "moon_de440_200625.tf" => ([NAIF * "fk/satellites/a_old_versions/moon_de440_200625.tf"],
        "a1352d73783c87ba7c452441ade92348ac7ac97dfd00e34daa4bbcae46086e8f"),
    "moon_080317.tf" => ([NAIF * "fk/satellites/moon_080317.tf"],
        "78732477b96f9863e7b0d65bcee3c22b8707ca5ed0db56d1173319cb2e8c7993"),
    "gm_de431.tpc" => ([NAIF * "pck/gm_de431.tpc"],
        "15756c162151853f329d473fa00c38d35e62d83c3e25ddceece6516bbd98738b"),
    "de440s.bsp" => ([NAIF * "spk/planets/de440s.bsp"],
        "c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2"),
    "de430.bsp" => ([NAIF * "spk/planets/de430.bsp"],
        "6e1b277c5f07135a84950604b83e56b736be696a7f3560bcddb1d4aeb944fca1"),
)

"Seconds of the waits before the second and the third attempt of a download."
const RETRY_WAITS = (10.0, 30.0)

"""
    download_any(urls, dest, io) -> Bool

Download the first of `urls` that answers into `dest`, with three attempts for each URL
and the waits `RETRY_WAITS` between them, and report each failure on `io`. Return
`false` if every attempt failed.
"""
function download_any(urls::Vector{String}, dest::AbstractString, io::IO)
    for url in urls, attempt in 1:length(RETRY_WAITS)+1
        # invariant: every attempt before this one failed
        attempt > 1 && sleep(RETRY_WAITS[attempt-1])
        try
            println(io, "Downloading ", url, attempt > 1 ? " (attempt $attempt)" : "")
            Downloads.download(url, dest)
            return true
        catch err
            println(io, "Download of ", url, " failed: ", sprint(showerror, err))
        end
    end
    return false
end

"Paths of all kernels of the tests on this machine."
all_kernels() = String[LSK, PCK11, PCK10, MOON_BPC, MOON_TF, MOON_TF_OLD, GM_TPC, DE440S, DE430]

"""
    fetch_kernels(io = stdout) -> Vector{String}

Download each kernel of [`all_kernels`](@ref) that is not on disk and whose file name
is in `SOURCES`, unless ALLSPICE_DOWNLOAD_KERNELS=0, and return the paths that are
still missing. Each URL of a kernel gets three attempts (see `download_any`). A download goes to a temporary file in the target directory and is
moved into place after its SHA-256 is checked, so that processes that fetch at the same
time never see a partial file. A failed download is reported on `io` and leaves the
kernel missing; a file with the wrong SHA-256 is deleted and throws an error.
"""
function fetch_kernels(io::IO=stdout)
    get(ENV, "ALLSPICE_DOWNLOAD_KERNELS", "1") == "0" && return missing_kernels(all_kernels()...)
    for path in all_kernels()     # invariant: the kernels before `path` were fetched or reported
        (isfile(path) || !haskey(SOURCES, basename(path))) && continue
        urls, sha = SOURCES[basename(path)]
        mkpath(dirname(path))
        tmp, tio = mktemp(dirname(path); cleanup=false)
        close(tio)
        if !download_any(urls, tmp, io)
            rm(tmp; force=true)
            continue
        end
        got = bytes2hex(open(SHA.sha256, tmp))
        if got != sha
            rm(tmp; force=true)
            error("SHA-256 of $(basename(path)) is $got, expected $sha")
        end
        mv(tmp, path; force=true)
    end
    gone = missing_kernels(all_kernels()...)
    @assert all(p -> !isfile(p), gone)
    return gone
end

"Paths of `paths` that are not regular files."
missing_kernels(paths...) = String[p for p in paths if !isfile(p)]

const ALIAS_PATH = Ref("")

"""
    moon_alias_kernel() -> String

Write (once per process) a frame kernel that defines `MOON_PA` and `MOON_ME` as
fixed-offset aliases of `MOON_PA_DE440` and `MOON_ME_DE440_ME421`, and return its path.

`moon_080317.tf` defines these two names as aliases of the DE421 frames, and the binary
PCK for DE421 is not among the test kernels, so CSPICE itself cannot evaluate them from the
kernels at hand. The alias kernel keeps the generic names, uses the DE440 orientation
of the test kernels, and exercises a chain of two fixed-offset frames (`MOON_ME`).
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
