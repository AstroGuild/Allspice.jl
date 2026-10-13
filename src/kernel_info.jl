# Information calls: the bodies and the coverage of an SPK file (`spkobj`, `spkcov`),
# the load list of a context (`ktotal`, `kdata`, `kinfo`) and the toolkit version.

"CSPICE version whose behaviour and file formats this package follows."
const CSPICE_VERSION = "CSPICE_N0067"

"""
    tkvrsn(item="TOOLKIT") -> String

Version string of the toolkit, as `tkvrsn_c`. For the item `"TOOLKIT"` (any case, or
the symbol `:TOOLKIT` as in SPICE.jl) the result is
`"Allspice_v<package version>/CSPICE_N0067"`: the version of this package and, after
the slash, the CSPICE version that the port follows (`tkvrsn_c` gives
`"CSPICE_N0067"`). Any other item gives `"No version found."`, as in CSPICE.
"""
function tkvrsn(item::Union{AbstractString,Symbol}="TOOLKIT")
    uppercase(strip(String(item))) == "TOOLKIT" || return "No version found."
    return string("Allspice_v", pkgversion(@__MODULE__), "/", CSPICE_VERSION)
end

# ---------------------------------------------------------------------------------------
# Load list
# ---------------------------------------------------------------------------------------

"Kernel types that `ktotal`, `kdata` and `kinfo` report, as in `keeper.f`."
const KERNEL_KINDS = ("SPK", "CK", "PCK", "DSK", "EK", "TEXT", "META")

"""
    kernel_filetype(entry::Loaded) -> String

Type of a load-list entry as CSPICE reports it: `"SPK"`, `"CK"`, `"PCK"`, `"DSK"` for a
binary kernel, `"META"` for a text kernel that assigns `KERNELS_TO_LOAD`, `"TEXT"` for every
other text kernel.
"""
function kernel_filetype(entry::Loaded)
    content = entry.content
    content isa TextKernel || return String(content.kind)
    return any(a -> a.name == "KERNELS_TO_LOAD", content.assignments) ? "META" : "TEXT"
end

"""
    _kind_filter(kind) -> Union{Nothing,Vector{String}}

The kernel types that the string `kind` selects: `nothing` for `"ALL"` (every type),
else the blank-separated words of `kind` in upper case. A word that is not a kernel
type selects nothing, as in CSPICE.
"""
function _kind_filter(kind::AbstractString)
    words = String[uppercase(w) for w in split(kind)]
    words == ["ALL"] && return nothing
    return words
end

"Load-list entries of `ks` whose type is selected by `kind`, in load order."
function _kernels_of_kind(ks::KernelSet, kind::AbstractString)
    words = _kind_filter(kind)
    words === nothing && return ks.loaded
    return Loaded[e for e in ks.loaded if kernel_filetype(e) in words]
end

"""
    ktotal(kind="ALL") -> Int

Number of loaded kernels of the types in `kind`, as `ktotal_c`: `kind` is `"ALL"` or a
blank-separated list of `SPK`, `CK`, `PCK`, `DSK`, `EK`, `TEXT`, `META` (any case). A
meta-kernel counts as one `META` entry plus one entry for each file that it lists.
"""
ktotal(ks::KernelSet, kind::AbstractString="ALL") = length(_kernels_of_kind(ks, kind))
ktotal(ctx::SpiceContext, kind::AbstractString) = ktotal(kernelset(ctx), kind)
ktotal(kind::AbstractString) = ktotal(default_context(), kind)

"""
    kdata(which, kind="ALL") -> Union{Nothing,Tuple{String,String,String,Int}}

Data of the loaded kernel number `which` (from 1, as in SPICE.jl) among the kernels
of the types in `kind`, in load order, as `kdata_c`: `(file, filtyp, source, handle)`,
where `source` is the meta-kernel that loaded the file or `""`. The handle is always
0: a loaded kernel has no DAF handle in this package. Returns `nothing` if there is
no such kernel.
"""
function kdata(ks::KernelSet, which::Integer, kind::AbstractString="ALL")
    list = _kernels_of_kind(ks, kind)
    (1 <= which <= length(list)) || return nothing
    e = list[which]
    return e.path, kernel_filetype(e), e.source, 0
end
kdata(ctx::SpiceContext, which::Integer, kind::AbstractString="ALL") =
    kdata(kernelset(ctx), which, kind)
kdata(which::Integer, kind::AbstractString="ALL") = kdata(default_context(), which, kind)

"""
    kinfo(file) -> Union{Nothing,Tuple{String,String,Int}}

Data of the loaded kernel with the path `file` (compared as a string, as in CSPICE),
as `kinfo_c`: `(filtyp, source, handle)` of the first load-list entry with that path,
or `nothing` if the file is not loaded. The handle is always 0 (see [`kdata`](@ref)).
"""
function kinfo(ks::KernelSet, file::AbstractString)
    i = findfirst(e -> e.path == file, ks.loaded)
    i === nothing && return nothing
    e = ks.loaded[i]
    return kernel_filetype(e), e.source, 0
end
kinfo(ctx::SpiceContext, file::AbstractString) = kinfo(kernelset(ctx), file)
kinfo(file::AbstractString) = kinfo(default_context(), file)

# ---------------------------------------------------------------------------------------
# Bodies and coverage of an SPK file
# ---------------------------------------------------------------------------------------

"""
    _require_binary_kernel(who, path, kind) -> Nothing

Check with `getfat` that `path` is a binary kernel of type `kind` (`"SPK"`, `"CK"`,
`"PCK"`), as `spkobj`, `spkcov`, `ckobj`, `ckcov` and `pckcov` do. Throws
`SpiceError(:invalidformat)` for a transfer file, `SpiceError(:invalidarchtype)` for a
file that is not a DAF and `SpiceError(:invalidfiletype)` for a DAF of another type.
"""
function _require_binary_kernel(who::AbstractString, path::AbstractString, kind::AbstractString)
    arch, typ = getfat(path)
    arch == "XFR" && throw(SpiceError(:invalidformat, "$who: input file $path has " *
        "architecture $arch; the file must be a binary $kind file"))
    arch == "DAF" || throw(SpiceError(:invalidarchtype, "$who: input file $path has " *
        "architecture $arch; the file must be a binary $kind file (DAF architecture)"))
    typ == kind || throw(SpiceError(:invalidfiletype, "$who: input file $path has file " *
        "type $typ; the file must be a binary $kind file"))
    return nothing
end

"""
    _each_summary(f, path) -> Nothing

Call `f(dc, ic, handle)` for the unpacked summary of every array of the DAF `path`,
in file order, on a private read handle that is closed afterwards. This is the loop
of `spkobj` and `spkcov`; a coverage routine for another kernel type (`ckcov`,
`pckcov`) has the same structure: `_require_binary_kernel`, then this loop with a
body that inserts the intervals of the wanted object into a `Window`.
"""
function _each_summary(f, path::AbstractString)
    h = dafopr(path)
    try
        dafbfs(h)
        while daffna(h)                  # invariant: f was called for every array before this one
            dc, ic = dafus(dafgs(h), h.nd, h.ni)
            f(dc, ic, h)
        end
    finally
        dafcls(h)
    end
    return nothing
end

"""
    spkobj!(ids, path) -> ids

Insert the ID codes of the bodies that have data in the SPK file `path` into the
sorted vector `ids`, which stays sorted and without duplicates (the set semantics of
the SPICE cell of `spkobj_c`; argument order of SPICE.jl). The file need not be
loaded. Throws `SpiceError` as [`spkobj`](@ref).
"""
function spkobj!(ids::AbstractVector{<:Integer}, path::AbstractString)
    _require_binary_kernel("spkobj", path, "SPK")
    _each_summary(path) do dc, ic, h
        k = searchsortedfirst(ids, ic[1])
        (k <= length(ids) && ids[k] == ic[1]) || insert!(ids, k, ic[1])
    end
    return ids
end

"""
    spkobj(path) -> Vector{Int}

ID codes of the bodies that have data in the SPK file `path`, sorted and without
duplicates, as `spkobj_c`. The file is read from disk and need not be loaded.
Preconditions (each throws `SpiceError`): the file exists (`:filenotfound`) and is a
binary SPK (`:invalidformat`, `:invalidarchtype`, `:invalidfiletype`).
"""
spkobj(path::AbstractString) = spkobj!(Int[], path)

"""
    spkcov!(cover, path, idcode) -> cover

Insert the coverage of the body `idcode` in the SPK file `path` into the window
`cover` (union with its old content), as `spkcov_c` (argument order of SPICE.jl).
"""
function spkcov!(cover::Window, path::AbstractString, idcode::Integer)
    _require_binary_kernel("spkcov", path, "SPK")
    _each_summary(path) do dc, ic, h
        ic[1] == idcode && wninsd!(cover, dc[1], dc[2])
    end
    return cover
end

"""
    spkcov(path, idcode) -> Window

Coverage window of the body `idcode` in the SPK file `path`, as `spkcov_c`: the union
of the descriptor intervals [start, stop] (seconds past J2000 TDB) of the segments of
that body, as sorted disjoint intervals. `collect(w)` gives the endpoints as a
`Vector{Float64}`. The file is read from disk and need not be loaded. The window is
empty if the body has no data. Throws `SpiceError` as [`spkobj`](@ref).
"""
spkcov(path::AbstractString, idcode::Integer) = spkcov!(Window(), path, idcode)
