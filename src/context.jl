# Kernel management: the load list, the immutable `KernelSet` snapshot that is built
# from it, and the `SpiceContext` that holds the current snapshot.

"Pool variables that `furnsh` uses for a meta-kernel and then deletes, as CSPICE does."
const METAKERNEL_VARIABLES = ("KERNELS_TO_LOAD", "PATH_SYMBOLS", "PATH_VALUES")

empty_kernelset() = build_kernelset(Loaded[])

"""
    build_kernelset(loaded::Vector{Loaded}) -> KernelSet

Build the snapshot for the load list `loaded` (load order). The function does not read
any file: it uses the `DAF` and `TextKernel` objects of the entries. Segment `file`
numbers index `daf`; `spk`, `pck` and `ck` are in CSPICE priority order (later file
first, and inside one file later segment first); the pool holds all assignments in load
order. The spacecraft clocks and the CK meta data are built from the pool.

    build_kernelset(loaded, edited, bodydefs) -> KernelSet

The same with the kernel pool `edited` in the place of the pool of the text kernels
(a pool with direct edits, see `pool_write.jl`; `nothing` for no edits) and with the
name and code pairs `bodydefs` of `boddef`.
"""
function build_kernelset(loaded::Vector{Loaded},
        edited::Union{Nothing,Dict{String,PoolValue}}=nothing,
        bodydefs::Vector{Tuple{String,Int32}}=Tuple{String,Int32}[])
    daf = Vector{Float64}[]
    spk = Dict{Int32,Vector{Segment}}()
    pck = Dict{Int32,Vector{Segment}}()
    ck = Dict{Int32,Vector{CKSegment}}()
    pool = Dict{String,PoolValue}()
    dskfiles = DSKFile[]
    ekfiles = EKFile[]
    filenumber = zeros(Int, length(loaded))
    for (i, entry) in enumerate(loaded)
        content = entry.content
        if content isa DAF
            push!(daf, content.data)
            filenumber[i] = length(daf)
        elseif content isa DSKFile
            push!(dskfiles, content)
        elseif content isa EKFile
            any(f -> f.path == content.path, ekfiles) || push!(ekfiles, content)
        elseif edited === nothing
            apply_assignments!(pool, content)
        end
    end
    for name in METAKERNEL_VARIABLES
        delete!(pool, name)
    end
    edited === nothing || (pool = edited)
    for i in length(loaded):-1:1
        content = loaded[i].content
        content isa DAF || continue
        if content.kind === :CK
            for j in length(content.cksegments):-1:1
                seg = withfile(content.cksegments[j], filenumber[i])
                push!(get!(() -> CKSegment[], ck, seg.inst), seg)
            end
            continue
        end
        table = content.kind === :PCK ? pck : spk
        for j in length(content.segments):-1:1
            seg = withfile(content.segments[j], filenumber[i])
            push!(get!(() -> Segment[], table, seg.target), seg)
        end
    end
    npool = _pool_x_namepool(pool, edited !== nothing)
    name2code, code2name = build_names(npool)
    _pool_x_bodydefs!(name2code, code2name, bodydefs)
    frames, framename2id = build_frames(npool)
    tk = build_tk(pool, frames)
    bodyrot = build_bodyrot(pool)
    lsk = build_lsk(pool)
    return _dyn_levels(KernelSet(copy(loaded), daf, spk, pck, pool, name2code, code2name, frames,
        framename2id, tk, bodyrot, lsk, ck, build_sclk(pool), build_ckmeta(pool),
        Dict{Int32,DynFrame}(), Dict{Int32,SwitchFrame}(), nothing, false,
        Tuple{String,Int32}[], dskfiles, _dsk_segment_table(dskfiles), ekfiles,
        _ek_build_tables(ekfiles)),
        edited !== nothing, bodydefs)
end

SpiceContext() = SpiceContext(empty_kernelset(), ReentrantLock())

const DEFAULT_CONTEXT = SpiceContext()

"""
    default_context() -> SpiceContext

Return the process-wide context that every call without a `ctx` argument uses.

The result is always the same object.

# Examples
```jldoctest
julia> default_context() === default_context()
true
```

See also [`SpiceContext`](@ref).
"""
default_context() = DEFAULT_CONTEXT

"""
    kernelset(ctx::SpiceContext) -> KernelSet

Return the current snapshot of `ctx` with one atomic read; the result never changes.
"""
kernelset(ctx::SpiceContext) = @atomic ctx.ks

"""
    read_kernel(path, yearbound=1969) -> Union{DAF,TextKernel}

Open the binary kernel (DAF, or DSK) or parse the text kernel at `path` (with the year
bound `yearbound` of `tsetyr` for its `@` dates); throws `SpiceError` if the file does
not exist.
"""
function read_kernel(path::AbstractString, yearbound::Int=1969)
    isfile(path) || throw(SpiceError(:nofile, "kernel file not found: $path"))
    arch = kernel_arch(path)
    arch === :DSK && return open_dsk(path)
    arch === :EK && return open_ek(path)
    return (arch === :SPK || arch === :PCK || arch === :CK) ? open_daf(path) : parse_text_kernel(path, yearbound)
end

"""
    load_entries(path, yearbound=1969) -> Vector{Loaded}

Read the kernel at `path` and return its load-list entries: one entry for the kernel
and, for a meta-kernel, one entry for each listed file (in order, `source == path`).
"""
function load_entries(path::AbstractString, yearbound::Int=1969)
    p = String(path)
    content = read_kernel(p, yearbound)
    entries = Loaded[Loaded(p, "", content)]
    if content isa TextKernel
        for file in metakernel_files(content)
            sub = read_kernel(file, yearbound)
            if sub isa TextKernel && !isempty(metakernel_files(sub))
                throw(SpiceError(:unsupported,
                    "meta-kernel $p lists the meta-kernel $file; CSPICE does not allow this"))
            end
            push!(entries, Loaded(String(file), p, sub))
        end
    end
    return entries
end

"""
    furnsh([ctx,] path) -> Nothing

Load the kernel file at `path` into `ctx`, as `furnsh_c`.

The type is read from the file: SPK, binary PCK, CK, DSK, EK, text kernel or
meta-kernel. A meta-kernel loads the files of its `KERNELS_TO_LOAD` list in order,
with `PATH_SYMBOLS` and `PATH_VALUES` applied as in CSPICE.
Without `ctx` the kernel goes into [`default_context`](@ref).

- A load that fails leaves `ctx` unchanged.
- A later kernel has priority over an earlier one, as in CSPICE.
- A path that is already loaded gets a second entry, as in CSPICE.
- Direct edits of the kernel pool ([`pdpool`](@ref), [`lmpool`](@ref), ...) stay; the
  assignments of a text kernel are applied on top of them.
- An SPK, PCK or CK file in the non-native byte order is read; it is converted in memory.

Throws `SpiceError(:badarg)` if `path` is blank, `SpiceError(:nofile)` if the file (or
a file that a meta-kernel lists) does not exist, and `SpiceError(:unsupported)` if a
meta-kernel lists a meta-kernel. A binary file that is damaged throws
`SpiceError(:baddaf)`, and a meta-kernel that is not consistent throws
`SpiceError(:badkernel)`.

# Examples
```julia
furnsh("naif0012.tls")      # leap seconds
furnsh("de440s.bsp")        # planetary ephemeris
et = str2et("2025-03-01T12:00:00")
state, lt = spkezr("MOON", et, "J2000", "LT+S", "EARTH")

ctx = SpiceContext()        # a second, independent kernel set
furnsh(ctx, "mission.tm")
```

See also [`unload`](@ref), [`kclear`](@ref), [`ktotal`](@ref), [`kdata`](@ref).
"""
function furnsh(ctx::SpiceContext, path::AbstractString)
    isempty(strip(path)) && throw(SpiceError(:badarg, "furnsh: the kernel path is blank"))
    lock(ctx.lock) do
        entries = load_entries(path, (@atomic ctx.timedefs).yearbound)
        list = copy(kernelset(ctx).loaded)
        append!(list, entries)
        _pool_x_swap_loaded!(ctx, list, entries)
    end
    return nothing
end

"""
    unload([ctx,] path) -> Nothing

Remove the kernel that was loaded from `path` out of `ctx`, as `unload_c`.

`path` is compared as a string with the path given to [`furnsh`](@ref). If the path
was loaded more than once, the last entry is removed. A meta-kernel takes the files
that it loaded with it. A path that is not loaded changes nothing, as in CSPICE.

If a text kernel is removed, the kernel pool is built again from the text kernels
that stay, so direct edits of the pool ([`pdpool`](@ref), ...) are lost, as in CSPICE.
The kernels that stay are not read from disk again.

See also [`kclear`](@ref).
"""
function unload(ctx::SpiceContext, path::AbstractString)
    p = String(path)
    lock(ctx.lock) do
        old = kernelset(ctx).loaded
        i = findlast(e -> e.path == p, old)
        i === nothing && return
        list = Loaded[]
        for (j, e) in enumerate(old)
            j == i && continue
            j > i && e.source == p && old[i].source == "" && continue
            push!(list, e)
        end
        _pool_x_swap_unloaded!(ctx, list)
    end
    return nothing
end

"""
    kclear([ctx]) -> Nothing

Unload every kernel of `ctx` and empty its kernel pool, as `kclear_c`.

Direct edits of the pool are removed and the pool watchers are told of the change.
The names of [`boddef`](@ref) stay, as in CSPICE. The time defaults of
[`timdef`](@ref) are not reset.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pdpool(ctx, "A", [1.0]);

julia> kclear(ctx)

julia> expool(ctx, "A")
false
```

See also [`unload`](@ref).
"""
function kclear(ctx::SpiceContext)
    lock(ctx.lock) do
        _pool_x_swap_cleared!(ctx)
    end
    return nothing
end

ktotal(ctx::SpiceContext) = length(kernelset(ctx).loaded)

furnsh(path::AbstractString) = furnsh(default_context(), path)
unload(path::AbstractString) = unload(default_context(), path)
kclear() = kclear(default_context())
ktotal() = ktotal(default_context())

const furnsh! = furnsh
const unload! = unload
const kclear! = kclear
