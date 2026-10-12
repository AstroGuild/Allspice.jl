# Kernel management: the load list, the immutable `KernelSet` snapshot that is built
# from it, and the `SpiceContext` that holds the current snapshot.

"Pool variables that `furnsh` uses for a meta-kernel and then deletes, as CSPICE does."
const METAKERNEL_VARIABLES = ("KERNELS_TO_LOAD", "PATH_SYMBOLS", "PATH_VALUES")

"""
    empty_kernelset() -> KernelSet

Return a `KernelSet` with no kernel loaded; every container in it is empty and
`lsk === nothing`.
"""
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
    filenumber = zeros(Int, length(loaded))
    for (i, entry) in enumerate(loaded)
        content = entry.content
        if content isa DAF
            push!(daf, content.data)
            filenumber[i] = length(daf)
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
        framename2id, tk, bodyrot, lsk, ck, build_sclk(pool), build_ckmeta(pool)),
        edited !== nothing, bodydefs)
end

"""
    SpiceContext()

Create a context whose current snapshot is the empty `KernelSet`.
"""
SpiceContext() = SpiceContext(empty_kernelset(), ReentrantLock())

"The process-wide context of the calls that take no `ctx` argument."
const DEFAULT_CONTEXT = SpiceContext()

"""
    default_context() -> SpiceContext

Return the one process-wide context that the calls without a `ctx` argument use.
"""
default_context() = DEFAULT_CONTEXT

"""
    kernelset(ctx::SpiceContext) -> KernelSet

Return the current snapshot of `ctx` with one atomic read; the result never changes.
"""
kernelset(ctx::SpiceContext) = @atomic ctx.ks

"""
    read_kernel(path) -> Union{DAF,TextKernel}

Open the binary kernel or parse the text kernel at `path`; throws `SpiceError` if the
file does not exist.
"""
function read_kernel(path::AbstractString)
    isfile(path) || throw(SpiceError(:nofile, "kernel file not found: $path"))
    arch = kernel_arch(path)
    return (arch === :SPK || arch === :PCK || arch === :CK) ? open_daf(path) : parse_text_kernel(path)
end

"""
    load_entries(path) -> Vector{Loaded}

Read the kernel at `path` and return its load-list entries: one entry for the kernel
and, for a meta-kernel, one entry for each listed file (in order, `source == path`).
"""
function load_entries(path::AbstractString)
    p = String(path)
    content = read_kernel(p)
    entries = Loaded[Loaded(p, "", content)]
    if content isa TextKernel
        for file in metakernel_files(content)
            sub = read_kernel(file)
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
    furnsh(ctx::SpiceContext, path) -> Nothing

Load the kernel at `path` (SPK, binary PCK, text kernel or meta-kernel) into `ctx`.
The new snapshot replaces the old one atomically; if the load fails, `ctx` is unchanged.
A path that is already loaded gets a second entry, as in CSPICE. Direct edits of the
kernel pool (`pdpool`, ...) stay; the assignments of a text kernel are applied on top.
"""
function furnsh(ctx::SpiceContext, path::AbstractString)
    isempty(strip(path)) && throw(SpiceError(:badarg, "furnsh: the kernel path is blank"))
    lock(ctx.lock) do
        entries = load_entries(path)
        list = copy(kernelset(ctx).loaded)
        append!(list, entries)
        _pool_x_swap_loaded!(ctx, list, entries)
    end
    return nothing
end

"""
    unload(ctx::SpiceContext, path) -> Nothing

Remove the last load-list entry with path `path` and, if it is a meta-kernel, the
entries that it caused. Kernels that stay loaded are not read again. A path that is not
loaded changes nothing, as in CSPICE. If a text kernel is removed, the kernel pool is
built again from the text kernels that stay, so direct edits of the pool (`pdpool`,
...) are lost, as in CSPICE.
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
    kclear(ctx::SpiceContext) -> Nothing

Replace the snapshot of `ctx` with the empty `KernelSet`. Direct edits of the kernel
pool are removed; the names of `boddef` stay, as in CSPICE.
"""
function kclear(ctx::SpiceContext)
    lock(ctx.lock) do
        _pool_x_swap_cleared!(ctx)
    end
    return nothing
end

"""
    ktotal(ctx::SpiceContext) -> Int

Return the number of load-list entries of `ctx` (a meta-kernel counts as one entry plus
one for each file that it lists).
"""
ktotal(ctx::SpiceContext) = length(kernelset(ctx).loaded)

"`furnsh(path)`: load the kernel at `path` into the default context."
furnsh(path::AbstractString) = furnsh(default_context(), path)
"`unload(path)`: unload the kernel at `path` from the default context."
unload(path::AbstractString) = unload(default_context(), path)
"`kclear()`: unload all kernels of the default context."
kclear() = kclear(default_context())
"`ktotal()`: number of load-list entries of the default context."
ktotal() = ktotal(default_context())

"Alias of [`furnsh`](@ref)."
const furnsh! = furnsh
"Alias of [`unload`](@ref)."
const unload! = unload
"Alias of [`kclear`](@ref)."
const kclear! = kclear
