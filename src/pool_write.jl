# Kernel-pool write interface (`pdpool`, `pcpool`, `pipool`, `dvpool`, `lmpool`,
# `ldpool`, `clpool`), the read calls that `pool_api.jl` does not have (`expool`,
# `dtpool`, `gnpool`, `szpool`, `stpool`, `badkpv`), body name definitions (`boddef`),
# and the watchers (`swpool`, `cvpool`, `poolversion`).
#
# An edit never changes a snapshot. It goes through the context as `furnsh` does: under
# the lock of the context a new pool is made from the current one, a new `KernelSet` is
# built with it (so the names, frames, clocks and every other table that is derived from
# the pool follow the edit), and the reference is swapped atomically. A query that runs
# at the same time works on the old or on the new snapshot, never on a mixture.
#
# Rule for the load list, which is the rule of CSPICE:
#
# - an edit stays through a later `furnsh` (the assignments of a text kernel are applied
#   on top of the edited pool) and through the `unload` of a binary kernel;
# - the `unload` of a text kernel or of a meta-kernel builds the pool again from the
#   text kernels that stay loaded, so every edit is lost;
# - `kclear` removes every edit;
# - `clpool` empties the pool and leaves the load list as it is (a later `unload` of a
#   text kernel brings back the variables of the other text kernels).
#
# Watchers. CSPICE keeps one global agent table. Here each `SpiceContext` has its own
# (`PoolWatch`), guarded by the lock of the context, so `swpool` and `cvpool` are thread
# safe and two contexts do not see the watches of each other. In addition every change
# of the pool of a context increments an atomic counter, `poolversion(ctx)`, which a
# caller can read without a lock and compare with a value that it kept.

"Sizes of `szpool`: the parameters of `pool.c`."
const _POOL_X_SIZES = ("MAXVAR" => 26003, "MAXVAL" => 400000, "MAXLIN" => 15000,
    "MAXCHR" => 80, "MXNOTE" => 130015, "MAXLEN" => 32, "MAXAGT" => 1000)

# ---------------------------------------------------------------------------------------
# Snapshot swaps (called under the lock of the context)
# ---------------------------------------------------------------------------------------

"""
    _pool_x_namepool(pool, edited::Bool) -> AbstractDict

Contract: the pool from which the body names and the frames are built: `pool` itself,
or, for a pool with direct edits whose pair `NAIF_BODY_NAME`, `NAIF_BODY_CODE` is not
complete or not consistent, a copy without the two variables. The two variables are
set by two calls (`pcpool`, then `pipool`), so the state between the calls must give a
valid snapshot; the pair is used as soon as it is consistent. A pool without edits is
returned as it is, and `build_names` reports its errors.
"""
function _pool_x_namepool(pool::Dict{String,PoolValue}, edited::Bool)
    edited || return pool
    try
        build_names(pool)
        return pool
    catch e
        (e isa SpiceError && e.kind === :badpool) || rethrow()
        reduced = copy(pool)
        delete!(reduced, "NAIF_BODY_NAME")
        delete!(reduced, "NAIF_BODY_CODE")
        return reduced
    end
end

"""
    _pool_x_bodydefs!(name2code, code2name, bodydefs) -> Nothing

Contract: add the name and code pairs of `boddef` to the name tables of the kernel
pool. A later pair for the same name replaces an earlier one; a name or a code that
the kernel pool defines keeps the definition of the pool. These tables are searched
before the built-in table, so a `boddef` pair masks a built-in name, as in CSPICE.
"""
function _pool_x_bodydefs!(name2code::Dict{String,Int32}, code2name::Dict{Int32,String},
        bodydefs::Vector{Tuple{String,Int32}})
    isempty(bodydefs) && return nothing
    final = Dict{String,Int32}()
    for (name, code) in bodydefs
        final[normalize_name(name)] = code
    end
    names = Dict{Int32,String}()
    for (name, code) in bodydefs    # invariant: `names` holds the latest valid name per code
        final[normalize_name(name)] == code && (names[code] = name)
    end
    for (nn, code) in final
        haskey(name2code, nn) || (name2code[nn] = code)
    end
    for (code, name) in names
        haskey(code2name, code) || (code2name[code] = name)
    end
    return nothing
end

"""
    _pool_x_notify!(ctx, names) -> Nothing

Contract: called under the lock of `ctx` after a change of the kernel pool. `names` is
the list of the variables that were assigned or deleted, or `:all`. Every agent that
watches one of them gets an update pending, and the pool version of the context
increases by one. Nothing happens for an empty list.
"""
function _pool_x_notify!(ctx::SpiceContext, names)
    names === :all || !isempty(names) || return nothing
    w = ctx.watch
    for (agent, watched) in w.agents
        if names === :all ? !isempty(watched) : any(in(names), watched)
            push!(w.notify, agent)
        end
    end
    @atomic w.version += 1
    return nothing
end

"Names that the text kernels among `entries` assign, with the meta-kernel variables left out."
function _pool_x_assigned(entries::Vector{Loaded})
    names = String[]
    for e in entries
        c = e.content
        c isa TextKernel || continue
        for a in c.assignments
            a.name in METAKERNEL_VARIABLES || push!(names, a.name)
        end
    end
    return names
end

"""
    _pool_x_swap_loaded!(ctx, list, entries) -> Nothing

Contract: called by `furnsh` under the lock of `ctx`. `list` is the new load list and
`entries` are its new entries. Swaps in the snapshot for `list`; if the pool of the old
snapshot has direct edits, the assignments of the new text kernels are applied to that
pool, otherwise the pool is the one of the text kernels of `list`. Notifies the
watchers of the assigned variables.
"""
function _pool_x_swap_loaded!(ctx::SpiceContext, list::Vector{Loaded}, entries::Vector{Loaded})
    old = kernelset(ctx)
    pool = nothing
    if old.pooledited
        pool = copy(old.pool)
        for e in entries
            c = e.content
            c isa TextKernel && apply_assignments!(pool, c)
        end
        for name in METAKERNEL_VARIABLES
            delete!(pool, name)
        end
    end
    ks = build_kernelset(list, pool, old.bodydefs)
    @atomic ctx.ks = ks
    _pool_x_notify!(ctx, _pool_x_assigned(entries))
    return nothing
end

"""
    _pool_x_swap_unloaded!(ctx, list) -> Nothing

Contract: called by `unload` under the lock of `ctx` with the new load list `list`. If
a text kernel was removed, the pool is built again from the text kernels of `list`
(direct edits are lost) and every watcher is notified, as `unload_c` does with
`clpool_` and a reload. Otherwise the pool of the old snapshot is kept.
"""
function _pool_x_swap_unloaded!(ctx::SpiceContext, list::Vector{Loaded})
    old = kernelset(ctx)
    ntext(l) = count(e -> e.content isa TextKernel, l)
    textgone = ntext(list) != ntext(old.loaded)
    pool = (old.pooledited && !textgone) ? old.pool : nothing
    ks = build_kernelset(list, pool, old.bodydefs)
    @atomic ctx.ks = ks
    textgone && _pool_x_notify!(ctx, :all)
    return nothing
end

"""
    _pool_x_swap_cleared!(ctx) -> Nothing

Contract: called by `kclear` under the lock of `ctx`. Swaps in the snapshot without
kernels and without pool edits; the pairs of `boddef` stay. Notifies every watcher.
"""
function _pool_x_swap_cleared!(ctx::SpiceContext)
    ks = build_kernelset(Loaded[], nothing, kernelset(ctx).bodydefs)
    @atomic ctx.ks = ks
    _pool_x_notify!(ctx, :all)
    return nothing
end

"""
    _pool_x_edit!(f!, ctx) -> Nothing

Contract: apply the edit `f!` to the kernel pool of `ctx`. Under the lock of `ctx`,
`f!(pool)` changes a private copy of the pool (it must replace values, never mutate a
stored vector) and returns the names that it touched, or `:all`, or an empty list for
no change. If there is a change, the snapshot with the new pool is swapped in
atomically and the watchers are notified. If `f!` throws, `ctx` is unchanged.
"""
function _pool_x_edit!(f!, ctx::SpiceContext)
    lock(ctx.lock) do
        old = kernelset(ctx)
        pool = copy(old.pool)
        names = f!(pool)
        (names === :all || !isempty(names)) || return
        ks = build_kernelset(old.loaded, pool, old.bodydefs)
        @atomic ctx.ks = ks
        _pool_x_notify!(ctx, names)
    end
    return nothing
end

"Contract: the variable name `name` without outer blanks; throws `SpiceError(:badvarname)` if it is empty or longer than 32 characters, as `pdpool_`."
function _pool_x_name(name::AbstractString)
    n = String(strip(name))
    (1 <= length(n) <= _TK_MAXNAME) || throw(SpiceError(:badvarname, "the kernel pool " *
        "variable name '$n' has length $(length(n)); the allowed length is 1 to $_TK_MAXNAME"))
    return n
end

# ---------------------------------------------------------------------------------------
# Edits
# ---------------------------------------------------------------------------------------

"""
    pdpool([ctx,] name, values) -> Nothing

Set the kernel-pool variable `name` to the numbers `values`, as `pdpool_c`.

A variable of that name is replaced, whatever its type. An empty `values` changes
nothing, as in CSPICE. Blanks around `name` are ignored.

Every table that is built from the pool (body names, frames, clocks) follows the
edit. The edit stays through a later
[`furnsh`](@ref) and through the [`unload`](@ref) of a binary kernel. It is lost at
[`kclear`](@ref) and at the `unload` of a text kernel or a meta-kernel, as in CSPICE.

Throws `SpiceError(:badvarname)` if `name` is blank or longer than 32 characters.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pdpool(ctx, "BODY999_GM", [42.0]);

julia> gdpool(ctx, "BODY999_GM")
1-element Vector{Float64}:
 42.0
```

See also [`pipool`](@ref), [`pcpool`](@ref), [`dvpool`](@ref), [`lmpool`](@ref).
"""
function pdpool(ctx::SpiceContext, name::AbstractString, values::AbstractVector{<:Real})
    n = _pool_x_name(name)
    isempty(values) && return nothing
    v = Float64[x for x in values]
    return _pool_x_edit!(pool -> (pool[n] = v; (n,)), ctx)
end

"""
    pipool([ctx,] name, ivals) -> Nothing

Set the kernel-pool variable `name` to the integers `ivals`, as `pipool_c`.

The pool stores them as `Float64` numbers; [`gipool`](@ref) reads them back as
integers. The rules and the error are those of [`pdpool`](@ref).
"""
pipool(ctx::SpiceContext, name::AbstractString, ivals::AbstractVector{<:Integer}) =
    pdpool(ctx, name, ivals)

"""
    pcpool([ctx,] name, cvals) -> Nothing

Set the kernel-pool variable `name` to the strings `cvals`, as `pcpool_c`.

A string longer than 80 characters is cut to 80, as in the kernel pool of CSPICE. The
rules and the error are those of [`pdpool`](@ref).
"""
function pcpool(ctx::SpiceContext, name::AbstractString, cvals::AbstractVector{<:AbstractString})
    n = _pool_x_name(name)
    isempty(cvals) && return nothing
    v = String[first(String(s), _TK_MAXSTR) for s in cvals]
    return _pool_x_edit!(pool -> (pool[n] = v; (n,)), ctx)
end

"""
    dvpool([ctx,] name) -> Nothing

Delete the kernel-pool variable `name`, as `dvpool_c`.

A name that is not in the pool changes nothing. The deletion is a direct edit of the
pool, with the rules of [`pdpool`](@ref).
"""
function dvpool(ctx::SpiceContext, name::AbstractString)
    n = String(strip(name))
    return _pool_x_edit!(ctx) do pool
        haskey(pool, n) || return ()
        delete!(pool, n)
        return (n,)
    end
end

"""
    clpool([ctx]) -> Nothing

Remove every variable from the kernel pool, as `clpool_c`.

The loaded kernels are not changed: binary kernels stay usable, and the text kernels
stay in the list of [`ktotal`](@ref) although their variables are gone. A later
[`unload`](@ref) of a text kernel builds the pool again from the text kernels that
stay loaded.

See also [`kclear`](@ref).
"""
clpool(ctx::SpiceContext) = _pool_x_edit!(pool -> (empty!(pool); :all), ctx)

function _pool_x_apply!(pool::Dict{String,PoolValue}, tk::TextKernel)
    apply_assignments!(pool, tk)
    return String[a.name for a in tk.assignments]
end

"""
    lmpool([ctx,] cvals) -> Nothing

Load the text lines `cvals` into the kernel pool as the data section of a text kernel,
as `lmpool_c`.

Each line holds text-kernel assignments (`NAME = value`, `NAME += value`); an
assignment can continue on the next line. The lines start in a data section, and
`\\begintext` and `\\begindata` work as in a file. The variables are direct edits of
the pool, with the rules of [`pdpool`](@ref).

Throws `SpiceError(:badkernel)` if an assignment is malformed; the pool is then
unchanged.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> lmpool(ctx, ["A = 1.5", "B = ( 'x', 'y' )"]);

julia> gdpool(ctx, "A"), gcpool(ctx, "B")
([1.5], ["x", "y"])
```

See also [`ldpool`](@ref), [`furnsh`](@ref).
"""
function lmpool(ctx::SpiceContext, cvals::AbstractVector{<:AbstractString})
    defs = timedefaults(ctx)
    tk = parse_text_lines("lmpool", vcat(["\\begindata"], String[String(s) for s in cvals]),
        defs.yearbound, defs.check)
    return _pool_x_edit!(pool -> _pool_x_apply!(pool, tk), ctx)
end

"""
    ldpool([ctx,] fname) -> Nothing

Load the variables of the text kernel file `fname` into the kernel pool, as
`ldpool_c`.

Different from [`furnsh`](@ref), the file does not enter the list of loaded kernels:
[`ktotal`](@ref) does not count it and [`unload`](@ref) cannot remove it. Its
variables are direct edits of the pool, with the rules of [`pdpool`](@ref). A
`KERNELS_TO_LOAD` list in the file is not followed.

Throws `SpiceError(:emptystring)` if `fname` is empty, `SpiceError(:fileopenfailed)` if the
file does not exist and `SpiceError(:badkernel)` if an assignment is malformed.
"""
function ldpool(ctx::SpiceContext, fname::AbstractString)
    isempty(fname) && throw(SpiceError(:emptystring, "ldpool: the file name is empty"))
    isfile(fname) || throw(SpiceError(:fileopenfailed, "ldpool: could not open $fname"))
    defs = timedefaults(ctx)
    tk = parse_text_kernel(String(fname), defs.yearbound, defs.check)
    return _pool_x_edit!(pool -> _pool_x_apply!(pool, tk), ctx)
end

"""
    boddef([ctx,] name, code) -> Nothing

Define `name` as a name of the body with the NAIF ID `code`, as `boddef_c`.

The pair has priority over the built-in names. A name or a code that the kernel pool
defines (`NAIF_BODY_NAME`, `NAIF_BODY_CODE`) keeps the definition of the pool. The
pair stays through [`furnsh`](@ref), [`unload`](@ref) and [`kclear`](@ref), as in
CSPICE. A lookup ignores letter case and blanks around the name;
[`bodc2n`](@ref) returns the name as given here.

Throws `SpiceError(:badarg)` if `name` is blank.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> boddef(ctx, "My Probe", -999);

julia> bodn2c(ctx, "my probe"), bodc2n(ctx, -999)
(-999, "My Probe")
```
"""
function boddef(ctx::SpiceContext, name::AbstractString, code::Integer)
    nm = String(strip(name))
    isempty(nm) && throw(SpiceError(:badarg, "boddef: the body name is blank"))
    lock(ctx.lock) do
        old = kernelset(ctx)
        defs = copy(old.bodydefs)
        push!(defs, (nm, Int32(code)))
        ks = build_kernelset(old.loaded, old.pooledited ? old.pool : nothing, defs)
        @atomic ctx.ks = ks
    end
    return nothing
end

# ---------------------------------------------------------------------------------------
# Reads
# ---------------------------------------------------------------------------------------

"""
    expool([ctx,] name) -> Bool

Return `true` if `name` is a numeric variable of the kernel pool, as `expool_c`.

A character variable gives `false`, as in CSPICE. Use [`dtpool`](@ref) to find a
variable of any type.
"""
expool(ks::KernelSet, name::AbstractString) =
    get(ks.pool, String(strip(name)), nothing) isa Vector{Float64}

"""
    dtpool([ctx,] name) -> (n::Int, type::Symbol)

Return the number of values and the type of the kernel-pool variable `name`, as
`dtpool_c`.

`type` is `:N` for numbers and `:C` for strings. A variable that is not in the pool
gives `(0, :X)`; this takes the place of the `found` flag of CSPICE, as in SPICE.jl.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pcpool(ctx, "NAMES", ["a", "b"]);

julia> dtpool(ctx, "NAMES"), dtpool(ctx, "ABSENT")
((2, :C), (0, :X))
```
"""
function dtpool(ks::KernelSet, name::AbstractString)
    v = get(ks.pool, String(strip(name)), nothing)
    v === nothing && return 0, :X
    return length(v), v isa Vector{Float64} ? :N : :C
end

"""
    badkpv([ctx,] caller, name, comp, size, divby, type) -> Bool

Check the number of values and the type of the kernel-pool variable `name`, as
`badkpv_c`, and return `false` if they are as expected.

The number of values `n` must satisfy `n comp size`, where `comp` is one of the
operators of CSPICE, `"="`, `"<"`, `">"`, `"<="` and `"=>"` (not `">="` or `"=<"`);
trailing blanks of `comp` are ignored. `n` must also be a
multiple of `divby`; `divby = 0` fails for every variable, as in CSPICE. `type` is
`'N'` for numbers or `'C'` for strings; letter case is ignored. `caller` is the name
that starts the error messages. Blanks around `name` are ignored, as in
[`dtpool`](@ref); CSPICE ignores only the trailing ones.

CSPICE signals an error and returns `true` for a bad variable; this call throws, so
the result is always `false`.

Throws `SpiceError`:
- `:emptystring` if `caller`, `name` or `comp` is empty;
- `:variablenotfound` if the variable is not in the pool;
- `:unknowncompare` if `comp` is not one of the above;
- `:badvariablesize` if the number of values fails the comparison or is not a multiple
  of `divby`;
- `:badvariabletype` if the variable has the other type;
- `:invalidtype` if `type` is not `'N'` or `'C'`.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pdpool(ctx, "BODY399_RADII", [6378.1366, 6378.1366, 6356.7519]);

julia> badkpv(ctx, "MYPROG", "BODY399_RADII", "=", 3, 1, 'N')
false

julia> try
           badkpv(ctx, "MYPROG", "BODY399_RADII", "=>", 6, 1, 'N')
       catch err
           err.kind
       end
:badvariablesize
```
"""
function badkpv(ks::KernelSet, caller::AbstractString, name::AbstractString,
        comp::AbstractString, size::Integer, divby::Integer, type::AbstractChar)
    for (arg, s) in (("caller", caller), ("name", name), ("comp", comp))
        isempty(s) && throw(SpiceError(:emptystring, "badkpv: the argument $arg is empty"))
    end
    dim, class = dtpool(ks, name)
    class === :X && throw(SpiceError(:variablenotfound, "$caller: the kernel pool " *
        "variable '$name' is not currently present in the kernel pool"))
    op = rstrip(comp, ' ')
    ok = op == "=" ? dim == size : op == "<" ? dim < size : op == ">" ? dim > size :
         op == "<=" ? dim <= size : op == "=>" ? dim >= size :
         throw(SpiceError(:unknowncompare, "$caller: the comparison operator '$comp' is " *
            "not a recognized value; the recognized values are '<', '<=', '=', '=>', '>'"))
    ok || throw(SpiceError(:badvariablesize, "$caller: the kernel pool variable '$name' " *
        "is expected to have a number of components DIM such that DIM $op $size, but it " *
        "has $dim"))
    # The rule of CSPICE: a division by zero is replaced by a ratio of 1.
    ratio = divby != 0 ? div(dim, divby) : 1
    divby * ratio == dim || throw(SpiceError(:badvariablesize, "$caller: the number of " *
        "components of the kernel pool variable '$name' is $dim, which is not divisible by $divby"))
    t = uppercase(type)
    (t == 'C' || t == 'N') || throw(SpiceError(:invalidtype, "$caller: the expected type " *
        "'$type' of the kernel pool variable '$name' is neither 'C' nor 'N'"))
    want = t == 'C' ? :C : :N
    class === want || throw(SpiceError(:badvariabletype, "$caller: the kernel pool " *
        "variable '$name' must be of type " * (t == 'C' ? "\"CHARACTER\"; it is numeric" :
        "\"NUMERIC\"; it is character")))
    return false
end

"""
    _pool_x_match(s, t) -> Bool

Contract: true if the string `s` matches the template `t`, as `matchi_` with the
wild cards of `gnpool_`: `*` stands for any substring, `%` for any one character, and
letter case does not matter.
"""
function _pool_x_match(s::AbstractString, t::AbstractString)
    a, b = collect(uppercase(s)), collect(uppercase(t))
    i, j, star, mark = 1, 1, 0, 0
    while i <= length(a)    # invariant: a[1:i-1] matches b[1:j-1]
        if j <= length(b) && (b[j] == '%' || (b[j] != '*' && b[j] == a[i]))
            i += 1
            j += 1
        elseif j <= length(b) && b[j] == '*'
            star, mark = j, i
            j += 1
        elseif star != 0
            j = star + 1
            mark += 1
            i = mark
        else
            return false
        end
    end
    while j <= length(b) && b[j] == '*'
        j += 1
    end
    return j > length(b)
end

"""
    gnpool([ctx,] name, start=1, room=typemax(Int)) -> Vector{String}

Return the names of the kernel-pool variables that match the template `name`, as
`gnpool_c`.

In the template `*` stands for any substring and `%` for any one character; letter
case is ignored. The result holds at most `room` names, from match number `start`.
`start` counts from 1, as in SPICE.jl; CSPICE counts from 0. The matches are in
ascending order of the name; CSPICE returns them in the order of its hash table. The
result is empty if nothing matches; there is no `found` flag.

Throws `SpiceError(:badarraysize)` if `room < 1`.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pdpool(ctx, "BODY399_GM", [398600.4]); pdpool(ctx, "BODY301_GM", [4902.8]);

julia> gnpool(ctx, "BODY%%%_GM")
2-element Vector{String}:
 "BODY301_GM"
 "BODY399_GM"
```
"""
function gnpool(ks::KernelSet, name::AbstractString, start::Integer=1, room::Integer=typemax(Int))
    room >= 1 || throw(SpiceError(:badarraysize, "the amount of room specified as " *
        "available for output was $room; the amount of room must be positive"))
    t = String(strip(name))
    hits = sort!(String[k for k in keys(ks.pool) if _pool_x_match(k, t)])
    first = max(1, Int(start))
    first > length(hits) && return String[]
    return hits[first:min(length(hits), first + min(Int(room), length(hits)) - 1)]
end

"""
    szpool(name) -> Union{Int,Nothing}

Return the kernel-pool size parameter `name` of CSPICE N0067, or `nothing` if the name
is not one, as `szpool_c`.

The names are `MAXVAR`, `MAXVAL`, `MAXLIN`, `MAXCHR`, `MXNOTE`, `MAXLEN` and `MAXAGT`;
letter case and blanks are ignored. Of these limits the package has only `MAXLEN`
(32 characters in a variable name) and `MAXCHR` (80 characters in a string value). The
other values are reported for programs that read them and limit nothing here.
"""
function szpool(name::AbstractString)
    key = uppercase(replace(name, ' ' => ""))
    for (k, v) in _POOL_X_SIZES
        k == key && return v
    end
    return nothing
end

"""
    stpool([ctx,] item, nth, contin) -> Union{String,Nothing}

Return string number `nth` of the character kernel-pool variable `item`, with
continued elements joined, as `stpool_c`.

An element that ends with the continuation marker `contin` is joined with the next
element; the marker is removed and the text before it is kept. `nth` counts the joined
strings from 1, as in SPICE.jl; CSPICE counts from 0. The result is `nothing` if the
variable is absent, is numeric, or has fewer strings; this takes the place of the
`found` flag.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pcpool(ctx, "PATH", ["/data/kernels/+", "spk/de440s.bsp", "other"]);

julia> stpool(ctx, "PATH", 1, "+"), stpool(ctx, "PATH", 2, "+"), stpool(ctx, "PATH", 3, "+")
("/data/kernels/spk/de440s.bsp", "other", nothing)
```
"""
function stpool(ks::KernelSet, item::AbstractString, nth::Integer, contin::AbstractString)
    v = get(ks.pool, String(strip(item)), nothing)
    v isa Vector{String} || return nothing
    mark = rstrip(contin, ' ')
    isempty(mark) && return 1 <= nth <= length(v) ? String(rstrip(v[nth], ' ')) : nothing
    return pool_string_joined(v, nth, mark)
end

# ---------------------------------------------------------------------------------------
# Watchers
# ---------------------------------------------------------------------------------------

"""
    swpool([ctx,] agent, names) -> Nothing

Watch the kernel-pool variables `names` for the agent `agent`, as `swpool_c`. An
agent is a name that the caller chooses for one watcher.

After this call, and after every later change of one of the variables,
[`cvpool`](@ref) returns `true` once for `agent`. The names are added to the ones
that the agent watches already. A watch belongs to one context: CSPICE keeps one
table for the process, and here two contexts do not see the watches of each other.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> swpool(ctx, "MY_AGENT", ["BODY399_GM"]);

julia> cvpool(ctx, "MY_AGENT"), cvpool(ctx, "MY_AGENT")
(true, false)

julia> pdpool(ctx, "BODY399_GM", [398600.4]);

julia> cvpool(ctx, "MY_AGENT")
true
```

See also [`poolversion`](@ref).
"""
function swpool(ctx::SpiceContext, agent::AbstractString, names::AbstractVector{<:AbstractString})
    a = String(strip(agent))
    lock(ctx.lock) do
        watched = get!(() -> String[], ctx.watch.agents, a)
        for n in names
            s = String(strip(n))
            s in watched || push!(watched, s)
        end
        push!(ctx.watch.notify, a)
    end
    return nothing
end

"""
    cvpool([ctx,] agent) -> Bool

Return `true` if a variable that `agent` watches changed since the last call for that
agent, or since [`swpool`](@ref), as `cvpool_c`.

The call clears the pending update. A change is an assignment by [`furnsh`](@ref),
[`ldpool`](@ref), [`lmpool`](@ref), [`pdpool`](@ref), [`pcpool`](@ref) or
[`pipool`](@ref), or a deletion by [`dvpool`](@ref). [`clpool`](@ref),
[`kclear`](@ref) and the [`unload`](@ref) of a text kernel count as a change of every
watched variable, as in CSPICE. An agent without a watch gives `false`.
"""
function cvpool(ctx::SpiceContext, agent::AbstractString)
    a = String(strip(agent))
    return lock(ctx.lock) do
        a in ctx.watch.notify || return false
        delete!(ctx.watch.notify, a)
        return true
    end
end

"""
    poolversion([ctx]) -> Int

Return the number of changes of the kernel pool of `ctx` since the context was made.
This call has no CSPICE counterpart.

Keep the value and compare it later to find out whether the pool changed. The call
takes no lock; [`cvpool`](@ref) takes the lock of the context, because it clears the
pending update. A change of the load list that assigns no variable (a binary kernel)
does not count.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> v = poolversion(ctx);

julia> pdpool(ctx, "A", [1.0]);

julia> poolversion(ctx) - v
1
```

See also [`swpool`](@ref), [`cvpool`](@ref).
"""
poolversion(ctx::SpiceContext) = @atomic ctx.watch.version

# ---------------------------------------------------------------------------------------
# Context and default-context forms
# ---------------------------------------------------------------------------------------

pdpool(name::AbstractString, values::AbstractVector{<:Real}) = pdpool(default_context(), name, values)
pipool(name::AbstractString, ivals::AbstractVector{<:Integer}) = pipool(default_context(), name, ivals)
pcpool(name::AbstractString, cvals::AbstractVector{<:AbstractString}) =
    pcpool(default_context(), name, cvals)
dvpool(name::AbstractString) = dvpool(default_context(), name)
clpool() = clpool(default_context())
lmpool(cvals::AbstractVector{<:AbstractString}) = lmpool(default_context(), cvals)
ldpool(fname::AbstractString) = ldpool(default_context(), fname)
boddef(name::AbstractString, code::Integer) = boddef(default_context(), name, code)
swpool(agent::AbstractString, names::AbstractVector{<:AbstractString}) =
    swpool(default_context(), agent, names)
cvpool(agent::AbstractString) = cvpool(default_context(), agent)
poolversion() = poolversion(default_context())
for f in (:expool, :dtpool)
    @eval begin
        $f(ctx::SpiceContext, name::AbstractString) = $f(kernelset(ctx), name)
        $f(name::AbstractString) = $f(default_context(), name)
    end
end
gnpool(ctx::SpiceContext, name::AbstractString, start::Integer=1, room::Integer=typemax(Int)) =
    gnpool(kernelset(ctx), name, start, room)
gnpool(name::AbstractString, start::Integer=1, room::Integer=typemax(Int)) =
    gnpool(default_context(), name, start, room)
stpool(ctx::SpiceContext, item::AbstractString, nth::Integer, contin::AbstractString) =
    stpool(kernelset(ctx), item, nth, contin)
stpool(item::AbstractString, nth::Integer, contin::AbstractString) =
    stpool(default_context(), item, nth, contin)
badkpv(ctx::SpiceContext, caller::AbstractString, name::AbstractString, comp::AbstractString,
    size::Integer, divby::Integer, type::AbstractChar) =
    badkpv(kernelset(ctx), caller, name, comp, size, divby, type)
badkpv(caller::AbstractString, name::AbstractString, comp::AbstractString, size::Integer,
    divby::Integer, type::AbstractChar) =
    badkpv(default_context(), caller, name, comp, size, divby, type)
