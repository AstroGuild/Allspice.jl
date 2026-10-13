# Kernel-pool, name, frame-name and time interface on a context.

"Return the value of pool variable `name`; throws `SpiceError(:notfound)` if it is absent."
function pool_value(ks::KernelSet, name::AbstractString)
    v = get(ks.pool, String(strip(name)), nothing)
    v === nothing && throw(SpiceError(:notfound, "kernel-pool variable $name is not loaded"))
    return v
end

"""
    gdpool([ctx,] name) -> Vector{Float64}

Return the numeric values of the kernel-pool variable `name`, as `gdpool_c`.

The result is a new vector with every value of the variable. Different from CSPICE,
the call has no `start` and `room` arguments and no `found` flag: an absent variable
is an error. Use [`expool`](@ref) or [`dtpool`](@ref) to test first. Blanks around
`name` are ignored; letter case is not.

Throws `SpiceError(:notfound)` if the variable is not in the pool and
`SpiceError(:badtype)` if it holds character values.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pdpool(ctx, "BODY399_RADII", [6378.1366, 6378.1366, 6356.7519]);

julia> gdpool(ctx, "BODY399_RADII")
3-element Vector{Float64}:
 6378.1366
 6378.1366
 6356.7519
```

See also [`gipool`](@ref), [`gcpool`](@ref), [`bodvrd`](@ref), [`pdpool`](@ref).
"""
function gdpool(ks::KernelSet, name::AbstractString)
    v = pool_value(ks, name)
    v isa Vector{Float64} || throw(SpiceError(:badtype, "kernel-pool variable $name is not numeric"))
    return copy(v)
end

"""
    gipool([ctx,] name) -> Vector{Int}

Return the values of the numeric kernel-pool variable `name` rounded to the nearest
integers, as `gipool_c`.

The result holds every value; there is no `found` flag (see [`gdpool`](@ref)).

Throws `SpiceError(:notfound)` if the variable is not in the pool and
`SpiceError(:badtype)` if it holds character values.
"""
gipool(ks::KernelSet, name::AbstractString) = Int[round(Int, x) for x in gdpool(ks, name)]

"""
    gcpool([ctx,] name) -> Vector{String}

Return the character values of the kernel-pool variable `name`, as `gcpool_c`.

The result is a new vector with every value; there is no `found` flag (see
[`gdpool`](@ref)).

Throws `SpiceError(:notfound)` if the variable is not in the pool and
`SpiceError(:badtype)` if it holds numeric values.
"""
function gcpool(ks::KernelSet, name::AbstractString)
    v = pool_value(ks, name)
    v isa Vector{String} || throw(SpiceError(:badtype, "kernel-pool variable $name is not character"))
    return copy(v)
end

"""
    bodvrd([ctx,] body, item) -> Vector{Float64}

Return the values of the kernel-pool variable `BODY<id>_<item>` for `body`, as
`bodvrd_c`.

`body` is a body name or a NAIF ID code. `item` is case sensitive, as in CSPICE;
blanks around it are ignored. The values have the units of the kernel: km for
`"RADII"`, km^3/s^2 for `"GM"` in the NAIF kernels. There is no `maxn` argument.

Throws `SpiceError(:unknownbody)` if the body name is not known,
`SpiceError(:notfound)` if the variable is not in the pool and `SpiceError(:badtype)`
if it holds character values.

# Examples
```julia
furnsh("pck00011.tpc")
a, b, c = bodvrd("EARTH", "RADII")    # km
```

See also [`bodvcd`](@ref), [`gdpool`](@ref).
"""
bodvrd(ks::KernelSet, body::NameOrID, item::AbstractString) =
    gdpool(ks, string("BODY", body_code(ks, body), "_", strip(item)))

"""
    bods2c([ctx,] name) -> Union{Int,Nothing}

Return the NAIF ID code of the body `name`, or `nothing` if the name is not known, as
`bods2c_c`.

`name` is a body name or a string that holds an integer, such as `"399"`. Letter case
and blanks around the name are ignored. Names from the kernel pool and from
[`boddef`](@ref) have priority over the built-in names. `nothing` takes the place of
the `found` flag of CSPICE.

# Examples
```jldoctest
julia> bods2c("Earth"), bods2c("-77"), bods2c("NO SUCH BODY")
(399, -77, nothing)
```

See also [`bodn2c`](@ref), [`bodc2n`](@ref).
"""
function bods2c(ks::KernelSet, name::AbstractString)
    try
        return Int(body_code(ks, name))
    catch e
        e isa SpiceError && e.kind === :unknownbody && return nothing
        rethrow()
    end
end

"""
    bodn2c([ctx,] name) -> Union{Int,Nothing}

Return the NAIF ID code of the body name `name`, or `nothing` if the name is not
known, as `bodn2c_c`.

Different from [`bods2c`](@ref), a string that holds an integer is not a name.
`nothing` takes the place of the `found` flag of CSPICE.

# Examples
```jldoctest
julia> bodn2c("MOON"), bodn2c("301")
(301, nothing)
```
"""
function bodn2c(ks::KernelSet, name::AbstractString)
    code = bods2c(ks, name)
    code === nothing && return nothing
    # `body_code` also accepts integer strings; accept one only if it is a defined name.
    if tryparse(Int, strip(name)) !== nothing
        haskey(ks.name2code, normalize_name(name)) || return nothing
    end
    return code
end

"""
    bodc2n([ctx,] code) -> Union{String,Nothing}

Return the name of the body with the NAIF ID `code`, or `nothing` if the code has no
name, as `bodc2n_c`.

A name from the kernel pool or from [`boddef`](@ref) has priority over the built-in
name.

# Examples
```jldoctest
julia> bodc2n(399), bodc2n(123456789)
("EARTH", nothing)
```

See also [`bodn2c`](@ref).
"""
bodc2n(ks::KernelSet, code::Integer) = body_name(ks, code)

"""
    namfrm([ctx,] name) -> Int

Return the ID of the frame `name`, or 0 if the frame is not known, as `namfrm_c`.

Built-in frames and frames defined in the kernel pool are known.

# Examples
```jldoctest
julia> namfrm("J2000"), namfrm("eclipj2000"), namfrm("NO_SUCH_FRAME")
(1, 17, 0)
```

See also [`frmnam`](@ref), [`frinfo`](@ref).
"""
function namfrm(ks::KernelSet, name::AbstractString)
    try
        return Int(frame_id(ks, name))
    catch e
        e isa SpiceError && e.kind === :unknownframe && return 0
        rethrow()
    end
end

"""
    frmnam([ctx,] id) -> String

Return the name of the frame with ID `id`, or `""` if the frame is not known, as
`frmnam_c`.

# Examples
```jldoctest
julia> frmnam(17), frmnam(-123456)
("ECLIPJ2000", "")
```

See also [`namfrm`](@ref).
"""
function frmnam(ks::KernelSet, id::Integer)
    name = frame_name(ks, id)
    return name === nothing ? "" : name
end

# Context and default-context forms of the functions above.
for f in (:gdpool, :gipool, :gcpool, :bods2c, :bodn2c, :namfrm)
    @eval begin
        $f(ctx::SpiceContext, name::AbstractString) = $f(kernelset(ctx), name)
        $f(name::AbstractString) = $f(default_context(), name)
    end
end
for f in (:bodc2n, :frmnam)
    @eval begin
        $f(ctx::SpiceContext, code::Integer) = $f(kernelset(ctx), code)
        $f(code::Integer) = $f(default_context(), code)
    end
end
bodvrd(ctx::SpiceContext, body::NameOrID, item::AbstractString) =
    bodvrd(kernelset(ctx), body, item)
bodvrd(body::NameOrID, item::AbstractString) = bodvrd(default_context(), body, item)
