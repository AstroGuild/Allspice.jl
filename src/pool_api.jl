# Kernel-pool, name, frame-name and time interface on a context.

"Return the value of pool variable `name`; throws `SpiceError(:notfound)` if it is absent."
function pool_value(ks::KernelSet, name::AbstractString)
    v = get(ks.pool, String(strip(name)), nothing)
    v === nothing && throw(SpiceError(:notfound, "kernel-pool variable $name is not loaded"))
    return v
end

"""
    gdpool(name) -> Vector{Float64}

Return a copy of the numeric values of pool variable `name`. Throws
`SpiceError(:notfound)` if it is absent and `SpiceError(:badtype)` if it is character.
"""
function gdpool(ks::KernelSet, name::AbstractString)
    v = pool_value(ks, name)
    v isa Vector{Float64} || throw(SpiceError(:badtype, "kernel-pool variable $name is not numeric"))
    return copy(v)
end

"""
    gipool(name) -> Vector{Int}

Return the values of numeric pool variable `name` rounded to integers, as `gipool_c`.
"""
gipool(ks::KernelSet, name::AbstractString) = Int[round(Int, x) for x in gdpool(ks, name)]

"""
    gcpool(name) -> Vector{String}

Return a copy of the character values of pool variable `name`. Throws
`SpiceError(:notfound)` if it is absent and `SpiceError(:badtype)` if it is numeric.
"""
function gcpool(ks::KernelSet, name::AbstractString)
    v = pool_value(ks, name)
    v isa Vector{String} || throw(SpiceError(:badtype, "kernel-pool variable $name is not character"))
    return copy(v)
end

"""
    bodvrd(body, item) -> Vector{Float64}

Return the values of pool variable `BODY<id>_<ITEM>` for `body` (name or ID), as
`bodvrd_c`; `item` is case sensitive, as in CSPICE.
"""
bodvrd(ks::KernelSet, body::NameOrID, item::AbstractString) =
    gdpool(ks, string("BODY", body_code(ks, body), "_", strip(item)))

"""
    bods2c(name) -> Union{Int,Nothing}

Return the ID code of the body `name` (a body name or a string that holds an integer),
or `nothing` if the name is unknown.
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
    bodn2c(name) -> Union{Int,Nothing}

Return the ID code of the body name `name`, or `nothing` if the name is unknown. A
string that holds an integer is not a name, as in `bodn2c_c`.
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
    bodc2n(code) -> Union{String,Nothing}

Return the name of the body with ID `code`, or `nothing` if it has no name.
"""
bodc2n(ks::KernelSet, code::Integer) = body_name(ks, code)

"""
    namfrm(name) -> Int

Return the ID of the frame `name`, or 0 if the frame is unknown, as `namfrm_c`.
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
    frmnam(id) -> String

Return the name of the frame with ID `id`, or `""` if the frame is unknown, as
`frmnam_c`.
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

# Time conversions (`time.jl` defines the forms on a `KernelSet`).
deltet(ctx::SpiceContext, epoch::Real, eptype::AbstractString) =
    deltet(kernelset(ctx), epoch, eptype)
deltet(epoch::Real, eptype::AbstractString) = deltet(default_context(), epoch, eptype)
utc2et(ctx::SpiceContext, s::AbstractString) = utc2et(kernelset(ctx), s)
utc2et(s::AbstractString) = utc2et(default_context(), s)
str2et(ctx::SpiceContext, s::AbstractString) = str2et(kernelset(ctx), s)
str2et(s::AbstractString) = str2et(default_context(), s)
et2utc(ctx::SpiceContext, et::Real, format::AbstractString, prec::Integer) =
    et2utc(kernelset(ctx), et, format, prec)
et2utc(et::Real, format::AbstractString, prec::Integer) =
    et2utc(default_context(), et, format, prec)
