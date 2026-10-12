# Frame information (`frinfo`, `cidfrm`, `cnmfrm`, `ccifrm` of `framex.c`), fixed-offset
# frame data (`tkfram`), body orientation (`tipbod`, `tisbod`), the contents of a binary
# PCK file (`pckfrm`, `pckcov`) and body constants (`bodfnd`, `bodvcd`). Every call has
# a form on a `KernelSet`, on a `SpiceContext` and on the default context.

"""
    frinfo(frcode) -> Union{Nothing,Tuple{Int,Int,Int}}

Center, class and class ID of the frame with ID `frcode`, as `frinfo_c`, or `nothing`
if the frame is not known. A built-in frame has priority over a frame of the kernel
pool, and an ID in 13000:13999 without a specification is a class 2 frame of the Earth
with the class ID `frcode - 10000`.
"""
function frinfo(ks::KernelSet, frcode::Integer)
    (typemin(Int32) <= frcode <= typemax(Int32)) || return nothing
    id = Int32(frcode)
    f = get(BUILTIN_FRAMES, id, nothing)
    f === nothing && (f = get(ks.frames, id, nothing))
    f === nothing || return Int(f.center), Int(f.class), Int(f.classid)
    13000 <= id < 14000 && return 399, 2, Int(id) - 10000
    return nothing
end

"""
    _frame_of_object(ks, value) -> Union{Nothing,Tuple{Int,String}}

Contract: frame ID and frame name for the value of an `OBJECT_<body>_FRAME` variable,
as `cidfrm_` and `cnmfrm_` read it: a frame name (built-in, or with a `FRAME_<name>`
variable) or a frame ID (built-in, or with a `FRAME_<id>_NAME` variable). `nothing` if
the value does not identify a frame.
"""
function _frame_of_object(ks::KernelSet, value::PoolValue)
    isempty(value) && return nothing
    if value isa Vector{String}
        name = String(strip(value[1]))
        code = get(BUILTIN_FRAMENAME2ID, normalize_name(name), nothing)
        code === nothing || return Int(code), name
        v = get(ks.pool, "FRAME_" * name, nothing)
        (v isa Vector{Float64} && length(v) == 1) || return nothing
        return round(Int, v[1]), name
    end
    code = round(Int, value[1])
    f = get(BUILTIN_FRAMES, Int32(clamp(code, typemin(Int32), typemax(Int32))), nothing)
    f === nothing || return code, f.name
    v = get(ks.pool, "FRAME_$(code)_NAME", nothing)
    (v isa Vector{String} && length(v) == 1) || return nothing
    return code, String(strip(v[1]))
end

"""
The built-in frames of `zzfdat_` as `(ID, name, center)`, in the order of the centers
and, for equal centers, of the table (the order vector that `orderi_` gives).
"""
const _FRAME_CENTER_ORDER = let all = Tuple{Int,String,Int}[]
    for (i, d) in enumerate(INERTIAL_FRAME_DEFS)
        push!(all, (i, String(d[1]), 0))
    end
    for f in BUILTIN_NONINERTIAL_FRAMES
        push!(all, (Int(f.id), f.name, Int(f.center)))
    end
    sort!(all; by=x -> x[3], alg=MergeSort)
end

"""
    cidfrm(cent) -> Union{Nothing,Tuple{Int,String}}

ID and name of the frame that is associated with the body with ID `cent`, as
`cidfrm_c`, or `nothing` if there is none. The variable `OBJECT_<cent>_FRAME` (or
`OBJECT_<name>_FRAME`) of the kernel pool has priority; without it the answer is
`IAU_EARTH` for the Earth and otherwise the built-in frame whose center is `cent` (for
a center with several built-in frames, the one that the binary search of CSPICE finds).
"""
function cidfrm(ks::KernelSet, cent::Integer)
    v = get(ks.pool, "OBJECT_$(cent)_FRAME", nothing)
    if v === nothing
        name = body_name(ks, cent)
        name === nothing || (v = get(ks.pool, uppercase("OBJECT_$(name)_FRAME"), nothing))
    end
    v === nothing || return _frame_of_object(ks, v)
    cent == 399 && return 10013, "IAU_EARTH"
    # The binary search of `bschoi_` over the built-in frames in the order of their
    # centers. For a center with several frames (all inertial frames have the center 0)
    # the answer is the frame that this search meets first, as in CSPICE.
    left, right = 1, length(_FRAME_CENTER_ORDER)
    while left <= right     # invariant: a frame with the center `cent` is in left:right, if any
        i = (left + right) ÷ 2
        id, name, center = _FRAME_CENTER_ORDER[i]
        cent == center && return id, name
        cent < center ? (right = i - 1) : (left = i + 1)
    end
    return nothing
end

"""
    cnmfrm(cname) -> Union{Nothing,Tuple{Int,String}}

ID and name of the frame that is associated with the body called `cname`, as
`cnmfrm_c`, or `nothing` if there is none. The variable `OBJECT_<cname>_FRAME` (or
`OBJECT_<id>_FRAME`) of the kernel pool has priority; without it the answer is the
built-in frame `IAU_<CNAME>`, if it exists.
"""
function cnmfrm(ks::KernelSet, cname::AbstractString)
    name = String(strip(cname))
    v = get(ks.pool, uppercase("OBJECT_$(name)_FRAME"), nothing)
    if v === nothing
        code = bodn2c(ks, name)
        code === nothing || (v = get(ks.pool, "OBJECT_$(code)_FRAME", nothing))
    end
    v === nothing || return _frame_of_object(ks, v)
    frname = uppercase("IAU_" * name)
    code = get(BUILTIN_FRAMENAME2ID, normalize_name(frname), nothing)
    return code === nothing ? nothing : (Int(code), frname)
end

"""
    ccifrm(frclss, clssid) -> Union{Nothing,Tuple{Int,String,Int}}

ID, name and center of the frame with the class `frclss` and the class ID `clssid`, as
`ccifrm_c`, or `nothing` if there is none. The built-in frames are searched first, in
the order of their table; then the frames of the kernel pool, of which the one with
the smallest ID is returned if several match (CSPICE returns the first one in the
order of its hash table).
"""
function ccifrm(ks::KernelSet, frclss::Integer, clssid::Integer)
    if frclss == 1 && 1 <= clssid <= N_INERTIAL_FRAMES
        return Int(clssid), INERTIAL_FRAME_DEFS[clssid][1], 0
    end
    for f in BUILTIN_NONINERTIAL_FRAMES
        f.class == frclss && f.classid == clssid && return Int(f.id), f.name, Int(f.center)
    end
    best = nothing
    for f in values(ks.frames)      # invariant: `best` is the match with the smallest ID so far
        (f.class == frclss && f.classid == clssid && !haskey(BUILTIN_FRAMES, f.id)) || continue
        (best === nothing || f.id < best.id) && (best = f)
    end
    return best === nothing ? nothing : (Int(best.id), best.name, Int(best.center))
end

"""
    tkfram(frcode) -> Union{Nothing,Tuple{SMatrix{3,3,Float64,9},Int}}

Data of the fixed-offset (class 4) frame with the class ID `frcode`, as `tkfram_`: the
matrix that rotates vectors from that frame to its relative frame, and the ID of the
relative frame. `nothing` if no valid `TKFRAME_<frcode>_*` definition is loaded (CSPICE
signals an error for an incomplete definition).
"""
function tkfram(ks::KernelSet, frcode::Integer)
    r = ccifrm(ks, 4, frcode)
    r === nothing && return nothing
    tk = get(ks.tk, Int32(r[1]), nothing)
    return tk === nothing ? nothing : (tk.rot, Int(tk.relative))
end

"ID of the built-in inertial frame `ref` (name or ID); throws `SpiceError(:irfnotrec)` otherwise, as `irfrot_`."
function _irf_id(ref::AbstractString)
    id = get(BUILTIN_FRAMENAME2ID, normalize_name(ref), Int32(0))
    1 <= id <= N_INERTIAL_FRAMES || throw(SpiceError(:irfnotrec,
        "the frame '$ref' is not a recognized inertial frame"))
    return id
end
function _irf_id(ref::Integer)
    1 <= ref <= N_INERTIAL_FRAMES || throw(SpiceError(:irfnotrec,
        "the frame with ID $ref is not a recognized inertial frame"))
    return Int32(ref)
end

"""
    _tis_jet(ks, ref, body, et) -> (R, dR)

Contract: the matrix from the inertial frame `ref` to the body-fixed frame of `body`
(a body ID, or the frame class ID of a binary PCK frame) at `et` and its time
derivative: `_body_rotation` (binary PCK first, text constants otherwise) times the
constant rotation from `ref` to J2000.
"""
function _tis_jet(ks::KernelSet, ref::NameOrID, body::Integer, et::Real)
    req = _irf_id(ref)
    m = _body_rotation(ks, Int32(body), taylor_var(check_epoch(et), Val(1)))
    r, dr = map(x -> x.v, m), map(x -> x.d[1], m)
    req == _J2000_FRAME && return r, dr
    c = transpose(inertial_rotation(req))
    return r * c, dr * c
end

"""
    tipbod(ref, body, et) -> SMatrix{3,3}

Matrix that rotates positions from the inertial frame `ref` (name or ID) to the
body-fixed frame of `body` (ID) at ET `et`, as `tipbod_c`. Throws
`SpiceError(:irfnotrec)` if `ref` is not a built-in inertial frame and
`SpiceError(:framedatanotfound)` without PCK data for the body.
"""
tipbod(ks::KernelSet, ref::NameOrID, body::Integer, et::Real) = _tis_jet(ks, ref, body, et)[1]

"""
    tisbod(ref, body, et) -> SMatrix{6,6}

State transformation `[R 0; dR R]` from the inertial frame `ref` (name or ID) to the
body-fixed frame of `body` (ID) at ET `et`, as `tisbod_c`. Throws as `tipbod`.
"""
function tisbod(ks::KernelSet, ref::NameOrID, body::Integer, et::Real)
    r, dr = _tis_jet(ks, ref, body, et)
    return vcat(hcat(r, zero(r)), hcat(dr, r))
end

"""
    pckfrm!(ids, pck) -> ids

Insert the frame class IDs of the frames that have data in the binary PCK file `pck`
into the vector `ids`, which is used as a set (sorted, without duplicates), as `pckfrm_c`
(argument order of SPICE.jl). The file need not be loaded. Throws `SpiceError` as
[`spkobj`](@ref), for a file that is not a binary PCK.
"""
function pckfrm!(ids::AbstractVector{<:Integer}, pck::AbstractString)
    _require_binary_kernel("pckfrm", pck, "PCK")
    sort!(unique!(ids))
    _each_summary(pck) do dc, ic, h
        k = searchsortedfirst(ids, ic[1])
        (k <= length(ids) && ids[k] == ic[1]) || insert!(ids, k, ic[1])
    end
    return ids
end

"""
    pckfrm(pck) -> Vector{Int}

Frame class IDs of the frames that have data in the binary PCK file `pck`, sorted and
without duplicates, as `pckfrm_c`.
"""
pckfrm(pck::AbstractString) = pckfrm!(Int[], pck)

"""
    pckcov!(cover, pck, idcode) -> cover

Insert the coverage of the frame with the class ID `idcode` in the binary PCK file
`pck` into the window `cover` (union with its old content), as `pckcov_c`.
"""
function pckcov!(cover::Window, pck::AbstractString, idcode::Integer)
    _require_binary_kernel("pckcov", pck, "PCK")
    _each_summary(pck) do dc, ic, h
        ic[1] == idcode && wninsd!(cover, dc[1], dc[2])
    end
    return cover
end

"""
    pckcov(pck, idcode) -> Window

Coverage window (ET seconds) of the frame with the class ID `idcode` in the binary PCK
file `pck`, as `pckcov_c`; empty if the frame has no data.
"""
pckcov(pck::AbstractString, idcode::Integer) = pckcov!(Window(), pck, idcode)

"""
    bodfnd(body, item) -> Bool

True if the kernel pool holds the variable `BODY<body>_<item>` (numeric or character),
as `bodfnd_c`. `item` is case sensitive.
"""
bodfnd(ks::KernelSet, body::Integer, item::AbstractString) =
    haskey(ks.pool, string("BODY", body, "_", strip(item)))

"""
    bodvcd(bodyid, item) -> Vector{Float64}

Values of the kernel-pool variable `BODY<bodyid>_<item>`, as `bodvcd_c`. Throws
`SpiceError(:kernelvarnotfound)` if the variable is absent and
`SpiceError(:typemismatch)` if it is character.
"""
function bodvcd(ks::KernelSet, bodyid::Integer, item::AbstractString)
    key = string("BODY", bodyid, "_", strip(item))
    v = get(ks.pool, key, nothing)
    v === nothing && throw(SpiceError(:kernelvarnotfound,
        "the variable $key could not be found in the kernel pool"))
    v isa Vector{Float64} || throw(SpiceError(:typemismatch,
        "the data associated with variable $key are not of numeric type"))
    return copy(v)
end

# Context and default-context forms.
for f in (:frinfo, :cidfrm, :tkfram)
    @eval begin
        $f(ctx::SpiceContext, code::Integer) = $f(kernelset(ctx), code)
        $f(code::Integer) = $f(default_context(), code)
    end
end
cnmfrm(ctx::SpiceContext, cname::AbstractString) = cnmfrm(kernelset(ctx), cname)
cnmfrm(cname::AbstractString) = cnmfrm(default_context(), cname)
ccifrm(ctx::SpiceContext, frclss::Integer, clssid::Integer) = ccifrm(kernelset(ctx), frclss, clssid)
ccifrm(frclss::Integer, clssid::Integer) = ccifrm(default_context(), frclss, clssid)
for f in (:tipbod, :tisbod)
    @eval begin
        $f(ctx::SpiceContext, ref::NameOrID, body::Integer, et::Real) =
            $f(kernelset(ctx), ref, body, et)
        $f(ref::NameOrID, body::Integer, et::Real) = $f(default_context(), ref, body, et)
    end
end
bodfnd(ctx::SpiceContext, body::Integer, item::AbstractString) = bodfnd(kernelset(ctx), body, item)
bodfnd(body::Integer, item::AbstractString) = bodfnd(default_context(), body, item)
bodvcd(ctx::SpiceContext, bodyid::Integer, item::AbstractString) =
    bodvcd(kernelset(ctx), bodyid, item)
bodvcd(bodyid::Integer, item::AbstractString) = bodvcd(default_context(), bodyid, item)
