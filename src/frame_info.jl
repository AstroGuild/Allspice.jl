# Frame information (`frinfo`, `cidfrm`, `cnmfrm`, `ccifrm` of `framex.c`), the lists of
# frames (`bltfrm`, `kplfrm`), fixed-offset frame data (`tkfram`), body orientation
# (`tipbod`, `tisbod`), the contents of a binary PCK file (`pckfrm`, `pckcov`) and body
# constants (`bodfnd`, `bodvcd`). Every call that reads kernel data has a form on a
# `KernelSet`, on a `SpiceContext` and on the default context.

"""
    frinfo([ctx,] frcode) -> Union{Nothing,Tuple{Int,Int,Int}}

Return `(center, class, classid)` of the frame with ID `frcode`, as `frinfo_c`, or
`nothing` if the frame is not known.

`center` is the NAIF ID of the body at the frame origin. `class` is 1 for inertial, 2
for PCK, 3 for CK, 4 for fixed offset, 5 for dynamic and 6 for switch frames;
`classid` is the ID of the frame inside its class. A built-in frame has priority over
a frame of the kernel pool. An ID in 13000:13999 without a definition is a PCK frame
of the Earth with the class ID `frcode - 10000`, as in CSPICE. `nothing` takes the
place of the `found` flag.

# Examples
```jldoctest
julia> frinfo(namfrm("IAU_EARTH")), frinfo(-999)
((399, 2, 399), nothing)
```

See also [`namfrm`](@ref), [`ccifrm`](@ref).
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
    cidfrm([ctx,] cent) -> Union{Nothing,Tuple{Int,String}}

Return `(frcode, frname)` of the frame that is associated with the body with ID
`cent`, as `cidfrm_c`, or `nothing` if there is none.

The kernel-pool variable `OBJECT_<cent>_FRAME` (or `OBJECT_<name>_FRAME`) has
priority. Without it the answer is `IAU_EARTH` for the Earth, and otherwise a built-in
frame whose center is `cent`. For a center with several built-in frames the answer is
the frame that CSPICE gives. `nothing` takes the place of the `found` flag.

# Examples
```jldoctest
julia> cidfrm(499)
(10014, "IAU_MARS")
```

See also [`cnmfrm`](@ref).
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
    cnmfrm([ctx,] cname) -> Union{Nothing,Tuple{Int,String}}

Return `(frcode, frname)` of the frame that is associated with the body called
`cname`, as `cnmfrm_c`, or `nothing` if there is none.

The kernel-pool variable `OBJECT_<cname>_FRAME` (or `OBJECT_<id>_FRAME`) has priority.
Without it the answer is the built-in frame `IAU_<CNAME>`, if it exists. `nothing`
takes the place of the `found` flag.

# Examples
```jldoctest
julia> cnmfrm("Mars"), cnmfrm("NO SUCH BODY")
((10014, "IAU_MARS"), nothing)
```

See also [`cidfrm`](@ref).
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
    ccifrm([ctx,] frclss, clssid) -> Union{Nothing,Tuple{Int,String,Int}}

Return `(frcode, frname, center)` of the frame with the class `frclss` and the class
ID `clssid`, as `ccifrm_c`, or `nothing` if there is none.

The built-in frames are searched first. If several frames of the kernel pool match,
the one with the smallest ID is returned; CSPICE returns the first one in the order
of its hash table. `nothing` takes the place of the `found` flag.

# Examples
```jldoctest
julia> ccifrm(2, 499)
(10014, "IAU_MARS", 499)
```

See also [`frinfo`](@ref).
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

"Throws `SpiceError(:badframeclass)` unless `frmcls` is -1 or a class 1 to 6, as `bltfrm_` and `kplfrm_`."
function _check_frame_class(who::AbstractString, frmcls::Integer)
    (frmcls == -1 || 1 <= frmcls <= 6) || throw(SpiceError(:badframeclass,
        "$who: the frame class $frmcls is not supported; it must be -1 (all) or 1 to 6"))
    return nothing
end

"The built-in frames of `zzfdat_` as `(ID, class)`, in the order of the IDs."
const _BUILTIN_FRAME_CLASSES = sort!(vcat(
    [(i, 1) for i in eachindex(INERTIAL_FRAME_DEFS)],
    [(Int(f.id), Int(f.class)) for f in BUILTIN_NONINERTIAL_FRAMES]))

"""
    bltfrm(frmcls) -> Vector{Int}

Return the IDs of the built-in frames of the class `frmcls`, in increasing order, as
`bltfrm_c` of CSPICE.

`frmcls` is a frame class of [`frinfo`](@ref) (1 to 6), or -1 for every class. The
frames that the kernel pool defines are not included; [`kplfrm`](@ref) gives those.

Throws `SpiceError(:badframeclass)` if `frmcls` is not -1 or a class 1 to 6.

# Examples
```jldoctest
julia> length(bltfrm(1)), bltfrm(1)[1:3], bltfrm(3)
(21, [1, 2, 3], Int64[])
```
"""
function bltfrm(frmcls::Integer)
    _check_frame_class("bltfrm", frmcls)
    return Int[id for (id, class) in _BUILTIN_FRAME_CLASSES if frmcls == -1 || class == frmcls]
end

"""
    kplfrm([ctx,] frmcls) -> Vector{Int}

Return the IDs of the frames of the class `frmcls` that the kernel pool defines, in
increasing order and without duplicates, as `kplfrm_c` of CSPICE.

`frmcls` is a frame class of [`frinfo`](@ref) (1 to 6), or -1 for every class. As in
CSPICE, a frame is listed if the pool has

- `FRAME_<name> = <id>`, for the value `<name>` of a `FRAME_*_NAME` variable;
- `FRAME_<id>_NAME`;
- `FRAME_<id>_CLASS` or `FRAME_<name>_CLASS`.

The other variables of a frame definition are not checked, so an ID can be listed that
[`frinfo`](@ref) does not know. The built-in frames are not included unless the pool
defines them; [`bltfrm`](@ref) gives those.

Throws `SpiceError(:badframeclass)` if `frmcls` is not -1 or a class 1 to 6.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> lmpool(ctx, ["FRAME_TOPO = 1500001", "FRAME_1500001_NAME = 'TOPO'",
                    "FRAME_1500001_CLASS = 4", "FRAME_1500001_CLASS_ID = 1500001",
                    "FRAME_1500001_CENTER = 399"]);

julia> kplfrm(ctx, 4), kplfrm(ctx, 3)
([1500001], Int64[])
```
"""
function kplfrm(ks::KernelSet, frmcls::Integer)
    _check_frame_class("kplfrm", frmcls)
    pool = ks.pool
    # `gipool_` gives the first value of a numeric variable, rounded; names are 32
    # characters long in CSPICE, so a longer name is cut.
    getint(key) = (v = get(pool, first(key, 32), nothing);
        v isa Vector{Float64} && !isempty(v) ? _nint(v[1]) : nothing)
    ids = Int[]
    for (key, val) in pool
        _pool_x_match(key, "FRAME_*_NAME") || continue
        (val isa Vector{String} && !isempty(val)) || continue
        name = strip(first(val[1], 32), ' ')
        isempty(name) && (name = " ")
        idcode = getint("FRAME_" * name)
        idcode === nothing && continue
        v = get(pool, "FRAME_$(idcode)_NAME", nothing)
        (v isa Vector{String} && !isempty(v)) || continue
        class = getint("FRAME_$(idcode)_CLASS")
        class === nothing && (class = getint("FRAME_$(name)_CLASS"))
        class === nothing && continue
        (frmcls == -1 || class == frmcls) && push!(ids, idcode)
    end
    return unique!(sort!(ids))
end

"""
    tkfram([ctx,] frcode) -> Tuple{SMatrix{3,3,Float64,9},Int}

Return `(rot, relative)` of the fixed-offset (class 4) frame with the class ID
`frcode`, as `tkfram_c`. `frcode` is the value of `FRAME_<id>_CLASS_ID`.

`rot` rotates a vector from the frame to its relative frame, the frame with ID
`relative`; this is the meaning of `TKFRAME_<frcode>_MATRIX`. The definition comes
from the `TKFRAME_<frcode>_*` variables of the kernel pool, or from the
`TKFRAME_<name>_*` variables, where `<name>` is the name of the frame with the ID
`frcode`, as in CSPICE. A TK frame whose class ID is not its frame ID therefore needs
another frame with that ID, in CSPICE as here.

Throws `SpiceError`:
- `:zeroframeid` if `frcode == 0`;
- `:incompleteframe` if no frame has the ID `frcode`;
- `:competingframespec` if a variable is present in both forms;
- `:variablenotfound`, `:badvariablesize` or `:badvariabletype` if a variable that the
  definition needs is missing, has the wrong number of values or the wrong type;
- `:badframespec` if the relative frame is not known, `:badframespec2` if it is the
  frame itself;
- `:unknownframespec` if `TKFRAME_<frcode>_SPEC` is not `MATRIX`, `ANGLES` or
  `QUATERNION`;
- `:unitsnotrec` or `:incompatibleunits` if `TKFRAME_<frcode>_UNITS` is not an angle
  unit, and `:badaxisnumbers` if an axis is not 1, 2 or 3.

[`pxform`](@ref) and the other calls that pass through such a frame throw the same.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> lmpool(ctx, ["FRAME_TOPO = -1000", "FRAME_-1000_NAME = 'TOPO'",
                    "FRAME_-1000_CLASS = 4", "FRAME_-1000_CLASS_ID = -1000",
                    "FRAME_-1000_CENTER = 399", "TKFRAME_-1000_RELATIVE = 'J2000'",
                    "TKFRAME_-1000_SPEC = 'ANGLES'", "TKFRAME_-1000_ANGLES = (0, 0, 90)",
                    "TKFRAME_-1000_AXES = (1, 2, 3)", "TKFRAME_-1000_UNITS = 'DEGREES'"])

julia> rot, relative = tkfram(ctx, -1000);

julia> relative, round.(rot; digits=12) == [0 1 0; -1 0 0; 0 0 1]
(1, true)

julia> pdpool(ctx, "TKFRAME_-1000_AXES", [1.0, 2.0]);

julia> try tkfram(ctx, -1000) catch err; err.kind end
:badvariablesize
```
"""
function tkfram(ks::KernelSet, frcode::Integer)
    r = frcode == 0 ? nothing : ccifrm(ks, 4, frcode)
    if r !== nothing
        tk = get(ks.tk, Int32(r[1]), nothing)
        tk === nothing || return tk.rot, Int(tk.relative)
    end
    tk = _tk_parse(ks.pool, ks.frames, frcode)    # throws the error of `tkfram_`
    return tk.rot, Int(tk.relative)
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
    return _m_mxm(r, c), _m_mxm(dr, c)
end

"""
    tipbod([ctx,] ref, body, et) -> SMatrix{3,3}

Return the matrix that rotates a position vector from the inertial frame `ref` to the
body-fixed frame of `body` at `et`, as `tipbod_c`.

`ref` is the name or the ID of a built-in inertial frame. `body` is a NAIF body ID, or
the frame class ID of a binary PCK frame, the ID that its segments carry (31008 for
`MOON_PA_DE440`). `et` is in seconds past J2000 TDB. A binary PCK has priority over the
constants of a text PCK.

Throws `SpiceError(:irfnotrec)` if `ref` is not a built-in inertial frame,
`SpiceError(:badarg)` if `et` is not finite, and `SpiceError(:framedatanotfound)` if
no PCK data are loaded for the body.

See also [`tisbod`](@ref), [`pxform`](@ref).
"""
tipbod(ks::KernelSet, ref::NameOrID, body::Integer, et::Real) = _tis_jet(ks, ref, body, et)[1]

"""
    tisbod([ctx,] ref, body, et) -> SMatrix{6,6}

Return the matrix that transforms a state from the inertial frame `ref` to the
body-fixed frame of `body` at `et`, as `tisbod_c`.

The matrix is `[R 0; dR R]`, with `R` the matrix of [`tipbod`](@ref) and `dR` its
time derivative in 1/s. The arguments and the errors are those of `tipbod`.

See also [`sxform`](@ref).
"""
function tisbod(ks::KernelSet, ref::NameOrID, body::Integer, et::Real)
    r, dr = _tis_jet(ks, ref, body, et)
    return vcat(hcat(r, zero(r)), hcat(dr, r))
end

"""
    pckfrm!(ids, pck) -> ids

Insert the frame class IDs of the frames that have data in the binary PCK file `pck`
into the vector `ids`, as `pckfrm_c` does with a SPICE cell.

`ids` is used as a set: on return it is sorted and without duplicates. The vector
comes first, as in SPICE.jl; CSPICE has the cell last. The file is read from disk and
need not be loaded.

Throws the errors of [`pckfrm`](@ref).
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

Return the frame class IDs of the frames that have data in the binary PCK file `pck`,
sorted and without duplicates, as `pckfrm_c`.

The file is read from disk and need not be loaded.

Throws `SpiceError(:filenotfound)` if the file does not exist, `SpiceError(:invalidformat)`
if it is a transfer file (a binary kernel in text form), `SpiceError(:invalidarchtype)` if
it is not a DAF, and `SpiceError(:invalidfiletype)` if it is a DAF of another kernel type.

See also [`pckcov`](@ref).
"""
pckfrm(pck::AbstractString) = pckfrm!(Int[], pck)

"""
    pckcov!(cover, pck, idcode) -> cover

Insert the coverage of the frame with the class ID `idcode` in the binary PCK file
`pck` into the window `cover`, as `pckcov_c`.

The result is the union of the old content of `cover` with the new intervals. The
window comes first, as in SPICE.jl; CSPICE has the cell last. See [`pckcov`](@ref) for
the intervals and the errors.
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

Return the coverage [`Window`](@ref) of the frame with the class ID `idcode` in the
binary PCK file `pck`, as `pckcov_c`.

The window is the union of the intervals of the segments of that frame, in seconds
past J2000 TDB. It is empty if the frame has no data. The file is read from disk and
need not be loaded.

Throws the errors of [`pckfrm`](@ref).

See also [`pckcov!`](@ref), [`spkcov`](@ref).
"""
pckcov(pck::AbstractString, idcode::Integer) = pckcov!(Window(), pck, idcode)

"""
    bodfnd([ctx,] body, item) -> Bool

Return `true` if the kernel pool holds the variable `BODY<body>_<item>`, as
`bodfnd_c`.

`body` is a NAIF ID code. `item` is case sensitive; blanks around it are ignored. The
variable can be numeric or character.

See also [`bodvcd`](@ref).
"""
bodfnd(ks::KernelSet, body::Integer, item::AbstractString) =
    haskey(ks.pool, string("BODY", body, "_", strip(item)))

"""
    bodvcd([ctx,] bodyid, item) -> Vector{Float64}

Return the values of the kernel-pool variable `BODY<bodyid>_<item>`, as `bodvcd_c`.

`bodyid` is a NAIF ID code; [`bodvrd`](@ref) also takes a body name. `item` is case
sensitive. The result is a new vector with every value, in the units of the kernel;
there is no `maxn` argument.

Throws `SpiceError(:kernelvarnotfound)` if the variable is not in the pool and
`SpiceError(:typemismatch)` if it holds character values.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> pdpool(ctx, "BODY399_RADII", [6378.1366, 6378.1366, 6356.7519]);

julia> bodvcd(ctx, 399, "RADII")
3-element Vector{Float64}:
 6378.1366
 6378.1366
 6356.7519
```
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
kplfrm(ctx::SpiceContext, frmcls::Integer) = kplfrm(kernelset(ctx), frmcls)
kplfrm(frmcls::Integer) = kplfrm(default_context(), frmcls)
