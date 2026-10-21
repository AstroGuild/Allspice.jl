# Surface names and codes (`srfc2s`, `srfs2c`, `srfcss`, `srfscc`), the computation
# method strings, and the surface point calls `latsrf` and `srfnrm`.
#
# The surface name map is defined by the kernel variables NAIF_SURFACE_NAME,
# NAIF_SURFACE_CODE and NAIF_SURFACE_BODY. It is read from the `KernelSet` snapshot at
# each query, so there is no hidden state and no cache.

const _DSK_SRF_VARS = ("NAIF_SURFACE_NAME", "NAIF_SURFACE_CODE", "NAIF_SURFACE_BODY")
const _DSK_MXNSRF = 2000                  # largest number of surface name assignments
const _DSK_MAXSRF = 100                   # largest number of surfaces in a method string

"""
    _dsk_surface_vars(ks) -> Union{Nothing,Tuple{Vector{String},Vector{Float64},Vector{Float64}}}

The three kernel variables of the surface name map after the checks of `zzsrfker`, or
`nothing` if none of them is in the pool. Throws `SpiceError` (`:badsurfacemap` if only
some are present, `:badvariabletype`, `:toomanysurfaces`, `:arraysizemismatch`,
`:blanknameassigned`).
"""
function _dsk_surface_vars(ks::KernelSet)
    names = get(ks.pool, _DSK_SRF_VARS[1], nothing)
    codes = get(ks.pool, _DSK_SRF_VARS[2], nothing)
    bodies = get(ks.pool, _DSK_SRF_VARS[3], nothing)
    (names === nothing && codes === nothing && bodies === nothing) && return nothing
    (names === nothing || codes === nothing || bodies === nothing) && throw(SpiceError(
        :badsurfacemap, "the surface name map needs the three kernel variables " *
        "NAIF_SURFACE_NAME, NAIF_SURFACE_CODE and NAIF_SURFACE_BODY, but only some are present"))
    (names isa Vector{String} && codes isa Vector{Float64} && bodies isa Vector{Float64}) ||
        throw(SpiceError(:badvariabletype, "NAIF_SURFACE_NAME must be character, and " *
            "NAIF_SURFACE_CODE and NAIF_SURFACE_BODY must be numeric"))
    max(length(names), length(codes), length(bodies)) > _DSK_MXNSRF && throw(SpiceError(
        :toomanysurfaces, "a surface name map variable has more than $_DSK_MXNSRF values"))
    (length(codes) == length(names) && length(bodies) == length(names)) || throw(SpiceError(
        :arraysizemismatch, "the surface name map variables have the sizes " *
        "$(length(names)), $(length(codes)) and $(length(bodies)); they must be equal"))
    for n in names
        isempty(strip(n, ' ')) && throw(SpiceError(:blanknameassigned,
            "a blank surface name is assigned in NAIF_SURFACE_NAME"))
    end
    return names, codes, bodies
end

"""
    _dsk_surface_n2c(ks, name, bodyid) -> Union{Int,Nothing}

`zzsrfn2c`: the code of the surface `name` of the body `bodyid`. Case and extra blanks
do not matter, and the last assignment of a name for a body has priority.
"""
function _dsk_surface_n2c(ks::KernelSet, name::AbstractString, bodyid::Integer)
    vars = _dsk_surface_vars(ks)
    vars === nothing && return nothing
    names, codes, bodies = vars
    nn = normalize_name(name)
    for i in length(names):-1:1          # invariant: no assignment after i matches
        (round(Int, bodies[i]) == bodyid && normalize_name(names[i]) == nn) &&
            return round(Int, codes[i])
    end
    return nothing
end

"""
    _dsk_surface_c2n(ks, code, bodyid) -> Union{String,Nothing}

`zzsrfc2n`: the name of the surface `code` of the body `bodyid`. As in `zzsrfini`, an
assignment whose name and body occur again later in the variables is masked and gives
no name for its code; among the others the last one for the code and the body has
priority.
"""
function _dsk_surface_c2n(ks::KernelSet, code::Integer, bodyid::Integer)
    vars = _dsk_surface_vars(ks)
    vars === nothing && return nothing
    names, codes, bodies = vars
    n = length(names)
    for i in n:-1:1                      # invariant: no unmasked assignment after i matches
        (round(Int, bodies[i]) == bodyid && round(Int, codes[i]) == code) || continue
        nn = normalize_name(names[i])
        masked = false
        for j in i+1:n
            if round(Int, bodies[j]) == bodyid && normalize_name(names[j]) == nn
                masked = true
                break
            end
        end
        masked || return String(rstrip(names[i], ' '))
    end
    return nothing
end

"Throws `SpiceError(:emptystring)` if the string argument `name` of `who` has no character, as the C interface of CSPICE."
function _dsk_nonempty(who::AbstractString, name::AbstractString, s::AbstractString)
    isempty(s) && throw(SpiceError(:emptystring, "$who: the string \"$name\" has length zero"))
    return nothing
end

"`nparsi`: the integer that the string holds, or `nothing` if it does not hold a number."
function _dsk_parse_int(s::AbstractString)
    t = strip(s, ' ')
    isempty(t) && return nothing
    v = tryparse(Int, t)
    v === nothing || return (typemin(Int32) <= v <= typemax(Int32) ? v : nothing)
    all(c -> isdigit(c) || c in "+-.eEdD", t) || return nothing
    x = tryparse(Float64, replace(t, r"[dD]" => "e"))
    (x === nothing || !isfinite(x)) && return nothing
    x = trunc(x)                         # `nparsi` drops the fraction
    (typemin(Int32) <= x <= typemax(Int32)) || return nothing
    return Int(x)
end

"""
    srfc2s([ctx | ks,] code, bodyid) -> (srfstr, isname)

Return the name of the surface `code` of the body `bodyid`, as `srfc2s_c`.

`isname` is `true` if the kernel pool assigns a name to the surface, with the variables
`NAIF_SURFACE_NAME`, `NAIF_SURFACE_CODE` and `NAIF_SURFACE_BODY`. Otherwise `srfstr` is
the code as a decimal string and `isname` is `false`.

Throws `SpiceError` if the three kernel variables are not consistent:
- `:badsurfacemap` if only some of the three are present;
- `:badvariabletype` if the names are not strings or the other two are not numbers;
- `:toomanysurfaces` if a variable has more than 2000 values;
- `:arraysizemismatch` if the three do not have the same number of values;
- `:blanknameassigned` if a name is blank.

See also [`srfcss`](@ref), [`srfscc`](@ref), [`srfs2c`](@ref).
"""
function srfc2s(ks::KernelSet, code::Integer, bodyid::Integer)
    name = _dsk_surface_c2n(ks, code, bodyid)
    return name === nothing ? (string(code), false) : (name, true)
end

"""
    srfcss([ctx | ks,] code, bodstr) -> (srfstr, isname)

Return the name of the surface `code` of the body `bodstr`, as `srfcss_c`. `bodstr` is a
body name or a string that holds a body ID; the result is that of [`srfc2s`](@ref).

A body that is not known gives the code as a decimal string and `isname = false`.

Throws `SpiceError(:emptystring)` if `bodstr` is empty, and the kinds of
[`srfc2s`](@ref).
"""
function srfcss(ks::KernelSet, code::Integer, bodstr::AbstractString)
    _dsk_nonempty("srfcss", "bodstr", bodstr)
    bodyid = bods2c(ks, bodstr)
    bodyid === nothing && return string(code), false
    return srfc2s(ks, code, bodyid)
end

"""
    srfscc([ctx | ks,] srfstr, bodyid) -> Union{Int,Nothing}

Return the code of the surface `srfstr` of the body `bodyid`, as `srfscc_c`, or `nothing`
if the string is not a surface name and does not hold an integer.

A name that the kernel pool assigns has priority; case and extra blanks do not matter.
Otherwise the result is the integer that the string holds. `nothing` replaces the `found`
flag of CSPICE.

Throws `SpiceError(:emptystring)` if `srfstr` is empty, and the kinds of
[`srfc2s`](@ref).
"""
function srfscc(ks::KernelSet, srfstr::AbstractString, bodyid::Integer)
    _dsk_nonempty("srfscc", "srfstr", srfstr)
    code = _dsk_surface_n2c(ks, srfstr, bodyid)
    return code === nothing ? _dsk_parse_int(srfstr) : code
end

"""
    srfs2c([ctx | ks,] srfstr, bodstr) -> Union{Int,Nothing}

Return the code of the surface `srfstr` of the body `bodstr`, as `srfs2c_c`. `bodstr` is
a body name or a string that holds a body ID; the result is that of [`srfscc`](@ref).

A body that is not known gives `nothing`.

Throws `SpiceError(:emptystring)` if `srfstr` or `bodstr` is empty, and the kinds of
[`srfc2s`](@ref).
"""
function srfs2c(ks::KernelSet, srfstr::AbstractString, bodstr::AbstractString)
    _dsk_nonempty("srfs2c", "srfstr", srfstr)
    _dsk_nonempty("srfs2c", "bodstr", bodstr)
    bodyid = bods2c(ks, bodstr)
    bodyid === nothing && return nothing
    return srfscc(ks, srfstr, bodyid)
end

srfc2s(ctx::SpiceContext, code::Integer, bodyid::Integer) = srfc2s(kernelset(ctx), code, bodyid)
srfc2s(code::Integer, bodyid::Integer) = srfc2s(default_context(), code, bodyid)
srfcss(ctx::SpiceContext, code::Integer, bodstr::AbstractString) = srfcss(kernelset(ctx), code, bodstr)
srfcss(code::Integer, bodstr::AbstractString) = srfcss(default_context(), code, bodstr)
srfscc(ctx::SpiceContext, srfstr::AbstractString, bodyid::Integer) = srfscc(kernelset(ctx), srfstr, bodyid)
srfscc(srfstr::AbstractString, bodyid::Integer) = srfscc(default_context(), srfstr, bodyid)
srfs2c(ctx::SpiceContext, srfstr::AbstractString, bodstr::AbstractString) =
    srfs2c(kernelset(ctx), srfstr, bodstr)
srfs2c(srfstr::AbstractString, bodstr::AbstractString) = srfs2c(default_context(), srfstr, bodstr)

# ---------------------------------------------------------------------------------------
# Computation method strings
# ---------------------------------------------------------------------------------------

"""
    _dsk_lex_method(method) -> Vector{String}

`zzlexmet`: the tokens of a method string. A token is one of the delimiters `/ , = :`,
a quoted string (with its quotes), or a run of other non-blank characters. Throws
`SpiceError(:syntaxerror)` for a quote that is not closed.
"""
function _dsk_lex_method(method::AbstractString)
    chars = collect(method)
    n = length(chars)
    tokens = String[]
    i = 1
    while i <= n                         # invariant: the tokens of chars[1:i-1] are in `tokens`
        c = chars[i]
        if c == ' '
            i += 1
        elseif c == '"'
            j = i + 1
            closed = false
            while j <= n
                if chars[j] == '"'
                    if j < n && chars[j+1] == '"'
                        j += 2
                        continue
                    end
                    closed = true
                    break
                end
                j += 1
            end
            closed || throw(SpiceError(:syntaxerror,
                "a quoted string in the method string is not terminated: $method"))
            push!(tokens, String(chars[i:j]))
            i = j + 1
        elseif c in "/,=:"
            push!(tokens, string(c))
            i += 1
        else
            j = i
            while j < n && !(chars[j+1] in " /,=:\"")
                j += 1
            end
            push!(tokens, String(chars[i:j]))
            i = j + 1
        end
    end
    length(tokens) > _DSK_MAXSRF && throw(SpiceError(:arraytoosmall,
        "the method string has more than $_DSK_MAXSRF tokens"))
    return tokens
end

@noinline _dsk_bad_method(method, why) = throw(SpiceError(:badmethodsyntax,
    "$why in the method string <$method>"))

"True if `a` and `b` are equal when case and blanks are ignored (`eqstr`)."
_dsk_eqstr(a::AbstractString, b::AbstractString) =
    uppercase(replace(a, ' ' => "")) == uppercase(replace(b, ' ' => ""))

"""
    _dsk_parse_method(ks, bodyid, method) -> (shape, subtyp, pri, srflst)

`zzprsmet`: the shape (`:ellipsoid` or `:dsk`), the sub-type (`""` if there is none),
the priority flag and the surface ID list of a method string such as `"ELLIPSOID"` or
`"DSK/UNPRIORITIZED/SURFACES = 1, \\"name\\""`. Throws `SpiceError`
(`:badmethodsyntax`, `:badpriorityspec` for a DSK method without `UNPRIORITIZED`,
`:idcodenotfound` for an unknown surface, `:toomanysurfaces`).
"""
function _dsk_parse_method(ks::KernelSet, bodyid::Integer, method::AbstractString)
    tok = pushfirst!(_dsk_lex_method(normalize_name(method)), "/")
    n = length(tok)
    nshape = count(t -> t == "ELLIPSOID" || t == "DSK", tok)
    nshape == 1 || _dsk_bad_method(method, "the shape (ELLIPSOID or DSK) must occur one time")
    subtyp = ""
    srflst = Int[]
    nsubtp = nlmbtp = ntrmtp = 0
    if "DSK" in tok
        pri = true
        npri = nsrfls = 0
        i = 1
        while i < n                      # invariant: the tokens before i are processed
            if tok[i] == "/"
                i += 1
                t = tok[i]
                if t == "SURFACES"
                    nsrfls += 1
                    nsrfls > 1 && _dsk_bad_method(method, "extra SURFACES keyword")
                    i > n - 2 && _dsk_bad_method(method, "the surface list is too short")
                    i += 1
                    tok[i] == "=" || _dsk_bad_method(method, "'=' is missing after SURFACES")
                    i += 1
                    while true
                        s = tok[i]
                        if first(s) == '"'
                            length(s) > 2 || _dsk_bad_method(method, "an empty surface name")
                            s = s[nextind(s, 1):prevind(s, lastindex(s))]
                        end
                        code = srfscc(ks, s, bodyid)
                        code === nothing && throw(SpiceError(:idcodenotfound, "could not " *
                            "map the surface name $s to an ID code for the body $bodyid"))
                        length(srflst) == _DSK_MAXSRF && throw(SpiceError(:toomanysurfaces,
                            "the method string names more than $_DSK_MAXSRF surfaces"))
                        push!(srflst, code)
                        i == n && break
                        i += 1
                        if tok[i] == "/"
                            i -= 1
                            break
                        elseif tok[i] == ","
                            i < n || _dsk_bad_method(method, "a surface is missing after ','")
                            i += 1
                        else
                            _dsk_bad_method(method, "unexpected token <$(tok[i])>")
                        end
                    end
                elseif t == "UNPRIORITIZED"
                    npri += 1
                    npri > 1 && _dsk_bad_method(method, "extra UNPRIORITIZED keyword")
                    pri = false
                elseif t == "INTERCEPT" || t == "NADIR"
                    nsubtp += 1
                    nsubtp > 1 && _dsk_bad_method(method, "extra sub-observer point type")
                    subtyp = t
                elseif t == "TANGENT" || t == "GUIDED"
                    nlmbtp += 1
                    nlmbtp > 1 && _dsk_bad_method(method, "extra limb type")
                elseif t == "UMBRAL" || t == "PENUMBRAL"
                    ntrmtp += 1
                    ntrmtp > 1 && _dsk_bad_method(method, "extra terminator type")
                elseif t == "/"
                    _dsk_bad_method(method, "an extra '/'")
                elseif t != "DSK"
                    _dsk_bad_method(method, "unexpected token <$t>")
                end
            end
            i += 1
        end
        pri && throw(SpiceError(:badpriorityspec, "the method string <$method> does not " *
            "contain UNPRIORITIZED; prioritized DSK data are not supported, as in CSPICE"))
        return :dsk, subtyp, false, srflst
    end
    _dsk_eqstr(method, "ELLIPSOID") && return :ellipsoid, "", true, srflst
    for (a, b) in (("NEARPOINT", "NEAR POINT"), ("INTERCEPT", "INTERCEPT")), sep in (":", "/")
        (_dsk_eqstr(method, a * sep * "ELLIPSOID") || _dsk_eqstr(method, "ELLIPSOID" * sep * a)) &&
            return :ellipsoid, b, true, srflst
    end
    i = 1
    while i <= n                         # invariant: the tokens before i are processed
        tok[i] == "/" || _dsk_bad_method(method, "unexpected token <$(tok[i])>")
        i == n && _dsk_bad_method(method, "a keyword is missing after '/'")
        i += 1
        t = tok[i]
        if t == "ELLIPSOID"
        elseif t == "INTERCEPT" || t == "NADIR" || t == "NEAR"
            nsubtp += 1
            nsubtp > 1 && _dsk_bad_method(method, "extra sub-observer point type")
            if t == "NEAR"
                (i < n && tok[i+1] == "POINT") && (subtyp = "NADIR")
            else
                subtyp = t
            end
        elseif t == "TANGENT" || t == "GUIDED"
            nlmbtp += 1
            nlmbtp > 1 && _dsk_bad_method(method, "extra limb type")
        elseif t == "UMBRAL" || t == "PENUMBRAL"
            ntrmtp += 1
            ntrmtp > 1 && _dsk_bad_method(method, "extra terminator type")
        else
            _dsk_bad_method(method, "unexpected token <$t>")
        end
        i += 1
    end
    return :ellipsoid, subtyp, true, srflst
end

# ---------------------------------------------------------------------------------------
# Surface points from longitude and latitude; normals at surface points
# ---------------------------------------------------------------------------------------

"Target code, frame ID, shape and surface list for `latsrf` and `srfnrm`, after their common checks."
function _dsk_surface_setup(ks::KernelSet, method::AbstractString, target::NameOrID,
        fixref::NameOrID)
    _dsk_nonempty("latsrf/srfnrm", "method", method)
    target isa AbstractString && _dsk_nonempty("latsrf/srfnrm", "target", target)
    fixref isa AbstractString && _dsk_nonempty("latsrf/srfnrm", "fixref", fixref)
    trgcde = _dsk_target(ks, target)
    fixfid = _dsk_frame_code(ks, fixref)
    center = fixfid == 0 ? nothing : _dsk_frame_center(ks, fixfid)
    center === nothing && throw(SpiceError(:noframe, "the reference frame $fixref is not recognized"))
    center == trgcde || throw(SpiceError(:invalidframe, "the reference frame $fixref is " *
        "not centered at the target body $target"))
    shape, subtyp, _, srflst = _dsk_parse_method(ks, trgcde, method)
    isempty(subtyp) || throw(SpiceError(:invalidmethod, "the method string <$method> " *
        "contains the sub-observer point type $subtyp, which is not valid here"))
    return trgcde, fixfid, shape, srflst
end

"The three radii of the body `code` from the kernel pool; throws `SpiceError(:badaxislength)` unless all are positive."
function _dsk_radii(ks::KernelSet, code::Integer)
    radii = bodvrd(ks, code, "RADII")
    length(radii) == 3 || throw(SpiceError(:badaxislength,
        "BODY$(code)_RADII has $(length(radii)) values; 3 are required"))
    min(min(radii[1], radii[2]), radii[3]) <= 0.0 && throw(SpiceError(:badaxislength,
        "the radii of the body $code must be positive but were $radii"))
    return radii[1], radii[2], radii[3]
end

@inline _dsk_get2(a::AbstractMatrix, i::Int) = (Float64(a[1, i]), Float64(a[2, i]))
@inline function _dsk_get2(a::AbstractVector, i::Int)
    e = a[begin+i-1]
    length(e) == 2 || throw(SpiceError(:badarg, "element $i does not have 2 components"))
    return Float64(e[begin]), Float64(e[begin+1])
end
function _dsk_count2(a::AbstractMatrix)
    size(a, 1) == 2 || throw(SpiceError(:badarg, "a 2 x n matrix is required, got $(size(a))"))
    return size(a, 2)
end
_dsk_count2(a::AbstractVector) = length(a)

"""
    latsrf([ctx,] method, target, et, fixref, lonlat) -> Vector{SVector{3,Float64}}

Return the surface points of `target` (km, in the body-fixed frame `fixref`) at the
planetocentric longitudes and latitudes `lonlat`, as `latsrf_c`.

`lonlat` holds the pairs (longitude, latitude) in radians: a 2 x n matrix, or a vector of
n pairs. `target` and `fixref` are names or IDs; the frame must be centered on the target.
`method` selects the surface:

- `"ELLIPSOID"`: the points are on the ellipsoid of the radii `BODY<id>_RADII` of the
  kernel pool. `et` is not used.
- `"DSK/UNPRIORITIZED[/SURFACES = <list>]"`: the points are on the loaded DSK segments
  at the epoch `et` (seconds past J2000 TDB). A ray starts outside of the body in the
  given direction and goes to the center; the point is its first intercept. The list
  names surfaces by name or by ID; without it every surface of the body is used.

A DSK method uses the DSK tolerances of the context ([`dskstl`](@ref)).

Throws `SpiceError`:
- `:emptystring` if `method`, `target` or `fixref` is an empty string;
- `:idcodenotfound` if the target or a surface of the list is not known;
- `:noframe` if the frame is not known, `:invalidframe` if it is not centered on the
  target;
- `:badmethodsyntax` if `method` is not valid, `:badpriorityspec` if a DSK method does not
  have `UNPRIORITIZED`, `:invalidmethod` if `method` has a sub-observer point type;
- `:kernelvarnotfound` if the radii are not loaded, `:badaxislength` if they are not three
  positive values (ellipsoid);
- `:dskdatanotfound` if no DSK segment of the target with a listed surface is loaded,
  `:pointnotfound` if a ray meets no surface, `:shapenotsupported` if an intercept is on
  the far side of the center (DSK).

See also [`srfnrm`](@ref), [`dskxsi`](@ref).
"""
function latsrf(ks::KernelSet, method::AbstractString, target::NameOrID, et::Real,
        fixref::NameOrID, lonlat)
    return _dsk_latsrf(ks, _DSK_DEFAULT_TOLS, method, target, et, fixref, lonlat)
end
latsrf(ctx::SpiceContext, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    lonlat) = _dsk_latsrf(kernelset(ctx), dsktolerances(ctx), method, target, et, fixref, lonlat)
latsrf(method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID, lonlat) =
    latsrf(default_context(), method, target, et, fixref, lonlat)

function _dsk_latsrf(ks::KernelSet, tols::DSKTolerances, method::AbstractString, target, et::Real,
        fixref, lonlat)
    trgcde, fixfid, shape, srflst = _dsk_surface_setup(ks, method, target, fixref)
    npts = _dsk_count2(lonlat)
    out = Vector{Vec3{Float64}}(undef, npts)
    if shape === :ellipsoid
        a, b, c = _dsk_radii(ks, trgcde)
        for i in 1:npts                  # invariant: out[1:i-1] holds the points 1 .. i-1
            lon, lat = _dsk_get2(lonlat, i)
            p = _dsk_latrec(1.0, lon, lat)
            d1 = p[1] / a
            d2 = p[2] / b
            d3 = p[3] / c
            level = d1 * d1 + d2 * d2 + d3 * d3
            level <= 0.0 && throw(SpiceError(:pointtoosmall, "the direction vector is too small"))
            sq = sqrt(level)
            out[i] = Vec3{Float64}(p[1] / sq, p[2] / sq, p[3] / sq)
        end
        return out
    end
    segs = get(ks.dsk, trgcde, nothing)
    segs === nothing && _dsk_no_segments(trgcde, srflst)
    r = _dsk_sph(ks, segs, trgcde, srflst)[2] * 2
    t = Float64(et)
    for i in 1:npts                      # invariant: out[1:i-1] holds the points 1 .. i-1
        lon, lat = _dsk_get2(lonlat, i)
        x = _dsk_latrec(r, lon, lat)
        xpt, _, _, found = _dsk_bux(ks, tols, trgcde, srflst, t, fixfid, x, _dsk_vminus(x))
        found || throw(SpiceError(:pointnotfound, "no surface point was found on body " *
            "$trgcde at planetocentric longitude $(rad2deg(lon)) deg, latitude $(rad2deg(lat)) deg"))
        _dsk_vdot(x, xpt) < 0.0 && throw(SpiceError(:shapenotsupported, "the intercept is " *
            "on the far side of the center of body $trgcde; such a shape is not supported"))
        out[i] = xpt
    end
    return out
end

"""
    srfnrm([ctx,] method, target, et, fixref, srfpts) -> Vector{SVector{3,Float64}}

Return the outward unit normals of `target` at the surface points `srfpts`, in the
body-fixed frame `fixref`, as `srfnrm_c`.

`srfpts` holds the points (km, in the frame `fixref`): a 3 x n matrix, or a vector of n
3-vectors. `method`, `target`, `et` and `fixref` are as for [`latsrf`](@ref). For a DSK
method the normal is that of the plate on which the point lies.

Throws `SpiceError`:
- the kinds of [`latsrf`](@ref) for `method`, `target` and `fixref`;
- `:kernelvarnotfound` if the radii are not loaded, `:badaxislength` if they are not three
  positive values, `:pointnotonsurface` if a
  point is not on the ellipsoid within the point-membership margin (item `PTMEMM` of
  [`dskgtl`](@ref)) (ellipsoid);
- `:nodsksegments` if no DSK segment of the target is loaded, `:pointnotinsegment` if no
  segment holds a point, `:pointoffsurface` if a point is not on a plate (DSK).
"""
function srfnrm(ks::KernelSet, method::AbstractString, target::NameOrID, et::Real,
        fixref::NameOrID, srfpts)
    return _dsk_srfnrm(ks, _DSK_DEFAULT_TOLS, method, target, et, fixref, srfpts)
end
srfnrm(ctx::SpiceContext, method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID,
    srfpts) = _dsk_srfnrm(kernelset(ctx), dsktolerances(ctx), method, target, et, fixref, srfpts)
srfnrm(method::AbstractString, target::NameOrID, et::Real, fixref::NameOrID, srfpts) =
    srfnrm(default_context(), method, target, et, fixref, srfpts)

function _dsk_srfnrm(ks::KernelSet, tols::DSKTolerances, method::AbstractString, target, et::Real,
        fixref, srfpts)
    trgcde, fixfid, shape, srflst = _dsk_surface_setup(ks, method, target, fixref)
    npts = _dsk_count3(srfpts)
    out = Vector{Vec3{Float64}}(undef, npts)
    if shape === :ellipsoid
        a, b, c = _dsk_radii(ks, trgcde)
        limit = tols.ptmemm * max(max(a, b), c)
        a2 = a * a
        b2 = b * b
        c2 = c * c
        for i in 1:npts                  # invariant: out[1:i-1] holds the normals 1 .. i-1
            p = _dsk_get3(srfpts, i)
            level = (p[1] * p[1] / a2 + p[2] * p[2] / b2 + p[3] * p[3] / c2)^0.5
            abs(level - 1.0) >= limit && throw(SpiceError(:pointnotonsurface,
                "the input point at index $i is not on the target body's reference ellipsoid"))
            out[i] = _dsk_surfnm(a, b, c, p)
        end
        return out
    end
    t = Float64(et)
    for i in 1:npts                      # invariant: out[1:i-1] holds the normals 1 .. i-1
        out[i] = _dsk_vhat(_dsk_bun(ks, tols, trgcde, srflst, t, fixfid, _dsk_get3(srfpts, i)))
    end
    return out
end
