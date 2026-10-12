# Reference frame definitions: the built-in inertial frames of `chgirf.c`, the
# built-in non-inertial frames of `zzfdat.c`, and the frames that a frame kernel
# defines in the kernel pool (`framex.c`).

# ---------------------------------------------------------------------------------
# Built-in inertial frames (CHGIRF)
# ---------------------------------------------------------------------------------

"""
    INERTIAL_FRAME_DEFS

The inertial frames of CSPICE `chgirf.c` in ID order (1..21) as
(name, base frame, definition). The definition is a list of `angle axis` pairs with
the angle in arcseconds; the rotation from the base frame to the frame is
`[angle1]_axis1 * [angle2]_axis2 * ...`.
"""
const INERTIAL_FRAME_DEFS = (
    ("J2000", "J2000", "0.0  1"),
    ("B1950", "J2000", "1152.84248596724 3  -1002.26108439117  2  1153.04066200330  3"),
    ("FK4", "B1950", "0.525  3"),
    ("DE-118", "B1950", "0.53155  3"),
    ("DE-96", "B1950", "0.4107  3"),
    ("DE-102", "B1950", "0.1359  3"),
    ("DE-108", "B1950", "0.4775  3"),
    ("DE-111", "B1950", "0.5880  3"),
    ("DE-114", "B1950", "0.5529  3"),
    ("DE-122", "B1950", "0.5316  3"),
    ("DE-125", "B1950", "0.5754  3"),
    ("DE-130", "B1950", "0.5247  3"),
    ("GALACTIC", "FK4", "1177200.0  3  225360.0  1  1016100.0  3"),
    ("DE-200", "J2000", "0.0  3"),
    ("DE-202", "J2000", "0.0  3"),
    ("MARSIAU", "J2000", "324000.0D0 3 133610.4D0 2 -152348.4D0 3"),
    ("ECLIPJ2000", "J2000", "84381.448 1"),
    ("ECLIPB1950", "B1950", "84404.836 1"),
    ("DE-140", "J2000", "1152.71013777252 3  -1002.25042010533  2  1153.75719544491  3"),
    ("DE-142", "J2000", "1152.72061453864 3  -1002.25052830351  2  1153.74663857521  3"),
    ("DE-143", "J2000", "1153.03919093833, 3, -1002.24822382286, 2, 1153.42900222357, 3"),
)

const N_INERTIAL_FRAMES = length(INERTIAL_FRAME_DEFS)

"Contract: `a` and `b` are 3x3; returns `a*b` with the summation order of CSPICE `mxm`."
function _mxm(a::Matrix{Float64}, b::Matrix{Float64})
    c = Matrix{Float64}(undef, 3, 3)
    for i in 1:3, j in 1:3
        c[i, j] = a[i, 1] * b[1, j] + a[i, 2] * b[2, j] + a[i, 3] * b[3, j]
    end
    return c
end

"Contract: `axis` is 1, 2 or 3; returns `[angle]_axis * m` with the arithmetic of CSPICE `rotmat`."
function _rotmat(m::Matrix{Float64}, angle::Float64, axis::Int)
    s, c = sin(angle), cos(angle)
    idx = (3, 1, 2, 3, 1)
    t = mod(axis, 3)
    i1, i2, i3 = idx[t+1], idx[t+2], idx[t+3]
    out = Matrix{Float64}(undef, 3, 3)
    for i in 1:3
        out[i1, i] = m[i1, i]
        out[i2, i] = c * m[i2, i] + s * m[i3, i]
        out[i3, i] = -s * m[i2, i] + c * m[i3, i]
    end
    return out
end

"Contract: returns the 21 matrices that rotate vectors from J2000 to each built-in inertial frame, composed from the base-frame chain as the initialization of `irfrot` composes them."
function _build_inertial_rotations()
    dpr = 180.0 / acos(-1.0)
    arcsec = 2.7777777777777778e-4
    trans = Matrix{Float64}[]
    for (i, (_, base, def)) in enumerate(INERTIAL_FRAME_DEFS)
        m = Matrix{Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
        words = split(def, ' '; keepempty=false)
        # Loop invariant: `m` is the product of the rotations of pairs j/2 .. end.
        for j in length(words):-2:2
            axis = parse(Int, rstrip(words[j], ','))
            angle = _nparsd(words[j-1])
            m = _rotmat(m, angle * arcsec / dpr, axis)
        end
        b = findfirst(k -> INERTIAL_FRAME_DEFS[k][1] == base, 1:i)
        b === nothing && error("base frame of inertial frame $i is not defined before it")
        push!(trans, b == i ? _mxm(m, Matrix{Float64}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])) :
                     _mxm(m, trans[b]))
    end
    return [Mat3{Float64}(t) for t in trans]
end

const INERTIAL_ROTATIONS = _build_inertial_rotations()

"""
    inertial_rotation(id) -> Mat3{Float64}

Contract: `1 <= id <= 21`; returns the constant matrix that rotates vectors from
J2000 (frame 1) to the built-in inertial frame `id`, equal to what CSPICE
`irfrot(1, id)` returns. Throws `SpiceError(:unknownframe)` for another ID.
"""
function inertial_rotation(id::Integer)::Mat3{Float64}
    1 <= id <= N_INERTIAL_FRAMES || throw(SpiceError(:unknownframe,
        "$id is not the ID of a built-in inertial frame"))
    return @inbounds INERTIAL_ROTATIONS[id]
end

# ---------------------------------------------------------------------------------
# Built-in non-inertial frames (ZZFDAT)
# ---------------------------------------------------------------------------------

"""
    BUILTIN_NONINERTIAL_FRAMES

The built-in non-inertial frames of CSPICE N0067 `zzfdat.c`, in table order:
`FrameDef(id, name, class, class ID, center)`.
"""
const BUILTIN_NONINERTIAL_FRAMES = FrameDef[
    FrameDef(10001, "IAU_MERCURY_BARYCENTER", 2, 1, 1),
    FrameDef(10002, "IAU_VENUS_BARYCENTER", 2, 2, 2),
    FrameDef(10003, "IAU_EARTH_BARYCENTER", 2, 3, 3),
    FrameDef(10004, "IAU_MARS_BARYCENTER", 2, 4, 4),
    FrameDef(10005, "IAU_JUPITER_BARYCENTER", 2, 5, 5),
    FrameDef(10006, "IAU_SATURN_BARYCENTER", 2, 6, 6),
    FrameDef(10007, "IAU_URANUS_BARYCENTER", 2, 7, 7),
    FrameDef(10008, "IAU_NEPTUNE_BARYCENTER", 2, 8, 8),
    FrameDef(10009, "IAU_PLUTO_BARYCENTER", 2, 9, 9),
    FrameDef(10010, "IAU_SUN", 2, 10, 10),
    FrameDef(10011, "IAU_MERCURY", 2, 199, 199),
    FrameDef(10012, "IAU_VENUS", 2, 299, 299),
    FrameDef(10013, "IAU_EARTH", 2, 399, 399),
    FrameDef(10014, "IAU_MARS", 2, 499, 499),
    FrameDef(10015, "IAU_JUPITER", 2, 599, 599),
    FrameDef(10016, "IAU_SATURN", 2, 699, 699),
    FrameDef(10017, "IAU_URANUS", 2, 799, 799),
    FrameDef(10018, "IAU_NEPTUNE", 2, 899, 899),
    FrameDef(10019, "IAU_PLUTO", 2, 999, 999),
    FrameDef(10020, "IAU_MOON", 2, 301, 301),
    FrameDef(10021, "IAU_PHOBOS", 2, 401, 401),
    FrameDef(10022, "IAU_DEIMOS", 2, 402, 402),
    FrameDef(10023, "IAU_IO", 2, 501, 501),
    FrameDef(10024, "IAU_EUROPA", 2, 502, 502),
    FrameDef(10025, "IAU_GANYMEDE", 2, 503, 503),
    FrameDef(10026, "IAU_CALLISTO", 2, 504, 504),
    FrameDef(10027, "IAU_AMALTHEA", 2, 505, 505),
    FrameDef(10028, "IAU_HIMALIA", 2, 506, 506),
    FrameDef(10029, "IAU_ELARA", 2, 507, 507),
    FrameDef(10030, "IAU_PASIPHAE", 2, 508, 508),
    FrameDef(10031, "IAU_SINOPE", 2, 509, 509),
    FrameDef(10032, "IAU_LYSITHEA", 2, 510, 510),
    FrameDef(10033, "IAU_CARME", 2, 511, 511),
    FrameDef(10034, "IAU_ANANKE", 2, 512, 512),
    FrameDef(10035, "IAU_LEDA", 2, 513, 513),
    FrameDef(10036, "IAU_THEBE", 2, 514, 514),
    FrameDef(10037, "IAU_ADRASTEA", 2, 515, 515),
    FrameDef(10038, "IAU_METIS", 2, 516, 516),
    FrameDef(10039, "IAU_MIMAS", 2, 601, 601),
    FrameDef(10040, "IAU_ENCELADUS", 2, 602, 602),
    FrameDef(10041, "IAU_TETHYS", 2, 603, 603),
    FrameDef(10042, "IAU_DIONE", 2, 604, 604),
    FrameDef(10043, "IAU_RHEA", 2, 605, 605),
    FrameDef(10044, "IAU_TITAN", 2, 606, 606),
    FrameDef(10045, "IAU_HYPERION", 2, 607, 607),
    FrameDef(10046, "IAU_IAPETUS", 2, 608, 608),
    FrameDef(10047, "IAU_PHOEBE", 2, 609, 609),
    FrameDef(10048, "IAU_JANUS", 2, 610, 610),
    FrameDef(10049, "IAU_EPIMETHEUS", 2, 611, 611),
    FrameDef(10050, "IAU_HELENE", 2, 612, 612),
    FrameDef(10051, "IAU_TELESTO", 2, 613, 613),
    FrameDef(10052, "IAU_CALYPSO", 2, 614, 614),
    FrameDef(10053, "IAU_ATLAS", 2, 615, 615),
    FrameDef(10054, "IAU_PROMETHEUS", 2, 616, 616),
    FrameDef(10055, "IAU_PANDORA", 2, 617, 617),
    FrameDef(10056, "IAU_ARIEL", 2, 701, 701),
    FrameDef(10057, "IAU_UMBRIEL", 2, 702, 702),
    FrameDef(10058, "IAU_TITANIA", 2, 703, 703),
    FrameDef(10059, "IAU_OBERON", 2, 704, 704),
    FrameDef(10060, "IAU_MIRANDA", 2, 705, 705),
    FrameDef(10061, "IAU_CORDELIA", 2, 706, 706),
    FrameDef(10062, "IAU_OPHELIA", 2, 707, 707),
    FrameDef(10063, "IAU_BIANCA", 2, 708, 708),
    FrameDef(10064, "IAU_CRESSIDA", 2, 709, 709),
    FrameDef(10065, "IAU_DESDEMONA", 2, 710, 710),
    FrameDef(10066, "IAU_JULIET", 2, 711, 711),
    FrameDef(10067, "IAU_PORTIA", 2, 712, 712),
    FrameDef(10068, "IAU_ROSALIND", 2, 713, 713),
    FrameDef(10069, "IAU_BELINDA", 2, 714, 714),
    FrameDef(10070, "IAU_PUCK", 2, 715, 715),
    FrameDef(10071, "IAU_TRITON", 2, 801, 801),
    FrameDef(10072, "IAU_NEREID", 2, 802, 802),
    FrameDef(10073, "IAU_NAIAD", 2, 803, 803),
    FrameDef(10074, "IAU_THALASSA", 2, 804, 804),
    FrameDef(10075, "IAU_DESPINA", 2, 805, 805),
    FrameDef(10076, "IAU_GALATEA", 2, 806, 806),
    FrameDef(10077, "IAU_LARISSA", 2, 807, 807),
    FrameDef(10078, "IAU_PROTEUS", 2, 808, 808),
    FrameDef(10079, "IAU_CHARON", 2, 901, 901),
    FrameDef(13000, "ITRF93", 2, 3000, 399),
    FrameDef(10081, "EARTH_FIXED", 4, 10081, 399),
    FrameDef(10082, "IAU_PAN", 2, 618, 618),
    FrameDef(10083, "IAU_GASPRA", 2, 9511010, 9511010),
    FrameDef(10084, "IAU_IDA", 2, 2431010, 2431010),
    FrameDef(10085, "IAU_EROS", 2, 2000433, 2000433),
    FrameDef(10086, "IAU_CALLIRRHOE", 2, 517, 517),
    FrameDef(10087, "IAU_THEMISTO", 2, 518, 518),
    FrameDef(10088, "IAU_MEGACLITE", 2, 519, 519),
    FrameDef(10089, "IAU_TAYGETE", 2, 520, 520),
    FrameDef(10090, "IAU_CHALDENE", 2, 521, 521),
    FrameDef(10091, "IAU_HARPALYKE", 2, 522, 522),
    FrameDef(10092, "IAU_KALYKE", 2, 523, 523),
    FrameDef(10093, "IAU_IOCASTE", 2, 524, 524),
    FrameDef(10094, "IAU_ERINOME", 2, 525, 525),
    FrameDef(10095, "IAU_ISONOE", 2, 526, 526),
    FrameDef(10096, "IAU_PRAXIDIKE", 2, 527, 527),
    FrameDef(10097, "IAU_BORRELLY", 2, 1000005, 1000005),
    FrameDef(10098, "IAU_TEMPEL_1", 2, 1000093, 1000093),
    FrameDef(10099, "IAU_VESTA", 2, 2000004, 2000004),
    FrameDef(10100, "IAU_ITOKAWA", 2, 2025143, 2025143),
    FrameDef(10101, "IAU_CERES", 2, 2000001, 2000001),
    FrameDef(10102, "IAU_PALLAS", 2, 2000002, 2000002),
    FrameDef(10103, "IAU_LUTETIA", 2, 2000021, 2000021),
    FrameDef(10104, "IAU_DAVIDA", 2, 2000511, 2000511),
    FrameDef(10105, "IAU_STEINS", 2, 2002867, 2002867),
    FrameDef(10106, "IAU_BENNU", 2, 2101955, 2101955),
    FrameDef(10107, "IAU_52_EUROPA", 2, 2000052, 2000052),
    FrameDef(10108, "IAU_NIX", 2, 902, 902),
    FrameDef(10109, "IAU_HYDRA", 2, 903, 903),
    FrameDef(10110, "IAU_RYUGU", 2, 2162173, 2162173),
    FrameDef(10111, "IAU_ARROKOTH", 2, 2486958, 2486958),
    FrameDef(10112, "IAU_DIDYMOS_BARYCENTER", 2, 20065803, 20065803),
    FrameDef(10113, "IAU_DIDYMOS", 2, 920065803, 920065803),
    FrameDef(10114, "IAU_DIMORPHOS", 2, 120065803, 120065803),
    FrameDef(10115, "IAU_DONALDJOHANSON", 2, 20052246, 20052246),
    FrameDef(10116, "IAU_EURYBATES", 2, 920003548, 920003548),
    FrameDef(10117, "IAU_EURYBATES_BARYCENTER", 2, 20003548, 20003548),
    FrameDef(10118, "IAU_QUETA", 2, 120003548, 120003548),
    FrameDef(10119, "IAU_POLYMELE", 2, 20015094, 20015094),
    FrameDef(10120, "IAU_LEUCUS", 2, 20011351, 20011351),
    FrameDef(10121, "IAU_ORUS", 2, 20021900, 20021900),
    FrameDef(10122, "IAU_PATROCLUS_BARYCENTER", 2, 20000617, 20000617),
    FrameDef(10123, "IAU_PATROCLUS", 2, 920000617, 920000617),
    FrameDef(10124, "IAU_MENOETIUS", 2, 120000617, 120000617),
]

"Contract: returns the (ID => definition, normalized name => ID) dictionaries of all built-in frames, inertial first; a later table entry overrides an earlier one."
function _builtin_frame_dicts()
    frames = Dict{Int32,FrameDef}()
    name2id = Dict{String,Int32}()
    for (i, (name, _, _)) in enumerate(INERTIAL_FRAME_DEFS)
        frames[Int32(i)] = FrameDef(i, name, 1, i, 0)
        name2id[normalize_name(name)] = i
    end
    for f in BUILTIN_NONINERTIAL_FRAMES
        frames[f.id] = f
        name2id[normalize_name(f.name)] = f.id
    end
    return frames, name2id
end

const (BUILTIN_FRAMES, BUILTIN_FRAMENAME2ID) = _builtin_frame_dicts()

# ---------------------------------------------------------------------------------
# Frames defined in the kernel pool
# ---------------------------------------------------------------------------------

"Contract: returns the single integer value of pool variable `key`, or `nothing` if the variable is missing, is not numeric, or has more than one element."
function _pool_int(pool::AbstractDict{String,PoolValue}, key::AbstractString)
    v = get(pool, key, nothing)
    (v isa Vector{Float64} && length(v) == 1 && isfinite(v[1]) &&
     abs(v[1]) <= typemax(Int32)) || return nothing
    return round(Int32, v[1])
end

"""
    _frame_item(pool, name2code, id, name, item) -> Union{Int32,Nothing}

Contract: value of the frame definition variable `FRAME_<id>_<item>`, or of
`FRAME_<name>_<item>` if the first is absent (the two forms of `zzdynbid`). A string
value is translated as a body name. `nothing` if neither form gives one integer.
"""
function _frame_item(pool::AbstractDict{String,PoolValue}, name2code, id::Int32,
    name::AbstractString, item::AbstractString)
    for key in ("FRAME_$(id)_$(item)", "FRAME_$(name)_$(item)")
        length(key) <= _TK_MAXNAME || continue
        v = get(pool, key, nothing)
        v === nothing && continue
        length(v) == 1 || return nothing
        v isa Vector{String} && return _body_code_lookup(name2code, v[1])
        return _pool_int(pool, key)
    end
    return nothing
end

"""
    build_frames(pool) -> (frames, framename2id)

Contract: reads the frame definitions of a frame kernel from `pool`. A frame exists
for each `FRAME_<id>_NAME` that has its class, class ID and center, each under the
`FRAME_<id>_<item>` key or the `FRAME_<name>_<item>` key; a center may be a body
name. Returns `frames::Dict{Int32,FrameDef}` and `framename2id::Dict{String,Int32}`
keyed by normalized name, the latter from the `FRAME_<NAME> = id` variables and from
the frame names. An incomplete definition is left out, so the frame is unknown.
"""
function build_frames(pool::AbstractDict{String,PoolValue})
    frames = Dict{Int32,FrameDef}()
    framename2id = Dict{String,Int32}()
    name2code, _ = build_names(pool)
    for (key, val) in pool
        (startswith(key, "FRAME_") && endswith(key, "_NAME") && length(key) > 11) || continue
        (val isa Vector{String} && length(val) == 1) || continue
        id = tryparse(Int32, SubString(key, 7, lastindex(key) - 5))
        id === nothing && continue
        name = String(strip(val[1], ' '))
        isempty(name) && continue
        class = _frame_item(pool, name2code, id, name, "CLASS")
        classid = _frame_item(pool, name2code, id, name, "CLASS_ID")
        center = _frame_item(pool, name2code, id, name, "CENTER")
        (class === nothing || classid === nothing || center === nothing) && continue
        frames[id] = FrameDef(id, name, class, classid, center)
    end
    # The `FRAME_<NAME> = id` form is the one `namfrm` reads; it has priority.
    for (id, f) in frames
        framename2id[normalize_name(f.name)] = id
    end
    for (key, _) in pool
        startswith(key, "FRAME_") || continue
        id = _pool_int(pool, key)
        (id === nothing || !haskey(frames, id)) && continue
        suffix = SubString(key, 7)
        (endswith(suffix, "_CLASS") || endswith(suffix, "_CLASS_ID") ||
         endswith(suffix, "_CENTER")) && continue
        framename2id[normalize_name(suffix)] = id
    end
    return frames, framename2id
end

# ---------------------------------------------------------------------------------
# Queries
# ---------------------------------------------------------------------------------

"""
    frame_def(ks, id) -> FrameDef

Contract: definition of frame `id`: a pool-defined frame first, then a built-in
frame, then the rule of `frinfo` that an ID in 13000..13999 is a class 2 Earth frame
with class ID `id - 10000`. Throws `SpiceError(:unknownframe)` otherwise. Reads
`ks.frames` only.
"""
function frame_def(ks::KernelSet, id::Integer)::FrameDef
    if typemin(Int32) <= id <= typemax(Int32)
        i = Int32(id)
        f = get(ks.frames, i, nothing)
        f === nothing || return f
        f = get(BUILTIN_FRAMES, i, nothing)
        f === nothing || return f
        13000 <= i < 14000 && return FrameDef(i, "", 2, i - Int32(10000), 399)
    end
    throw(SpiceError(:unknownframe, "there is no definition of the frame with ID $id"))
end

"""
    frame_id(ks, name) -> Int32

Contract: ID of the frame `name` (case and outer blanks do not matter): a
pool-defined frame first, then a built-in frame. Throws `SpiceError(:unknownframe)`
otherwise. An integer argument is returned as `Int32`. Reads `ks.framename2id` only.
"""
function frame_id(ks::KernelSet, name::AbstractString)::Int32
    nn = normalize_name(name)
    id = get(ks.framename2id, nn, nothing)
    id === nothing || return id
    id = get(BUILTIN_FRAMENAME2ID, nn, nothing)
    id === nothing || return id
    throw(SpiceError(:unknownframe, "the frame name '$name' is not known"))
end
frame_id(::KernelSet, id::Integer) = Int32(id)

"""
    frame_name(ks, id) -> Union{String,Nothing}

Contract: name of the frame `id`: a pool-defined frame first, then a built-in frame,
or `nothing` if the ID has no name. Reads `ks.frames` only.
"""
function frame_name(ks::KernelSet, id::Integer)::Union{String,Nothing}
    typemin(Int32) <= id <= typemax(Int32) || return nothing
    i = Int32(id)
    f = get(ks.frames, i, nothing)
    f === nothing || return f.name
    f = get(BUILTIN_FRAMES, i, nothing)
    return f === nothing ? nothing : f.name
end
