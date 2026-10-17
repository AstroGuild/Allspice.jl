# Ephemeris interface: geometric states from the SPK segments of a `KernelSet`, with
# time derivatives of any order. The aberration corrections are in `aberration.jl`.

"Speed of light in km/s, the value of `clight_c`."
const CLIGHT = 299792.458

"Number of bodies of a target chain that `spk_jet` holds; a longer chain goes to `spk_jet_long`."
const MAXCHAIN = 10

"""
    rotate_jet(Rjet::NTuple{M}, rjet::NTuple{M}) -> NTuple{M}

Return the jet of the product `R(t) * r(t)` from the jets of `R` and `r` (ordinary
derivatives, orders 0 to M-1), by the Leibniz rule.
"""
@inline function rotate_jet(Rjet::NTuple{M,Any}, rjet::NTuple{M,Any}) where {M}
    return ntuple(Val(M)) do i
        k = i - 1
        acc = _m_mxv(Rjet[1], rjet[i])
        for j in 1:k
            acc += binomial(k, j) * _m_mxv(Rjet[j+1], rjet[i-j])
        end
        acc
    end
end

"""
    frame_jet(ks, from::Int32, to::Int32, et, ::Val{N}) -> NTuple{N+1,Mat3}

Return the rotation matrix from frame `from` to frame `to` at `et` and its time
derivatives to order N, for any epoch number type (`rotation_jet` of `rotations.jl`
takes a plain epoch only, so a dual-number epoch goes through `frame_rotation`).
"""
@inline frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Float64, v::Val) =
    rotation_jet(ks, from, to, et, v)
@inline frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Float64, v::Val{0}) =
    rotation_jet(ks, from, to, et, v)
@inline frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Real, ::Val{0}) =
    (frame_rotation(ks, from, to, et),)
@inline function frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Real, v::Val{N}) where {N}
    d = map(derivatives, frame_rotation(ks, from, to, taylor_var(et, v)))
    return ntuple(k -> map(x -> x[k], d), Val(N + 1))
end

"""
    segment_jet(ks, seg::Segment, et, ::Val{N}) -> NTuple{N+1,Vec3}

Return the jet of `seg.target` relative to `seg.center` at `et` in J2000; a segment in
another frame is rotated with the rotation jet of that frame.
"""
@inline function segment_jet(ks::KernelSet, seg::Segment, et::T, ::Val{N}) where {T<:Real,N}
    # Type 10 (two-line elements) needs the leap seconds: spk_types/type10.jl.
    jet = seg.dtype == Int32(10) ? _t10_jet(_lsk(ks), ks.daf[seg.file], seg, et, Val(N)) :
          spk_eval(ks.daf[seg.file], seg, et, Val(N))
    seg.frame == Int32(1) && return jet
    return rotate_jet(frame_jet(ks, seg.frame, Int32(1), et, Val(N)), jet)
end

@inline function chain_index(chain::NTuple{MAXCHAIN,Int32}, n::Int, body::Int32)
    for k in 1:n
        chain[k] == body && return k
    end
    return 0
end

"Text that lists the SPK coverage intervals of `body` in `ks` (ET seconds)."
function coverage_text(ks::KernelSet, body::Integer)
    segs = get(ks.spk, Int32(body), nothing)
    segs === nothing && return "no SPK segment is loaded for this body"
    parts = ["[$(s.t0), $(s.t1)] (center $(s.center))" for s in Iterators.take(segs, 8)]
    more = length(segs) > 8 ? " and $(length(segs) - 8) more" : ""
    return "available coverage: " * join(parts, ", ") * more
end

"""
    throw_no_connection(ks, targ, obs, t, tend, oend) -> Union{}

Throw the [`SpkInsufficientDataError`](@ref) for two chains that do not meet at `t`;
`tend` and `oend` are the last bodies of the target chain and of the observer chain.
CSPICE signals NOLOADEDFILES with no SPK file and SPKINSUFFDATA for every other case.
"""
@noinline function throw_no_connection(ks::KernelSet, targ, obs, t::Float64, tend, oend)
    isempty(ks.spk) && throw(SpkInsufficientDataError(:noloadedfiles, "no SPK file is " *
        "loaded (needed for target $targ relative to observer $obs at ET $t s); load " *
        "one with furnsh", targ, obs, t, nothing, false))
    for body in (oend, tend)
        if body != 0 && haskey(ks.spk, Int32(body))
            throw(SpkInsufficientDataError(:spkinsuffdata, "no SPK segment covers body " *
                "$body at ET $t s (needed for target $targ relative to observer $obs); " *
                coverage_text(ks, body), targ, obs, t, body, true))
        end
    end
    # No end of a chain is in the files at another epoch: the end that is not the solar
    # system barycenter is the body without data.
    missing_body = tend != 0 ? tend : oend != 0 ? oend : nothing
    throw(SpkInsufficientDataError(:spkinsuffdata, "no chain of SPK segments connects " *
        "target $targ and observer $obs at ET $t s: the target chain ends at body $tend " *
        "and the observer chain ends at body $oend", targ, obs, t, missing_body, false))
end

# CSPICE does not return on a cycle; `spk_jet` stops after as many steps as there are
# bodies with SPK data, which no chain without a cycle exceeds.
@noinline throw_chain_cycle(targ, obs, t::Float64) = throw(SpkInsufficientDataError(
    :spkinsuffdata, "the SPK chain between target $targ and observer $obs at ET $t s " *
    "has a cycle", targ, obs, t, nothing, false))

"""
    spk_jet(ks, targ::Int32, et, ref::Int32, obs::Int32, ::Val{N}) -> NTuple{N+1,Vec3}

Return the geometric position of `targ` relative to `obs` in frame `ref` at `et` and
its time derivatives to order N (no aberration correction). The function does not
allocate, except for a target chain of more than MAXCHAIN bodies. Throws
[`SpkInsufficientDataError`](@ref).
"""
function spk_jet(ks::KernelSet, targ::Int32, et::T, ref::Int32, obs::Int32,
        ::Val{N}) where {T<:Real,N}
    t = primal(et)
    # Pass 1: bodies and segment indices of the target chain up to the observer or the
    # root, without evaluation.
    tbody = ntuple(_ -> Int32(0), Val(MAXCHAIN))
    tseg = ntuple(_ -> 0, Val(MAXCHAIN))
    tbody = Base.setindex(tbody, targ, 1)
    nt = 1
    body = targ
    while body != 0 && body != obs
        # invariant: tbody[1:nt] are the bodies reached, tseg[1:nt-1] their segments
        nt == MAXCHAIN && return spk_jet_long(ks, targ, et, ref, obs, Val(N))
        segs = get(ks.spk, body, nothing)
        segs === nothing && break
        i = find_segment(segs, t)
        i == 0 && break
        tseg = Base.setindex(tseg, i, nt)
        body = segs[i].center
        nt += 1
        tbody = Base.setindex(tbody, body, nt)
    end
    # Pass 2: walk the observer chain until it meets the target chain.
    ojet = ntuple(_ -> zero(Vec3{T}), Val(N + 1))
    body = obs
    depth = 0
    common = chain_index(tbody, nt, body)
    while common == 0
        segs = get(ks.spk, body, nothing)
        segs === nothing && throw_no_connection(ks, targ, obs, t, tbody[nt], body)
        i = find_segment(segs, t)
        i == 0 && throw_no_connection(ks, targ, obs, t, tbody[nt], body)
        depth += 1
        depth > length(ks.spk) && throw_chain_cycle(targ, obs, t)
        seg = segs[i]
        ojet = map(+, ojet, segment_jet(ks, seg, et, Val(N)))
        body = seg.center
        common = chain_index(tbody, nt, body)
    end
    # Pass 3: target chain up to the common body.
    tjet = ntuple(_ -> zero(Vec3{T}), Val(N + 1))
    for k in 1:common-1
        seg = ks.spk[tbody[k]][tseg[k]]
        tjet = map(+, tjet, segment_jet(ks, seg, et, Val(N)))
    end
    jet = map(-, tjet, ojet)
    ref == Int32(1) && return jet
    return rotate_jet(frame_jet(ks, Int32(1), ref, et, Val(N)), jet)
end

"""
    spk_jet_long(ks, targ::Int32, et, ref::Int32, obs::Int32, ::Val{N}) -> NTuple{N+1,Vec3}

`spk_jet` for a target chain of more than MAXCHAIN bodies; it allocates. CSPICE
(`SPKGEO`) has no limit on the length of a chain either.
"""
@noinline function spk_jet_long(ks::KernelSet, targ::Int32, et::T, ref::Int32, obs::Int32,
        ::Val{N}) where {T<:Real,N}
    t = primal(et)
    tbody = Int32[targ]
    tseg = Segment[]
    body = targ
    while body != 0 && body != obs
        # invariant: tseg[k] is the segment from tbody[k] to tbody[k+1]
        segs = get(ks.spk, body, nothing)
        segs === nothing && break
        i = find_segment(segs, t)
        i == 0 && break
        length(tseg) < length(ks.spk) || throw_chain_cycle(targ, obs, t)
        push!(tseg, segs[i])
        body = segs[i].center
        push!(tbody, body)
    end
    ojet = ntuple(_ -> zero(Vec3{T}), Val(N + 1))
    body = obs
    depth = 0
    common = findfirst(==(body), tbody)
    while common === nothing
        segs = get(ks.spk, body, nothing)
        segs === nothing && throw_no_connection(ks, targ, obs, t, last(tbody), body)
        i = find_segment(segs, t)
        i == 0 && throw_no_connection(ks, targ, obs, t, last(tbody), body)
        depth += 1
        depth > length(ks.spk) && throw_chain_cycle(targ, obs, t)
        ojet = map(+, ojet, segment_jet(ks, segs[i], et, Val(N)))
        body = segs[i].center
        common = findfirst(==(body), tbody)
    end
    tjet = ntuple(_ -> zero(Vec3{T}), Val(N + 1))
    for k in 1:common-1
        tjet = map(+, tjet, segment_jet(ks, tseg[k], et, Val(N)))
    end
    jet = map(-, tjet, ojet)
    ref == Int32(1) && return jet
    return rotate_jet(frame_jet(ks, Int32(1), ref, et, Val(N)), jet)
end

@inline function check_epoch(et::Real)
    isfinite(primal(et)) || throw(SpiceError(:badarg, "the epoch is not finite"))
    return float(et)
end

"A body or a frame given by name or by integer ID."
const NameOrID = Union{AbstractString,Integer}

function spkder(ks::KernelSet, targ::NameOrID, et::Real, ref::NameOrID, obs::NameOrID,
        ::Val{N}) where {N}
    N isa Int && N >= 0 || throw(SpiceError(:badarg, "spkder: the order must be an integer >= 0"))
    return spk_jet(ks, body_code(ks, targ), check_epoch(et), frame_id(ks, ref),
        body_code(ks, obs), Val(N))
end
spkder(ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID, obs::NameOrID, v::Val) =
    spkder(kernelset(ctx), targ, et, ref, obs, v)
spkder(targ::NameOrID, et::Real, ref::NameOrID, obs::NameOrID, v::Val) =
    spkder(default_context(), targ, et, ref, obs, v)

"""
    spkgeo([ctx,] targ::Integer, et, ref, obs::Integer) -> (state, lt)

Return the geometric state of body `targ` relative to body `obs` in frame `ref` at
`et`, and the one-way light time, as `spkgeo_c`.

`state` is an `SVector{6}`: position in km, then velocity in km/s. `lt` is the
distance divided by the speed of light, in seconds. `et` is in seconds past J2000
TDB. `targ` and `obs` are NAIF ID codes; `ref` is a frame name or a frame ID.

Throws `SpiceError(:unknownframe)` if `ref` is not a known frame name,
`SpiceError(:badarg)` if `et` is not finite, and [`SpkInsufficientDataError`](@ref) if
no SPK file is loaded (`:noloadedfiles`) or the loaded SPK files do not connect `targ`
and `obs` at `et` (`:spkinsuffdata`). As in CSPICE, `targ == obs` gives a zero state
and a zero light time before `ref` and `et` are checked; so does [`spkezr`](@ref) with
the correction `"NONE"`.

See also [`spkezr`](@ref), [`spkssb`](@ref).
"""
function spkgeo(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID, obs::Integer)
    # As in CSPICE, a body relative to itself is zero before `ref` and `et` are checked.
    targ == obs && return zero(SVector{6,typeof(float(et))}), zero(float(et))
    r, v = spk_jet(ks, Int32(targ), check_epoch(et), frame_id(ks, ref), Int32(obs), Val(1))
    return vcat(r, v), norm(r) / CLIGHT
end
spkgeo(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID, obs::Integer) =
    spkgeo(kernelset(ctx), targ, et, ref, obs)
spkgeo(targ::Integer, et::Real, ref::NameOrID, obs::Integer) =
    spkgeo(default_context(), targ, et, ref, obs)

"""
    spkgps([ctx,] targ::Integer, et, ref, obs::Integer) -> (pos, lt)

Return the geometric position of body `targ` relative to body `obs` in frame `ref` at
`et`, and the one-way light time, as `spkgps_c`.

`pos` is an `SVector{3}` in km and `lt` is in seconds. The arguments and the errors are
those of [`spkgeo`](@ref). As in CSPICE, `targ == obs` gives a zero position and a zero
light time before `ref` and `et` are checked.

See also [`spkezp`](@ref), [`spkpos`](@ref).
"""
function spkgps(ks::KernelSet, targ::Integer, et::Real, ref::NameOrID, obs::Integer)
    targ == obs && return zero(Vec3{typeof(float(et))}), zero(float(et))
    (r,) = spk_jet(ks, Int32(targ), check_epoch(et), frame_id(ks, ref), Int32(obs), Val(0))
    return r, norm(r) / CLIGHT
end
spkgps(ctx::SpiceContext, targ::Integer, et::Real, ref::NameOrID, obs::Integer) =
    spkgps(kernelset(ctx), targ, et, ref, obs)
spkgps(targ::Integer, et::Real, ref::NameOrID, obs::Integer) =
    spkgps(default_context(), targ, et, ref, obs)

"""
    spkezr([ctx,] targ, et, ref, abcorr, obs) -> (state, lt)

Return the state of `targ` relative to `obs` in frame `ref` at `et`, and the one-way
light time, as `spkezr_c`.

`state` is an `SVector{6}`: position in km, then velocity in km/s. `lt` is in seconds.

- `targ`, `obs`: body name or NAIF ID code.
- `et`: epoch of the observer in seconds past J2000 TDB (see [`str2et`](@ref)).
- `ref`: frame name or frame ID. For a frame that is not inertial and a correction
  other than `"NONE"`, the frame is evaluated at `et` corrected for the light time
  between the observer and the frame center, as in CSPICE.
- `abcorr`: `"NONE"` (geometric), `"LT"`, `"LT+S"`, `"CN"`, `"CN+S"` for reception, or
  the same with the prefix `X` for transmission. Letter case and blanks are ignored.
  An [`Aberration`](@ref) value is also accepted.

The call does not allocate when the bodies and the frame are integer IDs and `abcorr`
is an `Aberration`. With ForwardDiff loaded, `et` can be a dual number.

The SPK segments connect a body to its center, and that center to the next one. The
call needs such a chain from `targ` and from `obs` to a common body.

Throws `SpiceError`:
- `:idcodenotfound` if a body name is not known, `:unknownframe` if a frame name is
  not known;
- `:invalidoption` if `abcorr` is not valid, `:badarg` if `et` is not finite;
- `:noloadedfiles` if no SPK file is loaded;
- `:spkinsuffdata` if the loaded SPK files hold no chain between `targ` and `obs` at
  the epoch needed: a body is in no file, or the epoch is outside its segments. These
  two are a [`SpkInsufficientDataError`](@ref), which holds the body and the epoch;
- the kinds of [`pxform`](@ref) if `ref` is not inertial and its orientation data are
  not loaded.

# Examples
```julia
furnsh("naif0012.tls")
furnsh("de440s.bsp")
et = str2et("2025-03-01T12:00:00")
state, lt = spkezr("MARS BARYCENTER", et, "ECLIPJ2000", "LT+S", "EARTH")
r, v = state[1:3], state[4:6]                  # km, km/s
state, lt = spkezr(4, et, 17, Aberration("LT+S"), 399)   # the same with IDs
```

See also [`spkpos`](@ref), [`spkgeo`](@ref), [`spkder`](@ref), [`furnsh`](@ref).
"""
spkezr(ks::KernelSet, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString,
    obs::NameOrID) = spkezr(ks, targ, et, ref, Aberration(abcorr), obs)
spkezr(ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString,
    obs::NameOrID) = spkezr(kernelset(ctx), targ, et, ref, Aberration(abcorr), obs)
spkezr(targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString, obs::NameOrID) =
    spkezr(default_context(), targ, et, ref, abcorr, obs)

"""
    spkpos([ctx,] targ, et, ref, abcorr, obs) -> (pos, lt)

Return the position of `targ` relative to `obs` in frame `ref` at `et`, and the
one-way light time, as `spkpos_c`.

`pos` is an `SVector{3}` in km and `lt` is in seconds. The arguments and the errors
are those of [`spkezr`](@ref), except that an `abcorr` that is not valid throws
`SpiceError(:spkinvalidoption)`, as in CSPICE.

As in CSPICE, the `CN` corrections make exactly three light-time iterations here, and
[`spkezr`](@ref) iterates to convergence. The position can therefore differ from the
position part of `spkezr` at rounding level.

# Examples
```julia
et = str2et("2025-03-01T12:00:00")
pos, lt = spkpos("SUN", et, "J2000", "LT+S", "EARTH")    # km, s
```

See also [`spkezr`](@ref), [`spkder`](@ref).
"""
spkpos(ks::KernelSet, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString,
    obs::NameOrID) = spkpos(ks, targ, et, ref, _spk_abcorr(abcorr), obs)
spkpos(ctx::SpiceContext, targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString,
    obs::NameOrID) = spkpos(kernelset(ctx), targ, et, ref, _spk_abcorr(abcorr), obs)
spkpos(targ::NameOrID, et::Real, ref::NameOrID, abcorr::AbstractString, obs::NameOrID) =
    spkpos(default_context(), targ, et, ref, abcorr, obs)
