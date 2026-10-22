# CK interface: pointing queries, CK meta data, coverage of a CK file, and the rotation
# jets that the frame and state interfaces use when a chain holds a CK frame.

# --- rotation jets with CK frames ---------------------------------------------------------

"""
    frame_jet(ks, from::Int32, to::Int32, et::Float64, Val(1)) -> (R, dR)

Order 1 rotation jet at a plain epoch, the one behind `sxform`, `spkezr` and `spkgeo`.
If CK kernels are loaded and one of the two frame chains holds a CK frame, `dR` comes
from the angular velocity of the CK segment, as `sxform_c` of CSPICE gives it
(`_ck_state_rotation`). A chain with a dynamic or a switch frame goes the same way, and
`dR` is the derivative block of `zzdynfrm_`. Otherwise the result is `rotation_jet`,
unchanged.
"""
@inline function frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Float64, v::Val{1})
    (isempty(ks.ck) && isempty(ks.dyn) && isempty(ks.switch)) &&
        return rotation_jet(ks, from, to, et, v)
    return _ck_frame_jet(ks, from, to, et)
end

"""
    _ck_frame_jet(ks, from::Int32, to::Int32, et::Float64) -> (R, dR)

Contract: `frame_jet` of order 1 for a snapshot that has CK, dynamic or switch frames:
the state transformation of CSPICE if one of the two chains holds such a frame, and
`rotation_jet` otherwise. Not inlined, so that the callers of `frame_jet` stay small
for a snapshot without such frames.
"""
@noinline function _ck_frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Float64)
    (_ck_in_chain(ks, from) || _ck_in_chain(ks, to)) ||
        return rotation_jet(ks, from, to, et, Val(1))
    x = _ck_state_rotation(ks, from, to, et)
    return map(e -> e.v, x), map(e -> e.d[1], x)
end

"""
    _exact_frame_jet(ks, from::Int32, to::Int32, et, Val(N)) -> NTuple{N+1,Mat3}

Rotation matrix from frame `from` to frame `to` at `et` and its exact time derivatives
to order N: the derivatives of the function that `frame_rotation` evaluates, also for
a CK frame.
"""
@inline _exact_frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Float64, v::Val) =
    rotation_jet(ks, from, to, et, v)
@inline _exact_frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Real, v::Val) =
    frame_jet(ks, from, to, et, v)

# --- pointing -----------------------------------------------------------------------------

"Throw `SpiceError(:valueoutofrange)` unless the tolerance `tol` is a number `>= 0`."
@inline function _ck_check_tol(tol::Real)
    tol >= 0 || throw(SpiceError(:valueoutofrange,
        "the tolerance must be non-negative but was $tol"))
    return Float64(tol)
end

"ET at which `ckgp` changes frames: 0 for two inertial frames, else the ET of `clkout`."
@inline function _ck_change_epoch(ks::KernelSet, inst::Int32, refreq::Int32, refseg::Int32,
        clkout::Float64)
    (frame_def(ks, refreq).class == 1 && frame_def(ks, refseg).class == 1) && return 0.0
    return sct2e(ks, _ck_sclkid(ks, inst), clkout)
end

"""
    ckgp([ctx,] inst, sclkdp, tol, ref) -> Union{Nothing,Tuple{SMatrix{3,3},Float64}}

Return the pointing of the instrument or structure `inst` at the encoded spacecraft
clock value `sclkdp`, as `ckgp_c`.

- `inst`: CK ID of the instrument or structure.
- `sclkdp`: encoded spacecraft clock, in ticks (see [`sce2c`](@ref)).
- `tol`: time tolerance in ticks, `>= 0`. A pointing within `tol` of `sclkdp` is
  accepted.
- `ref`: name or ID of the reference frame of the result.

The result is `(cmat, clkout)`: `cmat` is an `SMatrix{3,3}` that rotates a vector from
`ref` to the instrument frame, and `clkout` is the clock value, in ticks, of the
pointing that was used. The result is `nothing` if no loaded CK segment gives a
pointing; `nothing` takes the place of the `found` flag of CSPICE. Segments with and
without angular velocity are searched in priority order.

If `ref` is not the frame of the segment, the matrix is changed with the rotation
between the two frames at the epoch of `clkout`. This needs the SCLK kernel of the
spacecraft, unless both frames are inertial.

Throws `SpiceError(:valueoutofrange)` if `tol < 0`. Throws
`SpiceError(:unknownframe)` if `ref` is not a known frame name; `ref` is read only
when a pointing is found. The change of frame throws the kinds of [`pxform`](@ref),
and the kinds of [`sct2e`](@ref) if the SCLK kernel is not loaded: [`ckfrot`](@ref)
returns `nothing` in this case, and `ckgp` throws, both as in CSPICE.

# Examples
```julia
furnsh("naif0012.tls"); furnsh("mission.tsc"); furnsh("mission.bc")
et = str2et("2026-10-16T12:00:00")
sclkdp = sce2c(-77, et)
r = ckgp(-77000, sclkdp, 0.0, "J2000")
if r !== nothing
    cmat, clkout = r
    boresight_j2000 = cmat' * [0.0, 0.0, 1.0]
end
```

See also [`ckgpav`](@ref), [`ckcov`](@ref), [`pxform`](@ref).
"""
function ckgp(ks::KernelSet, inst::Integer, sclkdp::Real, tol::Real, ref::NameOrID)
    i = Int32(inst)
    found, cmat, _, clkout, refseg = _ck_search(ks, i, Float64(sclkdp), _ck_check_tol(tol), false)
    found || return nothing
    refreq = frame_id(ks, ref)
    if refreq != refseg
        et = _ck_change_epoch(ks, i, refreq, refseg, clkout)
        cmat = _m_mxm(cmat, frame_rotation(ks, refreq, refseg, et))
    end
    return cmat, clkout
end

"""
    ckgpav([ctx,] inst, sclkdp, tol, ref) -> Union{Nothing,Tuple{SMatrix,SVector,Float64}}

Return the pointing and the angular velocity of the instrument or structure `inst` at
the encoded spacecraft clock value `sclkdp`, as `ckgpav_c`.

The result is `(cmat, av, clkout)`, where `av` is an `SVector{3}`: the angular velocity
of the instrument frame relative to `ref`, in rad/s, with components in `ref`. Only
segments with angular velocity are used, and the result is `nothing` if none gives a
pointing. The other arguments, results and errors are those of [`ckgp`](@ref).
"""
function ckgpav(ks::KernelSet, inst::Integer, sclkdp::Real, tol::Real, ref::NameOrID)
    i = Int32(inst)
    found, cmat, av, clkout, refseg = _ck_search(ks, i, Float64(sclkdp), _ck_check_tol(tol), true)
    found || return nothing
    refreq = frame_id(ks, ref)
    if refreq != refseg
        et = _ck_change_epoch(ks, i, refreq, refseg, clkout)
        x = _ck_state_rotation(ks, refreq, refseg, et)
        rot = map(e -> e.v, x)
        drdt = map(e -> e.d[1], x)
        omega = _m_mtxm(drdt, rot)           # `xf2rav_`
        cmat = _m_mxm(cmat, rot)
        av = Vec3{Float64}(omega[3, 2], omega[1, 3], omega[2, 1]) + _m_mtxv(rot, av)
    end
    return cmat, av, clkout
end

# --- one link of a CK frame -----------------------------------------------------------------

"Clock of the CK ID `inst`, or `nothing` if no CK is loaded (`ckhave_`) or no SCLK kernel of the clock is loaded (`zzsclk_`)."
@inline function _ck_link_clock(ks::KernelSet, inst::Int32)
    isempty(ks.ck) && return nothing
    return get(ks.sclk, _ck_sclkid(ks, inst), nothing)
end

"""
    ckfrot([ctx,] inst, et) -> Union{Nothing,Tuple{SMatrix{3,3},Int}}

Return the rotation from the CK frame of the instrument or structure `inst` to its base
frame, the reference frame of the CK segment, at `et`, and the frame ID of the base
frame, as `ckfrot_c`.

- `inst`: CK ID (frame class ID) of the CK frame.
- `et`: seconds past J2000 TDB.

For the result `(rot, ref)`, `rot * v` gives in frame `ref` the components of a position
`v` given in the CK frame. `ref` is the reference frame of the CK segment that gives the
pointing: the clock of `inst` (see [`ckmeta`](@ref)) converts `et` to ticks, and the
loaded segments of `inst` are searched in priority order with tolerance zero, with and
without angular velocity. This is the link of the frame chain that [`pxform`](@ref)
follows for a CK frame. As in `pxform`, `et` can be a dual number when ForwardDiff is
loaded; the elements of `rot` then are dual numbers.

`nothing` takes the place of the `found` flag of CSPICE. The checks come in this
order: `SpiceError(:badarg)` if `et` is not finite; `nothing` if no CK is loaded or the
SCLK kernel of the clock is not loaded, as in CSPICE; `SpiceError(:valueoutofrange)`
if `et` is outside the range of the clock; `nothing` if no segment covers `et`.

# Examples
```julia
furnsh("naif0012.tls"); furnsh("mission.tsc"); furnsh("mission.bc")
et = str2et("2026-10-16T12:00:00")
r = ckfrot(-77000, et)
if r !== nothing
    rot, ref = r
    boresight = rot * [0.0, 0.0, 1.0]        # in the frame `ref`
end
```

See also [`ckfxfm`](@ref), [`ckgp`](@ref).
"""
function ckfrot(ks::KernelSet, inst::Integer, et::Real)
    i = Int32(inst)
    t = check_epoch(et)
    clock = _ck_link_clock(ks, i)
    clock === nothing && return nothing
    found, cmat, _, _, ref = _ck_search(ks, i, _sce2c(ks, clock, t), 0.0, false)
    found || return nothing
    return Mat3{eltype(cmat)}(transpose(cmat)), Int(ref)
end

"""
    ckfxfm([ctx,] inst, et) -> Union{Nothing,Tuple{SMatrix{6,6},Int}}

Return the state transformation from the CK frame of the instrument or structure `inst`
to its base frame at `et`, and the frame ID of the base frame, as `ckfxfm_c`.

The matrix is `[R 0; dR R]`: `R` is the rotation of [`ckfrot`](@ref), and `dR` (1/s)
comes from the angular velocity of the CK segment, as [`sxform`](@ref) gives it for a CK
frame. Only segments with angular velocity are searched, so the result is `nothing` also
where only a segment without angular velocity covers `et`. The other `nothing` cases,
the arguments and the errors are those of [`ckfrot`](@ref), except that `et` is a plain
number.
"""
function ckfxfm(ks::KernelSet, inst::Integer, et::Real)
    i = Int32(inst)
    t = Float64(check_epoch(et))
    clock = _ck_link_clock(ks, i)
    clock === nothing && return nothing
    found, cmat, av, _, ref = _ck_search(ks, i, _sce2c(ks, clock, t), 0.0, true)
    found || return nothing
    # `rav2xf_` gives the transformation from the base frame to the CK frame, and
    # `invstm_` inverts it by transposing its blocks.
    x = rav2xf(cmat, av)
    r = Mat3{Float64}(transpose(cmat))
    dr = Mat3{Float64}(transpose(x[SVector(4, 5, 6), SVector(1, 2, 3)]))
    return vcat(hcat(r, zero(r)), hcat(dr, r)), Int(ref)
end

"""
    ckmeta([ctx,] ckid, meta) -> Int

Return the ID of the spacecraft clock or of the body that belongs to the CK ID
`ckid`, as `ckmeta_c`.

`meta` is `"SCLK"` for the clock ID or `"SPK"` for the body ID; letter case is
ignored. The kernel-pool variables `CK_<ckid>_SCLK` and `CK_<ckid>_SPK` give the
answer. Without the variable the answer is `ckid ÷ 1000` for `ckid <= -1000`, and 0
otherwise.

Throws `SpiceError(:unknownckmeta)` if `meta` is not one of the two.

# Examples
```jldoctest
julia> ckmeta(-77000, "SCLK"), ckmeta(5, "SPK")
(-77, 0)
```
"""
function ckmeta(ks::KernelSet, ckid::Integer, meta::AbstractString)
    item = uppercase(strip(meta))
    k = item == "SCLK" ? 1 : item == "SPK" ? 2 : throw(SpiceError(:unknownckmeta,
        "the CK meta data item \"$meta\" is not SCLK or SPK"))
    id = Int32(ckid)
    m = get(ks.ckmeta, id, nothing)
    return Int(m === nothing ? _ckmeta_default(id) : m[k])
end

# --- contents of a CK file ----------------------------------------------------------------

"""
    _ck_open(f, path, who) -> f(daf)

Open `path` as a CK file, with the checks of `_require_binary_kernel`, and return
`f(daf)`. The map of the file is removed when `f` returns, so the result of `f` must not
refer to `daf.data`. Throws `SpiceError(:invalidfiletype)` for another DAF kernel.
"""
function _ck_open(f, path::AbstractString, who::AbstractString)
    _require_binary_kernel(who, path, "CK")
    daf = open_daf(path)
    try
        daf.kind === :CK || throw(SpiceError(:invalidfiletype,
            "the file $path is a binary $(daf.kind) kernel; a CK file is needed"))
        return f(daf)
    finally
        _unmap!(daf.data)
    end
end

"""
    ckobj(ck) -> Vector{Int}

Return the CK IDs of the instruments and structures that have a segment in the CK file
`ck`, sorted and without duplicates, as `ckobj_c`.

The file is read from disk and need not be loaded.

Throws `SpiceError`:
- `:filenotfound` if the file does not exist;
- `:invalidformat` if it is a transfer file (a binary kernel in text form),
  `:invalidarchtype` if it is not a DAF;
- `:invalidfiletype` if it is an SPK or a binary PCK;
- `:baddaf` if the DAF is damaged.

See also [`ckcov`](@ref), [`ckobj!`](@ref).
"""
ckobj(ck::AbstractString) =
    _ck_open(daf -> sort!(unique!(Int[s.inst for s in daf.cksegments])), ck, "ckobj")

"""
    ckobj!(ids, ck) -> ids

Insert the CK IDs that have a segment in the CK file `ck` into the vector `ids`, as
`ckobj_c` does with a SPICE cell.

`ids` is used as a set: on return it is sorted and without duplicates. The vector
comes first, as in SPICE.jl; CSPICE has the cell last. A call that throws leaves `ids`
unchanged. Throws the errors of [`ckobj`](@ref).
"""
ckobj!(ids::AbstractVector{<:Integer}, ck::AbstractString) = _add_ids!(ids, ckobj(ck))

"""
    _ck_intervals!(out, data, seg) -> out

Append to `out` the intervals `(begin, finish)` in ticks on which the CK segment `seg`
can give a pointing with tolerance zero, as `zzckcv01_` to `zzckcv06_`: the single
pointing instances of type 1, the stored intervals of types 2 and 4, the
interpolation intervals of types 3 and 5 from their start to the last instance in
them, and the mini-segment intervals of type 6 up to their last packet.
"""
function _ck_intervals!(out::Vector{NTuple{2,Float64}}, data::Vector{Float64}, seg::CKSegment)
    beg, en = seg.baddr, seg.eaddr
    dtype = seg.dtype
    if dtype == 1
        nrec = trunc(Int, data[en])
        psiz = seg.avflag ? 7 : 4
        tbase = beg + nrec * psiz - 1
        for i in 1:nrec
            push!(out, (data[tbase+i], data[tbase+i]))
        end
    elseif dtype == 2
        nrec = round(Int, ((en - beg + 1) * 100.0 + 1.0) / 1001.0, RoundNearestTiesAway)
        for i in 1:nrec
            push!(out, (data[beg+8*nrec+i-1], data[beg+9*nrec+i-1]))
        end
    elseif dtype == 3 || dtype == 5
        if dtype == 3
            invls, nrec = round(Int, data[en-1]), round(Int, data[en])
            rsize = seg.avflag ? 7 : 4
        else
            subtyp = round(Int, data[en-3])
            0 <= subtyp <= 3 || throw(SpiceError(:notsupported,
                "unexpected subtype $subtyp of a CK type 5 segment"))
            invls, nrec = round(Int, data[en-1]), round(Int, data[en])
            rsize = _CK5_PACKSZ[subtyp+1]
        end
        ndir = (nrec - 1) ÷ 100
        tickat = beg + rsize * nrec
        lsttik = tickat + nrec - 1
        intat = beg + rsize * nrec + nrec + ndir
        lstint = intat + invls - 1
        start = data[intat]
        tick = data[tickat]
        while tick < start && tickat < lsttik
            tickat += 1
            tick = data[tickat]
        end
        tick < start && return out
        bail = false
        while intat <= lstint && tickat <= lsttik && !bail
            # invariant: `tick` is the first instance at or after the interval start `start`
            first = start
            finish = tick
            if intat == lstint
                finish = data[lsttik]
                bail = true
            else
                intat += 1
                start = data[intat]
                while tick < start && tickat < lsttik
                    finish = tick
                    tickat += 1
                    tick = data[tickat]
                end
                if tick < start
                    finish = tick
                    tickat += 1
                end
            end
            if dtype == 5
                first = max(first, seg.t0)
                finish = min(finish, seg.t1)
            end
            first <= finish && push!(out, (first, finish))
        end
    elseif dtype == 4
        _, _, pdrbas, pdrtyp, pktbas, npkt, pktsz, pktoff = _ck_meta04(data, seg)
        for k in 1:npkt
            a = pdrtyp == 0 ? pktbas + (k - 1) * (pktsz + pktoff) + 1 + pktoff :
                pktbas + trunc(Int, data[pdrbas+k])
            push!(out, (data[a] - data[a+1], data[a] + data[a+1]))
        end
    elseif dtype == 6
        n = round(Int, data[en])
        ptrbas = en - 2 - (n + 1)
        ivlbas = ptrbas - n ÷ 100 - (n + 1)
        for i in 1:n
            minie = beg - 1 + round(Int, data[ptrbas+i+1]) - 1
            npkt = round(Int, data[minie])
            lstepc = data[minie-4-(npkt-1)÷100]
            first = max(data[ivlbas+i], seg.t0)
            finish = min(lstepc, data[ivlbas+i+1], seg.t1)
            first <= finish && push!(out, (first, finish))
        end
    else
        throw(SpiceError(:ckunknowndatatype, "CK data type $dtype is not supported"))
    end
    return out
end

"""
    ckcov!(cover, [ctx,] ck, idcode, needav, level, tol, timsys) -> cover

Insert the coverage of the instrument or structure `idcode` in the CK file `ck` into the
window `cover`, as `ckcov_c`.

The result is the union of the old content of `cover` with the new intervals. The
window comes first, as in [`spkcov!`](@ref) and [`pckcov!`](@ref); CSPICE and SPICE.jl
have it last. See [`ckcov`](@ref) for the arguments, the intervals and the errors. A
call that throws for an argument or for the file leaves `cover` unchanged.
"""
function ckcov!(cover::Window, ks::KernelSet, ck::AbstractString, idcode::Integer,
        needav::Bool, level::AbstractString, tol::Real, timsys::AbstractString)
    tolf = _ck_check_tol(tol)
    lev = uppercase(strip(level))
    lev in ("SEGMENT", "INTERVAL") || throw(SpiceError(:invalidoption,
        "allowed values of the level are SEGMENT and INTERVAL; the value was '$level'"))
    sys = uppercase(strip(timsys))
    sys in ("SCLK", "TDB") || throw(SpiceError(:notsupported,
        "the time system was '$timsys'; allowed values are SCLK and TDB"))
    raw = NTuple{2,Float64}[]
    _ck_open(ck, "ckcov") do daf
        for seg in daf.cksegments
            (seg.inst == idcode && (seg.avflag || !needav)) || continue
            if lev == "SEGMENT"
                push!(raw, (seg.t0, seg.t1))
            else
                _ck_intervals!(raw, daf.data, seg)
            end
        end
    end
    clkid = _ck_sclkid(ks, Int32(idcode))
    # invariant: `cover` is the union of its old content and the intervals before (a, b)
    for (a, b) in raw
        if tolf > 0 || lev == "SEGMENT"
            a = max(a - tolf, 0.0)
            b += tolf
        end
        if sys == "TDB"
            a = sct2e(ks, clkid, a)
            b = sct2e(ks, clkid, b)
        end
        a <= b && wninsd!(cover, a, b)
    end
    return cover
end

"""
    ckcov([ctx,] ck, idcode, needav, level, tol, timsys) -> Window

Return the coverage [`Window`](@ref) of the instrument or structure `idcode` in the CK
file `ck`, as `ckcov_c`.

The window holds the intervals on which the file gives a pointing, as sorted disjoint
intervals; it is empty if the file has no data for `idcode`. `wnfetd(w, i)` gives
interval `i` as `(start, stop)`.

- `needav`: `true` restricts the answer to segments with angular velocity.
- `level`: `"SEGMENT"` gives the bounds in the segment descriptors; `"INTERVAL"` gives
  the intervals on which a pointing with tolerance zero exists.
- `tol`: `>= 0`, in ticks. Each interval is widened by `tol` on both sides; the start
  does not go below zero.
- `timsys`: `"SCLK"` gives encoded clock ticks; `"TDB"` gives seconds past J2000 TDB.
  `"TDB"` needs the SCLK kernel of the spacecraft in `ctx`.

The CK file is read from disk and need not be loaded.

Throws `SpiceError(:valueoutofrange)` if `tol < 0`, `SpiceError(:invalidoption)` if
`level` is not one of the two, `SpiceError(:notsupported)` if `timsys` is not one of
the two, the file errors of [`ckobj`](@ref), and for `"TDB"` the errors of
[`sct2e`](@ref).

# Examples
```julia
furnsh("naif0012.tls"); furnsh("mission.tsc")
for id in ckobj("mission.bc")
    w = ckcov("mission.bc", id, false, "INTERVAL", 0.0, "TDB")
    for i in 1:wncard(w)
        start, stop = wnfetd(w, i)
        println(id, ": ", et2utc(start, "ISOC", 0), " to ", et2utc(stop, "ISOC", 0))
    end
end
```

See also [`ckcov!`](@ref), [`spkcov`](@ref), [`pckcov`](@ref).
"""
ckcov(ks::KernelSet, ck::AbstractString, idcode::Integer, needav::Bool,
    level::AbstractString, tol::Real, timsys::AbstractString) =
    ckcov!(Window(), ks, ck, idcode, needav, level, tol, timsys)

# --- forms on a context and on the default context ---------------------------------------

ckgp(ctx::SpiceContext, inst::Integer, sclkdp::Real, tol::Real, ref::NameOrID) =
    ckgp(kernelset(ctx), inst, sclkdp, tol, ref)
ckgp(inst::Integer, sclkdp::Real, tol::Real, ref::NameOrID) =
    ckgp(default_context(), inst, sclkdp, tol, ref)
ckgpav(ctx::SpiceContext, inst::Integer, sclkdp::Real, tol::Real, ref::NameOrID) =
    ckgpav(kernelset(ctx), inst, sclkdp, tol, ref)
ckgpav(inst::Integer, sclkdp::Real, tol::Real, ref::NameOrID) =
    ckgpav(default_context(), inst, sclkdp, tol, ref)
ckmeta(ctx::SpiceContext, ckid::Integer, meta::AbstractString) =
    ckmeta(kernelset(ctx), ckid, meta)
ckmeta(ckid::Integer, meta::AbstractString) = ckmeta(default_context(), ckid, meta)
ckcov(ctx::SpiceContext, ck::AbstractString, idcode::Integer, needav::Bool,
    level::AbstractString, tol::Real, timsys::AbstractString) =
    ckcov(kernelset(ctx), ck, idcode, needav, level, tol, timsys)
ckcov(ck::AbstractString, idcode::Integer, needav::Bool, level::AbstractString, tol::Real,
    timsys::AbstractString) = ckcov(default_context(), ck, idcode, needav, level, tol, timsys)
ckcov!(cover::Window, ctx::SpiceContext, ck::AbstractString, idcode::Integer, needav::Bool,
    level::AbstractString, tol::Real, timsys::AbstractString) =
    ckcov!(cover, kernelset(ctx), ck, idcode, needav, level, tol, timsys)
ckcov!(cover::Window, ck::AbstractString, idcode::Integer, needav::Bool,
    level::AbstractString, tol::Real, timsys::AbstractString) =
    ckcov!(cover, default_context(), ck, idcode, needav, level, tol, timsys)
ckfrot(ctx::SpiceContext, inst::Integer, et::Real) = ckfrot(kernelset(ctx), inst, et)
ckfrot(inst::Integer, et::Real) = ckfrot(default_context(), inst, et)
ckfxfm(ctx::SpiceContext, inst::Integer, et::Real) = ckfxfm(kernelset(ctx), inst, et)
ckfxfm(inst::Integer, et::Real) = ckfxfm(default_context(), inst, et)
