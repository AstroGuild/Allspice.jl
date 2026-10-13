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
    ckgp([ctx,] inst, sclkdp, tol, ref) -> (cmat, clkout, found)

Return the pointing of the instrument or structure `inst` at the encoded spacecraft
clock value `sclkdp`, as `ckgp_c`.

- `inst`: CK ID of the instrument or structure.
- `sclkdp`: encoded spacecraft clock, in ticks (see [`sce2c`](@ref)).
- `tol`: time tolerance in ticks, `>= 0`. A pointing within `tol` of `sclkdp` is
  accepted.
- `ref`: name or ID of the reference frame of the result.

`cmat` is an `SMatrix{3,3}` that rotates a vector from `ref` to the instrument frame.
`clkout` is the clock value, in ticks, of the pointing that was used. If no loaded CK
segment gives a pointing, `found` is `false`, `cmat` is zero and `clkout` is 0.
Segments with and without angular velocity are searched in priority order.

If `ref` is not the frame of the segment, the matrix is changed with the rotation
between the two frames at the epoch of `clkout`. This needs the SCLK kernel of the
spacecraft, unless both frames are inertial.

Throws `SpiceError(:valueoutofrange)` if `tol < 0`. Throws
`SpiceError(:unknownframe)` if `ref` is not a known frame name; `ref` is read only
when a pointing is found.

# Examples
```julia
furnsh("naif0012.tls"); furnsh("mission.tsc"); furnsh("mission.bc")
sclkdp = sce2c(-77, et)
cmat, clkout, found = ckgp(-77000, sclkdp, 0.0, "J2000")
found && (boresight_j2000 = cmat' * [0.0, 0.0, 1.0])
```

See also [`ckgpav`](@ref), [`ckcov`](@ref), [`pxform`](@ref).
"""
function ckgp(ks::KernelSet, inst::Integer, sclkdp::Real, tol::Real, ref::NameOrID)
    i = Int32(inst)
    found, cmat, _, clkout, refseg = _ck_search(ks, i, Float64(sclkdp), _ck_check_tol(tol), false)
    found || return zero(Mat3{Float64}), 0.0, false
    refreq = frame_id(ks, ref)
    if refreq != refseg
        et = _ck_change_epoch(ks, i, refreq, refseg, clkout)
        cmat = _m_mxm(cmat, frame_rotation(ks, refreq, refseg, et))
    end
    return cmat, clkout, true
end

"""
    ckgpav([ctx,] inst, sclkdp, tol, ref) -> (cmat, av, clkout, found)

Return the pointing and the angular velocity of the instrument or structure `inst` at
the encoded spacecraft clock value `sclkdp`, as `ckgpav_c`.

`av` is an `SVector{3}`: the angular velocity of the instrument frame relative to
`ref`, in rad/s, with components in `ref`. Only segments with angular velocity are
used. If none gives a pointing, `found` is `false` and `cmat`, `av` and `clkout` are
zero. The other arguments, results and errors are those of [`ckgp`](@ref).
"""
function ckgpav(ks::KernelSet, inst::Integer, sclkdp::Real, tol::Real, ref::NameOrID)
    i = Int32(inst)
    found, cmat, av, clkout, refseg = _ck_search(ks, i, Float64(sclkdp), _ck_check_tol(tol), true)
    found || return zero(Mat3{Float64}), zero(Vec3{Float64}), 0.0, false
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
    return cmat, av, clkout, true
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

"Open `path` as a CK file; throws `SpiceError(:invalidfiletype)` for another DAF kernel."
function _ck_open(path::AbstractString)
    daf = open_daf(path)
    daf.kind === :CK || throw(SpiceError(:invalidfiletype,
        "the file $path is a binary $(daf.kind) kernel; a CK file is needed"))
    return daf
end

"""
    ckobj(ck) -> Vector{Int}

Return the CK IDs of the instruments and structures that have a segment in the CK file
`ck`, sorted and without duplicates, as `ckobj_c`.

The file is read from disk and need not be loaded.

Throws `SpiceError(:nofile)` if the file does not exist, `SpiceError(:baddaf)` if it
cannot be read as a DAF, and `SpiceError(:invalidfiletype)` if it is an SPK or a
binary PCK.

See also [`ckcov`](@ref).
"""
ckobj(ck::AbstractString) = sort!(unique!(Int[s.inst for s in _ck_open(ck).cksegments]))

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
            0 <= subtyp <= 3 || throw(SpiceError(:unsupported,
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
        throw(SpiceError(:unsupported, "CK data type $dtype is not supported"))
    end
    return out
end

"""
    _wnmerge(v) -> Vector{NTuple{2,Float64}}

The window (ascending disjoint intervals) that is the union of the intervals `v`, as
repeated `wninsd_`: intervals that overlap or touch are joined.
"""
function _wnmerge(v::Vector{NTuple{2,Float64}})
    out = NTuple{2,Float64}[]
    for (a, b) in sort(v)       # invariant: `out` is the window of the intervals before (a, b)
        if !isempty(out) && a <= out[end][2]
            out[end] = (out[end][1], max(out[end][2], b))
        else
            push!(out, (a, b))
        end
    end
    return out
end

"""
    ckcov([ctx,] ck, idcode, needav, level, tol, timsys) -> Vector{NTuple{2,Float64}}

Return the coverage of the instrument or structure `idcode` in the CK file `ck`, as
`ckcov_c`.

The result is the list of the intervals `(start, stop)` on which the file gives a
pointing, ascending and disjoint. Different from [`spkcov`](@ref), the result is a
vector of tuples and not a `Window`, and there is no form that adds to a window.

- `needav`: `true` restricts the answer to segments with angular velocity.
- `level`: `"SEGMENT"` gives the bounds in the segment descriptors; `"INTERVAL"` gives
  the intervals on which a pointing with tolerance zero exists.
- `tol`: `>= 0`, in ticks. Each interval is widened by `tol` on both sides; the start
  does not go below zero.
- `timsys`: `"SCLK"` gives encoded clock ticks; `"TDB"` gives seconds past J2000 TDB.
  `"TDB"` needs the SCLK kernel of the spacecraft in `ctx`.

The CK file is read from disk and need not be loaded.

Throws `SpiceError(:valueoutofrange)` if `tol < 0`, `SpiceError(:invalidoption)` if
`level` is not one of the two, `SpiceError(:unsupported)` if `timsys` is not one of
the two, the file errors of [`ckobj`](@ref), and for `"TDB"` the errors of
[`sct2e`](@ref).

# Examples
```julia
furnsh("naif0012.tls"); furnsh("mission.tsc")
for id in ckobj("mission.bc")
    for (start, stop) in ckcov("mission.bc", id, false, "INTERVAL", 0.0, "TDB")
        println(id, ": ", et2utc(start, "ISOC", 0), " to ", et2utc(stop, "ISOC", 0))
    end
end
```
"""
function ckcov(ks::KernelSet, ck::AbstractString, idcode::Integer, needav::Bool,
        level::AbstractString, tol::Real, timsys::AbstractString)
    tolf = _ck_check_tol(tol)
    lev = uppercase(strip(level))
    lev in ("SEGMENT", "INTERVAL") || throw(SpiceError(:invalidoption,
        "allowed values of the level are SEGMENT and INTERVAL; the value was '$level'"))
    sys = uppercase(strip(timsys))
    sys in ("SCLK", "TDB") || throw(SpiceError(:unsupported,
        "the time system was '$timsys'; allowed values are SCLK and TDB"))
    daf = _ck_open(ck)
    raw = NTuple{2,Float64}[]
    for seg in daf.cksegments
        (seg.inst == idcode && (seg.avflag || !needav)) || continue
        if lev == "SEGMENT"
            push!(raw, (seg.t0, seg.t1))
        else
            _ck_intervals!(raw, daf.data, seg)
        end
    end
    clkid = _ck_sclkid(ks, Int32(idcode))
    out = NTuple{2,Float64}[]
    for (a, b) in raw
        if tolf > 0 || lev == "SEGMENT"
            a = max(a - tolf, 0.0)
            b += tolf
        end
        if sys == "TDB"
            a = sct2e(ks, clkid, a)
            b = sct2e(ks, clkid, b)
        end
        a <= b && push!(out, (a, b))
    end
    return _wnmerge(out)
end

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
