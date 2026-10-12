# CK interface: pointing queries, CK meta data, coverage of a CK file, and the rotation
# jets that the frame and state interfaces use when a chain holds a CK frame.

# --- rotation jets with CK frames ---------------------------------------------------------

"""
    frame_jet(ks, from::Int32, to::Int32, et::Float64, Val(1)) -> (R, dR)

Order 1 rotation jet at a plain epoch, the one behind `sxform`, `spkezr` and `spkgeo`.
If CK kernels are loaded and one of the two frame chains holds a CK frame, `dR` comes
from the angular velocity of the CK segment, as `sxform_c` of CSPICE gives it
(`_ck_state_rotation`); otherwise the result is `rotation_jet`, unchanged.
"""
@inline function frame_jet(ks::KernelSet, from::Int32, to::Int32, et::Float64, v::Val{1})
    (isempty(ks.ck) || !(_ck_in_chain(ks, from) || _ck_in_chain(ks, to))) &&
        return rotation_jet(ks, from, to, et, v)
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
    ckgp(inst, sclkdp, tol, ref) -> (cmat::SMatrix{3,3}, clkout::Float64, found::Bool)

Pointing of the instrument or structure `inst` (CK ID) at the encoded spacecraft clock
`sclkdp` within `tol` ticks, as `ckgp_c`: `cmat` rotates vectors from the frame `ref`
(name or ID) to the instrument frame, and `clkout` is the clock time of the pointing
that was used. Segments with and without angular velocity are searched in priority
order. If no segment gives a pointing, `found` is false and `cmat` is zero. If `ref`
is not the frame of the segment, the matrix is changed with the rotation between the
two frames at the ET of `clkout` (which needs the SCLK kernel of the instrument unless
both frames are inertial). Throws `SpiceError(:valueoutofrange)` for `tol < 0`.
Forms: on a `KernelSet`, on a `SpiceContext`, and on the default context.
"""
function ckgp(ks::KernelSet, inst::Integer, sclkdp::Real, tol::Real, ref::NameOrID)
    i = Int32(inst)
    found, cmat, _, clkout, refseg = _ck_search(ks, i, Float64(sclkdp), _ck_check_tol(tol), false)
    found || return zero(Mat3{Float64}), 0.0, false
    refreq = frame_id(ks, ref)
    if refreq != refseg
        et = _ck_change_epoch(ks, i, refreq, refseg, clkout)
        cmat = cmat * frame_rotation(ks, refreq, refseg, et)
    end
    return cmat, clkout, true
end

"""
    ckgpav(inst, sclkdp, tol, ref) -> (cmat, av::SVector{3}, clkout, found)

Pointing and angular velocity of the instrument `inst` at the encoded spacecraft clock
`sclkdp` within `tol` ticks, as `ckgpav_c`. Only segments with angular velocity are
used. `av` is the angular velocity of the instrument frame relative to `ref`, in
radians per second, expressed in `ref`. For the rest see `ckgp`.
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
        omega = transpose(drdt) * rot           # `xf2rav_`
        cmat = cmat * rot
        av = Vec3{Float64}(omega[3, 2], omega[1, 3], omega[2, 1]) + transpose(rot) * av
    end
    return cmat, av, clkout, true
end

"""
    ckmeta(ckid, meta) -> Int

ID that belongs to the CK ID `ckid`, as `ckmeta_c`: for `meta = "SCLK"` the spacecraft
clock (`CK_<ckid>_SCLK`), for `meta = "SPK"` the body (`CK_<ckid>_SPK`); without the
variable `ckid ÷ 1000` if `ckid <= -1000`, else 0. Throws `SpiceError(:unknownckmeta)`
for another item.
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

CK IDs of the instruments and structures that have a segment in the CK file at the path
`ck`, ascending and without repeats, as the set that `ckobj_c` returns. The file does
not have to be loaded.
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
    ckcov(ck, idcode, needav, level, tol, timsys) -> Vector{NTuple{2,Float64}}

Coverage of the instrument `idcode` in the CK file at the path `ck`, as `ckcov_c`: the
ascending disjoint intervals `(start, stop)` on which the file gives a pointing.
`needav` restricts the answer to segments with angular velocity. `level` is
`"SEGMENT"` (the bounds in the segment descriptors) or `"INTERVAL"` (the intervals on
which a pointing with tolerance zero exists). Each interval is widened by `tol >= 0`
ticks on both sides (the start not below zero). `timsys` is `"SCLK"` (ticks) or `"TDB"`
(ET seconds, converted with the SCLK kernel of the instrument, which must be loaded in
the `KernelSet`). The CK file does not have to be loaded. Forms: on a `KernelSet`, on a
`SpiceContext`, and on the default context.
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
