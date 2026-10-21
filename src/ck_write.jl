# CK writers: `ckopn`, `ckcls` and the segment writers `ckw01`, `ckw02`, `ckw03`,
# `ckw05`, `ckw06`. Each writer makes the input checks of the SPICELIB routine in the
# same order and writes the descriptor, the pointing data, the directories and the
# trailer in the layout of CSPICE, so the data area is equal bit for bit to the one
# that CSPICE writes from the same input.
#
# Argument lists follow SPICE.jl, except that the record count is always the length
# of the input (SPICE.jl 0.2.5 passes one record less to `ckw02_c` and `ckw03_c`).
# Quaternions are SPICE quaternions (scalar first). A table of quaternions, angular
# velocities or packets is a vector of vectors, a matrix with one row of the table in
# each column, or a flat vector. Times are encoded spacecraft clock ticks. See
# `spk_write.jl` for the frame argument, the keyword `ctx` and the error kinds.

"""
    ckopn(path, ifname="CK_file", ncomch=0) -> DAFHandle

Create the CK file `path` and open it for writing, as `ckopn_c`.

`ifname` is the internal file name and `ncomch` is the number of comment characters to
make room for. The file is on disk when the call returns. An existing file is never
replaced. Add segments with the `ckwNN` writers and close the handle with
[`ckcls`](@ref). Load the file with `furnsh` only after the close.

Throws `SpiceError(:emptystring)` if `path` or `ifname` is empty, and `:fileopenfailed` if
the file exists or cannot be created.

# Examples
Write a type 3 segment with two pointing records and one interpolation interval:
```jldoctest
julia> path = joinpath(mktempdir(), "example.bc");

julia> h = ckopn(path, "example CK");

julia> quats = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]];

julia> ckw03(h, 0.0, 10.0, -77000, "J2000", "attitude", [0.0, 10.0], quats, [0.0])

julia> ckcls(h)

julia> getfat(path)
("DAF", "CK")
```

See also [`ckw01`](@ref), [`ckw02`](@ref), [`ckw03`](@ref), [`ckw05`](@ref),
[`ckw06`](@ref), [`dafac`](@ref) to add comments.
"""
function ckopn(path::AbstractString, ifname::AbstractString="CK_file", ncomch::Integer=0)
    _seg_chkfstr("ckopn", "path", path)
    _seg_chkfstr("ckopn", "ifname", ifname)
    ncomr = ncomch > 0 ? (Int(ncomch) - 1) ÷ _DAF_CRECLEN + 1 : 0
    return dafonw(path, "CK", 2, 6, ifname, ncomr)
end

"""
    ckcls(handle) -> Nothing

Close the CK file of `handle`, as `ckcls_c`. The file is then complete on disk and can be
loaded with `furnsh`. The call does nothing if the handle is closed already.

Throws `SpiceError(:nosegmentsfound)` if the file is open for writing and holds no
segment. The handle then stays open; close it with [`dafcls`](@ref).
"""
ckcls(h::DAFHandle) = _seg_close(h, "CK")

function _ck_add!(h::DAFHandle, inst::Integer, refcod::Int, dtype::Int, avflag::Bool,
        begtim::Real, endtim::Real, segid::AbstractString, data::Vector{Float64})
    _daf_add_array!(h, Float64[begtim, endtim], Int[inst, refcod, dtype, avflag ? 1 : 0, 0, 0],
        segid, data)
    return nothing
end

"""
    _ck_avvs(who, avvs, nrec) -> (flat::Vector{Float64}, avflag::Bool)

Angular velocities for `nrec` records: the segment has angular velocity if `avvs` holds
`nrec` vectors; `nothing` (the default argument) means none, and so does one vector for
more than one record (the default argument `[zeros(3)]` of SPICE.jl). Another count
throws `SpiceError(:badarg)`.
"""
function _ck_avvs(who::AbstractString, avvs, nrec::Int)
    avvs === nothing && return Float64[], false
    av, n, _ = _seg_flat(who, "avvs", avvs, 3)
    n == nrec && return av, true
    n == 1 || throw(SpiceError(:badarg,
        "$who: `avvs` must have $nrec elements (or be omitted) but has $n"))
    return Float64[], false
end

"Throws `SpiceError(:zeroquaternion)` if one of the `n` quaternions that start every `step` doubles in `flat` is zero."
function _ck_nonzero_quats(who::AbstractString, flat::Vector{Float64}, n::Int, step::Int)
    for i in 1:n
        base = (i - 1) * step
        (flat[base+1] == 0.0 && flat[base+2] == 0.0 && flat[base+3] == 0.0 &&
            flat[base+4] == 0.0) && throw(SpiceError(:zeroquaternion,
            "$who: the quaternion at index $i has magnitude zero"))
    end
    return nothing
end

function _ck_increasing(who::AbstractString, what::AbstractString, t::Vector{Float64})
    for i in 2:length(t)
        t[i] <= t[i-1] && throw(SpiceError(:timesoutoforder, "$who: the $what times are " *
            "not strictly increasing: $what($(i - 1)) = $(t[i-1]) and $what($i) = $(t[i])"))
    end
    return nothing
end

"Quaternions and, if `avflag`, angular velocities of `nrec` records interleaved as CK types 1 and 3 store them."
function _ck_pointing(q::Vector{Float64}, av::Vector{Float64}, nrec::Int, avflag::Bool)
    avflag || return copy(q)
    data = Vector{Float64}(undef, 7nrec)
    for i in 1:nrec                      # invariant: records 1 .. i-1 are in data
        copyto!(data, 7(i - 1) + 1, q, 4(i - 1) + 1, 4)
        copyto!(data, 7(i - 1) + 5, av, 3(i - 1) + 1, 3)
    end
    return data
end

"""
    ckw01(handle, begtim, endtim, inst, ref, segid, sclkdp, quats, avvs=nothing;
          ctx=default_context()) -> Nothing

Write a CK segment of type 1 (discrete pointing) to the file of `handle`, as `ckw01_c`.

# Arguments
The arguments `handle`, `begtim`, `endtim`, `inst`, `ref`, `segid` and `ctx` have the same
meaning in every `ckwNN` writer. [`ckw06`](@ref) names the coverage `first`, `last`.
- `handle`: a [`DAFHandle`](@ref) from [`ckopn`](@ref).
- `begtim`, `endtim`: the coverage of the segment, encoded spacecraft clock ticks.
- `inst`: the NAIF ID of the instrument or structure.
- `ref`: the reference frame, a name or a frame ID. A name is resolved in `ctx` (a
  `SpiceContext` or a `KernelSet`); built-in frames need no kernel.
- `segid`: the segment identifier, at most 40 printable characters.
- `sclkdp`: the times of the pointing records, encoded spacecraft clock ticks.
- `quats`: one SPICE quaternion (scalar first) for each record.
- `avvs`: one angular velocity vector (rad/s) for each record, or `nothing`.

The number of records is `length(sclkdp)`. A table such as `quats` is a vector of vectors
(the form of SPICE.jl), a matrix with one record in each column, or a flat vector.

The segment has angular velocity if `avvs` holds one vector for each record. It has none
if `avvs` is `nothing`, which is the default. One vector for more than one record also
gives none: that is the default argument `[zeros(3)]` of SPICE.jl.

The segment is on disk when the call returns. It holds the same data as the segment that
`ckw01_c` writes. A call that throws writes nothing.

Throws `SpiceError`. Every `ckwNN` writer throws:
- `:dafnosuchhandle` if the handle is closed, `:dafinvalidaccess` if it is open for
  reading;
- `:emptystring` if `ref` or `segid` is empty;
- `:badarg` if the shape of an input does not fit;
- `:invalidrefframe` if the frame is not known;
- `:nonprintablechars` or `:segidtoolong` if `segid` is not valid.

`ckw01` also throws:
- `:invalidnumrec` if there is no record;
- `:invaliddescrtime` unless `begtim <= sclkdp[1]` and `endtim >= sclkdp[end]`;
- `:invalidsclktime` if `sclkdp[1] < 0`;
- `:timesoutoforder` unless `sclkdp` is strictly increasing;
- `:zeroquaternion` if a quaternion is zero.

See also [`ckw02`](@ref), [`ckw03`](@ref), [`ckw05`](@ref), [`ckw06`](@ref), and
[`ckopn`](@ref) for an example.
"""
function ckw01(h::DAFHandle, begtim::Real, endtim::Real, inst::Integer, ref::SegFrame,
        segid::AbstractString, sclkdp, quats, avvs=nothing; ctx::SegContext=default_context())
    who = "ckw01"
    _daf_check_write(h, who)
    ref isa AbstractString && _seg_chkfstr(who, "ref", ref)
    _seg_chkfstr(who, "segid", segid)
    t = _seg_vec(who, "sclkdp", sclkdp, -1)
    nrec = length(t)
    q, nq, _ = _seg_flat(who, "quats", quats, 4)
    nq == nrec || throw(SpiceError(:badarg, "$who: `quats` must have $nrec elements but has $nq"))
    av, avflag = _ck_avvs(who, avvs, nrec)
    nrec <= 0 && throw(SpiceError(:invalidnumrec,
        "$who: $nrec is an invalid number of pointing instances for type 1"))
    begtim > t[1] && throw(SpiceError(:invaliddescrtime, "$who: the first d.p. component " *
        "of the descriptor is invalid: DCD(1) = $begtim and SCLKDP(1) = $(t[1])"))
    endtim < t[nrec] && throw(SpiceError(:invaliddescrtime, "$who: the second d.p. component " *
        "of the descriptor is invalid: DCD(2) = $endtim and SCLKDP(NREC) = $(t[nrec])"))
    refcod = _seg_refcode(who, ctx, ref)
    _seg_idchars(who, segid)
    _seg_idlen(who, segid)
    t[1] < 0.0 && throw(SpiceError(:invalidsclktime, "$who: the first SCLKDP time $(t[1]) is negative"))
    _ck_increasing(who, "SCLKDP", t)
    _ck_nonzero_quats(who, q, nrec, 4)
    data = _ck_pointing(q, av, nrec, avflag)
    append!(data, t)
    for i in 1:(nrec-1)÷100              # directory: midpoints between records 100i and 100i+1
        push!(data, (t[100i] + t[100i+1]) / 2.0)
    end
    push!(data, Float64(nrec))
    _ck_add!(h, inst, refcod, 1, avflag, begtim, endtim, segid, data)
    return nothing
end

"""
    ckw02(handle, begtim, endtim, inst, ref, segid, start, stop, quats, avvs, rates;
          ctx=default_context()) -> Nothing

Write a CK segment of type 2 (pointing with a constant angular velocity on intervals), as
`ckw02_c`.

Record `i` is valid from `start[i]` to `stop[i]` (encoded spacecraft clock ticks).
`quats[i]` is its quaternion, `avvs[i]` its angular velocity (rad/s) and `rates[i]` its
number of seconds per tick. The number of records is `length(start)`. SPICE.jl 0.2.5
writes one record less. See [`ckw01`](@ref) for the other arguments and for
the errors of every writer.

Throws `SpiceError`:
- `:invalidnumrec` if there is no record;
- `:invaliddescrtime` unless `begtim <= start[1]` and `endtim >= stop[end]`;
- `:invalidsclktime` if `start[1] < 0`;
- `:degenerateinterval` unless `stop[i] > start[i]`;
- `:timesoutoforder` unless `start` is strictly increasing;
- `:badstoptime` if `stop[i-1] > start[i]`;
- `:zeroquaternion` if a quaternion is zero.
"""
function ckw02(h::DAFHandle, begtim::Real, endtim::Real, inst::Integer, ref::SegFrame,
        segid::AbstractString, start, stop, quats, avvs, rates; ctx::SegContext=default_context())
    who = "ckw02"
    _daf_check_write(h, who)
    ref isa AbstractString && _seg_chkfstr(who, "ref", ref)
    _seg_chkfstr(who, "segid", segid)
    t0 = _seg_vec(who, "start", start, -1)
    nrec = length(t0)
    t1 = _seg_vec(who, "stop", stop, nrec)
    rt = _seg_vec(who, "rates", rates, nrec)
    q, nq, _ = _seg_flat(who, "quats", quats, 4)
    av, na, _ = _seg_flat(who, "avvs", avvs, 3)
    (nq == nrec && na == nrec) || throw(SpiceError(:badarg,
        "$who: `quats` and `avvs` must have $nrec elements but have $nq and $na"))
    nrec <= 0 && throw(SpiceError(:invalidnumrec,
        "$who: $nrec is an invalid number of pointing instances for type 2"))
    begtim > t0[1] && throw(SpiceError(:invaliddescrtime, "$who: the first d.p. component " *
        "of the descriptor is invalid: DCD(1) = $begtim and START(1) = $(t0[1])"))
    endtim < t1[nrec] && throw(SpiceError(:invaliddescrtime, "$who: the second d.p. component " *
        "of the descriptor is invalid: DCD(2) = $endtim and STOP(NREC) = $(t1[nrec])"))
    refcod = _seg_refcode(who, ctx, ref)
    _seg_idchars(who, segid)
    _seg_idlen(who, segid)
    t0[1] < 0.0 && throw(SpiceError(:invalidsclktime, "$who: the first START time $(t0[1]) is negative"))
    for i in 1:nrec                      # invariant: intervals 1 .. i-1 are valid and ordered
        t1[i] <= t0[i] && throw(SpiceError(:degenerateinterval, "$who: the STOP time is " *
            "less than or equal to the START time for interval number $i"))
        i == 1 && continue
        t0[i] <= t0[i-1] && throw(SpiceError(:timesoutoforder, "$who: the START times are " *
            "not strictly increasing: START($(i - 1)) = $(t0[i-1]) and START($i) = $(t0[i])"))
        t1[i-1] > t0[i] && throw(SpiceError(:badstoptime, "$who: the STOP time for " *
            "interval $(i - 1) is greater than the following START time"))
    end
    _ck_nonzero_quats(who, q, nrec, 4)
    data = Vector{Float64}(undef, 8nrec)
    for i in 1:nrec                      # invariant: records 1 .. i-1 are in data
        copyto!(data, 8(i - 1) + 1, q, 4(i - 1) + 1, 4)
        copyto!(data, 8(i - 1) + 5, av, 3(i - 1) + 1, 3)
        data[8i] = rt[i]
    end
    append!(data, t0)
    append!(data, t1)
    for i in 1:(nrec-1)÷100
        push!(data, (t1[100i] + t0[100i+1]) / 2.0)
    end
    _ck_add!(h, inst, refcod, 2, true, begtim, endtim, segid, data)
    return nothing
end

"""
    ckw03(handle, begtim, endtim, inst, ref, segid, sclkdp, quats, starts, avvs=nothing;
          ctx=default_context()) -> Nothing

Write a CK segment of type 3 (pointing with linear interpolation inside intervals), as
`ckw03_c`.

`sclkdp` holds the times of the `nrec` records (encoded spacecraft clock ticks), `quats`
their quaternions and `starts` the start times of the interpolation intervals. `avvs`
is as for [`ckw01`](@ref), which also gives the other arguments and the errors of every
writer. All `nrec` records are written. SPICE.jl 0.2.5 writes one record less.

Throws `SpiceError`:
- `:invalidnumrec` if there is no record, `:invalidnumint` if there is no interval;
- `:invaliddescrtime` unless `begtim <= sclkdp[1]` and `endtim >= sclkdp[end]`;
- `:invalidsclktime` if `sclkdp[1] < 0`;
- `:timesoutoforder` unless `sclkdp` and `starts` are strictly increasing;
- `:timesdontmatch` unless `starts[1] == sclkdp[1]`;
- `:invalidstarttime` unless every start time is equal to one of the `sclkdp`;
- `:zeroquaternion` if a quaternion is zero.

See [`ckopn`](@ref) for an example.
"""
function ckw03(h::DAFHandle, begtim::Real, endtim::Real, inst::Integer, ref::SegFrame,
        segid::AbstractString, sclkdp, quats, starts, avvs=nothing;
        ctx::SegContext=default_context())
    who = "ckw03"
    _daf_check_write(h, who)
    ref isa AbstractString && _seg_chkfstr(who, "ref", ref)
    _seg_chkfstr(who, "segid", segid)
    t = _seg_vec(who, "sclkdp", sclkdp, -1)
    nrec = length(t)
    st = _seg_vec(who, "starts", starts, -1)
    nints = length(st)
    q, nq, _ = _seg_flat(who, "quats", quats, 4)
    nq == nrec || throw(SpiceError(:badarg, "$who: `quats` must have $nrec elements but has $nq"))
    av, avflag = _ck_avvs(who, avvs, nrec)
    nrec <= 0 && throw(SpiceError(:invalidnumrec,
        "$who: $nrec is an invalid number of pointing instances for type 3"))
    nints <= 0 && throw(SpiceError(:invalidnumint,
        "$who: $nints is an invalid number of interpolation intervals for type 3"))
    begtim > t[1] && throw(SpiceError(:invaliddescrtime, "$who: the segment begin time is " *
        "greater than the time of the first pointing instance: DCD(1) = $begtim and " *
        "SCLKDP(1) = $(t[1])"))
    endtim < t[nrec] && throw(SpiceError(:invaliddescrtime, "$who: the segment end time is " *
        "less than the time of the last pointing instance: DCD(2) = $endtim and " *
        "SCLKDP($nrec) = $(t[nrec])"))
    refcod = _seg_refcode(who, ctx, ref)
    _seg_idchars(who, segid)
    _seg_idlen(who, segid)
    t[1] < 0.0 && throw(SpiceError(:invalidsclktime, "$who: the first SCLKDP time $(t[1]) is negative"))
    _ck_increasing(who, "SCLKDP", t)
    t[1] != st[1] && throw(SpiceError(:timesdontmatch, "$who: the start time of the first " *
        "interval $(st[1]) and the time of the first pointing instance $(t[1]) are not the same"))
    _ck_increasing(who, "STARTS", st)
    index = 0
    for i in 1:nints                     # invariant: starts 1 .. i-1 match records up to `index`
        match = false
        while !match && index < nrec
            index += 1
            match = st[i] == t[index]
        end
        match || throw(SpiceError(:invalidstarttime,
            "$who: interval start time number $i is invalid: STARTS($i) = $(st[i])"))
    end
    _ck_nonzero_quats(who, q, nrec, 4)
    data = _ck_pointing(q, av, nrec, avflag)
    append!(data, t)
    _seg_directory!(data, t, (nrec - 1) ÷ 100)
    append!(data, st)
    _seg_directory!(data, st, (nints - 1) ÷ 100)
    push!(data, Float64(nints), Float64(nrec))
    _ck_add!(h, inst, refcod, 3, avflag, begtim, endtim, segid, data)
    return nothing
end

"""
    ckw05(handle, subtyp, degree, begtim, endtim, inst, ref, avflag, segid, sclkdp, packts,
          rate, nints, starts; ctx=default_context()) -> Nothing

Write a CK segment of type 5, as `ckw05_c`. The type holds the interpolation methods of
the attitude files of ESOC, the operations center of ESA (format DDID).

`subtyp` selects the interpolation and the content of a packet. `q` is a SPICE
quaternion (scalar first), `dq` its derivative, `av` the angular velocity in `ref`
(rad/s) and `dav` its derivative (rad/s²):

| `subtyp` | Interpolation | Packet                 | Values |
|:---------|:--------------|:-----------------------|:-------|
| 0        | Hermite       | `q`, `dq`              | 8      |
| 1        | Lagrange      | `q`                    | 4      |
| 2        | Hermite       | `q`, `dq`, `av`, `dav` | 14     |
| 3        | Lagrange      | `q`, `av`              | 7      |

The other arguments:
- `degree`: the degree of the interpolating polynomial. The window size is the number of
  packets that one interpolation uses: `degree + 1` for Lagrange and `(degree + 1) ÷ 2`
  for Hermite.
- `avflag`: `true` if the segment has angular velocity.
- `sclkdp`: the times of the packets `packts`, encoded spacecraft clock ticks.
- `rate`: the nominal clock rate, seconds per tick.
- `nints`, `starts`: the first `nints` elements of `starts` are the start times of the
  interpolation intervals. No interpolation uses packets of two intervals.

See [`ckw01`](@ref) for the remaining arguments and for the errors of every writer.

Throws `SpiceError`:
- `:toofewpackets` if there is no packet, `:invalidnumints` unless `nints > 0`;
- `:invalidsclktime` if `sclkdp[1] < 0`;
- `:timesoutoforder` unless `sclkdp` and `starts` are strictly increasing;
- `:invalidstarttime` unless every start time is equal to one of the `sclkdp`;
- `:notsupported` if `subtyp` is not in `0:3`;
- `:zeroquaternion` if a quaternion is zero;
- `:invaliddegree` unless `1 <= degree <= 23` and the window size is even: the degree is
  3 modulo 4 for a Hermite subtype, and odd for a Lagrange subtype;
- `:baddescrtimes` if `endtim < begtim`;
- `:emptysegment` if no time of `sclkdp` is inside [`begtim`, `endtim`];
- `:invalidvalue` if `rate` is zero.
"""
function ckw05(h::DAFHandle, subtyp::Integer, degree::Integer, begtim::Real, endtim::Real,
        inst::Integer, ref::SegFrame, avflag::Bool, segid::AbstractString, sclkdp, packts,
        rate::Real, nints::Integer, starts; ctx::SegContext=default_context())
    who = "ckw05"
    _daf_check_write(h, who)
    ref isa AbstractString && _seg_chkfstr(who, "ref", ref)
    _seg_chkfstr(who, "segid", segid)
    t = _seg_vec(who, "sclkdp", sclkdp, -1)
    n = length(t)
    allst = _seg_vec(who, "starts", starts, -1)
    nints <= length(allst) || throw(SpiceError(:badarg,
        "$who: `starts` has $(length(allst)) elements, fewer than nints = $nints"))
    n < 1 && throw(SpiceError(:toofewpackets,
        "$who: at least 1 packet is required for CK type 5; number of packets supplied: $n"))
    nints <= 0 && throw(SpiceError(:invalidnumints,
        "$who: $nints is an invalid number of interpolation intervals for type 5"))
    st = allst[1:nints]
    refcod = _seg_refcode(who, ctx, ref)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    t[1] < 0.0 && throw(SpiceError(:invalidsclktime, "$who: the first SCLKDP time $(t[1]) is negative"))
    _ck_increasing(who, "SCLKDP", t)
    _ck_increasing(who, "STARTS", st)
    for i in 1:nints
        k = searchsortedfirst(t, st[i])
        (k <= n && t[k] == st[i]) || throw(SpiceError(:invalidstarttime,
            "$who: interval start time number $i is invalid: STARTS($i) = $(st[i])"))
    end
    (0 <= subtyp <= 3) || throw(SpiceError(:notsupported, "$who: CK type 5 subtype <$subtyp> is not supported"))
    packsz = (8, 4, 14, 7)[subtyp+1]
    winsiz = isodd(subtyp) ? degree + 1 : (degree + 1) ÷ 2
    pk, np, _ = _seg_flat(who, "packts", packts, packsz)
    np == n || throw(SpiceError(:badarg, "$who: `packts` must have $n elements but has $np"))
    _ck_nonzero_quats(who, pk, n, packsz)
    (degree < 1 || degree > 23) && throw(SpiceError(:invaliddegree, "$who: the " *
        "interpolating polynomials have degree $degree; the valid degree range is [1, 23]"))
    isodd(winsiz) && throw(SpiceError(:invaliddegree, "$who: the interpolating polynomials " *
        "have degree $degree; for CK type 5, the degree must be equivalent to 3 mod 4 " *
        "for Hermite interpolation and odd for Lagrange interpolation"))
    endtim < begtim && throw(SpiceError(:baddescrtimes,
        "$who: descriptor bounds are non-increasing: $begtim:$endtim"))
    i = searchsortedfirst(t, begtim) - 1     # number of time tags less than begtim
    i == n && throw(SpiceError(:emptysegment,
        "$who: all time tags are less than segment start time $begtim"))
    t[i+1] > endtim && throw(SpiceError(:emptysegment, "$who: no time tags lie between " *
        "the segment start time $begtim and segment end time $endtim"))
    rate == 0 && throw(SpiceError(:invalidvalue, "$who: the SCLK rate RATE was zero"))
    data = append!(pk, t)
    _seg_directory!(data, t, (n - 1) ÷ 100)
    append!(data, st)
    _seg_directory!(data, st, (nints - 1) ÷ 100)
    push!(data, Float64(rate), Float64(subtyp), Float64(winsiz), Float64(nints), Float64(n))
    _ck_add!(h, inst, refcod, 5, avflag, begtim, endtim, segid, data)
    return nothing
end

const _CK06_PKTSZS = (8, 4, 14, 7)

"""
    ckw06(handle, inst, ref, avflag, first, last, segid, npkts, subtps, degres, packts,
          rates, sclkdp, ivlbds, sellst; ctx=default_context()) -> Nothing

Write a CK segment of type 6 (mini-segments of the ESOC/DDID types), as `ckw06` of
SPICELIB. CSPICE has no C wrapper for it.

A mini-segment is a set of packets with its own subtype, degree and clock rate, valid on
its own interval. The number of mini-segments is `length(npkts)`. The subtypes, their
packets and the window size are those of [`ckw05`](@ref).

- `first`, `last`: the coverage of the segment, encoded spacecraft clock ticks. They
  take the place of `begtim` and `endtim`.
- `avflag`: `true` if the segment has angular velocity.
- `npkts[i]`, `subtps[i]`, `degres[i]`, `rates[i]`: the number of packets, the subtype,
  the degree and the clock rate (seconds per tick) of mini-segment `i`.
- `packts`: the flat vector of all packets, in order.
- `sclkdp`: the flat vector of the times of all packets, encoded spacecraft clock ticks.
- `ivlbds`: the `length(npkts) + 1` bounds of the intervals of the mini-segments, ticks.
- `sellst`: if `true`, a time on a bound of two intervals uses the later interval; else
  it uses the earlier one.

See [`ckw01`](@ref) for the other arguments and for the errors of every writer.

Throws `SpiceError`:
- `:baddescrtimes` if `first > last`;
- `:invalidcount` if there is no mini-segment;
- `:boundsoutoforder` unless `ivlbds` is strictly increasing;
- `:coveragegap` unless [`first`, `last`] is inside [`ivlbds[1]`, `ivlbds[end]`].

For a mini-segment it throws:
- `:toofewpackets` if it has fewer than two packets;
- `:invalidsubtype` if its subtype is not in `0:3`;
- `:invalidsclkrate` unless its rate is positive;
- `:invaliddegree` if its degree is not in `1:23`;
- `:badwindowsize` if its window size is odd;
- `:timesoutoforder` unless its times are strictly increasing;
- `:boundsdisagree` if its first time is after the start of its interval, or its last
  time is before that start. As in SPICELIB, the last time is not compared with the end
  of the interval, so the times can end before it;
- `:zeroquaternion` if a quaternion is zero;
- `:badquatsign`, for a Hermite subtype, if a quaternion is farther from its
  predecessor than its negative is.
"""
function ckw06(h::DAFHandle, inst::Integer, ref::SegFrame, avflag::Bool, first::Real,
        last::Real, segid::AbstractString, npkts::AbstractVector{<:Integer},
        subtps::AbstractVector{<:Integer}, degres::AbstractVector{<:Integer}, packts, rates,
        sclkdp, ivlbds, sellst::Bool; ctx::SegContext=default_context())
    who = "ckw06"
    _daf_check_write(h, who)
    ref isa AbstractString && _seg_chkfstr(who, "ref", ref)
    _seg_chkfstr(who, "segid", segid)
    nmini = length(npkts)
    (length(subtps) == nmini && length(degres) == nmini) || throw(SpiceError(:badarg,
        "$who: `npkts`, `subtps` and `degres` must have the same length"))
    pk = _seg_vec(who, "packts", packts, -1)
    t = _seg_vec(who, "sclkdp", sclkdp, -1)
    rt = _seg_vec(who, "rates", rates, nmini)
    bd = _seg_vec(who, "ivlbds", ivlbds, nmini < 1 ? -1 : nmini + 1)
    last < first && throw(SpiceError(:baddescrtimes, "$who: segment start time is " *
        "$first; stop time is $last; bounds must be in nondecreasing order"))
    refcod = _seg_refcode(who, ctx, ref)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    nmini < 1 && throw(SpiceError(:invalidcount,
        "$who: mini-segment count was $nmini; this count must be positive"))
    for i in 1:nmini
        bd[i] >= bd[i+1] && throw(SpiceError(:boundsoutoforder, "$who: mini-segment interval " *
            "bounds at indices $i and $(i + 1) are $(bd[i]) and $(bd[i+1]); they must be " *
            "strictly increasing"))
    end
    (first < bd[1] || last > bd[nmini+1]) && throw(SpiceError(:coveragegap, "$who: segment " *
        "coverage [$first, $last] must be contained within the union of the mini-segment " *
        "intervals [$(bd[1]), $(bd[nmini+1])]"))
    winsizes = zeros(Int, nmini)
    pktend = 0
    eepix = 0
    for i in 1:nmini                     # invariant: mini-segments 1 .. i-1 are valid
        n = Int(npkts[i])
        n < 2 && throw(SpiceError(:toofewpackets, "$who: at least 2 packets are required " *
            "for CK type 6; $n were supplied in mini-segment $i"))
        subtyp = Int(subtps[i])
        (subtyp < 0 || subtyp > 3) && throw(SpiceError(:invalidsubtype,
            "$who: unexpected CK type 6 subtype $subtyp found in mini-segment $i"))
        pktsiz = _CK06_PKTSZS[subtyp+1]
        winsiz = isodd(subtyp) ? degres[i] + 1 : (degres[i] + 1) ÷ 2
        rt[i] <= 0.0 && throw(SpiceError(:invalidsclkrate,
            "$who: SCLK rate in mini-segment $i was $(rt[i]); the rate must be positive"))
        pktbeg = pktend + 1
        pktend = pktbeg - 1 + n * pktsiz
        (degres[i] < 1 || degres[i] > 23) && throw(SpiceError(:invaliddegree, "$who: the " *
            "interpolating polynomials of mini-segment $i have degree $(degres[i]); the " *
            "valid degree range is [1, 23]"))
        isodd(winsiz) && throw(SpiceError(:badwindowsize, "$who: the interpolating " *
            "polynomials of mini-segment $i have window size $winsiz and degree $(degres[i])"))
        winsizes[i] = winsiz
        bepix = eepix + 1
        eepix = bepix - 1 + n
        (length(t) >= eepix && length(pk) >= pktend) || throw(SpiceError(:badarg,
            "$who: `sclkdp` or `packts` is shorter than the packet counts require"))
        for k in bepix:eepix-1
            t[k] >= t[k+1] && throw(SpiceError(:timesoutoforder, "$who: in mini-segment $i, " *
                "epoch $(t[k]) having index $k is greater than or equal to its successor"))
        end
        t[bepix] > bd[i] && throw(SpiceError(:boundsdisagree, "$who: mini-segment interval " *
            "$i start time $(bd[i]) precedes mini-segment's first epoch $(t[bepix])"))
        # SPICELIB compares the last epoch with the interval start here, not with its end.
        t[eepix] < bd[i] && throw(SpiceError(:boundsdisagree, "$who: mini-segment interval " *
            "$i end time exceeds mini-segment's last epoch $(t[eepix])"))
        for j in 1:n
            a = pktbeg + pktsiz * (j - 1)
            (pk[a] == 0.0 && pk[a+1] == 0.0 && pk[a+2] == 0.0 && pk[a+3] == 0.0) &&
                throw(SpiceError(:zeroquaternion, "$who: the quaternion at index $j within " *
                    "mini-segment $i has magnitude zero"))
            if j > 1 && iseven(subtyp)
                p = a - pktsiz
                dneg = 0.0
                dpos = 0.0
                for c in 0:3
                    dneg += (pk[p+c] + pk[a+c])^2
                    dpos += (pk[p+c] - pk[a+c])^2
                end
                sqrt(dneg) < sqrt(dpos) && throw(SpiceError(:badquatsign, "$who: the " *
                    "quaternion at index $j within mini-segment $i is farther than its " *
                    "negative from its predecessor at index $(j - 1)"))
            end
        end
    end
    data = Float64[]
    pktend = 0
    eepix = 0
    pointers = Float64[]
    segend = 0
    for i in 1:nmini                     # invariant: mini-segments 1 .. i-1 are in data
        pktsiz = _CK06_PKTSZS[subtps[i]+1]
        n = Int(npkts[i])
        append!(data, view(pk, pktend+1:pktend+n*pktsiz))
        pktend += n * pktsiz
        bepix = eepix + 1
        eepix = bepix - 1 + n
        append!(data, view(t, bepix:eepix))
        for j in 1:(n-1)÷100
            push!(data, t[bepix-1+100j])
        end
        push!(data, rt[i], Float64(subtps[i]), Float64(winsizes[i]), Float64(n))
        minisz = n * (pktsiz + 1) + (n - 1) ÷ 100 + 4
        push!(pointers, Float64(segend + 1))
        segend += minisz
    end
    append!(data, bd)
    for i in 1:nmini÷100
        push!(data, bd[100i])
    end
    append!(data, pointers)
    push!(data, Float64(segend + 1), sellst ? 1.0 : -1.0, Float64(nmini))
    _ck_add!(h, inst, refcod, 6, avflag, first, last, segid, data)
    return nothing
end
