# Records of CK segments of types 2 and 3 through a DAF handle: `cknr02`, `ckgr02`,
# `cknr03` and `ckgr03`, as the SPICELIB routines of the same names. The descriptor is
# the packed summary of `dafgs` (2 doubles, 6 integers); CSPICE does not check that it
# belongs to the file of the handle, and neither does this file.

"Integer components of the CK descriptor `descr`; throws `SpiceError(:ckwrongdatatype)` unless the data type is `dtype`."
function _ck_descr_ic(who::AbstractString, descr::AbstractVector{<:Real}, dtype::Int)
    _, ic = dafus(descr, 2, 6)
    ic[3] == dtype || throw(SpiceError(:ckwrongdatatype, "$who: the data type of the " *
        "segment should be $dtype; the descriptor shows type $(ic[3])"))
    return ic
end

function _ck_check_recno(who::AbstractString, recno::Integer, nrec::Int)
    1 <= recno <= nrec || throw(SpiceError(:cknonexistrec, "$who: the requested record " *
        "number ($recno) does not exist; there are $nrec records in the segment"))
    return nothing
end

"""
    cknr02(handle, descr) -> Int

Return the number of pointing records in the CK type 2 segment with the descriptor
`descr`, as `cknr02_c`.

`descr` is the packed descriptor of the segment (5 doubles), as [`dafgs`](@ref) returns
it during a search on `handle`. The count follows from the addresses in the descriptor,
so the file is not read, and a closed handle gives the count too.

Throws `SpiceError(:ckwrongdatatype)` if the segment is not of type 2, and
`:dafinvalidparams` if `descr` has fewer than 5 doubles.

# Examples
```julia
h = dafopr("mission.bc")
dafbfs(h)
while daffna(h)
    descr = dafgs(h)
    _, ic = dafus(descr, 2, 6)
    ic[3] == 2 || continue
    for i in 1:cknr02(h, descr)
        start, stop, rate = ckgr02(h, descr, i)[1:3]
    end
end
dafcls(h)
```

See also [`ckgr02`](@ref), [`cknr03`](@ref).
"""
function cknr02(h::DAFHandle, descr::AbstractVector{<:Real})
    ic = _ck_descr_ic("cknr02", descr, 2)
    return _nint(((ic[6] - ic[5] + 1) * 100.0 + 1.0) / 1001.0)
end

"""
    ckgr02(handle, descr, recno) -> Vector{Float64}

Return the pointing record `recno` of the CK type 2 segment with the descriptor `descr`
in the file of `handle`, as `ckgr02_c`. Records count from 1, as in CSPICE.

The record has 10 doubles: the start and the stop of its interval (ticks), the clock rate
(seconds per tick), the quaternion of the C-matrix at the start (SPICE form, scalar
first) and the constant angular velocity (rad/s, components in the reference frame of the
segment). `descr` is as in [`cknr02`](@ref).

Throws `SpiceError`:
- `:ckwrongdatatype` if the segment is not of type 2;
- `:cknonexistrec` unless `1 <= recno <= cknr02(handle, descr)`;
- `:dafnosuchhandle` if the handle is closed;
- `:dafinvalidparams` if `descr` has fewer than 5 doubles.
"""
function ckgr02(h::DAFHandle, descr::AbstractVector{<:Real}, recno::Integer)
    ic = _ck_descr_ic("ckgr02", descr, 2)
    nrec = cknr02(h, descr)
    _ck_check_recno("ckgr02", recno, nrec)
    beg = ic[5]
    addr = beg + 8 * (recno - 1)
    prec = dafgda(h, addr, addr + 7)             # quaternion, angular velocity, rate
    start = dafgda(h, beg + 8nrec + recno - 1, beg + 8nrec + recno - 1)[1]
    stop = dafgda(h, beg + 9nrec + recno - 1, beg + 9nrec + recno - 1)[1]
    record = Float64[start, stop, prec[8], prec[1], prec[2], prec[3], prec[4], prec[5],
        prec[6], prec[7]]
    @assert length(record) == 10
    return record
end

"""
    cknr03(handle, descr) -> Int

Return the number of pointing instances in the CK type 3 segment with the descriptor
`descr` in the file of `handle`, as `cknr03_c`. `descr` is as in [`cknr02`](@ref).

Throws `SpiceError(:ckwrongdatatype)` if the segment is not of type 3,
`:dafnosuchhandle` if the handle is closed, and `:dafinvalidparams` if `descr` has fewer
than 5 doubles.
"""
function cknr03(h::DAFHandle, descr::AbstractVector{<:Real})
    ic = _ck_descr_ic("cknr03", descr, 3)
    return _nint(dafgda(h, ic[6], ic[6])[1])
end

"""
    ckgr03(handle, descr, recno) -> Vector{Float64}

Return the pointing instance `recno` of the CK type 3 segment with the descriptor `descr`
in the file of `handle`, as `ckgr03_c`. Records count from 1, as in CSPICE.

The record is the time tag (ticks) and the quaternion of the C-matrix (SPICE form, scalar
first): 5 doubles for a segment without angular velocity, and 8 doubles, the last three
the angular velocity (rad/s, components in the reference frame of the segment), for a
segment with it. The fourth integer of the descriptor (see [`dafus`](@ref)) is 1 for a
segment with angular velocity. `descr` is as in [`cknr02`](@ref).

Throws `SpiceError`:
- `:ckwrongdatatype` if the segment is not of type 3;
- `:cknonexistrec` unless `1 <= recno <= cknr03(handle, descr)`;
- `:dafnosuchhandle` if the handle is closed;
- `:dafinvalidparams` if `descr` has fewer than 5 doubles.
"""
function ckgr03(h::DAFHandle, descr::AbstractVector{<:Real}, recno::Integer)
    ic = _ck_descr_ic("ckgr03", descr, 3)
    psiz = ic[4] == 1 ? 7 : 4
    nrec = _nint(dafgda(h, ic[6], ic[6])[1])
    _ck_check_recno("ckgr03", recno, nrec)
    addr = ic[5] + psiz * (recno - 1)
    tick = ic[5] + psiz * nrec + recno - 1
    record = [dafgda(h, tick, tick); dafgda(h, addr, addr + psiz - 1)]
    @assert length(record) == psiz + 1
    return record
end
