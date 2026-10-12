# SPK writers: `spkopn`, `spkopa`, `spkcls` and the segment writers `spkwNN`. Each
# writer makes the input checks of the SPICELIB routine in the same order, packs the
# descriptor as `spkpds` does and writes the data, the directory and the trailer of the
# segment in the layout of CSPICE, so the data area is equal bit for bit to the one
# that CSPICE writes from the same input.
#
# Argument lists follow SPICE.jl. A table of states or packets is a vector of vectors
# (the form of SPICE.jl), a matrix with one state in each column, or a flat vector.
# The reference frame is a frame name, which is resolved against the kernel set of
# the keyword `ctx` (default: the default context; built-in frames need no kernel), or
# a frame ID. Errors have the short message of CSPICE in lower case as kind, for
# example `SpiceError(:invaliddegree, …)` for SPICE(INVALIDDEGREE).

"Reference frame of a segment writer: a frame name or a frame ID."
const SegFrame = Union{AbstractString,Integer}
"Where a segment writer resolves a frame name."
const SegContext = Union{SpiceContext,KernelSet}

_seg_ks(ctx::SpiceContext) = kernelset(ctx)
_seg_ks(ks::KernelSet) = ks

"Throws `SpiceError(:emptystring)` if the string argument `what` is empty, as the macro CHKFSTR of CSPICE."
function _seg_chkfstr(who::AbstractString, what::AbstractString, s::AbstractString)
    isempty(s) && throw(SpiceError(:emptystring, "$who: the string `$what` is empty"))
    return nothing
end

"""
    _seg_refcode(who, ctx, frame) -> Int

Frame ID of `frame` (a name resolved in `ctx`, or an ID). Throws
`SpiceError(:invalidrefframe)` if the name is unknown or the ID is 0.
"""
function _seg_refcode(who::AbstractString, ctx::SegContext, frame::SegFrame)
    code = frame isa Integer ? Int(frame) : namfrm(_seg_ks(ctx), frame)
    code == 0 && throw(SpiceError(:invalidrefframe,
        "$who: the reference frame $frame is not supported"))
    return code
end

"Throws `SpiceError(:segidtoolong)` if `segid` has more than 40 characters without trailing blanks."
function _seg_idlen(who::AbstractString, segid::AbstractString)
    _daf_lastnb(segid) > 40 && throw(SpiceError(:segidtoolong,
        "$who: segment identifier contains more than 40 characters"))
    return nothing
end

"Throws `SpiceError(:nonprintablechars)` if `segid` has a character outside the codes 32 to 126."
function _seg_idchars(who::AbstractString, segid::AbstractString)
    cu = codeunits(String(segid))
    for i in 1:_daf_lastnb(segid)
        (cu[i] < 0x20 || cu[i] > 0x7e) && throw(SpiceError(:nonprintablechars,
            "$who: the segment identifier contains nonprintable characters"))
    end
    return nothing
end

"The checks of `chckid` for a segment identifier: `:idstringtoolong`, then `:nonprintablechars`."
function _seg_chckid(who::AbstractString, segid::AbstractString)
    _daf_lastnb(segid) > 40 && throw(SpiceError(:idstringtoolong,
        "$who: the SPK segment identifier is invalid; its last non-blank character " *
        "is at position $(_daf_lastnb(segid)); the maximum allowed length is 40"))
    _seg_idchars(who, segid)
    return nothing
end

"""
    _seg_flat(who, what, x, rowlen) -> (flat::Vector{Float64}, n::Int, rowlen::Int)

The table `x` of `n` rows of `rowlen` doubles as one vector, row after row (the memory
layout that the C routine takes). `x` is a matrix with one row of the table in each
column, a vector of vectors, or a flat vector. With `rowlen == 0` the row length is
taken from `x` (not possible for a flat vector). Throws `SpiceError(:badarg)` if the
shape does not fit.
"""
function _seg_flat(who::AbstractString, what::AbstractString, x, rowlen::Int)
    if x isa AbstractMatrix{<:Real}
        m = size(x, 1)
        (rowlen == 0 || m == rowlen) || throw(SpiceError(:badarg,
            "$who: `$what` must have $rowlen rows (one table row in each column), it has $m"))
        return Float64[Float64(v) for v in vec(x)], size(x, 2), m
    elseif x isa AbstractVector{<:Real}
        rowlen > 0 || throw(SpiceError(:badarg,
            "$who: `$what` must be a matrix or a vector of vectors"))
        length(x) % rowlen == 0 || throw(SpiceError(:badarg,
            "$who: the length $(length(x)) of `$what` is not a multiple of $rowlen"))
        return Float64[Float64(v) for v in x], length(x) ÷ rowlen, rowlen
    elseif x isa AbstractVector
        n = length(x)
        m = rowlen == 0 ? (n == 0 ? 0 : length(x[1])) : rowlen
        flat = Vector{Float64}(undef, m * n)
        for (k, row) in enumerate(x)     # invariant: rows 1 .. k-1 are copied
            length(row) == m || throw(SpiceError(:badarg,
                "$who: element $k of `$what` has $(length(row)) values, expected $m"))
            for j in 1:m
                flat[(k-1)*m+j] = Float64(row[j])
            end
        end
        return flat, n, m
    end
    throw(SpiceError(:badarg, "$who: `$what` must be a matrix or a vector of vectors"))
end

"The vector `x` as `Vector{Float64}`; throws `SpiceError(:badarg)` unless it has `n` elements (`n < 0`: any length)."
function _seg_vec(who::AbstractString, what::AbstractString, x, n::Int)
    x isa AbstractVector{<:Real} || throw(SpiceError(:badarg, "$who: `$what` must be a vector of numbers"))
    (n < 0 || length(x) == n) || throw(SpiceError(:badarg,
        "$who: `$what` must have $n elements but has $(length(x))"))
    return Float64[Float64(v) for v in x]
end

"The three doubles of `v`; throws `SpiceError(:badarg)` unless `v` has `n` elements."
function _seg_fixed(who::AbstractString, what::AbstractString, v, n::Int)
    length(v) == n || throw(SpiceError(:badarg, "$who: `$what` must have $n elements but has $(length(v))"))
    return Float64[Float64(x) for x in v]
end

"""
    _seg_increasing(who, epochs) -> Nothing

Throws `SpiceError(:timesoutoforder)` unless `epochs` is strictly increasing.
"""
function _seg_increasing(who::AbstractString, epochs::Vector{Float64})
    for i in 2:length(epochs)            # invariant: epochs[1:i-1] is strictly increasing
        epochs[i] <= epochs[i-1] && throw(SpiceError(:timesoutoforder,
            "$who: epoch $(epochs[i]) having index $i is not greater than its " *
            "predecessor $(epochs[i-1])"))
    end
    return nothing
end

"""
    spkpds(body, center, frame, type, first, last; ctx=default_context()) -> Vector{Float64}

Packed descriptor (five doubles) of an SPK segment, as `spkpds_c`; the two address
components are zero. Preconditions (each throws `SpiceError`): `body != 0`
(`:barycenterephem`), `body != center` (`:bodyandcentersame`), a known frame
(`:invalidrefframe`), `first < last` (`:baddescrtimes`), `1 <= type <= 1000`
(`:unknownspktype`).
"""
function spkpds(body::Integer, center::Integer, frame::SegFrame, type::Integer, first::Real,
        last::Real; ctx::SegContext=default_context())
    body == 0 && throw(SpiceError(:barycenterephem, "spkpds: you've attempted to create " *
        "a segment for the solar system barycenter; this is not supported"))
    body == center && throw(SpiceError(:bodyandcentersame, "spkpds: you've attempted to " *
        "create a segment for a body relative to itself; the body ID code was $body"))
    refcod = _seg_refcode("spkpds", ctx, frame)
    first >= last && throw(SpiceError(:baddescrtimes, "spkpds: the segment start time " *
        "$first is at or after the segment stop time $last"))
    (type <= 0 || type > 1000) && throw(SpiceError(:unknownspktype,
        "spkpds: the type $type is not supported within the SPK system"))
    return dafps(Float64[first, last], Int[body, center, refcod, type, 0, 0])
end

"Add an SPK segment of type `dtype` with the checked descriptor values and the complete `data`."
function _spk_add!(h::DAFHandle, body::Integer, center::Integer, refcod::Int, dtype::Int,
        first::Real, last::Real, segid::AbstractString, data::Vector{Float64})
    _daf_add_array!(h, Float64[first, last], Int[body, center, refcod, dtype, 0, 0], segid, data)
    return nothing
end

"Add an SPK segment whose descriptor is packed by `spkpds` (with its checks)."
function _spk_add_pds!(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame,
        dtype::Int, first::Real, last::Real, segid::AbstractString, data::Vector{Float64},
        ctx::SegContext)
    descr = spkpds(body, center, frame, dtype, first, last; ctx=ctx)
    dafbna(h, descr, segid)
    dafada(h, data)
    dafena(h)
    return nothing
end

"`data` followed by every 100th epoch from the first `ndir * 100` epochs (the directory of the epoch table)."
function _seg_directory!(data::Vector{Float64}, epochs::Vector{Float64}, ndir::Int)
    for i in 1:ndir
        push!(data, epochs[100i])
    end
    return data
end

"""
    _seg_coverage_gap(who, first, last, btime, ltime, tol) -> Nothing

Throws `SpiceError(:coveragegap)` if the descriptor interval [`first`, `last`] passes
the data interval [`btime`, `ltime`] by more than `tol`.
"""
function _seg_coverage_gap(who::AbstractString, first::Float64, last::Float64,
        btime::Float64, ltime::Float64, tol::Float64)
    first < btime - tol && throw(SpiceError(:coveragegap, "$who: the segment descriptor " *
        "start time $first is too much less than the beginning time of the segment " *
        "data $btime; the difference is $(btime - first) s, the tolerance is $tol s"))
    last > ltime + tol && throw(SpiceError(:coveragegap, "$who: the segment descriptor " *
        "end time $last is too much greater than the end time of the segment data " *
        "$ltime; the difference is $(last - ltime) s, the tolerance is $tol s"))
    return nothing
end

# ---------------------------------------------------------------------------------------
# Open and close
# ---------------------------------------------------------------------------------------

"""
    spkopn(path, ifname, ncomch=0) -> DAFHandle

Create the new SPK file `path` with the internal file name `ifname` and room for
`ncomch` comment characters, and open it for write, as `spkopn_c`. Preconditions
(each throws `SpiceError`): `path` and `ifname` are not empty (`:emptystring`), `path`
does not exist (`:fileopenfailed`). Postcondition: the file is byte for byte the file
that `spkopn_c` creates. Close the handle with `spkcls` (after at least one segment)
or `dafcls`; load the file with `furnsh` only after the close. The handle belongs to
the caller and must not be shared between tasks without a lock.
"""
function spkopn(path::AbstractString, ifname::AbstractString, ncomch::Integer=0)
    _seg_chkfstr("spkopn", "path", path)
    _seg_chkfstr("spkopn", "ifname", ifname)
    ncomr = ncomch > 0 ? (Int(ncomch) - 1) ÷ _DAF_CRECLEN + 1 : 0
    return dafonw(path, "SPK", 2, 6, ifname, ncomr)
end

"""
    spkopa(path) -> DAFHandle

Open the existing SPK file `path` to append segments, as `spkopa_c`. Preconditions
(each throws `SpiceError`): the file exists (`:filenotfound`), is a DAF of type SPK
(`:fileisnotspk`) and has the native byte order (`:unsupportedbff`). Unload the file
from every context before the call, and load it again after `spkcls`.
"""
function spkopa(path::AbstractString)
    _seg_chkfstr("spkopa", "path", path)
    isfile(path) || throw(SpiceError(:filenotfound,
        "spkopa: the file $path is not recognized as an existing file"))
    arch, kind = getfat(path)
    (arch == "DAF" && kind == "SPK") || throw(SpiceError(:fileisnotspk, "spkopa: the file " *
        "$path was not an SPK file; its architecture and type are $arch and $kind"))
    return dafopw(path)
end

"Close `h` as `spkcls`/`ckcls` do: a file open for write must hold a segment (`:nosegmentsfound`)."
function _seg_close(h::DAFHandle, what::AbstractString)
    h.isopen || return nothing
    if h.writable
        dafbfs(h)
        daffna(h) || throw(SpiceError(:nosegmentsfound, "no segments were found in the " *
            "$what file $(h.path); there must be at least one segment in the file when " *
            "it is closed"))
    end
    dafcls(h)
    return nothing
end

"""
    spkcls(handle) -> Nothing

Close the SPK file of `handle`, as `spkcls_c`. Precondition: a file that is open for
write holds at least one segment, else `SpiceError(:nosegmentsfound)` and the handle
stays open (close it with `dafcls`). Postcondition: the file is complete on disk and
can be loaded with `furnsh`.
"""
spkcls(h::DAFHandle) = _seg_close(h, "SPK")

# ---------------------------------------------------------------------------------------
# Types 2, 3 and 20: Chebyshev polynomials on equal intervals
# ---------------------------------------------------------------------------------------

"Checks and data that `spkw02` and `spkw03` share; `per` is 3 (positions) or 6 (states)."
function _spk_cheby(who::AbstractString, dtype::Int, per::Int, h::DAFHandle, body, center,
        frame, first, last, segid, intlen, cdata, btime, ctx)
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    cd, n, ninrec = _seg_flat(who, "cdata", cdata, 0)
    first, last, intlen, btime = Float64(first), Float64(last), Float64(intlen), Float64(btime)
    n <= 0 && throw(SpiceError(:numcoeffsnotpos,
        "$who: the number of sets of coordinate coefficients is not positive: $n"))
    ninrec % per == 0 && ninrec > 0 || throw(SpiceError(:badarg, "$who: every record of " *
        "`cdata` must hold $per * (degree + 1) coefficients, it holds $ninrec"))
    polydg = ninrec ÷ per - 1
    (polydg < 0 || polydg > 27) && throw(SpiceError(:invaliddegree, "$who: the " *
        "interpolating polynomials have degree $polydg; the valid degree range is [0, 27]"))
    intlen <= 0.0 && throw(SpiceError(:intlennotpos, "$who: the interval length is not positive: $intlen"))
    refcod = _seg_refcode(who, ctx, frame)
    first > last && throw(SpiceError(:baddescrtimes,
        "$who: the segment start time $first is greater than the segment end time $last"))
    tol = max(abs(first), abs(last)) * 1e-13
    ltime = btime + n * intlen
    _seg_coverage_gap(who, first, last, btime, ltime, tol)
    _seg_chckid(who, segid)
    data = Vector{Float64}(undef, n * (ninrec + 2) + 4)
    radius = intlen / 2.0
    for i in 1:n                         # invariant: records 1 .. i-1 are in data
        base = (i - 1) * (ninrec + 2)
        data[base+1] = btime + radius + (i - 1) * intlen
        data[base+2] = radius
        copyto!(data, base + 3, cd, (i - 1) * ninrec + 1, ninrec)
    end
    data[end-3] = btime
    data[end-2] = intlen
    data[end-1] = Float64(ninrec + 2)
    data[end] = Float64(n)
    _spk_add!(h, body, center, refcod, dtype, first, last, segid, data)
    return nothing
end

"""
    spkw02(handle, body, center, frame, first, last, segid, intlen, cdata, btime; ctx=default_context())

Write an SPK segment of type 2 (Chebyshev polynomials for the position on equal
intervals), as `spkw02_c`. `cdata` holds one record for each interval, each with
`3 * (polydg + 1)` coefficients (X, then Y, then Z); the number of records `n` and the
degree `polydg` follow from its shape. `btime` is the start of the first interval and
`intlen` the length of every interval in seconds.

Preconditions (each throws `SpiceError`): `n > 0` (`:numcoeffsnotpos`),
`0 <= polydg <= 27` (`:invaliddegree`), `intlen > 0` (`:intlennotpos`), a known frame
(`:invalidrefframe`), `first <= last` (`:baddescrtimes`), [`first`, `last`] inside
[`btime`, `btime + n*intlen`] with the relative tolerance 1e-13 (`:coveragegap`), a
printable `segid` of at most 40 characters (`:idstringtoolong`, `:nonprintablechars`).
Postcondition: the segment is in the file, with the data area of `spkw02_c`.
"""
spkw02(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real, last::Real,
    segid::AbstractString, intlen::Real, cdata, btime::Real; ctx::SegContext=default_context()) =
    _spk_cheby("spkw02", 2, 3, h, body, center, frame, first, last, segid, intlen, cdata, btime, ctx)

"""
    spkw03(handle, body, center, frame, first, last, segid, intlen, cdata, btime; ctx=default_context())

Write an SPK segment of type 3 (Chebyshev polynomials for position and velocity on
equal intervals), as `spkw03_c`. `cdata` holds one record for each interval, each
with `6 * (polydg + 1)` coefficients. The preconditions are those of [`spkw02`](@ref).
"""
spkw03(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real, last::Real,
    segid::AbstractString, intlen::Real, cdata, btime::Real; ctx::SegContext=default_context()) =
    _spk_cheby("spkw03", 3, 6, h, body, center, frame, first, last, segid, intlen, cdata, btime, ctx)

"""
    spkw20(handle, body, center, frame, first, last, segid, intlen, cdata, dscale, tscale,
           initjd, initfr; ctx=default_context())

Write an SPK segment of type 20 (Chebyshev polynomials for the velocity, with the
position at the interval midpoint), as `spkw20_c`. `cdata` holds one record for each
interval, each with `3 * (polydg + 2)` values; `intlen` is in days; `dscale` and
`tscale` are the distance scale (km) and the time scale (s); `initjd + initfr` is the
TDB Julian date of the start of the first interval.

Preconditions (each throws `SpiceError`): `n > 0` (`:invalidcount`),
`0 <= polydg <= 50` (`:invaliddegree`), `intlen > 0` (`:intlennotpos`), a known frame
(`:invalidrefframe`), `first <= last` (`:baddescrtimes`), positive scales
(`:nonpositivescale`), the descriptor interval inside the data interval
(`:coveragegap`), a valid `segid` (`:idstringtoolong`, `:nonprintablechars`).
"""
function spkw20(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, intlen::Real, cdata, dscale::Real, tscale::Real,
        initjd::Real, initfr::Real; ctx::SegContext=default_context())
    who = "spkw20"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    cd, n, ninrec = _seg_flat(who, "cdata", cdata, 0)
    first, last, intlen = Float64(first), Float64(last), Float64(intlen)
    dscale, tscale, initjd, initfr = Float64(dscale), Float64(tscale), Float64(initjd), Float64(initfr)
    n <= 0 && throw(SpiceError(:invalidcount,
        "$who: the number of sets of coordinate coefficients is not positive: $n"))
    (ninrec % 3 == 0 && ninrec >= 6) || throw(SpiceError(:badarg, "$who: every record of " *
        "`cdata` must hold 3 * (degree + 2) values, it holds $ninrec"))
    polydg = ninrec ÷ 3 - 2
    polydg > 50 && throw(SpiceError(:invaliddegree, "$who: the interpolating polynomials " *
        "have degree $polydg; the valid degree range is [0, 50]"))
    intlen <= 0.0 && throw(SpiceError(:intlennotpos, "$who: the interval length is not positive: $intlen"))
    refcod = _seg_refcode(who, ctx, frame)
    first > last && throw(SpiceError(:baddescrtimes,
        "$who: the segment start time $first is greater than the segment end time $last"))
    dscale <= 0.0 && throw(SpiceError(:nonpositivescale, "$who: the distance scale is not positive: $dscale"))
    tscale <= 0.0 && throw(SpiceError(:nonpositivescale, "$who: the time scale is not positive: $tscale"))
    btime = 86400.0 * (initjd - 2451545.0 + initfr)
    ltime = btime + n * intlen * 86400.0
    tol = max(abs(btime), abs(ltime)) * 1e-13
    _seg_coverage_gap(who, first, last, btime, ltime, tol)
    _seg_chckid(who, segid)
    data = append!(cd, (dscale, tscale, initjd, initfr, intlen, Float64(ninrec), Float64(n)))
    _spk_add!(h, body, center, refcod, 20, first, last, segid, data)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Types 5, 8, 9, 12, 13, 18: discrete states
# ---------------------------------------------------------------------------------------

"""
    spkw05(handle, body, center, frame, first, last, segid, gm, states, epochs; ctx=default_context())

Write an SPK segment of type 5 (two-body propagation between discrete states), as
`spkw05_c`. `states` holds `n` states (x, y, z, dx, dy, dz) in km and km/s, `epochs`
their `n` epochs, `gm` the gravitational parameter of the central body.

Preconditions (each throws `SpiceError`): `gm > 0` (`:nonpositivemass`), a known frame
(`:invalidrefframe`), `n > 0` (`:numstatesnotpos`), strictly increasing epochs
(`:unorderedtimes`), `first <= last` (`:baddescrtimes`), a valid `segid`
(`:nonprintablechars`, `:segidtoolong`).
"""
function spkw05(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, gm::Real, states, epochs;
        ctx::SegContext=default_context())
    who = "spkw05"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    st, n, _ = _seg_flat(who, "states", states, 6)
    ep = _seg_vec(who, "epochs", epochs, n)
    gm <= 0 && throw(SpiceError(:nonpositivemass, "$who: GM = $gm; non-positive gravitational parameter"))
    refcod = _seg_refcode(who, ctx, frame)
    n <= 0 && throw(SpiceError(:numstatesnotpos,
        "$who: the number of states and epochs is not positive: $n"))
    for i in 2:n
        ep[i] <= ep[i-1] && throw(SpiceError(:unorderedtimes, "$who: epoch $i is out of order"))
    end
    first > last && throw(SpiceError(:baddescrtimes,
        "$who: the segment start time $first is greater than the segment end time $last"))
    _seg_idchars(who, segid)
    _seg_idlen(who, segid)
    data = append!(st, ep)
    _seg_directory!(data, ep, n ÷ 100)
    push!(data, Float64(gm), Float64(n))
    _spk_add!(h, body, center, refcod, 5, first, last, segid, data)
    return nothing
end

"""
    spkw08(handle, body, center, frame, first, last, segid, degree, states, epoch1, step; ctx=default_context())

Write an SPK segment of type 8 (Lagrange interpolation of states at equal steps), as
`spkw08_c`. The state number `k` (from 1) is at `epoch1 + (k-1)*step`.

Preconditions (each throws `SpiceError`): a known frame (`:invalidrefframe`), a valid
`segid` (`:segidtoolong`, `:nonprintablechars`), `1 <= degree <= 27`
(`:invaliddegree`), `n > degree` (`:toofewstates`), `first < last`
(`:baddescrtimes`), `step > 0` (`:invalidstepsize`), the descriptor interval inside
the data interval with the relative tolerance 1e-13 (`:coveragegap`).
"""
function spkw08(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, degree::Integer, states, epoch1::Real, step::Real;
        ctx::SegContext=default_context())
    who = "spkw08"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    st, n, _ = _seg_flat(who, "states", states, 6)
    first, last, begtim, step = Float64(first), Float64(last), Float64(epoch1), Float64(step)
    refcod = _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    (degree < 1 || degree > 27) && throw(SpiceError(:invaliddegree, "$who: the " *
        "interpolating polynomials have degree $degree; the valid degree range is [1, 27]"))
    n <= degree && throw(SpiceError(:toofewstates, "$who: at least $(degree + 1) states " *
        "are required to define a polynomial of degree $degree; number of states supplied: $n"))
    first >= last && throw(SpiceError(:baddescrtimes, "$who: the segment start time " *
        "$first is greater than or equal to the segment end time $last"))
    step <= 0.0 && throw(SpiceError(:invalidstepsize, "$who: the step size must be > 0 but was $step"))
    tol = max(abs(first), abs(last)) * 1e-13
    ltime = begtim + (n - 1) * step
    _seg_coverage_gap(who, first, last, begtim, ltime, tol)
    data = push!(st, begtim, step, Float64(degree), Float64(n))
    _spk_add!(h, body, center, refcod, 8, first, last, segid, data)
    return nothing
end

"""
    spkw09(handle, body, center, frame, first, last, segid, degree, states, epochs; ctx=default_context())

Write an SPK segment of type 9 (Lagrange interpolation of states at unequal steps), as
`spkw09_c`.

Preconditions (each throws `SpiceError`): a known frame (`:invalidrefframe`),
`first < last` (`:baddescrtimes`), a valid `segid` (`:nonprintablechars`,
`:segidtoolong`), `1 <= degree <= 27` (`:invaliddegree`), `n > degree`
(`:toofewstates`), strictly increasing epochs (`:timesoutoforder`),
`epochs[1] <= first` and `epochs[n] >= last` (`:baddescrtimes`).
"""
function spkw09(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, degree::Integer, states, epochs;
        ctx::SegContext=default_context())
    who = "spkw09"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    st, n, _ = _seg_flat(who, "states", states, 6)
    ep = _seg_vec(who, "epochs", epochs, n)
    refcod = _seg_refcode(who, ctx, frame)
    first >= last && throw(SpiceError(:baddescrtimes,
        "$who: the segment start time $first is greater than the segment end time $last"))
    _seg_idchars(who, segid)
    _seg_idlen(who, segid)
    (degree < 1 || degree > 27) && throw(SpiceError(:invaliddegree, "$who: the " *
        "interpolating polynomials have degree $degree; the valid degree range is [1, 27]"))
    n <= degree && throw(SpiceError(:toofewstates, "$who: at least $(degree + 1) states " *
        "are required to define a polynomial of degree $degree; number of states supplied: $n"))
    _seg_increasing(who, ep)
    ep[1] > first && throw(SpiceError(:baddescrtimes,
        "$who: segment start time $first precedes first epoch $(ep[1])"))
    ep[n] < last && throw(SpiceError(:baddescrtimes,
        "$who: segment end time $last follows last epoch $(ep[n])"))
    data = append!(st, ep)
    _seg_directory!(data, ep, (n - 1) ÷ 100)
    push!(data, Float64(degree), Float64(n))
    _spk_add!(h, body, center, refcod, 9, first, last, segid, data)
    return nothing
end

"""
    spkw12(handle, body, center, frame, first, last, segid, degree, states, epoch1, step; ctx=default_context())

Write an SPK segment of type 12 (Hermite interpolation of states at equal steps), as
`spkw12_c`.

Preconditions (each throws `SpiceError`): a known frame (`:invalidrefframe`), a valid
`segid` (`:segidtoolong`, `:nonprintablechars`), an odd `degree` in 1:27
(`:invaliddegree`), `n >= (degree + 1) ÷ 2` (`:toofewstates`), `first < last`
(`:baddescrtimes`), `step > 0` (`:invalidstepsize`), the descriptor interval inside
the data interval (`:coveragegap`), and the checks of [`spkpds`](@ref)
(`:barycenterephem`, `:bodyandcentersame`).
"""
function spkw12(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, degree::Integer, states, epoch1::Real, step::Real;
        ctx::SegContext=default_context())
    who = "spkw12"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    st, n, _ = _seg_flat(who, "states", states, 6)
    first, last, begtim, step = Float64(first), Float64(last), Float64(epoch1), Float64(step)
    winsiz = (degree + 1) ÷ 2
    _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    (degree < 1 || degree > 27) && throw(SpiceError(:invaliddegree, "$who: the " *
        "interpolating polynomials have degree $degree; the valid degree range is [1, 27]"))
    iseven(degree) && throw(SpiceError(:invaliddegree, "$who: the interpolating " *
        "polynomials have degree $degree; for SPK type 12, the degree must be odd"))
    n < winsiz && throw(SpiceError(:toofewstates, "$who: at least $winsiz states are " *
        "required to define a Hermite polynomial of degree $degree; number of states supplied: $n"))
    first >= last && throw(SpiceError(:baddescrtimes,
        "$who: the segment start time $first is greater than the segment end time $last"))
    step <= 0.0 && throw(SpiceError(:invalidstepsize, "$who: the step size must be > 0 but was $step"))
    tol = max(abs(first), abs(last)) * 1e-13
    ltime = begtim + (n - 1) * step
    _seg_coverage_gap(who, first, last, begtim, ltime, tol)
    data = push!(st, begtim, step, Float64(winsiz - 1), Float64(n))
    _spk_add_pds!(h, body, center, frame, 12, first, last, segid, data, ctx)
    return nothing
end

"""
    spkw13(handle, body, center, frame, first, last, segid, degree, states, epochs; ctx=default_context())

Write an SPK segment of type 13 (Hermite interpolation of states at unequal steps),
as `spkw13_c`. `states` holds `n` states (x, y, z, dx, dy, dz) in km and km/s and
`epochs` their `n` epochs in seconds past J2000 TDB.

Preconditions (each throws `SpiceError`, in the order of SPICELIB): a known frame
(`:invalidrefframe`), a `segid` of at most 40 characters (`:segidtoolong`) that are
printable (`:nonprintablechars`), `1 <= degree <= 27` and `degree` odd
(`:invaliddegree`), `n >= (degree + 1) ÷ 2` (`:toofewstates`), `first < last`
(`:baddescrtimes`), strictly increasing epochs (`:timesoutoforder`),
`epochs[1] <= first` and `epochs[n] >= last` (`:baddescrtimes`), `body != 0`
(`:barycenterephem`), `body != center` (`:bodyandcentersame`). An empty `frame` or
`segid` throws `:emptystring`, as the C wrapper does.

Postcondition: the file has one more segment with the descriptor of `spkpds_c` and
the data area of `spkw13_c`: the states, the epochs, every 100th epoch for the first
`(n-1) ÷ 100` hundreds, the window size minus one and `n`.
"""
function spkw13(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, degree::Integer, states, epochs;
        ctx::SegContext=default_context())
    who = "spkw13"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    st, n, _ = _seg_flat(who, "states", states, 6)
    ep = _seg_vec(who, "epochs", epochs, n)
    winsiz = (degree + 1) ÷ 2
    _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    (degree < 1 || degree > 27) && throw(SpiceError(:invaliddegree, "$who: the " *
        "interpolating polynomials have degree $degree; the valid degree range is [1, 27]"))
    iseven(degree) && throw(SpiceError(:invaliddegree, "$who: the interpolating " *
        "polynomials have degree $degree; for SPK type 13, the degree must be odd"))
    n < winsiz && throw(SpiceError(:toofewstates, "$who: at least $winsiz states are " *
        "required to define a Hermite polynomial of degree $degree; number of states supplied: $n"))
    first >= last && throw(SpiceError(:baddescrtimes,
        "$who: the segment start time $first is greater than the segment end time $last"))
    _seg_increasing(who, ep)
    ep[1] > first && throw(SpiceError(:baddescrtimes,
        "$who: segment start time $first precedes first epoch $(ep[1])"))
    ep[n] < last && throw(SpiceError(:baddescrtimes,
        "$who: segment end time $last follows last epoch $(ep[n])"))
    data = append!(st, ep)
    _seg_directory!(data, ep, (n - 1) ÷ 100)
    push!(data, Float64(winsiz - 1), Float64(n))
    _spk_add_pds!(h, body, center, frame, 13, first, last, segid, data, ctx)
    return nothing
end

"Subtype code of an SPK type 18 segment from the code or the SPICE.jl symbol `:S18TP0`, `:S18TP1`."
_spk18_subtype(s::Integer) = Int(s)
_spk18_subtype(s::Symbol) = s === :S18TP0 ? 0 : s === :S18TP1 ? 1 : -1

"""
    spkw18(handle, subtyp, body, center, frame, first, last, segid, degree, packts, epochs; ctx=default_context())

Write an SPK segment of type 18 (ESOC/DDID interpolation), as `spkw18_c`. `subtyp` is
0 or `:S18TP0` (Hermite, packets of 12 values) or 1 or `:S18TP1` (Lagrange, packets
of 6 values).

Preconditions (each throws `SpiceError`): a known subtype (`:invalidvalue`), a known
frame (`:invalidrefframe`), a valid `segid` (`:segidtoolong`, `:nonprintablechars`),
`1 <= degree <= 15` and an even window size, that is a degree of 3 mod 4 for subtype
0 and an odd degree for subtype 1 (`:invaliddegree`), `n >= 2` (`:toofewstates`),
`first <= last` (`:baddescrtimes`), strictly increasing epochs (`:timesoutoforder`),
`epochs[1] <= first` and `epochs[n] >= last` (`:baddescrtimes`).
"""
function spkw18(h::DAFHandle, subtyp::Union{Integer,Symbol}, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, degree::Integer,
        packts, epochs; ctx::SegContext=default_context())
    who = "spkw18"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    sub = _spk18_subtype(subtyp)
    (sub == 0 || sub == 1) || throw(SpiceError(:invalidvalue,
        "$who: unexpected SPK type 18 subtype requested: $subtyp"))
    packsz = sub == 0 ? 12 : 6
    winsiz = sub == 0 ? (degree + 1) ÷ 2 : degree + 1
    pk, n, _ = _seg_flat(who, "packts", packts, packsz)
    ep = _seg_vec(who, "epochs", epochs, n)
    refcod = _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    (degree < 1 || degree > 15) && throw(SpiceError(:invaliddegree, "$who: the " *
        "interpolating polynomials have degree $degree; the valid degree range is [1, 15]"))
    isodd(winsiz) && throw(SpiceError(:invaliddegree, "$who: the interpolating " *
        "polynomials have degree $degree; for SPK type 18, the degree must be " *
        "equivalent to 3 mod 4 for Hermite interpolation and odd for Lagrange interpolation"))
    n < 2 && throw(SpiceError(:toofewstates,
        "$who: at least 2 packets are required for SPK type 18; number of packets supplied: $n"))
    first > last && throw(SpiceError(:baddescrtimes,
        "$who: the segment start time $first is greater than the segment end time $last"))
    _seg_increasing(who, ep)
    ep[1] > first && throw(SpiceError(:baddescrtimes,
        "$who: segment start time $first precedes first epoch $(ep[1])"))
    ep[n] < last && throw(SpiceError(:baddescrtimes,
        "$who: segment end time $last follows last epoch $(ep[n])"))
    data = append!(pk, ep)
    _seg_directory!(data, ep, (n - 1) ÷ 100)
    push!(data, Float64(sub), Float64(winsiz), Float64(n))
    _spk_add!(h, body, center, refcod, 18, first, last, segid, data)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Types 15 and 17: precessing conic elements
# ---------------------------------------------------------------------------------------

"Norm of the 3-vector `v`, with the scaling of `vnorm` of SPICELIB (same rounding)."
function _seg_vnorm(v)
    vmax = max(abs(v[1]), abs(v[2]), abs(v[3]))
    vmax == 0.0 && return 0.0
    a = v[1] / vmax
    b = v[2] / vmax
    c = v[3] / vmax
    return vmax * sqrt(a * a + b * b + c * c)
end

"Unit vector of the 3-vector `v`, or zero for the zero vector, as `vhat` of SPICELIB (same rounding)."
function _seg_vhat(v)
    m = _seg_vnorm(v)
    return m > 0.0 ? (v[1] / m, v[2] / m, v[3] / m) : (0.0, 0.0, 0.0)
end

"""
    spkw15(handle, body, center, frame, first, last, segid, epoch, tp, pa, p, ecc, j2flg,
           pv, gm, j2, radius; ctx=default_context())

Write an SPK segment of type 15 (precessing conic propagation), as `spkw15_c`: `epoch`
of periapsis, trajectory pole `tp`, periapsis vector `pa`, semi-latus rectum `p`,
eccentricity `ecc`, J2 processing flag `j2flg`, central body pole `pv`, `gm`, `j2`
and the equatorial `radius` of the central body. The three vectors are stored as unit
vectors.

Preconditions (each throws `SpiceError`): `p > 0` (`:badlatusrectum`), `ecc >= 0`
(`:badeccentricity`), `gm > 0` (`:nonpositivemass`), non-zero `tp`, `pa`, `pv`
(`:badvector`), `radius >= 0` (`:badradius`), `tp` orthogonal to `pa` within 1e-5
(`:badinitstate`), a valid `segid` (`:segidtoolong`, `:nonprintablechars`), and the
checks of [`spkpds`](@ref).
"""
function spkw15(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, epoch::Real, tp, pa, p::Real, ecc::Real,
        j2flg::Real, pv, gm::Real, j2::Real, radius::Real; ctx::SegContext=default_context())
    who = "spkw15"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    tp = _seg_fixed(who, "tp", tp, 3)
    pa = _seg_fixed(who, "pa", pa, 3)
    pv = _seg_fixed(who, "pv", pv, 3)
    mypa = _seg_vhat(pa)
    mytp = _seg_vhat(tp)
    record = Float64[epoch, mytp..., mypa..., p, ecc, j2flg, _seg_vhat(pv)..., gm, j2, radius]
    p <= 0 && throw(SpiceError(:badlatusrectum, "$who: the semi-latus rectum must be positive; it was $p"))
    ecc < 0 && throw(SpiceError(:badeccentricity, "$who: the eccentricity must be non-negative; it was $ecc"))
    gm <= 0 && throw(SpiceError(:nonpositivemass, "$who: the mass of the central body must be positive; GM was $gm"))
    all(iszero, tp) && throw(SpiceError(:badvector, "$who: the trajectory pole vector had length zero"))
    all(iszero, pa) && throw(SpiceError(:badvector, "$who: the periapse vector had length zero"))
    all(iszero, pv) && throw(SpiceError(:badvector, "$who: the central pole vector had length zero"))
    radius < 0 && throw(SpiceError(:badradius, "$who: the central body radius was negative: $radius"))
    dot = mypa[1] * mytp[1] + mypa[2] * mytp[2] + mypa[3] * mytp[3]
    abs(dot) > 1e-5 && throw(SpiceError(:badinitstate,
        "$who: the periapsis and trajectory pole vectors are not orthogonal (dot product $dot)"))
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    _spk_add_pds!(h, body, center, frame, 15, first, last, segid, record, ctx)
    return nothing
end

"""
    spkw17(handle, body, center, frame, first, last, segid, epoch, eqel, rapol, decpol; ctx=default_context())

Write an SPK segment of type 17 (equinoctial elements), as `spkw17_c`: `eqel` holds
the nine equinoctial elements at `epoch`, and `rapol`, `decpol` are the right
ascension and declination of the pole of the reference plane.

Preconditions (each throws `SpiceError`): `eqel[1] > 0` (`:badsemiaxis`), eccentricity
`sqrt(eqel[2]^2 + eqel[3]^2) <= 0.9` (`:badeccentricity`), a valid `segid`
(`:segidtoolong`, `:nonprintablechars`), and the checks of [`spkpds`](@ref).
"""
function spkw17(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, epoch::Real, eqel, rapol::Real, decpol::Real;
        ctx::SegContext=default_context())
    who = "spkw17"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    el = _seg_fixed(who, "eqel", eqel, 9)
    record = Float64[epoch, el..., rapol, decpol]
    a = record[2]
    ecc = sqrt(record[3] * record[3] + record[4] * record[4])
    a <= 0.0 && throw(SpiceError(:badsemiaxis, "$who: the semimajor axis must be positive; it was $a"))
    ecc > 0.9 && throw(SpiceError(:badeccentricity,
        "$who: the eccentricity $ecc is greater than 0.9"))
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    _spk_add_pds!(h, body, center, frame, 17, first, last, segid, record, ctx)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Types 1 and 21: modified difference arrays
# ---------------------------------------------------------------------------------------

"""
    spkw01(handle, body, center, frame, first, last, segid, dlines, epochs; ctx=default_context())

Write an SPK segment of type 1 (modified difference arrays), as `spkw01` of SPICELIB
(CSPICE has no C wrapper; the count `n` follows from the input). `dlines` holds `n`
difference lines of 71 doubles and `epochs` the final epoch of each line.

Preconditions (each throws `SpiceError`): a known frame (`:invalidrefframe`), a valid
`segid` (`:segidtoolong`, `:nonprintablechars`), `n >= 1` (`:invalidcount`),
`first < last` (`:baddescrtimes`), strictly increasing epochs (`:timesoutoforder`),
`epochs[n] >= last` (`:baddescrtimes`), and the checks of [`spkpds`](@ref).
"""
function spkw01(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, dlines, epochs; ctx::SegContext=default_context())
    who = "spkw01"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    dl, n, _ = _seg_flat(who, "dlines", dlines, 71)
    ep = _seg_vec(who, "epochs", epochs, n)
    _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    n < 1 && throw(SpiceError(:invalidcount,
        "$who: the difference line count was $n; the count must be at least one"))
    first >= last && throw(SpiceError(:baddescrtimes,
        "$who: the segment start time $first is greater than the segment end time $last"))
    _seg_increasing(who, ep)
    ep[n] < last && throw(SpiceError(:baddescrtimes,
        "$who: segment end time $last follows last epoch $(ep[n])"))
    data = append!(dl, ep)
    _seg_directory!(data, ep, n ÷ 100)
    push!(data, Float64(n))
    _spk_add_pds!(h, body, center, frame, 1, first, last, segid, data, ctx)
    return nothing
end

"""
    spkw21(handle, body, center, frame, first, last, segid, dlines, epochs; ctx=default_context())

Write an SPK segment of type 21 (extended modified difference arrays), as `spkw21` of
SPICELIB (CSPICE has no C wrapper). `dlines` holds `n` difference lines of equal size
`dlsize` (a matrix with one line in each column, or a vector of vectors); `n` and
`dlsize` follow from its shape.

Preconditions (each throws `SpiceError`): `dlsize <= 111` (`:difflinetoolarge`),
`dlsize >= 71` (`:difflinetoosmall`), a known frame (`:invalidrefframe`), a valid
`segid` (`:segidtoolong`, `:nonprintablechars`), `n >= 1` (`:invalidcount`),
`first <= last` (`:baddescrtimes`), strictly increasing epochs (`:timesoutoforder`),
`epochs[n] >= last` (`:coveragegap`), no zero step size inside a line (`:zerostep`),
and the checks of [`spkpds`](@ref).
"""
function spkw21(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, dlines, epochs; ctx::SegContext=default_context())
    who = "spkw21"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    dl, n, dlsize = _seg_flat(who, "dlines", dlines, 0)
    ep = _seg_vec(who, "epochs", epochs, n)
    dlsize > 111 && throw(SpiceError(:difflinetoolarge,
        "$who: the input difference line size is $dlsize, the maximum supported is 111"))
    dlsize < 71 && throw(SpiceError(:difflinetoosmall,
        "$who: the input difference line size is $dlsize, the minimum supported is 71"))
    _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    n < 1 && throw(SpiceError(:invalidcount,
        "$who: the difference line count was $n; the count must be at least one"))
    first > last && throw(SpiceError(:baddescrtimes,
        "$who: the segment start time $first is greater than the segment end time $last"))
    _seg_increasing(who, ep)
    ep[n] < last && throw(SpiceError(:coveragegap,
        "$who: segment has coverage gap: segment end time $last follows last epoch $(ep[n])"))
    maxdim = (dlsize - 11) ÷ 4
    kqmloc = 4 * maxdim + 8
    for i in 1:n                         # invariant: lines 1 .. i-1 have no zero step
        base = (i - 1) * dlsize
        kqmax1 = round(Int, dl[base+kqmloc], RoundNearestTiesAway)
        for j in 2:kqmax1-1
            dl[base+j] == 0.0 && throw(SpiceError(:zerostep, "$who: step size was zero " *
                "at step size vector index $(j - 1) within difference line $i"))
        end
    end
    data = append!(dl, ep)
    _seg_directory!(data, ep, n ÷ 100)
    push!(data, Float64(maxdim), Float64(n))
    _spk_add_pds!(h, body, center, frame, 21, first, last, segid, data, ctx)
    return nothing
end

# ---------------------------------------------------------------------------------------
# Type 19: mini-segments (ESOC/DDID, piecewise)
# ---------------------------------------------------------------------------------------

const _SPK19_PKTSZS = (12, 6, 6)

"""
    spkw19(handle, body, center, frame, first, last, segid, npkts, subtps, degres, packts,
           epochs, ivlbds, sellst; ctx=default_context())

Write an SPK segment of type 19 (mini-segments of the ESOC/DDID types), as `spkw19` of
SPICELIB (CSPICE has no C wrapper; the interval count is `length(npkts)`).
For mini-segment `i`, `npkts[i]` is the packet count, `subtps[i]` the subtype (0
Hermite with 12 values, 1 Lagrange with 6, 2 Hermite with 6) and `degres[i]` the
degree. `packts` is the flat vector of all packets and `epochs` of all epochs, in
order; `ivlbds` holds the `nintvl + 1` interval bounds; `sellst` selects the last
interval when an epoch is on a common bound.

Preconditions (each throws `SpiceError`): `first <= last` (`:baddescrtimes`), a known
frame (`:invalidrefframe`), a valid `segid` (`:segidtoolong`, `:nonprintablechars`),
at least one interval (`:invalidcount`), strictly increasing bounds
(`:boundsoutoforder`), [`first`, `last`] inside the bounds (`:coveragegap`); for each
mini-segment at least two packets (`:toofewpackets`), a subtype in 0:2
(`:invalidsubtype`), a degree in 1:27 (`:invaliddegree`), an even window size
(`:badwindowsize`), strictly increasing epochs (`:timesoutoforder`), and epochs that
cover the interval (`:boundsdisagree`).
"""
function spkw19(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
        last::Real, segid::AbstractString, npkts::AbstractVector{<:Integer},
        subtps::AbstractVector{<:Integer}, degres::AbstractVector{<:Integer}, packts, epochs,
        ivlbds, sellst::Bool; ctx::SegContext=default_context())
    who = "spkw19"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    nintvl = length(npkts)
    (length(subtps) == nintvl && length(degres) == nintvl) || throw(SpiceError(:badarg,
        "$who: `npkts`, `subtps` and `degres` must have the same length"))
    pk = _seg_vec(who, "packts", packts, -1)
    ep = _seg_vec(who, "epochs", epochs, -1)
    bd = _seg_vec(who, "ivlbds", ivlbds, nintvl < 1 ? -1 : nintvl + 1)
    last < first && throw(SpiceError(:baddescrtimes, "$who: segment start time is " *
        "$first; stop time is $last; bounds must be in nondecreasing order"))
    refcod = _seg_refcode(who, ctx, frame)
    _seg_idlen(who, segid)
    _seg_idchars(who, segid)
    nintvl < 1 && throw(SpiceError(:invalidcount,
        "$who: mini-segment/interval count was $nintvl; this count must be positive"))
    for i in 1:nintvl
        bd[i] >= bd[i+1] && throw(SpiceError(:boundsoutoforder, "$who: interval bounds at " *
            "indices $i and $(i + 1) are $(bd[i]) and $(bd[i+1]); they must be strictly increasing"))
    end
    (first < bd[1] || last > bd[nintvl+1]) && throw(SpiceError(:coveragegap, "$who: " *
        "segment coverage [$first, $last] must be contained within the union of the " *
        "interpolation intervals [$(bd[1]), $(bd[nintvl+1])]"))
    winsizes = zeros(Int, nintvl)
    eepix = 0
    npk = 0
    for i in 1:nintvl                    # invariant: mini-segments 1 .. i-1 are valid
        npkts[i] < 2 && throw(SpiceError(:toofewpackets, "$who: at least 2 packets are " *
            "required for SPK type 19; $(npkts[i]) were supplied in mini-segment $i"))
        subtyp = Int(subtps[i])
        (subtyp < 0 || subtyp > 2) && throw(SpiceError(:invalidsubtype,
            "$who: unexpected SPK type 19 subtype $subtyp found in mini-segment $i"))
        winsiz = subtyp == 1 ? degres[i] + 1 : (degres[i] + 1) ÷ 2
        (degres[i] < 1 || degres[i] > 27) && throw(SpiceError(:invaliddegree, "$who: the " *
            "interpolating polynomials of mini-segment $i have degree $(degres[i]); the " *
            "valid degree range is [1, 27]"))
        isodd(winsiz) && throw(SpiceError(:badwindowsize, "$who: the interpolating " *
            "polynomials of mini-segment $i have window size $winsiz and degree $(degres[i])"))
        winsizes[i] = winsiz
        bepix = eepix + 1
        eepix = bepix - 1 + Int(npkts[i])
        npk += Int(npkts[i]) * _SPK19_PKTSZS[subtyp+1]
        (length(ep) >= eepix && length(pk) >= npk) || throw(SpiceError(:badarg,
            "$who: `epochs` or `packts` is shorter than the packet counts require"))
        for k in bepix:eepix-1
            ep[k] >= ep[k+1] && throw(SpiceError(:timesoutoforder, "$who: in mini-segment " *
                "$i, epoch $(ep[k]) having index $k is greater than or equal to its successor"))
        end
        ep[bepix] > bd[i] && throw(SpiceError(:boundsdisagree, "$who: interpolation interval " *
            "$i start time $(bd[i]) precedes mini-segment's first epoch $(ep[bepix])"))
        ep[eepix] < bd[i+1] && throw(SpiceError(:boundsdisagree, "$who: interpolation interval " *
            "$i end time $(bd[i+1]) exceeds mini-segment's last epoch $(ep[eepix])"))
    end
    data = Float64[]
    pktend = 0
    eepix = 0
    pointers = Float64[]
    segend = 0
    for i in 1:nintvl                    # invariant: mini-segments 1 .. i-1 are in data
        pktsiz = _SPK19_PKTSZS[subtps[i]+1]
        n = Int(npkts[i])
        append!(data, view(pk, pktend+1:pktend+n*pktsiz))
        pktend += n * pktsiz
        bepix = eepix + 1
        eepix = bepix - 1 + n
        append!(data, view(ep, bepix:eepix))
        for j in 1:(n-1)÷100
            push!(data, ep[bepix-1+100j])
        end
        push!(data, Float64(subtps[i]), Float64(winsizes[i]), Float64(n))
        minisz = n * (pktsiz + 1) + (n - 1) ÷ 100 + 3
        push!(pointers, Float64(segend + 1))
        segend += minisz
    end
    append!(data, bd)
    for i in 1:nintvl÷100
        push!(data, bd[100i])
    end
    append!(data, pointers)
    push!(data, Float64(segend + 1), sellst ? 1.0 : -1.0, Float64(nintvl))
    _spk_add!(h, body, center, refcod, 19, first, last, segid, data)
    return nothing
end
