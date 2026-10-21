# Two-line element sets: `getelm` (parse a pair of lines, as `getelm.f` and
# `zzgetelm.f`), `evsgp4` (state in the TEME frame from one element set, as
# `evsgp4.f`) and `spkw10` (write an SPK segment of type 10, as `spkw10.f` with the
# generic-segment writer of `sgseqw.f`).
#
# The geophysical constants have the convention of CSPICE: a vector of eight numbers
# in the order J2, J3, J4, KE, QO, SO, ER, AE. The kernel `geophysical.ker` of the
# toolkit gives them as `BODY399_J2` ... `BODY399_AE`, so that
#
#     geophs = [bodvrd(399, name)[1] for name in ("J2", "J3", "J4", "KE", "QO", "SO", "ER", "AE")]
#
# reads them, as the example of `evsgp4_c` does with `bodvcd_c`.
#
# `getelm` and `evsgp4` need the leap seconds of a kernel set, as in CSPICE.

"Names of the eight geophysical constants in the order of `geophs` (the items of `BODY399_<name>`)."
const _TLE_GEOPHS_NAMES = ("J2", "J3", "J4", "KE", "QO", "SO", "ER", "AE")

"Powers of ten for the exponents of NDD60 and BSTAR, which have one digit."
const _TLE_POWER = (1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9)

"The terms of a two-line element set in the order in which `zzgetelm` parses them: name, line, first and last character."
const _TLE_TERMS = (("YEAR", 1, 19, 20), ("DAY", 1, 21, 32), ("NDT20", 1, 34, 43),
    ("NDD60", 1, 45, 45), ("IEXP", 1, 51, 52), ("BSTAR", 1, 54, 54), ("IBEXP", 1, 60, 61),
    ("INCL", 2, 9, 16), ("NODE0", 2, 18, 25), ("ECC", 2, 27, 33), ("OMEGA", 2, 35, 42),
    ("MO", 2, 44, 51), ("NO", 2, 53, 63))

@noinline _tle_error(msg::String) = throw(SpiceError(:badtle, "Error in TLE set. " * msg))

"`s` with blanks added at the end so that it has at least `n` code units (a Fortran string of length `n`)."
_tle_pad(s::AbstractString, n::Int) = ncodeunits(s) >= n ? String(s) : rpad(String(s), n)

"""
    _tle_number(text, term::Int, integer::Bool) -> Float64

Value of the text of term number `term` of `_TLE_TERMS`, parsed as `nparsd` does (as
`nparsi` if `integer`: the value is cut to a whole number and must fit a 32-bit
integer). Throws `SpiceError(:badtle)` with the name of the term if the text is not
a number.
"""
function _tle_number(text::AbstractString, term::Int, integer::Bool)
    name, line, first, last = _TLE_TERMS[term]
    where = "\"$name\" (characters $first to $last of the " *
            (line == 1 ? "first" : "second") * " line of a two-line element set)"
    x = try
        _nparsd(text)
    catch err
        err isa SpiceError || rethrow()
        _tle_error("An error occurred while trying to parse the term $where. The " *
                   "diagnostic was:  $(err.msg) ")
    end
    integer || return x
    abs(trunc(x)) <= 2147483647.0 || _tle_error("An error occurred while trying to parse " *
        "the term $where. The diagnostic was:  NPARSI: Value entered is beyond the " *
        "bounds of representable integers. ")
    return trunc(x)
end

"""
    getelm([ctx | ks,] frstyr, lines) -> (epoch::Float64, elems::SVector{10,Float64})

Parse the two lines `lines` of a two-line element set, as `getelm_c`, and return the epoch
of the elements (seconds past J2000 TDB) and the ten elements.

`frstyr` is the first year of the 100-year interval in which the two-digit year of the
set is placed. `elems` has the order and the units that [`evsgp4`](@ref) and
[`spkw10`](@ref) take:

| Index | Element | Unit            | Index | Element | Unit                    |
|:------|:--------|:----------------|:------|:--------|:------------------------|
| 1     | NDT20   | radians/minute² | 6     | ECC     |                         |
| 2     | NDD60   | radians/minute³ | 7     | OMEGA   | radians                 |
| 3     | BSTAR   |                 | 8     | M0      | radians                 |
| 4     | INCL    | radians         | 9     | N0      | radians/minute          |
| 5     | NODE0   | radians         | 10    | EPOCH   | seconds past J2000 TDB  |

The context `ctx` (default: the default context) or the kernel set `ks` must hold a
leap-second kernel. As in CSPICE N0067, the checksum in column 69 and the line numbers in
column 1 are not read.

Throws `SpiceError`:
- `:badarg` unless `lines` holds two lines;
- `:noleapseconds` if no leap-second kernel is loaded;
- `:badtle`, as CSPICE, in each of these cases:
  - the satellite numbers of the two lines differ;
  - a line does not have 68 or 69 characters (blanks at the end are not counted);
  - a term is not a number, or an exponent has more than one digit;
  - a value is out of range. The ranges are [0, 360) degrees for NODE0, OMEGA and M0,
    [0, 180] degrees for INCL, and [0, 20] revolutions per day for N0.

# Examples
The file names are those of the generic kernels of NAIF; `geophs` holds the constants of
WGS-72 that the kernel `geophysical.ker` gives.
```julia
furnsh("naif0012.tls")
lines = ["1 18123U 87 53  A 87324.61041692 -.00000023  00000-0 -75103-5 0 00675",
         "2 18123  98.8296 152.0074 0014950 168.7820 191.3688 14.12912554 21686"]
epoch, elems = getelm(1950, lines)
geophs = [1.082616e-3, -2.53881e-6, -1.65597e-6, 7.43669161e-2, 120.0, 78.0, 6378.135, 1.0]
state = evsgp4(epoch + 3600.0, geophs, elems)      # km and km/s, frame TEME (see evsgp4)
```

See also [`evsgp4`](@ref), [`spkw10`](@ref).
"""
function getelm(ks::KernelSet, frstyr::Integer, lines)
    length(lines) == 2 || throw(SpiceError(:badarg,
        "getelm: `lines` must hold two lines, it holds $(length(lines))"))
    l = _lsk(ks)
    raw = (String(lines[1]), String(lines[2]))
    (isascii(raw[1]) && isascii(raw[2])) ||
        _tle_error("A line of the TLE pair has a character that is not ASCII.")
    ln = (_tle_pad(raw[1], 69), _tle_pad(raw[2], 69))
    ln[1][2:7] == ln[2][2:7] || _tle_error("Line 1 of the TLE pair tagged with vehicle ID " *
        "$(ln[1][2:7]),  line 2 of TLE pair tagged with vehicle ID $(ln[2][2:7])")
    for k in 1:2          # invariant: the lines before `k` have 68 or 69 characters
        n = _daf_lastnb(raw[k])
        (n == 68 || n == 69) || _tle_error("Line $k of the TLE has incorrect data length. " *
            "Expected length 68 or 69 elements, actual length: $n. TLE line value: $(raw[k])")
    end
    a, b = ln
    # The pieces of the two lines; NDD60, BSTAR and ECC have an implied decimal point.
    text = (a[19:20], a[21:32], a[34:43], a[45:45] * "." * a[46:50], a[51:52],
        a[54:54] * "." * a[55:59], a[60:61], b[9:16], b[18:25], "0." * b[27:33], b[35:42],
        b[44:51], b[53:63])
    isint = (true, false, false, false, true, false, true, false, false, false, false, false,
        false)
    val = ntuple(i -> _tle_number(text[i], i, isint[i]), Val(13))
    yr, day, ndt20, ndd60, nexp, bstar, bexp, incl, node0, ecc, omega, mo, no = val
    abs(nexp) > 9 && _tle_error("NEXP (exponent) not a single digit. Actual value $(Int(nexp))")
    abs(bexp) > 9 && _tle_error("BEXP (exponent) not a single digit. Actual value $(Int(bexp))")
    (node0 < 0.0 || node0 >= 360.0) &&
        _tle_error("NODE0 (RA acend node) expected bounds [0,360). Actual value $node0")
    (omega < 0.0 || omega >= 360.0) &&
        _tle_error("OMEGA (arg periap) expected bounds [0,360). Actual value $omega")
    (mo < 0.0 || mo >= 360.0) &&
        _tle_error("MO (mean anomoly) expected bounds [0,360). Actual value $mo")
    (incl < 0.0 || incl > 180.0) &&
        _tle_error("INCL (inclination) expected bounds [0,180). Actual value $incl")
    (no > 20.0 || no < 0.0) &&
        _tle_error("NO (mean motion) expected bounds (0,20). Actual value $no")
    pow10(n) = n >= 0 ? _TLE_POWER[Int(n)+1] : 1.0 / _TLE_POWER[1-Int(n)]
    ndd60 *= pow10(nexp)
    bstar *= pow10(bexp)
    d2r = Float64(pi) / 180.0
    pi2 = 2.0 * Float64(pi)
    node0 *= d2r
    omega *= d2r
    mo *= d2r
    incl *= d2r
    no = no * pi2 / 1440.0
    ndt20 = ndt20 * pi2 / 1440.0 / 1440.0
    ndd60 = ndd60 * pi2 / 1440.0 / 1440.0 / 1440.0
    # The year of the 100-year interval that starts at FRSTYR.
    first = Int(frstyr)
    year = div(first, 100) * 100 + Int(yr)
    year < first && (year += 100)
    epoch = _ttrans(l, _TT_YDD, _TT_TDB, _tvec(year, day))[1]
    return epoch, SVector{10,Float64}(ndt20, ndd60, bstar, incl, node0, ecc, omega, mo, no, epoch)
end
getelm(ctx::SpiceContext, frstyr::Integer, lines) = getelm(kernelset(ctx), frstyr, lines)
getelm(frstyr::Integer, lines) = getelm(default_context(), frstyr, lines)

@inline function _tle_tuple(who::String, what::String, x, ::Val{N}) where {N}
    length(x) == N || throw(SpiceError(:badarg,
        "$who: `$what` must have $N elements but has $(length(x))"))
    return ntuple(i -> Float64(@inbounds x[firstindex(x)+i-1]), Val(N))
end

"""
    evsgp4([ctx | ks,] et, geophs, elems) -> SVector{6}

Return the state (km, km/s) in the TEME frame (true equator, mean equinox), at the epoch
`et` (seconds past J2000 TDB), of the two-line element set `elems`, as `evsgp4_c`.

`elems` holds the ten elements that [`getelm`](@ref) returns. `geophs` holds eight
geophysical constants of the Earth:

| Index | Name | Meaning                                                        |
|:------|:-----|:---------------------------------------------------------------|
| 1 – 3 | J2, J3, J4 | gravitational harmonics                                  |
| 4     | KE   | square root of GM, in Earth radii^1.5 per minute               |
| 5     | QO   | upper altitude bound of the atmospheric model, km              |
| 6     | SO   | lower altitude bound of the atmospheric model, km              |
| 7     | ER   | equatorial radius, km                                          |
| 8     | AE   | distance units per Earth radius                                |

As in CSPICE, QO, SO and AE are not used. The propagator is SGP4, or SDP4 for an orbital
period of 225 minutes or more.

The context `ctx` (default: the default context) or the kernel set `ks` must hold a
leap-second kernel. `et` can have any real number type; a `ForwardDiff.Dual` epoch gives
the time derivatives of the state. The call does not allocate for arguments that are
tuples, static vectors or vectors.

Throws `SpiceError`:
- `:badarg` unless `geophs` has 8 elements and `elems` has 10;
- `:noleapseconds` if no leap-second kernel is loaded;
- `:suborbital`, `:badmeanmotion`, `:badmeccentricity`, `:badmsemimajor`,
  `:badpeccentricity`, `:badsemilatus` or `:orbitdecay`, as CSPICE, for elements that
  cannot be propagated.

See [`getelm`](@ref) for an example.
"""
function evsgp4(ks::KernelSet, et::Real, geophs, elems)
    g = _tle_tuple("evsgp4", "geophs", geophs, Val(8))
    e = _tle_tuple("evsgp4", "elems", elems, Val(10))
    p, v = _t10_propagate(_lsk(ks), g, e, et)
    return vcat(p, v)
end
evsgp4(ctx::SpiceContext, et::Real, geophs, elems) = evsgp4(kernelset(ctx), et, geophs, elems)
evsgp4(et::Real, geophs, elems) = evsgp4(default_context(), et, geophs, elems)

"""
    spkw10([ctx,] handle, body, center, frame, first, last, segid, consts, elems, epochs)
        -> Nothing

Write an SPK segment of type 10 (two-line element sets), as `spkw10_c`.

`consts` holds the eight geophysical constants of [`evsgp4`](@ref). `elems`
holds `n` element sets of ten numbers each, in the form that [`getelm`](@ref) returns, and
`epochs` holds their `n` epochs (seconds past J2000 TDB). Each epoch must be equal to the
tenth number of its set; as in CSPICE, the call does not check this. See [`spkw02`](@ref)
for the other arguments and for the errors of every writer.

The call adds to each element set the nutation angles of the IAU 1980 series at its epoch
and their rates, as `spkw10_c` does, so a packet of the segment has 14 doubles.

Throws `SpiceError`:
- the error kinds of [`spkpds`](@ref);
- `:segidtoolong` or `:nonprintingchars` if `segid` is not valid;
- `:numpacketsnotpos` if there is no element set;
- `:unorderedrefs` unless the epochs are strictly increasing.
"""
function spkw10(ctx::SegContext, h::DAFHandle, body::Integer, center::Integer,
        frame::SegFrame, first::Real, last::Real, segid::AbstractString, consts, elems,
        epochs)
    who = "spkw10"
    _daf_check_write(h, who)
    frame isa AbstractString && _seg_chkfstr(who, "frame", frame)
    _seg_chkfstr(who, "segid", segid)
    con = _seg_fixed(who, "consts", consts, 8)
    el, n, _ = _seg_flat(who, "elems", elems, 10)
    ep = _seg_vec(who, "epochs", epochs, n)
    descr = spkpds(ctx, body, center, frame, 10, first, last)
    _seg_idlen(who, segid)
    cu = codeunits(String(segid))
    for i in 1:_daf_lastnb(segid)
        (cu[i] < 0x20 || cu[i] > 0x7e) && throw(SpiceError(:nonprintingchars,
            "$who: the segment identifier contains a nonprinting character at location $i"))
    end
    n <= 0 && throw(SpiceError(:numpacketsnotpos,
        "$who: the number of packets to store is not positive: $n"))
    for i in 2:n          # invariant: ep[1:i-1] is strictly increasing
        ep[i-1] >= ep[i] && throw(SpiceError(:unorderedrefs,
            "$who: reference values are out of order; the value $(ep[i-1]) was found to " *
            "be out of order for index $(i - 1)"))
    end
    packets = Vector{Float64}(undef, 14n)
    for i in 1:n          # invariant: the packets of the sets before `i` are in `packets`
        copyto!(packets, 14(i - 1) + 1, el, 10(i - 1) + 1, 10)
        psi, eps, dpsi, deps = _t10_wahr(ep[i])
        packets[14i-3:14i] .= (eps, psi, deps, dpsi)
    end
    dafbna(h, descr, segid)
    dafada(h, _t10_data(con, ep, packets))
    dafena(h)
    return nothing
end
spkw10(h::DAFHandle, body::Integer, center::Integer, frame::SegFrame, first::Real,
    last::Real, segid::AbstractString, consts, elems, epochs) =
    spkw10(default_context(), h, body, center, frame, first, last, segid, consts, elems,
        epochs)

"""
    _t10_data(consts, refs, packets) -> Vector{Float64}

The data of a type 10 segment as `sgbwfs`, `sgwfpk` and `sgwes` of CSPICE write it: the
eight constants, one record of the epoch and the 14 doubles of each packet, the epochs,
every 100th epoch, and the 17 meta data items. `packets` holds the packets one after the
other.
"""
function _t10_data(con::AbstractVector{Float64}, refs::AbstractVector{Float64},
        packets::AbstractVector{Float64})
    n = length(refs)
    @assert length(con) == 8 && length(packets) == 14n
    data = Vector{Float64}(undef, 0)
    sizehint!(data, 8 + 16 * n + (n - 1) ÷ 100 + 17)
    append!(data, con)
    for i in 1:n          # invariant: the records of the packets before `i` are in `data`
        push!(data, refs[i])
        append!(data, @view packets[14(i-1)+1:14i])
    end
    append!(data, refs)
    for i in 1:(n-1)÷100
        push!(data, refs[100i])
    end
    # Meta data: CONBAS, NCON, RDRBAS, NRDR, RDRTYP, REFBAS, NREF, PDRBAS, NPDR, PDRTYP,
    # PKTBAS, NPKT, RSVBAS, NRSV, PKTSZ, PKTOFF, NMETA.
    refbas = 8 + 15 * n
    append!(data, Float64[0, 8, refbas + n, (n - 1) ÷ 100, 4, refbas, n, 0, 0, 0, 8, n, 0, 0,
        14, 1, 17])
    return data
end
