# Chebyshev types with their own record layout: type 14 (unequal intervals, a generic
# segment) and type 20 (velocity series only, position by integration).

"Julian date of J2000 and seconds in one day (`j2000` and `spd` of CSPICE)."
const _J2000_JD = 2451545.0
const _SPD = 86400.0

"""
    _spk14(data, seg, et::T) -> (pos, vel)

State of a type 14 segment, as `spkr14` and `spke14`. The segment is a generic
segment: constants (one, the number of coefficients of a series), packets of fixed
size (`MID`, `RADIUS`, three position series, three velocity series), the start
epochs of the packets as reference values, a directory of them, and the meta data
that give the place and the count of each part. The packet is the last one whose
start epoch is not after the epoch.

Only the layout that `spk14b` of CSPICE writes is accepted (fixed packet size,
explicit reference values with the index type "last less than or equal"); another
one throws `SpiceError(:notsupported)`.
"""
function _spk14(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    lay = _spk14_layout(data, seg)
    b, ncoef, stride = seg.baddr, lay.ncoef, lay.stride
    k = clamp(_lstled(data, b - 1 + lay.refbas, lay.nref, primal(et)), 1, lay.npkt)
    a = b + lay.pktbas + (k - 1) * stride + lay.pktoff  # address of MID
    @inbounds mid, radius = data[a], data[a+1]
    tau = (et - mid) / radius
    (pos,) = _cheb3(data, a + 2, ncoef, tau, 1 / radius, Val(0))
    (vel,) = _cheb3(data, a + 2 + 3 * ncoef, ncoef, tau, 1 / radius, Val(0))
    return pos, vel
end

"""
    _spk14_layout(data, seg) -> NamedTuple

The parts of the type 14 segment `seg` from its meta data, as offsets from the address
before the segment: `conbas`, `refbas`, `pktbas`, the counts `nref` and `npkt`, the
packet size `pktsz`, the packet offset `pktoff`, `stride = pktsz + pktoff`, and the
coefficient count `ncoef` of a series. Throws `SpiceError(:notsupported)` for a layout
that `spk14b` does not write and `SpiceError(:baddaf)` for one that does not fit.
"""
@inline function _spk14_layout(data::Vector{Float64}, seg::Segment)
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    nmeta = _seg_count(seg, @inbounds(data[e]), 15, len, "meta data count")
    # A segment of 15 meta items was written as 16 doubles (`sgmeta`).
    mb = e - max(nmeta, 16)                     # meta item i is data[mb+i]
    mb >= b - 1 || _badsegment(seg, "the meta data do not fit the segment")
    @inbounds begin
        conbas = _seg_count(seg, data[mb+1], 0, len, "constants base")
        ncon = _seg_count(seg, data[mb+2], 1, len, "constants count")
        rdrtyp = _seg_count(seg, data[mb+5], 0, 4, "reference directory type")
        refbas = _seg_count(seg, data[mb+6], 0, len, "reference base")
        nref = _seg_count(seg, data[mb+7], 1, len, "reference count")
        pdrtyp = _seg_count(seg, data[mb+10], 0, 1, "packet directory type")
        pktbas = _seg_count(seg, data[mb+11], 0, len, "packet base")
        npkt = _seg_count(seg, data[mb+12], 1, len, "packet count")
        pktsz = _seg_count(seg, data[mb+15], 8, len, "packet size")
        pktoff = nmeta >= 16 ? _seg_count(seg, data[mb+16], 0, len, "packet offset") : 0
    end
    (rdrtyp == 3 && pdrtyp == 0) || throw(SpiceError(:notsupported,
        "SPK type 14 segment with reference directory type $rdrtyp and packet " *
        "directory type $pdrtyp is not supported"))
    stride = pktsz + pktoff
    (nref == npkt && refbas + nref <= len && Int128(pktbas) + Int128(npkt) * stride <= len) ||
        _badsegment(seg, "the meta data do not agree with the length")
    ncoef = _seg_count(seg, @inbounds(data[b-1+conbas+1]), 1, len, "coefficient count")
    6 * ncoef + 2 == pktsz || _badsegment(seg, "the packet size is not 2 + 6*NCOEF")
    return (; conbas, refbas, nref, pktbas, npkt, pktsz, pktoff, stride, ncoef)
end

"""
    _chbigr(data, a, ncof, scale, mid, radius, x::T) -> (p::T, itgrl::T)

Value at `x` of the Chebyshev series with the `ncof` coefficients
`data[a:a+ncof-1] * scale` on the interval `mid ± radius`, and its integral from
`mid` to `x`, by the recurrences of `chbigr`; also the kernel of the public `chbigr`.
Requires `ncof >= 1`, `radius > 0` and the coefficients inside `data` (not checked).
"""
@inline function _chbigr(data::AbstractVector, a::Int, ncof::Int, scale::Float64, mid::Real,
        radius::Real, x::T) where {T}
    cp(i) = @inbounds data[a+i-1] * scale       # coefficient i, 1-based
    degp = ncof - 1
    s = (x - mid) / radius
    s2 = s * 2.0
    # Coefficients of the integral: a_i = (c_{i-1} - c_{i+1}) / (2 (i-1)), with the
    # special cases of the ends.
    a2 = ncof >= 3 ? cp(1) - cp(3) * 0.5 : cp(1)
    adegp1 = degp >= 2 ? cp(degp) * 0.5 / degp : zero(a2)
    adegp2 = degp >= 1 ? cp(degp + 1) * 0.5 / (degp + 1) : zero(a2)
    f0 = zero(T) + (degp == 0 ? a2 : adegp2)
    f1 = zero(T)
    z0, z1 = (degp == 0 ? a2 : adegp2), zero(a2)    # the integral series at s = 0
    w0 = zero(T)
    w1 = zero(T)
    i = ncof
    # invariant: f, z and w hold the Clenshaw sums of the terms above i-1 of the
    # integral series at s, of the integral series at 0, and of the series itself
    while i > 1
        ai = i == 2 ? a2 : i < ncof ? (cp(i - 1) - cp(i + 1)) * 0.5 / (i - 1) : adegp1
        f0, f1 = ai + (s2 * f0 - f1), f0
        z0, z1 = ai - z1, z0
        w0, w1 = cp(i) + (s2 * w0 - w1), w0
        i -= 1
    end
    itgrl = z1 + s * f0 - f1
    p = cp(1) + (s * w0 - w1)
    return p, radius * itgrl
end

"""
    _spk20(data, seg, et::T) -> (pos, vel)

State of a type 20 segment, as `spkr20` and `spke20`. The segment ends with DSCALE,
TSCALE, INITJD, INITFR, INTLEN (days), RSIZE and N; `open_daf` checks RSIZE and N and
puts them in the descriptor. A record holds, for each component, the coefficients
of the velocity series (in units of DSCALE km per TSCALE s) and then the position at
the middle of the record interval (in units of DSCALE km). The velocity is the
series; the position is the middle position plus the integral of the series.
"""
function _spk20(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    b, e = seg.baddr, seg.eaddr
    recsiz, nrec = Int(seg.rsize), Int(seg.nrec)
    (nrec >= 1 && recsiz >= 1 && e - b + 1 == Int128(recsiz) * nrec + 7) ||
        _badsegment(seg, "the record size and count do not agree with the length")
    nterms = recsiz ÷ 3
    ncof = nterms - 1
    (ncof >= 1 && 3 * nterms == recsiz) || _badsegment(seg, "the record size $recsiz is not valid")
    @inbounds dscale, tscale, initjd, initfr, intlen = data[e-6], data[e-5], data[e-4],
        data[e-3], data[e-2]
    (isfinite(dscale) && isfinite(tscale) && tscale != 0.0 && isfinite(initjd) &&
     isfinite(initfr) && isfinite(intlen) && intlen > 0.0) ||
        _badsegment(seg, "the scales or the interval values are not valid")
    init = (initjd - _J2000_JD + initfr) * _SPD
    intrvl = intlen * _SPD
    recno = clamp(_trunc_index((primal(et) - init) / intrvl) + 1, 1, nrec)
    recbeg = (initjd - _J2000_JD + (recno - 1) * intlen) * _SPD
    radius = intrvl / 2.0
    mid = recbeg + initfr * _SPD + radius
    r = b + (recno - 1) * recsiz
    vscale = dscale / tscale
    vx, ix = _chbigr(data, r, ncof, vscale, mid, radius, et)
    vy, iy = _chbigr(data, r + nterms, ncof, vscale, mid, radius, et)
    vz, iz = _chbigr(data, r + 2 * nterms, ncof, vscale, mid, radius, et)
    @inbounds px, py, pz = dscale * data[r+nterms-1], dscale * data[r+2*nterms-1],
        dscale * data[r+3*nterms-1]
    return Vec3{T}(px + ix, py + iy, pz + iz), Vec3{T}(vx, vy, vz)
end
