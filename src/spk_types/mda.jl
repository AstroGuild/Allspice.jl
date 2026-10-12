# Modified difference arrays: type 1 (difference lines of 71 doubles, order up to 15)
# and type 21 (extended difference lines, order up to 25).

"Largest order of a difference line (`MAXTRM` of `spke21`); type 1 has 15."
const _MDA_MAXDIM = 25

"Sum of `data[a+j] * w[j+ks]` for `j = kqq, ..., 1` (in this order, as `spke01`); requires the indices in range (not checked)."
@inline function _mda_sum(data::Vector{Float64}, a::Int, kqq::Int, w::MVector{L,T}, ks::Int) where {L,T}
    sum = zero(T)
    @inbounds for j in kqq:-1:1       # invariant: sum holds the terms kqq .. j+1
        sum += data[a+j] * w[j+ks]
    end
    return sum
end

"""
    _spk_mda(data, seg, et::T) -> (pos, vel)

State of a type 1 or type 21 segment, as `spkr01`/`spke01` and `spkr21`/`spke21`.

The segment holds the difference lines, the final epoch of each line, an epoch
directory and the number of lines; a type 21 segment has the order `MAXDIM` of its
lines before that number. The line is the first one whose final epoch is not before
the epoch. A line of order `m` holds the reference epoch TL, `m` step sizes G, the
reference position and velocity (interleaved), the `m` by 3 difference table, the
order KQMAX1 and the three orders KQ; position and velocity are the reference state
plus sums of the differences times coefficients that depend on `et - TL` and G.
"""
function _spk_mda(data::Vector{Float64}, seg::Segment, et::T) where {T<:Real}
    b, e = seg.baddr, seg.eaddr
    len = e - b + 1
    trailer = seg.dtype == 1 ? 1 : 2
    len > trailer || _badsegment(seg, "the segment is too short")
    m = seg.dtype == 1 ? 15 :
        _seg_count(seg, @inbounds(data[e-1]), 1, _MDA_MAXDIM, "difference line order")
    nrec = _seg_count(seg, @inbounds(data[e]), 1, _MAXCOUNT, "difference line count")
    dlsize = 4 * m + 11
    nrec * (dlsize + 1) + nrec ÷ 100 + trailer == len ||
        _badsegment(seg, "the difference line count does not agree with the length")
    ebase = b - 1 + nrec * dlsize               # final epoch i is data[ebase+i]
    recno = min(_lstltd(data, ebase, nrec, primal(et)) + 1, nrec)
    r = b + (recno - 1) * dlsize                # address of TL
    @inbounds begin
        tl = data[r]
        # The arrays of `spke01` allow KQMAX1 up to 15, the ones of `spke21` up to m+1.
        kqmax1 = _seg_count(seg, data[r+4m+7], 1, seg.dtype == 1 ? m : m + 1, "KQMAX1")
        kq1 = _seg_count(seg, data[r+4m+8], 0, m, "KQ(1)")
        kq2 = _seg_count(seg, data[r+4m+9], 0, m, "KQ(2)")
        kq3 = _seg_count(seg, data[r+4m+10], 0, m, "KQ(3)")
    end
    (kqmax1 >= 2 || max(kq1, kq2, kq3) == 0) ||
        _badsegment(seg, "KQMAX1 = 1 with a difference order above 0")
    g = r                                       # G(j) is data[g+j]
    ref = r + m                                 # position i is data[ref+2i-1], velocity data[ref+2i]
    dtb = r + m + 6                             # DT(j,i) is data[dtb+(i-1)m+j]

    delta = et - tl
    fc = MVector{_MDA_MAXDIM,T}(undef)
    wc = MVector{_MDA_MAXDIM,T}(undef)
    w = MVector{_MDA_MAXDIM + 2,T}(undef)
    tp = delta
    mq2 = kqmax1 - 2
    ks = kqmax1 - 1
    @inbounds begin
        fc[1] = one(T)
        for j in 1:mq2
            gj = data[g+j]
            gj == 0.0 && _badsegment(seg, "a step size of difference line $recno is zero")
            fc[j+1] = tp / gj
            wc[j] = delta / gj
            tp = delta + gj
        end
        for j in 1:_MDA_MAXDIM+2
            w[j] = j <= kqmax1 ? one(T) / j : zero(T)
        end
        # Each pass lowers by one the order of the divided differences in `w`; after
        # the loop `w[j+1]` is the coefficient of difference j of the position sum.
        jx = 0
        ks1 = ks - 1
        while ks >= 2
            jx += 1
            for j in 1:jx
                w[j+ks] = fc[j+1] * w[j+ks1] - wc[j] * w[j+ks]
            end
            ks = ks1
            ks1 -= 1
        end
        s1 = _mda_sum(data, dtb, kq1, w, ks)
        s2 = _mda_sum(data, dtb + m, kq2, w, ks)
        s3 = _mda_sum(data, dtb + 2m, kq3, w, ks)
        pos = Vec3{T}(data[ref+1] + delta * (data[ref+2] + delta * s1),
            data[ref+3] + delta * (data[ref+4] + delta * s2),
            data[ref+5] + delta * (data[ref+6] + delta * s3))
        # One more pass gives the coefficients of the velocity sum in `w[j]`.
        for j in 1:jx
            w[j+ks] = fc[j+1] * w[j+ks1] - wc[j] * w[j+ks]
        end
        ks -= 1
        s1 = _mda_sum(data, dtb, kq1, w, ks)
        s2 = _mda_sum(data, dtb + m, kq2, w, ks)
        s3 = _mda_sum(data, dtb + 2m, kq3, w, ks)
        vel = Vec3{T}(data[ref+2] + delta * s1, data[ref+4] + delta * s2,
            data[ref+6] + delta * s3)
    end
    return pos, vel
end
