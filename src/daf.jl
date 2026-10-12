# Reader of the DAF container (binary SPK and binary PCK kernels).
#
# A DAF file is a sequence of 1024-byte records. Record 1 is the file record. Summary
# records form a doubly linked list that starts at record FWARD; each holds three
# control doubles (next, previous, count) and then the packed segment summaries. The
# data of a segment are doubles addressed from 1 at the start of the file.

const _DAF_RECLEN = 1024
const _DAF_RECWORDS = 128
# Characters that a text-mode transfer would damage; a valid file holds them unchanged.
const _DAF_FTPSTR = UInt8[codeunits("FTPSTR:\r:\n:\r\n:\r")..., 0x00, UInt8(':'), 0x81,
    UInt8(':'), 0x10, 0xce, codeunits(":ENDFTP")...]

"32-bit integer at 0-based byte offset `off` of `raw`, stored little-endian if `little`, else big-endian; requires `off + 4 <= length(raw)`."
@inline function _daf_i32(raw::AbstractVector{UInt8}, off::Int, little::Bool)
    v = UInt32(raw[off+1]) | (UInt32(raw[off+2]) << 8) | (UInt32(raw[off+3]) << 16) |
        (UInt32(raw[off+4]) << 24)
    return reinterpret(Int32, little ? v : bswap(v))
end

"Double at 0-based byte offset `off` of `raw`, stored little-endian if `little`, else big-endian; requires `off + 8 <= length(raw)`."
@inline function _daf_f64(raw::AbstractVector{UInt8}, off::Int, little::Bool)
    v = UInt64(0)
    for i in 1:8                      # invariant: v holds bytes off+1 .. off+i-1, little-endian
        v |= UInt64(raw[off+i]) << (8 * (i - 1))
    end
    return reinterpret(Float64, little ? v : bswap(v))
end

"Value of `x` as an `Int` if `x` is a finite whole number in `lo:hi`, else throws `SpiceError(:baddaf)` that names `what`."
function _daf_int(x::Float64, lo::Int, hi::Int, what::AbstractString, path::AbstractString)
    (isfinite(x) && x == trunc(x) && lo <= x <= hi) ||
        throw(SpiceError(:baddaf, "$what = $x is not a whole number in $lo:$hi in $path"))
    return Int(x)
end

"""
    _daf_directory(data, path, kind, dtype, baddr, eaddr) -> (init, intlen, rsize, nrec)

Directory values (the last four doubles) of a fixed-interval Chebyshev segment that
occupies `data[baddr:eaddr]`; throws `SpiceError(:baddaf)` unless RSIZE and N are
positive whole numbers that agree with the segment length and the data type.
"""
function _daf_directory(data::Vector{Float64}, path::AbstractString, kind::Symbol,
        dtype::Int32, baddr::Int, eaddr::Int)
    len = eaddr - baddr + 1
    # SPK type 20 ends with DSCALE, TSCALE, INITJD, INITFR, INTLEN, RSIZE, N.
    trailer = (kind === :SPK && dtype == 20) ? 7 : 4
    len > trailer || throw(SpiceError(:baddaf,
        "segment of type $dtype at $baddr:$eaddr is too short ($len doubles) in $path"))
    init = data[eaddr-3]
    intlen = data[eaddr-2]
    (isfinite(init) && isfinite(intlen) && intlen > 0) || throw(SpiceError(:baddaf,
        "segment at $baddr:$eaddr has INIT = $init, INTLEN = $intlen in $path"))
    rsize = _daf_int(data[eaddr-1], 1, Int(typemax(Int32)), "RSIZE", path)
    nrec = _daf_int(data[eaddr], 1, Int(typemax(Int32)), "N", path)
    Int128(rsize) * nrec + trailer == len || throw(SpiceError(:baddaf,
        "segment at $baddr:$eaddr: RSIZE = $rsize and N = $nrec do not agree with " *
        "the length $len in $path"))
    if dtype == 2 || dtype == 3
        per = dtype == 2 ? 3 : 6      # coefficient sets in one record
        (rsize > 2 && (rsize - 2) % per == 0) || throw(SpiceError(:baddaf,
            "segment at $baddr:$eaddr: RSIZE = $rsize is not 2 + $per*NCOEF in $path"))
    end
    return init, intlen, Int32(rsize), Int32(nrec)
end

"""
    open_daf(path::AbstractString) -> DAF

Open the binary SPK or binary PCK kernel at `path` and return its data array and its
segment descriptors in file order.

Contract: `path` names a readable DAF file. In the result `data[a]` is the double at
DAF address `a` (memory-mapped for a file of native byte order, an owned byte-swapped
copy otherwise), `kind` is `:SPK` or `:PCK`, and every segment has finite epochs and
`1 <= baddr <= eaddr <= length(data)`. For SPK types 2, 3, 20 and PCK type 2 the
fields `init`, `intlen`, `rsize`, `nrec` hold the last four doubles of the segment and
`rsize`, `nrec` agree with the segment length. Throws `SpiceError` with kind
`:nofile` or `:baddaf` if the file does not meet these conditions.

The ID words `DAF/SPK` and `DAF/PCK` are accepted, and also the old `NAIF/DAF` (as in
`de405.bsp`), for which the kind follows from NI (6 is SPK, 5 is PCK).
"""
function open_daf(path::AbstractString)
    isfile(path) || throw(SpiceError(:nofile, "no such file: $path"))
    nbytes = filesize(path)
    nbytes >= 2 * _DAF_RECLEN || throw(SpiceError(:baddaf,
        "file of $nbytes bytes is too short to be a DAF: $path"))
    nwords = Int(nbytes ÷ 8)
    nrecords = Int(nbytes ÷ _DAF_RECLEN)

    return open(path, "r") do io
        raw = Mmap.mmap(io, Vector{UInt8}, (8 * nwords,))

        # --- file record ---------------------------------------------------------
        idword = rstrip(String(raw[1:8]))
        idword in ("DAF/SPK", "DAF/PCK", "NAIF/DAF") || throw(SpiceError(:baddaf,
            "ID word $(repr(idword)) is not DAF/SPK, DAF/PCK or NAIF/DAF in $path"))
        bfmt = rstrip(String(raw[89:96]))
        little = if bfmt == "LTL-IEEE"
            true
        elseif bfmt == "BIG-IEEE"
            false
        else
            throw(SpiceError(:baddaf, "unknown byte-order string $(repr(bfmt)) in $path"))
        end
        swap = little != (ENDIAN_BOM == 0x04030201)
        if raw[700:705] == _DAF_FTPSTR[1:6]      # the string is absent in old files
            raw[700:727] == _DAF_FTPSTR || throw(SpiceError(:baddaf,
                "FTP validation string is damaged (text-mode transfer?) in $path"))
        end
        nd = Int(_daf_i32(raw, 8, little))
        ni = Int(_daf_i32(raw, 12, little))
        fward = Int(_daf_i32(raw, 76, little))
        kind = idword == "DAF/SPK" ? :SPK : idword == "DAF/PCK" ? :PCK :
               ni == 6 ? :SPK : ni == 5 ? :PCK : :UNKNOWN
        want_ni = kind === :SPK ? 6 : 5
        (kind !== :UNKNOWN && nd == 2 && ni == want_ni) || throw(SpiceError(:baddaf,
            "ND = $nd, NI = $ni do not fit ID word $(repr(idword)) in $path"))
        1 < fward <= nrecords || throw(SpiceError(:baddaf,
            "first summary record $fward is outside the file ($nrecords records): $path"))

        # --- data array -----------------------------------------------------------
        data = if swap
            d = Vector{Float64}(undef, nwords)
            for a in 1:nwords                    # invariant: d[1:a-1] is converted
                d[a] = _daf_f64(raw, 8 * (a - 1), little)
            end
            d
        else
            Mmap.mmap(io, Vector{Float64}, (nwords,))
        end

        # --- summary records ------------------------------------------------------
        ss = nd + (ni + 1) ÷ 2                   # doubles in one summary
        maxsum = (_DAF_RECWORDS - 3) ÷ ss
        segments = Segment[]
        rec = fward
        visited = 0
        while rec != 0
            # invariant: `visited` summary records were read, and none twice while
            # `visited <= nrecords`
            visited += 1
            visited <= nrecords || throw(SpiceError(:baddaf,
                "the chain of summary records has a cycle in $path"))
            base = (rec - 1) * _DAF_RECLEN       # 0-based byte offset of the record
            next = _daf_int(data[(rec-1)*_DAF_RECWORDS+1], 0, nrecords,
                "next summary record", path)
            nsum = _daf_int(data[(rec-1)*_DAF_RECWORDS+3], 0, maxsum,
                "summary count of record $rec", path)
            for i in 1:nsum
                off = base + 8 * (3 + (i - 1) * ss)
                t0 = _daf_f64(raw, off, little)
                t1 = _daf_f64(raw, off + 8, little)
                (isfinite(t0) && isfinite(t1)) || throw(SpiceError(:baddaf,
                    "segment $(length(segments) + 1) has non-finite epochs in $path"))
                ic = ntuple(j -> _daf_i32(raw, off + 16 + 4 * (j - 1), little), 6)
                target, center, frame, dtype, b, e = if kind === :SPK
                    ic[1], ic[2], ic[3], ic[4], Int(ic[5]), Int(ic[6])
                else
                    ic[1], Int32(0), ic[2], ic[3], Int(ic[4]), Int(ic[5])
                end
                (1 <= b <= e <= nwords) || throw(SpiceError(:baddaf,
                    "segment $(length(segments) + 1) has addresses $b:$e outside " *
                    "the file ($nwords doubles): $path"))
                chebyshev = kind === :SPK ? (dtype == 2 || dtype == 3 || dtype == 20) :
                            dtype == 2
                init, intlen, rsize, nrec = chebyshev ?
                    _daf_directory(data, path, kind, dtype, b, e) :
                    (0.0, 0.0, Int32(0), Int32(0))
                push!(segments, Segment(target, center, frame, dtype, t0, t1, Int32(0),
                    b, e, init, intlen, rsize, nrec))
            end
            rec = next
        end
        DAF(String(path), kind, data, segments)
    end
end
