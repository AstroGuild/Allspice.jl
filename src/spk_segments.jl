# SPK segments named by a file and a descriptor: the search of the loaded segments
# (`spksfs`), the evaluation of one segment (`spkpvn`), the unpacking of a descriptor
# (`spkuds`), and the subset of a segment written to another file (`spksub`, with the
# type-specific routines `spksNN` of CSPICE as `_spksNN`).
#
# A loaded file has no handle in this package (`kdata` gives 0), so `spksfs` returns the
# path of the loaded file in the place of the handle, and `spkpvn` takes that path with
# the context that holds the file. `spkpvn` and `spksub` also take a `DAFHandle`.
#
# The `_spksNN` functions follow the C code line by line: the same reads, the same
# integer conversions (`_cint` for a C cast, `_nint` for `i_dnnt`) and the same index
# arithmetic, so that the new segment is equal double for double to the one of CSPICE.
# The early returns of the C routines for `begin > end` are left out: `spksub` has
# checked `begin <= end` before.

# ---------------------------------------------------------------------------------------
# Descriptors
# ---------------------------------------------------------------------------------------

"""
    spkuds(descr) -> (body, center, frame, type, first, last, baddr, eaddr)

Unpack the SPK segment descriptor `descr`, as `spkuds_c`.

`body` moves relative to `center` in the frame with the ID `frame`, and `type` is the SPK
data type. `first` and `last` bound the coverage, in seconds past J2000 TDB. `baddr` and
`eaddr` are the first and the last DAF address of the segment data. As in CSPICE, the
values are not checked. The integers are `Int`.

Throws `SpiceError(:dafinvalidparams)` if `descr` has fewer than 5 doubles.

See also [`spksfs`](@ref), [`spkpds`](@ref), [`dafus`](@ref).
"""
function spkuds(descr::AbstractVector{<:Real})
    dc, ic = dafus(descr, 2, 6)
    return ic[1], ic[2], ic[3], ic[4], dc[1], dc[2], ic[5], ic[6]
end

"The file of a `Segment` with `file == k`: entry `k` among the binary kernels of `ks.loaded`."
function _daf_entry(ks::KernelSet, k::Integer)
    n = 0
    for entry in ks.loaded               # invariant: n binary kernels come before `entry`
        entry.content isa DAF || continue
        n += 1
        n == k && return entry
    end
    @assert false "the segment belongs to a file of the kernel set"
end

"""
    _daf_array_name(daf, k) -> String

The name of array `k`, in file order, of the loaded DAF `daf`, as `dafgn` reads it. The
data of a file of the other byte order hold each double converted, so the bytes of the
file record and of the name records are swapped back first.
"""
function _daf_array_name(daf::DAF, k::Int)
    d = daf.data
    swap = String(reinterpret(UInt8, [d[12]])) != _DAF_NATIVE_BFF
    rawrec(r) = collect(reinterpret(UInt8,
        [swap ? bswap(d[a]) : d[a] for a in _daf_rwa(r):_daf_rwa(r)+_DAF_RECWORDS-1]))
    little = swap ? !_DAF_NATIVE_LITTLE : _DAF_NATIVE_LITTLE
    rec = Int(_daf_i32(rawrec(1), 76, little))
    namsiz = 8 * _daf_sumsize(2, 6)
    count = 0
    while rec != 0                       # invariant: `count` arrays are in the records before `rec`
        base = _daf_rwa(rec) - 1
        nsum = Int(d[base+3])
        if k <= count + nsum
            off = (k - count - 1) * namsiz
            return String(rstrip(String(rawrec(rec + 1)[off+1:off+namsiz]), [' ', '\0']))
        end
        count += nsum
        rec = Int(d[base+1])
    end
    @assert false "array $k is in the summary records of $(daf.path)"
end

"""
    spksfs([ctx,] body::Integer, et) -> Union{Nothing,Tuple{String,Vector{Float64},String}}

Return `(file, descr, ident)` for the loaded SPK segment of the highest priority for
`body` at `et`, as `spksfs_c`, or `nothing` if no loaded segment of `body` covers `et`.

- `file`: the path of the loaded SPK file as it was passed to [`furnsh`](@ref) (the
  path that [`kdata`](@ref) gives). It takes the place of the handle of CSPICE: pass it
  to [`spkpvn`](@ref) with the same context, or open the file with [`dafopr`](@ref).
- `descr`: the descriptor of the segment, 5 doubles; [`spkuds`](@ref) unpacks it.
- `ident`: the segment identifier, without trailing blanks.

`et` is in seconds past J2000 TDB. A segment covers its first and its last epoch. The
priority is that of CSPICE: a file loaded later comes first, and in one file a later
segment comes first. `nothing` takes the place of the `found` flag of CSPICE.

Throws `SpiceError(:noloadedfiles)` if no SPK file is loaded.

# Examples
```jldoctest
julia> path = joinpath(mktempdir(), "example.bsp");

julia> h = spkopn(path, "example SPK");

julia> spkw02(h, 499, 4, "J2000", 0.0, 100.0, "mars", 100.0, [[1.0, 0.5, 2.0, 0.0, 3.0, 0.0]], 0.0)

julia> spkcls(h)

julia> ctx = SpiceContext();

julia> furnsh(ctx, path)

julia> file, descr, ident = spksfs(ctx, 499, 50.0);

julia> ident, spkuds(descr)
("mars", (499, 4, 1, 2, 0.0, 100.0, 385, 396))

julia> ref, state, center = spkpvn(ctx, file, descr, 100.0);

julia> ref, state == [1.5, 2.0, 3.0, 0.01, 0.0, 0.0], center
(1, true, 4)

julia> spksfs(ctx, 499, 200.0) === nothing
true
```

See also [`spksub`](@ref), [`spkobj`](@ref), [`spkcov`](@ref).
"""
function spksfs(ks::KernelSet, body::Integer, et::Real)
    segs = get(ks.spk, Int32(body), nothing)
    i = segs === nothing ? 0 : find_segment(segs, primal(et))
    if i == 0
        any(e -> e.content isa DAF && e.content.kind === :SPK, ks.loaded) ||
            throw(SpiceError(:noloadedfiles, "spksfs: at least one SPK file needs to be " *
                "loaded before a segment search"))
        return nothing
    end
    seg = segs[i]
    entry = _daf_entry(ks, seg.file)
    daf = entry.content::DAF
    # Of two segments with the same data, the later one has the priority.
    k = findlast(s -> s.baddr == seg.baddr && s.eaddr == seg.eaddr && s.target == seg.target,
        daf.segments)
    descr = dafps(Float64[seg.t0, seg.t1],
        Int[seg.target, seg.center, seg.frame, seg.dtype, seg.baddr, seg.eaddr])
    return entry.path, descr, _daf_array_name(daf, k)
end
spksfs(ctx::SpiceContext, body::Integer, et::Real) = spksfs(kernelset(ctx), body, et)
spksfs(body::Integer, et::Real) = spksfs(default_context(), body, et)

# ---------------------------------------------------------------------------------------
# Evaluation of one segment
# ---------------------------------------------------------------------------------------

"""
    _spk_descr_segment(data, descr, path, shift) -> Segment

The segment of the descriptor `descr` over `data`, where DAF address `a` is
`data[a - shift]`; the directory values of types 2, 3 and 20 are read and checked as
`open_daf` does (`SpiceError(:baddaf)`).
"""
function _spk_descr_segment(data::Vector{Float64}, descr, path::AbstractString, shift::Int)
    dc, ic = dafus(descr, 2, 6)
    b, e = ic[5] - shift, ic[6] - shift
    @assert 1 <= b <= e <= length(data)
    dtype = Int32(ic[4])
    init, intlen, rsize, nrec = (dtype == 2 || dtype == 3 || dtype == 20) ?
        _daf_directory(data, path, :SPK, dtype, b, e) : (0.0, 0.0, Int32(0), Int32(0))
    return Segment(Int32(ic[1]), Int32(ic[2]), Int32(ic[3]), dtype, dc[1], dc[2], Int32(0),
        b, e, init, intlen, rsize, nrec)
end

function _spkpvn(ks::KernelSet, data::Vector{Float64}, seg::Segment, et::Real)
    e = check_epoch(et)
    # The readers of CSPICE for the types 8, 9, 12, 13, 14, 18 and 19 make this check; for
    # the others CSPICE extrapolates, reads outside the segment or stops the process.
    seg.t0 <= primal(e) <= seg.t1 || throw(SpiceError(:timeoutofbounds, "spkpvn: request " *
        "time $(primal(e)) is outside of descriptor bounds $(seg.t0) : $(seg.t1)"))
    p, v = seg.dtype == Int32(10) ? _t10_jet(_lsk(ks), data, seg, e, Val(1)) :
           spk_eval(data, seg, e, Val(1))
    return Int(seg.frame), vcat(p, v), Int(seg.center)
end

"The loaded SPK file with the path `file`; throws `SpiceError(:nosuchhandle)` if there is none."
function _spk_loaded(ks::KernelSet, file::AbstractString)
    for entry in ks.loaded
        c = entry.content
        entry.path == file && c isa DAF && c.kind === :SPK && return c
    end
    throw(SpiceError(:nosuchhandle, "spkpvn: no SPK file with the path $file is loaded " *
        "in this kernel set"))
end

"""
    spkpvn([ctx,] handle::DAFHandle, descr, et) -> (ref, state, center)
    spkpvn([ctx,] file::AbstractString, descr, et) -> (ref, state, center)

Return the state of the target of one SPK segment relative to the center of the segment
at `et`, as `spkpvn_c`.

The segment has the descriptor `descr` (5 doubles, from [`dafgs`](@ref) or
[`spksfs`](@ref)). It is in the file of `handle`, or in the file with the path `file` that
is loaded in the context (the `file` of `spksfs`; it must be the same string, because
the path is not normalized). `ref` is the ID of the frame of the segment and `center`
the NAIF ID of its center, both `Int`. `state` is an `SVector{6}`: position in
km, then velocity in km/s, in the frame `ref`. `et` is in seconds past J2000 TDB.

The context gives the leap seconds that a segment of type 10 needs; the form with a
handle reads nothing else from it.

Throws `SpiceError`:
- `:timeoutofbounds` if `et` is outside the coverage of `descr`, for every type. CSPICE
  signals this for the types 8, 9, 12, 13, 14, 18 and 19 only; for the other types it
  extrapolates, reads data outside the segment, or stops the process (type 1);
- `:spktypenotsupp` if the data type is not one of CSPICE version N0067;
- `:nosuchhandle` if no SPK file with the path `file` is loaded in the context, and
  `:dafnosuchhandle` if `handle` is closed;
- `:dafnegaddr`, `:dafbeggtend` or `:filereadfailed` if the addresses of `descr` are not
  a range in the file, and `:baddaf` if the data of the segment are not consistent;
- `:badarg` if `et` is not finite;
- `:noleapseconds` for a segment of type 10 and a context without leap seconds.

See also [`spksfs`](@ref) for an example, [`spkuds`](@ref).
"""
function spkpvn(ks::KernelSet, h::DAFHandle, descr::AbstractVector{<:Real}, et::Real)
    _daf_check_open(h, "spkpvn")
    _, ic = dafus(descr, 2, 6)
    data = dafgda(h, ic[5], ic[6])
    return _spkpvn(ks, data, _spk_descr_segment(data, descr, h.path, ic[5] - 1), et)
end
function spkpvn(ks::KernelSet, file::AbstractString, descr::AbstractVector{<:Real}, et::Real)
    daf = _spk_loaded(ks, file)
    _, ic = dafus(descr, 2, 6)
    b, e = ic[5], ic[6]
    b <= 0 && throw(SpiceError(:dafnegaddr, "spkpvn: negative value for the first address: $b"))
    b > e && throw(SpiceError(:dafbeggtend,
        "spkpvn: beginning address ($b) greater than ending address ($e)"))
    e <= length(daf.data) || throw(SpiceError(:filereadfailed,
        "spkpvn: the addresses $b:$e are not inside the file $file"))
    return _spkpvn(ks, daf.data, _spk_descr_segment(daf.data, descr, daf.path, 0), et)
end
spkpvn(ctx::SpiceContext, src::Union{DAFHandle,AbstractString}, descr::AbstractVector{<:Real},
    et::Real) = spkpvn(kernelset(ctx), src, descr, et)
spkpvn(src::Union{DAFHandle,AbstractString}, descr::AbstractVector{<:Real}, et::Real) =
    spkpvn(default_context(), src, descr, et)

# ---------------------------------------------------------------------------------------
# Subset of a segment
# ---------------------------------------------------------------------------------------

"""
    _SegReader(handle, baddr, eaddr)

The data of the segment at `baddr:eaddr` of an open DAF, read once. `_rd` reads by DAF
address: an address inside the segment from the copy, another one from the file, as the
`dafgda` calls of the C code do.
"""
struct _SegReader
    h::DAFHandle
    baddr::Int
    data::Vector{Float64}
end
_SegReader(h::DAFHandle, baddr::Int, eaddr::Int) =
    _SegReader(h, baddr, baddr <= eaddr ? dafgda(h, baddr, eaddr) : Float64[])

function _rd(r::_SegReader, a::Int, b::Int)
    i, j = a - r.baddr + 1, b - r.baddr + 1
    (1 <= i <= j <= length(r.data)) && return r.data[i:j]
    return dafgda(r.h, a, b)
end
function _rd(r::_SegReader, a::Int)
    i = a - r.baddr + 1
    return 1 <= i <= length(r.data) ? r.data[i] : dafgda(r.h, a, a)[1]
end

# `max` and `min` of f2c (`a >= b ? a : b`), which differ from Julia for -0.0 and NaN.
_fmax(a::Float64, b::Float64) = a >= b ? a : b
_fmin(a::Float64, b::Float64) = a <= b ? a : b

"Types 1 and 21 (`spks01`, `spks21`): the difference lines whose final epochs cover `t0:t1`."
function _spks01(r::_SegReader, baddr::Int, eaddr::Int, t0::Float64, t1::Float64, t21::Bool)
    out = Float64[]
    if t21
        maxdim = _nint(_rd(r, eaddr - 1))
        nrec = _nint(_rd(r, eaddr))
        offe = eaddr - nrec ÷ 100 - nrec - 2
        dlsize = 4maxdim + 11
    else
        nrec = _cint(_rd(r, eaddr))
        offe = eaddr - nrec ÷ 100 - nrec - 1
        dlsize = 71
    end
    first = last = 0
    for i in 1:nrec      # invariant: `first` and `last` are those of the epochs before `i`
        ep = _rd(r, offe + i)
        first == 0 && ep >= t0 && (first = i)
        first != 0 && last == 0 && ep >= t1 && (last = i)
    end
    offset = baddr - 1 + (first - 1) * dlsize
    for _ in first:last
        append!(out, _rd(r, offset + 1, offset + dlsize))
        offset += dlsize
    end
    for i in first:last
        push!(out, _rd(r, offe + i))
    end
    for i in first+99:100:last
        push!(out, _rd(r, offe + i))
    end
    t21 && push!(out, Float64(maxdim))
    push!(out, Float64(last - first + 1))
    return out
end

"Types 2 and 3 (`spks02`, `spks03`): the records of equal length that cover `t0:t1`."
function _spks02(r::_SegReader, baddr::Int, eaddr::Int, t0::Float64, t1::Float64)
    init, intlen, rs, nr = _rd(r, eaddr - 3, eaddr)
    recsiz = _cint(rs)
    nrec = _cint(nr)
    first = min(_cint((t0 - init) / intlen) + 1, nrec)
    last = min(_cint((t1 - init) / intlen) + 1, nrec)
    nrec = last - first + 1
    remain = nrec * recsiz
    addr = baddr + (first - 1) * recsiz
    out = remain > 0 ? _rd(r, addr, addr + remain - 1) : Float64[]
    init += (first - 1) * intlen
    push!(out, init, intlen, Float64(recsiz), Float64(nrec))
    return out
end

"Type 5 (`spks05`): the states of the epochs that bracket `t0:t1`."
function _spks05(r::_SegReader, baddr::Int, eaddr::Int, t0::Float64, t1::Float64)
    gm, nr = _rd(r, eaddr - 1, eaddr)
    nrec = _cint(nr)
    offe = eaddr - nrec ÷ 100 - nrec - 2
    r2 = 1
    ep = _rd(r, offe + r2)
    while r2 < nrec && ep < t1           # invariant: the epochs before r2 are before t1
        r2 += 1
        ep = _rd(r, offe + r2)
    end
    r1 = nrec
    ep = _rd(r, offe + r1)
    while r1 > 1 && ep > t0              # invariant: the epochs after r1 are after t0
        r1 -= 1
        ep = _rd(r, offe + r1)
    end
    out = Float64[]
    for i in r1:r2
        offset = baddr - 1 + (i - 1) * 6
        append!(out, _rd(r, offset + 1, offset + 6))
    end
    for i in r1:r2
        push!(out, _rd(r, offe + i))
    end
    for i in r1+99:100:r2
        push!(out, _rd(r, offe + i))
    end
    push!(out, gm, Float64(r2 - r1 + 1))
    return out
end

"Types 8 and 12 (`spks08`): the states of equal spacing that the interpolation over `t0:t1` needs."
function _spks08(r::_SegReader, baddr::Int, eaddr::Int, t0::Float64, t1::Float64)
    start, step, dg, nr = _rd(r, eaddr - 3, eaddr)
    degree = _nint(dg)
    nrec = _nint(nr)
    (t1 < t0 || t1 < start || t0 > start + (nrec - 1) * step) && return Float64[]
    ratio = _fmax(0.0, (t0 - start) / step)
    r1 = min(_cint(ratio), nrec - 1) + 1
    r1 = max(1, r1 - degree ÷ 2)
    r1 = min(r1, nrec - degree)
    ratio = (t1 - start) / step
    r2 = ratio == trunc(ratio) ? min(_cint(ratio), nrec - 1) + 1 :
         min(_cint(ratio) + 1, nrec - 1) + 1
    r2 = min(nrec, r2 + degree ÷ 2)
    r2 = max(r2, degree + 1)
    out = Float64[]
    for i in r1:r2
        offset = baddr - 1 + (i - 1) * 6
        append!(out, _rd(r, offset + 1, offset + 6))
    end
    push!(out, start + (r1 - 1) * step, step, Float64(degree), Float64(r2 - r1 + 1))
    return out
end

"""
    _spks_window(r, offe, nrec, t0, t1, half) -> (r1, r2)

The first and the last record of the epochs at `offe + 1 : offe + nrec` that a subset to
`t0:t1` keeps, with `half` more records at each end, as `spks09` and `spks18` find them.
"""
function _spks_window(r::_SegReader, offe::Int, nrec::Int, t0::Float64, t1::Float64, half::Int)
    r2 = 1
    ep = _rd(r, offe + r2)
    while r2 < nrec && ep < t1           # invariant: the epochs before r2 are before t1
        r2 += 1
        ep = _rd(r, offe + r2)
    end
    r2 = min(nrec, r2 + half)
    r1 = nrec
    ep = _rd(r, offe + r1)
    while r1 > 1 && ep > t0              # invariant: the epochs after r1 are after t0
        r1 -= 1
        ep = _rd(r, offe + r1)
    end
    return max(1, r1 - half), r2
end

"Types 9 and 13 (`spks09`): the states of unequal spacing that the interpolation over `t0:t1` needs."
function _spks09(r::_SegReader, baddr::Int, eaddr::Int, t0::Float64, t1::Float64)
    dg, nr = _rd(r, eaddr - 1, eaddr)
    degree = _nint(dg)
    nrec = _nint(nr)
    offe = eaddr - (nrec - 1) ÷ 100 - nrec - 2
    r1, r2 = _spks_window(r, offe, nrec, t0, t1, degree ÷ 2)
    r2 = max(r2, degree + 1)
    r1 = min(r1, nrec - degree)
    out = Float64[]
    for i in r1:r2
        offset = baddr - 1 + (i - 1) * 6
        append!(out, _rd(r, offset + 1, offset + 6))
    end
    for i in r1:r2
        push!(out, _rd(r, offe + i))
    end
    for i in r1+99:100:r2-1
        push!(out, _rd(r, offe + i))
    end
    push!(out, Float64(degree), Float64(r2 - r1 + 1))
    return out
end

"Type 18 (`spks18`): the packets of the interpolation window over `t0:t1`."
function _spks18(r::_SegReader, baddr::Int, eaddr::Int, t0::Float64, t1::Float64)
    st, ws, nr = _rd(r, eaddr - 2, eaddr)
    subtyp = _nint(st)
    winsiz = _nint(ws)
    nrec = _nint(nr)
    packsz = subtyp == 0 ? 12 : subtyp == 1 ? 6 : throw(SpiceError(:invalidvalue,
        "spksub: unexpected SPK type 18 subtype $subtyp found in type 18 record"))
    offe = eaddr - (nrec - 1) ÷ 100 - nrec - 3
    r1, r2 = _spks_window(r, offe, nrec, t0, t1, (winsiz - 1) ÷ 2)
    out = Float64[]
    for i in r1:r2
        offset = baddr - 1 + (i - 1) * packsz
        append!(out, _rd(r, offset + 1, offset + packsz))
    end
    for i in r1:r2
        push!(out, _rd(r, offe + i))
    end
    for i in r1+99:100:r2-1
        push!(out, _rd(r, offe + i))
    end
    push!(out, Float64(subtyp), Float64(winsiz), Float64(r2 - r1 + 1))
    return out
end

"Type 20 (`spks20`): the records of equal length that cover `t0:t1`, with the new initial epoch."
function _spks20(r::_SegReader, baddr::Int, eaddr::Int, t0::Float64, t1::Float64)
    dscale, tscale, initjd, initfr, intlen, rs, nr = _rd(r, eaddr - 6, eaddr)
    recsiz = _cint(rs)
    nrec = _cint(nr)
    btime = (initjd - j2000() + initfr) * spd()
    intrvl = intlen * spd()
    first = max(1, min(_cint((t0 - btime) / intrvl) + 1, nrec))
    last = max(1, min(_cint((t1 - btime) / intrvl) + 1, nrec))
    nrec = last - first + 1
    remain = nrec * recsiz
    addr = baddr + (first - 1) * recsiz
    out = remain > 0 ? _rd(r, addr, addr + remain - 1) : Float64[]
    subbeg = j2000() + (btime + (first - 1) * intrvl) / spd()
    subijd = trunc(subbeg)
    push!(out, dscale, tscale, subijd, subbeg - subijd, intlen, Float64(recsiz), Float64(nrec))
    return out
end

@noinline _spks19_error(msg::String) = throw(SpiceError(:spkstructureerror, "spksub: " * msg))

"""
    _spks19_last(r, minbep, npkt, bufbas, nread, t1) -> (l, bufbas, nread, buf)

The search of `spks19` for the last epoch `<= t1` of a mini-segment, by blocks of 100
epochs from the block at `bufbas + 1 : bufbas + nread`. CSPICE stops the process with
a subscript error when the block in which the search ends starts after `t1` (`l == 0`);
here that case throws `SpiceError(:spkstructureerror)`.
"""
function _spks19_last(r::_SegReader, minbep::Int, npkt::Int, bufbas::Int, nread::Int, t1::Float64)
    buf = _rd(r, bufbas + 1, bufbas + nread)
    remain = npkt - (bufbas - minbep) - nread
    while remain > 0 && buf[nread] <= t1  # invariant: the epochs up to `bufbas + nread` are <= t1
        bufbas += nread
        nread = min(remain, 100)
        buf = _rd(r, bufbas + 1, bufbas + nread)
        remain -= nread
    end
    l = _lstled(buf, 0, nread, t1)
    return l, bufbas, nread, buf
end

"Type 19 (`spks19`): the mini-segments of the intervals that cover `t0:t1`, the first and the last one cut."
function _spks19(r::_SegReader, baddr::Int, eaddr::Int, t0::Float64, t1::Float64)
    out = Float64[]
    final = false
    sl, ni = _rd(r, eaddr - 1, eaddr)
    isel = _nint(sl)
    nintvl = _nint(ni)
    ndir = nintvl ÷ 100
    ptrbas = eaddr - (nintvl + 3)
    ivlbas = ptrbas - (ndir + nintvl + 1)
    # The interval that holds t0.
    nread = min(100, nintvl)
    bufbas = ivlbas
    buf = _rd(r, bufbas + 1, bufbas + nread)
    remain = nintvl - nread
    while remain > 0 && buf[nread] < t0  # invariant: the bounds up to `bufbas + nread` are < t0
        bufbas += nread
        nread = min(100, remain)
        buf = _rd(r, bufbas + 1, bufbas + nread)
        remain -= nread
    end
    begidx = max(1, bufbas - ivlbas + _lstltd(buf, 0, nread, t0))
    p1, p2 = _rd(r, ptrbas + begidx, ptrbas + begidx + 1)
    minib = baddr - 1 + _nint(p1)
    minie = baddr - 1 + _nint(p2) - 1
    subtyp, wndsiz, npkt = _nint.(_rd(r, minie - 2, minie))
    0 <= subtyp < 3 || throw(SpiceError(:notsupported, "spksub: unexpected SPK type 19 " *
        "subtype $subtyp found in type 19 segment within mini-segment $begidx"))
    pktsiz = _SPK19_PKTSZS[subtyp+1]
    minbep = minib - 1 + npkt * pktsiz
    # The first packet of the first mini-segment.
    nread = min(100, npkt)
    bufbas = minbep
    buf = _rd(r, bufbas + 1, bufbas + nread)
    remain = npkt - nread
    while remain > 0 && buf[nread] < t0  # invariant: the epochs up to `bufbas + nread` are < t0
        bufbas += nread
        nread = min(100, remain)
        buf = _rd(r, bufbas + 1, bufbas + nread)
        remain -= nread
    end
    bepidx = bufbas - minbep + _lstled(buf, 0, nread, t0)
    bepidx = min(max(1, bepidx - (wndsiz ÷ 2 - 1)), npkt - 1)
    # The last packet of the first mini-segment.
    ivlend = _rd(r, ivlbas + begidx + 1)
    if ivlend <= t1
        eepidx = npkt
    else
        nread = min(100, npkt - bepidx + 1)
        nread < 2 && _spks19_error("NREAD is $nread; end time of interval $begidx is $ivlend")
        l, bufbas, nread, buf = _spks19_last(r, minbep, npkt, minbep + bepidx - 1, nread, t1)
        eepidx = bufbas - minbep + l + 1
        (eepidx < 2 || eepidx > npkt) && _spks19_error("last epoch is $(buf[nread]); " *
            "end time of interval $begidx is $ivlend")
        l == 0 && _spks19_error("the epochs of interval $begidx from $(buf[1]) are after " *
            "the end $t1 (CSPICE stops with a subscript error)")
        eepidx = min(npkt, eepidx + (buf[l] == t1 ? wndsiz ÷ 2 - 2 : wndsiz ÷ 2 - 1))
    end
    ivlbeg = _rd(r, ivlbas + begidx)
    for i in bepidx:eepidx
        b = minib - 1 + (i - 1) * pktsiz
        append!(out, _rd(r, b + 1, b + pktsiz))
    end
    iv1beg = iv1end = 0.0
    for i in bepidx:eepidx
        ep = _rd(r, minib - 1 + npkt * pktsiz + i)
        push!(out, ep)
        i == bepidx && (iv1beg = _fmax(ivlbeg, ep))
        i == eepidx && (iv1end = _fmin(ivlend, ep))
    end
    minnpk = eepidx - bepidx + 1
    minndr = (minnpk - 1) ÷ 100
    for i in 1:minndr
        push!(out, _rd(r, minib - 1 + npkt * pktsiz + (bepidx - 1) + 100i))
    end
    push!(out, Float64(subtyp), Float64(wndsiz), Float64(minnpk))
    min1sz = minnpk * (pktsiz + 1) + minndr + 3
    # The mini-segments that follow, whole, and the last one, cut.
    ivfbeg = ivfend = 0.0
    minfsz = 0
    if ivlend > t1 || begidx == nintvl
        endidx = begidx
    else
        curivl = begidx + 1
        ivfbeg = ivlend
        ivlend = _rd(r, ivlbas + curivl + 1)
        while ivlend <= t1 && curivl <= nintvl
            # invariant: the mini-segments begidx+1 .. curivl-1 are copied whole
            q1, q2 = _rd(r, ptrbas + curivl, ptrbas + curivl + 1)
            minib = baddr - 1 + _nint(q1)
            minie = baddr - 1 + _nint(q2) - 1
            minie >= minib && append!(out, _rd(r, minib, minie))
            ivfbeg = ivlend
            curivl += 1
            curivl <= nintvl && (ivlend = _rd(r, ivlbas + curivl + 1))
        end
        if curivl > nintvl
            endidx = nintvl
        else
            final = true
            endidx = curivl
            q1, q2 = _rd(r, ptrbas + endidx, ptrbas + endidx + 1)
            minib = baddr - 1 + _nint(q1)
            minie = baddr - 1 + _nint(q2) - 1
            subtyp, wndsiz, npkt = _nint.(_rd(r, minie - 2, minie))
            0 <= subtyp < 3 || throw(SpiceError(:notsupported, "spksub: unexpected SPK " *
                "type 19 subtype $subtyp found in type 19 segment within mini-segment $curivl"))
            pktsiz = _SPK19_PKTSZS[subtyp+1]
            minbep = minib - 1 + npkt * pktsiz
            l, bufbas, nread, buf = _spks19_last(r, minbep, npkt, minbep, min(100, npkt), t1)
            eepidx = bufbas - minbep + l + 1
            eepidx < 2 && _spks19_error("no epochs in final input interval exceed END; " *
                "interval index is $endidx, END is $t1")
            l == 0 && _spks19_error("the epochs of interval $endidx from $(buf[1]) are " *
                "after the end $t1 (CSPICE stops with a subscript error)")
            eepidx = max(min(npkt, eepidx + (buf[l] == t1 ? wndsiz ÷ 2 - 2 : wndsiz ÷ 2 - 1)), 2)
            for i in 1:eepidx
                b = minib - 1 + (i - 1) * pktsiz
                append!(out, _rd(r, b + 1, b + pktsiz))
            end
            for i in 1:eepidx
                ep = _rd(r, minib - 1 + npkt * pktsiz + i)
                push!(out, ep)
                i == eepidx && (ivfend = ep)
            end
            minndr = (eepidx - 1) ÷ 100
            for i in 1:minndr
                push!(out, _rd(r, minib - 1 + npkt * pktsiz + 100i))
            end
            push!(out, Float64(subtyp), Float64(wndsiz), Float64(eepidx))
            minfsz = eepidx * (pktsiz + 1) + minndr + 3
        end
    end
    # Interval bounds, their directory, the mini-segment pointers and the control words.
    noivl = endidx - begidx + 1
    push!(out, iv1beg)
    if noivl == 1
        push!(out, iv1end)
    else
        for i in 2:(final ? noivl - 1 : noivl + 1)
            push!(out, _rd(r, ivlbas + (begidx - 1) + i))
        end
        final && push!(out, ivfbeg, ivfend)
    end
    for i in 1:noivl÷100
        push!(out, _rd(r, ivlbas + (begidx - 1) + 100i))
    end
    push!(out, 1.0)
    if noivl == 1
        push!(out, Float64(min1sz + 1))
    else
        start = min1sz + 1
        shift = 0
        for i in 2:(final ? noivl : noivl + 1)
            v = _nint(_rd(r, ptrbas + (begidx - 1) + i))
            i == 2 && (shift = min1sz + 1 - v)
            start = v + shift
            push!(out, Float64(start))
        end
        final && push!(out, Float64(start + minfsz))
    end
    push!(out, Float64(isel), Float64(noivl))
    return out
end

"A `Segment` of type `dtype` over the local copy of `r`, for the layout readers of types 10 and 14."
_spks_local(r::_SegReader, dtype::Integer) = Segment(Int32(0), Int32(0), Int32(0),
    Int32(dtype), 0.0, 0.0, Int32(0), 1, length(r.data), 0.0, 0.0, Int32(0), Int32(0))

"Type 10 (`spks10`): the element sets closest to `t0` and `t1` and those between them."
function _spks10(newh::DAFHandle, r::_SegReader, ndscr::Vector{Float64}, ident::AbstractString,
        t0::Float64, t1::Float64)
    d = r.data
    lay = _t10_layout(d, _spks_local(r, 10))
    lay.pksz == 14 || throw(SpiceError(:notsupported, "spksub: a type 10 segment with " *
        "packets of $(lay.pksz) doubles cannot be subset; spkw10 writes 14"))
    con = d[lay.con+1:lay.con+8]
    dafbna(newh, ndscr, ident)
    begidx = _t10_closest(d, lay, t0)
    d[lay.ref+begidx] > t0 && (begidx = max(1, begidx - 1))
    endidx = _t10_closest(d, lay, t1)
    d[lay.ref+endidx] < t1 && (endidx = min(lay.n, endidx + 1))
    refs = d[lay.ref+begidx:lay.ref+endidx]
    packets = Float64[]
    for i in begidx:endidx
        a = lay.pkt + (i - 1) * lay.stride
        append!(packets, @view d[a+1:a+14])
    end
    dafada(newh, _t10_data(con, refs, packets))
    dafena(newh)
    return nothing
end

"Type 14 (`spks14`): the records from the one of `t0` to the one of `t1`, through `spk14b`."
function _spks14(newh::DAFHandle, r::_SegReader, ndscr::Vector{Float64}, ident::AbstractString)
    dc, ic = dafus(ndscr, 2, 6)
    t0, t1 = dc
    d = r.data
    lay = _spk14_layout(d, _spks_local(r, 14))
    chbdeg = _cint(d[lay.conbas+1]) - 1
    recsiz = (chbdeg + 1) * 6 + 2
    recsiz > 128 && throw(SpiceError(:spkrectoolarge, "spksub: storage for $recsiz double " *
        "precision numbers is needed for an SPK data record and only 128 locations were available"))
    begidx = max(1, _lstled(d, lay.refbas, lay.nref, t0))
    endidx = max(1, _lstled(d, lay.refbas, lay.nref, t1))
    # The frame goes through its name in CSPICE (`irfnam`), which only an inertial frame
    # of the built-in table has: another frame is the blank name, which `spkpds` rejects.
    frame = 1 <= ic[3] <= 21 ? ic[3] : 0
    w = spk14b(newh, ident, ic[1], ic[2], frame, t0, t1, chbdeg)
    n = endidx - begidx + 1
    if n > 0
        records = Float64[]
        for i in begidx:endidx
            a = lay.pktbas + (i - 1) * lay.stride + lay.pktoff
            append!(records, @view d[a+1:a+lay.pktsz])
        end
        spk14a(w, n, records, d[lay.refbas+begidx:lay.refbas+endidx])
    end
    spk14e(w)
    return nothing
end

"""
    spksub(newh, handle, descr, ident, start, stop) -> Nothing

Write the part of an SPK segment that covers the interval `[start, stop]` as a new
segment of another SPK file, as `spksub_c`.

- `newh`: the [`DAFHandle`](@ref) of the new file, open for writing ([`spkopn`](@ref),
  [`spkopa`](@ref)).
- `handle`, `descr`: the file that holds the source segment ([`dafopr`](@ref)) and the
  descriptor of the segment ([`dafgs`](@ref), [`spksfs`](@ref)).
- `ident`: the identifier of the new segment.
- `start`, `stop`: the coverage of the new segment, seconds past J2000 TDB.

The new segment holds the records of the source that the evaluation in `[start, stop]`
needs: the data that `spksub_c` writes from the same source, double for double. It is on
disk when the call returns. `newh` comes first, as in `spksub!` of SPICE.jl; CSPICE has it
last.

In `[start, stop]` the new segment gives the states of the source segment, within
rounding for the types 8 and 12. For type 20 the states differ from those of the source
by up to about 3e-8 relative, in CSPICE as here, because the new segment stores its first
epoch as a Julian date.

For the types other than 14, an `ident` of more than 40 characters is cut to 40, as in
CSPICE. A segment of type 14 is written with [`spk14b`](@ref), so `ident` must be a
valid identifier, `start` must be before `stop`, and the frame of the source must be an
inertial frame of the built-in table (IDs 1 to 21).

Throws `SpiceError`:
- `:emptystring` if `ident` is empty;
- `:spknotasubset` unless `first <= start <= stop <= last`, where `[first, last]` is
  the coverage of `descr`;
- `:spktypenotsupp` if the data type is not one of CSPICE N0067;
- for type 14: `:spkrectoolarge` if the degree of the polynomials is above 20,
  `:invalidrefframe` if the frame is not inertial, `:baddescrtimes` if `start == stop`,
  and the error kinds of [`spk14b`](@ref) for the identifier;
- `:invalidvalue` (type 18), `:notsupported` (type 19) for an unknown subtype, and
  `:spkstructureerror` for a type 19 segment whose structure is damaged, or where the
  block of 100 epochs in which the search for `stop` ends begins after `stop`, a case
  in which CSPICE stops the process;
- the error kinds of [`dafbna`](@ref) and [`dafgda`](@ref) if a handle is closed, `newh`
  is not open for writing, or the addresses of `descr` are not in the file.

A call that throws writes nothing and leaves no array begun on `newh`.

# Examples
```jldoctest
julia> dir = mktempdir();

julia> h = spkopn(joinpath(dir, "source.bsp"), "source");

julia> cdata = [[25.0, 0.25, 50.0, 0.0, 75.0, 0.0] for _ in 1:4];   # four records of 25 s

julia> spkw02(h, -9, 399, "J2000", 0.0, 100.0, "four records", 25.0, cdata, 0.0)

julia> spkcls(h)

julia> src = dafopr(joinpath(dir, "source.bsp")); dafbfs(src); daffna(src);

julia> newh = spkopn(joinpath(dir, "part.bsp"), "part");

julia> spksub(newh, src, dafgs(src), "records 2 and 3", 30.0, 60.0)

julia> spkcls(newh); dafcls(src)

julia> h = dafopr(joinpath(dir, "part.bsp")); dafbfs(h); daffna(h);

julia> spkuds(dafgs(h))
(-9, 399, 1, 2, 30.0, 60.0, 385, 404)

julia> dafcls(h)
```
"""
function spksub(newh::DAFHandle, handle::DAFHandle, descr::AbstractVector{<:Real},
        ident::AbstractString, start::Real, stop::Real)
    _seg_chkfstr("spksub", "ident", ident)
    dc, ic = dafus(descr, 2, 6)
    alpha, omega = dc
    t0, t1 = Float64(start), Float64(stop)
    (alpha <= t0 && t0 <= t1 && t1 <= omega) || throw(SpiceError(:spknotasubset,
        "spksub: specified interval [$t0, $t1] is not a subset of segment interval " *
        "[$alpha, $omega]"))
    ndscr = dafps(Float64[t0, t1], ic)
    dtype, baddr, eaddr = ic[4], ic[5], ic[6]
    if dtype == 10 || dtype == 14
        # These two read the source before they begin the new segment, as in CSPICE.
        r = _SegReader(handle, baddr, eaddr)
        dtype == 10 ? _spks10(newh, r, ndscr, ident, t0, t1) : _spks14(newh, r, ndscr, ident)
        return nothing
    end
    dtype in (1, 2, 3, 5, 8, 9, 12, 13, 15, 17, 18, 19, 20, 21) || throw(SpiceError(
        :spktypenotsupp, "spksub: SPK data type $dtype is not supported"))
    # The data are read before the array is begun, so that a read error leaves no array
    # begun on `newh`; the checks of the handle come first, as in CSPICE.
    _daf_check_write(newh, "dafbna")
    r = _SegReader(handle, baddr, eaddr)
    data = if dtype == 1 || dtype == 21
        _spks01(r, baddr, eaddr, t0, t1, dtype == 21)
    elseif dtype == 2 || dtype == 3
        _spks02(r, baddr, eaddr, t0, t1)
    elseif dtype == 5
        _spks05(r, baddr, eaddr, t0, t1)
    elseif dtype == 8 || dtype == 12
        _spks08(r, baddr, eaddr, t0, t1)
    elseif dtype == 9 || dtype == 13
        _spks09(r, baddr, eaddr, t0, t1)
    elseif dtype == 15 || dtype == 17
        # `spks15`, `spks17`: the whole segment of 16 or 12 elements.
        _rd(r, baddr, baddr + (dtype == 15 ? 15 : 11))
    elseif dtype == 18
        _spks18(r, baddr, eaddr, t0, t1)
    elseif dtype == 19
        _spks19(r, baddr, eaddr, t0, t1)
    else
        _spks20(r, baddr, eaddr, t0, t1)
    end
    dafbna(newh, ndscr, ident)
    dafada(newh, data)
    dafena(newh)
    return nothing
end
