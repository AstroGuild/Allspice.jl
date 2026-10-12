# Spacecraft clocks (SCLK), type 1, as `sc01.f` and its callers in CSPICE.
#
# An encoded SCLK value ("ticks") counts the ticks of the least significant clock field
# from the start of the first partition. A type 1 clock relates ticks to a parallel time
# system (TDB or TDT) by a piecewise linear function, given as coefficient records
# (ticks, parallel time, rate). The data of every clock are read from the kernel pool
# once, when the `KernelSet` snapshot is built; the conversions are pure functions of
# that snapshot.

"Delimiters that separate the fields of a clock string on input (`sctk01`)."
const _SCLK_DELIMS = ".:-, "
"Output delimiters by code 1 to 5 (`SCLK01_OUTPUT_DELIM`)."
const _SCLK_OUTPUT_DELIM = ('.', ':', '-', ',', ' ')
"Largest number of clock fields (`MXNFLD`)."
const _SCLK_MXNFLD = 10

"`_dnint(x)`: nearest whole number, halves away from zero, as `d_nint` (Fortran `ANINT`)."
@inline _dnint(x::Float64) = round(x, RoundNearestTiesAway)

"""
    _sclk_build(pool, sc) -> Union{SclkType1,Nothing}

Data of the type 1 clock `sc` from the pool variables with the suffix `_<-sc>`, or
`nothing` if a required variable is absent or malformed. Required: `SCLK01_COEFFICIENTS`
(a positive multiple of three values), `SCLK_PARTITION_START` and `SCLK_PARTITION_END`
(equal counts, at least one), `SCLK01_MODULI`, `SCLK01_OFFSETS`, `SCLK01_N_FIELDS` (1 to
10, not more than the moduli and the offsets) and `SCLK01_OUTPUT_DELIM`.
`SCLK01_TIME_SYSTEM` is optional and defaults to 1 (TDB), as in `zzscad01`.
"""
function _sclk_build(pool::AbstractDict, sc::Integer)
    suffix = "_$(-Int(sc))"
    num(name) = _pool_numbers(pool, name * suffix)
    coef = num("SCLK01_COEFFICIENTS")
    pstart = num("SCLK_PARTITION_START")
    pstop = num("SCLK_PARTITION_END")
    moduli = num("SCLK01_MODULI")
    offsets = num("SCLK01_OFFSETS")
    nf = num("SCLK01_N_FIELDS")
    del = num("SCLK01_OUTPUT_DELIM")
    any(isnothing, (coef, pstart, pstop, moduli, offsets, nf, del)) && return nothing
    (length(coef) >= 3 && length(coef) % 3 == 0) || return nothing
    (length(pstart) == length(pstop) && !isempty(pstart)) || return nothing
    (length(nf) == 1 && length(del) == 1) || return nothing
    nfield = round(Int, nf[1])
    (1 <= nfield <= _SCLK_MXNFLD && length(moduli) >= nfield &&
     length(offsets) >= nfield) || return nothing
    ts = num("SCLK01_TIME_SYSTEM")
    timsys = (ts === nothing || isempty(ts)) ? 1 : round(Int, ts[1])
    tikmsc = 1.0
    for i in nfield:-1:2              # invariant: tikmsc = product of moduli[i+1:nfield]
        tikmsc *= moduli[i]
    end
    mxtick = 0.0
    for i in eachindex(pstart)        # invariant: mxtick = rounded length of partitions 1:i-1
        mxtick = _dnint(pstop[i] - pstart[i] + mxtick)
    end
    return SclkType1(Int32(sc), coef, pstart, pstop, moduli, offsets, Int32(nfield),
        Int32(round(Int, del[1])), Int32(timsys), tikmsc, mxtick)
end

"""
    build_sclk(pool) -> Dict{Int32,SclkType1}

Type 1 clocks of the pool, by clock ID: one for each variable `SCLK_DATA_TYPE_<n>` with
the value 1 (the clock ID is `-n`) whose `SCLK01_*` data are complete. A clock with
incomplete data is left out, so that loading a kernel never fails; a conversion with
that clock then throws `SpiceError`.
"""
function build_sclk(pool::AbstractDict)
    out = Dict{Int32,SclkType1}()
    prefix = "SCLK_DATA_TYPE_"
    for (name, value) in pool
        startswith(name, prefix) || continue
        n = tryparse(Int32, SubString(name, ncodeunits(prefix) + 1))
        n === nothing && continue
        (value isa Vector{Float64} && length(value) == 1 && value[1] == 1) || continue
        clock = _sclk_build(pool, -n)
        clock === nothing || (out[-n] = clock)
    end
    return out
end

"""
    sctype(ks, sc) -> Int

Data type of the spacecraft clock `sc` (`SCLK_DATA_TYPE_<-sc>`), as `sctype_`. Throws
`SpiceError(:kernelvarnotfound)` if no SCLK kernel for `sc` is loaded.
"""
function sctype(ks::KernelSet, sc::Integer)
    name = "SCLK_DATA_TYPE_$(-Int(sc))"
    v = _pool_numbers(ks.pool, name)
    (v === nothing || isempty(v)) && throw(SpiceError(:kernelvarnotfound,
        "kernel variable $name was not found in the kernel pool; no SCLK kernel is " *
        "loaded for the clock $sc"))
    return round(Int, v[1])
end

"Throw the error that tells why the clock `sc` has no data in `ks.sclk`."
@noinline function _sclk_missing(ks::KernelSet, sc::Integer)
    t = sctype(ks, sc)
    t == 1 || throw(SpiceError(:unsupported, "clock type $t of the clock $sc is not supported"))
    throw(SpiceError(:kernelvarnotfound, "the SCLK01_* or SCLK_PARTITION_* kernel " *
        "variables of the spacecraft clock $sc are absent or malformed"))
end

"""
    _sclk(ks, sc) -> SclkType1

Data of the type 1 clock `sc`. Throws `SpiceError` with kind `:kernelvarnotfound` (no or
incomplete data) or `:unsupported` (a clock type other than 1).
"""
@inline function _sclk(ks::KernelSet, sc::Integer)
    clock = get(ks.sclk, Int32(sc), nothing)
    clock === nothing && _sclk_missing(ks, sc)
    return clock
end

@noinline _sclk_range(sc, what) = throw(SpiceError(:valueoutofrange,
    "$what of the spacecraft clock $sc is outside the range of the SCLK kernel"))
@noinline _sclk_badrate(sc) = throw(SpiceError(:valueoutofrange,
    "invalid SCLK rate (not positive) for the spacecraft clock $sc"))
@noinline _sclk_badsys(sc, sys) = throw(SpiceError(:valueoutofrange,
    "invalid time system code $sys was found for the spacecraft clock $sc"))
@noinline _sclk_nolsk() = throw(SpiceError(:nolsk,
    "a spacecraft clock with parallel time system TDT needs a leap-second kernel"))

"""
    _sclk_record(coef, col, x) -> Int

Index of the coefficient record to use for the value `x` of column `col` (1 ticks,
2 parallel time): the last record whose value is `<= x`, and the first record if `x` is
before all, found with the bisection of `scte01`. Requires `length(coef) >= 3`.
"""
@inline function _sclk_record(coef::Vector{Float64}, col::Int, x::Float64)
    n = length(coef) ÷ 3
    @inbounds x < coef[3*(n-1)+col] || return n
    lower, upper = 1, n
    while lower < upper - 1     # invariant: record `lower` is <= x or lower == 1; record `upper` is > x
        middle = (lower + upper) ÷ 2
        if @inbounds x < coef[3*(middle-1)+col]
            upper = middle
        else
            lower = middle
        end
    end
    return lower
end

"TDB seconds past J2000 of the parallel time `partim` of `clock` (`unitim_`), generic over the number type."
@inline function _sclk_to_tdb(ks::KernelSet, clock::SclkType1, partim::T) where {T<:Real}
    clock.timsys == 1 && return partim
    clock.timsys == 2 || _sclk_badsys(clock.id, clock.timsys)
    l = ks.lsk
    l === nothing && _sclk_nolsk()
    m0, m1 = l.m
    return partim + l.k * sin(m0 + m1 * partim + l.eb * sin(m0 + m1 * partim))
end

"Parallel time of `clock` for the TDB epoch `et` (`unitim_`, three fixed-point steps), generic over the number type."
@inline function _sclk_from_tdb(ks::KernelSet, clock::SclkType1, et::T) where {T<:Real}
    clock.timsys == 1 && return et
    clock.timsys == 2 || _sclk_badsys(clock.id, clock.timsys)
    l = ks.lsk
    l === nothing && _sclk_nolsk()
    m0, m1 = l.m
    tdt = et
    for _ in 1:3
        tdt = et - l.k * sin(m0 + m1 * tdt + l.eb * sin(m0 + m1 * tdt))
    end
    return tdt
end

"""
    _sct2e(ks, clock, sclkdp::T) -> T

ET of the encoded SCLK `sclkdp` (`scte01`). The function is piecewise linear, so a
`Taylor` or dual number passes through; the record is chosen with `primal`. Throws
`SpiceError(:valueoutofrange)` for ticks outside the clock range.
"""
@inline function _sct2e(ks::KernelSet, clock::SclkType1, sclkdp::T) where {T<:Real}
    coef = clock.coef
    p = primal(sclkdp)
    (p < coef[1] || p > clock.mxtick || isnan(p)) && _sclk_range(clock.id, "the encoded SCLK $p")
    k = 3 * (_sclk_record(coef, 1, p) - 1)
    @inbounds c1, c2, c3 = coef[k+1], coef[k+2], coef[k+3]
    c3 <= 0.0 && _sclk_badrate(clock.id)
    rate = c3 / clock.tikmsc
    return _sclk_to_tdb(ks, clock, c2 + rate * (sclkdp - c1))
end

"""
    _sce2c(ks, clock, et::T) -> T

Continuous encoded SCLK of the ET `et` (`scec01`), piecewise linear and generic over the
number type as `_sct2e`. Throws `SpiceError(:valueoutofrange)` for an epoch outside the
clock range.
"""
@inline function _sce2c(ks::KernelSet, clock::SclkType1, et::T) where {T<:Real}
    coef = clock.coef
    partim = _sclk_from_tdb(ks, clock, et)
    p = primal(partim)
    (p < coef[2] || isnan(p)) && _sclk_range(clock.id, "the epoch $(primal(et))")
    k = 3 * (_sclk_record(coef, 2, p) - 1)
    @inbounds c1, c2, c3 = coef[k+1], coef[k+2], coef[k+3]
    c3 <= 0.0 && _sclk_badrate(clock.id)
    rate = 1.0 / (c3 / clock.tikmsc)
    sclkdp = c1 + rate * (partim - c2)
    primal(sclkdp) > clock.mxtick && _sclk_range(clock.id, "the epoch $(primal(et))")
    return sclkdp
end

"""
    _sce2t(ks, clock, et) -> Float64

Discrete encoded SCLK of the ET `et` (`scet01`): the value of the first or last record
is extended to epochs before the first record, and the result is rounded to whole ticks.
Throws `SpiceError(:valueoutofrange)` if the rounded value is outside the clock range.
"""
@inline function _sce2t(ks::KernelSet, clock::SclkType1, et::Float64)
    coef = clock.coef
    partim = _sclk_from_tdb(ks, clock, et)
    lower = partim < coef[2] ? 1 : _sclk_record(coef, 2, partim)
    k = 3 * (lower - 1)
    @inbounds c1, c2, c3 = coef[k+1], coef[k+2], coef[k+3]
    c3 <= 0.0 && _sclk_badrate(clock.id)
    rate = 1.0 / (c3 / clock.tikmsc)
    sclkdp = _dnint(c1 + rate * (partim - c2))
    (sclkdp < _dnint(coef[1]) || sclkdp > clock.mxtick || isnan(sclkdp)) &&
        _sclk_range(clock.id, "the epoch $et")
    return sclkdp
end

# --- numeric interface ------------------------------------------------------------------

"`x` as a floating-point number: an integer becomes `Float64`, every other real is kept."
@inline _sclk_num(x::Integer) = Float64(x)
@inline _sclk_num(x::Real) = x

"""
    sct2e(sc, sclkdp) -> ET

Convert the encoded spacecraft clock `sclkdp` (ticks) of the clock `sc` to ET, as
`sct2e_c`. Generic over the number type of `sclkdp`; does not allocate. Forms: on a
`KernelSet`, on a `SpiceContext`, and on the default context.
"""
sct2e(ks::KernelSet, sc::Integer, sclkdp::Real) =
    _sct2e(ks, _sclk(ks, sc), _sclk_num(sclkdp))

"""
    sce2c(sc, et) -> ticks

Convert ET `et` to the continuous encoded spacecraft clock of the clock `sc`, as
`sce2c_c`. Generic over the number type of `et`; does not allocate.
"""
sce2c(ks::KernelSet, sc::Integer, et::Real) = _sce2c(ks, _sclk(ks, sc), _sclk_num(et))

"""
    sce2t(sc, et) -> Float64

Convert ET `et` to the discrete encoded spacecraft clock of the clock `sc` (whole
ticks), as `sce2t_c`. A number that is not a `Float64` is taken at its `primal` value.
"""
sce2t(ks::KernelSet, sc::Integer, et::Real) = _sce2t(ks, _sclk(ks, sc), primal(et))

"""
    scpart(sc) -> (pstart::Vector{Float64}, pstop::Vector{Float64})

Partition start and stop times of the clock `sc` in ticks, as `scpart_c`. The vectors
are copies; they have one element for each partition.
"""
function scpart(ks::KernelSet, sc::Integer)
    clock = _sclk(ks, sc)
    return copy(clock.pstart), copy(clock.pstop)
end

# --- strings ----------------------------------------------------------------------------

"""
    _lparsm(list, delims, nmax) -> Vector{SubString{String}}

Items of `list` separated by any character of `delims`, at most `nmax`, as `lparsm_`:
blanks around an item are dropped, two delimiters in a row give a blank item, a run of
blanks counts as one delimiter and joins a non-blank delimiter next to it, and a list
that ends with a delimiter has a last blank item. A blank list has one blank item.
"""
function _lparsm(list::AbstractString, delims::AbstractString, nmax::Int)
    s = String(rstrip(list, ' '))
    items = SubString{String}[]
    isempty(s) && return push!(items, SubString(s, 1, 0))
    c = codeunits(s)
    isdelim(k) = UInt8(' ') == c[k] || Char(c[k]) in delims
    eol = length(c)
    b = 1
    while b <= eol
        # invariant: `items` holds the items that end before position b
        while b <= eol && c[b] == UInt8(' ')
            b += 1
        end
        e = b
        while e <= eol && !isdelim(e)
            e += 1
        end
        stop = e - 1                     # last character of the item
        if e <= eol && c[e] == UInt8(' ')
            while e <= eol && c[e] == UInt8(' ')
                e += 1
            end
            e <= eol && !isdelim(e) && (e -= 1)
        end
        push!(items, SubString(s, b, stop))
        length(items) < nmax || return items
        b = e + 1
    end
    Char(c[eol]) in delims && length(items) < nmax && push!(items, SubString(s, 1, 0))
    return items
end

"Ticks of one count of each of the `nfield` fields of `clock` (`CMPTKS` of `sc01`)."
function _sclk_cmptks(clock::SclkType1)
    nfield = Int(clock.nfield)
    cmptks = ones(nfield)
    for i in nfield-1:-1:1          # invariant: cmptks[i+1:nfield] are final
        cmptks[i] = cmptks[i+1] * clock.moduli[i+1]
    end
    return cmptks
end

"""
    sctiks(sc, clkstr) -> Float64

Ticks of the clock string `clkstr` (no partition number) of the clock `sc`, as
`sctiks_c`: the fields, most significant first and separated by `.`, `:`, `-`, `,` or
blanks, less their offsets, times the ticks of one count of the field. A blank field
counts as its offset. Throws `SpiceError(:invalidsclkstring)` for a blank string, too
many fields, a field that is not a number, or a field below its offset.
"""
function sctiks(ks::KernelSet, sc::Integer, clkstr::AbstractString)
    clock = _sclk(ks, sc)
    isempty(strip(clkstr)) && throw(SpiceError(:invalidsclkstring, "the clock string is blank"))
    cmptks = _sclk_cmptks(clock)
    items = _lparsm(clkstr, _SCLK_DELIMS, _SCLK_MXNFLD)
    length(items) > clock.nfield && throw(SpiceError(:invalidsclkstring,
        "the clock string '$clkstr' has $(length(items)) fields, which is too many"))
    ticks = 0.0
    for (i, item) in enumerate(items)       # invariant: ticks = sum over the fields 1:i-1
        value = if isempty(strip(item))
            clock.offsets[i]
        else
            try
                _nparsd(strip(item))
            catch err
                err isa SpiceError || rethrow()
                throw(SpiceError(:invalidsclkstring,
                    "could not parse SCLK component '$item' of '$clkstr' as a number"))
            end
        end
        value -= clock.offsets[i]
        _dnint(value) < 0 && throw(SpiceError(:invalidsclkstring,
            "component number $i of the SCLK string '$clkstr' is less than its field " *
            "offset $(clock.offsets[i])"))
        ticks += value * cmptks[i]
    end
    return ticks
end

"""
    scfmt(sc, ticks) -> String

Clock string (no partition number) of the tick count `ticks` of the clock `sc`, as
`scfmt_c`: `ticks` is rounded, split into the clock fields, and each field is written
with leading zeros to the width of its largest value, with the output delimiter of the
clock between the fields. Throws `SpiceError(:valueoutofrange)` for negative ticks.
"""
function scfmt(ks::KernelSet, sc::Integer, ticks::Real)
    clock = _sclk(ks, sc)
    nfield = Int(clock.nfield)
    cmptks = _sclk_cmptks(clock)
    rem = _dnint(Float64(ticks))
    (rem < 0 || !isfinite(rem)) && throw(SpiceError(:valueoutofrange,
        "negative or non-finite value for SCLK ticks: $ticks"))
    1 <= clock.delim <= 5 || throw(SpiceError(:valueoutofrange,
        "output delimiter code $(clock.delim) of the clock $sc is not in 1:5"))
    delim = _SCLK_OUTPUT_DELIM[clock.delim]
    io = IOBuffer()
    for i in 1:nfield       # invariant: fields 1:i-1 are written, rem is what they leave
        value = if i < nfield
            v = trunc(rem / cmptks[i]) + clock.offsets[i]
            rem = Base.rem(rem, cmptks[i])
            v
        else
            rem + clock.offsets[i]
        end
        width = trunc(Int, log10(clock.moduli[i] + clock.offsets[i] - 1.0 + 0.5)) + 1
        i > 1 && print(io, delim)
        print(io, lpad(string(trunc(BigInt, abs(value))), width, '0'))
    end
    return String(take!(io))
end

"Rounded partition bounds and rounded running totals of `clock` (`scencd_`, `scdecd_`)."
function _sclk_partitions(clock::SclkType1)
    pstart = _dnint.(clock.pstart)
    pstop = _dnint.(clock.pstop)
    ptotls = similar(pstart)
    acc = 0.0
    for i in eachindex(pstart)      # invariant: acc = rounded length of partitions 1:i-1
        acc = _dnint(acc + pstop[i] - pstart[i])
        ptotls[i] = acc
    end
    return pstart, pstop, ptotls
end

"""
    scencd(sc, sclkch) -> Float64

Encode the clock string `sclkch` of the clock `sc` (`"partition/fields"`, the partition
number optional) as ticks from the start of the first partition, as `scencd_c`. Without
a partition number the first partition that holds the count is taken. Throws
`SpiceError` with kind `:badpartnumber`, `:notinpart`, `:nopartition` or
`:invalidsclkstring`.
"""
function scencd(ks::KernelSet, sc::Integer, sclkch::AbstractString)
    clock = _sclk(ks, sc)
    pstart, pstop, ptotls = _sclk_partitions(clock)
    nparts = length(pstart)
    pos = findfirst('/', sclkch)
    pos !== nothing && isempty(strip(SubString(sclkch, 1, prevind(sclkch, pos)))) &&
        throw(SpiceError(:badpartnumber,
            "unable to parse the partition number from the SCLK string '$sclkch'"))
    fields = pos === nothing ? sclkch : SubString(sclkch, nextind(sclkch, pos))
    isempty(strip(fields)) && throw(SpiceError(:invalidsclkstring,
        "no SCLK components follow the slash in the SCLK string '$sclkch'"))
    ticks = _dnint(sctiks(ks, sc, fields))
    part = 1
    if pos !== nothing
        x = try
            _nparsd(strip(SubString(sclkch, 1, prevind(sclkch, pos))))
        catch err
            err isa SpiceError || rethrow()
            throw(SpiceError(:badpartnumber,
                "unable to parse the partition number from the SCLK string '$sclkch'"))
        end
        (1 <= trunc(x) <= nparts) || throw(SpiceError(:badpartnumber,
            "partition number $(trunc(x)) of the SCLK string '$sclkch' is not in 1:$nparts"))
        part = trunc(Int, x)
        (ticks < pstart[part] || ticks > pstop[part]) && throw(SpiceError(:notinpart,
            "SCLK count $ticks is not inside partition number $part"))
    else
        while part <= nparts && (ticks < pstart[part] || ticks > pstop[part])
            part += 1                # invariant: partitions 1:part-1 do not hold `ticks`
        end
        part > nparts && throw(SpiceError(:nopartition,
            "SCLK count $ticks is not inside any partition of the spacecraft clock $sc"))
    end
    return part > 1 ? ticks - pstart[part] + ptotls[part-1] : ticks - pstart[part]
end

"""
    scdecd(sc, sclkdp) -> String

Clock string `"partition/fields"` of the encoded SCLK `sclkdp` (rounded to whole ticks)
of the clock `sc`, as `scdecd_c`. Throws `SpiceError(:valueoutofrange)` if the ticks
are not inside a partition.
"""
function scdecd(ks::KernelSet, sc::Integer, sclkdp::Real)
    clock = _sclk(ks, sc)
    pstart, _, ptotls = _sclk_partitions(clock)
    nparts = length(pstart)
    ticks = _dnint(Float64(sclkdp))
    part = ticks == ptotls[nparts] ? nparts : searchsortedlast(ptotls, ticks) + 1
    (ticks < 0 || part > nparts || isnan(ticks)) && throw(SpiceError(:valueoutofrange,
        "value for ticks, $ticks, is not inside any partition of the spacecraft clock $sc"))
    ticks = part == 1 ? ticks + pstart[1] : ticks + pstart[part] - ptotls[part-1]
    return string(part, '/', scfmt(ks, sc, ticks))
end

"""
    sce2s(sc, et) -> String

Clock string of the ET `et` for the clock `sc`, as `sce2s_c`: `scdecd` of `sce2t`.
"""
sce2s(ks::KernelSet, sc::Integer, et::Real) = scdecd(ks, sc, sce2t(ks, sc, et))

"""
    scs2e(sc, sclkch) -> Float64

ET of the clock string `sclkch` of the clock `sc`, as `scs2e_c`: `sct2e` of `scencd`.
"""
scs2e(ks::KernelSet, sc::Integer, sclkch::AbstractString) =
    sct2e(ks, sc, scencd(ks, sc, sclkch))

# --- forms on a context and on the default context ---------------------------------------

for f in (:sct2e, :sce2c, :sce2t, :sce2s, :scdecd, :scfmt)
    @eval $f(ctx::SpiceContext, sc::Integer, x::Real) = $f(kernelset(ctx), sc, x)
    @eval $f(sc::Integer, x::Real) = $f(default_context(), sc, x)
end
for f in (:scs2e, :scencd, :sctiks)
    @eval $f(ctx::SpiceContext, sc::Integer, s::AbstractString) = $f(kernelset(ctx), sc, s)
    @eval $f(sc::Integer, s::AbstractString) = $f(default_context(), sc, s)
end
for f in (:scpart, :sctype)
    @eval $f(ctx::SpiceContext, sc::Integer) = $f(kernelset(ctx), sc)
    @eval $f(sc::Integer) = $f(default_context(), sc)
end
