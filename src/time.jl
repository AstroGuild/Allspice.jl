# Time conversions: leap seconds, ET - UTC (`deltet.c`), ISO 8601 calendar strings
# to ET and ET to UTC strings (the paths of `str2et.c`, `et2utc.c`, `ttrans.c` and
# `unitim.c` that those conversions use).

const _SPD = 86400.0
const _HALFD = 43200.0
const _JD2000 = 2451545.0

# ---------------------------------------------------------------------------------
# Leap-second kernel
# ---------------------------------------------------------------------------------

"""
    build_lsk(pool) -> Union{Nothing,LeapSeconds}

Contract: returns the leap-second data when `pool` has all of `DELTET/DELTA_T_A`,
`DELTET/K`, `DELTET/EB`, `DELTET/M` (two values) and `DELTET/DELTA_AT` (a non-empty
list of pairs) as numeric variables, and `nothing` otherwise.
"""
function build_lsk(pool::AbstractDict{String,PoolValue})::Union{Nothing,LeapSeconds}
    num(key) = begin
        v = get(pool, key, nothing)
        v isa Vector{Float64} ? v : nothing
    end
    dta, k, eb = num("DELTET/DELTA_T_A"), num("DELTET/K"), num("DELTET/EB")
    m, dat = num("DELTET/M"), num("DELTET/DELTA_AT")
    (dta === nothing || k === nothing || eb === nothing || m === nothing ||
     dat === nothing) && return nothing
    (isempty(dta) || isempty(k) || isempty(eb) || length(m) < 2 || length(dat) < 2 ||
     isodd(length(dat))) && return nothing
    return LeapSeconds(dta[1], k[1], eb[1], (m[1], m[2]), copy(dat))
end

"Contract: returns the leap-second data of `ks`; throws `SpiceError(:nolsk)` if no leap-second kernel is loaded."
function _lsk(ks::KernelSet)::LeapSeconds
    l = ks.lsk
    l === nothing && throw(SpiceError(:nolsk,
        "no leap-second kernel is loaded; load one (for example naif0012.tls) with furnsh"))
    return l
end

"Contract: `b` is upper case; returns whether `a` equals `b` when case is ignored, without allocation."
_eqci(a::AbstractString, b::AbstractString) =
    ncodeunits(a) == ncodeunits(b) && all(((x, y),) -> uppercase(x) == y, zip(a, b))

"""
    deltet(ks, epoch, eptype) -> Float64

Contract: `eptype` is `"UTC"` (`epoch` is UTC seconds past J2000) or `"ET"` (`epoch`
is TDB seconds past J2000), in any case; returns ET - UTC in seconds at `epoch` with
the arithmetic of CSPICE `deltet`. Throws `SpiceError(:nolsk)` without a leap-second
kernel and `SpiceError(:badepochtype)` for another `eptype`. Reads `ks.lsk` only.
"""
function deltet(ks::KernelSet, epoch::Real, eptype::AbstractString)::Float64
    l = _lsk(ks)
    typ = strip(eptype)
    isutc = _eqci(typ, "UTC")
    isutc || _eqci(typ, "ET") || throw(SpiceError(:badepochtype,
        "the epoch type '$eptype' is neither 'UTC' nor 'ET'"))
    t = primal(epoch)
    dleap = l.delta_at
    dta, k, eb = l.delta_t_a, l.k, l.eb
    m0, m1 = l.m
    leaps = dleap[1] - 1.0
    # Loop invariant: `leaps` is TAI - UTC of the last pair before `i` whose epoch is
    # not after `t`, or one second less than the first value if there is none.
    @inbounds for i in 1:2:length(dleap)-1
        if isutc
            t >= dleap[i+1] && (leaps = dleap[i])
        elseif t > dleap[i+1]
            aet = round(dleap[i+1] + dta + dleap[i], RoundNearestTiesAway)
            ma = m0 + m1 * aet
            ea = ma + eb * sin(ma)
            et = dleap[i+1] + dta + dleap[i] + k * sin(ea)
            t >= et && (leaps = dleap[i])
        end
    end
    aet = isutc ? round(t + dta + leaps, RoundNearestTiesAway) :
          round(t, RoundNearestTiesAway)
    ma = m0 + m1 * aet
    ea = ma + eb * sin(ma)
    return dta + leaps + k * sin(ea)
end

# ---------------------------------------------------------------------------------
# Calendar and the day/TAI tables of TTRANS
# ---------------------------------------------------------------------------------

const _DPJAN0 = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
const _DPBEGL = (0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)

"Contract: `year >= 1`; returns 1 for a Gregorian leap year and 0 otherwise."
_leap(year::Int) = (year % 4 == 0) - (year % 100 == 0) + (year % 400 == 0)

"Contract: `year >= 1`, `1 <= month <= 12`; returns the number of days from 0001-01-01 (day 0) to the given Gregorian date, as the `DAYNUM` of `ttrans`."
function _daynum(year::Int, month::Int, day::Int)::Int
    y = year - 1
    return y * 365 + div(y, 4) - div(y, 100) + div(y, 400) +
           (_DPJAN0[month] + (month > 2 ? _leap(year) : 0) + day) - 1
end

const _DN2000 = _daynum(2000, 1, 1)
const _JD1101 = _JD2000 - Float64(_DN2000) - 0.5

"Contract: returns (year, month, day, day of year) of day number `daynum` in the Gregorian calendar, as `ttrans` computes them."
function _calendar(daynum::Int)
    yr400 = div(daynum, 146097)
    rem = daynum - yr400 * 146097
    if rem < 0
        yr400 -= 1
        rem += 146097
    end
    yr100 = min(3, div(rem, 36524))
    rem -= yr100 * 36524
    yr4 = min(24, div(rem, 1461))
    rem -= yr4 * 1461
    yr1 = min(3, div(rem, 365))
    rem -= yr1 * 365
    dofyr = rem + 1
    year = yr400 * 400 + yr100 * 100 + yr4 * 4 + yr1 + 1
    tab = (mod(year, 4) == 0) - (mod(year, 100) == 0) + (mod(year, 400) == 0) == 0 ?
          _DPJAN0 : _DPBEGL
    month = 1
    while month < 12 && tab[month+1] < dofyr
        month += 1
    end
    return year, month, dofyr - tab[month], dofyr
end

"Contract: `1 <= i <= length(l.delta_at)`; entry `i` of the `TAITAB` of `ttrans`: TAI seconds past J2000 at the start of the day before (odd `i`) and of the day after (even `i`) each leap-second epoch."
@inline function _taitab(l::LeapSeconds, i::Int)::Float64
    d = l.delta_at
    j = (i + 1) >> 1
    @inbounds begin
        formal = d[2j]
        if isodd(i)
            lastdt = j == 1 ? d[1] - 1.0 : d[2j-3]
            return formal - _SPD + lastdt
        else
            return formal + d[2j-1]
        end
    end
end

"Contract: `1 <= i <= length(l.delta_at)`; entry `i` of the `DAYTAB` of `ttrans`: the day number that starts at `_taitab(l, i)`."
@inline function _daytab(l::LeapSeconds, i::Int)::Int
    j = (i + 1) >> 1
    daynum = trunc(Int, (@inbounds(l.delta_at[2j]) + _HALFD) / _SPD) + _DN2000
    return isodd(i) ? daynum - 1 : daynum
end

"Contract: returns the largest `i` with `_taitab(l, i) <= tai`, or 0 (`lstled`)."
function _last_tai(l::LeapSeconds, tai::Float64)::Int
    n = length(l.delta_at)
    i = 0
    while i < n && _taitab(l, i + 1) <= tai
        i += 1
    end
    return i
end

"Contract: returns the largest `i` with `_daytab(l, i) <= daynum`, or 0 (`lstlei`)."
function _last_day(l::LeapSeconds, daynum::Int)::Int
    n = length(l.delta_at)
    i = 0
    while i < n && _daytab(l, i + 1) <= daynum
        i += 1
    end
    return i
end

"Contract: `den > 0`; returns (q, r) with `num == q*den + r` and `0 <= r < den` up to rounding, as CSPICE `rmaind`."
@inline function _rmaind(num::Float64, den::Float64)
    q = trunc(num / den)
    r = num - q * den
    if r < 0.0
        q -= 1.0
        r += den
    end
    return q, r
end

"Contract: returns (day number, seconds into that UTC day) of the TAI epoch `tai`; the seconds reach 86400 and above only inside a leap second."
function _tai_to_day(l::LeapSeconds, tai::Float64)
    p = _last_tai(l, tai)
    if isodd(p)
        return _daytab(l, p), tai - _taitab(l, p)
    end
    p = max(p, 1)
    q, secs = _rmaind(tai - _taitab(l, p), _SPD)
    return trunc(Int, q) + _daytab(l, p), secs
end

"Contract: returns the TAI seconds past J2000 of `secs` seconds into UTC day `daynum`."
function _day_to_tai(l::LeapSeconds, daynum::Int, secs::Float64)::Float64
    p = max(1, _last_day(l, daynum))
    return _taitab(l, p) + (secs + Float64(daynum - _daytab(l, p)) * _SPD)
end

"Contract: the step of `ttrans` for a non-formal input: (day, seconds) with the seconds outside `[0, 86399]` are moved to the day that holds the instant."
function _regularize(l::LeapSeconds, daynum::Int, secs::Float64)
    (secs > _SPD - 1.0 || secs < 0.0) || return daynum, secs
    return _tai_to_day(l, _day_to_tai(l, daynum, secs))
end

"Contract: TDB seconds past J2000 of the TAI epoch `tai` (`unitim`)."
function _tai2tdb(l::LeapSeconds, tai::Float64)::Float64
    tdt = tai + l.delta_t_a
    return _tdt2tdb(l, tdt)
end

"Contract: TDB seconds past J2000 of the TDT epoch `tdt` (`unitim`)."
function _tdt2tdb(l::LeapSeconds, tdt::Float64)::Float64
    m0, m1 = l.m
    return tdt + l.k * sin(m0 + m1 * tdt + l.eb * sin(m0 + m1 * tdt))
end

"Contract: TAI seconds past J2000 of the TDB epoch `tdb`, by the three fixed-point steps of `unitim`."
function _tdb2tai(l::LeapSeconds, tdb::Float64)::Float64
    m0, m1 = l.m
    tdt = tdb
    for _ in 1:3
        tdt = tdb - l.k * sin(m0 + m1 * tdt + l.eb * sin(m0 + m1 * tdt))
    end
    return tdt - l.delta_t_a
end

# ---------------------------------------------------------------------------------
# Strings to ET
# ---------------------------------------------------------------------------------

const _ISO_TIME = raw"(?:(?:T|\s)\s*(\d{1,2})(?::(\d{1,2})(?::(\d{1,2}(?:\.\d*)?))?)?)?"
const _ISO_LABEL = raw"\s*(Z|UTC|TDB|ET|TDT|TT)?\s*$"
const _ISO_YMD = Regex(raw"^\s*(\d{4})-(\d{1,2})-(\d{1,2})" * _ISO_TIME * _ISO_LABEL, "i")
const _ISO_YD = Regex(raw"^\s*(\d{4})-(\d{3})" * _ISO_TIME * _ISO_LABEL, "i")

"""
    str2et(ks, s) -> Float64

Contract: `s` is an ISO 8601 calendar string `YYYY-MM-DD` or `YYYY-DDD` (day of
year), optionally followed by `T` or a blank and `hh[:mm[:ss[.fff]]]`, and optionally
by a label `UTC`, `Z`, `TDB`, `ET`, `TDT` or `TT` (default UTC); the year is 1 or
later. Returns TDB seconds past J2000 as CSPICE `str2et` does for these forms; a UTC
leap second `23:59:60` is accepted. Throws `SpiceError(:badtime)` for another form or
an impossible date, and `SpiceError(:nolsk)` when the conversion needs a leap-second
kernel (every system but TDB) and none is loaded. Reads `ks.lsk` only.
"""
function str2et(ks::KernelSet, s::AbstractString)::Float64
    bad(why) = throw(SpiceError(:badtime, "cannot convert the time string '$s': $why"))
    local year::Int, month::Int, day::Int
    m = match(_ISO_YD, s)
    if m !== nothing
        year, month, day = parse(Int, m[1]), 1, parse(Int, m[2])
        year >= 1 || bad("the year is before 1 A.D.")
        1 <= day <= 365 + _leap(year) || bad("the day of year is out of range")
    else
        m = match(_ISO_YMD, s)
        m === nothing && bad("the form is not ISO 8601 (YYYY-MM-DD or YYYY-DDD, " *
                             "optional Thh:mm:ss.fff, optional UTC or TDB label)")
        year, month, day = parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3])
        year >= 1 || bad("the year is before 1 A.D.")
        1 <= month <= 12 || bad("the month is out of range")
        dim = month == 12 ? 31 : (month == 2 ? 28 + _leap(year) :
                                  _DPJAN0[month+1] - _DPJAN0[month])
        1 <= day <= dim || bad("the day of month is out of range")
    end
    g = length(m.captures) - 4          # time and label are the last four groups
    hh = m[g+1] === nothing ? 0.0 : Float64(parse(Int, m[g+1]))
    mm = m[g+2] === nothing ? 0.0 : Float64(parse(Int, m[g+2]))
    ss = m[g+3] === nothing ? 0.0 : _nparsd(m[g+3])
    label = m[g+4] === nothing ? "UTC" : uppercase(m[g+4])
    label == "Z" && (label = "UTC")
    uniform = label in ("TDB", "ET", "TDT", "TT")
    hh < 24.0 || bad("the hour is out of range")
    mm < 60.0 || bad("the minute is out of range")
    ss < (uniform ? 60.0 : 61.0) || bad("the second is out of range")
    daynum = _daynum(year, month, day)
    secs = hh * 3600.0 + mm * 60.0 + ss
    if uniform
        q, r = _rmaind(secs, _SPD)
        formal = Float64(daynum + trunc(Int, q) - _DN2000) * _SPD - _HALFD + r
        return label in ("TDB", "ET") ? formal : _tdt2tdb(_lsk(ks), formal)
    end
    l = _lsk(ks)
    daynum, secs = _regularize(l, daynum, secs)
    return _tai2tdb(l, _day_to_tai(l, daynum, secs))
end

"""
    utc2et(ks, s) -> Float64

Contract: as `str2et`: `s` is an ISO 8601 calendar string, UTC unless it has another
label; returns TDB seconds past J2000.
"""
utc2et(ks::KernelSet, s::AbstractString)::Float64 = str2et(ks, s)

# ---------------------------------------------------------------------------------
# ET to strings
# ---------------------------------------------------------------------------------

const _POW10 = (1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13,
    1e14, 1e15, 1e16, 1e17)
const _IPOW10 = (1.0, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11,
    1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17)

"Contract: `x > 0`; returns the 14 significant decimal digits of `x` and its decimal exponent, generated with the arithmetic of CSPICE `dpstr`."
function _dpstr14(x::Float64)
    copy = x
    exp10 = log10(copy)
    if exp10 < 0.0
        expont = trunc(Int, exp10) - 1
        k = -expont
        while k > 16
            copy *= 1e16
            k -= 16
        end
        k != 0 && (copy *= _POW10[k+1])
    else
        expont = trunc(Int, exp10)
        k = expont
        while k > 16
            copy *= 1e-16
            k -= 16
        end
        k != 0 && (copy *= _IPOW10[k+1])
    end
    copy = (round(copy * _POW10[14], RoundNearestTiesAway) + 0.125) * _IPOW10[14]
    if copy >= 10.0
        copy *= 0.1
        expont += 1
    end
    digits = ntuple(_ -> 0, 14)
    # Loop invariant: digits 1..n-1 are final and `copy` holds the rest, scaled so
    # that digit n is its integer part.
    for n in 1:14
        d = clamp(trunc(Int, copy), 0, 9)
        digits = Base.setindex(digits, d, n)
        copy = (copy - d) * 10.0
    end
    return digits, expont
end

"Contract: `from <= to`; returns the characters of slots `from..to` of the virtual digit string of CSPICE `zzvsbstr` (negative slots are integer digits, slot 0 is the decimal point), rounded at slot `to`, and whether the rounding carried out of slot `from`."
function _vsbstr(digits::NTuple{14,Int}, expont::Int, from::Int, to::Int)
    index(slot) = slot < 0 ? expont + slot + 1 : expont + slot
    slot = to + 1 == 0 ? 1 : to + 1
    i = index(slot)
    incr = 0 <= i < 14 && digits[i+1] >= 5
    out = Vector{Char}(undef, to - from + 1)
    for slot in to:-1:from
        c = if slot == 0
            '.'
        else
            i = index(slot)
            if i < 0
                if incr
                    incr = false
                    '1'
                else
                    slot < -1 ? ' ' : '0'
                end
            elseif i < 14
                v = digits[i+1] + (incr ? 1 : 0)
                v == 10 ? '0' : (incr = false; Char('0' + v))
            else
                incr = false
                '0'
            end
        end
        out[slot-from+1] = c
    end
    return String(out), incr
end

"""
    _dpstrf(x, sigdig) -> String

Contract: returns `x` in fixed-point form with `clamp(sigdig, 1, 14)` significant
digits, as CSPICE `dpstrf(x, sigdig, 'F')` writes it: a blank or `-` for the sign,
then the digits; the rounding is the two-step rounding of CSPICE (to 14 digits, then
to the requested digit), so a tie rounds away from zero.
"""
function _dpstrf(x::Float64, sigdig::Integer)::String
    maxdig = Int(clamp(sigdig, 1, 14))
    sign = x < 0.0 ? "-" : " "
    x == 0.0 && return string(sign, "0.", '0'^maxdig)
    digits, expont = _dpstr14(abs(x))
    first = expont >= 0 ? -expont - 1 : -expont
    last = first + maxdig - 1
    first < 0 && last >= 0 && (last += 1)
    first = min(-1, first)
    str, overflow = _vsbstr(digits, expont, first, last)
    if overflow
        first -= 1
        str, _ = _vsbstr(digits, expont, first, last)
        last > 0 && (str = chop(str))
    end
    last < 0 && (str = string(str, '0'^(-1 - last), '.'))
    return string(sign, lstrip(str, ' '))
end

"Contract: `n >= 0` and `width >= 1`; returns the decimal digits of `n`, zero-padded on the left to `width`."
_pad(n::Integer, width::Int) = lpad(string(n), width, '0')

"""
    et2utc(ks, et, format, prec) -> String

Contract: `format` is `"C"` (`1986 APR 12 16:31:09.814`), `"D"`
(`1986-102 // 16:31:09.814`), `"ISOC"` (`1986-04-12T16:31:09.814`), `"ISOD"`
(`1986-102T16:31:09.814`) or `"J"` (`JD 2446533.18834276`), in any case; `prec` is
the number of decimal places and is limited to 0..14 (CSPICE writes at most 14
significant digits, so 13 places of seconds and 7 of a Julian date). Returns the UTC string of the
TDB epoch `et` as CSPICE `et2utc` builds it, with seconds 60 inside a leap second.
Throws `SpiceError(:badformat)` for another format, `SpiceError(:badtime)` for a year
before 1 A.D. in an ISO format, and `SpiceError(:nolsk)` without a leap-second
kernel. Reads `ks.lsk` only.
"""
function et2utc(ks::KernelSet, et::Real, format::AbstractString, prec::Integer)::String
    fmt = uppercase(strip(format))
    fmt in ("C", "D", "J", "ISOC", "ISOD") || throw(SpiceError(:badformat,
        "the format '$format' is not one of 'C', 'D', 'J', 'ISOC', 'ISOD'"))
    l = _lsk(ks)
    p = Int(clamp(prec, 0, 14))
    tai = _tdb2tai(l, primal(et))
    if fmt == "J"
        daynum, secs = _regularize(l, _tai_to_day(l, tai)...)
        if secs > _SPD
            daynum += 1
            secs = 0.0
        end
        jd = _JD1101 + Float64(daynum) + secs / _SPD
        return "JD" * _dpstrf(jd, p + 7)
    end
    whole = floor(tai)
    scale = round(10.0^p)
    frac = round(scale * (tai - whole), RoundNearestTiesAway)
    if frac == scale
        whole += 1.0
        frac = 0.0
    end
    daynum, secs = _regularize(l, _tai_to_day(l, whole)...)
    year, month, day, dofyr = _calendar(daynum)
    exsecs = max(0.0, secs - _SPD + 1.0)
    hours, rest = _rmaind(secs - exsecs, 3600.0)
    mins, tsecs = _rmaind(rest, 60.0)
    hh = round(Int, hours)
    mi = round(Int, mins)
    sc = round(Int, tsecs + exsecs)
    clock = string(_pad(hh, 2), ':', _pad(mi, 2), ':', _pad(sc, 2))
    if p > 0
        # CSPICE writes `1 + fraction` with `dpstrf`, so at most 13 places exist.
        fs = _dpstrf(frac / scale + 1.0, p + 1)
        dot = findfirst('.', fs)
        clock = string(clock, SubString(fs, dot, min(dot + p, lastindex(fs))))
    end
    calendar_style = fmt == "C" || fmt == "D"
    ystr = if year >= 1000
        string(year)
    elseif year > 0
        calendar_style ? string(year, " A.D.") : string(year)
    else
        calendar_style || throw(SpiceError(:badtime,
            "the year $(1 - year) B.C. cannot be written in an ISO format"))
        string(1 - year, " B.C.")
    end
    sep = (calendar_style && year < 1000) ? " " : "-"
    if fmt == "C"
        return string(ystr, ' ', _MONTH_ABBR[month], ' ', _pad(day, 2), ' ', clock)
    elseif fmt == "ISOC"
        return string(ystr, '-', _pad(month, 2), '-', _pad(day, 2), 'T', clock)
    elseif fmt == "D"
        return string(ystr, sep, _pad(dofyr, 3), " // ", clock)
    else
        return string(ystr, '-', _pad(dofyr, 3), 'T', clock)
    end
end
