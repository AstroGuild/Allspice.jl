# Time conversions, numeric part: leap seconds, ET - UTC (`deltet.c`), the uniform time
# scales (`unitim.c`), the time vectors (`ttrans.c`), the Julian and Gregorian
# calendars (`jul2gr.c`) and the constants. `timeparse.jl` has the strings to ET and
# `timeout.jl` the ET to strings.

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
    deltet([ctx,] epoch, eptype) -> Float64

Return the difference ET - UTC in seconds at `epoch`, as `deltet_c`.

`eptype` says how `epoch` is measured: `"UTC"` for UTC seconds past J2000 or `"ET"`
for seconds past J2000 TDB; letter case is ignored. A leap-second kernel must be
loaded.

Throws `SpiceError(:nolsk)` if no leap-second kernel is loaded and
`SpiceError(:badepochtype)` if `eptype` is not `"UTC"` or `"ET"`.

# Examples
```julia
furnsh("naif0012.tls")
deltet(0.0, "ET")       # seconds, at J2000
```

See also [`unitim`](@ref), [`str2et`](@ref).
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
    return _wrap32(_cint(q) + _daytab(l, p)), secs
end

"Contract: returns the TAI seconds past J2000 of `secs` seconds into UTC day `daynum`."
function _day_to_tai(l::LeapSeconds, daynum::Int, secs::Float64)::Float64
    p = max(1, _last_day(l, daynum))
    return _taitab(l, p) + (secs + Float64(_wrap32(daynum - _daytab(l, p))) * _SPD)
end

# ---------------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------------

"""
    j1900() -> Float64

Return the Julian date of 1899 DEC 31 12:00:00, 2415020.0, as `j1900_c`.
"""
j1900() = 2415020.0
"""
    j1950() -> Float64

Return the Julian date of 1950 JAN 01 00:00:00, 2433282.5, as `j1950_c`.
"""
j1950() = 2433282.5
"""
    j2000() -> Float64

Return the Julian date of 2000 JAN 01 12:00:00, 2451545.0, as `j2000_c`.
"""
j2000() = _JD2000
"""
    j2100() -> Float64

Return the Julian date of 2100 JAN 01 12:00:00, 2488070.0, as `j2100_c`.
"""
j2100() = 2488070.0
"""
    jyear() -> Float64

Return the number of seconds in a Julian year, 31557600.0, as `jyear_c`.
"""
jyear() = 31557600.0
"""
    tyear() -> Float64

Return the number of seconds in a tropical year, 31556925.9747, as `tyear_c`.
"""
tyear() = 31556925.9747
"""
    spd() -> Float64

Return the number of seconds in a day, 86400.0, as `spd_c`.
"""
spd() = _SPD
"""
    b1900() -> Float64

Return the Julian date of the Besselian epoch B1900.0, 2415020.31352, as `b1900_c`.
"""
b1900() = 2415020.31352
"""
    b1950() -> Float64

Return the Julian date of the Besselian epoch B1950.0, 2433282.42345905, as `b1950_c`.
"""
b1950() = 2433282.42345905
"""
    clight() -> Float64

Return the speed of light in vacuum in km/s, 299792.458, as `clight_c`.
"""
clight() = 299792.458
"""
    halfpi() -> Float64

Return `pi / 2`, as `halfpi_c`.
"""
halfpi() = Float64(pi) / 2
"""
    twopi() -> Float64

Return `2 pi`, as `twopi_c`.
"""
twopi() = 2 * Float64(pi)
"""
    rpd() -> Float64

Return the number of radians in a degree, `pi / 180`, as `rpd_c`.
"""
rpd() = Float64(pi) / 180.0
"""
    dpr() -> Float64

Return the number of degrees in a radian, `180 / pi`, as `dpr_c`.
"""
dpr() = 180.0 / Float64(pi)

# ---------------------------------------------------------------------------------
# Integer conversions of C and the calendar arithmetic of TTRANS and JUL2GR
# ---------------------------------------------------------------------------------

"Contract: returns `x` cut to an integer as the C cast `(integer) x` does on x86: toward zero, and the smallest 32-bit integer when `x` is NaN or out of range."
@inline _cint(x::Float64)::Int =
    (x > -2147483649.0 && x < 2147483648.0) ? Int(trunc(Int32, x)) : Int(typemin(Int32))

"Contract: returns `n` reduced to the range of a 32-bit integer, as integer arithmetic of the C code does."
@inline _wrap32(n::Int)::Int = Int(n % Int32)

"Contract: returns the integer nearest to `x`, a tie away from zero, as f2c `i_dnnt`."
@inline _nint(x::Float64)::Int = _cint(x >= 0.0 ? floor(x + 0.5) : -floor(0.5 - x))

"Contract: returns `x` rounded to a whole number, a tie away from zero, as f2c `d_nint`."
@inline _dnint(x::Float64)::Float64 = x >= 0.0 ? floor(x + 0.5) : -floor(0.5 - x)

"Contract: `1 <= month <= 12`, any `year`; returns the day number (0001-01-01 is day 0) of the date in the extended Gregorian calendar, as the `DAYNUM` statement function of `ttrans` with its shift of years before 1 by whole 400-year cycles."
function _tt_daynum(year::Int, month::Int, day::Int)::Int
    doffst = 0
    if year <= 0
        yr400, year = fldmod(year, 400)
        if year == 0
            year += 400
            yr400 -= 1
        end
        doffst = yr400 * 146097
    end
    y = year - 1
    n = y * 365 + div(y, 4) - div(y, 100) + div(y, 400) +
        (_DPJAN0[month] + (month > 2 ? _leap(year) : 0) + day) - 1 + doffst
    return Int(n % Int32)
end

"Day number of 1991 JAN 06, a Sunday (`SUNDAY` of `ttrans`)."
const _SUNDAY = _tt_daynum(1991, 1, 6)

"Contract: `1 <= dofyr <= 366`; returns (month, day of month) of day `dofyr` of a year that has (`leap`) or does not have a leap day."
@inline function _month_day(dofyr::Int, leap::Bool)
    tab = leap ? _DPBEGL : _DPJAN0
    month = 1
    # Loop invariant: `tab[month] < dofyr`.
    while month < 12 && tab[month+1] < dofyr
        month += 1
    end
    return month, dofyr - tab[month]
end

const _RDAYG = _tt_daynum(1582, 10, 15)
const _RDAYJ = (1582 - 1) * 365 + div(1582 - 1, 4) + (_DPJAN0[10] + 0 + 5) - 1
const _OFFSTJ = _RDAYJ - _RDAYG
const _OFFSTG = _RDAYG - _RDAYJ

"""
    jul2gr(year, month, day) -> (year, month, day, doy)

Return the Gregorian calendar date, and its day of year, of a date given in the
Julian calendar, as `jul2gr_` of CSPICE.

Years are astronomical: 1 B.C. is 0 and 2 B.C. is -1. A `month` outside 1:12 is
reduced into a year. The Gregorian calendar is extended to all dates. The result is a
tuple of integers; CSPICE changes its arguments in place.

# Examples
```jldoctest
julia> jul2gr(1582, 10, 5)    # the day after 1582 OCT 4 of the Julian calendar
(1582, 10, 15, 288)
```

See also [`gr2jul`](@ref).
"""
function jul2gr(year::Integer, month::Integer, day::Integer)
    yroff, mon = fldmod(Int(month) - 1, 12)
    yr = Int(year) + yroff
    mon += 1
    offset = 0
    if yr <= 0
        m4, yr = fldmod(yr, 4)
        if yr == 0
            yr += 4
            m4 -= 1
        end
        offset = m4 * 1461
    end
    jleap = yr % 4 == 0 ? 1 : 0
    dayj = (yr - 1) * 365 + div(yr - 1, 4) + (_DPJAN0[mon] + (mon > 2 ? jleap : 0) + Int(day)) -
           1 + offset
    m400, dayg = fldmod(dayj + _OFFSTG, 146097)
    m100 = min(3, div(dayg, 36524))
    dayg -= m100 * 36524
    m4 = min(24, div(dayg, 1461))
    dayg -= m4 * 1461
    m = min(3, div(dayg, 365))
    dayg -= m * 365
    dofyr = dayg + 1
    yr = m400 * 400 + m100 * 100 + m4 * 4 + m + 1
    mon, dy = _month_day(dofyr, _leap(yr) != 0)
    return yr, mon, dy, dofyr
end

"""
    gr2jul(year, month, day) -> (year, month, day, doy)

Return the Julian calendar date, and its day of year, of a date given in the extended
Gregorian calendar, as `gr2jul_` of CSPICE.

This is the inverse of [`jul2gr`](@ref), with the same conventions.

# Examples
```jldoctest
julia> gr2jul(1582, 10, 15)
(1582, 10, 5, 278)
```
"""
function gr2jul(year::Integer, month::Integer, day::Integer)
    yroff, mon = fldmod(Int(month) - 1, 12)
    yr = Int(year) + yroff
    mon += 1
    dayg = _tt_daynum(yr, mon, Int(day))
    m4, dayj = fldmod(dayg + _OFFSTJ, 1461)
    m = min(3, div(dayj, 365))
    dayj -= m * 365
    dofyr = dayj + 1
    yr = m4 * 4 + m + 1
    mon, dy = _month_day(dofyr, yr % 4 == 0)
    return yr, mon, dy, dofyr
end

# ---------------------------------------------------------------------------------
# Uniform time scales: UNITIM
# ---------------------------------------------------------------------------------

# Codes of the uniform scales; 1..5 are of the TDT kind and 6..9 of the TDB kind.
const _US_JDTDT = 1
const _US_TAI = 2
const _US_TDT = 3
const _US_GPS = 4
const _US_TT = 5
const _US_ET = 6
const _US_JDTDB = 7
const _US_JED = 8
const _US_TDB = 9

"Contract: `b` is upper case ASCII; returns whether `a` equals `b` when case and the trailing blanks of `a` are ignored, without allocation."
function _eqsys(a::AbstractString, b::String)::Bool
    ca = codeunits(a)
    cb = codeunits(b)
    n = length(cb)
    length(ca) >= n || return false
    @inbounds for i in 1:n
        c = ca[i]
        (UInt8('a') <= c <= UInt8('z')) && (c -= 0x20)
        c == cb[i] || return false
    end
    @inbounds for i in n+1:length(ca)
        ca[i] == UInt8(' ') || return false
    end
    return true
end

"Contract: returns the code of the uniform time scale `s` (any case, trailing blanks allowed), or 0 if `s` is not one of them. Does not allocate."
function _unitim_code(s::AbstractString)::Int
    _eqsys(s, "TDB") && return _US_TDB
    _eqsys(s, "ET") && return _US_ET
    _eqsys(s, "TDT") && return _US_TDT
    _eqsys(s, "TT") && return _US_TT
    _eqsys(s, "TAI") && return _US_TAI
    _eqsys(s, "JED") && return _US_JED
    _eqsys(s, "JDTDB") && return _US_JDTDB
    _eqsys(s, "JDTDT") && return _US_JDTDT
    _eqsys(s, "GPS") && return _US_GPS
    return 0
end

"Contract: `i`, `o` are codes of uniform scales; returns `epoch` of scale `i` in scale `o` with the arithmetic of CSPICE `unitim`. Generic over the number type; does not allocate."
@inline function _unitim(l::LeapSeconds, epoch::Real, i::Int, o::Int)
    t = epoch * 1.0
    i == o && return t
    dta, k, eb = l.delta_t_a, l.k, l.eb
    m0, m1 = l.m
    if i == _US_TAI
        t += dta
    elseif i == _US_GPS
        t += dta + 19.0
    elseif i == _US_JDTDT || i == _US_JED || i == _US_JDTDB
        t = (t - _JD2000) * _SPD
    end
    intdt = i <= _US_TT
    outtdt = o <= _US_TT
    if intdt && !outtdt
        t = t + k * sin(m0 + m1 * t + eb * sin(m0 + m1 * t))
    elseif !intdt && outtdt
        tdb = t
        # Loop invariant: `t` is the TDT estimate after the given number of fixed-point
        # steps; CSPICE takes exactly three.
        for _ in 1:3
            t = tdb - k * sin(m0 + m1 * t + eb * sin(m0 + m1 * t))
        end
    end
    if o == _US_TAI
        t -= dta
    elseif o == _US_GPS
        t -= dta + 19.0
    elseif o == _US_JDTDT || o == _US_JED || o == _US_JDTDB
        t = t / _SPD + _JD2000
    end
    return t
end

"""
    unitim([ctx,] epoch, insys, outsys)

Return the epoch `epoch` of the uniform time scale `insys` expressed in the scale
`outsys`, as `unitim_c`.

| Scale                    | Epoch is                                  |
|:-------------------------|:------------------------------------------|
| `"TAI"`                  | TAI seconds past J2000 TAI                |
| `"TDT"`, `"TT"`          | TDT seconds past J2000 TDT                |
| `"TDB"`, `"ET"`          | TDB seconds past J2000 TDB                |
| `"GPS"`                  | GPS seconds, on the scale TAI - 19 s      |
| `"JDTDT"`                | Julian date in TDT, days                  |
| `"JDTDB"`, `"JED"`       | Julian date in TDB, days                  |

Letter case is ignored. UTC is not a uniform scale: use [`str2et`](@ref) and
[`et2utc`](@ref) for it. A leap-second kernel must be loaded; the call reads its
constants for TDB - TDT and for TDT - TAI. The result has the number type of
`epoch * 1.0`, so a dual-number epoch keeps its derivatives. The call does not
allocate.

Throws `SpiceError(:badtimetype)` if a scale is not one of the above and
`SpiceError(:nolsk)` if no leap-second kernel is loaded.

# Examples
```julia
furnsh("naif0012.tls")
jd = unitim(et, "TDB", "JDTDB")     # Julian date, TDB
tt = unitim(et, "TDB", "TT")        # seconds past J2000, TT
```

See also [`deltet`](@ref).
"""
function unitim(ks::KernelSet, epoch::Real, insys::AbstractString, outsys::AbstractString)
    l = _lsk(ks)
    i = _unitim_code(insys)
    o = _unitim_code(outsys)
    (i == 0 || o == 0) && _unitim_error(insys, outsys)
    return _unitim(l, epoch, i, o)
end

@noinline _unitim_error(insys, outsys) = throw(SpiceError(:badtimetype,
    "one of the time scales '$insys' and '$outsys' is not a uniform scale; the scales " *
    "are TAI, TDT, TT, GPS, JDTDT, TDB, ET, JDTDB and JED"))

# ---------------------------------------------------------------------------------
# Time vectors: TTRANS
# ---------------------------------------------------------------------------------

# Codes of the representations of `ttrans`, in the order of its sorted name table.
const _TT_NAMES = ("DAYSEC", "DP2000", "ET", "FORMAL", "JDTDB", "JDTDT", "JDUTC", "JED",
    "TAI", "TDB", "TDT", "YD", "YD.D", "YD.DF", "YDF", "YMD", "YMDF", "YMWD", "YMWDF", "YWD",
    "YWDF")
const _TT_DAYSEC = 1
const _TT_DP2000 = 2
const _TT_ET = 3
const _TT_FORMAL = 4
const _TT_JDTDB = 5
const _TT_JDTDT = 6
const _TT_JDUTC = 7
const _TT_JED = 8
const _TT_TAI = 9
const _TT_TDB = 10
const _TT_TDT = 11
const _TT_YD = 12
const _TT_YDD = 13
const _TT_YDDF = 14
const _TT_YDF = 15
const _TT_YMD = 16
const _TT_YMDF = 17
const _TT_YMWD = 18
const _TT_YMWDF = 19
const _TT_YWD = 20
const _TT_YWDF = 21

"Contract: `1 <= p <= 21`; returns whether representation `p` is formal, that is, counts every day as 86400 seconds."
@inline _tt_formal(p::Int) = p == _TT_FORMAL || p == _TT_JDUTC || p == _TT_YDDF ||
                             p == _TT_YDF || p == _TT_YMDF || p == _TT_YMWDF || p == _TT_YWDF

"Contract: returns the `unitim` code of the uniform representation `p` of `ttrans`, or 0 if `p` is not uniform."
@inline function _tt_uniform(p::Int)::Int
    p == _TT_ET && return _US_ET
    p == _TT_JDTDB && return _US_JDTDB
    p == _TT_JDTDT && return _US_JDTDT
    p == _TT_JED && return _US_JED
    p == _TT_TAI && return _US_TAI
    p == _TT_TDB && return _US_TDB
    p == _TT_TDT && return _US_TDT
    return 0
end

"Contract: returns the code of the representation `name` of `ttrans` (first word, any case), or 0."
function _tt_code(name::AbstractString)::Int
    w = uppercase(first(split(name * " ")))
    i = findfirst(==(w), _TT_NAMES)
    return i === nothing ? 0 : i
end

const TVec = NTuple{7,Float64}

@inline _tvec(x::Vararg{Real,N}) where {N} = ntuple(i -> i <= N ? Float64(x[i]) : 0.0, Val(7))

"""
    _ttrans(l, pfrom, pto, tv) -> NTuple{7,Float64}

Contract: `pfrom` and `pto` are codes `_TT_*`; `tv` holds the components of an epoch
in representation `pfrom`. Returns the components of the same epoch in representation
`pto`, with the steps and the arithmetic of CSPICE `ttrans`; components that `pto`
does not have are copied from `tv` (as the in-place routine leaves them). Does not
allocate.
"""
function _ttrans(l::LeapSeconds, pfrom::Int, pto::Int, tv::TVec)::TVec
    nref = length(l.delta_at)
    daynum = 0
    secs = 0.0
    if pfrom == _TT_YMD || pfrom == _TT_YMDF
        year, month, day = _cint(tv[1]), _cint(tv[2]), _cint(tv[3])
        dyear, month = fldmod(month - 1, 12)
        daynum = _tt_daynum(year + dyear, month + 1, day)
        secs = tv[4] * 3600.0 + tv[5] * 60.0 + tv[6]
    elseif pfrom == _TT_YD || pfrom == _TT_YDF
        daynum = _tt_daynum(_cint(tv[1]), 1, _cint(tv[2]))
        secs = tv[3] * 3600.0 + tv[4] * 60.0 + tv[5]
    elseif pfrom == _TT_YDD || pfrom == _TT_YDDF
        day = _cint(tv[2])
        frac = tv[2] - Float64(day)
        daynum = _tt_daynum(_cint(tv[1]), 1, day)
        if pfrom == _TT_YDDF
            secs = frac * _SPD
        else
            daylen = _SPD
            p = _last_day(l, daynum)
            isodd(p) && (daylen = _taitab(l, p + 1) - _taitab(l, p))
            secs = frac * daylen
        end
    elseif pfrom == _TT_FORMAL
        dp2000, secs = _rmaind(tv[1] + _HALFD, _SPD)
        daynum = _wrap32(_cint(dp2000) + _DN2000)
    elseif pfrom == _TT_JDUTC
        daydp, secs = _rmaind((tv[1] - _JD1101) * _SPD, _SPD)
        daynum = _cint(daydp)
    elseif pfrom == _TT_DAYSEC
        daynum = _cint(tv[1])
        secs = tv[2]
    elseif pfrom == _TT_DP2000
        daynum = _wrap32(_cint(tv[1]) + _DN2000)
        secs = tv[2]
    elseif pfrom == _TT_YWD || pfrom == _TT_YWDF
        week = _cint(tv[2]) - 1
        wkday = _cint(tv[3])
        daynum = _tt_daynum(_cint(tv[1]), 1, 1)
        fyrday = mod(daynum - _SUNDAY, 7) + 1
        daynum += week * 7 + mod(wkday - fyrday, 7)
        secs = tv[4] * 3600.0 + tv[5] * 60.0 + tv[6]
    elseif pfrom == _TT_YMWD || pfrom == _TT_YMWDF
        month = _cint(tv[2])
        (1 <= month <= 12) || throw(SpiceError(:badtime,
            "the month $month of a year-month-week-day time vector is not in 1..12"))
        week = _cint(tv[3]) - 1
        day = _cint(tv[4])
        daynum = _tt_daynum(_cint(tv[1]), month, 1)
        fmday = mod(daynum - _SUNDAY, 7) + 1
        daynum += week * 7 + mod(day - fmday, 7)
        secs = tv[5] * 3600.0 + tv[6] * 60.0 + tv[7]
    else
        uin = _tt_uniform(pfrom)
        uout = _tt_uniform(pto)
        uout != 0 && return Base.setindex(tv, _unitim(l, tv[1], uin, uout), 1)
        daynum, secs = _tai_to_day(l, _unitim(l, tv[1], uin, _US_TAI))
    end
    if _tt_formal(pfrom)
        daydp, secs = _rmaind(secs, _SPD)
        daynum = _wrap32(daynum + _cint(daydp))
    end
    # A formal input that goes to a formal output keeps its seconds; every other pair
    # moves seconds outside [0, 86399] to the day that holds the instant.
    if !(_tt_formal(pto) && _tt_formal(pfrom)) && (secs > _SPD - 1.0 || secs < 0.0)
        daynum, secs = _tai_to_day(l, _day_to_tai(l, daynum, secs))
    end
    if _tt_formal(pto) && secs > _SPD
        daynum += 1
        secs = 0.0
    end
    year = month = day = dofyr = 0
    hours = mins = tsecs = 0.0
    if pto >= _TT_YD
        year, month, day, dofyr = _calendar(daynum)
        exsecs = max(0.0, secs - _SPD + 1.0)
        hours, tempd = _rmaind(secs - exsecs, 3600.0)
        mins, tsecs = _rmaind(tempd, 60.0)
        tsecs += exsecs
    end
    if pto == _TT_YMD || pto == _TT_YMDF
        return (Float64(year), Float64(month), Float64(day), hours, mins, tsecs, tv[7])
    elseif pto == _TT_YD || pto == _TT_YDF
        return (Float64(year), Float64(dofyr), hours, mins, tsecs, tv[6], tv[7])
    elseif pto == _TT_YDD
        daylen = _SPD
        p = _last_day(l, daynum)
        isodd(p) && (daylen = _taitab(l, p + 1) - _taitab(l, p))
        return (Float64(year), Float64(dofyr) + secs / daylen, tv[3], tv[4], tv[5], tv[6], tv[7])
    elseif pto == _TT_YDDF
        return (Float64(year), Float64(dofyr) + secs / _SPD, tv[3], tv[4], tv[5], tv[6], tv[7])
    elseif pto == _TT_FORMAL
        return Base.setindex(tv, Float64(_wrap32(daynum - _DN2000)) * _SPD - _HALFD + secs, 1)
    elseif pto == _TT_JDUTC
        return Base.setindex(tv, _JD1101 + Float64(daynum) + secs / _SPD, 1)
    elseif pto == _TT_DAYSEC
        return (Float64(daynum), secs, tv[3], tv[4], tv[5], tv[6], tv[7])
    elseif pto == _TT_DP2000
        return (Float64(_wrap32(daynum - _DN2000)), secs, tv[3], tv[4], tv[5], tv[6], tv[7])
    elseif pto == _TT_YWD || pto == _TT_YWDF
        wday = mod(daynum - _SUNDAY, 7)
        week = div(dofyr - 1, 7) + 1
        return (Float64(year), Float64(week), Float64(wday) + 1.0, hours, mins, tsecs, tv[7])
    elseif pto == _TT_YMWD || pto == _TT_YMWDF
        week = div(daynum - _tt_daynum(year, month, 1), 7) + 1
        wday = mod(daynum - _SUNDAY, 7)
        return (Float64(year), Float64(month), Float64(week), Float64(wday) + 1.0, hours, mins,
            tsecs)
    else
        tai = _day_to_tai(l, daynum, secs)
        return Base.setindex(tv, _unitim(l, tai, _US_TAI, _tt_uniform(pto)), 1)
    end
end

"""
    ttrans(ks, from, to, tvec) -> NTuple{7,Float64}

Contract: `from` and `to` are names of time-vector representations of CSPICE `ttrans`
(`"YMD"`, `"YMDF"`, `"YD"`, `"YDF"`, `"YD.D"`, `"YD.DF"`, `"YWD"`, `"YWDF"`, `"YMWD"`,
`"YMWDF"`, `"DAYSEC"`, `"DP2000"`, `"FORMAL"`, `"JDUTC"`, `"TAI"`, `"TDT"`, `"TDB"`,
`"ET"`, `"JED"`, `"JDTDB"`, `"JDTDT"`), in any case; `tvec` holds up to seven
components. Returns the seven components of the same epoch in representation `to`.
Throws `SpiceError(:badtimetype)` for another name and `SpiceError(:nolsk)` without a
leap-second kernel. Internal: not exported.
"""
function ttrans(ks::KernelSet, from::AbstractString, to::AbstractString, tv)
    l = _lsk(ks)
    pfrom = _tt_code(from)
    pto = _tt_code(to)
    pfrom == 0 && throw(SpiceError(:badtimetype,
        "the time-vector representation '$from' is not one that ttrans knows"))
    pto == 0 && throw(SpiceError(:badtimetype,
        "the time-vector representation '$to' is not one that ttrans knows"))
    return _ttrans(l, pfrom, pto, _tvec(tv...))
end
