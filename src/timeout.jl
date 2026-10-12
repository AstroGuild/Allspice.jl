# ET to time strings: ports of `dpstr.c`, `zzvstrng.c`, `dpstrf.c` and `dpfmt.c` (the
# number formats of CSPICE, with their two-step rounding), `scanit.c` (the picture
# scanner), `timout.c`, `et2utc.c` and `etcal.c`.

# ---------------------------------------------------------------------------------
# Numbers to strings
# ---------------------------------------------------------------------------------

const _POW10 = (1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12, 1e13,
    1e14, 1e15, 1e16, 1e17)
const _IPOW10 = (1.0, 0.1, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11,
    1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17)

"Contract: `x >= 0`, `1 <= maxsig <= 14`; returns the first `maxsig` significant decimal digits of `x` (in a tuple of 14, the rest 0) and its decimal exponent, generated with the arithmetic of CSPICE `dpstr`. Zero gives digits 0 and exponent 0."
function _dpstr_digits(x::Float64, maxsig::Int=14)
    digits = ntuple(_ -> 0, 14)
    x == 0.0 && return digits, 0
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
    copy = (round(copy * _POW10[maxsig], RoundNearestTiesAway) + 0.125) * _IPOW10[maxsig]
    if copy >= 10.0
        copy *= 0.1
        expont += 1
    end
    # Loop invariant: digits 1..n-1 are final and `copy` holds the rest, scaled so
    # that digit n is its integer part.
    for n in 1:maxsig
        d = clamp(trunc(Int, copy), 0, 9)
        digits = Base.setindex(digits, d, n)
        copy = (copy - d) * 10.0
    end
    return digits, expont
end

"""
    _dpstr(x, sigdig) -> String

Contract: returns `x` in the exponent form of CSPICE `dpstr` with `clamp(sigdig, 1,
14)` significant digits: a blank or `-`, one digit, the decimal point, the other
digits, `E+` or `E-` and the exponent with two digits at least.
"""
function _dpstr(x::Float64, sigdig::Integer)::String
    maxsig = Int(clamp(sigdig, 1, 14))
    x == 0.0 && return string(" 0.", '0'^(maxsig - 1), "E+00")
    digits, expont = _dpstr_digits(abs(x), maxsig)
    io = IOBuffer()
    print(io, x < 0.0 ? '-' : ' ', digits[1], '.')
    for n in 2:maxsig
        print(io, digits[n])
    end
    print(io, expont >= 0 ? "E+" : "E-", lpad(abs(expont), 2, '0'))
    return String(take!(io))
end

"""
    _vsbstr(digits, expont, from, to, fill=' ', rnd=true) -> (String, Bool)

Contract: `from <= to`; `digits`, `expont` come from `_dpstr_digits(x, 14)`. Returns
the characters of slots `from..to` of the virtual digit string of CSPICE `zzvsbstr`
(negative slots are integer digits, slot 0 is the decimal point, positive slots are
decimals), rounded at slot `to` if `rnd`, with `fill` in the integer slots to the left
of the number, and whether the rounding carried out of slot `from`.
"""
function _vsbstr(digits::NTuple{14,Int}, expont::Int, from::Int, to::Int, fill::Char=' ',
        rnd::Bool=true)
    index(slot) = slot < 0 ? expont + slot + 1 : expont + slot
    incr = false
    if rnd
        slot = to + 1 == 0 ? 1 : to + 1
        i = index(slot)
        incr = 0 <= i < 14 && digits[i+1] >= 5
    end
    out = Vector{Char}(undef, to - from + 1)
    # Loop invariant: slots `slot+1..to` are written and `incr` is the carry into `slot`.
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
                    slot < -1 ? fill : '0'
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
    digits, expont = _dpstr_digits(abs(x))
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

"""
    _dpfmt(x, pictur) -> String

Contract: `pictur` is a number picture of CSPICE `dpfmt` without blanks: an optional
`+` or `-`, then one character for each integer digit (a leading `0` asks for zero
fill), and optionally a decimal point and one character for each decimal. Returns `x`
written to the picture, rounded at the last decimal, in exactly `length(pictur)`
characters: right justified, the minus sign next to the first digit, the exponent form
of `dpstr` when the integer part does not fit and `*` characters when that does not
fit either.
"""
function _dpfmt(x::Float64, pictur::AbstractString)::String
    p = codeunits(pictur)
    lastch = length(p)
    stars = '*'^lastch
    lastch == 0 && return ""
    sign = x > 0.0 ? '+' : (x < 0.0 ? '-' : ' ')
    sprsiz = x < 0.0 ? 1 : 0
    if p[1] == UInt8('+')
        needsn, sgnlen = true, 1
    elseif p[1] == UInt8('-')
        needsn, sgnlen = true, 1
        x > 0.0 && (sign = ' ')
    else
        x > 0.0 && (sign = ' ')
        needsn, sgnlen = false, sprsiz
    end
    start = needsn ? 2 : 1
    fill = start <= lastch && p[start] == UInt8('0') ? '0' : ' '
    dpat = something(findfirst(==(UInt8('.')), p), 0)
    intlen, declen = dpat > 0 ? (dpat - start, lastch - dpat) : (lastch - start + 1, -1)
    shift = false
    if start == 1 && sgnlen == 1
        intlen -= 1
        shift = true
        if intlen == -1
            intlen = 0
            declen -= 1
            declen == 0 && return stars
        end
    end
    y = abs(x)
    digits, expont = _dpstr_digits(y)
    # Exponent form, right justified in the width of the picture.
    expform() = begin
        sigdig = lastch - sgnlen - (min(1, div(expont, 1000)) + 4 + min(1, div(expont, 100))) - 1
        sigdig < 1 && return stars
        e = _dpstr(x, sigdig)
        e = lstrip(string(sign, SubString(e, 2)), ' ')
        ncodeunits(e) >= lastch ? String(first(e, lastch)) : lpad(e, lastch)
    end
    (expont >= intlen && y != 0.0) && return expform()
    ((intlen == 0 && declen == 0) || -intlen > declen) && return stars
    body, ovflow = _vsbstr(digits, expont, -intlen, declen, fill, true)
    ovflow && return expform()
    str = collect(sgnlen == 1 ? string(sign, body) : body)
    if shift
        frstch = something(findfirst(c -> c != ' ' && c != '-', str), 0)
        if frstch > 2
            str[frstch-1] = str[1]
            str[1] = ' '
        end
    end
    return rpad(String(str), lastch)
end

"Contract: `n >= 0` and `width >= 1`; returns the decimal digits of `n`, zero-padded on the left to `width`."
_pad(n::Integer, width::Int) = lpad(string(n), width, '0')

# ---------------------------------------------------------------------------------
# TIMOUT
# ---------------------------------------------------------------------------------

# Marks of a `timout` picture. The number of a mark is its slot in the `ID` array of
# `timout.c`; slot 1 is text that is not a mark.
const _TO_MARKS = ("YYYY" => 2, "YR" => 3, "MON" => 4, "Mon" => 5, "mon" => 6, "MONTH" => 7,
    "Month" => 8, "month" => 9, "MM" => 10, "DOY" => 11, "WKD" => 12, "Wkd" => 13, "wkd" => 14,
    "WEEKDAY" => 15, "Weekday" => 16, "weekday" => 17, "DD" => 18, "MN" => 19, "HR" => 20,
    "SC" => 21, ".#" => 22, "#" => 23, "JULIAND" => 24, "::UTC" => 25, "::TDB" => 26,
    "::TDT" => 27, "::TT" => 28, "SP2000" => 29, "SP1950" => 30, "::RND" => 31, "::TRNC" => 32,
    "ERA" => 33, "era" => 34, "?ERA?" => 35, "?era?" => 36, "AMPM" => 37, "ampm" => 38,
    "::UTC+" => 39, "::UTC-" => 40, "::JCAL" => 41, "::GCAL" => 42, "::MCAL" => 43, "AP" => 46)

"The marks in descending ASCII order: the order in which `scan` tries the marks that start with one character, so that the longest of `::UTC+`, `::UTC` or `MONTH`, `MON` wins."
const _TO_SCAN = sort!(collect(_TO_MARKS); by=first, rev=true)

const _ID_TEXT = 1
const _ID_HASH = 23
const _ID_UTC = 25
const _ID_TDB = 26
const _ID_TDT = 27
const _ID_TT = 28
const _ID_RND = 31
const _ID_UTCP = 39
const _ID_UTCM = 40
const _ID_JCAL = 41
const _ID_GCAL = 42
const _ID_MCAL = 43

# Classes of the marks (`CLASS` of `timout.c`); a component class is also the index of
# its value in `values`.
const _TC_TEXT = 2
const _TC_YEAR = 3
const _TC_YR = 4
const _TC_MM = 11
const _TC_DOY = 12
const _TC_DD = 19
const _TC_MN = 20
const _TC_HR = 21
const _TC_SC = 22
const _TC_POINT = 23
const _TC_PLACE = 24
const _TC_JD = 25
const _TC_SP2000 = 30
const _TC_SP1950 = 31
const _TC_SYSTEM = 45
const _TC_CAL = 46
const _TC_AP = 47
const _TC_MONTH = 48
const _TC_WKDAY = 49
const _TC_ERA = 50
const _TC_AMPM = 51
const _TC_SIGNED_YEAR = 52

"Contract: `1 <= id <= 46` is the number of a mark; returns its class."
function _to_class(id::Int)::Int
    id == 1 && return _TC_TEXT
    4 <= id <= 9 && return _TC_MONTH
    12 <= id <= 17 && return _TC_WKDAY
    (25 <= id <= 28 || id == 39 || id == 40) && return _TC_SYSTEM
    33 <= id <= 36 && return _TC_ERA
    (id == 37 || id == 38) && return _TC_AMPM
    41 <= id <= 43 && return _TC_CAL
    id == 46 && return _TC_AP
    return id + 1
end

"Contract: returns half of the unit of the component class `type` in seconds, the amount that `::RND` adds before truncation (`PAD` of `timout.c`); 0 for a class that is not a number."
function _to_pad(type::Int)::Float64
    type == _TC_SC && return 0.5
    type == _TC_MN && return 0.5 * 60.0
    (type == _TC_HR || type == _TC_AP) && return 0.5 * 60.0 * 60.0
    (type == _TC_DD || type == _TC_DOY || type == _TC_JD) && return 0.5 * 60.0 * 60.0 * 24.0
    type == _TC_MM && return 0.5 * 60.0 * 60.0 * 24.0 * 30.0
    (type == _TC_YEAR || type == _TC_YR) && return 0.5 * 60.0 * 60.0 * 24.0 * 365.0
    (type == _TC_SP2000 || type == _TC_SP1950) && return 0.5
    return 0.0
end

"Contract: `type` is a numeric component class; returns its number picture and the width of that picture (`ORIGNL` and `LENGTH` of `timout.c`)."
function _to_format(type::Int)
    type == _TC_YEAR && return "YYYY", 4
    type == _TC_YR && return "0Y", 2
    type == _TC_DOY && return "0DD", 3
    type == _TC_DD && return "0D", 2
    type == _TC_MM && return "0M", 2
    (type == _TC_HR || type == _TC_AP) && return "0H", 2
    type == _TC_MN && return "0M", 2
    type == _TC_SC && return "0S", 2
    type == _TC_JD && return "XXXXXXX", 7
    return "XXXXXXXXXXX", 11
end

const _WEEKDAY_NAMES = ("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
    "Saturday")

"""
    _to_scan(mystr) -> (ident, beg, fin)

Contract: returns the tokens of the picture `mystr` in order: for each token the
number of its mark (`_ID_TEXT` for a run of characters that start no mark) and its
first and last position (`scan` with the marks of `timout`). At each position the
marks are tried in descending ASCII order.
"""
function _to_scan(mystr::Vector{UInt8})
    ident = Int[]
    beg = Int[]
    fin = Int[]
    known = true
    start = 1
    n = length(mystr)
    # Loop invariant: `mystr[1:start-1]` is covered by the tokens, and `known` says
    # whether the last token is a mark.
    while start <= n
        hit = 0
        stop = start
        for (mark, id) in _TO_SCAN
            m = codeunits(mark)
            stop = start + length(m) - 1
            if stop <= n && @views(mystr[start:stop] == m)
                hit = id
                break
            end
        end
        if hit != 0
            push!(ident, hit); push!(beg, start); push!(fin, stop)
            known = true
            start = stop + 1
        else
            if known
                push!(ident, _ID_TEXT); push!(beg, start); push!(fin, start)
                known = false
            else
                fin[end] = start
            end
            start += 1
        end
    end
    return ident, beg, fin
end

"""
    timout(ks, et, pictur, defs=TimeDefaults()) -> String

Contract: `et` is TDB seconds past J2000 and `pictur` a format picture of CSPICE
`timout`. The marks are `YYYY`, `YR`, `MM`, `MON`/`Mon`/`mon`, `MONTH`/`Month`/
`month`, `DD`, `DOY`, `WKD`/`Wkd`/`wkd`, `WEEKDAY`/`Weekday`/`weekday`, `HR`, `AP`
(hour 1..12), `MN`, `SC`, `AMPM`/`ampm`, `ERA`/`era`, `?ERA?`/`?era?` (blank from the
year 1000 A.D.), `JULIAND`, `SP2000`, `SP1950`, `.#` and `#` (decimals of the
component to the left), and the marks that write nothing: `::UTC`, `::TDB`, `::TDT`,
`::TT`, `::UTC+h:mm`, `::UTC-h:mm` (time system or zone), `::GCAL`, `::JCAL`, `::MCAL`
(calendar), `::RND`, `::TRNC` (round or truncate the last component; truncation is
the default). Every other character is copied. A picture without a system mark uses
the system or the zone of `defs`, and one without a calendar mark `defs.calendar`.
Returns the string that CSPICE builds, without trailing blanks. Throws
`SpiceError(:nolsk)` without a leap-second kernel and `SpiceError(:badformat)` for
a picture of length zero, which the C interface of CSPICE rejects. Reads `ks.lsk` only.
"""
function timout(ks::KernelSet, et::Real, pictur::AbstractString,
        defs::TimeDefaults=TimeDefaults())::String
    l = _lsk(ks)
    isempty(pictur) && throw(SpiceError(:badformat, "timout: the picture is an empty " *
        "string (a blank picture gives a blank string)"))
    # The picture, left justified, with one blank after it: the blank closes the last
    # numeric component.
    mystr = Vector{UInt8}(codeunits(first(strip(pictur, ' '), 255)))
    push!(mystr, _BLANK)
    ident, beg, fin = _to_scan(mystr)
    ntokns = length(ident)

    # Time system: the default, or the first system mark of the picture.
    hoff = moff = 0.0
    dozone = false
    tsys = defs.system
    timtyp = _ID_UTC
    if tsys == "TDB"
        timtyp = _ID_TDB
    elseif tsys == "TDT" || tsys == "TT"
        timtyp = _ID_TDT
    elseif tsys != "UTC"
        hoff, moff, _, ok = _zzutcpm(codeunits(string("::", defs.zone)), 1)
        dozone = ok
        timtyp = ok ? _ID_UTCP : _ID_UTC
    end
    i = 1
    unknwn = true
    while unknwn && i <= ntokns
        if _to_class(ident[i]) == _TC_SYSTEM
            timtyp = ident[i]
            unknwn = false
            dozone = false
            if ident[i] == _ID_UTCP || ident[i] == _ID_UTCM
                hoff, moff, last, ok = _zzutcpm(mystr, beg[i])
                if ok
                    # The digits of the offset belong to the mark, not to the text.
                    dozone = true
                    timtyp = _ID_UTCP
                    if i < ntokns && last == fin[i+1]
                        ident[i+1] = ident[i]
                        i += 1
                    else
                        fin[i] = last
                        i < ntokns && (beg[i+1] = last + 1)
                    end
                end
            end
        end
        i += 1
    end

    # Calendar: the default, or the first calendar mark of the picture.
    caltyp = defs.calendar == "GREGORIAN" ? _ID_GCAL :
             (defs.calendar == "JULIAN" ? _ID_JCAL : _ID_MCAL)
    k = findfirst(id -> _to_class(id) == _TC_CAL, ident)
    k === nothing || (caltyp = ident[k])

    pumpup = _ID_RND in ident
    doera = any(id -> 33 <= id <= 36, ident)
    vanish = false

    # Remove the marks that write nothing.
    keep = [!(_to_class(id) == _TC_SYSTEM || _to_class(id) == _TC_CAL || id == _ID_RND ||
              id == 32) for id in ident]
    ident, beg, fin = ident[keep], beg[keep], fin[keep]
    ntokns = length(ident)
    class(i) = i <= ntokns ? _to_class(ident[i]) : _TC_TEXT

    # With `::RND`, add half of the smallest unit that the picture shows.
    timpad = 0.0
    if pumpup
        i = 1
        while i <= ntokns
            type = class(i)
            if type == _TC_TEXT || type == _TC_POINT || type == _TC_PLACE || type == _TC_AMPM ||
               type == _TC_ERA || type == _TC_MONTH || type == _TC_WKDAY
                i += 1
            else
                factor = 1.0
                incr = _to_pad(type)
                i += 1
                if class(i) == _TC_POINT
                    factor *= 0.1
                    i += 1
                    while i <= ntokns && ident[i] == _ID_HASH
                        factor *= 0.1
                        i += 1
                    end
                end
                incr *= factor
                timpad = timpad != 0.0 ? min(timpad, incr) : incr
            end
        end
    end

    values = zeros(Float64, 52)
    have = falses(52)
    myet = primal(et)
    sp1950 = _SPD * (j2000() - j1950())
    if timtyp == _ID_TDB
        myet += timpad
        values[_TC_SP2000] = myet
        values[_TC_JD] = _unitim(l, myet, _US_TDB, _US_JDTDB)
        bastyp, ymdfmt, ywfmt = _TT_FORMAL, _TT_YMDF, _TT_YMWDF
    elseif timtyp == _ID_TDT || timtyp == _ID_TT
        myet = _unitim(l, myet, _US_TDB, _US_TDT) + timpad
        values[_TC_SP2000] = myet
        values[_TC_JD] = _unitim(l, myet, _US_TDT, _US_JDTDT)
        bastyp, ymdfmt, ywfmt = _TT_FORMAL, _TT_YMDF, _TT_YMWDF
    else
        myet = _unitim(l, myet, _US_TDB, _US_TDT) + timpad
        bastyp, ymdfmt, ywfmt = _TT_TDT, _TT_YMD, _TT_YMWD
    end
    if bastyp == _TT_FORMAL
        values[_TC_SP1950] = values[_TC_SP2000] + sp1950
        have[_TC_SP2000] = have[_TC_SP1950] = have[_TC_JD] = true
    end

    out = fill(_BLANK, 256)
    appnd = 1
    # Write `text` at `appnd`; the Fortran assignment also blanks the rest of the line.
    put(text) = appnd <= 256 && _assign!(out, appnd, 256, text)
    making = false
    go2jul = false
    fmt = ""
    width = 0
    numtyp = _TC_YEAR
    timfmt = 1
    intyp = ymdfmt
    value = 0.0
    ptvec = zeros(Float64, 7)
    ntvec = zeros(Float64, 7)
    tup(v) = (v[1], v[2], v[3], v[4], v[5], v[6], v[7])

    # Loop invariant: tokens `1:i-1` are written to `out[1:appnd-1]`, except a numeric
    # component that is still collecting its decimals (`making`, with picture `fmt`).
    for i in 1:ntokns
        type = class(i)
        text = @view mystr[beg[i]:fin[i]]
        len = fin[i] - beg[i] + 1
        if type != _TC_PLACE && making && (type != _TC_POINT || have[_TC_POINT])
            # The component is complete: truncate it to its decimals and write it. The
            # value is moved one eighth of the last place below the cut so that the
            # rounding of `dpfmt` gives the truncated digits.
            trncat = clamp(width - _to_format(numtyp)[2] - 1, 0, 14)
            frac = value - trunc(value)
            if frac < 0.0
                value -= 1.0
                frac += 1.0
            end
            intmed = (trunc(frac * _POW10[trncat+1]) - 0.125) / _POW10[trncat+1]
            value = trunc(value) + clamp(intmed, 0.0, 1.0)
            put(_dpfmt(value, fmt))
            appnd += width
            have[_TC_POINT] = false
            making = false
        end
        if type == _TC_TEXT
            put(text)
            appnd += len
        elseif type == _TC_PLACE
            if making
                fmt = string(fmt, String(text))
                width += len
            else
                appnd <= 256 && (out[appnd] = text[1])
                appnd += len
            end
        elseif type == _TC_POINT
            if !making
                put(text)
                appnd += len
                have[_TC_POINT] = false
            else
                fmt = string(fmt, String(text))
                width += len
                have[_TC_POINT] = true
                if timfmt != 2
                    # Fraction of the unit: where the epoch lies between the start of
                    # this unit and the start of the next one.
                    pt = _ttrans(l, intyp, bastyp, tup(ptvec))
                    nt = _ttrans(l, intyp, bastyp, tup(ntvec))
                    delta = max(1.0, nt[1] - pt[1])
                    value += clamp((myet - pt[1]) / delta, 0.0, 1.0)
                end
            end
        else
            making = true
            have[_TC_POINT] = false
            fmt, width = _to_format(type)
            numtyp = type
            if !have[type]
                if type == _TC_WKDAY
                    tv = _ttrans(l, bastyp, ywfmt, _tvec(myet))
                    if timtyp == _ID_UTCP
                        tv = _ttrans(l, _TT_YMWDF, _TT_YMWDF,
                            (tv[1], tv[2], tv[3], tv[4], tv[5] + hoff, tv[6] + moff, 0.0))
                    end
                    values[_TC_WKDAY] = tv[4]
                    have[_TC_WKDAY] = true
                elseif type == _TC_SP2000 || type == _TC_SP1950
                    values[_TC_SP2000] = _ttrans(l, bastyp, _TT_FORMAL, _tvec(myet))[1]
                    values[_TC_SP1950] = values[_TC_SP2000] + sp1950
                    have[_TC_SP2000] = have[_TC_SP1950] = true
                elseif type == _TC_JD
                    values[_TC_JD] = _ttrans(l, bastyp, _TT_JDUTC, _tvec(myet))[1]
                    have[_TC_JD] = true
                else
                    tv = _ttrans(l, bastyp, ymdfmt, _tvec(myet))
                    values[_TC_SC] = tv[6]
                    if timtyp == _ID_UTCP
                        tv = _ttrans(l, _TT_YMDF, _TT_YMDF,
                            (tv[1], tv[2], tv[3], tv[4] + hoff, tv[5] + moff, 0.0, tv[7]))
                    end
                    values[_TC_HR] = tv[4]
                    values[_TC_MN] = tv[5]
                    values[_TC_AP] = tv[4] == 0.0 ? 12.0 : (tv[4] > 12.0 ? tv[4] - 12.0 : tv[4])
                    jyear, jmonth, jday, jdoy = gr2jul(_nint(tv[1]), _nint(tv[2]), _nint(tv[3]))
                    gyear, gmonth, gday, gdoy = jul2gr(jyear, jmonth, jday)
                    if caltyp == _ID_GCAL
                        go2jul = false
                    elseif caltyp == _ID_JCAL
                        go2jul = true
                    else
                        go2jul = gyear < 1582 || (gyear == 1582 &&
                                  (gmonth < 10 || (gmonth == 10 && gday < 15)))
                    end
                    cy, cm, cd, cdoy = go2jul ? (jyear, jmonth, jday, jdoy) :
                                       (gyear, gmonth, gday, gdoy)
                    values[_TC_YEAR] = Float64(cy)
                    values[_TC_MM] = Float64(cm)
                    values[_TC_DD] = Float64(cd)
                    values[_TC_DOY] = Float64(cdoy)
                    values[_TC_SIGNED_YEAR] = values[_TC_YEAR]
                    if doera
                        if values[_TC_YEAR] < 1.0
                            values[_TC_YEAR] = 1.0 - values[_TC_YEAR]
                            values[_TC_ERA] = 1.0
                        else
                            values[_TC_ERA] = 2.0
                        end
                        vanish = values[_TC_YEAR] >= 1000.0
                    end
                    values[_TC_YR] = _rmaind(values[_TC_YEAR], 100.0)[2]
                    for t in (_TC_YEAR, _TC_YR, _TC_DOY, _TC_MM, _TC_MONTH, _TC_DD, _TC_HR,
                        _TC_MN, _TC_SC, _TC_AP, _TC_ERA)
                        have[t] = true
                    end
                end
            end
            if type == _TC_MONTH || type == _TC_WKDAY
                name = type == _TC_MONTH ? _MONTH_NAMES[clamp(_nint(values[_TC_MM]), 1, 12)] :
                       _WEEKDAY_NAMES[clamp(_nint(values[_TC_WKDAY]), 1, 7)]
                style = (ident[i] - (type == _TC_MONTH ? 4 : 12))   # 0..5, as the marks are listed
                style < 3 && (name = first(name, 3))
                name = style % 3 == 0 ? uppercase(name) : (style % 3 == 2 ? lowercase(name) : name)
                put(name)
                appnd += ncodeunits(name)
                making = false
            elseif type == _TC_ERA
                upper = ident[i] == 33 || ident[i] == 35
                ad = values[_TC_ERA] == 2.0
                put(ad ? (upper ? " A.D." : " a.d.") : (upper ? " B.C." : " b.c."))
                if ident[i] >= 35
                    if vanish
                        put(" ")
                        appnd += 1
                    else
                        appnd += 6
                    end
                else
                    appnd <= 256 && put(lstrip(String(out[appnd:256]), ' '))
                    appnd += 4
                end
                making = false
            elseif type == _TC_AMPM
                pm = values[_TC_HR] >= 12.0
                put(ident[i] == 37 ? (pm ? "P.M." : "A.M.") : (pm ? "p.m." : "a.m."))
                appnd += 4
                making = false
            else
                value = values[type]
            end
            if making
                # The time vectors of the start of this unit and of the next one, for
                # the fraction that a `.#` may ask for.
                fill!(ptvec, 0.0)
                fill!(ntvec, 0.0)
                incr = 1.0
                timfmt = 1
                intyp = ymdfmt
                if type == _TC_YEAR || type == _TC_YR
                    stopat = 1
                elseif type == _TC_MM
                    stopat = 2
                elseif type == _TC_DD || type == _TC_DOY
                    stopat = 3
                elseif type == _TC_HR || type == _TC_AP
                    stopat = 4
                elseif type == _TC_MN
                    stopat = 5
                elseif type == _TC_SC
                    stopat = 6
                    incr = 0.0
                else
                    stopat = 1
                    timfmt = 2
                    incr = 0.0
                    intyp = type != _TC_JD ? bastyp :
                            ((timtyp == _ID_TDT || timtyp == _ID_TT) ? _TT_JDTDT :
                             (timtyp == _ID_TDB ? _TT_JDTDB : _TT_JDUTC))
                end
                if timfmt == 1
                    comps = (_TC_SIGNED_YEAR, _TC_MM, _TC_DD, _TC_HR, _TC_MN, _TC_SC)
                    for j in 1:stopat
                        ptvec[j] = ntvec[j] = values[comps[j]]
                    end
                else
                    ptvec[1] = ntvec[1] = values[type]
                end
                ntvec[stopat] += incr
                if type == _TC_YEAR || type == _TC_YR
                    ptvec[2] = ntvec[2] = ptvec[3] = ntvec[3] = 1.0
                elseif type == _TC_MM
                    ptvec[3] = ntvec[3] = 1.0
                end
                if go2jul && timfmt != 2
                    for v in (ptvec, ntvec)
                        gy, gm, gd, _ = jul2gr(_nint(v[1]), _nint(v[2]), _nint(v[3]))
                        v[1], v[2], v[3] = Float64(gy), Float64(gm), Float64(gd)
                    end
                end
                if dozone && timfmt != 2
                    for v in (ptvec, ntvec)
                        v[4] -= hoff
                        v[5] -= moff
                        v[6] = 0.0
                        v .= _ttrans(l, _TT_YMDF, _TT_YMDF, tup(v))
                        type == _TC_SC && (v[6] = values[_TC_SC])
                    end
                end
            end
        end
    end
    return String(rstrip(String(out), ' '))
end

# ---------------------------------------------------------------------------------
# ET2UTC, ETCAL
# ---------------------------------------------------------------------------------

"""
    et2utc(ks, et, format, prec) -> String

Contract: `format` is `"C"` (`1986 APR 12 16:31:09.814`), `"D"`
(`1986-102 // 16:31:09.814`), `"ISOC"` (`1986-04-12T16:31:09.814`), `"ISOD"`
(`1986-102T16:31:09.814`) or `"J"` (`JD 2446533.18834276`), in any case; `prec` is
the number of decimal places and is limited to 0..14 (CSPICE writes at most 14
significant digits, so 13 places of seconds and 7 of a Julian date). Returns the UTC
string of the TDB epoch `et` as CSPICE `et2utc` builds it, with seconds 60 inside a
leap second; the fraction of the seconds is rounded, and a year before 1000 has
`A.D.` or `B.C.` in the formats `C` and `D`. Throws `SpiceError(:badformat)` for
another format, `SpiceError(:badtime)` for a year before 1 A.D. in an ISO format, and
`SpiceError(:nolsk)` without a leap-second kernel. Reads `ks.lsk` only.
"""
function et2utc(ks::KernelSet, et::Real, format::AbstractString, prec::Integer)::String
    fmt = uppercase(rstrip(first(format, 4), ' '))
    fmt in ("C", "D", "J", "ISOC", "ISOD") || throw(SpiceError(:badformat,
        "the format '$format' is not one of 'C', 'D', 'J', 'ISOC', 'ISOD'"))
    l = _lsk(ks)
    p = Int(clamp(prec, 0, 14))
    x = primal(et)
    if fmt == "J"
        jd = _ttrans(l, _TT_TDB, _TT_JDUTC, _tvec(x))[1]
        return "JD" * _dpstrf(jd, p + 7)
    end
    tai = _unitim(l, x, _US_TDB, _US_TAI)
    whole = floor(tai)
    scale = _POW10[p+1]
    frac = _dnint(scale * (tai - whole))
    if frac == scale
        whole += 1.0
        frac = 0.0
    end
    calendar_style = fmt == "C" || fmt == "D"
    if fmt == "C" || fmt == "ISOC"
        tv = _ttrans(l, _TT_TAI, _TT_YMD, _tvec(whole))
        year, month, day = _nint(tv[1]), _nint(tv[2]), _nint(tv[3])
        hh, mi, sc = _nint(tv[4]), _nint(tv[5]), _nint(tv[6])
    else
        tv = _ttrans(l, _TT_TAI, _TT_YD, _tvec(whole))
        year, month, day = _nint(tv[1]), 1, _nint(tv[2])
        hh, mi, sc = _nint(tv[3]), _nint(tv[4]), _nint(tv[5])
    end
    clock = string(_pad(hh, 2), ':', _pad(mi, 2), ':', _pad(sc, 2))
    if p > 0
        # CSPICE writes `1 + fraction` with `dpstrf`, so at most 13 places exist.
        fs = _dpstrf(frac / scale + 1.0, p + 1)
        dot = findfirst('.', fs)
        clock = string(clock, SubString(fs, dot, min(dot + p, lastindex(fs))))
    end
    ystr = if year >= 1000
        string(year)
    elseif year > 0
        calendar_style ? string(year, " A.D.") : string(year)
    else
        calendar_style || throw(SpiceError(:badtime,
            "the year $(1 - year) B.C. cannot be written in an ISO format"))
        string(1 - year, " B.C.")
    end
    if fmt == "C"
        return string(ystr, ' ', _MONTH_ABBR[month], ' ', _pad(day, 2), ' ', clock)
    elseif fmt == "ISOC"
        return string(ystr, '-', _pad(month, 2), '-', _pad(day, 2), 'T', clock)
    elseif fmt == "D"
        return string(ystr, year < 1000 ? ' ' : '-', _pad(day, 3), " // ", clock)
    else
        return string(ystr, '-', _pad(day, 3), 'T', clock)
    end
end

"""
    etcal(et) -> String

Contract: `et` is TDB seconds past J2000. Returns the calendar string of `et` on a
formal calendar in which every day has 86400 seconds, in the form
`1986 APR 12 16:31:09.814` (TDB, no leap seconds, no kernel needed), with `A.D.`
after a year before 1000 and `B.C.` for a year before 1, as CSPICE `etcal` builds it;
the seconds are rounded to milliseconds but never up to 60. An epoch whose day number
is outside the 32-bit integers gives `Epoch before ...` or `Epoch after ...` and the
first or last date that can be written.
"""
function etcal(et::Real)::String
    mynum = primal(et) + _HALFD
    q = trunc(mynum / _SPD)
    secs = mynum - q * _SPD
    if secs < 0.0
        q -= 1.0
        secs += _SPD
    end
    dmxint = 2147483647.0
    dmnint = -2147483648.0
    messge = ""
    if !(q + _DN2000 >= dmnint + 1.0)
        q = dmnint - _DN2000 + 1.0
        messge = "Epoch before"
        secs = 0.0
    elseif q + _DN2000 > dmxint - 1.0
        q = dmxint - _DN2000 - 1.0
        messge = "Epoch after"
        secs = 0.0
    end
    daynum = Int(q + _DN2000)
    offset, daynum = fldmod(daynum, 146097)
    adjust = offset < 0
    year, month, day, _ = _calendar(daynum)
    year += offset * 400
    era = ""
    if adjust
        year = 1 - year
        era = "B.C."
    elseif year < 1000
        era = "A.D."
    end
    tsecs = trunc(Int, secs * 1e3)
    hours = div(tsecs, 3600000)
    tsecs -= hours * 3600000
    mins = div(tsecs, 60000)
    tsecs -= mins * 60000
    secs = Float64(tsecs) / 1e3
    secs + 5e-4 < 60.0 && (secs += 5e-4)
    secs = trunc(Int, secs * 1e3) / 1e3
    sstr = lstrip(_dpstrf(secs, 6), ' ')
    (ncodeunits(sstr) >= 2 && sstr[2] == '.') && (sstr = "0" * sstr)
    sstr = first(sstr, 6)
    parts = (messge, string(year), era, _MONTH_ABBR[month], _pad(day, 2),
        string(_pad(hours, 2), ':', _pad(mins, 2), ':', sstr))
    return join(Iterators.filter(!isempty, parts), ' ')
end
