# Group 9: the time subsystem against CSPICE.
#
# Feature file: features/time_conversions.feature
#
# Strings to ET (`str2et`, `utc2et`, `tparse`), ET to strings (`timout`, `et2utc`,
# `etcal`), pictures (`tpictr`), uniform scales (`unitim`), time vectors (`ttrans`,
# compared with the Fortran-style entry `ttrans_` of the CSPICE library, which SPICE.jl
# does not wrap), the calendars (`jul2gr`, `gr2jul`), the constants and the defaults of
# `timdef`. CSPICE keeps the `timdef` defaults in process-wide state: the parts that
# change them restore `GREGORIAN` and `UTC` before they end.
#
# A conversion has one of two outcomes: a value, or a rejection (an exception of the
# engine). Both engines must have the same outcome for every input, and equal values.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff

"Value of `f()` from CSPICE, or `nothing` if CSPICE signals an error."
function oracle_outcome(f)
    try
        return f()
    catch e
        e isa SPICE.SpiceError || rethrow()
        return nothing
    end
end

"Value of `f()` from Allspice, or `nothing` if it throws a `SpiceError`."
function neo_outcome(f)
    try
        return f()
    catch e
        e isa Allspice.SpiceError || rethrow()
        return nothing
    end
end

"""
    compare_outcomes(what, inputs, oracle, neo; atol=0.0) -> NamedTuple

Run `oracle(x)` and `neo(x)` for every `x` of `inputs` and count the outcomes: `both`
(two values), `rejected` (two rejections), `mixed` (one value and one rejection),
`unequal` (two values that are not equal: `!==` for numbers within `atol`, `!=` for
strings) and the largest absolute difference of two numeric values. The first
disagreements are printed. Postcondition: `both + rejected + mixed == length(inputs)`.
"""
function compare_outcomes(what::AbstractString, inputs, oracle, neo; atol::Float64=0.0)
    both = rejected = mixed = unequal = 0
    worst = 0.0
    shown = 0
    for x in inputs
        ref = oracle_outcome(() -> oracle(x))
        got = neo_outcome(() -> neo(x))
        if ref === nothing && got === nothing
            rejected += 1
        elseif ref === nothing || got === nothing
            mixed += 1
            shown += 1
            shown <= 8 && println("    outcome differs for ", repr(x), ": CSPICE ",
                ref === nothing ? "rejects" : repr(ref), ", Allspice ",
                got === nothing ? "rejects" : repr(got))
        else
            both += 1
            same = if ref isa Real
                worst = max(worst, abs(Float64(got) - Float64(ref)))
                ref === got || abs(Float64(got) - Float64(ref)) <= atol
            else
                ref == got
            end
            if !same
                unequal += 1
                shown += 1
                shown <= 8 && println("    value differs for ", repr(x), ": CSPICE ",
                    repr(ref), ", Allspice ", repr(got))
            end
        end
    end
    @assert both + rejected + mixed == length(inputs)
    measured(@sprintf("  %-44s n=%-6d accepted=%-6d rejected=%-6d mixed=%d unequal=%d max|d|=%.3e",
        what, length(inputs), both, rejected, mixed, unequal, worst))
    return (; both, rejected, mixed, unequal, worst)
end

# ---------------------------------------------------------------------------------
# Generator of time strings
# ---------------------------------------------------------------------------------

const TS_MON = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
const TS_MONTH = ["January", "February", "March", "April", "May", "June", "July", "August",
    "September", "October", "November", "December"]
const TS_WKD = ["Sun", "Monday", "Tue", "Wednesday", "Thu", "Friday", "Sat", "Tues", "Thursday"]
const TS_ZONE = ["PDT", "PST", "MDT", "MST", "CDT", "CST", "EDT", "EST", "UTC+5", "UTC-7:30",
    "UTC+12:59", "UTC-0:15", "utc+3", "UTC+05:30", "pst"]
const TS_SYSTEM = ["UTC", "TDB", "TDT", "utc", "tdb", "(TDB)", "(UTC)"]

"Days of month `m` of the Gregorian year `y`."
ts_dim(y, m) = m == 2 ? (((y % 4 == 0) && (y % 100 != 0)) || (y % 400 == 0) ? 29 : 28) :
               (m in (4, 6, 9, 11) ? 30 : 31)

"Random valid Gregorian date (year, month, day) with the year in `years`."
function ts_date(rng, years=1000:2999)
    y = rand(rng, years)
    m = rand(rng, 1:12)
    return y, m, rand(rng, 1:ts_dim(y, m))
end

"Day of year of a Gregorian date."
ts_doy(y, m, d) = sum(ts_dim(y, k) for k in 1:m-1; init=0) + d

"Month name in a random spelling: abbreviation, full name, prefix of the full name, any case, optional period."
function ts_monthname(rng, m)
    full = TS_MONTH[m]
    name = rand(rng, (TS_MON[m], full, full[1:rand(rng, 3:length(full))]))
    name = rand(rng, (name, uppercase(name), lowercase(name)))
    return (length(name) == 3 && rand(rng) < 0.15) ? name * "." : name
end

"Seconds `0 <= s < 60` written with 0 to 9 decimals; `\"SS\"`, `\"SS.\"` or `\"SS.fff\"`."
function ts_seconds(rng)
    whole = @sprintf("%02d", rand(rng, 0:59))
    k = rand(rng, (0, 0, 1, 3, 3, 6, 9))
    k == 0 && return rand(rng) < 0.1 ? whole * "." : whole
    return whole * "." * join(rand(rng, '0':'9', k))
end

"Clock in a random form: `HH:MM:SS.fff`, `HH:MM`, `HH:MM.f`, or empty."
function ts_clock(rng)
    h = @sprintf("%02d", rand(rng, 0:23))
    mi = @sprintf("%02d", rand(rng, 0:59))
    r = rand(rng)
    r < 0.60 && return string(h, ':', mi, ':', ts_seconds(rng))
    r < 0.80 && return string(h, ':', mi)
    r < 0.88 && return string(h, ':', mi, '.', rand(rng, 0:9))
    return ""
end

"String of the ISO 8601 family: calendar or day of year, `T` separator, optional `Z`."
function ts_iso(rng)
    y, m, d = ts_date(rng, rand(rng) < 0.9 ? (1000:2999) : (1:99))
    date = rand(rng) < 0.7 ? @sprintf("%04d-%02d-%02d", y, m, d) : @sprintf("%04d-%03d", y, ts_doy(y, m, d))
    r = rand(rng)
    time = r < 0.6 ? @sprintf("%02d:%02d:", rand(rng, 0:23), rand(rng, 0:59)) * ts_seconds(rng) :
           r < 0.75 ? @sprintf("%02d:%02d", rand(rng, 0:23), rand(rng, 0:59)) :
           r < 0.85 ? @sprintf("%02d", rand(rng, 0:23)) :
           r < 0.92 ? @sprintf("%02d.%d", rand(rng, 0:23), rand(rng, 0:9)) : ""
    return string(date, 'T', time, rand(rng) < 0.3 ? "Z" : "")
end

"String with a month name, in one of the orders that `tpartv` knows, with random separators."
function ts_named(rng)
    y, m, d = ts_date(rng)
    mon = ts_monthname(rng, m)
    sep = rand(rng, (" ", " ", ", ", "-", "/", "  "))
    date = rand(rng, (string(mon, ' ', d, sep, y), string(d, ' ', mon, ' ', y),
        string(y, ' ', mon, ' ', d), string(mon, ' ', d, ", ", y), string(d, '-', mon, '-', y),
        string(y, sep, mon, sep, d), string(y, ' ', d, ' ', mon), string(mon, sep, d, sep, y)))
    clock = ts_clock(rng)
    s = clock == "" ? date : (rand(rng) < 0.8 ? string(date, ' ', clock) : string(clock, ' ', date))
    rand(rng) < 0.2 && (s = string(rand(rng, TS_WKD), rand(rng, (" ", ", ", ". ")), s))
    return s
end

"String with a numeric month: `MM/DD/YYYY`, `MM-DD-YYYY`, `YYYY MM DD`, two-digit and quoted years."
function ts_numeric(rng)
    y, m, d = ts_date(rng, 1969:2068)
    yy = y % 100
    date = rand(rng, (string(m, '/', d, '/', y), string(m, '-', d, '-', y), string(y, ' ', m, ' ', d),
        string(m, '/', d, '/', @sprintf("%02d", yy)), string(y, '-', m, '-', d),
        string(m, ' ', d, ' ', y), string(m, '/', d, "/'", @sprintf("%02d", yy)),
        string(y, ", ", m, ", ", d), string(y, '/', m, '/', d)))
    clock = ts_clock(rng)
    return clock == "" ? date : string(date, ' ', clock)
end

"String of the day-of-year family: `YYYY-DDD // HH:MM:SS`, `YYYY DDD::HH:MM`, `DDD-YYYY/...`."
function ts_dayofyear(rng)
    y, m, d = ts_date(rng)
    doy = ts_doy(y, m, d)
    ds = rand(rng, (string(doy), @sprintf("%03d", doy)))
    mark = rand(rng, ("//", "::", " // ", " :: ", "// ", "/"))
    date = rand(rng, (string(y, '-', ds), string(y, ' ', ds), string(ds, '-', y), string(ds, ' ', y)))
    clock = ts_clock(rng)
    return string(date, mark, clock)
end

"Julian date string: `JD 2451545.123`, `2451545.5 JD`, `(JD)`, a negative date, an optional system."
function ts_julian(rng)
    jd = rand(rng) < 0.9 ? 1.0e6 + 3.0e6 * rand(rng) : -5.0e5 * rand(rng)
    num = rand(rng, (@sprintf("%.5f", jd), @sprintf("%.1f", jd), @sprintf("%.0f", jd),
        @sprintf("%.9f", jd), @sprintf("%.0f.", jd)))
    label = rand(rng, ("JD", "jd", "(JD)", "Jd"))
    s = rand(rng) < 0.7 ? string(label, ' ', num) : string(num, ' ', label)
    r = rand(rng)
    return r < 0.6 ? s : string(s, ' ', rand(rng, TS_SYSTEM))
end

"String with an era: years 1 to 3000 with `B.C.`, `A.D.`, `BC`, `AD` in any case."
function ts_era(rng)
    y = rand(rng, (rand(rng, 1:99), rand(rng, 100:999), rand(rng, 1000:3000)))
    m = rand(rng, 1:12)
    d = rand(rng, 1:28)
    era = rand(rng, ("B.C.", "A.D.", "BC", "AD", "b.c.", "a.d.", "bc", "ad", "(B.C.)"))
    mon = ts_monthname(rng, m)
    date = rand(rng, (string(y, ' ', era, ' ', mon, ' ', d), string(mon, ' ', d, ", ", y, ' ', era),
        string(d, ' ', mon, ' ', y, ' ', era), string(era, ' ', y, ' ', mon, ' ', d)))
    clock = ts_clock(rng)
    return clock == "" ? date : string(date, ' ', clock)
end

"String with A.M./P.M.: hours 1 to 12 (sometimes 0 or 13, which CSPICE rejects)."
function ts_ampm(rng)
    y, m, d = ts_date(rng)
    h = rand(rng) < 0.9 ? rand(rng, 1:12) : rand(rng, (0, 13, 12))
    ap = rand(rng, ("A.M.", "P.M.", "AM", "PM", "am", "pm", "a.m.", "p.m."))
    clock = rand(rng, (string(h, ':', @sprintf("%02d", rand(rng, 0:59)), ':', ts_seconds(rng)),
        string(h, ':', @sprintf("%02d", rand(rng, 0:59)))))
    mon = ts_monthname(rng, m)
    return rand(rng, (string(mon, ' ', d, ", ", y, ' ', clock, ' ', ap),
        string(clock, ' ', ap, ' ', mon, ' ', d, ' ', y), string(d, ' ', mon, ' ', y, ' ', clock, ap)))
end

"String with a time zone or a time system (sometimes both, which CSPICE rejects)."
function ts_zoned(rng)
    y, m, d = ts_date(rng)
    mon = ts_monthname(rng, m)
    clock = string(@sprintf("%02d", rand(rng, 0:23)), ':', @sprintf("%02d", rand(rng, 0:59)), ':',
        ts_seconds(rng))
    r = rand(rng)
    tag = r < 0.55 ? rand(rng, TS_ZONE) : (r < 0.95 ? rand(rng, TS_SYSTEM) :
                                           string(rand(rng, TS_ZONE), ' ', rand(rng, TS_SYSTEM)))
    rand(rng) < 0.1 && !startswith(tag, "(") && (tag = string('(', tag, ')'))
    return rand(rng, (string(mon, ' ', d, ", ", y, ' ', clock, ' ', tag),
        string(y, ' ', mon, ' ', d, ' ', clock, ' ', tag), string(clock, ' ', tag, ' ', mon, ' ', d, ' ', y),
        string(y, '-', m, '-', d, ' ', clock, ' ', tag)))
end

"String at second 60 of a minute: at a real leap second, at another date, and in a zone."
function ts_leap(rng, leapdays)
    y, m, d = rand(rng) < 0.75 ? rand(rng, leapdays) : ts_date(rng, 1972:2030)
    frac = rand(rng, ("", ".0", ".5", ".999", ".123456"))
    r = rand(rng)
    r < 0.5 && return string(y, ' ', TS_MON[m], ' ', d, " 23:59:60", frac)
    r < 0.65 && return @sprintf("%04d-%02d-%02dT23:59:60%s", y, m, d, frac)
    r < 0.8 && return string(TS_MON[m], ' ', d, ", ", y, " 16:59:60", frac, " PDT")
    r < 0.9 && return string(TS_MON[m], ' ', d, ", ", y, " 11:59:60", frac, " PM")
    return string(y, ' ', TS_MON[m], ' ', d, " 23:59:60", frac, rand(rng, (" TDB", " TDT", " UTC+0", " UTC-7")))
end

"String with one component out of range: month 13, day 32, hour 24, minute 60, second 61, day 0."
function ts_outofrange(rng)
    y, m, d = ts_date(rng)
    h, mi, s = rand(rng, 0:23), rand(rng, 0:59), rand(rng, 0:59)
    k = rand(rng, 1:7)
    k == 1 && (m = rand(rng, (0, 13, 20)))
    k == 2 && (d = rand(rng, (0, 32, ts_dim(y, clamp(m, 1, 12)) + 1)))
    k == 3 && (h = rand(rng, (24, 25, 99)))
    k == 4 && (mi = rand(rng, (60, 61, 99)))
    k == 5 && (s = rand(rng, (60, 61, 99)))
    k == 6 && return @sprintf("%04d-%03dT%02d:%02d:%02d", y, rand(rng, (0, 366, 367, 400)), h, mi, s)
    k == 7 && return string(y, ' ', TS_MON[clamp(m, 1, 12)], ' ', d, ".5 ", h, ':', mi, ':', s)
    return rand(rng, (@sprintf("%04d-%02d-%02dT%02d:%02d:%02d", y, m, d, h, mi, s),
        string(m, '/', d, '/', y, ' ', h, ':', mi, ':', s),
        string(y, ' ', m, ' ', d, ' ', h, ':', mi, ':', s)))
end

const TS_ALPHABET = collect(" ,-/:.'()TZJDAPMUCBjtd+0123456789")
const TS_EXTRA = [" TDB", " PM", " B.C.", " UTC+3", " Monday", " JD", " A.D.", " PDT", " ::", " //",
    " 12:00", " 1999", " Jan", ",", "-", " (", ")", " '", " T", "Z", " UTC", " 7", ".5"]

"A valid string with one random edit: a character removed, doubled, inserted or replaced, two words exchanged, or a word added."
function ts_mutate(rng, s::String)
    c = collect(s)
    k = rand(rng, 1:6)
    i = rand(rng, 1:length(c))
    if k == 1 && length(c) > 1
        deleteat!(c, i)
    elseif k == 2
        insert!(c, i, c[i])
    elseif k == 3
        insert!(c, i, rand(rng, TS_ALPHABET))
    elseif k == 4
        c[i] = rand(rng, TS_ALPHABET)
    elseif k == 5
        w = split(s, ' ')
        if length(w) > 1
            a, b = rand(rng, 1:length(w)), rand(rng, 1:length(w))
            w[a], w[b] = w[b], w[a]
        end
        return join(w, ' ')
    else
        return rand(rng) < 0.5 ? s * rand(rng, TS_EXTRA) : lstrip(rand(rng, TS_EXTRA)) * " " * s
    end
    return String(c)
end

"Strings that both engines must reject."
const TS_MALFORMED = ["", " ", "\t", "hello", "Jan", "2000", "12:00:00", "Jan 2000", "2000 Jan",
    "2000-01-01T", "2000-01-01T12:00:00 TDB", "2000-1-1T12:00:00.0.0", "T12:00:00",
    "2000-01-01T25:00:00", "2000-02-30", "1999-02-29", "2000-13-01", "2000-00-10",
    "2000 Jan 1 12:00:00 PDT TDB", "2000 Jan 1 12:00:00 PDT EST", "2000 Jan 1 TDB TDT",
    "2000 Jan 1 13:00 PM", "2000 Jan 1 0:30 AM", "Jan 1 2000 B.C. A.D.", "'99 Jan 1 B.C.",
    "2000 Jan Feb 1", "2000 2001 Jan 1", "Jan 1, 2000 12:00:00:00:00", "2000 Jan 1 12::00",
    "2000,,Jan 1", "2000 Jan 1,", ",2000 Jan 1", "2000-Jan-1-", "2000 Jan 1 24:00:00",
    "2000 Jan 1 12:60:00", "2000 Jan 1 12:00:60", "1999 Jun 30 23:59:60", "2000 Jan 1 23:59:60 TDB",
    "2000 Jan 1.5 12:00:00", "2000 Jan 0", "2000 Jan 32", "2000 Apr 31", "2000-367T00:00:00",
    "JD", "JD JD 2451545", "JD 2451545 2451546", "JD 2451545 Jan 1", "JD 2451545 B.C.",
    "2451545.0", "1.5", "2000 Jan 1 UTC+13", "2000 Jan 1 UTC+", "2000 Jan 1 UTC+ 5",
    "2000 Jan 1 UTC+123", "2000 Jan 1 1.2E3", "2000 Jan 1 12:00 1D2", "2000 Jan 1 (", "2000 Jan 1 )",
    "2000 Jan 1 [PDT]", "2000 Jan 1 12:00:00 XYZ", "2000 Jann 1", "2000 Janu?ry 1", "Mon Tue Jan 1 2000",
    "2000 Jan 1 12:00 AM PM", "2000-001-01", "1/2/3/4", "2000 Jan 1 \x01", "2000 Jan 1 é",
    "2000 Jan 1 12.5:30", "12 12 12", "1 2", "2000 // 12", "Jan 1 2000 //", "2000-01-01T12:00:00+01:00",
    "99999999999 Jan 1", "2000 Jan 99999999999", "TDB", "UTC", "A.D.", "PM", "::", "//", ".", "-",
    repeat("1 ", 70), repeat("Wednesday", 40) * " 2000 Jan 1",
    # Fixed cases from the random edits: a Julian date with a calendar component, two
    # days of the month, and a day of the month without a month.
    "Jan JD 1246975.0", "23 Nov 2070 2:45:4:6.35 P.M.", "2370-251 :: 0:4:51:24"]

"Examples of the header of `str2et.c` (the forms that it lists as accepted)."
const TS_STR2ET_EXAMPLES = ["1996-12-18T12:28:28", "1986-01-18T12", "1986-01-18T12:19",
    "1986-01-18T12:19:52.18", "1986-01-18T12:19:52.18Z", "1995-08T18:28:12", "1995-08T18:28:12Z",
    "1995-18T", "0000-01-01T", "Tue Aug  6 11:10:57  1996", "1 DEC 1997 12:28:29.192",
    "2/3/1996 17:18:12.002", "Mon Apr 22 09:32:38 PDT 1996", "Dec 18, 11:59:59.9999 A.M. 1995",
    "1992 11:18:28  3 Jul", "June 12, 1989 01:21", "1978/3/12 23:28:59.29", "17JUN1982 18:28:28",
    "13:28:28.128 1992 27 Jun", "1972 27 jun 12:29", "'93 Jan 23 12:29:47.289",
    "27 Jan 3, 19:12:28.182", "23 A.D. APR 4, 18:28:29.29", "18 B.C. Jun 3, 12:29:28.291",
    "29 Jun  30 12:29:29.298", "29 Jun '30 12:29:29.298", "1997-162::12:18:28.827",
    "162-1996/12:28:28.287", "1993-321/12:28:28.287", "1992 183// 12 18 19", "17:28:01.287 1992-272//",
    "17:28:01.282 272-1994//", "'92-271/ 12:28:30.291", "92-182/ 18:28:28.281", "'92-271// 12:28:30.291",
    "182-'92/ 12:29:29.192", "182-92/ 12:29:29.192", "182-'92// 12:28:29.182", "92-182// 12:28:29.182",
    "jd 28272.291", "2451515.2981 (JD)", "2451515.2981 JD", "1988:18:21", "1/2/3 4:5:6",
    "1997 Jan 31 12:29:20 TDB", "1997 Jan 31 12:29:20 TDT", "1997 Jan 31 12:29:20 UTC",
    "1997 Jan 31 09:29:20.123456 PST", "1997 Jan 31 12:29:20.123 UTC-07:30", "1997-JAN-31 12:29:20.123456",
    "JD 2450480.02038194 TDB", "Jan 31, 1997 12:29:20 (TDB)", "1/1/1 A.D. 00:00:00",
    "1:1:1 1 Jan 1 B.C.", "12:28:28 A.M. Jan 1, 1 B.C.", "  2000-01-01T12:00:00  ",
    "December 31, 2016 23:59:60.500", "DEC 31 2016 3:59:60 PM PST", "1996-06-30T23:59:60",
    "0099-01-01T00:00:00", "99-01-01T00:00:00", "68-12-31T23:59:59", "69-01-01T00:00:00",
    "Thursday, April 4, 1996 7:23 p.m.", "Thu. April 4, 1996", "THURS apr 4 '96", "4.5 Apr 1996",
    "1996 Apr 4.75", "1996 Apr 4 12.5", "1996 Apr 4 12:30.5", "1996-095.5//", "5/6/7",
    "1996-12-18T12:28:28.9999999999999999", "2000 Jan 1 00:00:00.00000000000000000000001"]

"Pictures of the header of `timout.c` and other pictures that use every mark."
const TS_PICTURES = ["MON DD,YYYY  HR:MN:SC.#### (TDB) ::TDB", "MON DD,YYYY  HR:MN:SC.#### (TT) ::TT",
    "MON DD,YYYY  HR:MN:SC.#### (UTC)", "YYYY-MM-DD ::JCAL", "YYYY-MM-DD ::MCAL", "YYYY-MM-DD ::GCAL",
    "YYYY-DOYTHR:MN:SC ::RND", "YYYY-DOY // HR:MN:SC.### ::UTC", "MON DD, YYYY HR:MN:SC ::UTC",
    "Month DD, YYYY", "DD-MON-YYYY HR:MN", "Wkd Mon DD HR:MN:SC PDT YYYY ::UTC-7",
    "Weekday Month DD, YYYY AP:MN:SC AMPM", "WKD wkd Wkd WEEKDAY weekday MONTH month mon",
    "YYYY ERA Mon DD HR:MN:SC", "YYYY era Mon DD", "YYYY?ERA? Mon DD", "YYYY?era? MM/DD",
    "JULIAND.##### ::RND", "JULIAND.##### ::TDB", "JULIAND.### ::TDT", "JD JULIAND.######## ::UTC ::RND",
    "SP2000.#######", "SP1950.### ::TDB", "SP2000 ::TDT ::RND", "SP1950", "YR-MM-DD", "'YR Mon DD",
    "YYYY.########", "YYYY.#### ::TDB ::RND", "MM.#####", "DD.#####", "DOY.########", "DOY.### ::TDB",
    "HR.######", "MN.######", "SC.######", "SC.############## ::TDB", "AP:MN.### ampm",
    "YYYY-MM-DDTHR:MN:SC.### ::UTC+5:30", "YYYY-MM-DDTHR:MN:SC ::UTC-12", "HR:MN ::UTC+0",
    "HR:MN:SC ::UTC+13", "HR:MN:SC ::UTC+x", "DD.### Mon YYYY ::UTC-3:30 ::RND", "YYYY Mon DD ::UTC+5 ::JCAL",
    "DOY ::JCAL", "DOY ::MCAL", "DD/MM/YYYY ERA ::JCAL ::TDB", "YYYY.### ::JCAL", "MM.### ::MCAL ::RND",
    "Wkd ::UTC+9", "Weekday ::TDB", "wkd, DD ::JCAL", "", " ", "no marks at all", "#", ".#", "###.###",
    "YYYYMMDD", "YYYY#", "SC.", "SC..#", "SC.#.#", "HR:MN:SC.### xyz .## #", "::RND", "::TDB",
    "YYYY-MM-DD ::TRNC", "YYYY ::RND ::TRNC", "Mon DD, YYYY AP:MN:SC ampm ::RND",
    "DD Month YYYY, AP ampm ::UTC-8", "YYYY MM DD HR MN SC ::TT ::RND", "ERA era ?ERA? ?era? YYYY"]

"A random picture: marks and literal text in a random order, with optional system, calendar and rounding marks."
function ts_picture(rng)
    marks = ["YYYY", "YR", "MON", "Mon", "mon", "MONTH", "Month", "month", "MM", "DOY", "WKD", "Wkd",
        "wkd", "WEEKDAY", "Weekday", "weekday", "DD", "MN", "HR", "SC", "JULIAND", "SP2000", "SP1950",
        "ERA", "era", "?ERA?", "?era?", "AMPM", "ampm", "AP"]
    seps = [" ", "-", ":", "/", ", ", "T", " // ", "", " at "]
    io = IOBuffer()
    for _ in 1:rand(rng, 1:6)
        print(io, rand(rng, marks))
        rand(rng) < 0.35 && print(io, ".#", "#"^rand(rng, 0:5))
        print(io, rand(rng, seps))
    end
    rand(rng) < 0.5 && print(io, ' ', rand(rng, ("::UTC", "::TDB", "::TDT", "::TT", "::UTC+5:30",
        "::UTC-8", "::UTC+0:45", "::UTC-11:59")))
    rand(rng) < 0.4 && print(io, ' ', rand(rng, ("::GCAL", "::JCAL", "::MCAL")))
    rand(rng) < 0.4 && print(io, ' ', rand(rng, ("::RND", "::TRNC")))
    return String(take!(io))
end

"""
    oracle_et2utc(et, format, prec) -> String

`et2utc_c` with room for 200 characters. The wrapper of SPICE.jl gives the string a
room that depends on the precision, and CSPICE cuts a longer string (a year before
1000 A.D. with an era label) without an error.
"""
function oracle_et2utc(et::Float64, format::String, prec::Int)
    buf = zeros(UInt8, 200)
    # A wrapper of SPICE.jl that does not look at the error state (`tpictr`) can leave
    # an error of an earlier call in it.
    ccall((:reset_c, SPICE.libcspice), Cvoid, ())
    ccall((:et2utc_c, SPICE.libcspice), Cvoid, (Cdouble, Cstring, Cint, Cint, Ptr{UInt8}),
        et, format, prec, length(buf), buf)
    SPICE.handleerror()
    return GC.@preserve buf unsafe_string(pointer(buf))
end

"""
    oracle_unchecked(f, s)

`f(s)` of CSPICE with the range check of `tparch` off, which is the start state of
CSPICE. `str2et` of CSPICE turns the check on and does not turn it off again when it
converts a second 60 in a time zone or on the Julian calendar, so that `utc2et` and
`tparse` would depend on the strings converted before.
"""
function oracle_unchecked(f, s)
    ccall((:tparch_c, SPICE.libcspice), Cvoid, (Cstring,), "NO")
    return f(s)
end

"""
    oracle_tparse(s) -> Float64

`tparse_c` of CSPICE with buffers that are set to zero before the call; throws
`SPICE.SpiceError` if CSPICE signals an error or writes an error message. The wrapper
of SPICE.jl reads the message buffer without setting it first, and CSPICE does not
write the buffer when it rejects the input through its error system (the empty
string), so the answer of the wrapper for such an input depends on the content of the
memory: one run in about ten gave the value 0.0 for the empty string.
"""
function oracle_tparse(s::AbstractString)
    msg = zeros(UInt8, 128)
    out = Ref{Cdouble}(0.0)
    ccall((:reset_c, SPICE.libcspice), Cvoid, ())
    ccall((:tparse_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Ref{Cdouble}, Ptr{UInt8}),
        s, length(msg), out, msg)
    SPICE.handleerror()
    text = GC.@preserve msg unsafe_string(pointer(msg))
    isempty(text) || throw(SPICE.SpiceError(text))
    return out[]
end

"""
    oracle_tpictr(s) -> String

`tpictr_c` of CSPICE with buffers and the success flag set to zero before the call, for
the same reason as `oracle_tparse`; throws `SPICE.SpiceError` if CSPICE signals an error
or does not set the flag.
"""
function oracle_tpictr(s::AbstractString)
    pic = zeros(UInt8, 80)
    msg = zeros(UInt8, 128)
    ok = Ref{Cint}(0)
    ccall((:reset_c, SPICE.libcspice), Cvoid, ())
    ccall((:tpictr_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Cint, Ptr{UInt8}, Ref{Cint}, Ptr{UInt8}),
        s, length(pic), length(msg), pic, ok, msg)
    SPICE.handleerror()
    ok[] != 0 || throw(SPICE.SpiceError(GC.@preserve msg unsafe_string(pointer(msg))))
    return GC.@preserve pic unsafe_string(pointer(pic))
end

"""
    oracle_stale(s) -> Bool

True if `str2et` of CSPICE gives for `s` a value that depends on the conversion
before it. This occurs for a Julian date whose sign and number are apart (`- jd 1.5`):
the number is not read and the value of the earlier call stays in the time vector.
Allspice rejects such a string.
"""
function oracle_stale(s::AbstractString)
    occursin(r"jd"i, s) || return false
    a = oracle_outcome(() -> (SPICE.str2et("JD 2451545.0"); SPICE.str2et(s)))
    b = oracle_outcome(() -> (SPICE.str2et("JD 2460000.0"); SPICE.str2et(s)))
    return a !== b
end

"`ttrans_` of the CSPICE library on a copy of `tvec` (seven components), or `nothing` on an error."
function oracle_ttrans(from::String, to::String, tvec)
    v = zeros(Float64, 10)
    v[1:length(tvec)] .= tvec
    ccall((:reset_c, SPICE.libcspice), Cvoid, ())
    ccall((:ttrans_, SPICE.libcspice), Cint, (Cstring, Cstring, Ptr{Float64}, Clong, Clong),
        from, to, v, sizeof(from), sizeof(to))
    return oracle_outcome(() -> (SPICE.handleerror(); Tuple(v[1:7])))
end

"`jul2gr_` or `gr2jul_` of the CSPICE library: (year, month, day, doy)."
function oracle_calendar(entry::Symbol, y, m, d)
    yr, mo, dy, doy = Ref{Cint}(y), Ref{Cint}(m), Ref{Cint}(d), Ref{Cint}(0)
    if entry === :jul2gr
        ccall((:jul2gr_, SPICE.libcspice), Cint, (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}), yr, mo, dy, doy)
    else
        ccall((:gr2jul_, SPICE.libcspice), Cint, (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cint}), yr, mo, dy, doy)
    end
    return Int(yr[]), Int(mo[]), Int(dy[]), Int(doy[])
end

"Set the three `timdef` defaults of CSPICE and of `ctx` back to the start values."
function reset_timdef(ctx)
    SPICE.timdef("SET", "CALENDAR", "GREGORIAN")
    SPICE.timdef("SET", "SYSTEM", "UTC")
    timdef(ctx, "SET", "CALENDAR", "GREGORIAN")
    timdef(ctx, "SET", "SYSTEM", "UTC")
    @assert timdef(ctx, "GET", "ZONE") == "" && timdef(ctx, "GET", "SYSTEM") == "UTC"
    return nothing
end

"Random time vector of the representation `from` of `ttrans`; one in five has components outside the usual ranges."
function ts_timevector(rng, from::String)
    wild = rand(rng) < 0.2
    y = Float64(rand(rng, wild ? (-3000:6000) : (1960:2040)))
    h = Float64(wild ? rand(rng, -30:60) : rand(rng, 0:23))
    mi = Float64(wild ? rand(rng, -90:150) : rand(rng, 0:59))
    s = wild ? 200000.0 * (rand(rng) - 0.5) : 60.0 * rand(rng)
    rand(rng) < 0.2 && (s = Float64(rand(rng, (0, 59, 60, 61, 86400))))
    secs = 4.0e10 * (rand(rng) - 0.5) * (wild ? 5.0 : 1.0)
    from in ("YMD", "YMDF") && return (y, Float64(wild ? rand(rng, -20:40) : rand(rng, 1:12)),
        Float64(wild ? rand(rng, -40:400) : rand(rng, 1:28)), h, mi, s, 0.0)
    from in ("YD", "YDF") && return (y, Float64(wild ? rand(rng, -400:800) : rand(rng, 1:365)), h, mi, s, 0.0, 0.0)
    from in ("YD.D", "YD.DF") && return (y, 1.0 + 365.0 * rand(rng), 0.0, 0.0, 0.0, 0.0, 0.0)
    from in ("YWD", "YWDF") && return (y, Float64(rand(rng, 1:53)), Float64(rand(rng, 1:7)), h, mi, s, 0.0)
    from in ("YMWD", "YMWDF") && return (y, Float64(rand(rng, 1:12)), Float64(rand(rng, 1:5)),
        Float64(rand(rng, 1:7)), h, mi, s)
    from == "DAYSEC" && return (Float64(rand(rng, 700000:760000)), s, 0.0, 0.0, 0.0, 0.0, 0.0)
    from == "DP2000" && return (Float64(rand(rng, -20000:20000)), s, 0.0, 0.0, 0.0, 0.0, 0.0)
    from in ("JDUTC", "JDTDB", "JDTDT", "JED") && return (2451545.0 + secs / 86400.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    return (secs, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

group("9. time subsystem against CSPICE", Kernels.LSK) do
    rng = MersenneTwister(SEED + 9)
    oracle_load([Kernels.LSK])
    ctx = neo_context([Kernels.LSK])
    ks = Allspice.kernelset(ctx)
    reset_timdef(ctx)

    # Leap-second epochs of the kernel: UTC seconds past J2000 of the day start after
    # each leap second, the calendar day that ends with it, and the ET of that instant.
    at = SPICE.gdpool("DELTET/DELTA_AT"; room=400)[2:2:end]
    leap_et = [SPICE.str2et(SPICE.et2utc(t + 100.0, "ISOC", 0)[1:10] * "T00:00:00") for t in at]
    leapdays = map(leap_et[2:end]) do t
        day = SPICE.et2utc(t - 3600.0, "ISOC", 0)
        (parse(Int, day[1:4]), parse(Int, day[6:7]), parse(Int, day[9:10]))
    end
    @test length(leapdays) >= 27

    n = nsample(150, 400)     # strings of each family (the timdef checks take 150 of a family)
    families = ["ISO" => () -> ts_iso(rng), "month name" => () -> ts_named(rng),
        "numeric month" => () -> ts_numeric(rng), "day of year" => () -> ts_dayofyear(rng),
        "Julian date" => () -> ts_julian(rng), "era" => () -> ts_era(rng), "A.M./P.M." => () -> ts_ampm(rng),
        "zone or system" => () -> ts_zoned(rng), "second 60" => () -> ts_leap(rng, leapdays),
        "out of range" => () -> ts_outofrange(rng)]
    corpus = Dict{String,Vector{String}}(name => [gen() for _ in 1:n] for (name, gen) in families)
    valid = reduce(vcat, [corpus[name] for (name, _) in families])
    corpus["one random edit"] = [ts_mutate(rng, rand(rng, valid)) for _ in 1:10n]
    corpus["two random edits"] = [ts_mutate(rng, ts_mutate(rng, rand(rng, valid))) for _ in 1:5n]
    # Strings for which CSPICE has no defined value are taken out and tested apart.
    stale = String["- jd 2451545.0", "-jd 2451545.0", "- JD 1.5 TDB"]
    for name in ("one random edit", "two random edits")
        append!(stale, filter(oracle_stale, corpus[name]))
        filter!(!oracle_stale, corpus[name])
    end
    order = [first.(families)..., "one random edit", "two random edits"]

    @apart @testset "str2et, every string family" begin
        accepted = 0
        for name in order
            r = compare_outcomes("str2et, " * name, corpus[name], SPICE.str2et, s -> str2et(ctx, s))
            @test r.mixed == 0
            @test r.unequal == 0
            accepted += r.both
            name in ("out of range", "one random edit", "two random edits") || @test r.both > n ÷ 4
        end
        @test corpus["out of range"] != [] && accepted > 15n ÷ 2
        println("  strings with a value in CSPICE that depends on the call before: ", length(stale))
        @test all(oracle_stale, stale)
        @test all(s -> throws_kind(() -> str2et(ctx, s), :badtime), stale)
        @test all(s -> throws_kind(() -> utc2et(ctx, s), :badtime), stale)
        @test all(s -> throws_kind(() -> tparse(s), :badtime), filter(s -> !occursin("TDB", s), stale))
    end

    @apart @testset "str2et, header examples and malformed strings" begin
        r = compare_outcomes("str2et, examples of str2et.c", TS_STR2ET_EXAMPLES, SPICE.str2et,
            s -> str2et(ctx, s))
        @test r.mixed == 0 && r.unequal == 0
        @test r.both >= length(TS_STR2ET_EXAMPLES) - 12
        r = compare_outcomes("str2et, malformed strings", TS_MALFORMED, SPICE.str2et, s -> str2et(ctx, s))
        @test r.mixed == 0 && r.unequal == 0
        @test r.rejected >= length(TS_MALFORMED) - 6
        # A Julian date whose day number is outside the 32-bit integers: the day
        # arithmetic of CSPICE wraps, and the result has the sign of the wrapped number.
        r = compare_outcomes("str2et, Julian date above the integer range",
            ["238409829538 Jd", "238409829538.5 JD", "JD 5373484.5"], SPICE.str2et, s -> str2et(ctx, s))
        @test r.both == 3 && r.unequal == 0
        # The message says what is wrong with the string.
        msg(s) = try str2et(ctx, s); "" catch e; e.msg end
        @test occursin("<", msg("2000 Jan 1 12:00:00 XYZ")) && occursin("unrecognizable", msg("2000 Jann 1"))
        @test occursin("month", lowercase(msg("2000-13-01")))
        @test occursin("seconds", lowercase(msg("2000 Jan 1 23:59:60")))
        @test occursin("blank", msg("  "))
        @test throws_kind(() -> str2et(Allspice.empty_kernelset(), "2000 Jan 1"), :nolsk)
        @test throws_kind(() -> str2et(ctx, "not a time"), :badtime)
    end

    @apart @testset "utc2et, tparse, tpictr" begin
        all_strings = reduce(vcat, [corpus[name] for name in order])
        sample = vcat(all_strings[1:3:end], TS_STR2ET_EXAMPLES, TS_MALFORMED)
        r = compare_outcomes("utc2et", sample, s -> oracle_unchecked(SPICE.utc2et, s), s -> utc2et(ctx, s))
        @test r.mixed == 0 && r.unequal == 0 && r.both > nsample(200, 500)
        r = compare_outcomes("tparse", sample, s -> oracle_unchecked(oracle_tparse, s), tparse)
        @test r.mixed == 0 && r.unequal == 0 && r.both > nsample(200, 500)
        r = compare_outcomes("tpictr", sample, oracle_tpictr, tpictr)
        @test r.mixed == 0 && r.unequal == 0 && r.both > nsample(200, 500)
        @test tpictr(ctx, "Mon Apr 22 09:32:38 PDT 1996") == "Wkd Mon DD HR:MN:SC PDT YYYY ::UTC-7"
        @test tparse(ctx, "2000-01-01T12:00:00") == 0.0
        # A picture made from a string writes the epoch of that string in the same form.
        for s in ["1996 Apr 04 12:28:28.250", "Thu Apr 04 07:23:11 1996", "04-APR-1996 23:59", "1996-095 // 01:02:03"]
            @test timout(ctx, str2et(ctx, s), tpictr(s)) == s
        end
    end

    # Epochs: several thousand years, and a dense set around every leap second.
    wide = random_epochs(rng, nsample(300, 3000), -5000.0 * jyear(), 4500.0 * jyear())
    modern = random_epochs(rng, nsample(400, 2000), -1.3e9, 1.3e9)
    dense = Float64[]
    for t in leap_et, d in (-2.0, -1.5, -1.0, -0.999999, -0.5, -0.0005, -1e-6, 0.0, 1e-6, 0.0005, 0.5, 1.0, 2.0)
        push!(dense, t + d, t + d + 0.3 * (rand(rng) - 0.5))
    end
    epochs = vcat(wide, modern, dense, [0.0, -0.0, 1e-9, -43200.0, 43200.0, 0.5, -0.5])
    println("  epochs: ", length(wide), " over 9500 years, ", length(modern), " near the present, ",
        length(dense), " around the ", length(leap_et), " leap-second epochs")

    @apart @testset "et2utc, all formats and precisions" begin
        cases = [(et, fmt, prec) for et in epochs for (fmt, prec) in
                 (("C", rand(rng, 0:14)), ("D", rand(rng, 0:6)), ("J", rand(rng, 0:14)),
            ("ISOC", rand(rng, 0:14)), ("ISOD", rand(rng, 0:3)), (rand(rng, ("c", "isoc", "IsoD", "j")), 3))]
        r = compare_outcomes("et2utc", cases, c -> oracle_et2utc(c...), c -> et2utc(ctx, c...))
        @test r.mixed == 0 && r.unequal == 0 && r.both > nsample(5000, 20000)
        @test et2utc(ctx, leap_et[end] - 0.5, "C", 3) == "2016 DEC 31 23:59:60.500"
        @test et2utc(ctx, leap_et[end] - 0.5, "ISOD", 1) == "2016-366T23:59:60.5"
        for bad in ("", "X", "ISO", "CAL", " C")
            @test oracle_outcome(() -> SPICE.et2utc(0.0, bad, 3)) === nothing
            @test throws_kind(() -> et2utc(ctx, 0.0, bad, 3), :badformat)
        end
        @test throws_kind(() -> et2utc(ctx, -7.0e10, "ISOC", 0), :badtime)
        # Round trip through a string with nine decimals.
        @test all(abs(str2et(ctx, et2utc(ctx, et, "C", 9)) - et) < 2e-6 for et in modern[1:200])
    end

    @apart @testset "etcal" begin
        extreme = [1e12, -1e12, 1e15, -1e15, 1.8e14, -1.9e14, 1e17, -1e17, 59.9996, 86399.9996 - 43200.0]
        r = compare_outcomes("etcal", vcat(epochs, extreme), SPICE.etcal, etcal)
        @test r.mixed == 0 && r.unequal == 0
        @test etcal(ctx, 0.0) == "2000 JAN 01 12:00:00.000"
    end

    @apart @testset "timout, pictures of the header and random pictures" begin
        sub = vcat(epochs[1:40:end], dense[1:7:end])
        cases = [(et, pic) for pic in TS_PICTURES for et in sub]
        r = compare_outcomes("timout, listed pictures", cases, c -> SPICE.timout(c[1], c[2], 256),
            c -> timout(ctx, c...))
        @test r.mixed == 0 && r.unequal == 0 && r.rejected == length(sub)   # the picture ""
        cases = [(rand(rng, epochs), ts_picture(rng)) for _ in 1:nsample(600, 6000)]
        r = compare_outcomes("timout, random pictures", cases, c -> SPICE.timout(c[1], c[2], 256),
            c -> timout(ctx, c...))
        @test r.mixed == 0 && r.unequal == 0 && r.both == length(cases)
        # Examples of the header of `timout.c`, at the epoch that it uses.
        et = SPICE.str2et("24 Mar 2018  16:23:00 UTC")
        @test timout(ctx, et, "MON DD,YYYY  HR:MN:SC.#### (TDB) ::TDB") == "MAR 24,2018  16:24:09.1856 (TDB)"
        @test timout(ctx, et, "MON DD,YYYY  HR:MN:SC.#### (UTC)") == "MAR 24,2018  16:23:00.0000 (UTC)"
        @test timout(ctx, et, "Wkd Mon DD HR:MN:SC PDT YYYY ::UTC-7") == "Sat Mar 24 09:23:00 PDT 2018"
    end

    @apart @testset "unitim: every pair of scales, generic, no allocation" begin
        systems = ["TAI", "TDB", "TDT", "ET", "JDTDB", "JDTDT", "JED", "TT", "GPS"]
        worst = 0.0; nbit = 0; ncmp = 0
        for a in systems, b in systems, _ in 1:60
            x = startswith(a, "J") ? 2451545.0 + 4.0e5 * (rand(rng) - 0.5) : 3.0e10 * (rand(rng) - 0.5)
            ref = SPICE.unitim(x, a, b)
            got = unitim(ctx, x, lowercase(a), b)
            ncmp += 1
            nbit += ref === got
            worst = max(worst, abs(ref - got))
        end
        measured(@sprintf("  %-44s n=%-6d bit-equal=%-6d max|d|=%.3e", "unitim, 81 pairs of scales", ncmp, nbit, worst))
        @test nbit == ncmp
        @test throws_kind(() -> unitim(ctx, 0.0, "UTC", "TDB"), :badtimetype)
        @test throws_kind(() -> unitim(ctx, 0.0, "TDB", "XYZ"), :badtimetype)
        @test throws_kind(() -> unitim(Allspice.empty_kernelset(), 0.0, "TDB", "TAI"), :nolsk)
        @test oracle_outcome(() -> SPICE.unitim(0.0, "UTC", "TDB")) === nothing
        # No allocation, for `unitim` and for `deltet`.
        f1(k, t) = unitim(k, t, "TDB", "TAI")
        f2(k, t) = unitim(k, t, "TDT", "JED")
        f3(k, t) = deltet(k, t, "ET")
        f1(ks, 1.0); f2(ks, 1.0); f3(ks, 1.0)
        @test (@allocated f1(ks, 1.5e8)) == 0
        @test (@allocated f2(ks, 1.5e8)) == 0
        @test (@allocated f3(ks, 1.5e8)) == 0
        # Generic over the number type: a dual number and a Taylor series keep the value
        # and carry the derivative, which is 1 plus a periodic term of order 1e-9.
        for t in (0.0, 4.2e8, -9.0e8)
            d = ForwardDiff.derivative(x -> unitim(ks, x, "TDB", "TDT"), t)
            fd = (unitim(ks, t + 5.0e4, "TDB", "TDT") - unitim(ks, t - 5.0e4, "TDB", "TDT")) / 1.0e5
            @test abs(d - 1.0) < 1e-8 && d != 1.0
            @test abs(d - fd) < 5e-12
            @test ForwardDiff.derivative(x -> unitim(ks, x, "TAI", "JDTDT"), t) == 1.0 / 86400.0
            tay = unitim(ks, Allspice.taylor_var(t, Val(2)), "TDT", "TDB")
            @test Allspice.primal(tay) == unitim(ks, t, "TDT", "TDB")
            @test abs(Allspice.derivatives(tay)[2] - ForwardDiff.derivative(x -> unitim(ks, x, "TDT", "TDB"), t)) < 1e-14
            @test ForwardDiff.value(unitim(ks, ForwardDiff.Dual(t, 1.0), "TDB", "TAI")) === unitim(ks, t, "TDB", "TAI")
        end
    end

    @apart @testset "deltet agrees with unitim and ttrans" begin
        worst = 0.0; wref = 0.0
        # ET - UTC from the time vectors against `deltet` in both directions, away from
        # the leap seconds (the formal UTC count repeats one second there). `deltet`
        # evaluates the periodic term at a UTC epoch with the UTC value, so the two
        # differ by the change of that term in about one minute (below 1e-7 s).
        for et in modern[1:400]
            utc = Allspice.ttrans(ctx, "TDB", "FORMAL", (et,))[1]
            worst = max(worst, abs((et - utc) - deltet(ctx, et, "ET")), abs((et - utc) - deltet(ctx, utc, "UTC")))
            tai = unitim(ctx, et, "TDB", "TAI")
            leaps = deltet(ctx, et, "ET") - (et - tai)      # TAI - UTC, a whole number
            @test abs(leaps - round(leaps)) < 1e-6 && 10.0 <= round(leaps) + 1.0 <= 38.0
        end
        for et in vcat(modern, dense), sys in ("ET", "UTC", "et", "utc")
            wref = max(wref, abs(deltet(ctx, et, sys) - SPICE.deltet(et, uppercase(sys))))
        end
        measured(@sprintf("  %-44s max=%.3e", "deltet against ttrans, absolute s", worst))
        measured(@sprintf("  %-44s max=%.3e", "deltet against CSPICE, absolute s", wref))
        @test wref == 0.0
        @test worst < 1e-6
    end

    @apart @testset "ttrans against the CSPICE library" begin
        names = collect(Allspice._TT_NAMES)
        ncmp = 0; nbit = 0; nmixed = 0; worst = 0.0; shown = 0
        for from in names, to in names, _ in 1:nsample(3, 25)
            tv = ts_timevector(rng, from)
            ref = oracle_ttrans(from, to, tv)
            got = neo_outcome(() -> Allspice.ttrans(ctx, from, to, tv))
            if ref === nothing || got === nothing
                nmixed += (ref === nothing) != (got === nothing)
                continue
            end
            ncmp += 1
            same = all(ref .=== got)
            nbit += same
            worst = max(worst, maximum(abs.(ref .- got)))
            if !same
                shown += 1
                shown <= 6 && println("    ttrans ", from, " -> ", to, " of ", tv, ": CSPICE ", ref, ", Allspice ", got)
            end
        end
        measured(@sprintf("  %-44s n=%-6d bit-equal=%-6d mixed=%d max|d|=%.3e", "ttrans, 441 pairs of representations",
            ncmp, nbit, nmixed, worst))
        @test nbit == ncmp && nmixed == 0 && ncmp > nsample(1000, 10000)
        @test throws_kind(() -> Allspice.ttrans(ctx, "YMD", "XYZ", (2000.0, 1.0, 1.0)), :badtimetype)
        tv0 = Allspice._tvec(2016.0, 12.0, 31.0, 23.0, 59.0, 60.5)
        l = ks.lsk
        g(a, b) = Allspice._ttrans(a, Allspice._TT_YMD, Allspice._TT_TDB, b)
        g(l, tv0)
        @test (@allocated g(l, tv0)) == 0
    end

    @apart @testset "jul2gr, gr2jul and the constants" begin
        bad = 0
        ndate = nsample(400, 4000)
        for _ in 1:ndate
            y, m, d = rand(rng, -5000:5000), rand(rng, -5:20), rand(rng, 1:28)
            bad += jul2gr(y, m, d) != oracle_calendar(:jul2gr, y, m, d)
            bad += gr2jul(y, m, d) != oracle_calendar(:gr2jul, y, m, d)
            g = jul2gr(y, m, d)
            bad += gr2jul(g[1], g[2], g[3])[1:3] != jul2gr(gr2jul(g[1], g[2], g[3])[1:3]...)[1:3] &&
                   jul2gr(gr2jul(g[1], g[2], g[3])[1:3]...)[1:3] != g[1:3]
        end
        println("  jul2gr and gr2jul: ", 2ndate, " dates compared, different: ", bad)
        @test bad == 0
        @test jul2gr(1582, 10, 5) == (1582, 10, 15, 288)
        @test gr2jul(1582, 10, 15) == (1582, 10, 5, 278)
        for (mine, ref) in ((j1900, SPICE.j1900), (j1950, SPICE.j1950), (j2000, SPICE.j2000),
            (j2100, SPICE.j2100), (jyear, SPICE.jyear), (tyear, SPICE.tyear), (spd, SPICE.spd),
            (b1900, SPICE.b1900), (b1950, SPICE.b1950), (clight, SPICE.clight))
            @test mine() === ref()
        end
        @test halfpi() === SPICE.halfpi() && twopi() === SPICE.twopi()
        @test rpd() === SPICE.rpd() && dpr() === SPICE.dpr()
    end

    @apart @testset "timdef: defaults in the context" begin
        try
            @test timdef(ctx, "GET", "CALENDAR") == SPICE.timdef("GET", "CALENDAR") == "GREGORIAN"
            @test timdef(ctx, "get", " system") == SPICE.timdef("GET", "SYSTEM") == "UTC"
            @test timdef(ctx, "GET", "ZONE") == SPICE.timdef("GET", "ZONE") == ""
            strings = vcat(corpus["month name"][1:150], corpus["second 60"][1:150], corpus["era"][1:100],
                corpus["zone or system"][1:100], corpus["Julian date"][1:60], corpus["day of year"][1:60],
                ["1582 Oct 4", "1582 Oct 5", "1582 Oct 10", "1582 Oct 14", "1582 Oct 15", "1582-278//",
                    "1582-279//", "1500 Feb 29", "1900 Feb 29", "1579 Feb 29", "1581 Feb 29", "4 B.C. Feb 29",
                    "5 B.C. Feb 29", "1000 Feb 29 12:00", "2016 Dec 31 15:59:60 ", "2016 Dec 31 23:59:60"])
            pictures = ["YYYY Mon DD HR:MN:SC.###", "YYYY-DOYTHR:MN:SC ::RND", "Wkd DD/MM/YYYY ERA AP:MN ampm",
                "JULIAND.#####", "SP2000.###", "YYYY MM DD ::GCAL ::UTC", "YYYY.#### DOY.##"]
            cases = [(et, pic) for pic in pictures for et in vcat(wide[1:60], dense[1:9:end])]
            settings = [("SYSTEM", "TDB"), ("SYSTEM", "tdt"), ("SYSTEM", "TT"), ("ZONE", "PST"),
                ("ZONE", "UTC+5:30"), ("ZONE", "utc-11:59"), ("ZONE", "EDT"), ("SYSTEM", "UTC"),
                ("CALENDAR", "JULIAN"), ("ZONE", "UTC-3"), ("CALENDAR", "MIXED"), ("SYSTEM", "TDB"),
                ("SYSTEM", "UTC"), ("CALENDAR", "gregorian"), ("ZONE", "UTC+0")]
            for (item, value) in settings
                SPICE.timdef("SET", item, value)
                timdef(ctx, "SET", item, value)
                for it in ("CALENDAR", "SYSTEM", "ZONE")
                    @test timdef(ctx, "GET", it) == SPICE.timdef("GET", it)
                end
                label = string(item, '=', value, ", ", timdef(ctx, "GET", "CALENDAR"))
                r = compare_outcomes("str2et with " * label, strings, SPICE.str2et, s -> str2et(ctx, s))
                @test r.mixed == 0 && r.unequal == 0 && r.both > 100
                r = compare_outcomes("timout with " * label, cases, c -> SPICE.timout(c[1], c[2], 256),
                    c -> timout(ctx, c...))
                @test r.mixed == 0 && r.unequal == 0 && r.both == length(cases)
            end
            for (item, value) in (("SYSTEM", "TAI"), ("ZONE", "XYZ"), ("ZONE", "UTC+13"), ("ZONE", "UTC"),
                ("CALENDAR", "ROMAN"), ("EPOCH", "J2000"), ("SYSTEM", ""))
                before = (timdef(ctx, "GET", "CALENDAR"), timdef(ctx, "GET", "SYSTEM"), timdef(ctx, "GET", "ZONE"))
                @test oracle_outcome(() -> SPICE.timdef("SET", item, value)) === nothing
                @test throws_kind(() -> timdef(ctx, "SET", item, value), :badtimedefault)
                @test before == (timdef(ctx, "GET", "CALENDAR"), timdef(ctx, "GET", "SYSTEM"), timdef(ctx, "GET", "ZONE"))
            end
            @test throws_kind(() -> timdef(ctx, "PUT", "ZONE", "PST"), :badtimedefault)
            # The defaults belong to one context: another context keeps its own.
            other = neo_context([Kernels.LSK])
            timdef(ctx, "SET", "SYSTEM", "TDB")
            @test timdef(other, "GET", "SYSTEM") == "UTC"
            @test str2et(ctx, "2000 Jan 1 12:00:00") == 0.0
            @test str2et(other, "2000 Jan 1 12:00:00") == SPICE.str2et("2000 Jan 1 12:00:00 UTC")
            # The form on a `KernelSet` takes the defaults as an argument and has no state.
            @test str2et(ks, "2000 Jan 1 12:00:00") == str2et(other, "2000 Jan 1 12:00:00")
            @test str2et(ks, "2000 Jan 1 12:00:00", Allspice.TimeDefaults("GREGORIAN", "TDB", "")) == 0.0
            # `tsetyr`: the interval of the years of two digits, for one context only.
            yctx = neo_context([Kernels.LSK])
            try
                for bound in (1969, 1980, 2000, 1950, 2101, 1)
                    SPICE.tsetyr(bound)
                    tsetyr(yctx, bound)
                    for yy in (0, 1, 49, 50, 68, 69, 79, 80, 99)
                        str = "$(lpad(yy, 2, '0')) Mar 4 05:06:07"
                        @test str2et(yctx, str) == SPICE.str2et(str)
                        @test utc2et(yctx, str) == SPICE.utc2et(str)
                    end
                    @test str2et(yctx, "105 A.D. Jan 1") == SPICE.str2et("105 A.D. Jan 1")    # not abbreviated
                end
                # The `@` dates of a text kernel loaded afterwards take the same interval,
                # and the reader accepts and rejects the strings that `tparse_` does.
                kdates = ["60-JAN-1", "JAN-1-60", "99-JAN-1", "JAN-1-05", "1-JAN-49",
                    "15-MAR-87", "31-DEC-50/12:00", "1949-01-02", "01-02-1949", "1-2-49",
                    "1-2-999", "49-1-2", "1-2-0049", "13-2-1949", "1949-13-2", "1949-2-30",
                    "1-2-49 12:00", "JD2451545.5", "2000-JAN-1T12:00:00"]
                kdate(f) = try f() catch; :rejected end
                for bound in (1969, 1950, 2000), s in kdates
                    SPICE.tsetyr(bound)
                    tsetyr(yctx, bound)
                    line = ["YB_DATE = @$s"]
                    want = kdate(() -> (oracle_unchecked(SPICE.lmpool, line); SPICE.gdpool("YB_DATE")))
                    @test kdate(() -> (lmpool(yctx, line); gdpool(yctx, "YB_DATE"))) == want
                end
                mktempdir() do dir
                    tk = joinpath(dir, "yb.tk")
                    write(tk, "\\begindata\nYB_FILE = ( @60-JAN-1, @7-JUL-68 )\n\\begintext\n")
                    SPICE.tsetyr(1950); tsetyr(yctx, 1950)
                    oracle_unchecked(SPICE.furnsh, tk); furnsh(yctx, tk)
                    want = SPICE.gdpool("YB_FILE")
                    SPICE.unload(tk)
                    @test gdpool(yctx, "YB_FILE") == want
                    @test want[1] == tparse("1960 Jan 1")
                    unload(yctx, tk)
                    ldpool(yctx, tk)
                    @test gdpool(yctx, "YB_FILE") == want
                end
                tsetyr(yctx, 1980)
                @test tparse(yctx, "79 Jan 1") == tparse("2079 Jan 1")
                @test tparse(yctx, "80 Jan 1") == tparse("1980 Jan 1")
                @test tparse("80 Jan 1") == tparse("1980 Jan 1") && tparse("68 Jan 1") == tparse("2068 Jan 1")
                tsetyr(yctx, 0); tsetyr(yctx, -5)                   # no change, as in CSPICE
                @test str2et(yctx, "79 Jan 1") == str2et(yctx, "2079 Jan 1")
                timdef(yctx, "SET", "CALENDAR", "MIXED")            # keeps the interval
                @test str2et(yctx, "79 Jan 1") == str2et(yctx, "2079 Jan 1")
                @test str2et(other, "79 Jan 1 TDB") == str2et(other, "1979 Jan 1 TDB")
            finally
                SPICE.tsetyr(1969)
            end
            # `kclear` does not reset the defaults, as in CSPICE.
            kclear(other); timdef(other, "SET", "ZONE", "PST"); furnsh(other, Kernels.LSK)
            @test timdef(other, "GET", "ZONE") == "UTC-8" && timdef(other, "GET", "SYSTEM") == ""
        finally
            reset_timdef(ctx)
        end
    end

    @apart @testset "timdef: readers on other threads during a change" begin
        if Threads.nthreads() < 2
            skip_notice("9. timdef and threads", "Threads.nthreads() == 1")
        else
            c2 = neo_context([Kernels.LSK])
            s = "2000 Jan 1 12:00:00"
            a = str2et(c2, s)                       # UTC default
            timdef(c2, "SET", "ZONE", "UTC+5:30")
            b = str2et(c2, s)
            pa = timout(c2, 0.0, "HR:MN:SC")
            timdef(c2, "SET", "SYSTEM", "UTC")
            pb = timout(c2, 0.0, "HR:MN:SC")
            @test a != b && pa != pb
            stop = Threads.Atomic{Bool}(false)
            wrong = Threads.Atomic{Int}(0)
            count = Threads.Atomic{Int}(0)
            readers = [Threads.@spawn begin
                while !stop[]
                    v = str2et(c2, s)
                    p = timout(c2, 0.0, "HR:MN:SC")
                    (v == a || v == b) && (p == pa || p == pb) || Threads.atomic_add!(wrong, 1)
                    Threads.atomic_add!(count, 1)
                    yield()    # the root task stays on thread 1 and must get it back
                end
            end for _ in 1:max(2, Threads.nthreads() - 1)]
            changes = 0
            # Loop invariant: `changes` defaults were stored, each one complete.
            while changes < nsample(300, 3000) || (count[] < nsample(200, 2000) && changes < 3_000_000)
                changes += 1
                isodd(changes) ? timdef(c2, "SET", "ZONE", "UTC+5:30") : timdef(c2, "SET", "SYSTEM", "UTC")
            end
            stop[] = true
            foreach(wait, readers)
            println("  conversions during ", changes, " changes of the defaults: ", 2 * count[], ", wrong: ", wrong[])
            @test wrong[] == 0 && count[] > 0
        end
    end
    reset_timdef(ctx)
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
