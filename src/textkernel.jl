# Text kernels: architecture sniffing, the kernel-pool grammar of `pool.c`/`zzrvar.c`,
# meta-kernel expansion as `keeper.c` does it, and the number and date parsers
# (`nparsd.c`, `tparse.c`) that fix the exact value of every numeric pool variable.

const _TK_LINLEN = 132      # LINLEN of pool.c: longer lines are cut
const _TK_MAXNAME = 32      # MAXLEN: longest variable name
const _TK_MAXSTR = 80       # MAXCHR: longest string value

# ---------------------------------------------------------------------------------
# Architecture
# ---------------------------------------------------------------------------------

"""
    kernel_arch(path) -> Symbol

Contract: `path` names a readable file; returns `:SPK`, `:PCK` or `:CK` for a binary
DAF of that type (ID word `DAF/SPK`, `DAF/PCK`, `DAF/CK`, or the old `NAIF/DAF` with
summary format ND = 2 and NI = 6 or 5; with NI = 6 the first segment decides between
SPK and CK, as `open_daf` does), `:DSK` for a DAS file with the ID word `DAS/DSK`, and
`:TEXT` for everything that is not a DAF or DAS file.
Throws `SpiceError(:nofile)` for a missing file and `SpiceError(:unsupported)` for
another binary kernel type.
"""
function kernel_arch(path::AbstractString)::Symbol
    isfile(path) || throw(SpiceError(:nofile, "kernel file '$path' does not exist"))
    head = open(io -> read(io, 16), path)
    n = length(head)
    idword = String(head[1:min(n, 8)])
    if startswith(idword, "DAF/")
        typ = strip(idword[5:end])
        typ == "SPK" && return :SPK
        typ == "PCK" && return :PCK
        typ == "CK" && return :CK
        throw(SpiceError(:unsupported, "DAF kernel type '$typ' of '$path' is not supported"))
    elseif idword == "NAIF/DAF"
        n >= 16 || throw(SpiceError(:badfile, "DAF file '$path' is truncated"))
        nd = reinterpret(Int32, head[9:12])[1]
        ni = reinterpret(Int32, head[13:16])[1]
        for (d, i) in ((nd, ni), (bswap(nd), bswap(ni)))
            d == 2 && i == 6 && return open_daf(path).kind
            d == 2 && i == 5 && return :PCK
        end
        throw(SpiceError(:unsupported,
            "old-style DAF '$path' has a summary format that is neither SPK nor PCK"))
    elseif startswith(idword, "DAS/")
        strip(idword[5:end]) == "DSK" && return :DSK
        throw(SpiceError(:unsupported, "DAS kernel '$path' is not supported"))
    end
    return :TEXT
end

# ---------------------------------------------------------------------------------
# Numbers: port of NPARSD
# ---------------------------------------------------------------------------------

const _NPARSD_LOOKUP = (1.0, 10.0, 100.0, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10)

"""
    _nparsd(s) -> Float64

Contract: `s` is a number in the syntax of CSPICE `nparsd` (digits, one decimal
point, `E`/`D` exponent, signs, `PI`, embedded blanks and commas ignored); returns
the same double that `nparsd` computes, which is not always the correctly rounded
one. Throws `SpiceError(:badnumber)` otherwise.
"""
function _nparsd(s::Union{AbstractString,AbstractVector{UInt8}})::Float64
    bytes = s isa AbstractString ? codeunits(s) : s
    bad(i) = throw(SpiceError(:badnumber,
        "'$(String(copy(bytes)))' is not a number (character $i)"))
    intbnd = 1e15
    maxexp = 308.0
    minexp = -309.0
    dpmax = floatmax(Float64)
    smlbnd = dpmax / 1e10
    pntok = signok = expok = bpiok = doint = roundd = roundi = true
    dodec = doexp = mantsa = sigchr = epiok = zeroi = false
    intval = decval = expval = ecount = 0.0
    divisr = factor = 1.0
    signdx = 1
    sgn1 = sgn2 = 1.0
    l = length(bytes)
    while l > 0 && bytes[l] == UInt8(' ')
        l -= 1
    end
    l == 0 && throw(SpiceError(:badnumber, "a blank string is not a number"))
    @inbounds for i in 1:l
        c = bytes[i]
        if UInt8('0') <= c <= UInt8('9')
            v = Float64(c - UInt8('0'))
            bpiok = epiok = signok = false
            sigchr = true
            if doint
                mantsa = true
                if intval < intbnd
                    intval = intval * 10.0 + v
                else
                    ecount += 1.0
                    factor /= 10.0
                    if roundi
                        roundi = false
                        v > 5.0 && (intval += 1.0)
                    end
                end
            elseif dodec
                mantsa = true
                if zeroi
                    if decval < intbnd
                        decval = decval * 10.0 + v
                        ecount -= 1.0
                    elseif roundd
                        roundd = false
                        v >= 5.0 && (decval += 1.0)
                    end
                elseif divisr < intbnd
                    decval = decval * 10.0 + v
                    divisr *= 10.0
                end
            elseif doexp
                if expval + ecount > maxexp
                    throw(SpiceError(:badnumber,
                        "'$(String(copy(bytes)))' is too large for a double"))
                elseif expval + ecount < minexp
                else
                    expval = expval * 10.0 + sgn2 * v
                end
            else
                bad(i)
            end
        elseif c == UInt8('.')
            pntok || bad(i)
            bpiok = epiok = pntok = signok = false
            dodec = true
            doint = doexp = false
            zeroi = intval == 0.0
        elseif c == UInt8('E') || c == UInt8('D') || c == UInt8('e') || c == UInt8('d')
            sigchr = true
            expok || bad(i)
            bpiok = epiok = expok = pntok = dodec = doint = false
            doexp = signok = true
            signdx = 2
        elseif c == UInt8('+') || c == UInt8('-')
            signok || bad(i)
            v = c == UInt8('-') ? -1.0 : 1.0
            signdx == 1 ? (sgn1 = v) : (sgn2 = v)
            signok = false
        elseif c == UInt8('p') || c == UInt8('P')
            sigchr = true
            bpiok || bad(i)
            doint = dodec = doexp = expok = pntok = bpiok = signok = false
            epiok = true
        elseif c == UInt8('i') || c == UInt8('I')
            epiok || bad(i)
            doint = dodec = doexp = expok = pntok = bpiok = signok = epiok = false
            mantsa = true
            intval = Float64(pi)
        elseif c == UInt8(' ') || c == UInt8(',')
        else
            bad(i)
        end
    end
    epiok && bad(l)
    value = if mantsa
        intval + decval / divisr * factor
    elseif sigchr
        1.0
    else
        bad(l)
    end
    expval += ecount
    if expval < 0.0
        while expval < -10.0
            value /= 1e10
            expval += 10.0
        end
        value /= _NPARSD_LOOKUP[Int(-trunc(expval))+1]
    elseif expval > 0.0
        while expval > 10.0
            value >= smlbnd && throw(SpiceError(:badnumber,
                "'$(String(copy(bytes)))' is too large for a double"))
            value *= 1e10
            expval -= 10.0
        end
        f = _NPARSD_LOOKUP[Int(round(expval))+1]
        value < dpmax / f || throw(SpiceError(:badnumber,
            "'$(String(copy(bytes)))' is too large for a double"))
        value *= f
    end
    return sgn1 * value
end

# ---------------------------------------------------------------------------------
# Dates: the `@` forms, as TPARSE converts them
# ---------------------------------------------------------------------------------

const _MONTH_ABBR = ("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT",
    "NOV", "DEC")

"Contract: `y` is a two-digit year; returns the year in 1969..2068, as `texpyr` does."
_texpyr(y::Int) = y >= 69 ? y + 1900 : y + 2000

"""
    _tparse_seconds(year, month, day, h, m, s) -> Float64

Contract: a Gregorian calendar date with `year >= 1`; returns the seconds past J2000
on a uniform 86400 s day scale with the arithmetic of CSPICE `tparse`.
"""
function _tparse_seconds(year::Int, month::Int, day::Int, h::Float64, m::Float64, s::Float64)
    year >= 1 || throw(SpiceError(:badtime, "years before 1 A.D. are not supported"))
    d = 367 * year - div(7 * (year + div(month + 9, 12)), 4) -
        div(3 * (div(year + div(month - 9, 7), 100) + 1), 4) + div(275 * month, 9) + day -
        730516
    return (Float64(d) - 0.5) * 86400.0 + h * 3600.0 + m * 60.0 + s
end

"""
    _tparse(s) -> Float64

Contract: `s` is a kernel date without the leading `@` in one of the forms
`YYYY-MON-DD`, `DD-MON-YYYY`, `MON-DD-YYYY`, `YYYY-MM-DD`, each with an
optional time `hh:mm[:ss.fff]` after `-`, `/`, `T` or a blank, or `JD<number>`;
returns seconds past J2000 without leap seconds, as CSPICE `tparse` does. Throws
`SpiceError(:badtime)` for any other form.
"""
function _tparse(s::AbstractString)::Float64
    str = uppercase(strip(s))
    bad() = throw(SpiceError(:badtime, "cannot parse the kernel date '@$s'"))
    isempty(str) && bad()
    if startswith(str, "JD")
        return (_nparsd(strip(str[3:end])) - 2451545.0) * 86400.0
    end
    str = replace(str, r"(?<=\d)T(?=\d)" => " ")
    pieces = split(str, r"[-/,\s]+"; keepempty=false)
    h = mi = sec = 0.0
    datep = SubString{String}[]
    for p in pieces
        if occursin(':', p)
            t = split(p, ':')
            2 <= length(t) <= 3 || bad()
            h = _nparsd(t[1])
            mi = _nparsd(t[2])
            sec = length(t) == 3 ? _nparsd(t[3]) : 0.0
        else
            push!(datep, p)
        end
    end
    isnum(p) = all(isdigit, p)
    monthof(p) = begin
        length(p) >= 3 || bad()
        k = findfirst(==(p[1:3]), _MONTH_ABBR)
        k === nothing && bad()
        k
    end
    local year::Int, month::Int, day::Int
    length(datep) == 3 || bad()
    im = findfirst(!isnum, datep)
    if im === nothing
        a, b, c = parse(Int, datep[1]), parse(Int, datep[2]), parse(Int, datep[3])
        # A leading number of 1000 or more is a year; otherwise the order is M-D-Y.
        year, month, day = a >= 1000 ? (a, b, c) : (c, a, b)
    else
        month = monthof(datep[im])
        rest = [datep[i] for i in 1:3 if i != im]
        all(isnum, rest) || bad()
        a, b = parse(Int, rest[1]), parse(Int, rest[2])
        # `MON D Y` always; around a month name the first number is the year unless
        # the last one is 1000 or more (`15-MAR-87` is year 2015, day 87).
        yearlast = im == 1 || b >= 1000
        year, day = yearlast ? (b, a) : (a, b)
    end
    year < 100 && (year = _texpyr(year))
    (1 <= month <= 12 && day >= 1) || bad()
    return _tparse_seconds(year, month, day, h, mi, sec)
end

# ---------------------------------------------------------------------------------
# Grammar: port of the scanner of ZZRVAR
# ---------------------------------------------------------------------------------

# Token types of ZZRVAR.
const _TOK_QUOTED = 1
const _TOK_VALUE = 2
const _TOK_LPAREN = 3
const _TOK_RPAREN = 4
const _TOK_EQUAL = 5
const _TOK_PLUSEQ = 6

"""
    _tokenize(line) -> (begs, ends, types)

Contract: `line` is one data line as bytes; returns the token boundaries (1-based,
inclusive) and types exactly as the scanner loop of CSPICE `zzrvar` finds them.
"""
function _tokenize(line::Vector{UInt8})
    begs = Int[]
    ends = Int[]
    types = Int[]
    even = true
    intokn = false
    inquot = false
    insepf = true
    @inbounds for i in eachindex(line)
        c = line[i]
        if c == UInt8(' ') || c == UInt8(',') || c == UInt8('(') || c == UInt8(')') ||
           c == UInt8('=') || c == 0x09
            if intokn || (inquot && even)
                inquot = false
                intokn = false
                insepf = true
            end
            if insepf
                if c == UInt8('=')
                    push!(begs, i)
                    push!(ends, i)
                    push!(types, _TOK_EQUAL)
                    if i > 1 && line[i-1] == UInt8('+')
                        n = length(begs)
                        types[n] = _TOK_PLUSEQ
                        begs[n] = i - 1
                        if n >= 2 && begs[n-1] == ends[n-1]
                            pop!(begs)
                            pop!(ends)
                            pop!(types)
                            begs[n-1] = i - 1
                            ends[n-1] = i
                            types[n-1] = _TOK_PLUSEQ
                        elseif n >= 2
                            ends[n-1] -= 1
                        end
                    end
                elseif c == UInt8(')')
                    push!(begs, i)
                    push!(ends, i)
                    push!(types, _TOK_RPAREN)
                elseif c == UInt8('(')
                    push!(begs, i)
                    push!(ends, i)
                    push!(types, _TOK_LPAREN)
                end
            end
        elseif c == UInt8('\'')
            even = !even
            if !inquot
                insepf = false
                intokn = false
                inquot = true
                push!(begs, i)
                push!(ends, i)
                push!(types, _TOK_QUOTED)
            end
            ends[end] = i
        else
            if insepf || (inquot && even)
                inquot = false
                insepf = false
                intokn = true
                push!(begs, i)
                push!(ends, i)
                push!(types, _TOK_VALUE)
            end
            ends[end] = i
        end
    end
    return begs, ends, types
end

"Contract: `line[b:e]` is a quoted token; returns its value as `zzrvar` stores it (outer quotes removed, `''` reduced to `'`, at most 80 characters, trailing blanks removed)."
function _string_value(line::Vector{UInt8}, b::Int, e::Int)::String
    line[e] == UInt8('\'') || (e += 1)
    out = UInt8[]
    j = b + 1
    while j < e
        line[j] == UInt8('\'') && (j += 1)
        j <= length(line) || break
        length(out) < _TK_MAXSTR && push!(out, line[j])
        j += 1
    end
    while !isempty(out) && out[end] == UInt8(' ')
        pop!(out)
    end
    return String(out)
end

"Contract: `raw` is one physical line; returns the line as `rdkdat` hands it on (cut at 132 characters, tabs replaced by blanks, blanks removed at both ends)."
function _clean_line(raw::AbstractString)::Vector{UInt8}
    b = Vector{UInt8}(codeunits(raw))
    while !isempty(b) && (b[end] == UInt8('\r') || b[end] == UInt8('\n'))
        pop!(b)
    end
    length(b) > _TK_LINLEN && resize!(b, _TK_LINLEN)
    for i in eachindex(b)
        b[i] == 0x09 && (b[i] = UInt8(' '))
    end
    hi = length(b)
    while hi > 0 && b[hi] == UInt8(' ')
        hi -= 1
    end
    lo = 1
    while lo <= hi && b[lo] == UInt8(' ')
        lo += 1
    end
    return b[lo:hi]
end

const _BEGINDATA = Vector{UInt8}(codeunits("\\begindata"))
const _BEGINTEXT = Vector{UInt8}(codeunits("\\begintext"))

"""
    parse_text_lines(path, lines) -> TextKernel

Contract: `lines` iterates the physical lines of a text kernel; returns the
assignments in file order with the values CSPICE puts in the kernel pool. `path` is
used for the result and for error messages only. Throws `SpiceError(:badkernel)` for
a malformed assignment.
"""
function parse_text_lines(path::AbstractString, lines)::TextKernel
    out = Assignment[]
    indata = false
    invec = false
    name = ""
    append = false
    nums = Float64[]
    strs = String[]
    kind = 0                       # 0 unknown, 1 string, 2 numeric
    lineno = 0
    fail(msg) = throw(SpiceError(:badkernel, "$msg (line $lineno of '$path')"))
    finish() = begin
        kind == 1 && push!(out, Assignment(name, append, strs))
        kind == 2 && push!(out, Assignment(name, append, nums))
        invec = false
    end
    for raw in lines
        lineno += 1
        line = _clean_line(raw)
        if line == _BEGINTEXT
            indata = false
            continue
        elseif line == _BEGINDATA
            indata = true
            continue
        end
        (indata && !isempty(line)) || continue
        begs, ends, types = _tokenize(line)
        count = length(begs)
        nxt = 1
        if !invec
            count >= 3 || fail("a kernel variable assignment must have the form " *
                                "'<name> [+]= <values>'")
            ends[1] - begs[1] + 1 <= _TK_MAXNAME ||
                fail("a variable name is longer than $_TK_MAXNAME characters")
            name = String(line[begs[1]:ends[1]])
            types[2] == _TOK_EQUAL || types[2] == _TOK_PLUSEQ ||
                fail("the assignment operator of '$name' is neither '=' nor '+='")
            append = types[2] == _TOK_PLUSEQ
            nums = Float64[]
            strs = String[]
            kind = 0
            if types[3] == _TOK_LPAREN
                nxt = 4
            else
                nxt = 3
                push!(begs, 0)
                push!(ends, 0)
                push!(types, _TOK_RPAREN)
                count += 1
            end
        end
        # Loop invariant: tokens 1..nxt-1 of this line are consumed and every value
        # among them is in `nums` or `strs`.
        while nxt <= count && types[nxt] != _TOK_RPAREN
            t, b, e = types[nxt], begs[nxt], ends[nxt]
            if kind == 0
                kind = t == _TOK_QUOTED ? 1 : t == _TOK_VALUE ? 2 :
                       fail("the first item after the assignment operator of '$name' " *
                            "is not a value")
            end
            if kind == 1
                t == _TOK_QUOTED || fail("'$name' is a string variable but a value " *
                                          "is not a string")
                b + 1 >= e && fail("there is a quoted string with no characters")
                push!(strs, _string_value(line, b, e))
            else
                t == _TOK_VALUE || fail("'$name' is a numeric variable but a value " *
                                         "is not a number or a time")
                if line[b] == UInt8('@')
                    e == b && fail("'@' is not followed by a time")
                    push!(nums, _tparse(String(line[b+1:e])))
                else
                    push!(nums, _nparsd(view(line, b:e)))
                end
            end
            nxt += 1
        end
        nxt > count ? (invec = true) : finish()
    end
    invec && finish()
    return TextKernel(String(path), out)
end

"""
    parse_text_kernel(path) -> TextKernel

Contract: `path` names a readable text kernel; returns its assignments in file order.
Only text between `\\begindata` and `\\begintext` is data. Values are numbers (with
`D` or `E` exponents), `@` dates as TDB seconds past J2000, or quoted strings with
`''` as an escaped quote; string elements are kept raw, so a continuation marker is
still in place (see `pool_string_joined`). Variable names are case-sensitive.
"""
function parse_text_kernel(path::AbstractString)::TextKernel
    isfile(path) || throw(SpiceError(:nofile, "kernel file '$path' does not exist"))
    text = read(path, String)
    return parse_text_lines(path, eachsplit(text, r"\r\n|\n|\r"))
end

# ---------------------------------------------------------------------------------
# Pool updates and meta-kernels
# ---------------------------------------------------------------------------------

"""
    apply_assignments!(pool, tk) -> pool

Contract: applies the assignments of `tk` to `pool` in order; `=` replaces a variable
and `+=` appends to it (or creates it). A stored vector is never mutated, so vectors
shared with an older snapshot stay valid. Throws `SpiceError(:typemismatch)` when
`+=` mixes numeric and string values.
"""
function apply_assignments!(pool::Dict{String,PoolValue}, tk::TextKernel)
    for a in tk.assignments
        old = a.append ? get(pool, a.name, nothing) : nothing
        if old === nothing
            pool[a.name] = a.value
        elseif old isa Vector{Float64} && a.value isa Vector{Float64}
            pool[a.name] = vcat(old, a.value)
        elseif old isa Vector{String} && a.value isa Vector{String}
            pool[a.name] = vcat(old, a.value)
        else
            throw(SpiceError(:typemismatch, "'+=' on kernel variable '$(a.name)' in " *
                "'$(tk.path)' mixes numeric and string values"))
        end
    end
    return pool
end

"""
    pool_strings_joined(values, contin) -> Vector{String}

Contract: `values` are the raw elements of a string variable and `contin` a non-empty
continuation marker; returns the logical strings, where an element that ends in
`contin` is joined (marker removed) with the next element, as CSPICE `stpool` does.
"""
function pool_strings_joined(values::Vector{String}, contin::AbstractString)::Vector{String}
    isempty(contin) && throw(ArgumentError("the continuation marker is empty"))
    out = String[]
    buf = IOBuffer()
    open_string = false
    for v in values
        part = rstrip(v, ' ')
        if endswith(part, contin)
            print(buf, chop(part; tail=length(contin)))
            open_string = true
        else
            print(buf, part)
            push!(out, String(take!(buf)))
            open_string = false
        end
    end
    open_string && push!(out, String(take!(buf)))
    return out
end

"""
    pool_string_joined(values, nth, contin = "//") -> Union{String,Nothing}

Contract: the `nth` (1-based) logical string of a string variable whose elements may
be continued with `contin`, or `nothing` if there is none (`stpool`).
"""
function pool_string_joined(values::Vector{String}, nth::Integer,
    contin::AbstractString="//")
    all = pool_strings_joined(values, contin)
    return 1 <= nth <= length(all) ? all[nth] : nothing
end

"""
    metakernel_files(tk) -> Vector{String}

Contract: returns the files of `KERNELS_TO_LOAD` of `tk` in load order, with every
`\$SYMBOL` of `PATH_SYMBOLS` replaced by its `PATH_VALUES` entry (longest symbol
first, `+` as the continuation marker), and an empty vector if `tk` is not a
meta-kernel. As in CSPICE, a relative path is returned unchanged, so it is relative
to the current directory of the process and not to the directory of the meta-kernel.
"""
function metakernel_files(tk::TextKernel)::Vector{String}
    pool = apply_assignments!(Dict{String,PoolValue}(), tk)
    ktl = get(pool, "KERNELS_TO_LOAD", nothing)
    ktl isa Vector{String} || return String[]
    files = pool_strings_joined(ktl, "+")
    syms = get(pool, "PATH_SYMBOLS", nothing)
    syms isa Vector{String} || return files
    vals = get(pool, "PATH_VALUES", nothing)
    vals isa Vector{String} || throw(SpiceError(:badkernel,
        "meta-kernel '$(tk.path)' has PATH_SYMBOLS but no PATH_VALUES"))
    values = pool_strings_joined(vals, "+")
    length(values) == length(syms) || throw(SpiceError(:badkernel,
        "meta-kernel '$(tk.path)' has $(length(syms)) path symbols and " *
        "$(length(values)) path values"))
    return [_substitute_symbols(f, syms, values) for f in files]
end

"Contract: returns `file` with each `\$SYMBOL` replaced by its value; at each `\$` the longest matching symbol wins and the inserted text is not scanned again."
function _substitute_symbols(file::String, syms::Vector{String}, values::Vector{String})
    out = IOBuffer()
    i = firstindex(file)
    while i <= lastindex(file)
        c = file[i]
        if c == '$'
            rest = SubString(file, nextind(file, i))
            use = 0
            for (k, s) in enumerate(syms)
                !isempty(s) && startswith(rest, s) &&
                    (use == 0 || length(s) > length(syms[use])) && (use = k)
            end
            if use > 0
                print(out, values[use])
                i = nextind(file, i, length(syms[use]) + 1)
                continue
            end
        end
        print(out, c)
        i = nextind(file, i)
    end
    return String(take!(out))
end
