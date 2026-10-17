# Double precision numbers as strings in the base 16 "scientific notation" of SPICELIB
# (`dp2hx.c`, `hx2dp.c`, with `int2hx.c` and `hx2int.c` for the exponent). A string is
# a sign, hexadecimal digits of a mantissa in [1/16, 1), a caret and a hexadecimal
# exponent of 16: "2A^3" is (2/16 + 10/16²)·16³ = 672. The arithmetic follows the f2c
# source step by step, so a number and its string are those of CSPICE. Every private
# helper has the prefix `_hx_`.

const _HX_DIGITS = "0123456789ABCDEF"

"""
    dp2hx(number) -> String

Return the base 16 string of `number`, as `dp2hx_c` of CSPICE.

The string is an optional minus sign, the hexadecimal digits of the mantissa in
[1/16, 1), a caret `^` and the exponent of 16 in hexadecimal, as in `"2A^3"` for 672.
The digits are exact: every double has a finite string. Zero, and `NaN` as in CSPICE,
give `"0^0"`. [`hx2dp`](@ref) converts the string back.

Throws `SpiceError(:valueoutofrange)` for an infinite `number`; `dp2hx_c` does not
return for it.

# Examples
```jldoctest
julia> dp2hx(672.0), dp2hx(-11.0), dp2hx(0.1)
("2A^3", "-B^1", "1999999999999A^0")
```
"""
function dp2hx(number::Real)
    x = Float64(number)
    isinf(x) && throw(SpiceError(:valueoutofrange,
        "dp2hx: the number $x has no base 16 string; dp2hx_c does not return for it"))
    (x < 0.0 || x > 0.0) || return "0^0"
    # Normalize to 1/16 <= |x| < 1 with steps of 16^8, 16^4 and 16, as CSPICE does.
    # Every step is exact, because the factors are powers of 2.
    e = 0
    a = abs(x)
    if a < 1.0
        while a * 4294967296.0 < 1.0        # invariant: |x| = a * 16^e
            a *= 4294967296.0
            e -= 8
        end
        while a * 65536.0 < 1.0
            a *= 65536.0
            e -= 4
        end
        while a * 16.0 < 1.0
            a *= 16.0
            e -= 1
        end
    else
        while a * 2.3283064365386963e-10 >= 1.0
            a *= 2.3283064365386963e-10
            e += 8
        end
        while a * 1.52587890625e-5 >= 1.0
            a *= 1.52587890625e-5
            e += 4
        end
        while a >= 1.0
            a *= 0.0625
            e += 1
        end
    end
    @assert 0.0625 <= a < 1.0
    io = IOBuffer()
    x < 0.0 && print(io, '-')
    # The remainder becomes zero, because the mantissa is a sum of powers of 2.
    while a != 0.0                          # invariant: the digits written so far are exact
        t = a * 16.0
        digit = trunc(Int, t)
        a = t - Float64(digit)
        print(io, _HX_DIGITS[digit+1])
    end
    print(io, '^', e < 0 ? "-" : "", uppercase(string(abs(e); base=16)))
    return String(take!(io))
end

"Value of the hexadecimal digit `b` (a byte), or -1 if it is not one; letter case is ignored."
@inline function _hx_digit(b::UInt8)
    UInt8('0') <= b <= UInt8('9') && return Int(b - UInt8('0'))
    UInt8('A') <= b <= UInt8('F') && return Int(b - UInt8('A')) + 10
    UInt8('a') <= b <= UInt8('f') && return Int(b - UInt8('a')) + 10
    return -1
end

"The message of CSPICE for the byte at index `i` of `s`, which is not a digit."
function _hx_illegal(s::AbstractString, b::AbstractVector{UInt8}, i::Int)
    c = b[i] < 0x80 ? Char(b[i]) : (isvalid(s, i) ? s[i] : '�')
    return "ERROR: Illegal character '$c' encountered."
end

"""
    _hx_int(s, b, start) -> (value, errmsg)

The integer of the hexadecimal exponent in the bytes `b[start:end]` of `s`, as
`hx2int_`: blanks before and after are skipped, the sign is optional, and the value must
fit in 32 bits. `errmsg` is `nothing`, or the message of CSPICE and then `value` is 0.
"""
function _hx_int(s::AbstractString, b::AbstractVector{UInt8}, start::Int)
    n = length(b)
    lo = findnext(!=(UInt8(' ')), b, start)
    lo === nothing && return 0, "ERROR: A blank input string is not allowed."
    hi = findprev(!=(UInt8(' ')), b, n)
    pos = lo
    neg = b[pos] == UInt8('-')
    (neg || b[pos] == UInt8('+')) && (pos += 1)
    # The limits of `hx2int_` for 32-bit integers: INTMIN/16 with remainder 0 and
    # INTMAX/16 with remainder 15.
    mini, minmod = -134217728, 0
    maxi, maxmod = 134217727, 15
    value = 0
    while pos <= hi                         # invariant: value holds the digits before pos
        d = _hx_digit(b[pos])
        d < 0 && return 0, _hx_illegal(s, b, pos)
        if neg
            (value > mini || (value == mini && d <= minmod)) ||
                return 0, "ERROR: Integer too small to be represented."
            value = 16value - d
        else
            (value < maxi || (value == maxi && d <= maxmod)) ||
                return 0, "ERROR: Integer too large to be represented."
            value = 16value + d
        end
        pos += 1
    end
    return value, nothing
end

"""
    hx2dp(str) -> Float64

Return the double of the base 16 string `str`, as `hx2dp_c` of CSPICE.

The string has the form that [`dp2hx`](@ref) writes: an optional sign, at most 31
hexadecimal digits of the mantissa, a caret `^` and the exponent of 16 in hexadecimal
with an optional sign. Letter case is ignored, and blanks before and after the string
and after the caret are skipped. The digits are read as a fraction: `"2A^3"` is
(2/16 + 10/16²)·16³ = 672.

Throws `SpiceError(:emptystring)` if `str` is empty, and `SpiceError(:badnumber)` if
CSPICE sets its error flag; the message is then the one that `hx2dp_c` returns in
`errmsg`.

# Examples
```jldoctest
julia> hx2dp("2A^3"), hx2dp("-B^1"), hx2dp(dp2hx(0.1)) == 0.1
(672.0, -11.0, true)
```
"""
function hx2dp(str::AbstractString)
    isempty(str) && throw(SpiceError(:emptystring, "hx2dp: the string is empty"))
    bad(msg) = throw(SpiceError(:badnumber, msg))
    b = codeunits(str)
    n = length(b)
    lo = findfirst(!=(UInt8(' ')), b)
    lo === nothing && bad("ERROR: A blank input string is not allowed.")
    hi = findlast(!=(UInt8(' ')), b)
    pos = lo
    neg = b[pos] == UInt8('-')
    (neg || b[pos] == UInt8('+')) && (pos += 1)
    digits = Int[]
    fndexp = false
    while pos <= hi && !fndexp              # invariant: digits holds the mantissa before pos
        d = _hx_digit(b[pos])
        if d >= 0
            push!(digits, d)
        elseif b[pos] == UInt8('^')
            fndexp = true
        else
            bad(_hx_illegal(str, b, pos))
        end
        pos += 1
        length(digits) > 31 && bad("ERROR: Too many digits in the mantissa (> 31).")
    end
    pos <= hi || bad("ERROR: Missing exponent.")
    expon, msg = _hx_int(str, b, pos)
    msg === nothing || bad(msg)
    # The digits are added from the least significant one, as in CSPICE.
    scales = accumulate(*, fill(0.0625, length(digits)))
    t = 0.0
    for k in length(digits):-1:1
        t = neg ? t - Float64(digits[k]) * scales[k] : t + Float64(digits[k]) * scales[k]
    end
    if expon < 0
        # CSPICE negates the exponent in 32-bit arithmetic. For -2^31 that overflows, and
        # the compiled library then makes one step (measured with CSPICE N0067).
        steps = expon == -2147483648 ? 1 : -expon
        for _ in 1:steps                    # invariant: t is scaled by the steps made
            t == 0.0 && break               # no later step changes a zero
            t *= 0.0625
        end
    else
        mindp = -floatmax(Float64) * 0.0625
        maxdp = floatmax(Float64) * 0.0625
        for _ in 1:expon
            t == 0.0 && break
            if neg
                t >= mindp || bad("ERROR: Number is too small to be represented.")
            else
                t <= maxdp || bad("ERROR: Number is too large to be represented.")
            end
            t *= 16.0
        end
    end
    return t
end
