# Time interface on a context. `time.jl`, `timeparse.jl` and `timeout.jl` define the
# forms on a `KernelSet`; the forms here read the snapshot of the kernels and the
# snapshot of the time defaults of the context, each with one atomic read.

"""
    timedefaults(ctx::SpiceContext) -> TimeDefaults

Contract: returns the current time defaults of `ctx` with one atomic read; the result
never changes.
"""
timedefaults(ctx::SpiceContext) = @atomic ctx.timedefs

"""
    timdef([ctx], action, item, value="") -> String

Contract: `action` is `"SET"` or `"GET"` and `item` is `"CALENDAR"`, `"SYSTEM"` or
`"ZONE"`, in any case. `"GET"` returns the default of `ctx` for `item` (`""` when the
item is not set). `"SET"` stores `value` and returns the stored form: a calendar
`"GREGORIAN"`, `"JULIAN"` or `"MIXED"`; a system `"UTC"`, `"TDB"`, `"TDT"` or `"TT"`
(the zone becomes `""`); a zone `"EST"`, `"EDT"`, `"CST"`, `"CDT"`, `"MST"`, `"MDT"`,
`"PST"`, `"PDT"`, `"UTC+h[:mm]"` or `"UTC-h[:mm]"`, stored as an offset (the system
becomes `""`). The defaults are those of `str2et` and `timout` for a string or a
picture that names no system, zone or calendar, as in CSPICE `timdef`; here they
belong to the context. A `"SET"` builds a new immutable snapshot under the lock of the
context and swaps it atomically, so readers on other threads see the old or the new
defaults, never a mixture, and never wait. Throws `SpiceError(:badtimedefault)` for
another action, item or value; the defaults are then unchanged.
"""
function timdef(ctx::SpiceContext, action::AbstractString, item::AbstractString,
        value::AbstractString="")
    act = uppercase(strip(action, ' '))
    act == "GET" && return _timdef_get(timedefaults(ctx), item)
    act == "SET" || throw(SpiceError(:badtimedefault,
        "the action '$action' is neither 'SET' nor 'GET'"))
    return lock(ctx.lock) do
        defs = _timdef_set(timedefaults(ctx), item, value)
        @atomic ctx.timedefs = defs
        _timdef_get(defs, item)
    end
end
timdef(action::AbstractString, item::AbstractString, value::AbstractString="") =
    timdef(default_context(), action, item, value)

# Context and default-context forms.
deltet(ctx::SpiceContext, epoch::Real, eptype::AbstractString) =
    deltet(kernelset(ctx), epoch, eptype)
deltet(epoch::Real, eptype::AbstractString) = deltet(default_context(), epoch, eptype)

unitim(ctx::SpiceContext, epoch::Real, insys::AbstractString, outsys::AbstractString) =
    unitim(kernelset(ctx), epoch, insys, outsys)
unitim(epoch::Real, insys::AbstractString, outsys::AbstractString) =
    unitim(default_context(), epoch, insys, outsys)

ttrans(ctx::SpiceContext, from::AbstractString, to::AbstractString, tv) =
    ttrans(kernelset(ctx), from, to, tv)

utc2et(ctx::SpiceContext, s::AbstractString) = utc2et(kernelset(ctx), s, timedefaults(ctx))
utc2et(s::AbstractString) = utc2et(default_context(), s)

str2et(ctx::SpiceContext, s::AbstractString) = str2et(kernelset(ctx), s, timedefaults(ctx))
str2et(s::AbstractString) = str2et(default_context(), s)

et2utc(ctx::SpiceContext, et::Real, format::AbstractString, prec::Integer) =
    et2utc(kernelset(ctx), et, format, prec)
et2utc(et::Real, format::AbstractString, prec::Integer) =
    et2utc(default_context(), et, format, prec)

timout(ctx::SpiceContext, et::Real, pictur::AbstractString) =
    timout(kernelset(ctx), et, pictur, timedefaults(ctx))
timout(et::Real, pictur::AbstractString) = timout(default_context(), et, pictur)

"""
    tsetyr([ctx], year) -> Nothing

Set the first year of the 100 years that a year of two digits abbreviates in `str2et`,
`utc2et` and `tparse` of the context, as `tsetyr_c`: after `tsetyr(ctx, 1980)` the year
80 is 1980 and the year 79 is 2079. A `year` that is not positive changes nothing, as in
CSPICE. The value belongs to the context (CSPICE keeps one for the process), and a
reader sees the old or the new value, never a mixture.
"""
function tsetyr(ctx::SpiceContext, year::Integer)
    year > 0 || return nothing
    lock(ctx.lock) do
        d = timedefaults(ctx)
        @atomic ctx.timedefs = TimeDefaults(d.calendar, d.system, d.zone, Int(year))
    end
    @assert timedefaults(ctx).yearbound > 0
    return nothing
end
tsetyr(year::Integer) = tsetyr(default_context(), year)

# These three need no kernel; the context forms exist so that every call has them. The
# context form of `tparse` uses the interval of abbreviated years of the context.
etcal(::Union{SpiceContext,KernelSet}, et::Real) = etcal(et)
tparse(ctx::SpiceContext, s::AbstractString) = tparse(s, timedefaults(ctx))
tparse(::KernelSet, s::AbstractString) = tparse(s)
tpictr(::Union{SpiceContext,KernelSet}, sample::AbstractString) = tpictr(sample)
