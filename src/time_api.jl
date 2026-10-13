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
    timdef([ctx,] action, item, value="") -> String

Set or get a default of the time conversions of `ctx`, as `timdef_c`.

`action` is `"SET"` or `"GET"` and `item` is `"CALENDAR"`, `"SYSTEM"` or `"ZONE"`;
letter case is ignored. `"GET"` returns the default (`""` if the item is not set).
`"SET"` stores `value` and returns the stored form.

| `item`       | Values                                                                 | At start      |
|:-------------|:-----------------------------------------------------------------------|:--------------|
| `"CALENDAR"` | `"GREGORIAN"`, `"JULIAN"`, `"MIXED"` (Julian before 1582 OCT 15)       | `"GREGORIAN"` |
| `"SYSTEM"`   | `"UTC"`, `"TDB"`, `"TDT"`, `"TT"`                                      | `"UTC"`       |
| `"ZONE"`     | `"EST"`, `"EDT"`, `"CST"`, `"CDT"`, `"MST"`, `"MDT"`, `"PST"`, `"PDT"`, `"UTC+h[:mm]"`, `"UTC-h[:mm]"` | `""` |

A zone is stored as an offset, such as `"UTC-8"`. A system clears the zone and a zone
clears the system.

The defaults apply to [`str2et`](@ref) for a string that names no system or zone, and
to [`timout`](@ref) for a picture that names no system, zone or calendar. They belong
to the context; CSPICE keeps one set for the process. A `"SET"` is atomic for readers
on other threads.

Throws `SpiceError(:badtimedefault)` if the action, the item or the value is not one
of the above; the defaults are then unchanged.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> timdef(ctx, "GET", "SYSTEM")
"UTC"

julia> timdef(ctx, "SET", "ZONE", "PST")
"UTC-8"

julia> timdef(ctx, "GET", "SYSTEM")
""
```
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
    tsetyr([ctx,] year) -> Nothing

Set the first year of the 100 years that a year of two digits stands for, as
`tsetyr_c`.

At start the value is 1969: the year 69 is 1969 and the year 68 is 2068. After
`tsetyr(ctx, 1980)` the year 80 is 1980 and the year 79 is 2079. The value applies to
[`str2et`](@ref), [`utc2et`](@ref) and [`tparse`](@ref) on the context, and to the `@`
dates of the text kernels that the context loads afterwards ([`furnsh`](@ref),
[`ldpool`](@ref), [`lmpool`](@ref)). A `year` that is not positive changes nothing, as
in CSPICE. The value belongs to the context; CSPICE keeps one for the process.

# Examples
```jldoctest
julia> ctx = SpiceContext();

julia> etcal(tparse(ctx, "79 JAN 1"))
"1979 JAN 01 00:00:00.000"

julia> tsetyr(ctx, 1980)

julia> etcal(tparse(ctx, "79 JAN 1"))
"2079 JAN 01 00:00:00.000"
```
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
