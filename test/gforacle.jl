# CSPICE oracle of the geometry-finder groups (`test_gf_*.jl`).
#
# The functions call the `_c` routines of the C library directly, as `mathoracle.jl`
# does, for two reasons: the wrappers of SPICE.jl give the long error message only, and
# the tests need the short one (the `kind` of the `SpiceError` of NeoSpice); and the
# wrappers of `gfuds` and `gfudb` pass the scalar function to the callbacks, a form that
# NeoSpice does not use. A CSPICE error becomes a `GFO.CError` with the short message as
# a lower-case symbol.
#
# Every search returns the result window as a `Vector{Float64}` of end points. CSPICE is
# not thread safe, so these functions are called from the main task only.

module GFO

import SPICE
import Libdl
using SPICE: SpiceDoubleCell, SpiceDouble, SpiceInt, SpiceBoolean, Cell

struct CError <: Exception
    kind::Symbol
end

const HANDLE = Libdl.dlopen(SPICE.libcspice)

"Throw `CError` if CSPICE has signalled an error, and reset the error state."
function check()
    ccall(Libdl.dlsym(HANDLE, :failed_c), Cint, ()) == 0 && return nothing
    buf = zeros(UInt8, 64)
    ccall(Libdl.dlsym(HANDLE, :getmsg_c), Cvoid, (Cstring, Cint, Ptr{UInt8}), "SHORT", 64, buf)
    ccall(Libdl.dlsym(HANDLE, :reset_c), Cvoid, ())
    msg = unsafe_string(pointer(buf))
    m = match(r"SPICE\((.*)\)", msg)
    throw(CError(Symbol(lowercase(m === nothing ? msg : m.captures[1]))))
end

"Number of intervals of the workspace and of the result cell."
const NINTVLS = 20000

"SPICE cell that holds the end points `w` (a valid window)."
function cell(w::AbstractVector{<:Real}, room::Integer=max(length(w), 2))
    c = SpiceDoubleCell(room)
    for i in 1:2:length(w)
        SPICE.wninsd!(c, w[i], w[i+1])
    end
    return c
end

"Run `f(cnfine_cell, result_cell)`, check the error state, return the end points of the result."
function search(f, cnfine::AbstractVector{<:Real})
    cn = cell(cnfine)
    res = SpiceDoubleCell(2 * NINTVLS)
    cd, rd = cn.data, res.data
    GC.@preserve cd rd f(cn.cell, res.cell)
    check()
    return collect(Float64, res)
end

gfdist(target, abcorr, obsrvr, relate, refval, adjust, step, cnfine) = search(cnfine) do cn, res
    ccall(Libdl.dlsym(HANDLE, :gfdist_c), Cvoid,
        (Cstring, Cstring, Cstring, Cstring, SpiceDouble, SpiceDouble, SpiceDouble, SpiceInt,
            Ref{Cell{SpiceDouble}}, Ref{Cell{SpiceDouble}}),
        target, abcorr, obsrvr, relate, refval, adjust, step, NINTVLS, cn, res)
end

gfrr(target, abcorr, obsrvr, relate, refval, adjust, step, cnfine) = search(cnfine) do cn, res
    ccall(Libdl.dlsym(HANDLE, :gfrr_c), Cvoid,
        (Cstring, Cstring, Cstring, Cstring, SpiceDouble, SpiceDouble, SpiceDouble, SpiceInt,
            Ref{Cell{SpiceDouble}}, Ref{Cell{SpiceDouble}}),
        target, abcorr, obsrvr, relate, refval, adjust, step, NINTVLS, cn, res)
end

gfpa(target, illmn, abcorr, obsrvr, relate, refval, adjust, step, cnfine) = search(cnfine) do cn, res
    ccall(Libdl.dlsym(HANDLE, :gfpa_c), Cvoid,
        (Cstring, Cstring, Cstring, Cstring, Cstring, SpiceDouble, SpiceDouble, SpiceDouble,
            SpiceInt, Ref{Cell{SpiceDouble}}, Ref{Cell{SpiceDouble}}),
        target, illmn, abcorr, obsrvr, relate, refval, adjust, step, NINTVLS, cn, res)
end

gfsep(targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, relate, refval, adjust,
    step, cnfine) = search(cnfine) do cn, res
    ccall(Libdl.dlsym(HANDLE, :gfsep_c), Cvoid,
        (Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring,
            SpiceDouble, SpiceDouble, SpiceDouble, SpiceInt, Ref{Cell{SpiceDouble}},
            Ref{Cell{SpiceDouble}}),
        targ1, shape1, frame1, targ2, shape2, frame2, abcorr, obsrvr, relate, refval, adjust,
        step, NINTVLS, cn, res)
end

gfposc(target, frame, abcorr, obsrvr, crdsys, coord, relate, refval, adjust, step, cnfine) =
    search(cnfine) do cn, res
        ccall(Libdl.dlsym(HANDLE, :gfposc_c), Cvoid,
            (Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, SpiceDouble,
                SpiceDouble, SpiceDouble, SpiceInt, Ref{Cell{SpiceDouble}},
                Ref{Cell{SpiceDouble}}),
            target, frame, abcorr, obsrvr, crdsys, coord, relate, refval, adjust, step,
            NINTVLS, cn, res)
    end

gfsubc(target, fixref, method, abcorr, obsrvr, crdsys, coord, relate, refval, adjust, step,
    cnfine) = search(cnfine) do cn, res
    ccall(Libdl.dlsym(HANDLE, :gfsubc_c), Cvoid,
        (Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, SpiceDouble,
            SpiceDouble, SpiceDouble, SpiceInt, Ref{Cell{SpiceDouble}}, Ref{Cell{SpiceDouble}}),
        target, fixref, method, abcorr, obsrvr, crdsys, coord, relate, refval, adjust, step,
        NINTVLS, cn, res)
end

gfsntc(target, fixref, method, abcorr, obsrvr, dref, dvec, crdsys, coord, relate, refval,
    adjust, step, cnfine) = search(cnfine) do cn, res
    ccall(Libdl.dlsym(HANDLE, :gfsntc_c), Cvoid,
        (Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Ptr{SpiceDouble}, Cstring,
            Cstring, Cstring, SpiceDouble, SpiceDouble, SpiceDouble, SpiceInt,
            Ref{Cell{SpiceDouble}}, Ref{Cell{SpiceDouble}}),
        target, fixref, method, abcorr, obsrvr, dref, collect(Float64, dvec), crdsys, coord,
        relate, refval, adjust, step, NINTVLS, cn, res)
end

gfilum(method, angtyp, target, illmn, fixref, abcorr, obsrvr, spoint, relate, refval, adjust,
    step, cnfine) = search(cnfine) do cn, res
    ccall(Libdl.dlsym(HANDLE, :gfilum_c), Cvoid,
        (Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Ptr{SpiceDouble},
            Cstring, SpiceDouble, SpiceDouble, SpiceDouble, SpiceInt, Ref{Cell{SpiceDouble}},
            Ref{Cell{SpiceDouble}}),
        method, angtyp, target, illmn, fixref, abcorr, obsrvr, collect(Float64, spoint),
        relate, refval, adjust, step, NINTVLS, cn, res)
end

gfoclt(occtyp, front, fshape, fframe, back, bshape, bframe, abcorr, obsrvr, step, cnfine) =
    search(cnfine) do cn, res
        ccall(Libdl.dlsym(HANDLE, :gfoclt_c), Cvoid,
            (Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, Cstring,
                SpiceDouble, Ref{Cell{SpiceDouble}}, Ref{Cell{SpiceDouble}}),
            occtyp, front, fshape, fframe, back, bshape, bframe, abcorr, obsrvr, step, cn, res)
    end

gftfov(inst, target, tshape, tframe, abcorr, obsrvr, step, cnfine) = search(cnfine) do cn, res
    ccall(Libdl.dlsym(HANDLE, :gftfov_c), Cvoid,
        (Cstring, Cstring, Cstring, Cstring, Cstring, Cstring, SpiceDouble,
            Ref{Cell{SpiceDouble}}, Ref{Cell{SpiceDouble}}),
        inst, target, tshape, tframe, abcorr, obsrvr, step, cn, res)
end

gfrfov(inst, raydir, rframe, abcorr, obsrvr, step, cnfine) = search(cnfine) do cn, res
    ccall(Libdl.dlsym(HANDLE, :gfrfov_c), Cvoid,
        (Cstring, Ptr{SpiceDouble}, Cstring, Cstring, Cstring, SpiceDouble,
            Ref{Cell{SpiceDouble}}, Ref{Cell{SpiceDouble}}),
        inst, collect(Float64, raydir), rframe, abcorr, obsrvr, step, cn, res)
end

# The user functions of the two calls below. A C callback cannot be a closure on every
# platform, so the Julia functions are kept in global references for the time of a call.
const UDFUNS = Ref{Any}(nothing)
const UDQDEC = Ref{Any}(nothing)
const UDFUNB = Ref{Any}(nothing)

function c_udfuns(et::SpiceDouble, value::Ptr{SpiceDouble})
    unsafe_store!(value, Float64(UDFUNS[](et)))
    return nothing
end
function c_udqdec(::Ptr{Cvoid}, et::SpiceDouble, value::Ptr{SpiceBoolean})
    unsafe_store!(value, SpiceBoolean(UDQDEC[](et) ? 1 : 0))
    return nothing
end
function c_udfunb(::Ptr{Cvoid}, et::SpiceDouble, value::Ptr{SpiceBoolean})
    unsafe_store!(value, SpiceBoolean(UDFUNB[](et) ? 1 : 0))
    return nothing
end

"`gfuds_c` with the Julia functions `udfuns(et)` and `udqdec(et)`."
function gfuds(udfuns, udqdec, relate, refval, adjust, step, cnfine)
    UDFUNS[] = udfuns
    UDQDEC[] = udqdec
    pf = @cfunction(c_udfuns, Cvoid, (SpiceDouble, Ptr{SpiceDouble}))
    pd = @cfunction(c_udqdec, Cvoid, (Ptr{Cvoid}, SpiceDouble, Ptr{SpiceBoolean}))
    return search(cnfine) do cn, res
        ccall(Libdl.dlsym(HANDLE, :gfuds_c), Cvoid,
            (Ptr{Cvoid}, Ptr{Cvoid}, Cstring, SpiceDouble, SpiceDouble, SpiceDouble, SpiceInt,
                Ref{Cell{SpiceDouble}}, Ref{Cell{SpiceDouble}}),
            pf, pd, relate, refval, adjust, step, NINTVLS, cn, res)
    end
end

"`gfudb_c` with the Julia function `udfunb(et)`."
function gfudb(udfunb, step, cnfine)
    UDFUNS[] = et -> 0.0
    UDFUNB[] = udfunb
    pf = @cfunction(c_udfuns, Cvoid, (SpiceDouble, Ptr{SpiceDouble}))
    pb = @cfunction(c_udfunb, Cvoid, (Ptr{Cvoid}, SpiceDouble, Ptr{SpiceBoolean}))
    return search(cnfine) do cn, res
        ccall(Libdl.dlsym(HANDLE, :gfudb_c), Cvoid,
            (Ptr{Cvoid}, Ptr{Cvoid}, SpiceDouble, Ref{Cell{SpiceDouble}}, Ref{Cell{SpiceDouble}}),
            pf, pb, step, cn, res)
    end
end

"`gfstol_c`."
function gfstol(value)
    ccall(Libdl.dlsym(HANDLE, :gfstol_c), Cvoid, (SpiceDouble,), value)
    check()
    return nothing
end

"`uddf_c` with the Julia function `udfuns(et)`."
function uddf(udfuns, x, dx)
    UDFUNS[] = udfuns
    pf = @cfunction(c_udfuns, Cvoid, (SpiceDouble, Ptr{SpiceDouble}))
    out = Ref{SpiceDouble}(0.0)
    ccall(Libdl.dlsym(HANDLE, :uddf_c), Cvoid, (Ptr{Cvoid}, SpiceDouble, SpiceDouble,
        Ref{SpiceDouble}), pf, x, dx, out)
    check()
    return out[]
end

"Short CSPICE error of `f()`, or `:none`."
function kind(f)
    try
        f()
    catch e
        e isa CError || rethrow()
        return e.kind
    end
    return :none
end

end # module GFO

# ---------------------------------------------------------------------------------------
# Run time of the groups: light evaluation of a test body, and the cases of the default run
# ---------------------------------------------------------------------------------------

"""
Module in which `@gf_light` evaluates the body of a group. Its compiler options ask for
the interpreter, so Julia does not infer and optimize a test body of some hundred lines
that runs once. Measured on the fifteen geometry groups: the compilation of the bodies
takes about 40 percent of the compilation time of the groups. The code of NeoSpice, the
calls into CSPICE and the functions that the test files define at top level belong to
other modules and are compiled as usual; allocation and inference checks are therefore
not changed.
"""
module GFLight
Base.Experimental.@compiler_options compile=min optimize=0 infer=false
using Test, Random, LinearAlgebra, Printf, StaticArrays, NeoSpice
end

"""
    gf_light_eval(ex)

Evaluate the expression `ex` in `GFLight`, after every global name of the test module
is imported there. Postcondition: each name that is defined in the test module is
defined in `GFLight`.
"""
function gf_light_eval(ex)
    host = parentmodule(GFLight)
    for n in names(host; all=true, imported=true)
        (startswith(String(n), '#') || n in (:Base, :Core, :eval, :include, :GFLight, nameof(host))) && continue
        (isdefined(host, n) && !isdefined(GFLight, n)) || continue
        Core.eval(GFLight, Expr(:import, Expr(:(:), Expr(:., fullname(host)...), Expr(:., n))))
    end
    # (the imports are visible in the newest world only)
    @assert all(n -> Base.invokelatest(isdefined, GFLight, n), (:group, :SEED, :FULL))
    return Core.eval(GFLight, ex)
end

"Evaluate the expression (a call of `group` with its body) with the interpreter; see `GFLight`."
macro gf_light(ex)
    return :(gf_light_eval($(QuoteNode(ex))))
end

"""
    gf_pick(i, stride) -> Bool

True if case `i` (counted from 1) of a loop over a product of options belongs to this
run: every case in the exhaustive run (`NEOSPICE_TEST_FULL=1`), each `stride`-th case
in the default run. A stride without a common factor with the sizes of the inner loops
spreads the cases of the default run over all the options. Requires `stride >= 1`.
"""
function gf_pick(i::Integer, stride::Integer)
    @assert i >= 1 && stride >= 1
    return FULL || i % stride == 0
end

"Number of cases `1:n` that `gf_pick` selects with `stride`."
gf_npick(n::Integer, stride::Integer) = FULL ? Int(n) : n ÷ stride

"Largest end-point difference in seconds between NeoSpice and CSPICE, by label."
const GF_MAXDIFF = Dict{String,Float64}()
"Number of searches compared, by label."
const GF_COUNT = Dict{String,Int}()

"""
    gf_noise(value, et) -> Float64

Rounding noise of the quantity `value` at the epoch `et`: the largest change of the
quantity over the 16 nearest epochs after `et`. These epochs span about 6e-8 s at 1e8 s
from the reference epoch, so the real change of a geometric quantity is below the
rounding of its inputs, and a root cannot be found more exactly than this noise divided
by the rate of the quantity. Postcondition: the result is not negative.
"""
function gf_noise(value, et::Float64)
    v0 = value(et)
    noise = 0.0
    t = et
    for _ in 1:16
        t = nextfloat(t)
        noise = max(noise, abs(value(t) - v0))
    end
    @assert noise >= 0.0
    return noise
end

"""
    gf_compare(label, neo::Window, ref::Vector{Float64}; atol=2.0e-6, value=nothing,
               vtol=0.0, wide=1.0) -> Bool

True if the two windows have the same number of end points and each end point of `neo`
agrees with the one of CSPICE. Two end points agree if they differ by no more than
`atol` seconds. With a function `value(et)` they also agree if they differ by no more
than `wide` seconds and the quantity has the same value at both, within `vtol`: this is
the case of an ill-conditioned event (a root where the quantity is almost stationary,
or an extremum whose derivative test is rounding noise in an interval), where the two
libraries cannot agree better than the rounding of their states allows. The values
also agree if they differ by no more than 4 times the rounding noise of the quantity
at the end point (`gf_noise`), which is above `vtol` for a corrected quantity far from
the reference epoch. For an absolute extremum give `wide=Inf`: each epoch at which the
quantity has the extreme value to rounding is a correct answer. Records the
largest difference under `label`. Postcondition: `GF_COUNT[label]` is one more than
before.
"""
function gf_compare(label::AbstractString, neo::Window, ref::Vector{Float64};
        atol::Float64=2.0e-6, value=nothing, vtol::Float64=0.0, wide::Float64=1.0)
    GF_COUNT[label] = get(GF_COUNT, label, 0) + 1
    if length(neo) != length(ref)
        @info "window sizes differ" label neo = length(neo) cspice = length(ref)
        return false
    end
    ok = true
    for i in eachindex(ref)
        d = abs(neo[i] - ref[i])
        GF_MAXDIFF[label] = max(get(GF_MAXDIFF, label, 0.0), d)
        d <= atol && continue
        if value !== nothing && d <= wide
            dv = abs(value(neo[i]) - value(ref[i]))
            dv <= vtol && continue
            noise = gf_noise(value, neo[i])
            dv <= 4 * noise && continue
            @info "end points differ, and the values too" label d dv vtol noise
        else
            @info "end points differ" label d
        end
        ok = false
    end
    haskey(GF_MAXDIFF, label) || (GF_MAXDIFF[label] = 0.0)
    return ok
end

"Print the recorded differences whose label starts with `prefix`."
function gf_report(prefix::AbstractString)
    for k in sort!(collect(keys(GF_COUNT)))
        startswith(k, prefix) || continue
        measured(@sprintf("  %-46s searches=%-5d largest end-point difference %.3e s", k,
            GF_COUNT[k], get(GF_MAXDIFF, k, 0.0)))
    end
    return nothing
end

"Family of a relation, for the report: `ineq`, `eq`, `local`, `abs`, `abs+adjust`."
function gf_family(relate::AbstractString, adjust::Real)
    r = uppercase(strip(relate))
    r in ("<", ">") && return "inequality"
    r == "=" && return "equality"
    r in ("LOCMIN", "LOCMAX") && return "local extremum"
    return adjust == 0 ? "absolute extremum" : "absolute extremum, adjust"
end

"Random valid window of `n` intervals in [`a`, `b`]; some intervals are singletons."
function gf_random_window(rng::AbstractRNG, n::Integer, a::Real, b::Real)
    @assert n >= 1 && b > a
    x = sort!([a + (b - a) * rand(rng) for _ in 1:2n])
    for i in 1:2:2n
        rand(rng) < 0.1 && (x[i+1] = x[i])
    end
    w = Window(x)
    @assert wncard(w) >= 1
    return w
end

"The seven relations of the geometry finder."
const GF_RELATIONS = ("<", "=", ">", "LOCMIN", "LOCMAX", "ABSMIN", "ABSMAX")

"The seven relations and the two absolute extrema with a positive adjust value (flag true)."
const GF_VARIANTS = ((("<", false), ("=", false), (">", false), ("LOCMIN", false),
    ("LOCMAX", false), ("ABSMIN", false), ("ABSMAX", false), ("ABSMIN", true), ("ABSMAX", true)))

"The coordinate systems of `gfposc` and their coordinates."
const GF_SYSTEMS = (("RECTANGULAR", ("X", "Y", "Z")),
    ("LATITUDINAL", ("RADIUS", "LONGITUDE", "LATITUDE")),
    ("RA/DEC", ("RANGE", "RIGHT ASCENSION", "DECLINATION")),
    ("SPHERICAL", ("RADIUS", "COLATITUDE", "LONGITUDE")),
    ("CYLINDRICAL", ("RADIUS", "LONGITUDE", "Z")),
    ("GEODETIC", ("LONGITUDE", "LATITUDE", "ALTITUDE")),
    ("PLANETOGRAPHIC", ("LONGITUDE", "LATITUDE", "ALTITUDE")))

"The angular coordinates with a branch cut."
const GF_ANGULAR = (("LATITUDINAL", "LONGITUDE"), ("RA/DEC", "RIGHT ASCENSION"),
    ("SPHERICAL", "LONGITUDE"), ("CYLINDRICAL", "LONGITUDE"), ("GEODETIC", "LONGITUDE"),
    ("PLANETOGRAPHIC", "LONGITUDE"))

"The nine aberration corrections of `spkezr`."
const GF_ABCORRS = ("NONE", "LT", "LT+S", "CN", "CN+S", "XLT", "XLT+S", "XCN", "XCN+S")

"""
    gf_reference(rng, q, cn, rel, withadjust) -> (refval, adjust)

Reference value and adjust value for a search of the quantity `q` (a `GFQuantity` of
NeoSpice) on the window `cn`: for `<`, `=`, `>` a value between the quartiles of the
quantity on the window, so that the search has roots; for an absolute extremum with
`withadjust` a tenth to a third of the range of the quantity.
"""
function gf_reference(rng::AbstractRNG, q, cn::Window, rel::AbstractString, withadjust::Bool)
    ts = range(cn[1], cn[end]; length=60)
    vals = [NeoSpice._gf_value(q, Float64(t)) for t in ts]
    lo, hi = quantile_of(vals, 0.25), quantile_of(vals, 0.75)
    refval = rel in ("<", "=", ">") ? lo + (hi - lo) * rand(rng) : 0.0
    adjust = withadjust ? (0.1 + 0.2 * rand(rng)) * (maximum(vals) - minimum(vals)) : 0.0
    return refval, adjust
end

# User functions and callbacks of the tests behind fixed types. NeoSpice compiles a search
# anew for each type of the user function, of the step function and of each callback
# (0.2 s to 1.2 s for each), and an anonymous function has a type of its own. With the
# types below the tests compile each search once; the values are those of the wrapped
# function, so the results are not changed.

"A scalar user function of `gfuds` and `uddf`; `f(et)` gives a `Float64`."
struct GFScalar <: Function
    f::Any
end
(g::GFScalar)(et) = g.f(et)::Float64

"A user function with a true or false value: `udqdec` of `gfuds`, `udfunb` of `gfudb`."
struct GFState <: Function
    f::Any
end
(g::GFState)(et) = g.f(et)::Bool

"`GFState` that is true where the derivative `d` is negative (the test \"is decreasing\")."
gf_decreasing(d) = GFState(et -> d(et) < 0)

"A step function with the constant value `step`."
struct GFConstStep <: Function
    step::Float64
end
(s::GFConstStep)(et) = s.step

"Interrupt callback: the search stops at the first call if `stop` is true."
struct GFStop <: Function
    stop::Bool
end
(s::GFStop)() = s.stop

"Report callback that keeps the message of the start of each pass in `log`."
struct GFLog <: Function
    log::Vector{String}
end
(l::GFLog)(event, args...) = (event === :init && push!(l.log, args[2]); nothing)

"""
    gf_reference_skip(rng, rel, withadjust)

Draw the random numbers that `gf_reference` draws, without an evaluation of the
quantity. A loop calls it for a case that the default run skips, so that the cases of
the default run are cases of the exhaustive run (same seed, same random numbers).
"""
function gf_reference_skip(rng::AbstractRNG, rel::AbstractString, withadjust::Bool)
    rel in ("<", "=", ">") && rand(rng)
    withadjust && rand(rng)
    return nothing
end

"Kind of the `SpiceError` that `f()` throws, or `:none`."
function gf_kind(f)
    try
        f()
    catch e
        e isa NeoSpice.SpiceError || rethrow()
        return e.kind
    end
    return :none
end

"Shortest of three run times of `f()` in seconds (after one warm-up call)."
function gf_time(f)
    f()
    return minimum(@elapsed(f()) for _ in 1:3)
end
