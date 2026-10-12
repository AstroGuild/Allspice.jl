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

"Largest end-point difference in seconds between NeoSpice and CSPICE, by label."
const GF_MAXDIFF = Dict{String,Float64}()
"Number of searches compared, by label."
const GF_COUNT = Dict{String,Int}()

"""
    gf_compare(label, neo::Window, ref::Vector{Float64}; atol=2.0e-6, value=nothing,
               vtol=0.0, wide=1.0) -> Bool

True if the two windows have the same number of end points and each end point of `neo`
agrees with the one of CSPICE. Two end points agree if they differ by no more than
`atol` seconds. With a function `value(et)` they also agree if they differ by no more
than `wide` seconds and the quantity has the same value at both, within `vtol`: this is
the case of an ill-conditioned event (a root where the quantity is almost stationary,
or an extremum whose derivative test is rounding noise in an interval), where the two
libraries cannot agree better than the rounding of their states allows. Records the
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
            @info "end points differ, and the values too" label d dv vtol
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
