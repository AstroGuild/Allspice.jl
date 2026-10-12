# CSPICE oracle of the mathematics groups (`test_math_*.jl`).
#
# SPICE.jl replaces many of the vector and matrix calls by `LinearAlgebra` code and does
# not wrap some others, so the oracle calls the C library itself: `CS.name(args...)` is
# the function `name_c` of CSPICE and `CS.name_(args...)` the f2c routine `name_`. A
# CSPICE error becomes a `CS.CError` with the short message as a lower-case symbol,
# which is the `kind` of the `SpiceError` that NeoSpice throws for the same input.
#
# A wrapper is described by the return type and one code for each C argument:
#
#   d, i, b    input double, integer, boolean by value      s    input string
#   rd, ri     input double, integer by reference (f2c)
#   vN         input vector of N doubles                     mN   input N by N matrix
#   fN         input N by N matrix in Fortran order (f2c)
#   od, oi, ob output double, integer, boolean
#   oN         output vector of N doubles                    omN  output N by N matrix
#   ofN        output N by N matrix in Fortran order (f2c)
#   pl, el     input plane (4 doubles), ellipse (9 doubles); opl, oel the outputs
#
# The Julia function takes the inputs in order and returns the C return value (if any)
# followed by the outputs in order, as one value or a tuple.

module CS

using StaticArrays
import SPICE
import Libdl

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

rowmajor(m) = collect(Float64, permutedims(m))
unrow(buf, n) = SMatrix{n,n,Float64}(permutedims(reshape(buf, n, n)))

"Split a code such as `om3` into its kind and its size."
function parsecode(c::Symbol)
    m = match(r"^([a-z]+?)(\d*)$", String(c))
    return Symbol(m.captures[1]), isempty(m.captures[2]) ? 0 : parse(Int, m.captures[2])
end

"Define `CS.jname` as the wrapper of the C function `cname` (see the file header)."
function defc(jname::Symbol, cname::Symbol, ret, codes::Symbol...)
    params = Symbol[]
    pre = Expr[]
    ctypes = Any[]
    cargs = Any[]
    outs = Any[]
    for (k, code) in enumerate(codes)
        kind, n = parsecode(code)
        a = Symbol(:a, k)
        if kind === :d
            push!(params, a); push!(ctypes, Cdouble); push!(cargs, :(Float64($a)))
        elseif kind === :i || kind === :b
            push!(params, a); push!(ctypes, Cint); push!(cargs, :(Cint($a)))
        elseif kind === :s
            push!(params, a); push!(ctypes, Cstring); push!(cargs, :(String($a)))
        elseif kind === :rd
            push!(params, a); push!(ctypes, Ref{Cdouble}); push!(cargs, :(Float64($a)))
        elseif kind === :ri
            push!(params, a); push!(ctypes, Ref{Cint}); push!(cargs, :(Cint($a)))
        elseif kind === :v || kind === :pl || kind === :el
            push!(params, a); push!(ctypes, Ptr{Cdouble}); push!(cargs, :(collect(Float64, $a)))
        elseif kind === :m
            push!(params, a); push!(ctypes, Ptr{Cdouble}); push!(cargs, :(rowmajor($a)))
        elseif kind === :f
            push!(params, a); push!(ctypes, Ptr{Cdouble}); push!(cargs, :(collect(Float64, $a)))
        elseif kind === :od
            push!(pre, :($a = Ref{Cdouble}(NaN))); push!(ctypes, Ref{Cdouble}); push!(cargs, a)
            push!(outs, :($a[]))
        elseif kind === :oi
            push!(pre, :($a = Ref{Cint}(0))); push!(ctypes, Ref{Cint}); push!(cargs, a)
            push!(outs, :(Int($a[])))
        elseif kind === :ob
            push!(pre, :($a = Ref{Cint}(0))); push!(ctypes, Ref{Cint}); push!(cargs, a)
            push!(outs, :($a[] != 0))
        elseif kind === :o || kind === :opl || kind === :oel
            len = kind === :o ? n : kind === :opl ? 4 : 9
            push!(pre, :($a = fill(NaN, $len))); push!(ctypes, Ptr{Cdouble}); push!(cargs, a)
            push!(outs, :(SVector{$len,Float64}($a)))
        elseif kind === :om
            push!(pre, :($a = fill(NaN, $(n * n)))); push!(ctypes, Ptr{Cdouble}); push!(cargs, a)
            push!(outs, :(unrow($a, $n)))
        elseif kind === :of
            push!(pre, :($a = fill(NaN, $(n * n)))); push!(ctypes, Ptr{Cdouble}); push!(cargs, a)
            push!(outs, :(SMatrix{$n,$n,Float64}($a)))
        else
            error("unknown argument code $code")
        end
    end
    call = :(ccall(Libdl.dlsym(HANDLE, $(QuoteNode(cname))), $ret, ($(ctypes...),), $(cargs...)))
    body = if ret === Cvoid
        result = length(outs) == 1 ? outs[1] : Expr(:tuple, outs...)
        quote
            $(pre...)
            $call
            check()
            return $result
        end
    else
        conv = ret === Cint ? :(r) : :(r)
        result = isempty(outs) ? conv : Expr(:tuple, conv, outs...)
        quote
            $(pre...)
            r = $call
            check()
            return $result
        end
    end
    @eval $jname($(params...)) = $body
    return nothing
end

"`defc` for the `_c` function of `name`."
c(name::Symbol, ret, codes::Symbol...) = defc(name, Symbol(name, :_c), ret, codes...)
"`defc` for the f2c routine of `name`; the Julia name ends with an underscore."
f(name::Symbol, ret, codes::Symbol...) = defc(Symbol(name, :_), Symbol(name, :_), ret, codes...)

const V, D, I = Cvoid, Cdouble, Cint

# vectors
c(:vadd, V, :v3, :v3, :o3); c(:vsub, V, :v3, :v3, :o3); c(:vscl, V, :d, :v3, :o3)
c(:vminus, V, :v3, :o3); c(:vdot, D, :v3, :v3); c(:vcrss, V, :v3, :v3, :o3)
c(:vnorm, D, :v3); c(:vhat, V, :v3, :o3); c(:vdist, D, :v3, :v3); c(:vsep, D, :v3, :v3)
c(:vrel, D, :v3, :v3); c(:vproj, V, :v3, :v3, :o3); c(:vperp, V, :v3, :v3, :o3)
c(:vlcom, V, :d, :v3, :d, :v3, :o3); c(:vlcom3, V, :d, :v3, :d, :v3, :d, :v3, :o3)
c(:vrotv, V, :v3, :v3, :d, :o3); c(:vzero, I, :v3); c(:vequ, V, :v3, :o3)
c(:vpack, V, :d, :d, :d, :o3); c(:vupack, V, :v3, :od, :od, :od)
c(:ucrss, V, :v3, :v3, :o3); c(:unorm, V, :v3, :o3, :od); c(:vtmv, D, :v3, :m3, :v3)
f(:vperp, I, :v3, :v3, :o3); f(:vproj, I, :v3, :v3, :o3)
# matrices
c(:mxv, V, :m3, :v3, :o3); c(:mtxv, V, :m3, :v3, :o3); c(:mxm, V, :m3, :m3, :om3)
c(:mtxm, V, :m3, :m3, :om3); c(:mxmt, V, :m3, :m3, :om3); c(:xpose, V, :m3, :om3)
c(:xpose6, V, :m6, :om6); c(:invert, V, :m3, :om3); c(:invort, V, :m3, :om3)
c(:det, D, :m3); c(:trace, D, :m3); c(:isrot, I, :m3, :d, :d); c(:ident, V, :om3)
c(:rotate, V, :d, :i, :om3); c(:rotmat, V, :m3, :d, :i, :om3); c(:rotvec, V, :v3, :d, :i, :o3)
c(:twovec, V, :v3, :i, :v3, :i, :om3); c(:twovxf, V, :v6, :i, :v6, :i, :om6)
c(:dvhat, V, :v6, :o6); c(:dvdot, D, :v6, :v6); c(:dvcrss, V, :v6, :v6, :o6)
c(:ducrss, V, :v6, :v6, :o6); c(:dvnorm, D, :v6); c(:dvsep, D, :v6, :v6)
# rotations
c(:eul2m, V, :d, :d, :d, :i, :i, :i, :om3); c(:m2eul, V, :m3, :i, :i, :i, :od, :od, :od)
c(:eul2xf, V, :v6, :i, :i, :i, :om6); c(:xf2eul, V, :m6, :i, :i, :i, :o6, :ob)
c(:axisar, V, :v3, :d, :om3); c(:raxisa, V, :m3, :o3, :od); c(:m2q, V, :m3, :o4)
c(:q2m, V, :v4, :om3); c(:qxq, V, :v4, :v4, :o4); c(:qdq2av, V, :v4, :v4, :o3)
c(:rav2xf, V, :m3, :v3, :om6); c(:xf2rav, V, :m6, :om3, :o3)

"`sharpr_` works in place on a matrix in Fortran order."
function sharpr_(m)
    buf = collect(Float64, m)
    ccall(Libdl.dlsym(HANDLE, :sharpr_), Cint, (Ptr{Cdouble},), buf)
    check()
    return SMatrix{3,3,Float64}(buf)
end

# general-dimension forms: the dimension arguments are written out
vaddg(a, b) = (o = similar(a, Float64); ccall(Libdl.dlsym(HANDLE, :vaddg_c), Cvoid,
    (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), collect(Float64, a), collect(Float64, b), length(a), o); o)
vsubg(a, b) = (o = similar(a, Float64); ccall(Libdl.dlsym(HANDLE, :vsubg_c), Cvoid,
    (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), collect(Float64, a), collect(Float64, b), length(a), o); o)
vsclg(s, a) = (o = similar(a, Float64); ccall(Libdl.dlsym(HANDLE, :vsclg_c), Cvoid,
    (Cdouble, Ptr{Cdouble}, Cint, Ptr{Cdouble}), s, collect(Float64, a), length(a), o); o)
vminug(a) = (o = similar(a, Float64); ccall(Libdl.dlsym(HANDLE, :vminug_c), Cvoid,
    (Ptr{Cdouble}, Cint, Ptr{Cdouble}), collect(Float64, a), length(a), o); o)
vhatg(a) = (o = similar(a, Float64); ccall(Libdl.dlsym(HANDLE, :vhatg_c), Cvoid,
    (Ptr{Cdouble}, Cint, Ptr{Cdouble}), collect(Float64, a), length(a), o); o)
function unormg(a)
    o = similar(a, Float64)
    mag = Ref{Cdouble}(NaN)
    ccall(Libdl.dlsym(HANDLE, :unormg_c), Cvoid, (Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ref{Cdouble}),
        collect(Float64, a), length(a), o, mag)
    return o, mag[]
end
vdotg(a, b) = ccall(Libdl.dlsym(HANDLE, :vdotg_c), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint),
    collect(Float64, a), collect(Float64, b), length(a))
vnormg(a) = ccall(Libdl.dlsym(HANDLE, :vnormg_c), Cdouble, (Ptr{Cdouble}, Cint),
    collect(Float64, a), length(a))
vdistg(a, b) = ccall(Libdl.dlsym(HANDLE, :vdistg_c), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint),
    collect(Float64, a), collect(Float64, b), length(a))
vsepg(a, b) = ccall(Libdl.dlsym(HANDLE, :vsepg_c), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint),
    collect(Float64, a), collect(Float64, b), length(a))
vsepg_(a, b) = ccall(Libdl.dlsym(HANDLE, :vsepg_), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}),
    collect(Float64, a), collect(Float64, b), Cint(length(a)))
vrelg(a, b) = ccall(Libdl.dlsym(HANDLE, :vrelg_c), Cdouble, (Ptr{Cdouble}, Ptr{Cdouble}, Cint),
    collect(Float64, a), collect(Float64, b), length(a))
vzerog(a) = ccall(Libdl.dlsym(HANDLE, :vzerog_c), Cint, (Ptr{Cdouble}, Cint),
    collect(Float64, a), length(a)) != 0
vlcomg(x, a, y, b) = (o = similar(a, Float64); ccall(Libdl.dlsym(HANDLE, :vlcomg_c), Cvoid,
    (Cint, Cdouble, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), length(a), x,
    collect(Float64, a), y, collect(Float64, b), o); o)
vprojg(a, b) = (o = similar(a, Float64); ccall(Libdl.dlsym(HANDLE, :vprojg_c), Cvoid,
    (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), collect(Float64, a), collect(Float64, b), length(a), o); o)
vtmvg(a, m, b) = ccall(Libdl.dlsym(HANDLE, :vtmvg_c), Cdouble,
    (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint), collect(Float64, a), rowmajor(m),
    collect(Float64, b), size(m, 1), size(m, 2))
function mxvg(m, v)
    o = zeros(size(m, 1))
    ccall(Libdl.dlsym(HANDLE, :mxvg_c), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}),
        rowmajor(m), collect(Float64, v), size(m, 1), size(m, 2), o)
    return o
end
function mtxvg(m, v)
    o = zeros(size(m, 2))
    ccall(Libdl.dlsym(HANDLE, :mtxvg_c), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}),
        rowmajor(m), collect(Float64, v), size(m, 2), size(m, 1), o)
    return o
end
"Matrix of size `r` by `c` from a row-major buffer."
unrowg(buf, r, c) = collect(permutedims(reshape(buf, c, r)))
function mxmg(a, b)
    o = zeros(size(a, 1) * size(b, 2))
    ccall(Libdl.dlsym(HANDLE, :mxmg_c), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}),
        rowmajor(a), rowmajor(b), size(a, 1), size(a, 2), size(b, 2), o)
    return unrowg(o, size(a, 1), size(b, 2))
end
function mtxmg(a, b)
    o = zeros(size(a, 2) * size(b, 2))
    ccall(Libdl.dlsym(HANDLE, :mtxmg_c), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}),
        rowmajor(a), rowmajor(b), size(a, 2), size(a, 1), size(b, 2), o)
    return unrowg(o, size(a, 2), size(b, 2))
end
function mxmtg(a, b)
    o = zeros(size(a, 1) * size(b, 1))
    ccall(Libdl.dlsym(HANDLE, :mxmtg_c), Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cint, Ptr{Cdouble}),
        rowmajor(a), rowmajor(b), size(a, 1), size(a, 2), size(b, 1), o)
    return unrowg(o, size(a, 1), size(b, 1))
end
function xposeg(a)
    o = zeros(length(a))
    ccall(Libdl.dlsym(HANDLE, :xposeg_c), Cvoid, (Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}),
        rowmajor(a), size(a, 1), size(a, 2), o)
    return unrowg(o, size(a, 2), size(a, 1))
end

# coordinates
for n in (:reclat, :recrad, :recsph, :reccyl)
    c(n, V, :v3, :od, :od, :od)
end
for n in (:latrec, :radrec, :sphrec, :cylrec)
    c(n, V, :d, :d, :d, :o3)
end
for n in (:cyllat, :latcyl, :cylsph, :sphcyl, :sphlat, :latsph)
    c(n, V, :d, :d, :d, :od, :od, :od)
end
c(:recgeo, V, :v3, :d, :d, :od, :od, :od); c(:georec, V, :d, :d, :d, :d, :d, :o3)
c(:recpgr, V, :s, :v3, :d, :d, :od, :od, :od); c(:pgrrec, V, :s, :d, :d, :d, :d, :d, :o3)
for n in (:drdlat, :dlatdr, :drdsph, :dsphdr, :drdcyl, :dcyldr)
    c(n, V, :d, :d, :d, :om3)
end
c(:drdgeo, V, :d, :d, :d, :d, :d, :om3); c(:dgeodr, V, :d, :d, :d, :d, :d, :om3)
c(:drdpgr, V, :s, :d, :d, :d, :d, :d, :om3); c(:dpgrdr, V, :s, :d, :d, :d, :d, :d, :om3)
c(:xfmsta, V, :v6, :s, :s, :s, :o6); c(:convrt, V, :d, :s, :s, :od)
# ellipsoids, planes, ellipses
c(:nearpt, V, :v3, :d, :d, :d, :o3, :od); c(:dnearp, V, :v6, :d, :d, :d, :o6, :o2, :ob)
c(:surfpt, V, :v3, :v3, :d, :d, :d, :o3, :ob); c(:surfnm, V, :d, :d, :d, :v3, :o3)
c(:surfpv, V, :v6, :v6, :d, :d, :d, :o6, :ob); c(:npedln, V, :d, :d, :d, :v3, :v3, :o3, :od)
c(:npelpt, V, :v3, :el, :o3, :od); c(:nplnpt, V, :v3, :v3, :v3, :o3, :od)
c(:edlimb, V, :d, :d, :d, :v3, :oel); c(:inedpl, V, :d, :d, :d, :pl, :oel, :ob)
c(:inelpl, V, :el, :pl, :oi, :o3, :o3); c(:inrypl, V, :v3, :v3, :pl, :oi, :o3)
c(:cgv2el, V, :v3, :v3, :v3, :oel); c(:el2cgv, V, :el, :o3, :o3, :o3); c(:saelgv, V, :v3, :v3, :o3, :o3)
c(:nvc2pl, V, :v3, :d, :opl); c(:nvp2pl, V, :v3, :v3, :opl); c(:psv2pl, V, :v3, :v3, :v3, :opl)
c(:pl2nvc, V, :pl, :o3, :od); c(:pl2nvp, V, :pl, :o3, :o3); c(:pl2psv, V, :pl, :o3, :o3, :o3)
c(:pjelpl, V, :el, :pl, :oel); c(:ednmpt, V, :d, :d, :d, :v3, :o3); c(:edpnt, V, :v3, :d, :d, :d, :o3)
c(:vprjp, V, :v3, :pl, :o3); c(:vprjpi, V, :v3, :pl, :pl, :o3, :ob)

"`edterm_c`: the epoch, the observer position and the terminator points as rows."
function edterm(trmtyp, source, target, et, fixref, abcorr, obsrvr, npts)
    trgepc = Ref{Cdouble}(NaN)
    obspos = fill(NaN, 3)
    trmpts = fill(NaN, 3 * max(npts, 1))
    ccall(Libdl.dlsym(HANDLE, :edterm_c), Cvoid, (Cstring, Cstring, Cstring, Cdouble, Cstring,
            Cstring, Cstring, Cint, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        trmtyp, source, target, et, fixref, abcorr, obsrvr, npts, trgepc, obspos, trmpts)
    check()
    return trgepc[], SVector{3,Float64}(obspos), [SVector{3,Float64}(trmpts[3k-2:3k]) for k in 1:npts]
end

# conics
c(:conics, V, :v8, :d, :o6); c(:oscelt, V, :v6, :d, :d, :o8); c(:oscltx, V, :v6, :d, :d, :o20)
c(:prop2b, V, :d, :v6, :d, :o6); f(:stmp03, I, :rd, :od, :od, :od, :od)

end # module CS

# --- comparison helpers ----------------------------------------------------------------

"All numbers of a result (number, Bool, array, tuple, `nothing`) as one `Float64` vector."
flat(x::Real) = Float64[x]
flat(::Nothing) = Float64[]
flat(x::AbstractArray{<:Real}) = Float64[Float64(v) for v in vec(collect(x))]
flat(x::AbstractArray) = reduce(vcat, map(flat, x); init=Float64[])
flat(x::Tuple) = reduce(vcat, map(flat, x); init=Float64[])
flat(x::SpicePlane) = Float64[x.normal..., x.constant]
flat(x::SpiceEllipse) = Float64[x.center..., x.semimajor..., x.semiminor...]

"""
    reldiff(a, b) -> Float64

Largest absolute difference of two results divided by the largest absolute value of
`b` (not divided if that is zero). `Inf` if the shapes differ or a NaN is on one side
only.
"""
function reldiff(a, b)
    x = flat(a)
    y = flat(b)
    length(x) == length(y) || return Inf
    isempty(x) && return 0.0
    any(isnan.(x) .!= isnan.(y)) && return Inf
    d = 0.0
    s = 0.0
    for (p, q) in zip(x, y)
        isnan(p) && continue
        d = max(d, p == q ? 0.0 : abs(p - q))
        isfinite(q) && (s = max(s, abs(q)))
    end
    return s == 0.0 ? d : d / s
end

"Result of `f()`, or the error kind as a symbol for a NeoSpice or CSPICE error."
function outcome(f)
    try
        return f()
    catch e
        e isa NeoSpice.SpiceError && return e.kind
        e isa CS.CError && return e.kind
        rethrow()
    end
end

"""
    Tally

Largest difference against CSPICE and the count of acceptance mismatches (one side
throws, or the two sides throw different kinds) for each function name.
"""
struct Tally
    maxdiff::Dict{Symbol,Float64}
    count::Dict{Symbol,Int}
    mismatch::Dict{Symbol,Int}
    order::Vector{Symbol}
end
Tally() = Tally(Dict{Symbol,Float64}(), Dict{Symbol,Int}(), Dict{Symbol,Int}(), Symbol[])

"""
    versus!(t, name, neo, c) -> Float64

Run the two closures, record the difference of their results under `name` and return
it. Equal error kinds count as difference 0; an error on one side only is a mismatch.
"""
function versus!(t::Tally, name::Symbol, neo, c)
    a = outcome(neo)
    b = outcome(c)
    haskey(t.count, name) || (push!(t.order, name); t.maxdiff[name] = 0.0; t.count[name] = 0;
        t.mismatch[name] = 0)
    t.count[name] += 1
    if a isa Symbol || b isa Symbol
        a === b && return 0.0
        t.mismatch[name] += 1
        t.mismatch[name] <= 3 && @info "acceptance mismatch" name neo = a cspice = b
        return Inf
    end
    d = reldiff(a, b)
    t.maxdiff[name] = max(t.maxdiff[name], d)
    return d
end

"Print the table of `t` and test each function against `tol` (default) or `tols[name]`."
function report!(t::Tally, family::AbstractString; tol::Float64=0.0, tols=Dict{Symbol,Float64}())
    worst = 0.0
    for name in t.order
        limit = get(tols, name, tol)
        measured(@sprintf("  %-10s %-10s n=%-6d maxdiff=%.3e (limit %.1e) mismatches=%d", family,
            name, t.count[name], t.maxdiff[name], limit, t.mismatch[name]))
        @test t.mismatch[name] == 0
        @test t.maxdiff[name] <= limit
        worst = max(worst, t.maxdiff[name])
    end
    measured(@sprintf("  %-10s largest difference against CSPICE: %.3e", family, worst))
    return worst
end

# --- random inputs ---------------------------------------------------------------------

"Random magnitude, log-uniform in [1e-10, 1e10]."
rmag(rng) = 10.0^(20 * rand(rng) - 10)

const SPECIAL_VECTORS = (SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0),
    SVector(0.0, 0.0, 1.0), SVector(-1.0, 0.0, 0.0), SVector(0.0, 0.0, -1.0),
    SVector(1.0, 1.0, 0.0), SVector(1e-300, 0.0, 0.0), SVector(1e300, 1e300, 1e300))

"Random vector with a random magnitude; a special vector with probability `special`."
function rvec(rng; special::Float64=0.1)
    rand(rng) < special && return rand(rng, SPECIAL_VECTORS)
    return rmag(rng) * SVector(randn(rng), randn(rng), randn(rng))
end

"Pair of vectors; with probability `special` they are equal, parallel or opposite."
function rpair(rng; special::Float64=0.15)
    a = rvec(rng)
    u = rand(rng)
    u < special / 3 && return a, a
    u < 2special / 3 && return a, 3.0 * a
    u < special && return a, -0.5 * a
    return a, rvec(rng)
end

"Random rotation matrix (exact up to rounding) from three Euler angles."
function rrot(rng)
    a = 2pi * rand(rng, 3) .- pi
    return CS.eul2m(a[1], a[2], a[3], rand(rng, 1:3), rand(rng, 1:3), rand(rng, 1:3))
end

"Random state: a position and a velocity with independent magnitudes."
rstate(rng; special::Float64=0.1) = vcat(rvec(rng; special), rvec(rng; special))

"Allocation of the second call of `f(args...)`; the first call compiles."
function alloc_of(f::F, args::Vararg{Any,N}) where {F,N}
    f(args...)
    return @allocated f(args...)
end
