# Geometry finder, solver core: the step-and-bracket search for the transitions of a
# binary state (`zzgfsolv`/`zzgfsolvx` of CSPICE N0067), the relation driver that turns
# a scalar quantity and a relation into a window (`zzgfrel`/`zzgfrelx`), the window
# sifter `zzgfwsts`, the refinement `gfrefn` and the numerical derivative `uddf`/`uddc`.
#
# The core knows nothing about geometry. A search is described by two functions of the
# epoch, the value of the quantity and the test "is the quantity decreasing", so a new
# quantity needs no change here (`quantities.jl` builds these functions).
#
# CSPICE keeps the step, the tolerance, the reference value and the progress reporter
# in saved variables. Here they are fields of a `_GFRun` value that the caller builds
# for one search, so a search is a pure function of its arguments and many searches can
# run at the same time.

"Default convergence tolerance of the geometry finder in seconds (`CNVTOL` of CSPICE)."
const _GF_CNVTOL = 1.0e-6

"Largest number of refinement steps for one transition (`MXLOOP` of `zzgfsolv`)."
const _GF_MXLOOP = 1000

"Expansion of the confinement window for an equality search, in seconds (`ADDWIN`)."
const _GF_ADDWIN = 0.5

"""
    gfrefn(t1, t2, s1=false, s2=true) -> Float64

Return the midpoint of the bracket [`t1`, `t2`], as `gfrefn_c`. This is the default
refinement function of the geometry finder.

The result is `t1 * 0.5 + t2 * 0.5`, clamped between `t1` and `t2`; the two ends can be in
either order. The states `s1` and `s2` at the two ends are not used.

See also [`gfevnt`](@ref), [`gfocce`](@ref), [`gffove`](@ref).
"""
@inline function gfrefn(t1::Real, t2::Real, s1::Bool=false, s2::Bool=true)
    a = Float64(t1)
    b = Float64(t2)
    return _gf_brcktd(a * 0.5 + b * 0.5, a, b)
end

"`brcktd`: `x` bracketed between the two ends `a` and `b`, in either order."
@inline _gf_brcktd(x::Float64, a::Float64, b::Float64) =
    a < b ? max(a, min(b, x)) : max(b, min(a, x))

"""
    _GFRun(step, refine, bail, report, tol, threads)

Controls of one search. `step` is the search step in seconds, a `Float64` or a function
of the epoch (`udstep`); `refine` is the refinement function `(t1, t2, s1, s2) -> t`
(`udrefn`); `bail` is `nothing` or a function without arguments that returns `true` to
stop the search (`udbail`); `report` is `nothing` or the progress callback (`udrepi`,
`udrepu`, `udrepf`); `tol` is the convergence tolerance in seconds; `threads` asks for
one task for each interval that is searched.

Invariant: `tol > 0`, and a constant `step` is positive.
"""
struct _GFRun{S,R,B,P}
    step::S
    refine::R
    bail::B
    report::P
    tol::Float64
    threads::Bool
end

"""
    _gf_run(step, tol; refine=gfrefn, bail=nothing, report=nothing, threads=false) -> _GFRun

Check the step and the tolerance as `gfsstp` and `zzgfsolv` do and build the controls.
Throws `SpiceError(:invalidstep)` for a constant step that is not positive and
`SpiceError(:invalidtolerance)` for a tolerance that is not positive.
"""
function _gf_run(step, tol::Real; refine=gfrefn, bail=nothing, report=nothing,
        threads::Bool=false)
    st = step isa Real ? Float64(step) : step
    st isa Float64 && !(st > 0.0) && throw(SpiceError(:invalidstep,
        "Step has value $st; step size must be positive."))
    tol > 0 || throw(SpiceError(:invalidtolerance,
        "Tolerance must be positive but was $tol."))
    return _GFRun(st, refine, bail, report, Float64(tol), threads)
end

_gf_withtol(run::_GFRun, tol::Float64) =
    _GFRun(run.step, run.refine, run.bail, run.report, tol, run.threads)

@inline _gf_step(run::_GFRun{Float64}, et::Float64) = run.step
@inline function _gf_step(run::_GFRun, et::Float64)
    st = Float64(run.step(et))::Float64
    st > 0.0 || throw(SpiceError(:invalidstep,
        "Step has value $st at epoch $et; step size must be positive."))
    return st
end

@inline function _gf_bail(run::_GFRun)
    run.bail === nothing && return nothing
    (run.bail()::Bool) && throw(InterruptException())
    return nothing
end

# Progress report: the three calls of CSPICE become one callback with a first argument
# that names the event. The test on `nothing` is decided by the type of `run`.
@inline _gf_repi(run::_GFRun, w::Window, pre::String, suf::String) =
    (run.report === nothing || run.report(:init, w, pre, suf); nothing)
@inline _gf_repu(run::_GFRun, a::Float64, b::Float64, t::Float64) =
    (run.report === nothing || run.report(:update, a, b, t); nothing)
@inline _gf_repf(run::_GFRun) = (run.report === nothing || run.report(:final); nothing)

"""
    _gf_solve!(result, cond, run, cstep, step, start, finish) -> result

The search of `zzgfsolvx`: insert into `result` the intervals of [`start`, `finish`] on
which the binary state `cond(et)::Bool` is true. The state is sampled from `start` in
steps (`step` if `cstep`, else the step of `run`), the last sample is at `finish`, and
each change of state between two samples is refined with `run.refine` until the bracket
is not longer than `run.tol`; the transition is the midpoint of the last bracket.

Contract: precondition `start <= finish`, else `SpiceError(:badtimecase)`. The sequence
of epochs at which `cond` is evaluated is the sequence of CSPICE, but for the second
evaluation at `finish` that CSPICE makes when the last sample changes nothing.
Postcondition: `result` is a valid window that holds its old value and the new
intervals, which lie in [`start`, `finish`]. Throws `SpiceError(:noconverg)` after
$(_GF_MXLOOP) refinement steps of one transition.
"""
function _gf_solve!(result::Window, cond::F, run::_GFRun, cstep::Bool, step::Float64,
        start::Float64, finish::Float64) where {F}
    tol = run.tol
    start > finish && throw(SpiceError(:badtimecase,
        "Bad time interval result, START > FINISH."))
    _gf_repu(run, start, finish, start)
    curtim = start
    curste = cond(curtim)::Bool
    instat = curste
    begtim = curtim
    timest = step
    svdtim = curtim
    savst = curste
    # invariant: every transition before `svdtim` is recorded; `instat` is the state at
    # `svdtim`, and `begtim` is the start of the open interval if `instat`
    while svdtim < finish
        cstep || (timest = _gf_step(run, curtim))
        curtim = min(curtim + timest, finish)
        curste = cond(curtim)::Bool
        # invariant: the state is `savst` at every sample from the last transition to
        # `svdtim`; `curtim` is the sample after `svdtim`
        while savst == curste && svdtim < finish
            _gf_bail(run)
            _gf_repu(run, start, finish, svdtim)
            svdtim = curtim
            savst = curste
            cstep || (timest = _gf_step(run, curtim))
            nxt = min(curtim + timest, finish)
            # CSPICE evaluates the state at `finish` a second time here; the state is a
            # function of the epoch, so the value in hand is used.
            nxt == curtim || (curste = cond(nxt)::Bool)
            curtim = nxt
        end
        if savst != curste
            state1 = savst
            t1 = svdtim
            t2 = curtim
            l1 = savst
            l2 = curste
            t1 > t2 && throw(SpiceError(:badtimecase, "Bad time interval result, T1 > T2."))
            svdtim = curtim
            savst = curste
            prvdif = floatmax(Float64)
            diff = abs(t2 - t1)
            nloop = 0
            # invariant: the state is `state1` at `t1` and the other state at `t2`, and
            # the bracket became shorter in every step
            while diff > tol && diff < prvdif
                nloop += 1
                nloop >= _GF_MXLOOP && throw(SpiceError(:noconverg,
                    "Loop run exceeds maximum loop count. Unable to converge to TOL " *
                    "value $tol within MXLOOP value $(_GF_MXLOOP) iterations."))
                _gf_bail(run)
                t = _gf_brcktd(Float64(run.refine(t1, t2, l1, l2))::Float64, t1, t2)
                if t == t1
                    t2 = t
                elseif t == t2
                    t1 = t
                else
                    s = cond(t)::Bool
                    if s == state1
                        t1 = t
                        l1 = s
                    else
                        t2 = t
                        l2 = s
                    end
                end
                prvdif = diff
                diff = abs(t2 - t1)
            end
            trnstn = _gf_brcktd((t1 + t2) * 0.5, t1, t2)
            if instat
                wninsd!(result, begtim, trnstn)
            else
                begtim = trnstn
            end
            instat = !instat
        end
    end
    instat && wninsd!(result, begtim, finish)
    _gf_repu(run, start, finish, finish)
    return result
end

"""
    _gf_scan!(result, cond, w, run, cstep) -> result

Run `_gf_solve!` on each interval of the window `w` and collect the intervals on which
`cond` is true in `result`. With `cstep` the step of each search is the length of its
interval (the second pass of `zzgfrelx`, on intervals where the quantity is monotone).

Contract: with `run.threads` each interval of `w` is searched in its own task and the
results are inserted in the order of `w`; the searches of different intervals share no
state, so the window is the one of the serial run, bit for bit. `cond` and the
callbacks of `run` must then be safe to call from several threads.
"""
function _gf_scan!(result::Window, cond::F, w::Window, run::_GFRun, cstep::Bool) where {F}
    n = wncard(w)
    if run.threads && n > 1 && Threads.nthreads() > 1
        parts = Vector{Window}(undef, n)
        Threads.@threads for i in 1:n
            a = w.data[2i-1]
            b = w.data[2i]
            parts[i] = _gf_solve!(Window(), cond, run, cstep, b - a, a, b)
        end
        for p in parts                   # invariant: result holds the parts before p
            d = p.data
            for k in 1:2:length(d)
                wninsd!(result, d[k], d[k+1])
            end
        end
        return result
    end
    for i in 1:n                         # invariant: intervals 1 .. i-1 of w are searched
        a = w.data[2i-1]
        b = w.data[2i]
        _gf_solve!(result, cond, run, cstep, b - a, a, b)
    end
    return result
end

"""
    _gf_wnexpd(w, left, right) -> Window

Copy of `w` with each interval expanded by `left` at its left end and by `right` at its
right end, and overlapping intervals merged, as `wnexpd_c` (the amounts are not
negative here).
"""
function _gf_wnexpd(w::Window, left::Float64, right::Float64)
    d = copy(w.data)
    card = length(d)
    card == 0 && return Window()
    for i in 1:2:card                    # invariant: intervals before i are expanded
        d[i] -= left
        d[i+1] += right
    end
    i = 2
    j = 2
    while j < card                       # invariant: d[1:i-1] is final but for its last end
        if d[j] < d[j+1]
            d[i] = d[j]
            d[i+1] = d[j+1]
            i += 2
        end
        j += 2
    end
    d[i] = d[j]
    resize!(d, i)
    return _window_of(d)
end

"Window of the singleton intervals at the left (`side == 1`) or right (`side == 2`) ends of `w` (`wnextd_c`)."
function _gf_wnextd(w::Window, side::Int)
    d = Float64[]
    for i in 1:wncard(w)                 # invariant: d holds the ends of intervals before i
        x = w.data[2i-2+side]
        push!(d, x, x)
    end
    return _window_of(d)
end

"""
    _gf_wsts(w1, w2) -> Window

The intervals of `w1` that lie in the interior of an interval of `w2`, as `zzgfwsts`
with the inclusion flag `"()"`: both ends of the interval of `w1` are strictly inside.
"""
function _gf_wsts(w1::Window, w2::Window)
    a, b = w1.data, w2.data
    na, nb = length(a), length(b)
    out = Float64[]
    p1 = 1
    p2 = 1
    while p1 < na && p2 < nb             # invariant: intervals before p1 and p2 are sifted
        if a[p1+1] < b[p2]
            p1 += 2
        elseif b[p2+1] < a[p1]
            p2 += 2
        else
            (a[p1] > b[p2] && a[p1+1] < b[p2+1]) && push!(out, a[p1], a[p1+1])
            if a[p1+1] < b[p2+1]
                p1 += 2
            elseif b[p2+1] < a[p1+1]
                p2 += 2
            else
                p1 += 2
                p2 += 2
            end
        end
    end
    return _window_of(out)
end

# Relations, as small integers.
const _GF_LT = 1
const _GF_EQ = 2
const _GF_GT = 3
const _GF_LOCMIN = 4
const _GF_ABSMIN = 5
const _GF_LOCMAX = 6
const _GF_ABSMAX = 7

_gf_relname(rel::Int) = ("<", "=", ">", "LOCMIN", "ABSMIN", "LOCMAX", "ABSMAX")[rel]

function _gf_relcode(relate::AbstractString)
    s = uppercase(strip(relate, ' '))
    s == "<" && return _GF_LT
    s == "=" && return _GF_EQ
    s == ">" && return _GF_GT
    s == "LOCMIN" && return _GF_LOCMIN
    s == "ABSMIN" && return _GF_ABSMIN
    s == "LOCMAX" && return _GF_LOCMAX
    s == "ABSMAX" && return _GF_ABSMAX
    throw(SpiceError(:notrecognized, "The comparison operator, $relate is not " *
        "recognized.  Supported quantities are: <, =, >, LOCMIN, ABSMIN, LOCMAX, ABSMAX."))
end

"""
    _gf_checkadjust(rel, adjust)

The checks of `zzgfrelx` on the adjustment value: `SpiceError(:valueoutofrange)` if it
is negative, `SpiceError(:invalidvalue)` if it is not zero and the relation is not an
absolute extremum.
"""
function _gf_checkadjust(rel::Int, adjust::Real)
    adjust < 0 && throw(SpiceError(:valueoutofrange,
        "ADJUST was $adjust; must be non-negative."))
    (rel == _GF_ABSMIN || rel == _GF_ABSMAX || adjust == 0) || throw(SpiceError(:invalidvalue,
        "ADJUST should have value zero for all comparison operators except ABSMAX and ABSMIN"))
    return nothing
end

"The binary state \"the quantity is less than the reference value\" (`zzgfudlt`)."
struct _GFLess{V}
    value::V
    ref::Float64
end
@inline (c::_GFLess)(et::Float64) = c.value(et) < c.ref

"""
    _gf_relx(value, decreasing, relate, refval, adjust, cnfine, run, prefix, suffix) -> Window

The relation driver `zzgfrelx`: the window, inside the confinement window `cnfine`, on
which the scalar quantity `value(et)::Float64` satisfies the relation `relate` (`"<"`,
`"="`, `">"`, `"LOCMIN"`, `"LOCMAX"`, `"ABSMIN"`, `"ABSMAX"`) with the reference value
`refval`; `decreasing(et)::Bool` says if the quantity decreases at `et`.

The first pass finds the window on which the quantity decreases, with the step of
`run`; its complement is the window on which the quantity does not decrease, and the
ends of the two are the local extrema. The second pass searches each interval of the
two windows, where the quantity is monotone, for the single change of the state
`value(et) < reference`, with a step equal to the length of the interval. For `"="` the
confinement window is expanded by $(_GF_ADDWIN) s at each end before the search, so
that a root at an end of `cnfine` is found, and the result is the set of ends of the
"less than" window that lie in `cnfine`. An absolute extremum is the most extreme
value at the ends of the monotone intervals in `cnfine`; with `adjust > 0` the result
is the window on which the quantity is within `adjust` of that value.

Contract: `adjust >= 0`, and `adjust == 0` unless the relation is an absolute
extremum; the errors are those of `_gf_relcode` and `_gf_checkadjust`. An empty
`cnfine` gives an empty window. `cnfine` is not modified. The result is a new valid
window that is a subset of `cnfine`. `prefix` holds the progress messages of the two
passes.
"""
function _gf_relx(value::V, decreasing::D, relate::AbstractString, refval::Real,
        adjust::Real, cnfine::Window, run::_GFRun,
        prefix::NTuple{2,String}=("", ""), suffix::String="") where {V,D}
    rel = _gf_relcode(relate)
    _gf_checkadjust(rel, adjust)
    card(cnfine) == 0 && return Window()
    # CSPICE compares the relation as it was given with "=", so an equality relation
    # with a leading blank gets no expansion; this is kept.
    addw = rstrip(relate, ' ') == "=" ? _GF_ADDWIN : 0.0
    confin = _gf_wnexpd(cnfine, addw, addw)
    refer2 = Float64(refval)

    # Pass 1: the window on which the quantity decreases.
    _gf_repi(run, confin, prefix[1], suffix)
    decres = _gf_scan!(Window(), decreasing, confin, run, false)
    _gf_repf(run)

    # The right ends of `decres` are local minima and the left ends local maxima; the
    # ones in the interior of the confinement window are kept.
    rel == _GF_LOCMIN && return _gf_wsts(_gf_wnextd(decres, 2), cnfine)
    rel == _GF_LOCMAX && return _gf_wsts(_gf_wnextd(decres, 1), cnfine)

    if rel == _GF_ABSMIN || rel == _GF_ABSMAX
        dec = wnintd(cnfine, decres)
        inc = wndifd(cnfine, dec)
        result = Window()
        extrem = 0.0
        need = true
        # On a decreasing interval the maximum is at the left end and the minimum at
        # the right end; on an increasing interval the ends change roles.
        minat, maxat = 2, 1
        for w in (dec, inc)
            want = rel == _GF_ABSMIN ? minat : maxat
            # invariant: `extrem` is the most extreme value at the candidate epochs
            # seen so far, and `result` holds the last epoch that reached it
            for i in 1:wncard(w)
                t = w.data[2i-2+want]
                v = Float64(value(t))::Float64
                if need
                    need = false
                    extrem = v
                end
                if rel == _GF_ABSMIN
                    if adjust == 0 && v <= extrem
                        empty!(result)
                        wninsd!(result, t, t)
                    end
                    extrem = extrem <= v ? extrem : v
                else
                    if adjust == 0 && v >= extrem
                        empty!(result)
                        wninsd!(result, t, t)
                    end
                    extrem = extrem >= v ? extrem : v
                end
            end
            minat, maxat = maxat, minat
        end
        adjust == 0 && return result
        refer2 = rel == _GF_ABSMIN ? extrem + adjust : extrem - adjust
    end

    incres = wndifd(confin, decres)

    # Pass 2: the window on which the quantity is less than the reference value.
    _gf_repi(run, confin, prefix[2], suffix)
    cond = _GFLess(value, refer2)
    less = Window()
    _gf_scan!(less, cond, decres, run, true)
    _gf_scan!(less, cond, incres, run, true)
    _gf_repf(run)

    (rel == _GF_LT || rel == _GF_ABSMIN) && return wnintd(cnfine, less)
    (rel == _GF_GT || rel == _GF_ABSMAX) && return wndifd(cnfine, less)
    # Equality: the ends of `less` are the roots and the ends of the expanded
    # confinement window; the ones in `cnfine` are the roots.
    temp = Window()
    for x in less.data                   # invariant: temp holds the ends before x
        wninsd!(temp, x, x)
    end
    return wnintd(cnfine, temp)
end

"""
    uddf(udfunc, x, dx) -> Float64

Return the numerical derivative of the scalar function `udfunc` at `x`, as `uddf_c`.

The result is the central difference `(udfunc(x + dx) - udfunc(x - dx)) / (2 dx)`.
`udfunc` is called once at `x - dx` and then once at `x + dx`. Unlike CSPICE, `udfunc`
takes the abscissa alone.

Throws `SpiceError(:dividebyzero)` if `dx == 0`.

# Examples
```jldoctest
julia> uddf(x -> x^2, 3.0, 0.5)
6.0
```

See also [`uddc`](@ref), [`gfuds`](@ref).
"""
function uddf(udfunc::F, x::Real, dx::Real) where {F}
    h = Float64(dx)
    f0 = Float64(udfunc(Float64(x) - h))
    f2 = Float64(udfunc(Float64(x) + h))
    h == 0.0 && throw(SpiceError(:dividebyzero,
        "Delta abscissa value is zero; a non-zero value is required."))
    return (0.5 / h) * f2 + (-0.5 / h) * f0
end

"""
    uddc(udfunc, x, dx) -> Bool

Return `true` if the numerical derivative `uddf(udfunc, x, dx)` is negative, as `uddc_c`.

This is the "is decreasing" test that [`gfuds`](@ref) needs, for a function without an
analytical derivative: pass `et -> uddc(udfuns, et, dx)` as `udqdec`.

Throws `SpiceError(:dividebyzero)` if `dx == 0`.

# Examples
```jldoctest
julia> uddc(cos, 1.0, 1e-3)
true
```

See also [`uddf`](@ref), [`gfuds`](@ref).
"""
uddc(udfunc::F, x::Real, dx::Real) where {F} = uddf(udfunc, x, dx) < 0.0
