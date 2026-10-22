# Helpers shared by the programs.

_paths(x::AbstractString) = [String(x)]
_paths(x) = String[String(p) for p in x]

"A new context with the files of each list loaded in order; the last file has priority."
function _context(lists...)
    ctx = SpiceContext()
    for list in lists, p in _paths(list)
        furnsh(ctx, p)
    end
    return ctx
end

"Epoch in seconds past J2000 TDB: a number as is, a string through `str2et` on `ctx`."
_et(ctx, t::Real) = Float64(t)
_et(ctx, t::AbstractString) = str2et(ctx, t)

"NAIF ID code of a body given by its code or its name."
_bodyid(ctx, b::Integer) = Int(b)
function _bodyid(ctx, b::AbstractString)
    code = bodn2c(ctx, b)
    code === nothing && throw(ArgumentError("the body name \"$b\" is not known"))
    return code
end

"Name of the body `code`, or `\"\"` if it has none."
_bodyname(ctx, code::Integer) = something(bodc2n(ctx, code), "")

"Frame ID of a frame given by its ID or its name."
_frameid(ctx, f::Integer) = Int(f)
function _frameid(ctx, f::AbstractString)
    id = namfrm(ctx, f)
    id == 0 && throw(ArgumentError("the frame name \"$f\" is not known"))
    return id
end

"""
    _sample(window, n, step) -> Vector{Float64}

Epochs in `window`: with `step === nothing`, `n` epochs at equal distance along the
measure of the window, the first and the last at its ends; otherwise the epochs at
`step` from the left end of each interval, and the right end of each interval.
"""
function _sample(window::Window, n::Integer, step)
    @assert wncard(window) > 0
    m = sum(((a, b),) -> b - a, (wnfetd(window, i) for i in 1:wncard(window)))
    if step === nothing
        n >= 2 || throw(ArgumentError("the number of epochs must be at least 2, got $n"))
        if m == 0.0
            return fill(wnfetd(window, 1)[1], 1)
        end
        h = m / (n - 1)
        times = Float64[]
        for k in 0:(n - 1)
            s = k == n - 1 ? m : k * h     # distance along the measure
            for i in 1:wncard(window)    # invariant: s is past the intervals before i
                a, b = wnfetd(window, i)
                if s <= b - a || i == wncard(window)
                    push!(times, a + min(s, b - a))
                    break
                end
                s -= b - a
            end
        end
    else
        step > 0 || throw(ArgumentError("the step must be positive, got $step"))
        times = Float64[]
        for i in 1:wncard(window)
            a, b = wnfetd(window, i)
            k = 0
            while a + k * step < b      # invariant: the epochs a + j step, j < k, are in times
                push!(times, a + k * step)
                k += 1
            end
            push!(times, b)
        end
    end
    @assert issorted(times)
    return times
end

"Calendar TDB format of an epoch, as the NAIF programs print it."
_cal(et::Real) = etcal(Float64(et))
