# Group 9: the SPK data types other than 2 against CSPICE.
#
# Type 2 is the type of `de440s.bsp` and `de405.bsp` (groups 1 to 4). No kernel of the
# other types, type 3 included, is on this machine, so the group writes them: one SPK file
# with one segment per case, made with the segment writers of CSPICE (`spkwNN_c`, the
# `spk14b_c`/`spk14a_c`/`spk14e_c` sequence, and the Fortran-interface routines
# `spkw01_`, `spkw19_`, `spkw21_` for the types that have no C wrapper). Each case has
# its own target body, relative to the Earth (399), so a query needs no other kernel.
#
# The segment writers are called directly with `ccall`: the wrappers of SPICE.jl for
# some of them take the state count from the wrong dimension.
#
# Checks of each case:
#   (a) the `spkgeo` state against CSPICE at random epochs, at the two ends of the
#       segment and at the special epochs of the case (knots, record bounds), where the
#       record or window selection of CSPICE is decided;
#   (b) no allocation and inferred result types of `spkgeo` and `spkder`;
#   (c) orders 2 and 3 of `spkder` against the nested `ForwardDiff.derivative` of the
#       velocity part of `spkgeo` (the generic path with a `Dual` epoch), and order 2
#       against a central finite difference of the velocity;
#   (d) the orders of `spk_eval_position` against the nested derivative of `spkpos`.
#
# Tolerance of (a). The two engines make the same operations in the same order on the
# same numbers, so the states differ by a few units in the last place. The bound is
# `TYPE_RTOL` on |difference| / |vector|.
#
# Tolerance of (c). The Taylor path and the dual path differentiate the same function,
# so they differ in rounding only: `AD_VEL_RTOL` per order. The finite difference has
# the step `h` of the case and the bound `FD_VEL_RTOL`; it is the independent check that
# both paths carry the derivative through the Kepler solutions of types 5, 15 and 17.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff

const LIBCSPICE = SPICE.libcspice
const CInt = Cint                  # `SpiceInt`, and `integer`, `logical`, `ftnlen` of f2c

const TYPE_RTOL = 5e-13
const AD_VEL_RTOL = (1e-10, 1e-8)   # orders 2 and 3 of spkder
const AD_POS_RTOL = (1e-13, 1e-11, 1e-9, 1e-7)   # orders 0 to 3 of spk_eval_position
const FD_VEL_RTOL = 1e-5

const EARTH = 399
const GM_EARTH = 398600.4418

"""
One generated segment. `epochs` are the special epochs (knots, record bounds) inside
`[t0, t1]`; `h` is the step of the finite difference of the velocity; `rtol` is the
bound of the state comparison with CSPICE.
"""
struct TypeCase
    label::String
    dtype::Int
    body::Int
    t0::Float64
    t1::Float64
    epochs::Vector{Float64}
    h::Float64
    rtol::Float64
end

"Raise the pending CSPICE error, if there is one."
cspice_check() = SPICE.handleerror()

# --- reference trajectories -------------------------------------------------------------

"State (6-vector) at `dt` seconds after `state0` on the two-body orbit about the Earth."
orbit_state(state0, dt) = collect(SPICE.prop2b(GM_EARTH, state0, dt))

"Two-body acceleration at the position `r`."
orbit_acc(r) = -GM_EARTH .* r ./ norm(r)^3

"Ellipse: perigee radius 7000 km, eccentricity 0.2, inclined; period about 2.3 hours."
const ELLIPSE0 = [7000.0, 0.0, 0.0, 0.0, 7.0 * cosd(40.0) * 1.1783, 7.0 * sind(40.0) * 1.1783]
"Hyperbola: perigee radius 8000 km, speed above the escape speed."
const HYPERBOLA0 = [8000.0, 0.0, 0.0, 0.0, 11.0, 4.0]

"6 x n matrix of the states of the orbit from `state0` at `epochs` (time zero is `epochs[1]`)."
function state_matrix(state0, epochs)
    m = Matrix{Float64}(undef, 6, length(epochs))
    for (j, t) in enumerate(epochs)
        m[:, j] = orbit_state(state0, t - epochs[1])
    end
    return m
end

"12 x n matrix of type 18 subtype 0 packets: position, velocity, velocity, acceleration."
function packets12(states)
    p = Matrix{Float64}(undef, 12, size(states, 2))
    for j in axes(states, 2)
        p[1:6, j] = states[:, j]
        p[7:9, j] = states[4:6, j]
        p[10:12, j] = orbit_acc(states[1:3, j])
    end
    return p
end

"`n` ascending epochs from `t0` with steps of `step` seconds that vary by ±30 percent."
function uneven_epochs(rng, n, t0, step)
    out = Vector{Float64}(undef, n)
    out[1] = t0
    for i in 2:n
        out[i] = out[i-1] + step * (0.7 + 0.6 * rand(rng))
    end
    return out
end

# --- segment writers ------------------------------------------------------------------

function write_type5!(h, body, frame, states, epochs)
    ccall((:spkw05_c, LIBCSPICE), Cvoid, (CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring,
            Cdouble, CInt, Ref{Cdouble}, Ref{Cdouble}), h, body, EARTH, frame, epochs[1],
        epochs[end], "type 5", GM_EARTH, length(epochs), states, epochs)
    cspice_check()
end

"Types 8 and 12: states at equal steps."
function write_equal!(h, dtype, body, frame, degree, states, epoch1, step)
    n = size(states, 2)
    last = epoch1 + (n - 1) * step
    types = (CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring, CInt, CInt, Ref{Cdouble},
        Cdouble, Cdouble)
    if dtype == 8
        ccall((:spkw08_c, LIBCSPICE), Cvoid, (CInt, CInt, CInt, Cstring, Cdouble, Cdouble,
                Cstring, CInt, CInt, Ref{Cdouble}, Cdouble, Cdouble), h, body, EARTH, frame,
            epoch1, last, "type 8", degree, n, states, epoch1, step)
    else
        ccall((:spkw12_c, LIBCSPICE), Cvoid, (CInt, CInt, CInt, Cstring, Cdouble, Cdouble,
                Cstring, CInt, CInt, Ref{Cdouble}, Cdouble, Cdouble), h, body, EARTH, frame,
            epoch1, last, "type 12", degree, n, states, epoch1, step)
    end
    cspice_check()
    return last
end

"Types 9 and 13: states at arbitrary epochs."
function write_unequal!(h, dtype, body, frame, degree, states, epochs)
    n = length(epochs)
    if dtype == 9
        ccall((:spkw09_c, LIBCSPICE), Cvoid, (CInt, CInt, CInt, Cstring, Cdouble, Cdouble,
                Cstring, CInt, CInt, Ref{Cdouble}, Ref{Cdouble}), h, body, EARTH, frame,
            epochs[1], epochs[end], "type 9", degree, n, states, epochs)
    else
        ccall((:spkw13_c, LIBCSPICE), Cvoid, (CInt, CInt, CInt, Cstring, Cdouble, Cdouble,
                Cstring, CInt, CInt, Ref{Cdouble}, Ref{Cdouble}), h, body, EARTH, frame,
            epochs[1], epochs[end], "type 13", degree, n, states, epochs)
    end
    cspice_check()
end

function write_type18!(h, subtyp, body, frame, degree, packets, epochs)
    ccall((:spkw18_c, LIBCSPICE), Cvoid, (CInt, CInt, CInt, CInt, Cstring, Cdouble, Cdouble,
            Cstring, CInt, CInt, Ref{Cdouble}, Ref{Cdouble}), h, subtyp, body, EARTH, frame,
        epochs[1], epochs[end], "type 18", degree, length(epochs), packets, epochs)
    cspice_check()
end

"""
Type 19 through the Fortran-interface routine. `minis` is a vector of
`(subtype, degree, packets, epochs)`; `bounds` are the interval bounds.
"""
function write_type19!(h, body, frame, minis, bounds, later::Bool)
    @assert length(bounds) == length(minis) + 1
    nintvl = Ref{CInt}(length(minis))
    npkts = CInt[length(m[4]) for m in minis]
    subtps = CInt[m[1] for m in minis]
    degres = CInt[m[2] for m in minis]
    packts = reduce(vcat, [vec(m[3]) for m in minis])
    epochs = reduce(vcat, [m[4] for m in minis])
    seg = "type 19"
    ccall((:spkw19_, LIBCSPICE), CInt, (Ref{CInt}, Ref{CInt}, Ref{CInt}, Ptr{UInt8},
            Ref{Cdouble}, Ref{Cdouble}, Ptr{UInt8}, Ref{CInt}, Ref{CInt}, Ref{CInt}, Ref{CInt},
            Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{CInt}, CInt, CInt),
        h, body, EARTH, frame, bounds[1], bounds[end], seg, nintvl, npkts, subtps, degres,
        packts, epochs, bounds, CInt(later), ncodeunits(frame), ncodeunits(seg))
    cspice_check()
end

"Types 1 and 21 through the Fortran-interface routines; `dlines` has one line per column."
function write_mda!(h, dtype, body, frame, first, dlines, epochs)
    seg = "type $dtype"
    n = Ref{CInt}(length(epochs))
    if dtype == 1
        @assert size(dlines, 1) == 71
        ccall((:spkw01_, LIBCSPICE), CInt, (Ref{CInt}, Ref{CInt}, Ref{CInt}, Ptr{UInt8},
                Ref{Cdouble}, Ref{Cdouble}, Ptr{UInt8}, Ref{CInt}, Ref{Cdouble}, Ref{Cdouble},
                CInt, CInt), h, body, EARTH, frame, first, epochs[end], seg, n, dlines, epochs,
            ncodeunits(frame), ncodeunits(seg))
    else
        ccall((:spkw21_, LIBCSPICE), CInt, (Ref{CInt}, Ref{CInt}, Ref{CInt}, Ptr{UInt8},
                Ref{Cdouble}, Ref{Cdouble}, Ptr{UInt8}, Ref{CInt}, Ref{CInt}, Ref{Cdouble},
                Ref{Cdouble}, CInt, CInt), h, body, EARTH, frame, first, epochs[end], seg, n,
            size(dlines, 1), dlines, epochs, ncodeunits(frame), ncodeunits(seg))
    end
    cspice_check()
end

"""
Difference lines of order `m` (15 for type 1) with final epochs `epochs`: reference
epoch inside the interval of the line, step sizes of the size of the interval, a
reference state of an orbit and a random difference table that decreases with the
order, as the table of an integrator does.
"""
function difference_lines(rng, m, t0, epochs)
    dl = zeros(4m + 11, length(epochs))
    prev = t0
    for (i, tf) in enumerate(epochs)
        span = tf - prev
        dl[1, i] = prev + span * (0.3 + 0.4 * rand(rng))                 # TL
        dl[2:m+1, i] = span .* (0.5 .+ rand(rng, m))                     # G
        s = orbit_state(ELLIPSE0, dl[1, i] - t0)
        dl[m+2:2:m+6, i] = s[1:3]
        dl[m+3:2:m+7, i] = s[4:6]
        kq = rand(rng, m-5:m-1, 3)      # KQMAX1 = max(KQ) + 1 is not above m
        for c in 1:3, j in 1:m
            dl[m+7+(c-1)*m+j, i] = 1e-3 * randn(rng) / 2.0^j             # DT(j, c)
        end
        dl[4m+8, i] = maximum(kq) + 1                                    # KQMAX1
        dl[4m+9:4m+11, i] = kq
        prev = tf
    end
    return dl
end

"Chebyshev coefficients that decrease with the degree; `scale` is the size of term 0."
cheb_coefs(rng, ncoef, scale) = [scale * randn(rng) / 3.0^(k - 1) for k in 1:ncoef]

function write_type14!(rng, h, body, frame, degree, bounds)
    n = length(bounds) - 1
    ncoef = degree + 1
    recs = Matrix{Float64}(undef, 2 + 6 * ncoef, n)
    for i in 1:n
        recs[1, i] = (bounds[i] + bounds[i+1]) / 2
        recs[2, i] = (bounds[i+1] - bounds[i]) / 2
        recs[3:end, i] = vcat((cheb_coefs(rng, ncoef, c <= 3 ? 7000.0 : 7.0) for c in 1:6)...)
    end
    ccall((:spk14b_c, LIBCSPICE), Cvoid, (CInt, Cstring, CInt, CInt, Cstring, Cdouble, Cdouble,
            CInt), h, "type 14", body, EARTH, frame, bounds[1], bounds[end], degree)
    cspice_check()
    # Two calls, to cover the addition of records to a segment that is open.
    half = n ÷ 2
    ccall((:spk14a_c, LIBCSPICE), Cvoid, (CInt, CInt, Ref{Cdouble}, Ref{Cdouble}), h, half,
        recs[:, 1:half], bounds[1:half])
    cspice_check()
    ccall((:spk14a_c, LIBCSPICE), Cvoid, (CInt, CInt, Ref{Cdouble}, Ref{Cdouble}), h, n - half,
        recs[:, half+1:n], bounds[half+1:n])
    cspice_check()
    ccall((:spk14e_c, LIBCSPICE), Cvoid, (CInt,), h)
    cspice_check()
end

"Type 20: `n` records of `intlen` days from the Julian date `initjd + initfr`; returns the coverage."
function write_type20!(rng, h, body, frame, polydg, n, intlen, dscale, tscale, initjd, initfr)
    recs = Matrix{Float64}(undef, 3 * (polydg + 2), n)
    for i in 1:n, c in 1:3
        rows = (c - 1) * (polydg + 2) .+ (1:polydg+2)
        recs[rows, i] = vcat(cheb_coefs(rng, polydg + 1, 0.01), 0.5 * randn(rng))
    end
    first = (initjd - 2451545.0 + initfr) * 86400.0
    last = first + n * intlen * 86400.0
    ccall((:spkw20_c, LIBCSPICE), Cvoid, (CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring,
            Cdouble, CInt, CInt, Ref{Cdouble}, Cdouble, Cdouble, Cdouble, Cdouble), h, body,
        EARTH, frame, first, last, "type 20", intlen, n, polydg, recs, dscale, tscale, initjd,
        initfr)
    cspice_check()
    return first, last
end

"Type 3: `n` records of `intlen` seconds from `t0`; the velocity series is not the derivative of the position series."
function write_type3!(rng, h, body, frame, degree, n, t0, intlen)
    @assert n >= 1 && intlen > 0
    ncoef = degree + 1
    recs = Matrix{Float64}(undef, 6 * ncoef, n)
    for i in 1:n
        recs[:, i] = vcat((cheb_coefs(rng, ncoef, c <= 3 ? 7000.0 : 7.0) for c in 1:6)...)
    end
    last = t0 + n * intlen
    ccall((:spkw03_c, LIBCSPICE), Cvoid, (CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring,
            Cdouble, CInt, CInt, Ref{Cdouble}, Cdouble), h, body, EARTH, frame, t0, last, "type 3",
        intlen, n, degree, recs, t0)
    cspice_check()
    return last
end

function write_type15!(h, body, frame, t0, t1, epoch, tp, pa, p, ecc, j2flg, pv, j2, radius)
    ccall((:spkw15_c, LIBCSPICE), Cvoid, (CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring,
            Cdouble, Ref{Cdouble}, Ref{Cdouble}, Cdouble, Cdouble, Cdouble, Ref{Cdouble},
            Cdouble, Cdouble, Cdouble), h, body, EARTH, frame, t0, t1, "type 15", epoch, tp, pa,
        p, ecc, Float64(j2flg), pv, GM_EARTH, j2, radius)
    cspice_check()
end

function write_type17!(h, body, frame, t0, t1, epoch, eqel, rapol, decpol)
    ccall((:spkw17_c, LIBCSPICE), Cvoid, (CInt, CInt, CInt, Cstring, Cdouble, Cdouble, Cstring,
            Cdouble, Ref{Cdouble}, Cdouble, Cdouble), h, body, EARTH, frame, t0, t1, "type 17",
        epoch, eqel, rapol, decpol)
    cspice_check()
end

# --- the test kernel ------------------------------------------------------------------

"""
    write_type_kernel(path, rng) -> Vector{TypeCase}

Write the SPK file `path` with one segment per case and return the cases.
Postcondition: the file is closed and every case has `t0 < t1`.
"""
function write_type_kernel(path::AbstractString, rng::AbstractRNG)
    cases = TypeCase[]
    h = SPICE.spkopn(path, "Allspice test kernel", 0)
    body = -9000
    add!(label, dtype, t0, t1, epochs, hfd; rtol=TYPE_RTOL) = push!(cases,
        TypeCase(label, dtype, body, t0, t1, [t for t in epochs if t0 <= t <= t1], hfd, rtol))
    J = "J2000"

    # Type 5: more than 100 states, so that the epoch directory is in use.
    ep = uneven_epochs(rng, 230, 1000.0, 60.0)
    body -= 1; write_type5!(h, body, J, state_matrix(ELLIPSE0, ep), ep)
    add!("type 5, ellipse", 5, ep[1], ep[end], ep, 1e-2)
    ep = uneven_epochs(rng, 40, -5000.0, 3000.0)
    body -= 1; write_type5!(h, body, J, state_matrix(HYPERBOLA0, ep), ep)
    add!("type 5, hyperbola", 5, ep[1], ep[end], ep, 1.0)
    ep = [0.0, 40000.0]
    body -= 1; write_type5!(h, body, J, state_matrix(ELLIPSE0, ep), ep)
    add!("type 5, two states, many revolutions", 5, ep[1], ep[end], ep, 1e-2)

    # Types 8 and 12: even and odd window sizes.
    for (dtype, degree) in ((8, 7), (8, 4), (8, 1), (12, 7), (12, 5), (12, 1))
        n, t0, step = 150, 2000.0, 60.0
        ep = [t0 + (i - 1) * step for i in 1:n]
        body -= 1
        t1 = write_equal!(h, dtype, body, J, degree, state_matrix(ELLIPSE0, ep), t0, step)
        # Knots, and the middles between knots, where an odd window moves.
        add!("type $dtype, degree $degree", dtype, t0, t1, vcat(ep, ep[1:end-1] .+ step / 2), 1e-3)
    end

    # Types 9 and 13, with the largest degree and one segment in another frame.
    for (dtype, degree, frame) in ((9, 9, J), (9, 6, J), (9, 27, J), (9, 3, "ECLIPJ2000"),
            (13, 7, J), (13, 5, J), (13, 27, J), (13, 3, "ECLIPJ2000"))
        ep = uneven_epochs(rng, 230, 500.0, 60.0)
        body -= 1; write_unequal!(h, dtype, body, frame, degree, state_matrix(ELLIPSE0, ep), ep)
        mids = (ep[1:end-1] .+ ep[2:end]) ./ 2
        add!("type $dtype, degree $degree, $frame", dtype, ep[1], ep[end], vcat(ep, mids), 1e-3)
    end
    # The smallest segment that the writer accepts: as many states as the window.
    ep = uneven_epochs(rng, 4, 0.0, 60.0)
    body -= 1; write_unequal!(h, 9, body, J, 3, state_matrix(ELLIPSE0, ep), ep)
    add!("type 9, degree 3, 4 states", 9, ep[1], ep[end], ep, 1e-3)

    # Type 18: both subtypes, a small and the largest window.
    for (subtyp, degree) in ((0, 7), (0, 15), (1, 7), (1, 15))
        ep = uneven_epochs(rng, 230, 500.0, 60.0)
        st = state_matrix(ELLIPSE0, ep)
        body -= 1; write_type18!(h, subtyp, body, J, degree, subtyp == 0 ? packets12(st) : st, ep)
        add!("type 18, subtype $subtyp, degree $degree", 18, ep[1], ep[end], ep, 1e-3)
    end

    # Type 19: mini-segments of the three subtypes, with packets outside the interval
    # bounds (padding) and one mini-segment with a packet directory; both values of
    # the selection flag; and a segment with more than 100 intervals.
    for later in (true, false)
        minis = Any[]
        bounds = Float64[0.0]
        for (subtyp, degree, npkt) in ((0, 7, 30), (1, 7, 130), (2, 11, 25), (1, 3, 5), (2, 3, 2))
            ep = uneven_epochs(rng, npkt, bounds[end] - (npkt > 6 ? 150.0 : 0.0), 60.0)
            push!(bounds, npkt > 6 ? ep[end-3] + 10.0 : ep[end])
            st = state_matrix(ELLIPSE0, ep)
            push!(minis, (subtyp, degree, subtyp == 0 ? packets12(st) : st, ep))
        end
        body -= 1; write_type19!(h, body, J, minis, bounds, later)
        knots = vcat(bounds, (m[4] for m in minis)...)
        add!("type 19, mixed subtypes, select " * (later ? "later" : "earlier"), 19, bounds[1],
            bounds[end], knots, 1e-3)
    end
    minis = Any[]
    bounds = Float64[0.0]
    for i in 1:130
        ep = [bounds[end] + 50.0 * (j - 1) for j in 1:6]
        push!(bounds, ep[end])
        push!(minis, (isodd(i) ? 1 : 2, 3, state_matrix(ELLIPSE0, ep), ep))
    end
    body -= 1; write_type19!(h, body, J, minis, bounds, true)
    add!("type 19, 130 intervals", 19, bounds[1], bounds[end], bounds, 1e-3)

    # Types 1 and 21: more than 100 difference lines.
    for (dtype, m) in ((1, 15), (21, 15), (21, 25))
        t0 = 100.0
        ep = uneven_epochs(rng, 130, t0 + 80.0, 80.0)
        body -= 1; write_mda!(h, dtype, body, J, t0, difference_lines(rng, m, t0, ep), ep)
        add!("type $dtype, order $m", dtype, t0, ep[end], ep, 1e-3)
    end

    # Type 14: records of different lengths, more than 100 of them.
    bounds = uneven_epochs(rng, 131, -3000.0, 900.0)
    body -= 1; write_type14!(rng, h, body, J, 9, bounds)
    add!("type 14, degree 9", 14, bounds[1], bounds[end], bounds, 1e-2)
    bounds = uneven_epochs(rng, 4, 0.0, 900.0)
    body -= 1; write_type14!(rng, h, body, "ECLIPJ2000", 2, bounds)
    add!("type 14, degree 2, ECLIPJ2000", 14, bounds[1], bounds[end], bounds, 1e-2)

    # Type 20: scaled units (au and days) and a start that is a fraction of a day.
    for (polydg, n, intlen, dscale, tscale) in ((8, 40, 0.5, 149597870.7, 86400.0),
            (1, 3, 2.0, 1.0, 1.0), (12, 7, 1.25, 1000.0, 3600.0))
        body -= 1
        t0, t1 = write_type20!(rng, h, body, J, polydg, n, intlen, dscale, tscale, 2451545.0, 0.25)
        add!("type 20, degree $polydg", 20, t0, t1, [t0 + k * intlen * 86400.0 for k in 0:n], 1.0)
    end

    # Type 15: the four values of the J2 flag, many revolutions on both sides of the
    # periapsis epoch, and a hyperbola (no J2 drift).
    tp = [sind(30.0) * sind(50.0), -cosd(30.0) * sind(50.0), cosd(50.0)]     # orbit normal
    pa = normalize(cross([0.0, 0.0, 1.0], tp))                               # in the plane
    pv = [0.0, 0.0, 1.0]
    for (j2flg, ecc, p) in ((0, 0.2, 9000.0), (1, 0.2, 9000.0), (2, 0.05, 8000.0),
            (3, 0.2, 9000.0), (0, 0.6, 12000.0), (0, 1.5, 20000.0), (0, 0.0, 7500.0))
        epoch = 1234.5
        t0, t1 = ecc > 1 ? (epoch - 4.0e4, epoch + 6.0e4) : (epoch - 8.0e4, epoch + 1.2e5)
        body -= 1
        write_type15!(h, body, J, t0, t1, epoch, tp, pa, p, ecc, j2flg, pv, 1.0826e-3, 6378.14)
        add!("type 15, flag $j2flg, ecc $ecc", 15, t0, t1, [epoch], 1e-2)
    end

    # Type 17: an ellipse about a pole that is not the frame pole, and a near-circle.
    for (hh, kk, a) in ((0.12, -0.09, 12000.0), (1e-4, 2e-4, 42164.0), (0.6, 0.55, 30000.0))
        epoch = -777.0
        n = sqrt(GM_EARTH / a^3)
        eqel = [a, hh, kk, 0.7, 0.2, -0.15, 3.0e-7, n * (1 + 1e-3), -5.0e-7]
        body -= 1
        write_type17!(h, body, J, epoch - 1.0e5, epoch + 2.0e5, epoch, eqel, deg2rad(20.0), deg2rad(70.0))
        add!("type 17, h $hh, k $kk", 17, epoch - 1.0e5, epoch + 2.0e5, [epoch], 1e-2)
    end

    # Type 3: position and velocity series apart. It is the last case, so that the cases
    # before it keep their random data.
    t0 = -2000.0
    body -= 1
    t1 = write_type3!(rng, h, body, J, 9, 12, t0, 900.0)
    add!("type 3, degree 9", 3, t0, t1, [t0 + k * 900.0 for k in 0:12], 1e-2)

    SPICE.spkcls(h)
    @assert all(c -> c.t0 < c.t1, cases)
    return cases
end

# --- comparisons ----------------------------------------------------------------------

"Epochs of the state comparison of `c`: ends, special epochs, points just inside the ends, random points."
function case_epochs(rng, c::TypeCase, nrandom::Int)
    span = c.t1 - c.t0
    special = length(c.epochs) > 400 ? c.epochs[1:(length(c.epochs)÷400+1):end] : c.epochs
    near = [c.t0 + span * 1e-9, c.t1 - span * 1e-9, c.t0 + span * 1e-3, c.t1 - span * 1e-3]
    return vcat(c.t0, c.t1, special, near, random_epochs(rng, nrandom, c.t0, c.t1))
end

"Largest relative difference of position and of velocity between Allspice and CSPICE for `c`."
function state_errors(ctx, c::TypeCase, epochs)
    worst_r = 0.0
    worst_v = 0.0
    for et in epochs
        ref, _ = SPICE.spkgeo(c.body, et, "J2000", EARTH)
        state, _ = spkgeo(ctx, c.body, et, "J2000", EARTH)
        worst_r = max(worst_r, norm(state[1:3] - ref[1:3]) / norm(ref[1:3]))
        worst_v = max(worst_v, norm(state[4:6] - ref[4:6]) / norm(ref[4:6]))
    end
    return worst_r, worst_v
end

"Allocation of the second call of `f(args...)` (the first call compiles)."
function warm_alloc(f, args...)
    f(args...)
    return f(args...)
end

"Relative difference `|a - b| / |b|`; for `b == 0` the absolute difference."
function type_reldiff(a, b)
    nb = norm(b)
    return nb == 0 ? norm(a - b) : norm(a - b) / nb
end

"""
Nesting depth of ForwardDiff through `spk_eval` in check (d). Each depth is a new
compilation of the evaluators of all the data types (about 0.7 s, 1 s and 1.4 s), so
the default run nests to depth 1 and the exhaustive run to depth 3.
"""
const TYPE_AD_ORDER = nsample(1, 3)

"The function `t -> d^N f / dt^N` by nested forward mode; `f` itself for `N == 0`."
type_ad_nth(f, ::Val{0}) = f
function type_ad_nth(f, ::Val{N}) where {N}
    g = type_ad_nth(f, Val(N - 1))
    return s -> ForwardDiff.derivative(g, s)
end

"Nested forward-mode derivatives of `f` at `t`, orders 0 to `N`, as a tuple. Requires `N >= 0`."
function type_ad_jet(f, t::Float64, ::Val{N}) where {N}
    @assert N isa Int && N >= 0
    return ntuple(k -> type_ad_nth(f, Val(k - 1))(t), Val(N + 1))
end

alloc_type_geo(ctx, body, et) = @allocated spkgeo(ctx, body, et, 1, EARTH)
alloc_type_der(ctx, body, et, ::Val{N}) where {N} = @allocated spkder(ctx, body, et, 1, EARTH, Val(N))

"Epochs of `c` for the derivative checks: random, and not within `100 h` of a special epoch."
function smooth_epochs(rng, c::TypeCase, n::Int)
    out = Float64[]
    while length(out) < n
        t = c.t0 + (c.t1 - c.t0) * (0.02 + 0.96 * rand(rng))
        all(s -> abs(s - t) > 100 * c.h, c.epochs) && push!(out, t)
    end
    return out
end

group("9. SPK data types against CSPICE") do
    rng = MersenneTwister(SEED + 9)
    path = joinpath(mktempdir(), "allspice_types.bsp")
    cases = write_type_kernel(path, rng)
    oracle_load([path])
    ctx = neo_context([path])
    ks = Allspice.kernelset(ctx)
    println("  ", length(cases), " segments in ", path)
    @test sort(unique(c.dtype for c in cases)) == [1, 3, 5, 8, 9, 12, 13, 14, 15, 17, 18, 19, 20, 21]

    @apart @testset "states against CSPICE" begin
        by_type = Dict{Int,Tuple{Float64,Float64}}()
        for c in cases
            wr, wv = state_errors(ctx, c, case_epochs(rng, c, 300))
            measured(@sprintf("  %-44s position %.2e  velocity %.2e", c.label, wr, wv))
            @test wr <= c.rtol
            @test wv <= c.rtol
            old = get(by_type, c.dtype, (0.0, 0.0))
            by_type[c.dtype] = (max(old[1], wr), max(old[2], wv))
        end
        for dtype in sort(collect(keys(by_type)))
            measured(@sprintf("  type %-2d  max relative difference: position %.2e  velocity %.2e",
                dtype, by_type[dtype]...))
        end
    end

    @apart @testset "name forms and spkezr give the spkgeo state" begin
        c = cases[1]
        et = (c.t0 + c.t1) / 2
        byid = spkgeo(ctx, c.body, et, "J2000", EARTH)
        @test spkezr(ctx, string(c.body), et, "J2000", "NONE", "EARTH") == byid
        @test spkpos(ctx, c.body, et, "J2000", "NONE", EARTH)[1] == byid[1][SOneTo(3)]
    end

    @apart @testset "zero allocation and inferred types" begin
        for c in cases
            et = c.t0 + 0.37 * (c.t1 - c.t0)
            @testset "$(c.label)" begin
                @test warm_alloc(alloc_type_geo, ctx, c.body, et) == 0 skip = STACK_BUFFERS_ESCAPE
                for N in 0:3
                    @test warm_alloc(alloc_type_der, ctx, c.body, et, Val(N)) == 0 skip = STACK_BUFFERS_ESCAPE
                end
                @test (@inferred spkgeo(ctx, c.body, et, 1, EARTH)) isa Tuple{SVector{6,Float64},Float64}
                @test (@inferred spkder(ctx, c.body, et, 1, EARTH, Val(3))) isa NTuple{4,SVector{3,Float64}}
            end
        end
    end

    @apart @testset "spkder: lower orders, and higher orders against ForwardDiff and differences" begin
        for c in cases
            vel(t) = spkgeo(ctx, c.body, t, 1, EARTH)[1][SVector(4, 5, 6)]
            acc(t) = ForwardDiff.derivative(vel, t)
            jerk(t) = ForwardDiff.derivative(acc, t)
            worst = zeros(3)            # order 2, order 3 against AD; order 2 against FD
            for et in smooth_epochs(rng, c, nsample(3, 8))
                jet = spkder(ctx, c.body, et, 1, EARTH, Val(3))
                state, _ = spkgeo(ctx, c.body, et, 1, EARTH)
                # Equal to the last bit in a normal run; with coverage instrumentation the
                # compiler arranges the sums differently, so one rounding is allowed.
                @test jet[1] ≈ state[SOneTo(3)] rtol = 4eps()
                @test jet[2] ≈ state[SVector(4, 5, 6)] rtol = 4eps()
                @test spkder(ctx, c.body, et, 1, EARTH, Val(2))[3] ≈ jet[3] rtol = 1e-12
                # The dual path gives the same plain values.
                dstate = spkgeo(ctx, c.body, ForwardDiff.Dual(et, 1.0), 1, EARTH)[1]
                @test ForwardDiff.value.(dstate) ≈ state rtol = 1e-14
                fd = (vel(et + c.h) - vel(et - c.h)) / (2 * c.h)
                worst .= max.(worst, (type_reldiff(jet[3], acc(et)), type_reldiff(jet[4], jerk(et)),
                    type_reldiff(jet[3], fd)))
            end
            measured(@sprintf("  %-44s AD order 2 %.1e, 3 %.1e; FD order 2 %.1e", c.label, worst...))
            @test worst[1] <= AD_VEL_RTOL[1]
            @test worst[2] <= AD_VEL_RTOL[2]
            @test worst[3] <= FD_VEL_RTOL
        end
    end

    @apart @testset "spk_eval_position: every order is a derivative of the position" begin
        for c in cases
            seg = ks.spk[Int32(c.body)][1]
            data = ks.daf[seg.file]
            pos(t) = Allspice.spk_eval(data, seg, t, Val(0))[1]
            worst = zeros(TYPE_AD_ORDER + 1)
            for et in smooth_epochs(rng, c, nsample(2, 4))
                jet = Allspice.spk_eval_position(data, seg, et, Val(3))
                ad = type_ad_jet(pos, et, Val(TYPE_AD_ORDER))
                worst .= max.(worst, [type_reldiff(jet[k], ad[k]) for k in 1:TYPE_AD_ORDER+1])
                @test (@inferred Allspice.spk_eval_position(data, seg, et, Val(3))) isa
                      NTuple{4,SVector{3,Float64}}
                @test Allspice.spk_eval_position(data, seg, et, Val(0))[1] == jet[1]
            end
            @test all(worst[k] <= AD_POS_RTOL[k] for k in 1:TYPE_AD_ORDER+1)
        end
    end

    @apart @testset "epoch number types: Taylor epoch, and the Kepler solutions at their special epochs" begin
        for c in cases
            seg = ks.spk[Int32(c.body)][1]
            data = ks.daf[seg.file]
            et = c.t0 + 0.61 * (c.t1 - c.t0)
            jet = Allspice.spk_eval(data, seg, et, Val(3))
            # A Taylor epoch of order 1 with Val(2): series of series. The plain parts
            # are the jet of order 2 and the first-order parts the next orders.
            nested = Allspice.spk_eval(data, seg, Allspice.taylor_var(et, Val(1)), Val(2))
            @test all(map(x -> x.v, nested[k]) ≈ jet[k] for k in 1:3)
            @test map(x -> x.d[1], nested[2]) ≈ jet[3] rtol = 1e-9
            @test map(x -> x.d[1], nested[3]) ≈ jet[4] rtol = 1e-7
        end
        # Type 15 at the periapsis epoch (true anomaly zero, time interval zero) and type
        # 5 at a state epoch: the derivatives are finite and agree with the differences.
        for c in cases
            (c.dtype == 15 || c.dtype == 17 || c.dtype == 5) || continue
            et = c.dtype == 5 ? c.epochs[2] : c.epochs[1]
            vel(t) = spkgeo(ctx, c.body, t, 1, EARTH)[1][SVector(4, 5, 6)]
            jet = spkder(ctx, c.body, et, 1, EARTH, Val(3))
            @test all(all(isfinite, j) for j in jet)
            if c.dtype != 5       # at a type 5 state epoch the pair of states changes
                fd = (vel(et + c.h) - vel(et - c.h)) / (2 * c.h)
                @test type_reldiff(jet[3], fd) <= FD_VEL_RTOL
            end
            @test type_reldiff(jet[3], ForwardDiff.derivative(vel, et)) <= AD_VEL_RTOL[1]
        end
    end

    @apart @testset "numeric helpers" begin
        tv = Allspice.taylor_var
        der = Allspice.derivatives
        # sinh and cosh: the derivatives alternate.
        x = 0.7
        @test collect(der(sinh(tv(x, Val(4))))) ≈ [sinh(x), cosh(x), sinh(x), cosh(x), sinh(x)]
        @test collect(der(cosh(tv(x, Val(4))))) ≈ [cosh(x), sinh(x), cosh(x), sinh(x), cosh(x)]
        @test der(sinh(2 * tv(x, Val(2))))[3] ≈ 4 * sinh(2x)
        # atan of a point that turns at unit rate while its radius changes: the angle is t.
        for t in (0.3, 1.5707963, 2.5, -2.9)
            tt = tv(t, Val(4))
            radius = 2.0 + 0.5 * tt * tt
            ang = atan(radius * sin(tt), radius * cos(tt))
            @test collect(der(ang)) ≈ [t, 1.0, 0.0, 0.0, 0.0] atol = 1e-12
        end
        # Orders that a number type carries, and the replacement of the plain value.
        D1 = ForwardDiff.Dual{Nothing,Float64,1}
        @test Allspice._ad_order(Float64) == 0
        @test Allspice._ad_order(Allspice.Taylor{3,Float64}) == 3
        @test Allspice._ad_order(D1) == 1
        @test Allspice._ad_order(ForwardDiff.Dual{Nothing,D1,1}) == 2
        @test Allspice._ad_order(Allspice.Taylor{2,D1}) == 3
        @test Allspice._ad_order(Rational{Int}) == 8
        @test Allspice._with_primal(3.0, 1.5) === 1.5
        y = Allspice._with_primal(tv(3.0, Val(2)), 1.5)
        @test y.c == (1.5, 1.0, 0.0)
        d = Allspice._with_primal(ForwardDiff.Dual(3.0, 2.0), 1.5)
        @test ForwardDiff.value(d) == 1.5 && ForwardDiff.partials(d, 1) == 2.0
        @test Allspice._with_primal(1 // 2, 0.25) == 0.25
        # Stumpff functions: the three branches meet the definitions.
        for z in (-9.0, -0.5, 0.0, 0.5, 9.0)
            c0, c1, c2, c3 = Allspice._stumpff(z)
            r = sqrt(abs(z))
            e0, e1 = z < 0 ? (cosh(r), sinh(r) / r) : z > 0 ? (cos(r), sin(r) / r) : (1.0, 1.0)
            @test c0 ≈ e0 && c1 ≈ e1
            @test 1 - z * c2 ≈ c0 && 1 - z * c3 ≈ c1
        end
        @test throws_kind(() -> Allspice._stumpff(-1.0e6), :badorbit)
        # Angle between vectors and rotation about an axis.
        ex, ey, ez = SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0)
        @test Allspice._vsep(ex, ey) == pi / 2
        @test Allspice._vsep(ex, -ex) ≈ pi
        @test Allspice._vsep(ex, SVector(1.0, 1.0, 0.0)) ≈ pi / 4
        @test Allspice._vsep(ex, zero(ex)) == 0.0
        @test Allspice._vrotv(ex, ez, pi / 2) ≈ ey
        @test Allspice._vrotv(ex, zero(ez), 1.0) == ex
        @test Allspice._vrotv(zero(ex), ez, 1.0) == zero(ex)
        @test Allspice._dsign(2.0, -0.0) == 2.0 && Allspice._dsign(-2.0, -1.0) == -2.0
        # Two-body propagation: zero time gives the state back; bad input is refused.
        r0, v0 = SVector(7000.0, 0.0, 0.0), SVector(0.0, 7.5, 1.0)
        @test Allspice._prop2b(GM_EARTH, r0, v0, 0.0) == (r0, v0)
        pd, vd = Allspice._prop2b(GM_EARTH, r0, v0, ForwardDiff.Dual(0.0, 1.0))
        @test ForwardDiff.value.(pd) ≈ r0 && ForwardDiff.partials.(pd, 1) ≈ v0
        @test ForwardDiff.partials.(vd, 1) ≈ orbit_acc(r0)
        for dt in (-5.0e4, 300.0, 9.0e5)
            ref = orbit_state(collect(vcat(r0, v0)), dt)
            p, v = Allspice._prop2b(GM_EARTH, r0, v0, dt)
            @test p ≈ ref[1:3] rtol = 1e-13
            @test v ≈ ref[4:6] rtol = 1e-13
        end
        @test throws_kind(() -> Allspice._prop2b(-1.0, r0, v0, 1.0), :badorbit)
        @test throws_kind(() -> Allspice._prop2b(GM_EARTH, zero(r0), v0, 1.0), :badorbit)
        @test throws_kind(() -> Allspice._prop2b(GM_EARTH, r0, 2.0 * r0, 1.0), :badorbit)
        @test throws_kind(() -> Allspice._prop2b(GM_EARTH, r0, 3.0 * v0, 1.0e308), :badorbit)
        @test Allspice._kepler_equinoctial(0.0, 0.3) == 0.0
        @test throws_kind(() -> Allspice._kepler_equinoctial(0.8, 0.8), :badorbit)
    end

    @apart @testset "errors" begin
        # Type 10 (two-line elements) is not available.
        seg10 = Allspice.Segment(Int32(-1), Int32(399), Int32(1), Int32(10), 0.0, 1.0, Int32(0),
            1, 20, 0.0, 0.0, Int32(0), Int32(0))
        @test throws_kind(() -> Allspice.spk_eval(zeros(20), seg10, 0.5, Val(1)), :notsupported)
        @test throws_kind(() -> Allspice.spk_eval_position(zeros(20), seg10, 0.5, Val(2)), :notsupported)
        # A count at the end of a segment that does not agree with its length. (For
        # type 3 the count is checked when the file is opened, as for type 2.)
        for c in cases
            c.dtype in (3, 15, 17, 20) && continue
            seg = ks.spk[Int32(c.body)][1]
            bad = copy(ks.daf[seg.file])
            bad[seg.eaddr] += 1.0
            @test throws_kind(() -> Allspice.spk_eval(bad, seg, (c.t0 + c.t1) / 2, Val(1)), :baddaf)
        end
        # One wrong value in the segment: `change(data, seg)` makes it, and the query
        # in the middle of the segment must throw an error of `kind`.
        function refuses(label, kind, change)
            c = cases[findfirst(x -> x.label == label, cases)]
            seg = ks.spk[Int32(c.body)][1]
            bad = copy(ks.daf[seg.file])
            change(bad, seg)
            return throws_kind(() -> Allspice.spk_eval(bad, seg, (c.t0 + c.t1) / 2, Val(1)), kind)
        end
        e(seg) = seg.eaddr
        b(seg) = seg.baddr
        @test refuses("type 8, degree 7", :baddaf, (d, s) -> d[e(s)-2] = 0.0)             # step
        @test refuses("type 8, degree 7", :baddaf, (d, s) -> d[e(s)-1] = 40.0)            # window
        @test refuses("type 8, degree 7", :baddaf, (d, s) -> d[e(s)] = NaN)
        @test refuses("type 9, degree 9, J2000", :baddaf, (d, s) -> d[e(s)-1] = 500.0)
        @test refuses("type 9, degree 9, J2000", :baddaf,                                 # equal epochs
            (d, s) -> d[b(s)+6*230:b(s)+7*230-1] .= 0.0)
        @test refuses("type 13, degree 7, J2000", :baddaf,
            (d, s) -> d[b(s)+6*230:b(s)+7*230-1] .= 0.0)
        @test refuses("type 18, subtype 0, degree 7", :baddaf, (d, s) -> d[e(s)-1] = 3.0)  # odd window
        @test refuses("type 18, subtype 0, degree 7", :baddaf, (d, s) -> d[e(s)-1] = 10.0)
        @test refuses("type 18, subtype 0, degree 7", :baddaf, (d, s) -> d[e(s)-2] = 2.0)  # subtype
        @test refuses("type 19, 130 intervals", :baddaf, (d, s) -> d[e(s)-2-131:e(s)-2] .= 1.0)
        @test refuses("type 1, order 15", :baddaf, (d, s) -> d[b(s):b(s)+71*130-1] .= 0.0)
        @test refuses("type 21, order 25", :baddaf, (d, s) -> d[e(s)-1] = 26.0)
        @test refuses("type 14, degree 9", :notsupported, (d, s) -> d[e(s)-12] = 2.0)       # RDRTYP
        @test refuses("type 14, degree 9", :baddaf, (d, s) -> d[e(s)-2] = 7.0)             # PKTSZ
        @test refuses("type 20, degree 8", :baddaf, (d, s) -> d[e(s)-5] = 0.0)             # TSCALE
        @test refuses("type 15, flag 0, ecc 0.2", :baddaf, (d, s) -> d[b(s)+7] = -1.0)     # P
        @test refuses("type 15, flag 0, ecc 0.2", :baddaf, (d, s) -> d[b(s)+1:b(s)+3] .= 0.0)
        @test refuses("type 15, flag 0, ecc 0.2", :baddaf,                                # PA along TP
            (d, s) -> d[b(s)+4:b(s)+6] .= d[b(s)+1:b(s)+3])
        @test refuses("type 17, h 0.12, k -0.09", :baddaf, (d, s) -> d[b(s)+2] = 0.95)     # H
        @test refuses("type 17, h 0.12, k -0.09", :baddaf, (d, s) -> d[b(s)+1] = 0.0)      # A
        @test refuses("type 5, ellipse", :badorbit, (d, s) -> d[e(s)-1] = 0.0)             # GM
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
