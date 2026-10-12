# Kernels that the CK and SCLK groups write for themselves: no CK or SCLK kernel is on
# this machine, so the tests write a type 1 SCLK text kernel and a frame kernel by hand
# and CK files with the CSPICE writers, and then load them with both engines.
#
# The CK writers are called through `ccall` on `SPICE.libcspice`: `ckw01_c`, `ckw02_c`,
# `ckw03_c` and `ckw05_c` directly (the SPICE.jl wrappers of types 2 and 3 pass one
# record fewer than they are given), and the Fortran-style `ckw04b_`/`ckw04a_`/`ckw04e_`
# and `ckw06_`, which have no C wrapper.

module CKFix

import SPICE
using NeoSpice: taylor_var, derivatives

const LIB = SPICE.libcspice

# --- clocks ---------------------------------------------------------------------------------

"Clock with two fields, two partitions, four coefficient records, parallel time TDB."
const SC_TDB = -77
"Clock with three fields (the first counts from 1), three partitions, parallel time TDT."
const SC_TDT = -78
"Ticks in one second of clock `SC_TDB` (nominal), and nominal seconds per tick."
const TICKS_PER_SEC = 256.0
const RATE = 1 / TICKS_PER_SEC
"Ticks in one count of the first field of clock `SC_TDT` (60 * 50), one count is about 60 s."
const TDT_TIKMSC = 3000.0
const TDT_RATE = 60.0 / TDT_TIKMSC

"Partition bounds of the two clocks in ticks, as written to the kernel."
const TDB_PSTART = [0.0, 512000.0]
const TDB_PSTOP = [2560000.0, 1.0e12]
const TDT_PSTART = [0.0, 6000.0, 0.0]
const TDT_PSTOP = [3000000.0, 9000000.0, 3.0e11]
"Coefficient records (ticks, parallel time, rate) of the two clocks."
const TDB_COEF = [0.0 0.0 1.0; 1280000.0 5000.0 1.00001; 2560000.0 10050.0 0.99999;
    28160000.0 110049.0 1.0]
const TDT_COEF = [0.0 -5000.0 60.0; 1500000.0 25000.0 60.001; 3000000.0 55000.5 59.999]

const PATHS = Dict{Symbol,String}()
const DIR = Ref("")

"Directory of the generated kernels (made once per process)."
function dir()
    isdir(DIR[]) || (DIR[] = mktempdir())
    return DIR[]
end

fmt(v) = join((repr(x) for x in v), "\n      ")

"""
    sclk_kernel() -> String

Write (once per process) the SCLK kernel of the clocks `SC_TDB` and `SC_TDT` and return
its path.
"""
function sclk_kernel()
    haskey(PATHS, :sclk) && return PATHS[:sclk]
    path = joinpath(dir(), "test_clocks.tsc")
    write(path, """
    KPL/SCLK

    Two type 1 spacecraft clocks for the tests (written by test/ckfixtures.jl).

    \\begindata

    SCLK_KERNEL_ID            = ( @2026-10-11 )

    SCLK_DATA_TYPE_77         = ( 1 )
    SCLK01_TIME_SYSTEM_77     = ( 1 )
    SCLK01_N_FIELDS_77        = ( 2 )
    SCLK01_MODULI_77          = ( 4294967296 256 )
    SCLK01_OFFSETS_77         = ( 0 0 )
    SCLK01_OUTPUT_DELIM_77    = ( 1 )
    SCLK_PARTITION_START_77   = ( $(fmt(TDB_PSTART)) )
    SCLK_PARTITION_END_77     = ( $(fmt(TDB_PSTOP)) )
    SCLK01_COEFFICIENTS_77    = (
          $(fmt(vec(permutedims(TDB_COEF)))) )

    SCLK_DATA_TYPE_78         = ( 1 )
    SCLK01_TIME_SYSTEM_78     = ( 2 )
    SCLK01_N_FIELDS_78        = ( 3 )
    SCLK01_MODULI_78          = ( 100000000 60 50 )
    SCLK01_OFFSETS_78         = ( 1 0 0 )
    SCLK01_OUTPUT_DELIM_78    = ( 2 )
    SCLK_PARTITION_START_78   = ( $(fmt(TDT_PSTART)) )
    SCLK_PARTITION_END_78     = ( $(fmt(TDT_PSTOP)) )
    SCLK01_COEFFICIENTS_78    = (
          $(fmt(vec(permutedims(TDT_COEF)))) )

    \\begintext
    """)
    return PATHS[:sclk] = path
end

# --- frames ---------------------------------------------------------------------------------

"""
CK frames of the tests: name => (CK ID, clock). The frame ID is the CK ID. Frames with
the clock `nothing` have no `CK_<id>_SCLK` variable and use the default `id ÷ 1000`.
"""
const CK_FRAMES = [
    "TESTSC_T1" => (-77001, SC_TDB), "TESTSC_T1_NOAV" => (-77011, SC_TDB),
    "TESTSC_T2" => (-77002, SC_TDB), "TESTSC_T3" => (-77003, nothing),
    "TESTSC_T3_NOAV" => (-77013, nothing), "TESTSC_T3_EARTH" => (-77030, SC_TDB),
    "TESTSC_T4" => (-77004, SC_TDB),
    "TESTSC_T5S0" => (-77050, SC_TDB), "TESTSC_T5S1" => (-77051, SC_TDB),
    "TESTSC_T5S2" => (-77052, SC_TDB), "TESTSC_T5S3" => (-77053, nothing),
    "TESTSC_T6" => (-77006, SC_TDB), "TESTSC_TDT" => (-78100, nothing),
    "TESTSC_OTHERCLOCK" => (-99001, SC_TDB), "TESTSC_NODATA" => (-77099, SC_TDB),
]
"Fixed-offset frame defined relative to the CK frame `TESTSC_T3`."
const TK_FRAME = "TESTSC_CAM"

"Write (once per process) the frame kernel of the CK frames and return its path."
function frame_kernel()
    haskey(PATHS, :fk) && return PATHS[:fk]
    path = joinpath(dir(), "test_ck_frames.tf")
    io = IOBuffer()
    println(io, "KPL/FK\n\nCK frames for the tests (written by test/ckfixtures.jl).\n")
    println(io, "\\begindata\n")
    for (name, (id, clock)) in CK_FRAMES
        println(io, "   FRAME_$name = $id")
        println(io, "   FRAME_$(id)_NAME = '$name'")
        println(io, "   FRAME_$(id)_CLASS = 3")
        println(io, "   FRAME_$(id)_CLASS_ID = $id")
        println(io, "   FRAME_$(id)_CENTER = $(id ÷ 1000)")
        if clock !== nothing
            println(io, "   CK_$(id)_SCLK = $clock")
            println(io, "   CK_$(id)_SPK = $(id ÷ 1000)")
        end
        println(io)
    end
    println(io, """
       FRAME_$TK_FRAME = -77900
       FRAME_-77900_NAME = '$TK_FRAME'
       FRAME_-77900_CLASS = 4
       FRAME_-77900_CLASS_ID = -77900
       FRAME_-77900_CENTER = -77
       TKFRAME_-77900_SPEC = 'ANGLES'
       TKFRAME_-77900_RELATIVE = 'TESTSC_T3'
       TKFRAME_-77900_ANGLES = ( 10.0 -25.0 40.0 )
       TKFRAME_-77900_AXES = ( 3 1 2 )
       TKFRAME_-77900_UNITS = 'DEGREES'

    \\begintext
    """)
    write(path, take!(io))
    return PATHS[:fk] = path
end

# --- attitude profile -----------------------------------------------------------------------

"""
    quat(t, k) -> NTuple{4}

Smooth unit quaternion (SPICE style, scalar first) at `t` seconds for the profile
number `k`; rotation rates are about 2e-3 rad/s. Generic over the number type of `t`.
"""
function quat(t, k::Integer)
    theta = 0.7 + 2e-3 * (1 + 0.1 * k) * t + 0.3 * sin(3e-4 * t)
    a1 = sin(1e-4 * t + k)
    a2 = cos(2e-4 * t)
    a3 = 0.5 + 0.3 * sin(1.5e-4 * t)
    n = sqrt(a1 * a1 + a2 * a2 + a3 * a3)
    s, c = sincos(theta / 2)
    return (c, s * a1 / n, s * a2 / n, s * a3 / n)
end

"""
    attitude(t, k) -> (q, dq, av, dav)

Quaternion, its time derivative (per second), the angular velocity that CSPICE
associates with them (`qdq2av`: rad/s, in the base frame), and the time derivative of
the angular velocity (by a central difference), all as `Vector{Float64}`.
"""
function attitude(t::Float64, k::Integer)
    jets = map(derivatives, quat(taylor_var(t, Val(1)), k))
    q = Float64[j[1] for j in jets]
    dq = Float64[j[2] for j in jets]
    av = collect(SPICE.qdq2av(q, dq))
    h = 1e-2
    avat(x) = (j = map(derivatives, quat(taylor_var(x, Val(1)), k));
        collect(SPICE.qdq2av(Float64[a[1] for a in j], Float64[a[2] for a in j])))
    dav = (avat(t + h) - avat(t - h)) / (2h)
    return q, dq, av, dav
end

"Packet of subtype `subtyp` of CK types 5 and 6 at `t` seconds for profile `k`."
function packet(t::Float64, k::Integer, subtyp::Integer)
    q, dq, av, dav = attitude(t, k)
    return subtyp == 0 ? vcat(q, dq) : subtyp == 1 ? q : subtyp == 2 ? vcat(q, dq, av, dav) :
           vcat(q, av)
end

# --- time tags --------------------------------------------------------------------------------

"Number of pointing instances of the large segments (more than two directory groups)."
const NREC = 260

"""
    ticks_tdb() -> Vector{Float64}

The `NREC` time tags of the segments on clock `SC_TDB`: whole ticks, 20 to 60 nominal
seconds apart, from 200 s to past the partition boundary at 10000 s. The spacing is a
fixed pattern, so the kernels do not depend on the random seed.
"""
function ticks_tdb()
    t = Float64[]
    sec = 200
    for i in 1:NREC
        push!(t, sec * TICKS_PER_SEC)
        sec += 20 + (i * 37) % 41
    end
    return t
end

"Nominal seconds of the tick value `tick` of clock `SC_TDB`."
seconds(tick::Real) = tick * RATE

# --- thin wrappers of the CSPICE writers ------------------------------------------------------

"Flat C array (record after record) of the vectors `v`."
flat(v) = Float64[x for rec in v for x in rec]

function ckw01(h, t0, t1, inst, ref, segid, sclkdp, quats, avvs)
    avflag = avvs !== nothing
    ccall((:ckw01_c, LIB), Cvoid, (Cint, Cdouble, Cdouble, Cint, Cstring, Cint, Cstring, Cint,
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        h, t0, t1, inst, ref, avflag, segid, length(sclkdp), sclkdp, flat(quats),
        avflag ? flat(avvs) : zeros(3))
    SPICE.handleerror()
end

function ckw02(h, t0, t1, inst, ref, segid, start, stop, quats, avvs, rates)
    ccall((:ckw02_c, LIB), Cvoid, (Cint, Cdouble, Cdouble, Cint, Cstring, Cstring, Cint,
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
        h, t0, t1, inst, ref, segid, length(start), start, stop, flat(quats), flat(avvs), rates)
    SPICE.handleerror()
end

function ckw03(h, t0, t1, inst, ref, segid, sclkdp, quats, avvs, starts)
    avflag = avvs !== nothing
    ccall((:ckw03_c, LIB), Cvoid, (Cint, Cdouble, Cdouble, Cint, Cstring, Cint, Cstring, Cint,
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}),
        h, t0, t1, inst, ref, avflag, segid, length(sclkdp), sclkdp, flat(quats),
        avflag ? flat(avvs) : zeros(3), length(starts), starts)
    SPICE.handleerror()
end

function ckw05(h, subtyp, degree, t0, t1, inst, ref, avflag, segid, sclkdp, packts, rate, starts)
    ccall((:ckw05_c, LIB), Cvoid, (Cint, Cint, Cint, Cdouble, Cdouble, Cint, Cstring, Cint,
            Cstring, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cint, Ptr{Cdouble}),
        h, subtyp, degree, t0, t1, inst, ref, avflag, segid, length(sclkdp), sclkdp,
        flat(packts), rate, length(starts), starts)
    SPICE.handleerror()
end

"""
    ckw04(h, t0, t1, inst, ref, avflag, segid, starts, records)

Write a type 4 segment with the Fortran-style writers of CSPICE. `records[i]` is one
record in the layout of `ckw04a_` (midpoint, radius, seven coefficient counts, the
coefficients) and `starts[i]` the start of its interval.
"""
function ckw04(h, t0, t1, inst, ref, avflag, segid, starts, records)
    ccall((:ckw04b_, LIB), Cint, (Ref{Cint}, Ref{Cdouble}, Ref{Cint}, Ptr{UInt8}, Ref{Cint},
            Ptr{UInt8}, Cint, Cint),
        h, t0, inst, ref, avflag, segid, ncodeunits(ref), ncodeunits(segid))
    SPICE.handleerror()
    # One call for each record: with several records in one call `ckw04a_` packs the
    # array in place and, for the first records after the first one, overwrites the
    # midpoint and the radius of a record before it moves them.
    for (rec, start) in zip(records, starts)
        ccall((:ckw04a_, LIB), Cint, (Ref{Cint}, Ref{Cint}, Ptr{Cint}, Ptr{Cdouble},
                Ptr{Cdouble}), h, 1, Cint[length(rec)], copy(rec), Float64[start])
        SPICE.handleerror()
    end
    ccall((:ckw04e_, LIB), Cint, (Ref{Cint}, Ref{Cdouble}), h, t1)
    SPICE.handleerror()
end

"""
    ckw06(h, inst, ref, avflag, t0, t1, segid, minis, bounds, sellst)

Write a type 6 segment with `ckw06_`. `minis[i]` is a named tuple `(subtyp, degree,
rate, epochs, packets)` of one mini-segment and `bounds` holds the `length(minis) + 1`
interval bounds.
"""
function ckw06(h, inst, ref, avflag, t0, t1, segid, minis, bounds, sellst)
    ccall((:ckw06_, LIB), Cint, (Ref{Cint}, Ref{Cint}, Ptr{UInt8}, Ref{Cint}, Ref{Cdouble},
            Ref{Cdouble}, Ptr{UInt8}, Ref{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble},
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}, Cint, Cint),
        h, inst, ref, avflag, t0, t1, segid, length(minis),
        Cint[length(m.epochs) for m in minis], Cint[m.subtyp for m in minis],
        Cint[m.degree for m in minis], Float64[x for m in minis for p in m.packets for x in p],
        Float64[m.rate for m in minis], Float64[x for m in minis for x in m.epochs],
        bounds, sellst, ncodeunits(ref), ncodeunits(segid))
    SPICE.handleerror()
end

"Open a new CK file `name` in the fixture directory, call `body(handle)`, and close it."
function write_ck(body, name::AbstractString)
    path = joinpath(dir(), name)
    isfile(path) && rm(path)
    h = SPICE.ckopn(path, name, 0)
    try
        body(h)
    finally
        SPICE.ckcls(h)
    end
    return path
end

# --- the CK files -----------------------------------------------------------------------------

"Interpolation interval starts of the type 3 and type 5 segments: indices into the time tags."
const START_INDEX = [1, 40, 41, 100, 101, 180, 230]

"""
    chebyshev_record(ta, tb, k, ncf) -> Vector{Float64}

Type 4 record for the clock interval `[ta, tb]` (ticks of `SC_TDB`): Chebyshev
coefficients of the four quaternion components and the three angular velocity
components of profile `k`, `ncf[i]` coefficients for component `i`, from the values at
the Chebyshev nodes.
"""
function chebyshev_record(ta::Float64, tb::Float64, k::Integer, ncf::NTuple{7,Int})
    mid, rad = (ta + tb) / 2, (tb - ta) / 2
    rec = Float64[mid, rad, ncf...]
    for comp in 1:7
        n = ncf[comp]
        f = Float64[]
        for j in 0:n-1
            x = cos(pi * (j + 0.5) / n)
            q, _, av, _ = attitude(seconds(mid + rad * x), k)
            push!(f, comp <= 4 ? q[comp] : av[comp-4])
        end
        for m in 0:n-1
            c = sum(f[j+1] * cos(pi * m * (j + 0.5) / n) for j in 0:n-1) * 2 / n
            push!(rec, m == 0 ? c / 2 : c)
        end
    end
    return rec
end

"Clock intervals of the type 4 segment: 30 intervals of 300 s, a gap of 20 s after every third one."
function type4_intervals()
    out = NTuple{2,Float64}[]
    ta = ticks_tdb()[1]
    for i in 1:30
        tb = ta + 300 * TICKS_PER_SEC
        push!(out, (ta, tb))
        ta = tb + (i % 3 == 0 ? 20 * TICKS_PER_SEC : 0.0)
    end
    return out
end

"Clock intervals of the type 2 segment: every second tag starts one that ends 5 s before the next tag."
function type2_intervals()
    tk = ticks_tdb()
    return [(tk[i], tk[i+1] - 5 * TICKS_PER_SEC) for i in 1:2:length(tk)-1]
end

"""
    ck_files() -> Vector{String}

Write (once per process) the CK files of the tests and return their paths in load
order. Every segment is relative to J2000 unless its name says otherwise.
"""
function ck_files()
    haskey(PATHS, :ck1) && return [PATHS[:ck1], PATHS[:ck2], PATHS[:ck3]]
    tk = ticks_tdb()
    att = [attitude(seconds(t), 0) for t in tk]
    quats = [a[1] for a in att]
    avvs = [a[3] for a in att]
    starts = tk[START_INDEX]

    # File 1: types 1, 2 and 3.
    PATHS[:ck1] = write_ck("test_types123.bc") do h
        ckw01(h, tk[1], tk[end], -77001, "J2000", "type 1 with av", tk, quats, avvs)
        ckw01(h, tk[1], tk[end], -77011, "J2000", "type 1 without av", tk, quats, nothing)
        # Type 2: every second tag starts an interval that ends 5 s before the next tag.
        s2 = tk[1:2:end-1]
        e2 = tk[2:2:end] .- 5 * TICKS_PER_SEC
        ckw02(h, s2[1], e2[end], -77002, "J2000", "type 2", s2, e2, quats[1:2:end-1],
            avvs[1:2:end-1], fill(RATE, length(s2)))
        ckw03(h, tk[1], tk[end], -77003, "J2000", "type 3 with av", tk, quats, avvs, starts)
        ckw03(h, tk[1], tk[end], -77013, "J2000", "type 3 without av", tk, quats, nothing, starts)
        # Relative to a body-fixed frame, and for the instrument with another clock ID.
        ckw03(h, tk[1], tk[end], -77030, "IAU_EARTH", "type 3 in IAU_EARTH", tk, quats, avvs, starts)
        ckw03(h, tk[1], tk[end], -99001, "ECLIPJ2000", "type 3, clock by CK_*_SCLK", tk, quats,
            avvs, starts)
        # On the TDT clock: 80 tags, one count (3000 ticks, about 60 s) apart.
        tt = [3000.0 * (10 + i) for i in 1:80]
        at = [attitude(TDT_RATE * t, 2) for t in tt]
        ckw03(h, tt[1], tt[end], -78100, "J2000", "type 3 on the TDT clock", tt,
            [a[1] for a in at], [a[3] for a in at], [tt[1], tt[30]])
    end

    # File 2: types 4, 5 and 6.
    PATHS[:ck2] = write_ck("test_types456.bc") do h
        iv = type4_intervals()
        recs = [chebyshev_record(a, b, 0, (9, 8, 9, 10, 7, 7, 6)) for (a, b) in iv]
        ckw04(h, iv[1][1], iv[end][2], -77004, "J2000", true, "type 4", first.(iv), recs)
        for (subtyp, degree) in ((0, 7), (1, 9), (2, 3), (3, 5))
            ckw05(h, subtyp, degree, tk[1], tk[end], -77050 - subtyp, "J2000", true,
                "type 5 subtype $subtyp", tk, [packet(seconds(t), 0, subtyp) for t in tk],
                RATE, starts)
        end
        # Type 6: four mini-segments, one of each subtype; the second and the fourth end
        # before their interval does (a gap), and the rates differ.
        bounds = [tk[1], tk[60], tk[130], tk[200], tk[260]]
        spec = ((0, 7, 1:60, RATE), (1, 9, 60:120, RATE), (2, 3, 130:200, RATE * 1.0001),
            (3, 5, 200:250, RATE))
        minis = [(subtyp=s, degree=d, rate=r, epochs=tk[idx],
            packets=[packet(seconds(t), 0, s) for t in tk[idx]]) for (s, d, idx, r) in spec]
        ckw06(h, -77006, "J2000", true, tk[1], tk[260], "type 6", minis, bounds, false)
    end

    # File 3, loaded last: a second attitude profile for part of the range of the type 3
    # instrument, so that the priority rule shows; first a segment without angular
    # velocity, then (higher priority) one with it on a shorter range.
    PATHS[:ck3] = write_ck("test_priority.bc") do h
        r1 = 120:170
        a1 = [attitude(seconds(t), 5) for t in tk[r1]]
        ckw03(h, tk[120], tk[170], -77003, "J2000", "priority, without av", tk[r1],
            [a[1] for a in a1], nothing, [tk[120]])
        r2 = 140:160
        a2 = [attitude(seconds(t), 7) for t in tk[r2]]
        ckw03(h, tk[140], tk[160], -77003, "J2000", "priority, with av", tk[r2],
            [a[1] for a in a2], [a[3] for a in a2], [tk[140]])
    end
    return [PATHS[:ck1], PATHS[:ck2], PATHS[:ck3]]
end

"Allocation that `f(args...)` reports on its second call (the first call compiles)."
function warm(f, args...)
    f(args...)
    return f(args...)
end

"All generated kernels in load order: SCLK kernel, frame kernel, CK files."
all_kernels() = String[sclk_kernel(), frame_kernel(), ck_files()...]

end # module
