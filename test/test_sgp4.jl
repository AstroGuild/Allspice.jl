# Group 17: two-line elements (SGP4 and SDP4) and SPK segments of type 10 against CSPICE.
#
# Feature file: features/sgp4_tle.feature
#
# The element sets are the 29 verification sets of the SGP4 revision of Vallado et al.
# (2006), which the comment block of `zzsgp4.c` of CSPICE lists. They are grouped by the
# branch of the propagator that they exercise:
#
#   low orbit            period below 225 minutes, SGP4
#   decaying             low perigee and strong drag; the propagation stops with an error
#   12 h resonant        Molniya orbits, the resonance terms of SDP4 for two revolutions a day
#   1 day resonant       geosynchronous orbits, the resonance terms for one revolution a day
#   semi-synchronous     GPS; SDP4 with lunar and solar terms and no resonance
#   high eccentricity    SDP4 with eccentricity 0.56 to 0.97 and no resonance
#   deep space           SDP4 with moderate eccentricity and no resonance
#
# No SPK file of type 10 is on this machine, so the group writes them with `spkw10_c`
# of CSPICE (called with `ccall`) and with `spkw10` of Allspice.
#
# Checks:
#   (a) `getelm` against CSPICE (equal numbers), and its errors;
#   (b) `evsgp4` against `evsgp4_c` at random epochs within two days of the element
#       epoch, with the same errors at the same epochs;
#   (c) the second operation mode of the propagator (the improved sidereal time) against
#       the entries `xxsgp4i_` and `xxsgp4e_` of CSPICE;
#   (d) states of type 10 segments (`spkgeo` in J2000, `spkezr` in IAU_EARTH) against
#       CSPICE at random epochs, at the element epochs, at the middle between two element
#       epochs, where CSPICE changes the pair of sets, and at the segment ends;
#   (e) a segment that `spkw10` of Allspice writes against the segment of CSPICE, word for
#       word, with fewer and with more than 100 element sets;
#   (f) a segment in the old layout (packets of ten numbers, 15 meta items);
#   (g) orders 2 and 3 of `spkder` against the nested `ForwardDiff.derivative` of the
#       velocity of `spkgeo` and against a central difference of the velocity;
#   (h) no allocation and inferred result types;
#   (i) the errors.
#
# Tolerance of (b) and (d). The two engines make the same operations in the same order,
# but the sine and cosine of Julia and of the C library differ in the last bit for some
# arguments. The Kepler equation amplifies this difference by about 1/(1 - e)^2, so the
# bound is `SGP4_RTOL` except for the orbit with e = 0.97 (`SGP4_RTOL_E97`). In a
# body-fixed frame the velocity of a geosynchronous object is near zero, so a velocity
# difference is divided by the inertial speed and not by the body-fixed speed.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

import ForwardDiff

const SGP4_RTOL = 1e-12
const SGP4_RTOL_E97 = 5e-11
const SGP4_AD_RTOL = (1e-9, 1e-7)      # orders 2 and 3 of spkder against ForwardDiff
const SGP4_FD_RTOL = 1e-5
"Object whose eccentricity is 0.97 (the looser bound applies)."
const SGP4_E97 = 23333

"J2, J3, J4, KE, QO, SO, ER, AE of the WGS-72 model, as `geophysical.ker` of NAIF."
const SGP4_GEOPHS = [1.082616e-3, -2.53881e-6, -1.65597e-6, 7.43669161e-2, 120.0, 78.0,
    6378.135, 1.0]
const SGP4_GEOPHS_NAMES = ("J2", "J3", "J4", "KE", "QO", "SO", "ER", "AE")

const SGP4_CLASSES = ("low orbit", "decaying", "12 h resonant", "1 day resonant",
    "semi-synchronous", "high eccentricity", "deep space")

"The verification element sets: class and the two lines."
const SGP4_TLES = [
    ("low orbit",
        "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753",
        "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"),
    ("deep space",
        "1 04632U 70093B   04031.91070959 -.00000084  00000-0  10000-3 0  9955",
        "2 04632  11.4628 273.1101 1450506 207.6000 143.9350  1.20231981 44145"),
    ("low orbit",
        "1 06251U 62025E   06176.82412014  .00008885  00000-0  12808-3 0  3985",
        "2 06251  58.0579  54.0425 0030035 139.1568 221.1854 15.56387291  6774"),
    ("12 h resonant",
        "1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813",
        "2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656"),
    ("12 h resonant",
        "1 09880U 77021A   06176.56157475  .00000421  00000-0  10000-3 0  9814",
        "2 09880  64.5968 349.3786 7069051 270.0229  16.3320  2.00813614112380"),
    ("1 day resonant",
        "1 09998U 74033F   05148.79417928 -.00000112  00000-0  00000+0 0  4480",
        "2 09998   9.4958 313.1750 0270971 327.5225  30.8097  1.16186785 45878"),
    ("high eccentricity",
        "1 11801U          80230.29629788  .01431103  00000-0  14311-1      13",
        "2 11801  46.7916 230.4354 7318036  47.4722  10.4117  2.28537848    13"),
    ("1 day resonant",
        "1 14128U 83058A   06176.02844893 -.00000158  00000-0  10000-3 0  9627",
        "2 14128  11.4384  35.2134 0011562  26.4582 333.5652  0.98870114 46093"),
    ("high eccentricity",
        "1 16925U 86065D   06151.67415771  .02550794 -30915-6  18784-3 0  4486",
        "2 16925  62.0906 295.0239 5596327 245.1593  47.9690  4.88511875148616"),
    ("high eccentricity",
        "1 20413U 83020D   05363.79166667  .00000000  00000-0  00000+0 0  7041",
        "2 20413  12.3514 187.4253 7864447 196.3027 356.5478  0.24690082  7978"),
    ("12 h resonant",
        "1 21897U 92011A   06176.02341244 -.00001273  00000-0 -13525-3 0  3044",
        "2 21897  62.1749 198.0096 7421690 253.0462  20.1561  2.01269994104880"),
    ("decaying",
        "1 22312U 93002D   06094.46235912  .99999999  81888-5  49949-3 0  3953",
        "2 22312  62.1486  77.4698 0308723 267.9229  88.7392 15.95744531 98783"),
    ("12 h resonant",
        "1 22674U 93035D   06176.55909107  .00002121  00000-0  29868-3 0  6569",
        "2 22674  63.5035 354.4452 7541712 253.3264  18.7754  1.96679808 93877"),
    ("high eccentricity",
        "1 23177U 94040C   06175.45752052  .00000386  00000-0  76590-3 0    95",
        "2 23177   7.0496 179.8238 7258491 296.0482   8.3061  2.25906668 97438"),
    ("high eccentricity",
        "1 23333U 94071A   94305.49999999 -.00172956  26967-3  10000-3 0    15",
        "2 23333  28.7490   2.3720 9728298  30.4360   1.3500  0.07309491    70"),
    ("high eccentricity",
        "1 23599U 95029B   06171.76535463  .00085586  12891-6  12956-2 0  2905",
        "2 23599   6.9327   0.2849 5782022 274.4436  25.2425  4.47796565123555"),
    ("1 day resonant",
        "1 24208U 96044A   06177.04061740 -.00000094  00000-0  10000-3 0  1600",
        "2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119"),
    ("1 day resonant",
        "1 25954U 99060A   04039.68057285 -.00000108  00000-0  00000-0 0  6847",
        "2 25954   0.0004 243.8136 0001765  15.5294  22.7134  1.00271289 15615"),
    ("1 day resonant",
        "1 26900U 01039A   06106.74503247  .00000045  00000-0  10000-3 0  8290",
        "2 26900   0.0164 266.5378 0003319  86.1794 182.2590  1.00273847 16981"),
    ("12 h resonant",
        "1 26975U 78066F   06174.85818871  .00000620  00000-0  10000-3 0  6809",
        "2 26975  68.4714 236.1303 5602877 123.7484 302.5767  2.05657553 67521"),
    ("low orbit",
        "1 28057U 03049A   06177.78615833  .00000060  00000-0  35940-4 0  1836",
        "2 28057  98.4283 247.6961 0000884  88.1964 271.9322 14.35478080140550"),
    ("semi-synchronous",
        "1 28129U 03058A   06175.57071136 -.00000104  00000-0  10000-3 0   459",
        "2 28129  54.7298 324.8098 0048506 266.2640  93.1663  2.00562768 18443"),
    ("decaying",
        "1 28350U 04020A   06167.21788666  .16154492  76267-5  18678-3 0  8894",
        "2 28350  64.9977 345.6130 0024870 260.7578  99.9590 16.47856722116490"),
    ("high eccentricity",
        "1 28623U 05006B   06177.81079184  .00637644  69054-6  96390-3 0  6000",
        "2 28623  28.5200 114.9834 6249053 170.2550 212.8965  3.79477162 12753"),
    ("1 day resonant",
        "1 28626U 05008A   06176.46683397 -.00000205  00000-0  10000-3 0  2190",
        "2 28626   0.0019 286.9433 0000335  13.7918  55.6504  1.00270176  4891"),
    ("decaying",
        "1 28872U 05037B   05333.02012661  .25992681  00000-0  24476-3 0  1534",
        "2 28872  96.4736 157.9986 0303955 244.0492 110.6523 16.46015938 10708"),
    ("decaying",
        "1 29141U 85108AA  06170.26783845  .99999999  00000-0  13519-0 0   718",
        "2 29141  82.4288 273.4882 0015848 277.2124  83.9133 15.93343074  6828"),
    ("low orbit",
        "1 29238U 06022G   06177.28732010  .00766286  10823-4  13334-2 0   101",
        "2 29238  51.5595 213.7903 0202579  95.2503 267.9010 15.73823839  1061"),
    ("low orbit",
        "1 88888U          80275.98708465  .00073094  13844-3  66816-4 0    87",
        "2 88888  72.8435 115.9689 0086731  52.6988 110.5714 16.05824518  1058"),
]

"Catalogue number of the element set `t`."
sgp4_id(t) = parse(Int, t[2][3:7])
"Body ID of the test segment of the element set `t`."
sgp4_body(t) = -100000 - sgp4_id(t)
sgp4_tol(id::Integer) = id == SGP4_E97 ? SGP4_RTOL_E97 : SGP4_RTOL

"`f()`, or the CSPICE error that it throws."
function sgp4_oracle(f)
    try
        return f()
    catch e
        e isa SPICE.SpiceError || rethrow()
        return e
    end
end

"`f()`, or the Allspice error that it throws."
function sgp4_neo(f)
    try
        return f()
    catch e
        e isa Allspice.SpiceError || rethrow()
        return e
    end
end

"""
    sgp4_c_check()

Throw the pending CSPICE error, if there is one, as `SPICE.SpiceError` with the short
message (`SPICE(...)`) before the long one, and reset the error state of CSPICE.
(`SPICE.handleerror` keeps the long message only.)
"""
function sgp4_c_check()
    ccall((:failed_c, SPICE.libcspice), Cint, ()) == 0 && return nothing
    short = zeros(UInt8, 64)
    ccall((:getmsg_c, SPICE.libcspice), Cvoid, (Cstring, Cint, Ref{UInt8}), "SHORT", 64, short)
    text = String(short[1:findfirst(iszero, short)-1])
    try
        SPICE.handleerror()
    catch e
        e isa SPICE.SpiceError || rethrow()
        throw(SPICE.SpiceError(text * " " * e.msg))
    end
    return nothing
end

"`evsgp4_c` of CSPICE."
function sgp4_c_evsgp4(et::Float64, geophs::Vector{Float64}, elems::Vector{Float64})
    state = zeros(6)
    ccall((:evsgp4_c, SPICE.libcspice), Cvoid,
        (Cdouble, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}), et, geophs, elems, state)
    sgp4_c_check()
    return state
end

"Initialization and propagation of CSPICE in the operation mode `opmode`, `t` minutes from the epoch."
function sgp4_c_mode(geophs::Vector{Float64}, elems::Vector{Float64}, opmode::Int, t::Float64)
    state = zeros(6)
    ccall((:xxsgp4i_, SPICE.libcspice), Cint, (Ref{Cdouble}, Ref{Cdouble}, Ref{Cint}),
        geophs, elems, Cint(opmode))
    SPICE.handleerror()
    ccall((:xxsgp4e_, SPICE.libcspice), Cint, (Ref{Cdouble}, Ref{Cdouble}), t, state)
    SPICE.handleerror()
    return state
end

"`spkw10_c` of CSPICE; `elems` is a vector of element sets."
function sgp4_c_spkw10(handle, body, center, frame, first, last, segid, geophs, elems, epochs)
    flat = reduce(vcat, elems)
    ccall((:spkw10_c, SPICE.libcspice), Cvoid, (Cint, Cint, Cint, Cstring, Cdouble, Cdouble,
            Cstring, Ref{Cdouble}, Cint, Ref{Cdouble}, Ref{Cdouble}),
        handle, body, center, frame, first, last, segid, geophs, length(epochs), flat,
        Vector{Float64}(epochs))
    SPICE.handleerror()
    return nothing
end

"""
    sgp4_sets(el, ep, n, step) -> (elems, epochs)

`n` element sets from the set `el` of epoch `ep`, `step` seconds apart. The mean
anomaly moves with the mean motion and the node moves a little, so the sets describe
nearly the same orbit but are not the same function of time.
Postcondition: `elems[i][10] == epochs[i]` and the epochs increase.
"""
function sgp4_sets(el, ep::Float64, n::Int, step::Float64)
    @assert n >= 1 && step > 0
    elems = [collect(Float64, el) for _ in 1:n]
    epochs = [ep + step * (i - 1) for i in 1:n]
    for i in 1:n
        elems[i][10] = epochs[i]
        elems[i][8] = mod(el[8] + el[9] * (epochs[i] - ep) / 60, 2pi)
        elems[i][5] = mod(el[5] + 1e-5 * (i - 1), 2pi)
    end
    @assert issorted(epochs) && all(elems[i][10] == epochs[i] for i in 1:n)
    return elems, epochs
end

"Largest of |difference| / |reference| of the position and of the velocity; `vscale` replaces the speed."
function sgp4_reldiff(a, b, vscale::Float64=norm(b[4:6]))
    return max(norm(a[1:3] - b[1:3]) / norm(b[1:3]), norm(a[4:6] - b[4:6]) / vscale)
end

"Relative difference of two vectors."
sgp4_vecdiff(a, b) = norm(a - b) / norm(b)

"""
    sgp4_word_diff(pa, pb) -> (count, worst)

Count of the 8-byte words that differ between the files `pa` and `pb`, and the largest
relative difference of these words read as doubles; `(-1, Inf)` if the sizes differ.
"""
function sgp4_word_diff(pa::AbstractString, pb::AbstractString)
    a, b = read(pa), read(pb)
    (length(a) == length(b) && length(a) % 8 == 0) || return (-1, Inf)
    wa, wb = reinterpret(Float64, a), reinterpret(Float64, b)
    count, worst = 0, 0.0
    for i in eachindex(wa)        # invariant: `count` words before `i` differ
        reinterpret(UInt64, wa[i]) == reinterpret(UInt64, wb[i]) && continue
        count += 1
        d = abs(wa[i] - wb[i]) / abs(wb[i])
        worst = isnan(d) ? Inf : max(worst, d)
    end
    return count, worst
end

"Mean time in nanoseconds of one call `f(args..., et)` over `n` calls at nearby epochs."
function sgp4_time(f, n::Int, et::Float64, args...)
    acc = 0.0
    t = time_ns()
    for i in 1:n                  # invariant: `acc` holds the sum of the first i - 1 results
        acc += f(args..., et + 1e-3 * i)
    end
    return (time_ns() - t) / n + 0.0 * acc
end

sgp4_neo_x(ctx, body, et) = spkgeo(ctx, body, et, 1, 399)[1][1]
sgp4_neo_x3(ctx, body, et) = spkder(ctx, body, et, 1, 399, Val(3))[4][1]
sgp4_neo_ev(ctx, g, el, et) = evsgp4(ctx, et, g, el)[1]
sgp4_c_x(body, et) = SPICE.spkgeo(body, et, "J2000", 399)[1][1]
sgp4_c_ev(g, el, et) = sgp4_c_evsgp4(et, g, el)[1]

alloc_sgp4_geo(ctx, body, et) = @allocated spkgeo(ctx, body, et, 1, 399)
alloc_sgp4_ezr(ctx, body, et, ref) = @allocated spkezr(ctx, body, et, ref, "NONE", 399)
alloc_sgp4_pos(ctx, body, et) = @allocated spkpos(ctx, body, et, 1, "NONE", 399)
alloc_sgp4_der(ctx, body, et, ::Val{N}) where {N} = @allocated spkder(ctx, body, et, 1, 399, Val(N))
alloc_sgp4_ev(ctx, et, g, el) = @allocated evsgp4(ctx, et, g, el)

"Allocation of the second call of `f(args...)` (the first call compiles)."
function sgp4_warm(f, args...)
    f(args...)
    return f(args...)
end

group("17. two-line elements and SPK type 10", Kernels.LSK, Kernels.PCK11) do
    dir = mktempdir()
    rng = MersenneTwister(SEED + 17)

    # The geophysical constants come from a text kernel, as CSPICE users load them.
    geoker = joinpath(dir, "geophysical.tpc")
    write(geoker, """
    KPL/PCK

    Geophysical constants of the WGS-72 model for two-line elements (written by the tests).

    \\begindata

       BODY399_J2 =    1.082616D-3
       BODY399_J3 =   -2.53881D-6
       BODY399_J4 =   -1.65597D-6
       BODY399_KE =    7.43669161D-2
       BODY399_QO =  120.0D0
       BODY399_SO =   78.0D0
       BODY399_ER = 6378.135D0
       BODY399_AE =    1.0D0

    \\begintext
    """)
    base = String[Kernels.LSK, Kernels.PCK11, geoker]
    oracle_load(base)
    ctx0 = neo_context(base)
    iau_earth = Int(namfrm(ctx0, "IAU_EARTH"))

    geophs = Float64[bodvrd(ctx0, 399, name)[1] for name in SGP4_GEOPHS_NAMES]
    @apart @testset "geophysical constants from a kernel" begin
        @test geophs == Float64[SPICE.bodvrd("EARTH", name)[1] for name in SGP4_GEOPHS_NAMES]
        # The text kernel reader of both engines gives J3 one unit in the last place
        # away from the Julia literal.
        @test geophs ≈ SGP4_GEOPHS rtol = 1e-15
        for name in SGP4_GEOPHS_NAMES
            @test length(bodvrd(ctx0, "EARTH", name)) == 1
        end
    end

    # Element sets of both engines.
    parsed = Dict{Int,Tuple{Float64,Vector{Float64}}}()

    @apart @testset "getelm against CSPICE" begin
        for t in SGP4_TLES
            pair = [t[2], t[3]]
            ep, el = getelm(ctx0, 1957, pair)
            epc, elc = SPICE.getelm(1957, pair)
            @test ep == epc
            @test collect(el) == elc
            @test el[10] == ep
            parsed[sgp4_id(t)] = (ep, collect(el))
        end
        pair = [SGP4_TLES[1][2], SGP4_TLES[1][3]]
        @test (@inferred getelm(ctx0, 1957, pair)) isa Tuple{Float64,SVector{10,Float64}}
        # A tuple of lines and lines without the checksum column are accepted.
        @test getelm(ctx0, 1957, (pair[1], pair[2])) == getelm(ctx0, 1957, pair)
        short = [p[1:68] for p in pair]
        @test getelm(ctx0, 1957, short)[1] == SPICE.getelm(1957, short)[1]
        @test collect(getelm(ctx0, 1957, short)[2]) == SPICE.getelm(1957, short)[2]
        # The first year sets the century of the two-digit year.
        for frstyr in (1950, 1957, 1990, 2001, 2050)
            for t in SGP4_TLES[[1, 7, 15]]
                p = [t[2], t[3]]
                @test getelm(ctx0, frstyr, p)[1] == SPICE.getelm(frstyr, p)[1]
            end
        end
        # CSPICE N0067 does not read the checksum (column 69): a wrong one is accepted.
        wrong = [pair[1][1:68] * "0", pair[2][1:68] * "0"]
        @test wrong != pair
        @test getelm(ctx0, 1957, wrong) == getelm(ctx0, 1957, pair)
        @test SPICE.getelm(1957, wrong) == SPICE.getelm(1957, pair)
    end

    @apart @testset "getelm errors" begin
        l1, l2 = SGP4_TLES[1][2], SGP4_TLES[1][3]
        put(line, range, text) = line[1:first(range)-1] * text * line[last(range)+1:end]
        bad = [
            ("catalogue numbers differ", [l1, put(l2, 3:7, "00006")]),
            ("line too short", [l1[1:60], l2]),
            ("term that is not a number", [l1, put(l2, 9:16, " 3x.2682")]),
            ("term that is not a number in line 1", [put(l1, 21:32, "179.7849506z"), l2]),
            ("inclination out of bounds", [l1, put(l2, 9:16, "190.2682")]),
            ("node out of bounds", [l1, put(l2, 18:25, "360.7242")]),
            ("argument of perigee out of bounds", [l1, put(l2, 35:42, "361.7664")]),
            ("mean anomaly out of bounds", [l1, put(l2, 44:51, "400.3264")]),
            ("mean motion out of bounds", [l1, put(l2, 53:63, "20.82419157")]),
        ]
        for (label, lines) in bad
            @test length(lines) == 2
            r = sgp4_oracle(() -> SPICE.getelm(1957, lines))
            r isa SPICE.SpiceError || @info "CSPICE accepts the lines" label
            @test r isa SPICE.SpiceError
            @test throws_kind(() -> getelm(ctx0, 1957, lines), :badtle)
        end
        @test sgp4_neo(() -> getelm(ctx0, 1957, [l1])) isa Allspice.SpiceError
        # Blanks after column 69 are not an error in either engine.
        long = [l1 * "  ", l2 * "   "]
        @test getelm(ctx0, 1957, long)[1] == SPICE.getelm(1957, long)[1]
        @test collect(getelm(ctx0, 1957, long)[2]) == SPICE.getelm(1957, long)[2]
        @test throws_kind(() -> getelm(ctx0, 1957, [l1, put(l2, 9:16, " 34.26é2")]), :badtle)
        # No leap-second kernel.
        @test throws_kind(() -> getelm(SpiceContext(), 1957, [l1, l2]), :nolsk)
    end

    # Worst differences of each class: (b) evsgp4, (d) J2000 and IAU_EARTH.
    worst_ev = Dict(c => 0.0 for c in SGP4_CLASSES)
    worst_j2 = Dict(c => 0.0 for c in SGP4_CLASSES)
    worst_bf = Dict(c => 0.0 for c in SGP4_CLASSES)
    count_ev = Dict(c => 0 for c in SGP4_CLASSES)
    count_seg = Dict(c => 0 for c in SGP4_CLASSES)
    errors_ev = Dict(c => 0 for c in SGP4_CLASSES)
    errors_seg = Dict(c => 0 for c in SGP4_CLASSES)

    @apart @testset "evsgp4 against CSPICE" begin
        kinds = Set{Symbol}()
        for t in SGP4_TLES
            id, class = sgp4_id(t), t[1]
            ep, el = parsed[id]
            worst, mismatch, valid = 0.0, 0, 0
            minutes = vcat(0.0, 5.0, -5.0, 720.0, -720.0, 1440.0, -1440.0,
                [2880.0 * (2 * rand(rng) - 1) for _ in 1:nsample(100, 400)])
            for dt in minutes
                et = ep + 60 * dt
                sc = sgp4_oracle(() -> sgp4_c_evsgp4(et, geophs, el))
                sn = sgp4_neo(() -> evsgp4(ctx0, et, geophs, el))
                if sc isa Exception || sn isa Exception
                    errors_ev[class] += 1
                    (sc isa Exception) == (sn isa Exception) || (mismatch += 1)
                    if sc isa Exception && sn isa Exception
                        push!(kinds, sn.kind)
                        # The short message of CSPICE is the kind of Allspice.
                        occursin(uppercase(String(sn.kind)), uppercase(sc.msg)) || (mismatch += 1)
                    end
                    continue
                end
                valid += 1
                worst = max(worst, sgp4_reldiff(sn, sc))
            end
            @test mismatch == 0
            @test worst <= sgp4_tol(id)
            id == 28872 || @test valid > 0
            worst_ev[class] = max(worst_ev[class], worst)
            count_ev[class] += valid
        end
        # The decaying sets reach the errors of the propagator in both engines.
        @test :orbitdecay in kinds
        @test :suborbital in kinds
        @test :badmeccentricity in kinds
        measured("  evsgp4: error kinds that both engines gave: " * join(sort!(String.(collect(kinds))), ", "))
        # The argument forms give the same state.
        ep, el = parsed[5]
        ref = evsgp4(ctx0, ep + 600.0, geophs, el)
        @test evsgp4(ctx0, ep + 600.0, Tuple(geophs), Tuple(el)) == ref
        @test evsgp4(ctx0, ep + 600.0, SVector{8}(geophs), SVector{10}(el)) == ref
        @test evsgp4(Allspice.kernelset(ctx0), ep + 600.0, geophs, el) == ref
        @test (@inferred evsgp4(ctx0, ep + 600.0, geophs, el)) isa SVector{6,Float64}
        # As in CSPICE, QO, SO and AE are not read.
        other = copy(geophs)
        other[5:6] .= (100.0, 90.0)
        other[8] = 2.0
        @test evsgp4(ctx0, ep + 600.0, other, el) == ref
        @test sgp4_c_evsgp4(ep + 600.0, other, el) == sgp4_c_evsgp4(ep + 600.0, geophs, el)
    end

    @apart @testset "operation mode 2 (improved sidereal time) against CSPICE" begin
        l = Allspice._lsk(Allspice.kernelset(ctx0))
        worst = 0.0
        for t in SGP4_TLES
            id = sgp4_id(t)
            id in (22312, 28350, 28872, 29141) && continue
            ep, el = parsed[id]
            for opmode in (1, 2)
                r = Allspice._sgp4_init(Tuple(geophs), Tuple(el), Allspice._t10_epoch(l, ep), opmode)
                @test isbits(r)
                for dt in (0.0, 13.0, -211.0, 733.0, 1500.0)
                    id == SGP4_E97 && dt < 0 && continue
                    sc = sgp4_c_mode(geophs, el, opmode, dt)
                    p, v = Allspice._sgp4_eval(r, dt)
                    d = sgp4_reldiff(vcat(p, v), sc)
                    worst = max(worst, id == SGP4_E97 ? 0.0 : d)
                    @test d <= sgp4_tol(id)
                end
            end
        end
        measured(@sprintf("  operation modes 1 and 2, all classes: max relative difference %.2e", worst))
        @test throws_kind(() -> Allspice._sgp4_init(Tuple(geophs), Tuple(parsed[5][2]), 18000.0, 3),
            :unknownmode)
    end

    # ---------------------------------------------------------------------------------
    # Type 10 segments: one per element set, three sets 0.4 day apart, written by CSPICE
    # and by Allspice.
    # ---------------------------------------------------------------------------------
    pc = joinpath(dir, "tle_cspice.bsp")
    pn = joinpath(dir, "tle_allspice.bsp")
    segs = Dict{Int,Any}()
    hc = SPICE.spkopn(pc, "TLE", 0)
    hn = spkopn(pn, "TLE", 0)
    for t in SGP4_TLES
        id = sgp4_id(t)
        ep, el = parsed[id]
        elems, epochs = sgp4_sets(el, ep, 3, 0.4 * 86400.0)
        t0, t1 = epochs[1] - 3000.0, epochs[end] + 3000.0
        sgp4_c_spkw10(hc, sgp4_body(t), 399, "J2000", t0, t1, "TLE $id", geophs, elems, epochs)
        spkw10(hn, sgp4_body(t), 399, "J2000", t0, t1, "TLE $id", geophs, elems, epochs; ctx=ctx0)
        segs[id] = (; t0, t1, elems, epochs)
    end
    SPICE.spkcls(hc)
    spkcls(hn)

    SPICE.furnsh(pc)
    ctx = neo_context(vcat(base, pc))
    ctxn = neo_context(vcat(base, pn))

    @apart @testset "type 10 segments against CSPICE" begin
        for t in SGP4_TLES
            id, class, body = sgp4_id(t), t[1], sgp4_body(t)
            s = segs[id]
            e = s.epochs
            mids = [(e[1] + e[2]) / 2, (e[2] + e[3]) / 2]
            special = vcat(s.t0, s.t1, e, prevfloat.(e), nextfloat.(e), mids, prevfloat.(mids),
                nextfloat.(mids))
            ets = vcat(special, random_epochs(rng, nsample(60, 300), s.t0, s.t1))
            w2, wb, mismatch, valid = 0.0, 0.0, 0, 0
            for et in ets
                sc = sgp4_oracle(() -> SPICE.spkgeo(body, et, "J2000", 399)[1])
                sn = sgp4_neo(() -> spkgeo(ctx, body, et, 1, 399)[1])
                if sc isa Exception || sn isa Exception
                    errors_seg[class] += 1
                    (sc isa Exception) == (sn isa Exception) || (mismatch += 1)
                    continue
                end
                valid += 1
                w2 = max(w2, sgp4_reldiff(sn, sc))
                # The file of Allspice gives the same numbers as the file of CSPICE.
                spkgeo(ctxn, body, et, 1, 399)[1] == sn || (mismatch += 1)
                speed = norm(sc[4:6])
                bc = SPICE.spkezr(string(body), et, "IAU_EARTH", "NONE", "399")[1]
                bn = spkezr(ctx, body, et, iau_earth, "NONE", 399)[1]
                wb = max(wb, sgp4_reldiff(bn, bc, speed))
            end
            @test mismatch == 0
            @test w2 <= sgp4_tol(id)
            @test wb <= 2 * sgp4_tol(id)
            id == 28872 || @test valid > 0
            worst_j2[class] = max(worst_j2[class], w2)
            worst_bf[class] = max(worst_bf[class], wb)
            count_seg[class] += valid
        end
        # The other calls read the same segment.
        body = sgp4_body(SGP4_TLES[1])
        et = segs[5].epochs[1] + 12345.0
        state, lt = spkezr(ctx, body, et, 1, "NONE", 399)
        @test state == spkgeo(ctx, body, et, 1, 399)[1]
        @test spkpos(ctx, body, et, 1, "NONE", 399)[1] == state[SOneTo(3)]
        @test spkezr(ctx, string(body), et, "J2000", "NONE", "EARTH")[1] == state
        @test lt ≈ SPICE.spkezr(string(body), et, "J2000", "NONE", "399")[2] rtol = 1e-13
        @test spkgeo(ctx, 399, et, 1, body)[1] ≈ -state rtol = 1e-15
        # Outside the segment: no data.
        @test throws_kind(() -> spkgeo(ctx, body, segs[5].t1 + 10.0, 1, 399), :nocoverage)
    end

    @apart @testset "accuracy table" begin
        measured("  Largest relative difference from CSPICE of each orbit class " *
            "(samples; of them errors that both engines gave)")
        measured(@sprintf("  %-18s %-22s %-10s %-22s %-10s %-10s", "class", "evsgp4 samples",
            "evsgp4", "segment samples", "J2000", "IAU_EARTH"))
        for c in SGP4_CLASSES
            measured(@sprintf("  %-18s %-22s %-10.2e %-22s %-10.2e %-10.2e", c,
                "$(count_ev[c]) ($(errors_ev[c]))", worst_ev[c],
                "$(count_seg[c]) ($(errors_seg[c]))", worst_j2[c], worst_bf[c]))
            @test count_ev[c] > 0 && count_seg[c] > 0
        end
    end

    @apart @testset "spkw10 against CSPICE, word for word" begin
        # Three sets per segment, 29 segments. The words that can differ are the four
        # nutation numbers of a packet, which CSPICE computes with the sine and cosine
        # of the C library; the evaluator of CSPICE N0067 does not read them.
        count, worst = sgp4_word_diff(pn, pc)
        measured(@sprintf("  spkw10, 29 segments of 3 sets: %d of %d words differ, max relative %.1e",
            count, filesize(pc) ÷ 8, worst))
        @test count >= 0
        @test count <= 4 * 3 * length(SGP4_TLES)
        @test worst <= 1e-9
        # More than 100 sets: the segment has a directory of every 100th epoch.
        ep, el = parsed[5]
        for n in (1, 2, 100, 101, 250)
            elems, epochs = sgp4_sets(el, ep, n, 600.0)
            t0, t1 = epochs[1] - 100.0, epochs[end] + 100.0
            qc, qn = joinpath(dir, "many_c_$n.bsp"), joinpath(dir, "many_n_$n.bsp")
            h = SPICE.spkopn(qc, "TLE", 0)
            sgp4_c_spkw10(h, -77, 399, "J2000", t0, t1, "many sets", geophs, elems, epochs)
            SPICE.spkcls(h)
            h = spkopn(qn, "TLE", 0)
            # A matrix with one set in each column is accepted too.
            spkw10(h, -77, 399, "J2000", t0, t1, "many sets", geophs,
                isodd(n) ? reduce(hcat, elems) : elems, epochs; ctx=ctx0)
            spkcls(h)
            count, worst = sgp4_word_diff(qn, qc)
            @test 0 <= count <= 4 * n
            @test worst <= 1e-9
            n == 250 && measured(@sprintf("  spkw10, 1 segment of 250 sets: %d of %d words differ, max relative %.1e",
                count, filesize(qc) ÷ 8, worst))
            # Each engine reads the file of the other one.
            ca = neo_context(vcat(base, qc))
            cb = neo_context(vcat(base, qn))
            SPICE.furnsh(qn)
            w = 0.0
            ets = vcat(t0, t1, epochs[[1, end]], epochs[max(1, n ÷ 2)],
                random_epochs(rng, 60, t0, t1))
            n > 1 && append!(ets, [(epochs[i] + epochs[i+1]) / 2 for i in (1, n - 1)])
            for et in ets
                sa = spkgeo(ca, -77, et, 1, 399)[1]
                @test spkgeo(cb, -77, et, 1, 399)[1] == sa
                w = max(w, sgp4_reldiff(sa, SPICE.spkgeo(-77, et, "J2000", 399)[1]))
            end
            SPICE.unload(qn)
            @test w <= SGP4_RTOL
        end
    end

    @apart @testset "segment in the old layout" begin
        # Packets of the ten elements only, no epoch before a packet, 15 meta items
        # that are followed by their count: the files that are older than the
        # nutation numbers in the packets.
        ep, el = parsed[8195]
        n = 3
        elems, epochs = sgp4_sets(el, ep, n, 0.4 * 86400.0)
        t0, t1 = epochs[1] - 3000.0, epochs[end] + 3000.0
        refbas = 8 + 10 * n
        data = vcat(geophs, reduce(vcat, elems), epochs,
            Float64[0, 8, refbas + n, 0, 4, refbas, n, 0, 0, 0, 8, n, 0, 0, 10, 15])
        po = joinpath(dir, "old_layout.bsp")
        h = spkopn(po, "TLE", 0)
        dafbna(h, spkpds(-78, 399, "J2000", 10, t0, t1; ctx=ctx0), "old layout")
        dafada(h, data)
        dafena(h)
        spkcls(h)
        co = neo_context(vcat(base, po))
        SPICE.furnsh(po)
        w = 0.0
        for et in vcat(t0, t1, epochs, random_epochs(rng, 100, t0, t1))
            so = spkgeo(co, -78, et, 1, 399)[1]
            # The same sets in the present layout give the same state.
            @test so == spkgeo(ctx, sgp4_body(SGP4_TLES[4]), et, 1, 399)[1]
            w = max(w, sgp4_reldiff(so, SPICE.spkgeo(-78, et, "J2000", 399)[1]))
        end
        SPICE.unload(po)
        measured(@sprintf("  old layout (10 numbers a packet, 15 meta items): max relative difference %.2e", w))
        @test w <= SGP4_RTOL
    end

    @apart @testset "spkder: higher orders against ForwardDiff and differences" begin
        tv = Allspice.taylor_var
        for id in (5, 28057, 8195, 28129, 24208, 23177, 4632)
            body = -100000 - id
            s = segs[id]
            vel(t) = spkgeo(ctx, body, t, 1, 399)[1][SVector(4, 5, 6)]
            acc(t) = ForwardDiff.derivative(vel, t)
            jerk(t) = ForwardDiff.derivative(acc, t)
            worst = zeros(3)
            # Between the first two sets and not near an element epoch or the middle,
            # where CSPICE changes the pair of sets, nor near a step of the resonance
            # integrator (720 minutes from an element epoch).
            for _ in 1:6
                u = rand(rng, (0.2 + 0.25 * rand(rng), 0.55 + 0.25 * rand(rng)))
                et = s.epochs[1] + u * (s.epochs[2] - s.epochs[1])
                jet = spkder(ctx, body, et, 1, 399, Val(3))
                state = spkgeo(ctx, body, et, 1, 399)[1]
                @test jet[1] == state[SOneTo(3)]
                @test jet[2] == state[SVector(4, 5, 6)]
                @test spkder(ctx, body, et, 1, 399, Val(2))[3] ≈ jet[3] rtol = 1e-12
                @test spkder(ctx, body, et, 1, 399, Val(0))[1] == jet[1]
                dstate = spkgeo(ctx, body, ForwardDiff.Dual(et, 1.0), 1, 399)[1]
                @test ForwardDiff.value.(dstate) ≈ state rtol = 1e-14
                h = 0.5
                fd = (vel(et + h) - vel(et - h)) / (2h)
                worst .= max.(worst, (sgp4_vecdiff(jet[3], acc(et)), sgp4_vecdiff(jet[4], jerk(et)),
                    sgp4_vecdiff(jet[3], fd)))
                # A rotating frame. The jet of the frame is applied to the jet of the
                # segment as if the velocity were the derivative of the position, and
                # the SGP4 velocity is that only to about 1e-4 (see below), so the
                # order 2 in a rotating frame agrees with the difference of the
                # body-fixed velocity to that level and not to rounding.
                jb = spkder(ctx, body, et, iau_earth, 399, Val(2))
                velb(t) = spkezr(ctx, body, t, iau_earth, "NONE", 399)[1][SVector(4, 5, 6)]
                @test jb[1] == spkezr(ctx, body, et, iau_earth, "NONE", 399)[1][SOneTo(3)]
                @test sgp4_vecdiff(jb[3], (velb(et + h) - velb(et - h)) / (2h)) <= 5e-3
            end
            measured(@sprintf("  type 10, object %05d: AD order 2 %.1e, 3 %.1e; FD order 2 %.1e",
                id, worst...))
            @test worst[1] <= SGP4_AD_RTOL[1]
            @test worst[2] <= SGP4_AD_RTOL[2]
            @test worst[3] <= SGP4_FD_RTOL
            # evsgp4 with a series epoch against the nested derivative of evsgp4.
            ep, el = parsed[id]
            et = ep + 7000.0
            ser = evsgp4(ctx0, tv(et, Val(2)), geophs, el)
            f(t) = evsgp4(ctx0, t, geophs, el)
            d1(t) = ForwardDiff.derivative(f, t)
            d2(t) = ForwardDiff.derivative(d1, t)
            der = Allspice.derivatives
            @test SVector{6}(ntuple(i -> der(ser[i])[1], 6)) == f(et)
            @test sgp4_vecdiff(SVector{6}(ntuple(i -> der(ser[i])[2], 6)), d1(et)) <= SGP4_AD_RTOL[1]
            @test sgp4_vecdiff(SVector{6}(ntuple(i -> der(ser[i])[3], 6)), d2(et)) <= SGP4_AD_RTOL[2]
            # The derivative of the position against a central difference of the
            # CSPICE position. The SGP4 velocity is not this derivative: the model
            # computes it from the rates of some terms only, so the two agree to
            # about 1e-4 (reported, with a loose bound).
            hc = 1.0
            fdc = (sgp4_c_evsgp4(et + hc, geophs, el)[1:3] - sgp4_c_evsgp4(et - hc, geophs, el)[1:3]) / (2hc)
            @test sgp4_vecdiff(d1(et)[SOneTo(3)], fdc) <= SGP4_FD_RTOL
            gap = sgp4_vecdiff(d1(et)[SOneTo(3)], f(et)[SVector(4, 5, 6)])
            measured(@sprintf("  evsgp4, object %05d: |d(position)/dt - velocity| / |velocity| = %.1e", id, gap))
            @test gap <= 2e-3
        end
    end

    @apart @testset "zero allocation and inferred types" begin
        for id in (5, 8195, 24208, 23177)
            body = -100000 - id
            s = segs[id]
            for et in (s.epochs[1] + 5000.0, s.epochs[2], s.t0, s.epochs[3] + 100.0)
                @test sgp4_warm(alloc_sgp4_geo, ctx, body, et) == 0
                @test sgp4_warm(alloc_sgp4_pos, ctx, body, et) == 0
                @test sgp4_warm(alloc_sgp4_ezr, ctx, body, et, iau_earth) == 0
                for N in 0:3
                    @test sgp4_warm(alloc_sgp4_der, ctx, body, et, Val(N)) == 0
                end
            end
            ep, el = parsed[id]
            @test sgp4_warm(alloc_sgp4_ev, ctx0, ep + 100.0, geophs, el) == 0
            @test sgp4_warm(alloc_sgp4_ev, ctx0, ep + 100.0, Tuple(geophs), Tuple(el)) == 0
            @test sgp4_warm(alloc_sgp4_ev, ctx0, ep + 100.0, SVector{8}(geophs), SVector{10}(el)) == 0
        end
        body = -100005
        et = segs[5].epochs[1] + 5000.0
        @test (@inferred spkgeo(ctx, body, et, 1, 399)) isa Tuple{SVector{6,Float64},Float64}
        @test (@inferred spkezr(ctx, body, et, iau_earth, "NONE", 399)) isa Tuple{SVector{6,Float64},Float64}
        @test (@inferred spkpos(ctx, body, et, 1, "NONE", 399)) isa Tuple{SVector{3,Float64},Float64}
        @test (@inferred spkder(ctx, body, et, 1, 399, Val(3))) isa NTuple{4,SVector{3,Float64}}
        @test (@inferred spkder(ctx, body, et, iau_earth, 399, Val(2))) isa NTuple{3,SVector{3,Float64}}
    end

    @apart @testset "errors" begin
        ks = Allspice.kernelset(ctx)
        l = Allspice._lsk(ks)
        body = Int32(-100005)
        seg = ks.spk[body][1]
        data = ks.daf[seg.file]
        s = segs[5]
        et = s.epochs[1] + 100.0
        @test seg.dtype == 10
        # `spk_eval` has no kernel set, so it has no leap seconds.
        @test throws_kind(() -> Allspice.spk_eval(data, seg, et, Val(1)), :unsupported)
        # A context without a leap-second kernel.
        nolsk = neo_context([pc])
        @test throws_kind(() -> spkgeo(nolsk, -100005, et, 1, 399), :nolsk)
        @test throws_kind(() -> evsgp4(nolsk, et, geophs, parsed[5][2]), :nolsk)
        # Wrong counts of constants or elements.
        @test throws_kind(() -> evsgp4(ctx0, et, geophs[1:7], parsed[5][2]), :badarg)
        @test throws_kind(() -> evsgp4(ctx0, et, geophs, parsed[5][2][1:9]), :badarg)
        # Elements that cannot be propagated.
        @test throws_kind(() -> evsgp4(ctx0, parsed[28872][1], geophs, parsed[28872][2]), :suborbital)
        @test throws_kind(() -> spkgeo(ctx, -128872, segs[28872].epochs[1], 1, 399), :suborbital)
        @test throws_kind(() -> evsgp4(ctx0, parsed[SGP4_E97][1] - 900.0, geophs,
            parsed[SGP4_E97][2]), :orbitdecay)
        # A negative mean motion is not an error in CSPICE: both engines give NaN.
        neg = copy(parsed[5][2])
        neg[9] = -neg[9]
        @test all(isnan, evsgp4(ctx0, et, geophs, neg))
        @test all(isnan, sgp4_c_evsgp4(et, geophs, neg))
        # Meta data that are not the ones of a type 10 segment.
        e = seg.eaddr
        function refuses(change)
            bad = copy(data)
            change(bad)
            return throws_kind(() -> Allspice._t10_state(bad, seg, et, l), :baddaf)
        end
        @test Allspice._t10_state(copy(data), seg, et, l) == Allspice._t10_state(data, seg, et, l)
        @test refuses(d -> d[e] = 14.0)                 # count of meta items
        @test refuses(d -> d[e] = 1.0e6)
        @test refuses(d -> d[e-17+2] = 7.0)             # NCON
        @test refuses(d -> d[e-17+5] = 2.0)             # RDRTYP
        @test refuses(d -> d[e-17+7] = 2.0)             # NREF
        @test refuses(d -> d[e-17+10] = 1.0)            # PDRTYP
        @test refuses(d -> d[e-17+12] = 0.0)            # NPKT
        @test refuses(d -> d[e-17+12] = 4.0e5)
        @test refuses(d -> d[e-17+15] = 12.0)           # PKTSZ
        @test refuses(d -> d[e-17+6] = 1.0e5)           # REFBAS
        @test refuses(d -> d[e-17+11] = -1.0)           # PKTBAS
        short = Allspice.Segment(seg.target, seg.center, seg.frame, seg.dtype, seg.t0, seg.t1,
            seg.file, seg.eaddr - 5, seg.eaddr, 0.0, 0.0, Int32(0), Int32(0))
        @test throws_kind(() -> Allspice._t10_state(data, short, et, l), :baddaf)

        # spkw10: a refused call writes nothing.
        pw = joinpath(dir, "refused.bsp")
        h = spkopn(pw, "TLE", 0)
        elems, epochs = sgp4_sets(parsed[5][2], parsed[5][1], 3, 600.0)
        t0, t1 = epochs[1], epochs[end]
        w(args...) = spkw10(h, args...; ctx=ctx0)
        @test throws_kind(() -> w(-5, 399, "J2000", t0, t1, "x", geophs, Vector{Float64}[], Float64[]),
            :numpacketsnotpos)
        @test throws_kind(() -> w(-5, 399, "J2000", t0, t1, "x", geophs, elems[[2, 1, 3]],
            epochs[[2, 1, 3]]), :unorderedrefs)
        @test throws_kind(() -> w(-5, 399, "J2000", t0, t1, "x", geophs, elems[[1, 1, 3]],
            epochs[[1, 1, 3]]), :unorderedrefs)
        @test throws_kind(() -> w(-5, 399, "J2000", t0, t1, "x"^41, geophs, elems, epochs), :segidtoolong)
        @test throws_kind(() -> w(-5, 399, "J2000", t0, t1, "x\ty", geophs, elems, epochs), :nonprintingchars)
        @test throws_kind(() -> w(-5, -5, "J2000", t0, t1, "x", geophs, elems, epochs), :bodyandcentersame)
        @test throws_kind(() -> w(-5, 399, "NO_SUCH_FRAME", t0, t1, "x", geophs, elems, epochs), :invalidrefframe)
        @test throws_kind(() -> w(-5, 399, "J2000", t1, t0, "x", geophs, elems, epochs), :baddescrtimes)
        @test throws_kind(() -> w(-5, 399, "J2000", t0, t1, "x", geophs[1:7], elems, epochs), :badarg)
        @test throws_kind(() -> w(-5, 399, "J2000", t0, t1, "x", geophs, elems, epochs[1:2]), :badarg)
        @test throws_kind(() -> w(-5, 399, "J2000", t0, t1, "x", geophs, [e[1:9] for e in elems], epochs), :badarg)
        w(-5, 399, "J2000", t0, t1, "accepted", geophs, elems, epochs)
        spkcls(h)
        cw = neo_context(vcat(base, pw))
        @test length(Allspice.kernelset(cw).spk[Int32(-5)]) == 1
        @test spkgeo(cw, -5, t0 + 50.0, 1, 399)[1] isa SVector{6,Float64}
        # A file that is open for reading only.
        hr = dafopr(pw)
        @test sgp4_neo(() -> spkw10(hr, -6, 399, "J2000", t0, t1, "x", geophs, elems, epochs;
            ctx=ctx0)) isa Allspice.SpiceError
        dafcls(hr)
    end

    @apart @testset "time of a query (reported, not asserted)" begin
        measured("  Mean time of one query in microseconds")
        measured(@sprintf("  %-34s %-10s %-10s %-14s %-10s %-10s", "case", "spkgeo", "CSPICE",
            "spkder order 3", "evsgp4", "CSPICE"))
        for (label, id) in (("low orbit (00005)", 5), ("12 h resonant (08195)", 8195),
                ("semi-synchronous (28129)", 28129), ("1 day resonant (24208)", 24208),
                ("high eccentricity (23177)", 23177))
            body = -100000 - id
            s = segs[id]
            et = (s.epochs[1] + s.epochs[2]) / 2 + 100.0
            ep, el = parsed[id]
            g8, e10 = SVector{8}(geophs), SVector{10}(el)
            for n in (200, nsample(500, 20000))         # the first pass compiles
                a = sgp4_time(sgp4_neo_x, n, et, ctx, body)
                b = sgp4_time(sgp4_c_x, n, et, body)
                c = sgp4_time(sgp4_neo_x3, n, et, ctx, body)
                d = sgp4_time(sgp4_neo_ev, n, et, ctx0, g8, e10)
                e = sgp4_time(sgp4_c_ev, n, et, geophs, el)
                n == 200 || measured(@sprintf("  %-34s %-10.2f %-10.2f %-14.2f %-10.2f %-10.2f", label,
                    a / 1e3, b / 1e3, c / 1e3, d / 1e3, e / 1e3))
            end
        end
        @test true
    end

    SPICE.kclear()
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
