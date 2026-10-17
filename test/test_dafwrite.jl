# Group 9: DAF access, SPK, CK and binary PCK writers, kernel information.
#
# Feature file: features/daf_write.feature
#
# CSPICE is the oracle. Every writer gets the same random, valid segment with both
# engines, into two temporary files, and the files are compared:
#
#   (a) data area, double by double (bit equality), over the whole address range of
#       every array, read with `dafgda` of both engines from both files;
#   (b) descriptors and segment names;
#   (c) what CSPICE reads from the Allspice file (`spkgeo`, `ckgpav`, `ckgp`) against
#       what it reads from its own file;
#   (d) comments, written with `dafac` and read with `dafec`, in both directions;
#   (e) many segments in one file, past one summary record, and appending;
#   (f) every input check.
#
# C routines that SPICE.jl does not wrap (`spkw01_`, `spkw19_`, `spkw21_`, `ckw06_`) or
# wraps with a wrong count (`spkw02_c`, `spkw03_c`, `ckw02_c`, `ckw03_c`, `pckw02_c`) are
# called on `SPICE.libcspice` directly, as SPICE.jl does. CSPICE is called from the main
# task only.

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of Allspice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`ALLSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0
isdefined(@__MODULE__, Symbol("@apart")) || include("apart.jl")

using Random: randstring

const LIBCSPICE = SPICE.libcspice

# ---------------------------------------------------------------------------------------
# Direct calls into CSPICE
# ---------------------------------------------------------------------------------------

"Columns of the matrix `m` as a vector of vectors (the form that SPICE.jl takes)."
cols(m::AbstractMatrix) = [m[:, j] for j in 1:size(m, 2)]

function c_spkw0203(sym::Symbol, h, a)
    ninrec, n = size(a.cdata)
    per = sym === :spkw02_c ? 3 : 6
    if sym === :spkw02_c
        ccall((:spkw02_c, LIBCSPICE), Cvoid, (Cint, Cint, Cint, Cstring, Cdouble, Cdouble,
            Cstring, Cdouble, Cint, Cint, Ptr{Cdouble}, Cdouble), h, a.body, a.center,
            a.frame, a.first, a.last, a.segid, a.intlen, n, ninrec ÷ per - 1, a.cdata, a.btime)
    else
        ccall((:spkw03_c, LIBCSPICE), Cvoid, (Cint, Cint, Cint, Cstring, Cdouble, Cdouble,
            Cstring, Cdouble, Cint, Cint, Ptr{Cdouble}, Cdouble), h, a.body, a.center,
            a.frame, a.first, a.last, a.segid, a.intlen, n, ninrec ÷ per - 1, a.cdata, a.btime)
    end
    SPICE.handleerror()
end

function c_spkw01(h, a)
    ccall((:spkw01_, LIBCSPICE), Cint, (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ptr{UInt8},
        Ref{Cdouble}, Ref{Cdouble}, Ptr{UInt8}, Ref{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint),
        h, a.body, a.center, a.frame, a.first, a.last, a.segid, length(a.epochs), a.dlines,
        a.epochs, sizeof(a.frame), sizeof(a.segid))
    SPICE.handleerror()
end

"Type 14 through the three calls of CSPICE; the records go in two additions when `a.cut > 0`."
function c_spk14(h, a)
    ccall((:spk14b_c, LIBCSPICE), Cvoid, (Cint, Cstring, Cint, Cint, Cstring, Cdouble, Cdouble,
        Cint), h, a.segid, a.body, a.center, a.frame, a.first, a.last, a.degree)
    SPICE.handleerror()
    n = length(a.epochs)
    for r in (a.cut > 0 ? (1:a.cut, a.cut+1:n) : (1:n,))
        ccall((:spk14a_c, LIBCSPICE), Cvoid, (Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), h,
            length(r), a.coeffs[:, r], a.epochs[r])
        SPICE.handleerror()
    end
    ccall((:spk14e_c, LIBCSPICE), Cvoid, (Cint,), h)
    SPICE.handleerror()
end

"The same segment with the writer object of the package."
function n_spk14(h, a)
    w = spk14b(h, a.segid, a.body, a.center, a.frame, a.first, a.last, a.degree)
    n = length(a.epochs)
    if a.cut > 0
        spk14a(w, a.cut, a.coeffs[:, 1:a.cut], a.epochs[1:a.cut])
        spk14a(w, vec(a.coeffs[:, a.cut+1:n]), a.epochs[a.cut+1:n])
    else
        spk14a(w, n, a.coeffs, a.epochs)
    end
    spk14e(w)
end

function c_spkw21(h, a)
    ccall((:spkw21_, LIBCSPICE), Cint, (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ptr{UInt8},
        Ref{Cdouble}, Ref{Cdouble}, Ptr{UInt8}, Ref{Cint}, Ref{Cint}, Ptr{Cdouble},
        Ptr{Cdouble}, Cint, Cint), h, a.body, a.center, a.frame, a.first, a.last, a.segid,
        length(a.epochs), size(a.dlines, 1), a.dlines, a.epochs, sizeof(a.frame), sizeof(a.segid))
    SPICE.handleerror()
end

function c_spkw19(h, a)
    ccall((:spkw19_, LIBCSPICE), Cint, (Ref{Cint}, Ref{Cint}, Ref{Cint}, Ptr{UInt8},
        Ref{Cdouble}, Ref{Cdouble}, Ptr{UInt8}, Ref{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}, Cint, Cint), h, a.body, a.center,
        a.frame, a.first, a.last, a.segid, length(a.npkts), Cint.(a.npkts), Cint.(a.subtps),
        Cint.(a.degres), a.packts, a.epochs, a.ivlbds, Cint(a.sellst), sizeof(a.frame),
        sizeof(a.segid))
    SPICE.handleerror()
end

function c_ckw06(h, a)
    ccall((:ckw06_, LIBCSPICE), Cint, (Ref{Cint}, Ref{Cint}, Ptr{UInt8}, Ref{Cint},
        Ref{Cdouble}, Ref{Cdouble}, Ptr{UInt8}, Ref{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cint}, Cint, Cint), h,
        a.inst, a.ref, Cint(a.avflag), a.first, a.last, a.segid, length(a.npkts),
        Cint.(a.npkts), Cint.(a.subtps), Cint.(a.degres), a.packts, a.rates, a.sclkdp,
        a.ivlbds, Cint(a.sellst), sizeof(a.ref), sizeof(a.segid))
    SPICE.handleerror()
end

function c_ckw01(h, a)
    avflag = a.avvs !== nothing
    ccall((:ckw01_c, LIBCSPICE), Cvoid, (Cint, Cdouble, Cdouble, Cint, Cstring, Cint, Cstring,
        Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), h, a.begtim, a.endtim, a.inst, a.ref,
        avflag, a.segid, length(a.sclkdp), a.sclkdp, a.quats, avflag ? a.avvs : zeros(3, 1))
    SPICE.handleerror()
end

function c_ckw02(h, a)
    ccall((:ckw02_c, LIBCSPICE), Cvoid, (Cint, Cdouble, Cdouble, Cint, Cstring, Cstring, Cint,
        Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), h, a.begtim,
        a.endtim, a.inst, a.ref, a.segid, length(a.start), a.start, a.stop, a.quats, a.avvs,
        a.rates)
    SPICE.handleerror()
end

function c_ckw03(h, a)
    avflag = a.avvs !== nothing
    ccall((:ckw03_c, LIBCSPICE), Cvoid, (Cint, Cdouble, Cdouble, Cint, Cstring, Cint, Cstring,
        Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}), h, a.begtim,
        a.endtim, a.inst, a.ref, avflag, a.segid, length(a.sclkdp), a.sclkdp, a.quats,
        avflag ? a.avvs : zeros(3, 1), length(a.starts), a.starts)
    SPICE.handleerror()
end

"""
    c_pckw02(h, a)

`pckw02_c` with the state of a fresh CSPICE process (SPICE.jl 0.2.5 passes the degree
`length(cdata[1]) - 1`, three times too large, so the wrapper is not used).

A binary PCK summary has 5 integers, so the last 4 bytes of its third double are
padding. CSPICE packs the summary in a static buffer of `dafps_`/`dafus_` and leaves
there the sixth integer of the last summary of 6 or more integers that it packed or
unpacked, often the end address of an SPK or CK segment of an earlier call; Allspice
writes zero, the value of a fresh process. A `dafps_c` with 6 zero integers sets the
buffer to that state.
"""
function c_pckw02(h, a)
    ccall((:dafps_c, LIBCSPICE), Cvoid, (Cint, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}),
        2, 6, zeros(2), zeros(Cint, 6), zeros(5))
    ninrec, n = size(a.cdata)
    ccall((:pckw02_c, LIBCSPICE), Cvoid, (Cint, Cint, Cstring, Cdouble, Cdouble, Cstring,
        Cdouble, Cint, Cint, Ptr{Cdouble}, Cdouble), h, a.clssid, a.frame, a.first, a.last,
        a.segid, a.intlen, n, ninrec ÷ 3 - 1, a.cdata, a.btime)
    SPICE.handleerror()
end

function c_ckw05(h, a)
    ccall((:ckw05_c, LIBCSPICE), Cvoid, (Cint, Cint, Cint, Cdouble, Cdouble, Cint, Cstring,
        Cint, Cstring, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cint, Ptr{Cdouble}), h,
        a.subtyp, a.degree, a.begtim, a.endtim, a.inst, a.ref, a.avflag, a.segid,
        length(a.sclkdp), a.sclkdp, a.packts, a.rate, a.nints, a.starts)
    SPICE.handleerror()
end

# ---------------------------------------------------------------------------------------
# File comparison
# ---------------------------------------------------------------------------------------

"Everything in the DAF `path` as Allspice reads it: summaries, names, arrays, comments, file record."
function neo_dump(path::AbstractString)
    h = dafopr(path)
    try
        nd, ni, ifname, fward, bward, free = dafrfr(h)
        sums = Tuple{Vector{Float64},Vector{Int}}[]
        names = String[]
        data = Vector{Float64}[]
        dafbfs(h)
        while daffna(h)
            dc, ic = dafus(dafgs(h), nd, ni)
            push!(sums, (dc, ic))
            push!(names, dafgn(h))
            push!(data, ic[ni] >= ic[ni-1] ? dafgda(h, ic[ni-1], ic[ni]) : Float64[])
        end
        return (; sums, names, data, comments=dafec(h), record=(nd, ni, ifname, fward, bward, free))
    finally
        dafcls(h)
    end
end

"Everything in the DAF `path` as CSPICE reads it (same layout as `neo_dump`)."
function c_dump(path::AbstractString)
    h = SPICE.dafopr(path)
    try
        nd, ni, ifname, fward, bward, free = SPICE.dafrfr(h)
        sums = Tuple{Vector{Float64},Vector{Int}}[]
        names = String[]
        data = Vector{Float64}[]
        SPICE.dafbfs(h)
        while SPICE.daffna()
            dc, ic = SPICE.dafus(SPICE.dafgs(nd + (ni + 1) ÷ 2), nd, ni)
            push!(sums, (dc, ic))
            push!(names, SPICE.dafgn())
            push!(data, ic[ni] >= ic[ni-1] ? SPICE.dafgda(h, ic[ni-1], ic[ni]) : Float64[])
        end
        return (; sums, names, data, comments=String.(SPICE.dafec(h)),
            record=(nd, ni, ifname, fward, bward, free))
    finally
        SPICE.dafcls(h)
    end
end

"Number of doubles of `a` and `b` that differ as bit patterns (a length difference counts)."
function bit_mismatches(a::Vector{Vector{Float64}}, b::Vector{Vector{Float64}})
    bad = abs(sum(length, a; init=0) - sum(length, b; init=0))
    for (x, y) in zip(a, b)
        n = min(length(x), length(y))
        bad += count(i -> reinterpret(UInt64, x[i]) != reinterpret(UInt64, y[i]), 1:n)
    end
    return bad
end

"Bit equality of the descriptors (doubles as bit patterns, integers as values)."
same_sums(a, b) = length(a) == length(b) && all(zip(a, b)) do (x, y)
    reinterpret(UInt64, x[1]) == reinterpret(UInt64, y[1]) && x[2] == y[2]
end

"""
    compare_files(pc, pn) -> NamedTuple

Compare the file `pc` that CSPICE wrote with the file `pn` that Allspice wrote. The
reference is CSPICE reading its own file. `words` is the number of data doubles,
`badwords` the number that differ in any of the three other readings (CSPICE reads the
Allspice file; Allspice reads each file). `bytes` is true if the two files are equal
byte for byte.
"""
function compare_files(pc::AbstractString, pn::AbstractString)
    ref = c_dump(pc)
    others = (c_dump(pn), neo_dump(pc), neo_dump(pn))
    return (nseg=length(ref.sums), words=sum(length, ref.data; init=0),
        badwords=sum(o -> bit_mismatches(ref.data, o.data), others),
        sums=all(o -> same_sums(ref.sums, o.sums), others),
        names=all(o -> o.names == ref.names, others),
        comments=all(o -> o.comments == ref.comments, others),
        record=all(o -> o.record == ref.record, others),
        bytes=read(pc) == read(pn))
end

"True if all structural parts of a `compare_files` result agree."
files_agree(r) = r.badwords == 0 && r.sums && r.names && r.comments && r.record

# ---------------------------------------------------------------------------------------
# Random valid input
# ---------------------------------------------------------------------------------------

"`n` strictly increasing epochs from `t0` with random steps in [`lo`, `hi`]."
function inc_epochs(rng::AbstractRNG, n::Integer; t0::Float64=1000.0 * rand(rng), lo=1.0, hi=90.0)
    out = Vector{Float64}(undef, n)
    t = t0
    for i in 1:n
        out[i] = t
        t += lo + (hi - lo) * rand(rng)
    end
    return out
end

"Random segment identifier of printable characters, at most 40."
rand_segid(rng::AbstractRNG, tag::AbstractString) =
    tag * " " * randstring(rng, ['A':'Z'; '0':'9'; ' '; '_'; '-'], rand(rng, 0:(38 - length(tag)))) * "X"

"`n` random states (km, km/s) as a 6 x n matrix."
rand_states(rng::AbstractRNG, n::Integer) = randn(rng, 6, n) .* [7.0e3, 7.0e3, 7.0e3, 5.0, 5.0, 5.0]

"`n` unit quaternions as a 4 x n matrix; each is close to its predecessor (same hemisphere)."
function rand_quats(rng::AbstractRNG, n::Integer)
    q = Matrix{Float64}(undef, 4, n)
    cur = normalize(randn(rng, 4))
    for j in 1:n
        q[:, j] = cur
        cur = normalize(cur + 0.05 * randn(rng, 4))
    end
    return q
end

"Descriptor bounds for the epochs `ep`: the ends, or points inside the first and last step."
function rand_bounds(rng::AbstractRNG, ep::Vector{Float64})
    (length(ep) < 2 || rand(rng) < 0.5) && return ep[1], ep[end]
    return ep[1] + 0.25 * (ep[2] - ep[1]), ep[end] - 0.25 * (ep[end] - ep[end-1])
end

const SIZES = (2, 3, 7, 40, 99, 100, 101, 199, 200, 201, 257)
rand_n(rng::AbstractRNG, lo::Integer) = max(lo, rand(rng, SIZES))

"One difference line of `dlsize = 4*maxdim + 11` doubles with plausible content (final epoch `tl`)."
function rand_dline(rng::AbstractRNG, tl::Float64, maxdim::Integer)
    line = zeros(4 * maxdim + 11)
    line[1] = tl
    line[2:maxdim+1] .= 10.0 .+ rand(rng, maxdim)                   # step sizes G
    line[maxdim+2:maxdim+7] = [7000.0, 1.0, -3000.0, 6.0, 2000.0, 3.0] .* (1 .+ 0.1 .* randn(rng, 6))
    line[maxdim+8:4*maxdim+7] .= 1.0e-6 .* randn(rng, 3 * maxdim)   # differences DT
    line[4*maxdim+8] = 4.0                                          # KQMAX1
    line[4*maxdim+9:4*maxdim+11] .= 3.0                             # KQ
    return line
end

# ---------------------------------------------------------------------------------------
# Writer cases
# ---------------------------------------------------------------------------------------

"""
    SegArgs(nt)

Arguments of one segment for a writer: the named tuple `nt` behind a field without a
concrete type. Each writer has its own tuple type; with this wrapper the helpers that
take the arguments of any writer (`write_both`, `c_read_spk`, the `bad_*` lists) are
compiled once and not once for each writer. Fields are read as `a.name`, and
`merge(a, (name=value,))` gives the arguments with a field replaced.
"""
struct SegArgs
    nt::NamedTuple
end
Base.getproperty(a::SegArgs, s::Symbol) = getproperty(getfield(a, :nt), s)
Base.hasproperty(a::SegArgs, s::Symbol) = hasproperty(getfield(a, :nt), s)
Base.propertynames(a::SegArgs) = propertynames(getfield(a, :nt))
Base.values(a::SegArgs) = values(getfield(a, :nt))
Base.merge(a::SegArgs, b::NamedTuple) = SegArgs(merge(getfield(a, :nt), b))

"""
One writer under test. `gen(rng, k)` gives valid arguments for segment number `k` of a
file as a named tuple; `neo` and `cspice` write them; `bad(a)` lists invalid variants
of the arguments `a` as `(label, arguments, kind, both)`, where `both` says that CSPICE
rejects the variant too (false for a check that only the Julia interface has).
"""
struct WriterCase
    name::String
    kind::Symbol
    gen::Function
    neo::Function
    cspice::Function
    bad::Function
    # `gen` gives a named tuple; the case hands it out as `SegArgs`.
    WriterCase(name, kind, gen, neo, cspice, bad) =
        new(name, kind, (rng, k) -> SegArgs(gen(rng, k)), neo, cspice, bad)
end

"The frame and identifier checks that all SPK writers share."
function bad_ids(a; long::Symbol=:segidtoolong, fname::Symbol=:frame, empty_both::Bool=true)
    out = Any[
        ("unknown frame", merge(a, NamedTuple{(fname,)}(("NO_SUCH_FRAME",))), :invalidrefframe, true),
        ("identifier of 41 characters", merge(a, (segid="A"^41,)), long, true),
        ("nonprintable identifier", merge(a, (segid="BAD\tID",)), :nonprintablechars, true),
        ("empty frame", merge(a, NamedTuple{(fname,)}(("",))), :emptystring, empty_both),
        ("empty identifier", merge(a, (segid="",)), :emptystring, empty_both),
    ]
    return out
end

"The checks of `spkpds`."
bad_pds(a) = Any[
    ("body is the barycenter", merge(a, (body=0,)), :barycenterephem, true),
    ("body equal to center", merge(a, (center=a.body,)), :bodyandcentersame, true),
]

spk_ids(k) = (body=-1000 - k, center=399, frame="J2000")

# --- discrete states ---------------------------------------------------------------------

function gen_w13(rng, k)
    n = rand_n(rng, 2)
    degree = rand(rng, [d for d in 1:2:27 if (d + 1) ÷ 2 <= n])
    epochs = inc_epochs(rng, n)
    first, last = rand_bounds(rng, epochs)
    return (; spk_ids(k)..., first, last, segid=rand_segid(rng, "T13"), degree, states=rand_states(rng, n), epochs)
end
function gen_w14(rng, k)
    n = rand_n(rng, 1)
    degree = rand(rng, 0:12)
    bounds = inc_epochs(rng, n + 1)
    coeffs = randn(rng, 2 + 6 * (degree + 1), n)
    for i in 1:n
        coeffs[1, i] = (bounds[i] + bounds[i+1]) / 2
        coeffs[2, i] = (bounds[i+1] - bounds[i]) / 2
    end
    return (; spk_ids(k)..., first=bounds[1], last=bounds[end], segid=rand_segid(rng, "T14"), degree,
        coeffs, epochs=bounds[1:n], cut=n > 1 ? rand(rng, 0:n-1) : 0)
end
function bad_w14(a)
    out = Any[bad_ids(a)..., bad_pds(a)...,
        ("negative degree", merge(a, (degree=-1,)), :invalidargument, true)]
    if length(a.epochs) >= 3
        dup = copy(a.epochs)
        dup[3] = dup[2]
        # Not sent to CSPICE (SPICE(UNORDEREDREFS) in `sgwfpk_`): a rejected `spk14a_c`
        # leaves the segment that `spk14b_c` began in the file, and it cannot be removed.
        push!(out, ("repeated epoch", merge(a, (epochs=dup, cut=0)), :unorderedrefs, false))
        push!(out, ("second addition not after the first", merge(a, (epochs=dup, cut=2)), :unorderedrefs, false))
        push!(out, ("epochs shorter than records", merge(a, (epochs=a.epochs[1:end-1], cut=0)), :badarg, false))
    end
    return out
end
function bad_w13(a)
    dup = copy(a.epochs)
    dup[3] = dup[2]
    return Any[bad_ids(a)..., bad_pds(a)...,
        ("degree 0", merge(a, (degree=0,)), :invaliddegree, true),
        ("degree 29", merge(a, (degree=29,)), :invaliddegree, true),
        ("even degree", merge(a, (degree=4,)), :invaliddegree, true),
        ("too few states", merge(a, (degree=5, states=a.states[:, 1:2], epochs=a.epochs[1:2],
            last=a.epochs[2])), :toofewstates, true),
        ("first equal to last", merge(a, (last=a.first,)), :baddescrtimes, true),
        ("repeated epoch", merge(a, (epochs=dup,)), :timesoutoforder, true),
        ("first before the first epoch", merge(a, (first=a.epochs[1] - 1.0,)), :baddescrtimes, true),
        ("last after the last epoch", merge(a, (last=a.epochs[end] + 1.0,)), :baddescrtimes, true),
        ("epochs shorter than states", merge(a, (epochs=a.epochs[1:end-1],)), :badarg, false),
        ("states with 5 rows", merge(a, (states=a.states[1:5, :],)), :badarg, false),
    ]
end

function gen_w09(rng, k)
    n = rand_n(rng, 2)
    degree = rand(rng, 1:min(27, n - 1))
    epochs = inc_epochs(rng, n)
    first, last = rand_bounds(rng, epochs)
    return (; spk_ids(k)..., first, last, segid=rand_segid(rng, "T09"), degree, states=rand_states(rng, n), epochs)
end
function bad_w09(a)
    dup = copy(a.epochs)
    dup[3] = dup[2]
    return Any[bad_ids(a)...,
        ("degree 0", merge(a, (degree=0,)), :invaliddegree, true),
        ("degree 28", merge(a, (degree=28,)), :invaliddegree, true),
        ("too few states", merge(a, (degree=3, states=a.states[:, 1:3], epochs=a.epochs[1:3],
            last=a.epochs[3])), :toofewstates, true),
        ("first equal to last", merge(a, (last=a.first,)), :baddescrtimes, true),
        ("repeated epoch", merge(a, (epochs=dup,)), :timesoutoforder, true),
        ("first before the first epoch", merge(a, (first=a.epochs[1] - 1.0,)), :baddescrtimes, true),
        ("last after the last epoch", merge(a, (last=a.epochs[end] + 1.0,)), :baddescrtimes, true),
    ]
end

function gen_w05(rng, k)
    n = rand_n(rng, 1)
    epochs = inc_epochs(rng, n)
    first, last = n == 1 ? (epochs[1] - 10.0, epochs[1] + 10.0) : rand_bounds(rng, epochs)
    return (; spk_ids(k)..., first, last, segid=rand_segid(rng, "T05"), gm=398600.4418 * (1 + rand(rng)),
        states=rand_states(rng, n), epochs)
end
function bad_w05(a)
    dup = copy(a.epochs)
    dup[3] = dup[2]
    return Any[bad_ids(a)...,
        ("zero GM", merge(a, (gm=0.0,)), :nonpositivemass, true),
        ("repeated epoch", merge(a, (epochs=dup,)), :unorderedtimes, true),
        ("first after last", merge(a, (first=a.last + 1.0,)), :baddescrtimes, true),
        ("no state", merge(a, (states=zeros(6, 0), epochs=Float64[])), :numstatesnotpos, false),
    ]
end

function gen_w0812(rng, k, hermite::Bool)
    n = rand_n(rng, 2)
    degree = hermite ? rand(rng, [d for d in 1:2:27 if (d + 1) ÷ 2 <= n]) : rand(rng, 1:min(27, n - 1))
    epoch1 = 1000.0 * rand(rng)
    step = 10.0 + 50.0 * rand(rng)
    last = epoch1 + (n - 1) * step
    rand(rng) < 0.5 && (last -= 0.3 * step)
    return (; spk_ids(k)..., first=epoch1, last, segid=rand_segid(rng, hermite ? "T12" : "T08"),
        degree, states=rand_states(rng, n), epoch1, step)
end
function bad_w0812(a, hermite::Bool)
    out = Any[bad_ids(a)...,
        ("degree 0", merge(a, (degree=0,)), :invaliddegree, true),
        ("degree 28", merge(a, (degree=28,)), :invaliddegree, true),
        ("first equal to last", merge(a, (last=a.first,)), :baddescrtimes, true),
        ("zero step", merge(a, (step=0.0,)), :invalidstepsize, true),
        ("first before the data", merge(a, (first=a.epoch1 - 1.0,)), :coveragegap, true),
        ("last after the data", merge(a, (last=a.last + 2 * a.step,)), :coveragegap, true),
    ]
    if hermite
        push!(out, ("even degree", merge(a, (degree=4,)), :invaliddegree, true))
        push!(out, ("too few states", merge(a, (degree=7, states=a.states[:, 1:3],
            last=a.epoch1 + 2 * a.step)), :toofewstates, true))
        append!(out, bad_pds(a))
    else
        push!(out, ("too few states", merge(a, (degree=3, states=a.states[:, 1:3],
            last=a.epoch1 + 2 * a.step)), :toofewstates, true))
    end
    return out
end

function gen_w18(rng, k)
    subtyp = rand(rng, 0:1)
    degree = subtyp == 0 ? rand(rng, (3, 7, 11, 15)) : rand(rng, 1:2:15)
    n = rand_n(rng, 2)
    epochs = inc_epochs(rng, n)
    first, last = rand_bounds(rng, epochs)
    states = rand_states(rng, n)
    packts = subtyp == 0 ? vcat(states, states[4:6, :] .* (1 .+ 1e-3 .* randn(rng, 3, n)),
        1e-3 .* randn(rng, 3, n)) : states
    return (; subtyp, spk_ids(k)..., first, last, segid=rand_segid(rng, "T18"), degree, packts, epochs)
end
function bad_w18(a)
    dup = copy(a.epochs)
    dup[3] = dup[2]
    packsz = size(a.packts, 1)
    return Any[bad_ids(a)...,
        ("subtype 2", merge(a, (subtyp=2,)), :invalidvalue, false),
        ("degree 0", merge(a, (degree=0,)), :invaliddegree, true),
        ("degree 17", merge(a, (degree=17,)), :invaliddegree, true),
        ("odd window size", merge(a, (degree=a.subtyp == 0 ? 5 : 4,)), :invaliddegree, true),
        ("one packet", merge(a, (packts=a.packts[:, 1:1], epochs=a.epochs[1:1],
            first=a.epochs[1], last=a.epochs[1])), :toofewstates, true),
        ("first after last", merge(a, (first=a.last + 1.0,)), :baddescrtimes, true),
        ("repeated epoch", merge(a, (epochs=dup,)), :timesoutoforder, true),
        ("first before the first epoch", merge(a, (first=a.epochs[1] - 1.0,)), :baddescrtimes, true),
        ("last after the last epoch", merge(a, (last=a.epochs[end] + 1.0,)), :baddescrtimes, true),
        ("packets of the other size", merge(a, (packts=a.packts[1:packsz-1, :],)), :badarg, false),
    ]
end

# --- Chebyshev types ---------------------------------------------------------------------

function gen_w0203(rng, k, per::Int)
    n = rand(rng, (1, 2, 9, 40, 130))
    polydg = rand(rng, (0, 1, 5, 12, 27))
    intlen = 3600.0 * (1 + 23 * rand(rng))
    btime = 1.0e5 * randn(rng)
    cdata = randn(rng, per * (polydg + 1), n) .* 1.0e4
    first, last = rand(rng) < 0.5 ? (btime, btime + n * intlen) : (btime + 0.1 * intlen, btime + (n - 0.1) * intlen)
    return (; spk_ids(k)..., first, last, segid=rand_segid(rng, per == 3 ? "T02" : "T03"), intlen, cdata, btime)
end
function bad_w0203(a, per::Int)
    return Any[bad_ids(a; long=:idstringtoolong)...,
        ("no record", merge(a, (cdata=zeros(size(a.cdata, 1), 0),)), :numcoeffsnotpos, true),
        ("degree 28", merge(a, (cdata=zeros(per * 29, size(a.cdata, 2)),)), :invaliddegree, true),
        ("zero interval length", merge(a, (intlen=0.0,)), :intlennotpos, true),
        ("first after last", merge(a, (first=a.last + 1.0,)), :baddescrtimes, true),
        ("first before the data", merge(a, (first=a.btime - 10.0,)), :coveragegap, true),
        ("last after the data", merge(a, (last=a.btime + size(a.cdata, 2) * a.intlen + 10.0,)), :coveragegap, true),
        ("record size not a multiple", merge(a, (cdata=zeros(per + 1, size(a.cdata, 2)),)), :badarg, false),
    ]
end

function gen_w20(rng, k)
    n = rand(rng, (1, 2, 9, 40, 130))
    polydg = rand(rng, (0, 1, 5, 12, 27, 50))
    intlen = 1.0 + 15.0 * rand(rng)                 # days
    initjd = 2451545.0 + rand(rng, -3000:3000)
    initfr = rand(rng, (0.0, 0.25, 0.5, 0.123456789))
    btime = 86400.0 * (initjd - 2451545.0 + initfr)
    ltime = btime + n * intlen * 86400.0
    cdata = randn(rng, 3 * (polydg + 2), n)
    return (; spk_ids(k)..., first=btime, last=ltime, segid=rand_segid(rng, "T20"), intlen, cdata,
        dscale=rand(rng, (1.0, 149597870.7)), tscale=rand(rng, (1.0, 86400.0)), initjd, initfr)
end
function bad_w20(a)
    return Any[bad_ids(a; long=:idstringtoolong)...,
        ("degree 51", merge(a, (cdata=zeros(3 * 53, size(a.cdata, 2)),)), :invaliddegree, true),
        ("zero interval length", merge(a, (intlen=0.0,)), :intlennotpos, true),
        ("first after last", merge(a, (first=a.last + 1.0,)), :baddescrtimes, true),
        ("zero distance scale", merge(a, (dscale=0.0,)), :nonpositivescale, true),
        ("negative time scale", merge(a, (tscale=-1.0,)), :nonpositivescale, true),
        ("first before the data", merge(a, (first=a.first - 100.0,)), :coveragegap, true),
        ("last after the data", merge(a, (last=a.last + 100.0,)), :coveragegap, true),
        ("no record", merge(a, (cdata=zeros(size(a.cdata, 1), 0),)), :invalidcount, false),
    ]
end

# --- conic elements ----------------------------------------------------------------------

function gen_w15(rng, k)
    tp = normalize(randn(rng, 3))
    pa = normalize(cross(tp, randn(rng, 3)))
    scale = 10.0^rand(rng, -3:3)                    # the writer stores unit vectors
    return (; spk_ids(k)..., first=-1.0e5 * rand(rng), last=1.0e5 * rand(rng) + 1.0,
        segid=rand_segid(rng, "T15"), epoch=100.0 * randn(rng), tp=scale .* tp, pa=pa ./ scale,
        p=7000.0 + 3000.0 * rand(rng), ecc=0.5 * rand(rng), j2flg=Float64(rand(rng, 0:3)),
        pv=randn(rng, 3), gm=398600.4418, j2=1.08263e-3, radius=6378.137)
end
bad_w15(a) = Any[bad_ids(a)..., bad_pds(a)...,
    ("zero semi-latus rectum", merge(a, (p=0.0,)), :badlatusrectum, true),
    ("negative eccentricity", merge(a, (ecc=-0.1,)), :badeccentricity, true),
    ("zero GM", merge(a, (gm=0.0,)), :nonpositivemass, true),
    ("zero trajectory pole", merge(a, (tp=zeros(3),)), :badvector, true),
    ("zero periapsis vector", merge(a, (pa=zeros(3),)), :badvector, true),
    ("zero central pole", merge(a, (pv=zeros(3),)), :badvector, true),
    ("negative radius", merge(a, (radius=-1.0,)), :badradius, true),
    # `spkw15` compares the dot product of the unit vectors with 1e-5 (spkw15.c, line 463).
    # The generator scales `tp` and `pa` by reciprocal factors, so the variant is built from
    # the unit vectors: with the scaled vectors and the factor 1e-3 the dot product is 1e-5,
    # which is on the limit and both engines accept it.
    ("pole not orthogonal to periapsis", merge(a, (pa=normalize(a.tp) .+ 0.1 .* normalize(a.pa),)),
        :badinitstate, true),
    ("first equal to last", merge(a, (last=a.first,)), :baddescrtimes, true),
    ("pole with two components", merge(a, (tp=[1.0, 0.0],)), :badarg, false),
]

function gen_w17(rng, k)
    eqel = [8000.0 + 4000.0 * rand(rng), 0.2 * randn(rng), 0.2 * randn(rng), 2pi * rand(rng),
        0.1 * randn(rng), 0.1 * randn(rng), 1e-7 * randn(rng), 1e-3 * (1 + rand(rng)), 1e-7 * randn(rng)]
    return (; spk_ids(k)..., first=-1.0e5 * rand(rng), last=1.0e5 * rand(rng) + 1.0,
        segid=rand_segid(rng, "T17"), epoch=100.0 * randn(rng), eqel, rapol=2pi * rand(rng),
        decpol=pi / 2 * rand(rng))
end
bad_w17(a) = Any[bad_ids(a)..., bad_pds(a)...,
    ("zero semimajor axis", merge(a, (eqel=[0.0; a.eqel[2:end]],)), :badsemiaxis, true),
    ("eccentricity above 0.9", merge(a, (eqel=[a.eqel[1]; 0.8; 0.8; a.eqel[4:end]],)), :badeccentricity, true),
    ("first equal to last", merge(a, (last=a.first,)), :baddescrtimes, true),
    ("eight elements", merge(a, (eqel=a.eqel[1:8],)), :badarg, false),
]

# --- difference lines --------------------------------------------------------------------

function gen_w0121(rng, k, extended::Bool)
    n = rand_n(rng, 1)
    maxdim = extended ? rand(rng, (15, 18, 25)) : 15
    epochs = inc_epochs(rng, n; t0=1000.0 * rand(rng) + 50.0)
    dlines = reduce(hcat, [rand_dline(rng, t, maxdim) for t in epochs])
    last = rand(rng) < 0.5 ? epochs[end] : epochs[end] - 0.5
    return (; spk_ids(k)..., first=epochs[1] - 40.0, last, segid=rand_segid(rng, extended ? "T21" : "T01"),
        dlines, epochs)
end
function bad_w0121(a, extended::Bool)
    dup = copy(a.epochs)
    length(dup) >= 3 && (dup[3] = dup[2])
    out = Any[bad_ids(a; empty_both=false)..., bad_pds(a)...,
        ("last after the last epoch", merge(a, (last=a.epochs[end] + 1.0,)),
            extended ? :coveragegap : :baddescrtimes, true),
        ("no line", merge(a, (dlines=zeros(size(a.dlines, 1), 0), epochs=Float64[])), :invalidcount, false),
    ]
    length(dup) >= 3 && push!(out, ("repeated epoch", merge(a, (epochs=dup,)), :timesoutoforder, true))
    if extended
        zstep = copy(a.dlines)
        zstep[3, 1] = 0.0
        push!(out, ("zero step", merge(a, (dlines=zstep,)), :zerostep, true))
        push!(out, ("line of 115 doubles", merge(a, (dlines=zeros(115, size(a.dlines, 2)),)), :difflinetoolarge, true))
        push!(out, ("line of 67 doubles", merge(a, (dlines=zeros(67, size(a.dlines, 2)),)), :difflinetoosmall, true))
        push!(out, ("first after last", merge(a, (first=a.last + 1.0,)), :baddescrtimes, true))
    else
        push!(out, ("first equal to last", merge(a, (first=a.last,)), :baddescrtimes, true))
        push!(out, ("line of 70 doubles", merge(a, (dlines=a.dlines[1:70, :],)), :badarg, false))
    end
    return out
end

# --- mini-segments -----------------------------------------------------------------------

function gen_w19(rng, k; nintvl::Integer=rand(rng, 1:4))
    npkts = Int[]
    subtps = Int[]
    degres = Int[]
    packts = Float64[]
    epochs = Float64[]
    ivlbds = Float64[1000.0 * rand(rng)]
    for _ in 1:nintvl
        n = nintvl > 10 ? rand(rng, 2:4) : rand(rng, (2, 5, 100, 101, 230))
        sub = rand(rng, 0:2)
        deg = sub == 1 ? rand(rng, 1:2:27) : rand(rng, (3, 7, 11, 15, 19, 23, 27))
        ep = inc_epochs(rng, n; t0=ivlbds[end] - (rand(rng) < 0.5 ? 0.0 : 5.0))
        # The first epoch can be 5 s before the start of the interval and a step can be as
        # short as 1 s, so a short mini-segment (2 packets) can end before its interval
        # starts. `spkw19` requires strictly increasing interval bounds
        # (SPICE(BOUNDSOUTOFORDER)), thus such epochs are drawn again; 0.5 s is the largest
        # distance of the stop bound from the last epoch.
        while ep[end] - 0.5 <= ivlbds[end]
            ep = inc_epochs(rng, n; t0=ivlbds[end] - (rand(rng) < 0.5 ? 0.0 : 5.0))
        end
        push!(npkts, n)
        push!(subtps, sub)
        push!(degres, deg)
        append!(packts, randn(rng, (12, 6, 6)[sub+1] * n) .* 1.0e3)
        append!(epochs, ep)
        push!(ivlbds, rand(rng) < 0.5 ? ep[end] : ep[end] - 0.5)
    end
    @assert all(>(0.0), diff(ivlbds)) "the interval bounds of a type 19 segment increase strictly"
    return (; spk_ids(k)..., first=ivlbds[1], last=ivlbds[end], segid=rand_segid(rng, "T19"), npkts,
        subtps, degres, packts, epochs, ivlbds, sellst=rand(rng, Bool))
end
function bad_w19(a)
    swapped = copy(a.ivlbds)
    swapped[2] = swapped[1]
    dup = copy(a.epochs)
    dup[2] = dup[1]
    one = merge(a, (npkts=[1; a.npkts[2:end]],))
    return Any[bad_ids(a; empty_both=false)...,
        ("first after last", merge(a, (first=a.last + 1.0,)), :baddescrtimes, true),
        ("bounds out of order", merge(a, (ivlbds=swapped,)), :boundsoutoforder, true),
        ("first before the bounds", merge(a, (first=a.ivlbds[1] - 1.0,)), :coveragegap, true),
        ("last after the bounds", merge(a, (last=a.ivlbds[end] + 1.0,)), :coveragegap, true),
        ("one packet", one, :toofewpackets, true),
        ("subtype 3", merge(a, (subtps=[3; a.subtps[2:end]],)), :invalidsubtype, true),
        ("degree 29", merge(a, (subtps=[1; a.subtps[2:end]], degres=[29; a.degres[2:end]])), :invaliddegree, true),
        ("odd window size", merge(a, (subtps=[1; a.subtps[2:end]], degres=[2; a.degres[2:end]])), :badwindowsize, true),
        ("repeated epoch", merge(a, (epochs=dup,)), :timesoutoforder, true),
        ("no interval", merge(a, (npkts=Int[], subtps=Int[], degres=Int[])), :invalidcount, false),
    ]
end

# --- CK ----------------------------------------------------------------------------------

ck_ids(k) = (inst=-77000 - k, ref="J2000")
"Encoded clock times: whole ticks, strictly increasing, from a non-negative start."
ck_times(rng, n) = cumsum(Float64[rand(rng, 0:5000); rand(rng, 1:4000, n - 1)])

function gen_ckw01(rng, k)
    nrec = rand_n(rng, 1)
    sclkdp = ck_times(rng, nrec)
    begtim, endtim = rand(rng) < 0.5 ? (sclkdp[1], sclkdp[end]) : (max(sclkdp[1] - 3, 0.0), sclkdp[end] + 7)
    return (; begtim, endtim, ck_ids(k)..., segid=rand_segid(rng, "C01"), sclkdp, quats=rand_quats(rng, nrec),
        avvs=rand(rng, Bool) ? 1e-3 .* randn(rng, 3, nrec) : nothing)
end
function bad_ck_common(a, tname::Symbol)
    t = getproperty(a, tname)
    dup = copy(t)
    length(dup) >= 3 && (dup[3] = dup[2])
    zq = copy(a.quats)
    zq[:, end] .= 0.0
    out = Any[bad_ids(a; fname=:ref)...,
        ("begin after the first time", merge(a, (begtim=t[1] + 0.5,)), :invaliddescrtime, true),
        ("zero quaternion", merge(a, (quats=zq,)), :zeroquaternion, true),
        ("quaternions with 3 rows", merge(a, (quats=a.quats[1:3, :],)), :badarg, false),
    ]
    length(dup) >= 3 && push!(out, ("repeated time", merge(a, NamedTuple{(tname,)}((dup,))), :timesoutoforder, true))
    return out
end
function bad_ckw01(a)
    neg = copy(a.sclkdp)
    neg[1] = -1.0
    return Any[bad_ck_common(a, :sclkdp)...,
        ("end before the last time", merge(a, (endtim=a.sclkdp[end] - 0.5,)), :invaliddescrtime, true),
        ("negative first time", merge(a, (sclkdp=neg, begtim=-2.0)), :invalidsclktime, true),
        ("no record", merge(a, (sclkdp=Float64[], quats=zeros(4, 0), avvs=nothing)), :invalidnumrec, true),
        ("two angular velocities for more records", merge(a, (sclkdp=a.sclkdp[1:min(end, 3)],
            quats=a.quats[:, 1:min(end, 3)], endtim=a.endtim, avvs=zeros(3, 2))),
            length(a.sclkdp) == 2 ? :none : :badarg, false),
    ]
end

function gen_ckw02(rng, k)
    nrec = rand_n(rng, 1)
    edges = ck_times(rng, 2nrec)
    start = edges[1:2:end]
    stop = edges[2:2:end]
    for i in 1:nrec-1                               # some intervals touch the next one
        rand(rng) < 0.3 && (stop[i] = start[i+1])
    end
    return (; begtim=start[1], endtim=stop[end], ck_ids(k)..., segid=rand_segid(rng, "C02"), start, stop,
        quats=rand_quats(rng, nrec), avvs=1e-3 .* randn(rng, 3, nrec), rates=fill(1 / 256, nrec) .* (1 .+ 1e-6 .* randn(rng, nrec)))
end
function bad_ckw02(a)
    neg = copy(a.start)
    neg[1] = -1.0
    deg = copy(a.stop)
    deg[1] = a.start[1]
    out = Any[bad_ck_common(a, :start)[[1, 2, 3, 4, 5, 6, 7, 8]]...,
        ("end before the last stop", merge(a, (endtim=a.stop[end] - 0.5,)), :invaliddescrtime, true),
        ("negative first start", merge(a, (start=neg, begtim=-2.0)), :invalidsclktime, true),
        ("degenerate interval", merge(a, (stop=deg,)), :degenerateinterval, true),
        ("no record", merge(a, (start=Float64[], stop=Float64[], rates=Float64[], quats=zeros(4, 0),
            avvs=zeros(3, 0))), :invalidnumrec, true),
    ]
    if length(a.start) >= 3
        late = copy(a.stop)
        late[1] = a.start[2] + 0.5
        push!(out, ("stop after the next start", merge(a, (stop=late,)), :badstoptime, true))
        back = copy(a.start)
        back[3] = back[2]
        push!(out, ("repeated start", merge(a, (start=back,)), :timesoutoforder, true))
    end
    return out
end

function gen_ckw03(rng, k)
    nrec = rand_n(rng, 1)
    sclkdp = ck_times(rng, nrec)
    pick = sort!(unique!([1; rand(rng, 1:nrec, rand(rng, (0, 1, 3, 120, 230)))]))
    begtim, endtim = rand(rng) < 0.5 ? (sclkdp[1], sclkdp[end]) : (max(sclkdp[1] - 3, 0.0), sclkdp[end] + 7)
    return (; begtim, endtim, ck_ids(k)..., segid=rand_segid(rng, "C03"), sclkdp, quats=rand_quats(rng, nrec),
        starts=sclkdp[pick], avvs=rand(rng, Bool) ? 1e-3 .* randn(rng, 3, nrec) : nothing)
end
function bad_ckw03(a)
    neg = copy(a.sclkdp)
    neg[1] = -1.0
    out = Any[bad_ck_common(a, :sclkdp)...,
        ("end before the last time", merge(a, (endtim=a.sclkdp[end] - 0.5,)), :invaliddescrtime, true),
        ("negative first time", merge(a, (sclkdp=neg, begtim=-2.0, starts=[-1.0; a.starts[2:end]])), :invalidsclktime, true),
        ("no record", merge(a, (sclkdp=Float64[], quats=zeros(4, 0), avvs=nothing)), :invalidnumrec, true),
        ("no interval", merge(a, (starts=Float64[],)), :invalidnumint, true),
        ("first start is not the first time", merge(a, (starts=[a.sclkdp[1] + 0.5],)), :timesdontmatch, true),
    ]
    if length(a.sclkdp) >= 3
        push!(out, ("start that is not a time", merge(a, (starts=[a.sclkdp[1], a.sclkdp[2] + 0.5],)), :invalidstarttime, true))
        push!(out, ("repeated start", merge(a, (starts=[a.sclkdp[1], a.sclkdp[2], a.sclkdp[2]],)), :timesoutoforder, true))
    end
    return out
end

function gen_ckw05(rng, k)
    subtyp = rand(rng, 0:3)
    degree = isodd(subtyp) ? rand(rng, 1:2:23) : rand(rng, (3, 7, 11, 15, 19, 23))
    n = rand_n(rng, 1)
    sclkdp = ck_times(rng, n)
    q = rand_quats(rng, n)
    packts = subtyp == 0 ? vcat(q, 1e-4 .* randn(rng, 4, n)) :
             subtyp == 1 ? q :
             subtyp == 2 ? vcat(q, 1e-4 .* randn(rng, 4, n), 1e-3 .* randn(rng, 3, n), 1e-6 .* randn(rng, 3, n)) :
             vcat(q, 1e-3 .* randn(rng, 3, n))
    pick = sort!(unique!([1; rand(rng, 1:n, rand(rng, (0, 1, 3, 120)))]))
    return (; subtyp, degree, begtim=sclkdp[1], endtim=sclkdp[end], ck_ids(k)..., avflag=rand(rng, Bool),
        segid=rand_segid(rng, "C05"), sclkdp, packts, rate=1 / 256, nints=length(pick), starts=sclkdp[pick])
end
function bad_ckw05(a)
    neg = copy(a.sclkdp)
    neg[1] = -1.0
    zq = copy(a.packts)
    zq[1:4, end] .= 0.0
    out = Any[bad_ids(a; fname=:ref)...,
        ("no interval", merge(a, (nints=0,)), :invalidnumints, true),
        ("negative first time", merge(a, (sclkdp=neg, starts=[-1.0; a.starts[2:end]], begtim=-1.0)), :invalidsclktime, true),
        ("start that is not a time", merge(a, (starts=[a.sclkdp[1] + 0.5], nints=1)), :invalidstarttime, true),
        ("subtype 4", merge(a, (subtyp=4,)), :notsupported, true),
        ("zero quaternion", merge(a, (packts=zq,)), :zeroquaternion, true),
        ("degree 0", merge(a, (degree=0,)), :invaliddegree, true),
        ("degree 25", merge(a, (degree=25,)), :invaliddegree, true),
        ("odd window size", merge(a, (degree=isodd(a.subtyp) ? 2 : 5,)), :invaliddegree, true),
        ("end before begin", merge(a, (endtim=a.begtim - 1.0,)), :baddescrtimes, true),
        ("all times before begin", merge(a, (begtim=a.sclkdp[end] + 1.0, endtim=a.sclkdp[end] + 2.0)), :emptysegment, true),
        ("no time in the bounds", merge(a, (begtim=a.sclkdp[1] + 0.25, endtim=a.sclkdp[1] + 0.5)), :emptysegment, true),
        ("zero rate", merge(a, (rate=0.0,)), :invalidvalue, true),
        ("nints above the starts", merge(a, (nints=length(a.starts) + 1,)), :badarg, false),
        ("no packet", merge(a, (sclkdp=Float64[], packts=zeros(size(a.packts, 1), 0))), :toofewpackets, false),
    ]
    if length(a.sclkdp) >= 3
        dup = copy(a.sclkdp)
        dup[3] = dup[2]
        push!(out, ("repeated time", merge(a, (sclkdp=dup,)), :timesoutoforder, true))
        push!(out, ("repeated start", merge(a, (starts=a.sclkdp[[1, 2, 2]], nints=3)), :timesoutoforder, true))
    end
    return out
end

function gen_ckw06(rng, k; nmini::Integer=rand(rng, 1:4))
    npkts = Int[]
    subtps = Int[]
    degres = Int[]
    packts = Float64[]
    rates = Float64[]
    sclkdp = Float64[]
    ivlbds = Float64[Float64(rand(rng, 10:5000))]
    for _ in 1:nmini
        n = nmini > 10 ? rand(rng, 2:4) : rand(rng, (2, 5, 100, 101, 230))
        sub = rand(rng, 0:3)
        deg = isodd(sub) ? rand(rng, 1:2:23) : rand(rng, (3, 7, 11, 15, 19, 23))
        t = ivlbds[end] - rand(rng, (0.0, 4.0)) .+ cumsum(Float64[0; rand(rng, 1:400, n - 1)])
        q = rand_quats(rng, n)
        pk = sub == 0 ? vcat(q, 1e-4 .* randn(rng, 4, n)) :
             sub == 1 ? q :
             sub == 2 ? vcat(q, 1e-4 .* randn(rng, 4, n), 1e-3 .* randn(rng, 3, n), 1e-6 .* randn(rng, 3, n)) :
             vcat(q, 1e-3 .* randn(rng, 3, n))
        push!(npkts, n)
        push!(subtps, sub)
        push!(degres, deg)
        push!(rates, 1 / 256)
        append!(packts, vec(pk))
        append!(sclkdp, t)
        push!(ivlbds, t[end])
    end
    return (; ck_ids(k)..., avflag=rand(rng, Bool), first=ivlbds[1], last=ivlbds[end],
        segid=rand_segid(rng, "C06"), npkts, subtps, degres, packts, rates, sclkdp, ivlbds,
        sellst=rand(rng, Bool))
end
function bad_ckw06(a)
    swapped = copy(a.ivlbds)
    swapped[2] = swapped[1]
    dup = copy(a.sclkdp)
    dup[2] = dup[1]
    zq = copy(a.packts)
    zq[1:4] .= 0.0
    late = copy(a.sclkdp)
    late[1] = a.ivlbds[1] + 0.5
    out = Any[bad_ids(a; fname=:ref, empty_both=false)...,
        ("first after last", merge(a, (first=a.last + 1.0,)), :baddescrtimes, true),
        ("bounds out of order", merge(a, (ivlbds=swapped,)), :boundsoutoforder, true),
        ("first before the bounds", merge(a, (first=a.ivlbds[1] - 1.0,)), :coveragegap, true),
        ("last after the bounds", merge(a, (last=a.ivlbds[end] + 1.0,)), :coveragegap, true),
        ("one packet", merge(a, (npkts=[1; a.npkts[2:end]],)), :toofewpackets, true),
        ("subtype 4", merge(a, (subtps=[4; a.subtps[2:end]],)), :invalidsubtype, true),
        ("zero rate", merge(a, (rates=[0.0; a.rates[2:end]],)), :invalidsclkrate, true),
        ("degree 25", merge(a, (subtps=[a.subtps[1] | 1; a.subtps[2:end]], degres=[25; a.degres[2:end]])), :invaliddegree, true),
        ("odd window size", merge(a, (degres=[isodd(a.subtps[1]) ? 2 : 5; a.degres[2:end]],)), :badwindowsize, true),
        ("repeated time", merge(a, (sclkdp=dup,)), :timesoutoforder, true),
        ("no mini-segment", merge(a, (npkts=Int[], subtps=Int[], degres=Int[], rates=Float64[])), :invalidcount, false),
    ]
    if a.sclkdp[1] == a.ivlbds[1] && a.sclkdp[2] > late[1]
        push!(out, ("first time after the interval start", merge(a, (sclkdp=late,)), :boundsdisagree, true))
    end
    return out
end

# --- binary PCK type 2 ----------------------------------------------------------------------

"Body-fixed frames of the PCK segments: frame class ID => frame name (built-in class 2 frames)."
const PCK_FRAMES = Dict(301 => "IAU_MOON", 499 => "IAU_MARS", 3000 => "ITRF93")

function gen_pckw02(rng, k)
    n = rand(rng, (1, 2, 9, 40))
    polydg = rand(rng, (0, 1, 5, 12))
    intlen = 86400.0 * (1 + 9 * rand(rng))
    btime = 1.0e7 * randn(rng)
    cdata = randn(rng, 3 * (polydg + 1), n)         # angles in rad
    first, last = rand(rng) < 0.5 ? (btime, btime + n * intlen) : (btime + 0.1 * intlen, btime + (n - 0.1) * intlen)
    return (clssid=rand(rng, sort!(collect(keys(PCK_FRAMES)))),
        frame=rand(rng, ("J2000", "ECLIPJ2000", "B1950", "GALACTIC")), first, last,
        segid=rand_segid(rng, "P02"), intlen, cdata, btime)
end
function bad_pckw02(a)
    return Any[bad_ids(a; long=:idstringtoolong)...,
        ("frame that is not inertial", merge(a, (frame="IAU_EARTH",)), :invalidrefframe, true),
        ("no record", merge(a, (cdata=zeros(size(a.cdata, 1), 0),)), :numcoeffsnotpos, true),
        ("zero interval length", merge(a, (intlen=0.0,)), :intlennotpos, true),
        ("first after last", merge(a, (first=a.last + 1.0,)), :baddescrtimes, true),
        ("first before the data", merge(a, (first=a.btime - 1.0e-3,)), :baddescrtimes, true),
        ("last after the data", merge(a, (last=a.btime + size(a.cdata, 2) * a.intlen + 1.0e-3,)), :baddescrtimes, true),
        ("record size not a multiple of 3", merge(a, (cdata=zeros(4, size(a.cdata, 2)),)), :badarg, false),
    ]
end

const WRITER_CASES = WriterCase[
    WriterCase("spkw13", :spk, gen_w13,
        (h, a) -> spkw13(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.degree, a.states, a.epochs),
        (h, a) -> SPICE.spkw13(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.degree, cols(a.states), a.epochs),
        bad_w13),
    WriterCase("spkw09", :spk, gen_w09,
        (h, a) -> spkw09(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.degree, a.states, a.epochs),
        (h, a) -> SPICE.spkw09(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.degree, cols(a.states), a.epochs),
        bad_w09),
    WriterCase("spkw05", :spk, gen_w05,
        (h, a) -> spkw05(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.gm, a.states, a.epochs),
        (h, a) -> SPICE.spkw05(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.gm, cols(a.states), a.epochs),
        bad_w05),
    WriterCase("spkw08", :spk, (rng, k) -> gen_w0812(rng, k, false),
        (h, a) -> spkw08(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.degree, a.states, a.epoch1, a.step),
        (h, a) -> SPICE.spkw08(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.degree, cols(a.states), a.epoch1, a.step),
        a -> bad_w0812(a, false)),
    WriterCase("spkw12", :spk, (rng, k) -> gen_w0812(rng, k, true),
        (h, a) -> spkw12(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.degree, a.states, a.epoch1, a.step),
        (h, a) -> SPICE.spkw12(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.degree, cols(a.states), a.epoch1, a.step),
        a -> bad_w0812(a, true)),
    WriterCase("spkw18", :spk, gen_w18,
        (h, a) -> spkw18(h, a.subtyp, a.body, a.center, a.frame, a.first, a.last, a.segid, a.degree, a.packts, a.epochs),
        (h, a) -> SPICE.spkw18(h, (:S18TP0, :S18TP1)[a.subtyp+1], a.body, a.center, a.frame, a.first, a.last, a.segid, a.degree, cols(a.packts), a.epochs),
        bad_w18),
    WriterCase("spkw02", :spk, (rng, k) -> gen_w0203(rng, k, 3),
        (h, a) -> spkw02(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.intlen, a.cdata, a.btime),
        (h, a) -> c_spkw0203(:spkw02_c, h, a), a -> bad_w0203(a, 3)),
    WriterCase("spkw03", :spk, (rng, k) -> gen_w0203(rng, k, 6),
        (h, a) -> spkw03(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.intlen, a.cdata, a.btime),
        (h, a) -> c_spkw0203(:spkw03_c, h, a), a -> bad_w0203(a, 6)),
    WriterCase("spkw20", :spk, gen_w20,
        (h, a) -> spkw20(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.intlen, a.cdata, a.dscale, a.tscale, a.initjd, a.initfr),
        (h, a) -> SPICE.spkw20(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.intlen, cols(a.cdata), a.dscale, a.tscale, a.initjd, a.initfr),
        bad_w20),
    WriterCase("spkw15", :spk, gen_w15,
        (h, a) -> spkw15(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.epoch, a.tp, a.pa, a.p, a.ecc, a.j2flg, a.pv, a.gm, a.j2, a.radius),
        (h, a) -> SPICE.spkw15(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.epoch, a.tp, a.pa, a.p, a.ecc, a.j2flg, a.pv, a.gm, a.j2, a.radius),
        bad_w15),
    WriterCase("spkw17", :spk, gen_w17,
        (h, a) -> spkw17(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.epoch, a.eqel, a.rapol, a.decpol),
        (h, a) -> SPICE.spkw17(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.epoch, a.eqel, a.rapol, a.decpol),
        bad_w17),
    WriterCase("spkw01", :spk, (rng, k) -> gen_w0121(rng, k, false),
        (h, a) -> spkw01(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.dlines, a.epochs),
        c_spkw01, a -> bad_w0121(a, false)),
    WriterCase("spkw21", :spk, (rng, k) -> gen_w0121(rng, k, true),
        (h, a) -> spkw21(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.dlines, a.epochs),
        c_spkw21, a -> bad_w0121(a, true)),
    WriterCase("spk14", :spk, gen_w14, n_spk14, c_spk14, bad_w14),
    WriterCase("spkw19", :spk, gen_w19,
        (h, a) -> spkw19(h, a.body, a.center, a.frame, a.first, a.last, a.segid, a.npkts, a.subtps, a.degres, a.packts, a.epochs, a.ivlbds, a.sellst),
        c_spkw19, bad_w19),
    WriterCase("ckw01", :ck, gen_ckw01,
        (h, a) -> a.avvs === nothing ? ckw01(h, a.begtim, a.endtim, a.inst, a.ref, a.segid, a.sclkdp, a.quats) :
                  ckw01(h, a.begtim, a.endtim, a.inst, a.ref, a.segid, a.sclkdp, a.quats, a.avvs),
        c_ckw01, bad_ckw01),
    WriterCase("ckw02", :ck, gen_ckw02,
        (h, a) -> ckw02(h, a.begtim, a.endtim, a.inst, a.ref, a.segid, a.start, a.stop, a.quats, a.avvs, a.rates),
        c_ckw02, bad_ckw02),
    WriterCase("ckw03", :ck, gen_ckw03,
        (h, a) -> a.avvs === nothing ? ckw03(h, a.begtim, a.endtim, a.inst, a.ref, a.segid, a.sclkdp, a.quats, a.starts) :
                  ckw03(h, a.begtim, a.endtim, a.inst, a.ref, a.segid, a.sclkdp, a.quats, a.starts, a.avvs),
        c_ckw03, bad_ckw03),
    WriterCase("ckw05", :ck, gen_ckw05,
        (h, a) -> ckw05(h, a.subtyp, a.degree, a.begtim, a.endtim, a.inst, a.ref, a.avflag, a.segid, a.sclkdp, a.packts, a.rate, a.nints, a.starts),
        c_ckw05, bad_ckw05),
    WriterCase("ckw06", :ck, gen_ckw06,
        (h, a) -> ckw06(h, a.inst, a.ref, a.avflag, a.first, a.last, a.segid, a.npkts, a.subtps, a.degres, a.packts, a.rates, a.sclkdp, a.ivlbds, a.sellst),
        c_ckw06, bad_ckw06),
    WriterCase("pckw02", :pck, gen_pckw02,
        (h, a) -> pckw02(h, a.clssid, a.frame, a.first, a.last, a.segid, a.intlen, a.cdata, a.btime),
        c_pckw02, bad_pckw02),
]

"The writer case with the name `name`."
writer_case(name::AbstractString) = WRITER_CASES[findfirst(c -> c.name == name, WRITER_CASES)]

neo_open(kind::Symbol, path, ifname, ncomch=0) = kind === :spk ? spkopn(path, ifname, ncomch) :
    kind === :ck ? ckopn(path, ifname, ncomch) : pckopn(path, ifname, ncomch)
c_open(kind::Symbol, path, ifname, ncomch=0) = kind === :spk ? SPICE.spkopn(path, ifname, ncomch) :
    kind === :ck ? SPICE.ckopn(path, ifname, ncomch) : SPICE.pckopn(path, ifname, ncomch)
neo_close(kind::Symbol, h) = kind === :spk ? spkcls(h) : kind === :ck ? ckcls(h) : pckcls(h)
c_close(kind::Symbol, h) = kind === :spk ? SPICE.spkcls(h) : kind === :ck ? SPICE.ckcls(h) : SPICE.pckcls(h)
"File name extension and file type (as `getfat` gives it) of the kernel kind `kind`."
kind_ext(kind::Symbol) = kind === :spk ? ".bsp" : kind === :ck ? ".bc" : ".bpc"
kind_type(kind::Symbol) = kind === :spk ? "SPK" : kind === :ck ? "CK" : "PCK"

"Write the segments `segs` of the case `c` with both engines; returns the paths `(pc, pn)`."
function write_both(c::WriterCase, dir::AbstractString, tag::AbstractString, segs; ncomch::Integer=0,
        comments=nothing)
    ext = kind_ext(c.kind)
    pc = joinpath(dir, "c_" * tag * ext)
    pn = joinpath(dir, "n_" * tag * ext)
    hc = c_open(c.kind, pc, "IFNAME " * tag, ncomch)
    comments === nothing || SPICE.dafac(hc, comments)
    foreach(a -> c.cspice(hc, a), segs)
    c_close(c.kind, hc)
    hn = neo_open(c.kind, pn, "IFNAME " * tag, ncomch)
    comments === nothing || dafac(hn, comments)
    foreach(a -> c.neo(hn, a), segs)
    neo_close(c.kind, hn)
    return pc, pn
end

# ---------------------------------------------------------------------------------------
# What CSPICE reads back
# ---------------------------------------------------------------------------------------

"Run `f()`; a CSPICE error becomes the symbol `:error` (the same call must fail on both files)."
function c_value(f)
    try
        return f()
    catch e
        e isa SPICE.SpiceError || rethrow()
        return :error
    end
end

"States that CSPICE computes from the SPK file `path` for the segments `segs` at the epochs `ets[k]`."
function c_read_spk(path::AbstractString, segs, ets)
    SPICE.furnsh(path)
    try
        return [c_value(() -> SPICE.spkgeo(a.body, et, "J2000", a.center)) for (a, e) in zip(segs, ets) for et in e]
    finally
        SPICE.unload(path)
    end
end

"Pointing that CSPICE computes from the CK file `path` (with the SCLK kernel `sclk`) at the ticks `ticks[k]`."
function c_read_ck(path::AbstractString, sclk::AbstractString, segs, ticks)
    SPICE.furnsh(sclk)
    SPICE.furnsh(path)
    try
        return [(c_value(() -> SPICE.ckgpav(a.inst, t, 0.0, "J2000")), c_value(() -> SPICE.ckgp(a.inst, t, 0.0, "J2000")))
                for (a, ts) in zip(segs, ticks) for t in ts]
    finally
        SPICE.unload(path)
        SPICE.unload(sclk)
    end
end

"Rotations from J2000 to the body-fixed frames that CSPICE computes from the binary PCK `path` for the segments `segs` at the epochs `ets[k]`."
function c_read_pck(path::AbstractString, segs, ets)
    SPICE.furnsh(path)
    try
        return [c_value(() -> SPICE.pxform("J2000", PCK_FRAMES[a.clssid], et)) for (a, e) in zip(segs, ets) for et in e]
    finally
        SPICE.unload(path)
    end
end

"Epochs at which the SPK segment `a` is read back: the descriptor bounds and random epochs between them."
spk_probe_epochs(rng, a) = Float64[a.first; a.last; [a.first + (a.last - a.first) * rand(rng) for _ in 1:6]]

"Ticks at which the CK segment `a` is read back: time tags of the segment and points between them."
function ck_probe_ticks(rng, a)
    tags = hasproperty(a, :sclkdp) ? a.sclkdp : a.start
    pick = tags[unique!([1; length(tags); rand(rng, 1:length(tags), 4)])]
    return Float64[pick; [tags[1] + (tags[end] - tags[1]) * rand(rng) for _ in 1:3]]
end

"Hand-written type 1 SCLK kernel for the spacecraft -77: 256 ticks in a second, one partition, TDB."
function write_sclk(dir::AbstractString)
    path = joinpath(dir, "allspice_test.tsc")
    write(path, """
    KPL/SCLK

    Spacecraft clock of the test spacecraft -77 (written by test_dafwrite.jl).

    \\begindata

       SCLK_KERNEL_ID            = ( @2026-10-11/00:00:00 )
       SCLK_DATA_TYPE_77         = ( 1 )
       SCLK01_TIME_SYSTEM_77     = ( 1 )
       SCLK01_N_FIELDS_77        = ( 2 )
       SCLK01_MODULI_77          = ( 4294967296 256 )
       SCLK01_OFFSETS_77         = ( 0 0 )
       SCLK01_OUTPUT_DELIM_77    = ( 1 )
       SCLK_PARTITION_START_77   = ( 0.0000000000000E+00 )
       SCLK_PARTITION_END_77     = ( 1.0995116277750E+12 )
       SCLK01_COEFFICIENTS_77    = ( 0.0000000000000E+00  0.0000000000000E+00  1.0000000000000E+00 )

    \\begintext
    """)
    return path
end

"Random comment lines: printable ASCII, some empty, some with leading or trailing blanks."
function rand_comments(rng::AbstractRNG, n::Integer; maxlen::Integer=120)
    chars = Char.(32:126)
    return String[rand(rng) < 0.15 ? "" : (" "^rand(rng, 0:3)) * randstring(rng, chars, rand(rng, 0:maxlen)) *
                  (" "^rand(rng, 0:2)) for _ in 1:n]
end

"What `dafec` gives for the lines `lines` that `dafac` stored: trailing blanks removed."
stored(lines) = String[String(rstrip(l, ' ')) for l in lines]

"""
Copy of the native DAF `src` in the other byte order at `dst`: the integers of the file
record, the control words and summaries of the summary records and all data doubles are
reversed; the comment records and the name records are bytes and stay.
"""
function write_other_endian(src::AbstractString, dst::AbstractString)
    raw = read(src)
    out = copy(raw)
    rev!(off, n) = reverse!(out, off + 1, off + n)
    foreach(off -> rev!(off, 4), (8, 12, 76, 80, 84))
    out[89:96] = codeunits(ENDIAN_BOM == 0x04030201 ? "BIG-IEEE" : "LTL-IEEE")
    nd, ni, fward = Int.(reinterpret(Int32, raw[[9:12; 13:16; 77:80]]))
    ss = nd + (ni + 1) ÷ 2
    chars = Set{Int}(2:fward-1)
    rec = fward
    while rec != 0
        base = (rec - 1) * 1024
        ctl = reinterpret(Float64, raw[base+1:base+24])
        foreach(w -> rev!(base + 8w, 8), 0:2)
        for s in 0:Int(ctl[3])-1
            off = base + 24 + 8 * ss * s
            foreach(j -> rev!(off + 8j, 8), 0:nd-1)
            foreach(j -> rev!(off + 8nd + 4j, 4), 0:ni-1)
        end
        push!(chars, rec, rec + 1)
        rec = Int(ctl[1])
    end
    for r in 2:length(raw)÷1024
        r in chars && continue
        foreach(w -> rev!((r - 1) * 1024 + 8w, 8), 0:127)
    end
    write(dst, out)
    return dst
end

# ---------------------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------------------

group("9. DAF access, SPK and CK writers") do
    rng = MersenneTwister(SEED + 9)
    dir = mktempdir()
    SPICE.kclear()
    sclk = write_sclk(dir)

    @apart @testset "writers: data area, descriptors, names and read-back against CSPICE" begin
        # Files of each writer. Three files of one to three segments reach each layout
        # branch of each writer in the default run (measured with line coverage).
        ntrial = nsample(3, 6)
        measured(@sprintf("  %-8s %5s %9s %9s %5s %5s %6s %8s", "writer", "segs", "doubles", "bit diff",
            "descr", "names", "bytes", "readback"))
        for c in WRITER_CASES
            nseg = 0
            words = 0
            badwords = 0
            descr = true
            names = true
            bytes = true
            nread = 0
            readok = true
            for trial in 1:ntrial
                segs = [c.gen(rng, k) for k in 1:rand(rng, 1:3)]
                pc, pn = write_both(c, dir, "$(c.name)_$trial", segs)
                r = compare_files(pc, pn)
                nseg += r.nseg
                words += r.words
                badwords += r.badwords
                descr &= r.sums && r.record
                names &= r.names
                bytes &= r.bytes
                @test r.nseg == length(segs)
                if c.kind === :spk
                    ets = [spk_probe_epochs(rng, a) for a in segs]
                    rc = c_read_spk(pc, segs, ets)
                    rn = c_read_spk(pn, segs, ets)
                elseif c.kind === :pck
                    ets = [spk_probe_epochs(rng, a) for a in segs]
                    rc = c_read_pck(pc, segs, ets)
                    rn = c_read_pck(pn, segs, ets)
                else
                    ticks = [ck_probe_ticks(rng, a) for a in segs]
                    rc = c_read_ck(pc, sclk, segs, ticks)
                    rn = c_read_ck(pn, sclk, segs, ticks)
                end
                readok &= isequal(rc, rn)
                nread += count(x -> c.kind !== :ck ? x !== :error : (x[2] !== :error && x[2] !== nothing), rc)
            end
            @testset "$(c.name)" begin
                @test badwords == 0
                @test descr
                @test names
                @test bytes
                @test readok
                @test nread > 0
            end
            measured(@sprintf("  %-8s %5d %9d %9d %5s %5s %6s %8s", c.name, nseg, words, badwords,
                descr ? "same" : "DIFF", names ? "same" : "DIFF", bytes ? "same" : "DIFF",
                readok ? "same/$nread" : "DIFF"))
        end
    end

    @apart @testset "writers: forms of the table arguments" begin
        c = writer_case("spkw13")
        a = c.gen(rng, 1)
        paths = String[]
        for (i, states) in enumerate((a.states, cols(a.states), vec(a.states), [Tuple(col) for col in cols(a.states)],
            [SVector{6}(col) for col in cols(a.states)], view(a.states, :, :)))
            p = joinpath(dir, "forms_$i.bsp")
            h = spkopn(p, "FORMS")
            c.neo(h, merge(a, (states=states, epochs=i == 2 ? view(a.epochs, :) : a.epochs, frame=i == 3 ? 1 : "J2000")))
            spkcls(h)
            push!(paths, p)
        end
        @test all(p -> read(p) == read(paths[1]), paths)
        # A frame that only the kernel pool of an explicit context defines.
        fk = joinpath(dir, "writer_frame.tf")
        write(fk, """
        KPL/FK
        \\begindata
           FRAME_WRITER_TEST_FRAME   = 1977001
           FRAME_1977001_NAME        = 'WRITER_TEST_FRAME'
           FRAME_1977001_CLASS       = 4
           FRAME_1977001_CLASS_ID    = 1977001
           FRAME_1977001_CENTER      = 399
           TKFRAME_1977001_SPEC      = 'MATRIX'
           TKFRAME_1977001_RELATIVE  = 'J2000'
           TKFRAME_1977001_MATRIX    = ( 1 0 0 0 1 0 0 0 1 )
        \\begintext
        """)
        ctx = neo_context([fk])
        p = joinpath(dir, "frame_ctx.bsp")
        h = spkopn(p, "FRAME")
        b = merge(a, (frame="WRITER_TEST_FRAME",))
        @test throws_kind(() -> c.neo(h, b), :invalidrefframe)      # not in the default context
        spkw13(h, b.body, b.center, b.frame, b.first, b.last, b.segid, b.degree, b.states, b.epochs; ctx=ctx)
        spkw13(h, b.body, b.center, b.frame, b.first, b.last, b.segid, b.degree, b.states, b.epochs;
            ctx=Allspice.kernelset(ctx))
        spkcls(h)
        d = neo_dump(p)
        @test [s[2][3] for s in d.sums] == [1977001, 1977001]
        SPICE.furnsh(fk)
        pc = joinpath(dir, "frame_ctx_c.bsp")
        hc = SPICE.spkopn(pc, "FRAME", 0)
        c.cspice(hc, b)
        c.cspice(hc, b)
        SPICE.spkcls(hc)
        SPICE.unload(fk)
        @test read(pc) == read(p)
        # (fault fixed: without the `avvs` argument a segment of one record got the angular
        # velocity zero, and a segment of more records got none; the default is now no
        # angular velocity for every count. An explicit argument keeps its meaning.)
        p1 = joinpath(dir, "ck_default_avvs.bc")
        h = ckopn(p1, "AVVS", 0)
        q1 = [1.0, 0.0, 0.0, 0.0]
        ckw01(h, 5.0, 5.0, -77001, "J2000", "one record, default", [5.0], [q1])
        ckw03(h, 5.0, 5.0, -77001, "J2000", "one record, default", [5.0], [q1], [5.0])
        ckw01(h, 5.0, 5.0, -77001, "J2000", "one record, one vector", [5.0], [q1], [zeros(3)])
        ckw01(h, 5.0, 9.0, -77001, "J2000", "two records, one vector", [5.0, 9.0], [q1, q1], [zeros(3)])
        ckcls(h)
        @test [s[2][4] for s in neo_dump(p1).sums] == [0, 0, 1, 0]
        # Only the writers that call `spkpds` in CSPICE reject the barycenter as body and
        # a body equal to its center (`bad_pds`); `spkw02` is one of those that accept them.
        c2 = writer_case("spkw02")
        pc2, pn2 = write_both(c2, dir, "body_is_center", [merge(c2.gen(rng, 1), (body=0, center=0))])
        @test read(pc2) == read(pn2)
    end

    @apart @testset "writers: every input check" begin
        for c in WRITER_CASES
            @testset "$(c.name)" begin
                a = c.gen(rng, 1)
                # A case with at least three records, so that every variant can be built. The
                # writers of types 8 and 12 have no epoch vector, thus the count of their
                # records is the count of the columns of `states`.
                while !all(v -> !(v isa AbstractVector{Float64}) || length(v) >= 3, values(a)) ||
                    (hasproperty(a, :states) && size(a.states, 2) < 3) ||
                    (hasproperty(a, :starts) && length(a.starts) < 2)
                    a = c.gen(rng, 1)
                end
                ext = kind_ext(c.kind)
                pn = joinpath(dir, "bad_n_$(c.name)$ext")
                pc = joinpath(dir, "bad_c_$(c.name)$ext")
                pref = joinpath(dir, "bad_ref_$(c.name)$ext")
                hn = neo_open(c.kind, pn, "BAD")
                hc = c_open(c.kind, pc, "BAD")
                variants = c.bad(a)
                @test length(variants) >= 8
                for (label, b, kind, both) in variants
                    kind === :none && continue
                    ok = throws_kind(() -> c.neo(hn, b), kind)
                    ok || @info "input check failed" writer = c.name variant = label expected = kind
                    @test ok
                    if both
                        cok = c_value(() -> (c.cspice(hc, b); :written)) === :error
                        cok || @info "CSPICE accepts the variant" writer = c.name variant = label
                        @test cok
                    end
                end
                # A rejected call leaves nothing in the file: the file with the rejected calls
                # and one valid segment is the file with the valid segment only.
                c.neo(hn, a)
                neo_close(c.kind, hn)
                c.cspice(hc, a)
                c_close(c.kind, hc)
                href = neo_open(c.kind, pref, "BAD")
                c.neo(href, a)
                neo_close(c.kind, href)
                @test read(pn) == read(pref)
                @test read(pn) == read(pc)
                # A read handle and a closed handle are rejected.
                hr = dafopr(pn)
                @test throws_kind(() -> c.neo(hr, a), :dafinvalidaccess)
                dafcls(hr)
                @test throws_kind(() -> c.neo(hr, a), :dafnosuchhandle)
                @test read(pn) == read(pref)
            end
        end
    end

    @apart @testset "spk14: the writer object" begin
        a = writer_case("spk14").gen(rng, 3)
        p = joinpath(dir, "w14_object.bsp")
        h = neo_open(:spk, p, "W14")
        w = spk14b(h, a.segid, a.body, a.center, a.frame, a.first, a.last, a.degree)
        @test throws_kind(() -> spk14a(w, 0, Float64[], Float64[]), :invalidargument)
        @test occursin("0 records", sprint(show, w))
        spk14a(w, a.coeffs, a.epochs)
        spk14e(w)
        @test occursin("ended", sprint(show, w))
        @test throws_kind(() -> spk14e(w), :segmentnotopen)
        @test throws_kind(() -> spk14a(w, a.coeffs, a.epochs), :segmentnotopen)
        neo_close(:spk, h)
        # The segment that was written gives the series of its records.
        ctx = neo_context([p])
        i = length(a.epochs)
        mid, rad = a.coeffs[1, i], a.coeffs[2, i]
        st = spkgeo(ctx, a.body, mid, a.frame, a.center)[1]
        nc = a.degree + 1
        # At the midpoint the series is the sum of its even coefficients with changing sign.
        expect(k) = sum(a.coeffs[2+(k-1)*nc+j+1, i] * (iseven(j ÷ 2) ? 1.0 : -1.0) for j in 0:2:a.degree)
        @test all(isapprox(st[k], expect(k); rtol=1e-12, atol=1e-12) for k in 1:6)
        @test rad > 0
    end

    @apart @testset "many segments in one file: past one summary record, and appending" begin
        c13 = writer_case("spkw13")
        c03 = writer_case("ckw03")
        # A summary record holds 25 summaries. Default run: a full record, one summary in
        # the second record and one in the third record for SPK; the CK files have the
        # same summary format, so one count shows that `ckcls` closes such a file.
        plan = FULL ? ((c13, (24, 25, 26, 50, 51, 77)), (c03, (25, 26, 60))) : ((c13, (25, 26, 51)), (c03, (26,)))
        for (c, counts) in plan
            for nseg in counts
                segs = [merge(c.gen(rng, k), c.kind === :spk ? (states=rand_states(rng, 4), epochs=[0.0, 1.0, 2.0, 3.0],
                    first=0.0, last=3.0, degree=3) : (;)) for k in 1:nseg]
                c.kind === :ck && (segs = [a for a in segs if length(a.sclkdp) <= 101])
                pc, pn = write_both(c, dir, "many_$(c.name)_$nseg", segs)
                r = compare_files(pc, pn)
                @test r.nseg == length(segs)
                @test files_agree(r)
                @test r.bytes
                if c.kind === :spk
                    @test spkobj(pn) == sort!(unique!([a.body for a in segs]))
                    ets = [[a.first, 1.5, a.last] for a in segs]
                    @test isequal(c_read_spk(pc, segs, ets), c_read_spk(pn, segs, ets))
                end
            end
        end
        # Appending in the other engine gives the file that one engine writes alone.
        segs = [c13.gen(rng, k) for k in 1:57]
        pall, nall = write_both(c13, dir, "append_all", segs)
        pcn, pnc = write_both(c13, dir, "append_part", segs[1:23])
        hn = spkopa(pcn)                             # CSPICE began the file, Allspice appends
        foreach(a -> c13.neo(hn, a), segs[24:40])
        spkcls(hn)
        hn = dafopw(pcn)
        foreach(a -> c13.neo(hn, a), segs[41:57])
        dafcls(hn)
        hc = SPICE.spkopa(pnc)                       # Allspice began the file, CSPICE appends
        foreach(a -> c13.cspice(hc, a), segs[24:57])
        SPICE.spkcls(hc)
        for p in (pcn, pnc)
            r = compare_files(pall, p)
            @test r.nseg == 57 && r.badwords == 0 && r.sums && r.names && r.comments
        end
        # The internal file names differ by the tag only.
        strip_ifname(p) = (b = read(p); b[17:76] .= 0x20; b)
        @test strip_ifname(pcn) == strip_ifname(pall)
        @test strip_ifname(pnc) == strip_ifname(pall)
        @test read(nall) == read(pall)
        @test throws_kind(() -> spkopa(joinpath(dir, "no_such_file.bsp")), :filenotfound)
        @test throws_kind(() -> spkopa(sclk), :fileisnotspk)
        ckfile = write_both(c03, dir, "notspk", [c03.gen(rng, 1)])[2]
        @test throws_kind(() -> spkopa(ckfile), :fileisnotspk)
        @test_throws SPICE.SpiceError SPICE.spkopa(ckfile)
    end

    @apart @testset "comments: dafac and dafec in both directions" begin
        c = writer_case("spkw13")
        seg(k) = merge(c.gen(rng, k), (states=rand_states(rng, 30), epochs=collect(0.0:29.0), first=0.0, last=29.0, degree=3))
        nbytes = 0
        nfiles = 0
        # (reserved characters, comments before the segments, lines of 1000 characters).
        # In the last case the first line and its end mark fill one comment record exactly.
        # A comment record holds 1000 characters: no record, one full record, one character
        # in the second record, and several records. The exhaustive run adds more sizes.
        cases = [(0, true, false), (1000, true, false), (1001, true, false), (2500, true, false),
            (0, false, false), (0, false, true)]
        FULL && append!(cases, [(1, true, false), (999, true, false), (6000, true, false), (3000, true, false)])
        for (i, (ncomch, early, longlines)) in enumerate(cases)
            before = early ? rand_comments(rng, rand(rng, 1:30)) : String[]
            after = longlines ? ["x"^999, "y"^1000, "", "z"^1000, "w"] : rand_comments(rng, rand(rng, 1:60); maxlen=400)
            last = ["one more line", ""]
            segs = [seg(k) for k in 1:rand(rng, (1, 3, 30))]
            pc = joinpath(dir, "com_c_$i.bsp")
            pn = joinpath(dir, "com_n_$i.bsp")
            # The same sequence with both engines: comments, segments, comments (the comment
            # area grows and the data move), close, open again, comments.
            hc = SPICE.spkopn(pc, "COMMENTS", ncomch)
            isempty(before) || SPICE.dafac(hc, before)
            foreach(a -> c.cspice(hc, a), segs)
            SPICE.dafac(hc, after)
            SPICE.spkcls(hc)
            hc = SPICE.dafopw(pc)
            SPICE.dafac(hc, last)
            SPICE.dafcls(hc)
            hn = spkopn(pn, "COMMENTS", ncomch)
            isempty(before) || dafac(hn, before)
            foreach(a -> c.neo(hn, a), segs)
            dafac(hn, after)
            spkcls(hn)
            hn = dafopw(pn)
            dafac(hn, last)
            @test dafec(hn) == stored([before; after; last])       # readable on the write handle
            dafcls(hn)
            r = compare_files(pc, pn)
            @test files_agree(r)
            @test r.nseg == length(segs)
            @test neo_dump(pn).comments == stored([before; after; last])
            nfiles += 1
            nbytes += r.bytes
            # Cross direction: each engine appends to the file of the other.
            hc = SPICE.dafopw(pn)
            SPICE.dafac(hc, ["appended by CSPICE", " to the Allspice file"])
            SPICE.dafcls(hc)
            hn = dafopw(pc)
            dafac(hn, ["appended by CSPICE", " to the Allspice file"])
            dafcls(hn)
            r2 = compare_files(pc, pn)
            @test files_agree(r2)
            @test length(r2.nseg == length(segs) ? neo_dump(pc).comments : String[]) ==
                  length(before) + length(after) + length(last) + 2
            ets = [[0.0, 14.5, 29.0] for _ in segs]
            @test isequal(c_read_spk(pc, segs, ets), c_read_spk(pn, segs, ets))
        end
        measured(@sprintf("  %-58s files=%d byte-identical=%d", "comments: files with a grown comment area", nfiles, nbytes))
        # Delete the comment area with both engines.
        segs = [seg(k) for k in 1:27]
        pc, pn = write_both(c, dir, "dafdc", segs; ncomch=2300, comments=rand_comments(rng, 20))
        hc = SPICE.dafopw(pc)
        SPICE.dafdc(hc)
        SPICE.dafcls(hc)
        hn = dafopw(pn)
        dafdc(hn)
        @test dafec(hn) == String[]
        @test dafrfr(hn)[4] == 2
        dafcls(hn)
        r = compare_files(pc, pn)
        @test files_agree(r) && r.nseg == 27
        ets = [[0.0, 14.5, 29.0] for _ in segs]
        @test isequal(c_read_spk(pc, segs, ets), c_read_spk(pn, segs, ets))
        # Input checks.
        p = joinpath(dir, "dafac_checks.bsp")
        h = spkopn(p, "CHECKS", 100)
        @test throws_kind(() -> dafac(h, String[]), :invalidargument)
        @test throws_kind(() -> dafac(h, ["tab\there"]), :illegalcharacter)
        @test throws_kind(() -> dafac(h, ["ok", "é"]), :illegalcharacter)
        @test_throws SPICE.SpiceError (hc = SPICE.dafopw(pc); try SPICE.dafac(hc, ["tab\there"]) finally SPICE.dafcls(hc) end)
        @test dafec(h) == String[]
        dafbna(h, dafps([0.0, 1.0], [-5, 399, 1, 9, 0, 0]), "OPEN ARRAY")
        @test throws_kind(() -> dafac(h, ["while an array is open"]), :dafnewconflict)
        @test throws_kind(() -> dafdc(h), :dafnewconflict)
        dafcls(h)
        hr = dafopr(p)
        @test throws_kind(() -> dafac(hr, ["read only"]), :dafinvalidaccess)
        @test throws_kind(() -> dafdc(hr), :dafinvalidaccess)
        dafcls(hr)
    end

    @apart @testset "empty files, open and close" begin
        # The reserved comment records are the layout branch: none, one, two. The CK file
        # differs from the SPK file in the identification word only, the binary PCK file
        # also in the summary size.
        plan = FULL ? [(k, n) for k in (:spk, :ck, :pck) for n in (0, 1, 1000, 1001, 5000)] :
               [(:spk, 0), (:spk, 1000), (:spk, 1001), (:ck, 0), (:pck, 1001)]
        for (kind, ncomch) in plan
            ext = kind_ext(kind)
            pc = joinpath(dir, "empty_c_$(kind)_$ncomch$ext")
            pn = joinpath(dir, "empty_n_$(kind)_$ncomch$ext")
            hc = c_open(kind, pc, "EMPTY FILE", ncomch)
            SPICE.dafcls(hc)
            hn = neo_open(kind, pn, "EMPTY FILE", ncomch)
            @test isopen(hn)
            @test throws_kind(() -> neo_close(kind, hn), :nosegmentsfound)
            @test isopen(hn)
            dafcls(hn)
            @test !isopen(hn)
            dafcls(hn)                                         # a second close does nothing
            @test read(pc) == read(pn)
            @test getfat(pn) == SPICE.getfat(pn) == ("DAF", kind_type(kind))
            @test throws_kind(() -> neo_open(kind, pn, "EXISTS", 0), :fileopenfailed)
            @test read(pc) == read(pn)                         # the existing file is not touched
            @test_throws SPICE.SpiceError c_open(kind, pc, "EXISTS", 0)
        end
        @test throws_kind(() -> spkopn(joinpath(dir, "x.bsp"), ""), :emptystring)
        @test throws_kind(() -> spkopn("", "NAME"), :emptystring)
        @test throws_kind(() -> ckopn("", "NAME"), :emptystring)
        @test throws_kind(() -> pckopn(joinpath(dir, "x.bpc"), ""), :emptystring)
        @test throws_kind(() -> spkopn(joinpath(dir, "no_such_dir", "x.bsp"), "NAME"), :fileopenfailed)
        @test throws_kind(() -> dafopr(joinpath(dir, "no_such_file.bsp")), :filenotfound)
        @test throws_kind(() -> dafopw(joinpath(dir, "no_such_file.bsp")), :filenotfound)
        @test throws_kind(() -> dafopr(" "), :blankfilename)
        @test throws_kind(() -> dafopr(sclk), :filereadfailed)       # shorter than one record
        longtext = joinpath(dir, "long.tpc")
        write(longtext, "KPL/PCK\n" * "a line of text\n"^200)
        @test throws_kind(() -> dafopr(longtext), :unknownfilarc)   # as dafopr_c
        @test_throws SPICE.SpiceError SPICE.dafopr(longtext)
        # The binary file format: VAX is known but not supported; blanks are not known.
        for (bff, kind) in (("VAX-GFLT", :unsupportedbff), (" "^8, :unknownbff))
            bytes = read(joinpath(dir, "empty_c_spk_0.bsp"))
            bytes[89:96] = codeunits(bff)
            pb = joinpath(dir, "bff.bsp")
            write(pb, bytes)
            @test throws_kind(() -> dafopr(pb), kind)
        end
        junk = joinpath(dir, "junk.bin")
        write(junk, rand(rng, UInt8, 3000) .| 0x80)
        @test throws_kind(() -> dafopr(junk), :unknownfilarc)
        short = joinpath(dir, "short.bsp")
        write(short, "DAF/SPK ")
        @test throws_kind(() -> dafopr(short), :filereadfailed)
        @test throws_kind(() -> Allspice.dafonw(joinpath(dir, "p.daf"), " ", 2, 6, "X", 0), :blankfiletype)
        @test throws_kind(() -> Allspice.dafonw(joinpath(dir, "p.daf"), "SPK", 125, 6, "X", 0), :dafinvalidparams)
        @test throws_kind(() -> Allspice.dafonw(joinpath(dir, "p.daf"), "SPK", 2, 1, "X", 0), :dafinvalidparams)
        @test throws_kind(() -> Allspice.dafonw(joinpath(dir, "p.daf"), "SPK", 124, 6, "X", 0), :dafinvalidparams)
        @test throws_kind(() -> Allspice.dafonw(joinpath(dir, "p.daf"), "SPK", 2, 6, "X", -1), :dafnoresv)
        @test !ispath(joinpath(dir, "p.daf"))
    end

    @apart @testset "low-level DAF access" begin
        c = writer_case("spkw09")
        segs = [c.gen(rng, k) for k in 1:31]
        comments = rand_comments(rng, 12)
        pc, pn = write_both(c, dir, "lowlevel", segs; ncomch=1500, comments)
        h = dafopr(pn)
        hc = SPICE.dafopr(pc)
        @test dafrfr(h) == SPICE.dafrfr(hc)
        @test dafhsf(h) == (2, 6)
        @test dafec(h) == String.(SPICE.dafec(hc)) == stored(comments)
        # Forward and backward search give the same arrays in opposite order.
        fwd = Vector{Float64}[]
        fnames = String[]
        @test throws_kind(() -> daffna(h), :dafnosearch)
        @test throws_kind(() -> dafgs(h), :dafnosearch)
        dafbfs(h)
        @test throws_kind(() -> dafgs(h), :nocurrentarray)
        while daffna(h)
            push!(fwd, dafgs(h))
            push!(fnames, dafgn(h))
        end
        @test throws_kind(() -> dafgs(h), :nocurrentarray)
        @test !daffna(h)
        @test daffpa(h) && dafgs(h) == fwd[end]               # back from the end
        bwd = Vector{Float64}[]
        bnames = String[]
        dafbbs(h)
        while daffpa(h)
            push!(bwd, dafgs(h))
            push!(bnames, dafgn(h))
        end
        @test !daffpa(h)
        @test length(fwd) == 31 && bwd == reverse(fwd) && bnames == reverse(fnames)
        @test daffna(h) && dafgs(h) == fwd[1]
        @test fnames == [String(rstrip(a.segid)) for a in segs]
        cback = Vector{Float64}[]
        SPICE.dafbbs(hc)
        while SPICE.daffpa()
            push!(cback, SPICE.dafgs(5))
        end
        @test reinterpret(UInt64, reduce(vcat, cback)) == reinterpret(UInt64, reduce(vcat, bwd))
        # dafps and dafus against CSPICE.
        for s in fwd
            dc, ic = dafus(s, 2, 6)
            cdc, cic = SPICE.dafus(s, 2, 6)
            @test dc == cdc && ic == cic
            @test reinterpret(UInt64, dafps(dc, ic)) == reinterpret(UInt64, SPICE.dafps(dc, ic)) == reinterpret(UInt64, s)
        end
        @test dafus(dafps([1.5], [7, -8, 9]), 1, 3) == ([1.5], [7, -8, 9])        # odd integer count
        @test SPICE.dafus(dafps([1.5], [7, -8, 9]), 1, 3) == ([1.5], [7, -8, 9])
        @test throws_kind(() -> dafps([1.0], [2^40]), :dafinvalidparams)
        @test throws_kind(() -> dafus([1.0, 2.0], 2, 6), :dafinvalidparams)
        # dafgda and dafgsr.
        _, ic = dafus(fwd[7], 2, 6)
        for (b, e) in ((ic[5], ic[6]), (ic[5] + 3, ic[5] + 3), (ic[5] + 100, ic[5] + 420), (1, 128))
            @test reinterpret(UInt64, dafgda(h, b, e)) == reinterpret(UInt64, SPICE.dafgda(hc, b, e))
        end
        @test throws_kind(() -> dafgda(h, 0, 5), :dafnegaddr)
        @test throws_kind(() -> dafgda(h, 9, 8), :dafbeggtend)
        @test throws_kind(() -> dafgda(h, 1, 10^9), :filereadfailed)
        # A file that ends inside a record: the words that it holds are returned (CSPICE
        # gives zeros for the whole record), and a range past the end is an error.
        last = ic[5] + 10 + (iszero((ic[5] + 10) % 128) ? 1 : 0)
        cutfile = tempname()
        write(cutfile, read(h.path)[1:8*last])
        hx = dafopr(cutfile)
        @test reinterpret(UInt64, dafgda(hx, ic[5], last)) == reinterpret(UInt64, dafgda(h, ic[5], last))
        @test throws_kind(() -> dafgda(hx, last, last + 1), :filereadfailed)
        dafcls(hx)
        rm(cutfile; force=true)
        @test_throws SPICE.SpiceError SPICE.dafgda(hc, 0, 5)
        @test_throws SPICE.SpiceError ccall((:dafgda_c, LIBCSPICE), Cvoid, (Cint, Cint, Cint, Ptr{Cdouble}), hc, 9, 8, zeros(4)) === SPICE.handleerror()
        fward = dafrfr(h)[4]
        @test reinterpret(UInt64, dafgsr(h, fward, 1, 128)) == reinterpret(UInt64, SPICE.dafgsr(hc, fward, 1, 128))
        @test dafgsr(h, fward, 4, 8) == fwd[1]
        @test dafgsr(h, 10^6, 1, 128) === nothing
        @test throws_kind(() -> dafgsr(h, fward, 0, 128), :dafinvalidparams)
        # A read handle cannot write.
        @test throws_kind(() -> dafbna(h, fwd[1], "X"), :dafinvalidaccess)
        @test throws_kind(() -> dafws(h, fwd[1]), :dafinvalidaccess)
        dafcls(h)
        SPICE.dafcls(hc)
        @test throws_kind(() -> dafbfs(h), :dafnosuchhandle)
        @test throws_kind(() -> dafgda(h, 1, 2), :dafnosuchhandle)
        @test throws_kind(() -> dafec(h), :dafnosuchhandle)
        @test occursin("closed", sprint(show, h))

        # Change summaries and names with both engines (dafrs keeps the addresses, dafws does not).
        for (engine, p) in ((:c, pc), (:n, pn))
            if engine === :c
                hw = SPICE.dafopw(p)
                SPICE.dafbfs(hw)
                k = 0
                while SPICE.daffna()
                    k += 1
                    dc, ic = SPICE.dafus(SPICE.dafgs(5), 2, 6)
                    k % 3 == 0 && SPICE.dafrs(SPICE.dafps(dc .+ 0.5, [ic[1] - 5; ic[2:4]; 1; 2]))
                    k % 3 == 1 && ccall((:dafrn_, LIBCSPICE), Cint, (Ptr{UInt8}, Cint), "RENAMED $k", sizeof("RENAMED $k"))
                    k % 5 == 0 && ccall((:dafws_, LIBCSPICE), Cint, (Ptr{Cdouble},), SPICE.dafps(dc, [ic[1:3]; 13; ic[5:6]]))
                    SPICE.handleerror()
                end
                SPICE.dafcls(hw)
            else
                hw = dafopw(p)
                dafbfs(hw)
                k = 0
                while daffna(hw)
                    k += 1
                    dc, ic = dafus(dafgs(hw), 2, 6)
                    k % 3 == 0 && dafrs(hw, dafps(dc .+ 0.5, [ic[1] - 5; ic[2:4]; 1; 2]))
                    k % 3 == 1 && dafrn(hw, "RENAMED $k")
                    k % 5 == 0 && dafws(hw, dafps(dc, [ic[1:3]; 13; ic[5:6]]))
                end
                @test throws_kind(() -> dafws(hw, [1.0]), :nocurrentarray)
                dafcls(hw)
            end
        end
        r = compare_files(pc, pn)
        @test files_agree(r) && r.bytes && r.nseg == 31
        @test count(startswith("RENAMED"), neo_dump(pn).names) == 11

        # The generic array calls with another summary format (ND = 3, NI = 5).
        pg = joinpath(dir, "generic.daf")
        hg = Allspice.dafonw(pg, "TEST", 3, 5, "GENERIC DAF", 1)
        arrays = [randn(rng, n) for n in (1, 127, 128, 129, 1000, 0, 5)]
        for round in 1:4, (j, v) in enumerate(arrays)
            dafbna(hg, dafps([j, round, 0.5], [j, round, 7, 0, 0]), "ARRAY $j OF ROUND $round")
            dafada(hg, v[1:length(v)÷2])                           # the array grows in two steps
            dafada(hg, view(v, length(v)÷2+1:length(v)))
            dafada(hg, Float64[])
            dafena(hg)
        end
        @test throws_kind(() -> dafada(hg, [1.0]), :dafnewconflict)
        @test throws_kind(() -> dafena(hg), :dafnewconflict)
        @test throws_kind(() -> dafbna(hg, [1.0], "SHORT SUMMARY"), :dafinvalidparams)
        dafbna(hg, dafps([9, 9, 9.0], [9, 9, 9, 0, 0]), "ABANDONED")   # begun, never ended
        dafada(hg, [1.0, 2.0, 3.0])
        dafcls(hg)
        d = neo_dump(pg)
        cd = c_dump(pg)
        @test length(d.sums) == 28 && d.names[9] == "ARRAY 2 OF ROUND 2"
        @test d.data == repeat(arrays, 4)
        @test same_sums(d.sums, cd.sums) && d.names == cd.names && bit_mismatches(d.data, cd.data) == 0
        @test getfat(pg) == SPICE.getfat(pg) == ("DAF", "TEST")
        @test d.record[1:3] == (3, 5, "GENERIC DAF")

        # A file of the other byte order can be read and cannot be written.
        po = write_other_endian(pn, joinpath(dir, "other_endian.bsp"))
        @test getfat(po) == SPICE.getfat(po) == ("DAF", "SPK")
        dn = neo_dump(pn)
        dother = neo_dump(po)
        cother = c_dump(po)
        @test same_sums(dn.sums, dother.sums) && dn.names == dother.names && dn.comments == dother.comments
        @test bit_mismatches(dn.data, dother.data) == 0 && dn.record == dother.record
        @test same_sums(dn.sums, cother.sums) && bit_mismatches(dn.data, cother.data) == 0
        @test spkobj(po) == spkobj(pn) && collect(spkcov(po, dn.sums[1][2][1])) == collect(spkcov(pn, dn.sums[1][2][1]))
        @test throws_kind(() -> dafopw(po), :unsupportedbff)
        @test throws_kind(() -> spkopa(po), :unsupportedbff)
        @test_throws SPICE.SpiceError SPICE.dafopw(po)
        ho = dafopr(po)
        hn = dafopr(pn)
        @test dafgsr(ho, dafrfr(ho)[4], 1, 128) == dafgsr(hn, dafrfr(hn)[4], 1, 128)
        dafcls(ho)
        dafcls(hn)
        # A damaged FTP validation string is detected.
        damaged = read(pn)
        damaged[710] = 0x0a
        pd = joinpath(dir, "ftp_damaged.bsp")
        write(pd, damaged)
        @test throws_kind(() -> dafopr(pd), :ftpxfererror)
        @test_throws SPICE.SpiceError SPICE.dafopr(pd)
    end

    @apart @testset "getfat, spkobj, spkcov, tkvrsn" begin
        c = writer_case("spkw13")
        segs = [merge(c.gen(rng, k), (body=rand(rng, (-41, -42, -43, 501, 2000001)),)) for k in 1:70]
        pc, pn = write_both(c, dir, "coverage", segs)
        ids = spkobj(pn)
        @test ids == spk_bodies(pn) == spk_bodies(pc) == sort!(unique!([a.body for a in segs]))
        @test ids isa Vector{Int}
        for body in [ids; 999]
            w = spkcov(pn, body)
            cw = SPICE.spkcov!(SPICE.SpiceDoubleCell(400), pn, body)
            @test w isa Window
            @test collect(w) == Float64[cw[i] for i in 1:SPICE.card(cw)]
            @test card(w) == SPICE.card(cw)
        end
        # The mutating forms take the union with the old content.
        pc2, pn2 = write_both(c, dir, "coverage2", [merge(c.gen(rng, k), (body=-41,)) for k in 1:3])
        w = spkcov(pn, -41)
        @test spkcov!(w, pn2, -41) === w
        cw = SPICE.SpiceDoubleCell(400)
        SPICE.spkcov!(cw, pn, -41)
        SPICE.spkcov!(cw, pn2, -41)
        @test collect(w) == Float64[cw[i] for i in 1:SPICE.card(cw)]
        @test spkobj!([-9000, 7], pn2) == [-9000, -41, 7]
        ckfile = write_both(writer_case("ckw03"), dir, "cov_ck", [gen_ckw03(rng, 1)])[2]
        text = joinpath(dir, "plain.txt")
        write(text, "no kernel header here\n")
        for f in (spkobj, p -> spkcov(p, 1))
            @test throws_kind(() -> f(ckfile), :invalidfiletype)
            @test throws_kind(() -> f(sclk), :invalidarchtype)
            @test throws_kind(() -> f(text), :invalidarchtype)
            @test throws_kind(() -> f(joinpath(dir, "missing.bsp")), :filenotfound)
        end
        @test_throws SPICE.SpiceError SPICE.spkobj!(SPICE.SpiceIntCell(10), ckfile)
        @test_throws SPICE.SpiceError SPICE.spkcov!(SPICE.SpiceDoubleCell(10), sclk, 1)
        # getfat on the files of the tests and on the kernels of this machine.
        mk = joinpath(dir, "fat.tm")
        write(mk, "KPL/MK\n\\begindata\nKERNELS_TO_LOAD = ( '$sclk' )\n\\begintext\n")
        oneline = joinpath(dir, "oneline.txt")
        write(oneline, "\n")
        files = [pn, ckfile, sclk, text, mk, oneline, Kernels.DE440S, Kernels.DE405, Kernels.LSK,
            Kernels.PCK11, Kernels.MOON_BPC, Kernels.MOON_TF]
        seen = 0
        for f in files
            isfile(f) || continue
            seen += 1
            @test getfat(f) == SPICE.getfat(f)
        end
        @test seen >= 6
        @test getfat(text) == ("?", "?") && getfat(sclk) == ("KPL", "SCLK")
        @test throws_kind(() -> getfat(joinpath(dir, "missing.bsp")), :filenotfound)
        @test throws_kind(() -> getfat("  "), :blankfilename)
        for content in ("", "K", "DAF/SPK abc")                    # no 12 bytes and no line
            write(oneline, content)
            @test throws_kind(() -> getfat(oneline), :filereadfailed)
            @test_throws SPICE.SpiceError SPICE.getfat(oneline)
        end
        @test_throws SPICE.SpiceError SPICE.getfat(joinpath(dir, "missing.bsp"))
        @test endswith(tkvrsn("TOOLKIT"), "/" * SPICE.tkvrsn("TOOLKIT"))
        @test tkvrsn() == tkvrsn(:toolkit) && startswith(tkvrsn(), "Allspice_v")
        @test tkvrsn("SPKW13") == "No version found."
    end

    @apart @testset "furnsh of a written file; kinfo, kdata, ktotal" begin
        c2 = writer_case("spkw02")
        c3 = writer_case("spkw03")
        seg2 = merge(gen_w0203(rng, 1, 3), (body=-501, center=399))
        # The type 3 segment has the time span of the type 2 segment.
        seg3 = merge(gen_w0203(rng, 2, 6), (body=-502, center=-501, btime=seg2.btime, intlen=seg2.intlen,
            first=seg2.first, last=seg2.last, cdata=randn(rng, 24, size(seg2.cdata, 2)) .* 1.0e4))
        spk1 = joinpath(dir, "loaded_1.bsp")
        h = spkopn(spk1, "LOADED 1", 500)
        dafac(h, ["type 2 and type 3 segments for furnsh"])
        c2.neo(h, seg2)
        c3.neo(h, seg3)
        spkcls(h)
        spk2 = joinpath(dir, "loaded_2.bsp")
        h = spkopn(spk2, "LOADED 2")
        c2.neo(h, merge(seg2, (body=-503,)))
        spkcls(h)
        tk1 = joinpath(dir, "vars_1.tpc")
        write(tk1, "KPL/PCK\n\\begindata\nBODY-501_TEST = ( 1.0 2.0 )\n\\begintext\n")
        tk2 = joinpath(dir, "vars_2.tf")
        write(tk2, "KPL/FK\n\\begindata\nALLSPICE_TEST_VALUE = ( 3.0 )\n\\begintext\n")
        mk = joinpath(dir, "set.tm")
        write(mk, "KPL/MK\n\\begindata\nKERNELS_TO_LOAD = ( '$tk2',\n '$spk2' )\n\\begintext\n")
        order = [spk1, tk1, mk, spk1]
        isfile(Kernels.MOON_BPC) && push!(order, Kernels.MOON_BPC)
        ctx = neo_context(order)
        oracle_load(order)
        # The written file gives the states of CSPICE through the existing reader.
        worst = 0.0
        for (targ, obs, a) in ((-501, 399, seg2), (-502, -501, seg3), (-502, 399, seg2), (-503, 399, seg2))
            lo = max(seg2.first, seg3.first)
            hi = min(seg2.last, seg3.last)
            for et in (targ == -502 ? random_epochs(rng, 20, lo, hi) : random_epochs(rng, 20, a.first, a.last))
                s, _ = spkgeo(ctx, targ, et, "J2000", obs)
                cs, _ = SPICE.spkgeo(targ, et, "J2000", obs)
                worst = max(worst, maximum(abs.(s .- cs) ./ max.(abs.(cs), 1.0)))
            end
        end
        @test worst < 1e-11
        measured(@sprintf("  %-58s max rel=%.3e", "furnsh of a Allspice file: spkgeo vs CSPICE", worst))
        # Load-list information.
        kinds = ("ALL", "all", "SPK", "spk", "CK", "PCK", "EK", "TEXT", "META", "SPK TEXT", "text  meta ",
            "SPK PCK CK EK TEXT META", "FOO", "SPK FOO", " ")
        for kind in kinds
            n = ktotal(ctx, kind)
            @test n == SPICE.ktotal(kind)
            for which in 0:n+1
                got = kdata(ctx, which, kind)
                want = which == 0 ? nothing : SPICE.kdata(which, kind)
                @test (got === nothing) == (want === nothing)
                got === nothing || @test got[1:3] == want[1:3]
            end
        end
        @test ktotal(ctx) == ktotal(ctx, "ALL") == length(order) + 2
        @test ktotal(Allspice.kernelset(ctx), "META") == 1
        for file in [order; tk2; spk2; joinpath(dir, "not_loaded.bsp")]
            got = kinfo(ctx, file)
            want = SPICE.kinfo(file)
            @test (got === nothing) == (want === nothing)
            got === nothing || @test got[1:2] == want[1:2]
        end
        @test kinfo(ctx, spk2) == ("SPK", mk, 0)
        @test kdata(ctx, 1, "META") == (mk, "META", "", 0)
        # The same calls on the default context.
        kclear()
        @test ktotal("ALL") == 0 && kdata(1, "ALL") === nothing && kinfo(spk1) === nothing
        furnsh(spk1)
        @test ktotal("SPK") == 1 && ktotal("TEXT") == 0 && kdata(1, "SPK") == (spk1, "SPK", "", 0)
        @test kinfo(spk1) == ("SPK", "", 0) && kdata(1) == kdata(1, "ALL")
        kclear()

        # A loaded file must be unloaded before it is appended to (the documented rule).
        # A context that loaded the file before does not see the new segment.
        ctx2 = neo_context([spk1])
        before = [spkgeo(ctx2, -501, et, "J2000", 399) for et in random_epochs(MersenneTwister(1), 10, seg2.first, seg2.last)]
        unload(ctx2, spk1)
        snapshot = neo_context([spk1])
        h = spkopa(spk1)
        dafac(h, rand_comments(rng, 40; maxlen=300))               # the data move down
        c2.neo(h, merge(seg2, (body=-504,)))
        spkcls(h)
        furnsh(ctx2, spk1)
        after = [spkgeo(ctx2, -501, et, "J2000", 399) for et in random_epochs(MersenneTwister(1), 10, seg2.first, seg2.last)]
        @test after == before
        @test spkgeo(ctx2, -504, seg2.first, "J2000", 399) == spkgeo(ctx2, -501, seg2.first, "J2000", 399)
        @test spkobj(spk1) == [-504, -502, -501]
        @test throws_kind(() -> spkgeo(snapshot, -504, seg2.first, "J2000", 399), :spkinsuffdata)
        SPICE.kclear()
    end

    @apart @testset "binary PCK: Allspice reads the file that it writes" begin
        c = writer_case("pckw02")
        # One segment for each body-fixed frame, each with another base frame.
        segs = [merge(c.gen(rng, k), (clssid=id, frame=base)) for (k, (id, base)) in
                enumerate(zip(sort!(collect(keys(PCK_FRAMES))), ("ECLIPJ2000", "J2000", "B1950")))]
        pc, pn = write_both(c, dir, "pck_readback", segs)
        @test read(pc) == read(pn)
        ctx = neo_context([pn])
        oracle_load([pc])
        drot = 0.0
        dxf = 0.0
        try
            for a in segs, et in random_epochs(rng, 20, a.first, a.last)
                name = PCK_FRAMES[a.clssid]
                drot = max(drot, maximum(abs.(pxform(ctx, "J2000", name, et) .- SPICE.pxform("J2000", name, et))))
                X = sxform(ctx, "J2000", name, et)
                Xc = SPICE.sxform("J2000", name, et)
                dxf = max(dxf, maximum(abs.(X .- Xc) ./ max.(abs.(Xc), 1.0)))
            end
        finally
            SPICE.kclear()
        end
        measured(@sprintf("  %-58s rot=%.3e xform rel=%.3e", "Allspice reads its binary PCK: pxform, sxform", drot, dxf))
        # Measured: 1.1e-16 (rotation) and 1.1e-16 relative (state transformation), a
        # last-bit difference of one element.
        @test drot <= 1e-15
        @test dxf <= 1e-15
        # The base frame by ID, and by a name with blanks and in lower case, as `irfnum_`
        # compares names; a frame ID outside the inertial frames is rejected.
        a = segs[1]
        paths = String[]
        for (i, frame) in enumerate((17, " eclip J2000 ", "ECLIPJ2000"))
            p = joinpath(dir, "pck_frame_$i.bpc")
            h = pckopn(p, "IFNAME pck_readback")
            c.neo(h, merge(a, (frame=frame,)))
            pckcls(h)
            push!(paths, p)
        end
        @test allequal(read.(paths))
        h = pckopn(joinpath(dir, "pck_frame_bad.bpc"), "BAD")
        @test throws_kind(() -> c.neo(h, merge(a, (frame=22,))), :invalidrefframe)
        @test throws_kind(() -> c.neo(h, merge(a, (frame=0,))), :invalidrefframe)
        dafcls(h)
    end

    @apart @testset "handles on different tasks" begin
        c = writer_case("spkw13")
        jobs = [(seed=SEED + 900 + i, nseg=rand(rng, (3, 26, 40))) for i in 1:nsample(6, 16)]
        function job(tag, j)
            local r = MersenneTwister(j.seed)   # `local`: the test set has a variable `r` too
            p = joinpath(dir, "task_$(tag)_$(j.seed).bsp")
            h = spkopn(p, "TASK", 300)
            dafac(h, rand_comments(r, 5))
            for k in 1:j.nseg
                c.neo(h, c.gen(r, k))
            end
            dafac(h, rand_comments(r, 30; maxlen=200))
            spkcls(h)
            hr = dafopr(p)
            n = 0
            dafbfs(hr)
            while daffna(hr)
                _, ic = dafus(dafgs(hr), 2, 6)
                n += length(dafgda(hr, ic[5], ic[6]))
            end
            dafcls(hr)
            return p, n, spkobj(p), collect(spkcov(p, -1001))
        end
        serial = [job("serial", j) for j in jobs]
        tasks = [Threads.@spawn job("spawned", j) for j in jobs]
        spawned = fetch.(tasks)
        @test all(read(s[1]) == read(t[1]) for (s, t) in zip(serial, spawned))
        @test all(s[2:4] == t[2:4] for (s, t) in zip(serial, spawned))
        r = compare_files(serial[1][1], spawned[1][1])
        @test files_agree(r)
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
