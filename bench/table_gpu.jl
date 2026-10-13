# Ephemeris tables on a GPU: build two tables on the CPU, move them to the GPU and
# evaluate them at many epochs, by broadcast and in a hand-written kernel.
#
# CUDA.jl is not a dependency of the package or of `bench/Project.toml`. Run with a
# temporary environment:
#
#     julia --startup-file=no -e 'using Pkg; Pkg.activate(temp=true);
#         Pkg.develop(path="."); Pkg.add(["CUDA", "Adapt", "StaticArrays", "ForwardDiff"]);
#         include("bench/table_gpu.jl")'

using Allspice, CUDA, Adapt, StaticArrays, LinearAlgebra, Printf
import ForwardDiff

include(joinpath(@__DIR__, "..", "test", "kernels.jl"))

position_gpu(tab, et) = @inbounds table_position(tab, et)
state_gpu(tab, et) = @inbounds table_state(tab, et)
speed_gpu(tab, et) = ForwardDiff.derivative(t -> norm(position_gpu(tab, t)), et)

"Kernel: distance (km) between the positions of two tables at each epoch."
function distance_kernel!(out, a, b, ets)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if i <= length(ets)
        @inbounds out[i] = norm(table_position(a, ets[i]) - table_position(b, ets[i]))
    end
    return nothing
end

function main()
    ctx = SpiceContext()
    foreach(p -> furnsh(ctx, p), Kernels.full_set())
    t0, t1 = 7.0e8, 7.0e8 + 400 * 86400.0
    earth = ephemeris_table(ctx, 399, t0, t1, "J2000", "NONE", 10)
    moon = ephemeris_table(ctx, 301, t0, t1, "J2000", "NONE", 10)
    n = 4_000_000
    ets = collect(range(t0, t1, length=n))
    gets = CuArray(ets)
    println(CUDA.name(CUDA.device()), ", ", n, " epochs, table of ", earth.nint, " intervals, ",
        earth.ncoef, " coefficients")

    rows = Tuple{String,Float64,Float64}[]
    for (name, T) in (("Float64", Float64), ("Float32", Float32))
        cpu = T === Float64 ? earth : Allspice.EphemerisTable{T}(earth)
        gpu = adapt(CuArray, cpu)
        @assert isbits(cudaconvert(gpu))
        f = let gpu = gpu; et -> state_gpu(gpu, et); end
        out = Array(map(f, gets))                      # compile
        best = minimum(@elapsed(CUDA.@sync map(f, gets)) for _ in 1:5)
        exact = [spkezr(ctx, 399, et, "J2000", "NONE", 10)[1] for et in ets[1:997:end]]
        err = maximum(norm(Float64.(out[1+997*(k-1)][SVector(1, 2, 3)]) - exact[k][SVector(1, 2, 3)]) for k in eachindex(exact))
        push!(rows, ("GPU broadcast, table_state, " * name, best / n * 1e9, err))
    end

    # CPU reference, one thread
    cpu_table = minimum(@elapsed(map(et -> state_gpu(earth, et), ets)) for _ in 1:3)
    cpu_spk = minimum(@elapsed(map(et -> spkezr(ctx, 399, et, "J2000", "NONE", 10)[1], ets)) for _ in 1:3)
    push!(rows, ("CPU one thread, table_state", cpu_table / n * 1e9, earth.maxerr))
    push!(rows, ("CPU one thread, spkezr", cpu_spk / n * 1e9, 0.0))

    # Hand-written kernel with two tables
    ge, gm = adapt(CuArray, earth), adapt(CuArray, moon)
    gout = CUDA.zeros(Float64, n)
    kernel = @cuda launch=false distance_kernel!(gout, ge, gm, gets)
    threads = 256
    kernel(gout, ge, gm, gets; threads, blocks=cld(n, threads))
    best = minimum(@elapsed(CUDA.@sync kernel(gout, ge, gm, gets; threads, blocks=cld(n, threads))) for _ in 1:5)
    d = Array(gout)
    err = maximum(abs(d[i] - norm(spkpos(ctx, 301, ets[i], "J2000", "NONE", 399)[1])) for i in 1:997:n)
    push!(rows, ("GPU kernel, Earth-Moon distance, 2 tables", best / n * 1e9, err))

    # Dual-number epoch on the GPU: range rate of the Earth relative to the Sun
    g = let ge = ge; et -> speed_gpu(ge, et); end
    rr = Array(map(g, gets))
    best = minimum(@elapsed(CUDA.@sync map(g, gets)) for _ in 1:5)
    err = maximum(1:997:n) do i
        x = spkezr(ctx, 399, ets[i], "J2000", "NONE", 10)[1]
        abs(rr[i] - dot(x[SVector(1, 2, 3)], x[SVector(4, 5, 6)]) / norm(x[SVector(1, 2, 3)]))
    end
    push!(rows, ("GPU broadcast, range rate by ForwardDiff", best / n * 1e9, err))

    println("\n| Case | ns per epoch | largest difference from the kernels |")
    println("|---|---:|---:|")
    for (name, ns, err) in rows
        @printf("| %s | %.3f | %.2e |\n", name, ns, err)
    end
end

main()
