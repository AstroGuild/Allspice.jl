# Allspice.jl verification suite. CSPICE (through SPICE.jl) is the oracle.
#
# Run all groups, with threads (group 7 needs more than one):
#
#     julia --startup-file=no --threads=8 --project=test test/runtests.jl
#
# Run some groups only (names are the file names without `test_` and `.jl`):
#
#     julia --startup-file=no --threads=8 --project=test test/runtests.jl states threads
#
# Run the groups in several processes, which is the fast way to run the whole suite
# (about 35 s against 160 s on the development machine):
#
#     julia --startup-file=no --threads=8 test/parallel.jl
#
# ALLSPICE_TEST_FULL=1 selects the exhaustive run (large samples, all derivative orders).
#
# Kernel paths are in `test/kernels.jl`; the environment variables ALLSPICE_KERNEL_DIR,
# ALLSPICE_DE440S and ALLSPICE_DE405 override them. A group whose kernels are missing is
# skipped with a visible notice. ALLSPICE_TEST_SEED changes the random samples.
#
# Behaviour specifications (Given / When / Then) are in `test/features/*.feature`.

include("helpers.jl")

const GROUPS = ["states", "rotations", "derivatives", "identities", "pool", "management",
    "threads", "quality", "aberration", "spktypes", "time", "sclk", "ck", "dafwrite",
    "windows", "roundtrip", "sgp4", "math_vectors", "math_rotations", "math_coordinates",
    "math_ellipsoid", "math_conics", "dynframes", "poolwrite", "frameinfo", "das", "dsk",
    "gf_solver", "gf_quantities", "gf_posc", "gf_threads", "gf_surface", "gf_illum", "gf_occult",
    "gf_fov",
    "geom_subpoints", "geom_illum", "geom_occult", "geom_limbterm", "geom_fov", "geom_constobs",
    "geom_azel", "table",
    "ek_file", "ek_parse", "ek_query", "ek_threads", "ek_write"]

"Wall time of each group of this run, in the order of the run."
const GROUP_SECONDS = Pair{String,Float64}[]

selected = isempty(ARGS) ? GROUPS : ARGS
unknown = setdiff(selected, GROUPS)
isempty(unknown) || error("unknown test group(s) $(unknown); the groups are $(GROUPS)")

println("Allspice tests: Julia ", VERSION, ", ", Threads.nthreads(), " thread(s), seed ", SEED)
Threads.nthreads() == 1 && printstyled("NOTE  one thread only: group 7 will be skipped; run " *
    "with `julia --threads=8 --project=test test/runtests.jl`\n"; color=:yellow, bold=true)

try
    @testset "Allspice" begin
        for g in selected
            push!(GROUP_SECONDS, g => @elapsed include(joinpath(@__DIR__, "test_$g.jl")))
        end
    end
finally
    SPICE.kclear()
    kclear()
    println("\nMeasured values (for tightening the thresholds):")
    foreach(println, MEASURED)
    println("\nTime of each group (s), longest first; the whole suite has a budget of 120 s:")
    for (g, t) in sort(GROUP_SECONDS; by=last, rev=true)
        println("  ", rpad(g, 18), lpad(round(t; digits=1), 7))
    end
    println("  ", rpad("all", 18), lpad(round(sum(last, GROUP_SECONDS; init=0.0); digits=1), 7))
    if !isempty(SKIPPED)
        printstyled("\nSkipped:\n"; color=:yellow, bold=true)
        for (name, reason) in SKIPPED
            println("  ", name, ": ", reason)
        end
    end
end
