# NeoSpice.jl verification suite. CSPICE (through SPICE.jl) is the oracle.
#
# Run all groups, with threads (group 7 needs more than one):
#
#     julia --startup-file=no --threads=8 --project=test test/runtests.jl
#
# Run some groups only (names are the file names without `test_` and `.jl`):
#
#     julia --startup-file=no --threads=8 --project=test test/runtests.jl states threads
#
# Kernel paths are in `test/kernels.jl`; the environment variables NEOSPICE_KERNEL_DIR,
# NEOSPICE_DE440S and NEOSPICE_DE405 override them. A group whose kernels are missing is
# skipped with a visible notice. NEOSPICE_TEST_SEED changes the random samples.
#
# Behaviour specifications (Given / When / Then) are in `test/features/*.feature`.

include("helpers.jl")

const GROUPS = ["states", "rotations", "derivatives", "identities", "pool", "management",
    "threads", "quality", "aberration", "spktypes", "time", "sclk", "ck", "dafwrite",
    "windows", "roundtrip", "sgp4"]

selected = isempty(ARGS) ? GROUPS : ARGS
unknown = setdiff(selected, GROUPS)
isempty(unknown) || error("unknown test group(s) $(unknown); the groups are $(GROUPS)")

println("NeoSpice tests: Julia ", VERSION, ", ", Threads.nthreads(), " thread(s), seed ", SEED)
Threads.nthreads() == 1 && printstyled("NOTE  one thread only: group 7 will be skipped; run " *
    "with `julia --threads=8 --project=test test/runtests.jl`\n"; color=:yellow, bold=true)

try
    @testset "NeoSpice" begin
        for g in selected
            include(joinpath(@__DIR__, "test_$g.jl"))
        end
    end
finally
    SPICE.kclear()
    kclear()
    println("\nMeasured values (for tightening the thresholds):")
    foreach(println, MEASURED)
    if !isempty(SKIPPED)
        printstyled("\nSkipped:\n"; color=:yellow, bold=true)
        for (name, reason) in SKIPPED
            println("  ", name, ": ", reason)
        end
    end
end
