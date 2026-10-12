# Group 8: quality gates. Zero allocation and inferred return types for the integer-ID
# forms on an explicit context, after a warm-up call.
#
# Each measurement is inside a function, so that the arguments are not global variables
# (a call on untyped globals allocates on its own account).

# Compile time. In the default run the test code of this group is compiled without
# optimization (`@optlevel 0` for `Main`, until the end of this file), because the
# compilation of the test code is the larger part of the run time of the group. The
# methods of NeoSpice keep their own optimization level, except where the compiler
# inlines one into test code. The exhaustive run (`NEOSPICE_TEST_FULL=1`) is not changed.
FULL || @eval Base.Experimental.@optlevel 0

alloc_spkezr(ctx, et) = @allocated spkezr(ctx, 301, et, 1, "NONE", 399)
alloc_spkpos(ctx, et) = @allocated spkpos(ctx, 301, et, 1, "NONE", 399)
alloc_spkgeo(ctx, et) = @allocated spkgeo(ctx, 301, et, 1, 399)
alloc_spkezr_chain(ctx, et) = @allocated spkezr(ctx, 301, et, 1, "NONE", 10)
alloc_spkezr_rot(ctx, et, ref) = @allocated spkezr(ctx, 301, et, ref, "NONE", 399)
alloc_spkder(ctx, et, ::Val{N}) where {N} = @allocated spkder(ctx, 301, et, 1, 399, Val(N))
alloc_spkder_rot(ctx, et, ref, ::Val{N}) where {N} = @allocated spkder(ctx, 301, et, ref, 399, Val(N))
alloc_pxform(ctx, from, to, et) = @allocated pxform(ctx, from, to, et)
alloc_sxform(ctx, from, to, et) = @allocated sxform(ctx, from, to, et)
alloc_pxder(ctx, from, to, et, ::Val{N}) where {N} = @allocated pxder(ctx, from, to, et, Val(N))
alloc_names(ctx, et) = @allocated spkezr(ctx, "MOON", et, "J2000", "NONE", "EARTH")
alloc_pxform_names(ctx, et) = @allocated pxform(ctx, "J2000", "IAU_EARTH", et)

"Allocation of the second call of `f(args...)` (the first call compiles)."
function warm(f, args...)
    f(args...)
    return f(args...)
end

group("8. quality gates", Kernels.FULL_SET_FILES...) do
    ctx = neo_context(Kernels.full_set())
    et = 1.2345e8
    iau_earth = Int(namfrm(ctx, "IAU_EARTH"))
    moon_pa = Int(namfrm(ctx, "MOON_PA"))
    moon_me = Int(namfrm(ctx, "MOON_ME"))
    eclip = 17
    @test iau_earth == 10013 && moon_pa == 31000 && moon_me == 31001

    @testset "zero allocation, integer-ID forms" begin
        @test warm(alloc_spkezr, ctx, et) == 0
        @test warm(alloc_spkpos, ctx, et) == 0
        @test warm(alloc_spkgeo, ctx, et) == 0
        @test warm(alloc_spkezr_chain, ctx, et) == 0
        for ref in (eclip, iau_earth, moon_pa, moon_me)
            @test warm(alloc_spkezr_rot, ctx, et, ref) == 0
        end
        for N in 0:3
            @test warm(alloc_spkder, ctx, et, Val(N)) == 0
            @test warm(alloc_spkder_rot, ctx, et, iau_earth, Val(N)) == 0
            @test warm(alloc_spkder_rot, ctx, et, moon_pa, Val(N)) == 0
        end
        for (from, to) in ((1, eclip), (1, iau_earth), (iau_earth, 1), (1, moon_pa),
                (1, moon_me), (iau_earth, moon_pa))
            @test warm(alloc_pxform, ctx, from, to, et) == 0
            @test warm(alloc_sxform, ctx, from, to, et) == 0
            for N in 0:3
                @test warm(alloc_pxder, ctx, from, to, et, Val(N)) == 0
            end
        end
    end

    @testset "allocation of the name forms (reported, not asserted)" begin
        measured("  name form spkezr(ctx, \"MOON\", et, \"J2000\", \"NONE\", \"EARTH\"): " *
            "$(warm(alloc_names, ctx, et)) bytes")
        measured("  name form pxform(ctx, \"J2000\", \"IAU_EARTH\", et): " *
            "$(warm(alloc_pxform_names, ctx, et)) bytes")
        @test true
    end

    @testset "inferred return types" begin
        @test (@inferred spkezr(ctx, 301, et, 1, "NONE", 399)) isa Tuple{SVector{6,Float64},Float64}
        @test (@inferred spkpos(ctx, 301, et, 1, "NONE", 399)) isa Tuple{SVector{3,Float64},Float64}
        @test (@inferred spkgeo(ctx, 301, et, 1, 399)) isa Tuple{SVector{6,Float64},Float64}
        @test (@inferred spkezr(ctx, 301, et, iau_earth, "NONE", 399)) isa Tuple{SVector{6,Float64},Float64}
        @test (@inferred spkezr(ctx, "MOON", et, "J2000", "NONE", "EARTH")) isa Tuple{SVector{6,Float64},Float64}
        @test (@inferred spkder(ctx, 301, et, 1, 399, Val(0))) isa NTuple{1,SVector{3,Float64}}
        @test (@inferred spkder(ctx, 301, et, 1, 399, Val(1))) isa NTuple{2,SVector{3,Float64}}
        @test (@inferred spkder(ctx, 301, et, 1, 399, Val(2))) isa NTuple{3,SVector{3,Float64}}
        @test (@inferred spkder(ctx, 301, et, 1, 399, Val(3))) isa NTuple{4,SVector{3,Float64}}
        @test (@inferred spkder(ctx, 301, et, moon_pa, 399, Val(3))) isa NTuple{4,SVector{3,Float64}}
        @test (@inferred pxform(ctx, 1, iau_earth, et)) isa SMatrix{3,3,Float64,9}
        @test (@inferred pxform(ctx, 1, moon_me, et)) isa SMatrix{3,3,Float64,9}
        @test (@inferred pxform(ctx, "J2000", "IAU_EARTH", et)) isa SMatrix{3,3,Float64,9}
        @test (@inferred sxform(ctx, 1, iau_earth, et)) isa SMatrix{6,6,Float64,36}
        @test (@inferred pxder(ctx, 1, iau_earth, et, Val(0))) isa NTuple{1,SMatrix{3,3,Float64,9}}
        @test (@inferred pxder(ctx, 1, iau_earth, et, Val(1))) isa NTuple{2,SMatrix{3,3,Float64,9}}
        @test (@inferred pxder(ctx, 1, iau_earth, et, Val(2))) isa NTuple{3,SMatrix{3,3,Float64,9}}
        @test (@inferred pxder(ctx, 1, moon_pa, et, Val(3))) isa NTuple{4,SMatrix{3,3,Float64,9}}
        # An integer epoch and a Float32 epoch give Float64 results.
        @test (@inferred spkezr(ctx, 301, 0, 1, "NONE", 399)) isa Tuple{SVector{6,Float64},Float64}
    end

    @testset "package image and method table" begin
        # A method defined twice makes precompilation fail, and the package then loads
        # without its image (slow first calls) although every other test passes.
        project = dirname(Base.active_project())
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$project -e "using NeoSpice"`
        err = IOBuffer()
        ok = success(pipeline(cmd; stderr=err, stdout=devnull))
        text = String(take!(err))
        @test ok
        @test !occursin("overwritten", text)
        @test !occursin("precompilation", text)
        @test isempty(Test.detect_ambiguities(NeoSpice))
    end
end

FULL || @eval Base.Experimental.@optlevel -1    # the level of the session again
