# Group: static proof that the numeric kernels do not allocate (AllocCheck).
#
# The other groups measure the allocation of one call (`@allocated`), which covers the
# path that the call took. AllocCheck examines the compiled code, so a finding of zero
# covers every path that does not throw. The calls that take names of bodies and frames
# (`spkezr`, `pxform`, ...) are not in the list: they have paths that allocate by design
# (the load of a dynamic frame, an error message), and the groups `states`, `rotations`
# and `dynframes` measure them on the paths that matter.

import AllocCheck

group("allocation proof: numeric kernels (AllocCheck)") do
    V = SVector{3,Float64}; V4 = SVector{4,Float64}; V6 = SVector{6,Float64}
    M = SMatrix{3,3,Float64,9}; F = Float64
    ctx = neo_context([Kernels.LSK, Kernels.DE440S])
    tab = ephemeris_table(ctx, "EARTH", 0.0, 86400.0 * 30, "J2000", "NONE", "SUN")
    cases = Any[
        (vhat, (V,)), (vnorm, (V,)), (unorm, (V,)), (vcrss, (V, V)), (vdot, (V, V)), (vsep, (V, V)),
        (vproj, (V, V)), (vperp, (V, V)), (ucrss, (V, V)), (vrotv, (V, V, F)), (vlcom, (F, V, F, V)),
        (mxv, (M, V)), (mtxv, (M, V)), (mxm, (M, M)), (mtxm, (M, M)), (mxmt, (M, M)), (xpose, (M,)),
        (invert, (M,)), (Allspice.det, (M,)), (rotate, (F, Int)), (rotmat, (M, F, Int)),
        (twovec, (V, Int, V, Int)), (mxvg, (M, V)), (mxmg, (M, M)), (mxmtg, (M, M)),
        (eul2m, (F, F, F, Int, Int, Int)), (m2eul, (M, Int, Int, Int)), (m2q, (M,)), (q2m, (V4,)),
        (qxq, (V4, V4)), (axisar, (V, F)), (raxisa, (M,)),
        (reclat, (V,)), (latrec, (F, F, F)), (recsph, (V,)), (sphrec, (F, F, F)), (reccyl, (V,)),
        (cylrec, (F, F, F)), (recrad, (V,)), (radrec, (F, F, F)), (recgeo, (V, F, F)),
        (georec, (F, F, F, F, F)),
        (nearpt, (V, F, F, F)), (surfpt, (V, V, F, F, F)), (surfnm, (F, F, F, V)),
        (npedln, (F, F, F, V, V)), (dvhat, (V6,)), (dvdot, (V6, V6)), (dvcrss, (V6, V6)),
        (dvsep, (V6, V6)), (ducrss, (V6, V6)),
        (prop2b, (F, V6, F)), (conics, (SVector{8,Float64}, F)), (oscelt, (V6, F, F)),
        (stelab, (V, V)), (deltet, (SpiceContext, F, String)), (unitim, (SpiceContext, F, String, String)),
        (table_state, (typeof(tab), F)), (table_position, (typeof(tab), F)),
        (table_jet, (typeof(tab), F, Val{2})),
        (Allspice.spk_eval, (Vector{Float64}, Allspice.Segment, F, Val{1})),
    ]
    bad = String[]
    for (f, tt) in cases        # invariant: `bad` names the functions before this one with a finding
        isempty(AllocCheck.check_allocs(f, tt)) || push!(bad, string(f))
    end
    isempty(bad) || @info "AllocCheck findings" functions = bad
    # `spk_eval` reaches the MVector buffers of the interpolating SPK types.
    @test isempty(setdiff(bad, ["spk_eval"]))
    @test "spk_eval" ∉ bad skip = STACK_BUFFERS_ESCAPE
    # The analysis does report a function that allocates.
    @test !isempty(AllocCheck.check_allocs(x -> [x, x], (Float64,)))
    measured("alloccheck: $(length(cases)) numeric kernels without a possible allocation")
end
