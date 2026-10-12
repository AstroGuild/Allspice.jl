# Precompile directives for the common call signatures. No kernel is read here: the
# directives name argument types only, so the package image holds native code for the
# first `furnsh` and the first queries without any data file at build time.

let
    S, F, I = String, Float64, Int
    K, C = KernelSet, SpiceContext
    precompile(furnsh, (S,))
    precompile(furnsh, (C, S))
    precompile(unload, (S,))
    precompile(unload, (C, S))
    precompile(kclear, ())
    precompile(kclear, (C,))
    precompile(open_daf, (S,))
    precompile(parse_text_kernel, (S,))
    precompile(build_kernelset, (Vector{Loaded},))
    for B in (S, I), R in (S, I)
        for head in ((), (C,), (K,))
            precompile(spkezr, (head..., B, F, R, S, B))
            precompile(spkpos, (head..., B, F, R, S, B))
            for N in 0:3
                precompile(spkder, (head..., B, F, R, B, Val{N}))
                precompile(spkder, (head..., B, F, R, S, B, Val{N}))
            end
        end
    end
    for R in (S, I), head in ((), (C,), (K,))
        precompile(spkgeo, (head..., I, F, R, I))
        precompile(pxform, (head..., R, R, F))
        precompile(sxform, (head..., R, R, F))
        for N in 0:3
            precompile(pxder, (head..., R, R, F, Val{N}))
        end
    end
    for head in ((), (C,), (K,))
        precompile(utc2et, (head..., S))
        precompile(str2et, (head..., S))
        precompile(et2utc, (head..., F, S, I))
        precompile(deltet, (head..., F, S))
        precompile(timout, (head..., F, S))
        precompile(etcal, (head..., F))
        precompile(unitim, (head..., F, S, S))
        precompile(gdpool, (head..., S))
        precompile(gcpool, (head..., S))
        precompile(gipool, (head..., S))
        precompile(bodvrd, (head..., S, S))
        precompile(bodn2c, (head..., S))
        precompile(bodc2n, (head..., I))
        precompile(bods2c, (head..., S))
        precompile(namfrm, (head..., S))
        precompile(frmnam, (head..., I))
    end
end

# Geometry calls and geometry finder searches with names as strings, for a context and
# for a kernel set. A signature that does not exist gives `false` and no error.
let
    S, F, V, W = String, Float64, Vector{Float64}, Window
    for head in ((SpiceContext,), (KernelSet,))
        precompile(subpnt, (head..., S, S, F, S, S, S))
        precompile(subslr, (head..., S, S, F, S, S, S))
        precompile(sincpt, (head..., S, S, F, S, S, S, S, V))
        precompile(ilumin, (head..., S, S, F, S, S, S, V))
        precompile(illumg, (head..., S, S, S, F, S, S, S, V))
        precompile(illumf, (head..., S, S, S, F, S, S, S, V))
        precompile(phaseq, (head..., F, S, S, S, S))
        precompile(occult, (head..., S, S, S, S, S, S, S, S, F))
        precompile(limbpt, (head..., S, S, F, S, S, S, S, V, F, Int, F, F, Int))
        precompile(termpt, (head..., S, S, S, F, S, S, S, S, V, F, Int, F, F, Int))
        precompile(tangpt, (head..., S, S, F, S, S, S, S, S, V))
        precompile(fovtrg, (head..., S, S, S, S, S, S, F))
        precompile(fovray, (head..., S, V, S, S, S, F))
        precompile(getfov, (head..., Int))
        precompile(spkcpo, (head..., S, F, S, S, S, V, S, S))
        precompile(spkcpt, (head..., V, S, S, F, S, S, S, S))
        precompile(gfdist, (head..., S, S, S, S, F, F, F, W))
        precompile(gfrr, (head..., S, S, S, S, F, F, F, W))
        precompile(gfsep, (head..., S, S, S, S, S, S, S, S, S, F, F, F, W))
        precompile(gfpa, (head..., S, S, S, S, S, F, F, F, W))
        precompile(gfposc, (head..., S, S, S, S, S, S, S, F, F, F, W))
        precompile(gfsubc, (head..., S, S, S, S, S, S, S, S, F, F, F, W))
        precompile(gfsntc, (head..., S, S, S, S, S, S, V, S, S, S, F, F, F, W))
        precompile(gfilum, (head..., S, S, S, S, S, S, S, V, S, F, F, F, W))
        precompile(gfoclt, (head..., S, S, S, S, S, S, S, S, S, F, W))
        precompile(gftfov, (head..., S, S, S, S, S, S, F, W))
        precompile(gfrfov, (head..., S, V, S, S, S, F, W))
    end
end

# EK: open, read, query and the record writers.
let
    S, I, H = String, Int, DASHandle
    precompile(open_ek, (S,))
    precompile(ekopr, (S,))
    precompile(eknseg, (H,))
    precompile(ekssum, (H, I))
    for f in (ekrcec, ekrced, ekrcei)
        precompile(f, (H, I, I, S))
    end
    for head in ((SpiceContext,), (KernelSet, TimeDefaults))
        precompile(ekfind, (head..., S))
        precompile(ekpsel, (head..., S))
    end
    for f in (ekgc, ekgd, ekgi)
        precompile(f, (EKResult, I, I))
        precompile(f, (EKResult, I, I, I))
    end
    precompile(eknelt, (EKResult, I, I))
    precompile(ekopn, (S, S, I))
    precompile(ekopw, (S,))
    precompile(ekbseg, (H, S, Vector{S}, Vector{S}))
    precompile(ekappr, (H, I))
    precompile(ekinsr, (H, I, I))
    precompile(ekdelr, (H, I, I))
    precompile(ekifld, (H, S, I, Vector{S}, Vector{S}))
    precompile(ekffld, (H, I, EKFastLoad))
    for (f, T) in ((ekacec, S), (ekaced, Float64), (ekacei, I), (ekucec, S), (ekuced, Float64), (ekucei, I))
        precompile(f, (H, I, I, S, Vector{T}, Bool))
    end
    for (f, T) in ((ekaclc, S), (ekacld, Float64), (ekacli, I))
        precompile(f, (H, I, S, Vector{T}, Vector{I}, Vector{Bool}, EKFastLoad))
    end
end
