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
