# Rotation benchmark: Allspice against CSPICE (SPICE.jl).
#
#     julia --startup-file=no --project=bench bench/rotations.jl
#
# `pxform` from J2000 to IAU_EARTH (text PCK: polynomials and trigonometric terms) and to
# MOON_PA (fixed-offset frame on a binary PCK Chebyshev segment), `sxform`, and `pxder`
# for orders 0 to 3, at the 100 epochs of the state benchmark.

include("common.jl")
using Allspice
import SPICE

require_kernels(Kernels.TEXT_SET_FILES..., Kernels.MOON_BPC)

function loop_neo_pxform(ctx, epochs, from, to)
    acc = 0.0
    for et in epochs
        acc += pxform(ctx, from, to, et)[1]
    end
    return acc
end
function loop_neo_sxform(ctx, epochs, from, to)
    acc = 0.0
    for et in epochs
        acc += sxform(ctx, from, to, et)[1]
    end
    return acc
end
function loop_neo_pxder(ctx, epochs, from, to, v::Val{N}) where {N}
    acc = 0.0
    for et in epochs
        acc += pxder(ctx, from, to, et, v)[N+1][1]
    end
    return acc
end
function loop_cspice_pxform(epochs, from, to)
    acc = 0.0
    for et in epochs
        acc += SPICE.pxform(from, to, et)[1]
    end
    return acc
end
function loop_cspice_sxform(epochs, from, to)
    acc = 0.0
    for et in epochs
        acc += SPICE.sxform(from, to, et)[1]
    end
    return acc
end

function main()
    epochs = anise_epochs()
    n = length(epochs)
    set = [Kernels.text_set()..., Kernels.MOON_BPC]
    SPICE.kclear(); foreach(SPICE.furnsh, set)
    ctx = SpiceContext(); foreach(p -> furnsh(ctx, p), set)

    rows = Row[]
    agreement = String[]
    for to in ("IAU_EARTH", "MOON_PA")
        toid = Int(namfrm(ctx, to))
        @assert toid != 0
        # The engines must compute the same thing before they are timed.
        d = maximum(maximum(abs.(pxform(ctx, 1, toid, et) .- SPICE.pxform("J2000", to, et)))
                    for et in epochs)
        @assert d <= 1e-8 "Allspice differs from CSPICE for J2000 to $to: $d"
        push!(agreement, @sprintf("| J2000 to %s | %.2e |", to, d))
        g = "J2000 to $to"
        push!(rows, bench_row(g, "Allspice pxform, integer IDs", loop_neo_pxform, ctx, epochs, 1, toid; n))
        push!(rows, bench_row(g, "Allspice pxform, names", loop_neo_pxform, ctx, epochs, "J2000", to; n))
        push!(rows, bench_row(g, "Allspice sxform, integer IDs", loop_neo_sxform, ctx, epochs, 1, toid; n))
        for N in 0:3
            push!(rows, bench_row(g, "Allspice pxder Val($N)", loop_neo_pxder, ctx, epochs, 1, toid, Val(N); n))
        end
        push!(rows, bench_row(g, "CSPICE pxform (SPICE.jl)", loop_cspice_pxform, epochs, "J2000", to; n))
        push!(rows, bench_row(g, "CSPICE sxform (SPICE.jl)", loop_cspice_sxform, epochs, "J2000", to; n))
    end

    text = run_header("Rotation benchmark: Allspice, CSPICE") * """
    Kernels: naif0012.tls, pck00011.tpc, moon_de440_200625.tf, the MOON_PA alias kernel of
    the tests and moon_pa_de440_200625.bpc. Epochs: 100, evenly spaced from 1900-01-01
    12:00 ET to 2099-01-01 12:00 ET (end excluded).

    Largest element difference of `pxform` from CSPICE on the timed queries:

    | Rotation | max abs difference |
    |---|---:|
    $(join(agreement, "\n"))

    ## Time per query

    """ * markdown_table(rows)
    println()
    println(text)
    write_results("rotations", text)
    SPICE.kclear()
    return rows
end

main()
