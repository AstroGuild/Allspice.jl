# Contexts and threads

All loaded kernels belong to a [`SpiceContext`](@ref), and this page calls the loaded kernels of one context its kernel set. Every call that reads the loaded kernels, the kernel pool or a setting of a context exists in two forms: with a context as the first argument, and without one, in which case the call uses the context that [`default_context`](@ref) returns, which is one object for the whole process and is as safe to share as any other. A call that depends on its arguments alone, such as `vnorm`, has no context argument. A program that loads one set of kernels never needs to name a context.

```@setup ctx
using Allspice
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
kernels = DocKernels.write(mktempdir())
```

## Several kernel sets

A second context is an independent kernel set, with its own kernel pool (the variables that text kernels define) and its own settings, such as the defaults of [`timdef`](@ref). This is how a program compares two versions of an ephemeris, or keeps the kernels of two missions apart.

```@example ctx
using Allspice

mission = SpiceContext()
foreach(k -> furnsh(mission, k), kernels)

(ktotal(mission, "ALL"), ktotal("ALL"))     # the default context is still empty
```

## Queries from many threads

A query reads one snapshot of the kernel set, which is an immutable version that a later load or unload does not change, and takes no lock, so any number of threads can call `spkezr`, `pxform` and the other queries on the same context at once. The loop below computes a distance at 10,000 epochs on all the threads that Julia has.

```@example ctx
et = str2et(mission, "2025 MAR 01 12:00:00 UTC")
epochs = et .+ range(0, 86400; length=10_000)
distance = Vector{Float64}(undef, length(epochs))

Threads.@threads for i in eachindex(epochs)
    position, _ = spkpos(mission, "PROBE", epochs[i], "J2000", "LT+S", "EARTH")
    distance[i] = vnorm(position)
end

extrema(distance)                           # km
```

A thread can also load or unload a kernel while other threads query. [`furnsh`](@ref), [`unload`](@ref) and the calls that change the kernel pool build a new snapshot and replace the old one in one step, so that a reader sees the set before the change or after it and never a mixture. Two changes at the same time wait for each other. Two consecutive queries of one thread can therefore see different kernel sets if another thread changes the set between them.

## What must not be shared

A file that is open for direct access (a handle from `dafopr`, `dasopr`, `spkopn` and similar calls) belongs to the code that opened it, and one handle must not be used from two threads at once. The same holds for a [`Window`](@ref) that a thread changes. An [`EphemerisTable`](@ref) or a [`RotationTable`](@ref) is immutable and can be shared.
