# Event searches

The geometry finder answers questions of the form "when is this quantity equal to, above or below a value, or at an extremum?". Each search takes a confinement window, which is the set of times to search, and returns the window of times at which the condition holds.

```@setup gf
using Allspice
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Windows

A [`Window`](@ref) is a sorted list of disjoint intervals of time. [`wninsd!`](@ref) adds an interval, [`wncard`](@ref) counts the intervals, and [`wnfetd`](@ref) returns one of them, with indices from 1.

```@example gf
using Allspice

start = str2et("2025 MAR 18 00:00:00 TDB")
search = Window()
wninsd!(search, start, start + 2 * 86400.0)
wncard(search)
```

## Eclipses

[`gfoclt`](@ref) finds the times at which one body hides another from an observer. The spacecraft of the examples is in the shadow of the Earth when the Earth occults the Sun as seen from the spacecraft.

```@example gf
eclipses = gfoclt("ANY", "EARTH", "ELLIPSOID", "IAU_EARTH",
    "SUN", "ELLIPSOID", "IAU_SUN", "LT", "PROBE", 60.0, search)

for i in 1:wncard(eclipses)
    t0, t1 = wnfetd(eclipses, i)
    println(timout(t0, "YYYY Mon DD HR:MN:SC ::TDB"), "   ", round(t1 - t0; digits=1), " s")
end
```

## A condition on a coordinate

[`gfposc`](@ref) searches on one coordinate of a position. The search below finds when the spacecraft is north of 20 degrees of latitude, and [`wnsumd`](@ref) gives the total duration.

```@example gf
north = gfposc("PROBE", "IAU_EARTH", "NONE", "EARTH", "LATITUDINAL", "LATITUDE",
    ">", 20 * rpd(), 0.0, 600.0, search)
wnsumd(north)[1] / 3600                     # hours
```

## A function of your own

[`gfuds`](@ref) searches on any scalar function of time. It needs the function and a second function that says whether the first one decreases at an epoch; [`uddc`](@ref) builds that answer from a numerical derivative.

```@example gf
speed(t) = vnorm(spkezr("PROBE", t, "J2000", "NONE", "EARTH")[1][4:6])
fast = gfuds(speed, t -> uddc(speed, t, 10.0), ">", 8.0, 0.0, 600.0, search)
wncard(fast)                                # intervals with a speed above 8 km/s
```

## The step

Every search samples the condition at the step that the caller gives and then refines each change that it finds. The step must be shorter than the shortest interval on which the condition holds and shorter than the shortest gap between two such intervals; an event that begins and ends between two samples is not found. A long step is fast and a short step is safe, and the right value follows from the dynamics. The first eclipse above lasts less than two minutes, so the search uses a step of 60 seconds; with a step of 300 seconds the samples can fall on both sides of that eclipse, and the search then does not report it.

The relations are `"="`, `"<"`, `">"`, `"LOCMIN"`, `"LOCMAX"`, `"ABSMIN"` and `"ABSMAX"`. The convergence tolerance is one microsecond unless the keyword `tol` changes it.
