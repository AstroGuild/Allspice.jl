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

## The step

Every search samples its quantity at the step that the caller gives and then refines each change that it finds, to the tolerance in seconds of the keyword `tol`; without `tol`, the tolerance is the one that [`gfstol`](@ref) sets for the context, by default one microsecond. A long step is fast, and a step that is too long loses events:

| Search | The step must be shorter than |
|---|---|
| A condition that holds or not (`"="`, `"<"`, `">"`, an occultation, a field of view) | the shortest interval on which the condition holds, and the shortest gap between two such intervals |
| An extremum (`"LOCMIN"`, `"LOCMAX"`, `"ABSMIN"`, `"ABSMAX"`) | the shortest interval on which the quantity only increases or only decreases |

An event that begins and ends between two samples is not found, and the search gives no warning.

## Eclipses

[`gfoclt`](@ref) finds the times at which one body hides another from an observer. The spacecraft of the examples is in the shadow of the Earth when the Earth occults the Sun as seen from the spacecraft. The arguments are the kind of occultation (`"ANY"` accepts partial and full), then the name, the shape and the body-fixed frame of the front body and of the back body, the aberration correction, the observer, the step and the window. The first eclipse below lasts less than two minutes, so the step is 60 seconds; with 300 seconds the samples can fall on both sides of it, and the search then does not report it.

```@example gf
eclipses = gfoclt("ANY", "EARTH", "ELLIPSOID", "IAU_EARTH",
    "SUN", "ELLIPSOID", "IAU_SUN", "LT", "PROBE", 60.0, search)

for i in 1:wncard(eclipses)
    t0, t1 = wnfetd(eclipses, i)
    println(timout(t0, "YYYY Mon DD HR:MN:SC ::TDB"), "   ", round(t1 - t0; digits=1), " s")
end
```

## A condition on a coordinate

[`gfposc`](@ref) searches on one coordinate of a position. The search below finds when the spacecraft is north of 20 degrees of latitude; `rpd()` is the number of radians in a degree. After the relation come the reference value, an adjustment that only the absolute extrema use, and the step. The first result of [`wnsumd`](@ref) is the total duration.

```@example gf
north = gfposc("PROBE", "IAU_EARTH", "NONE", "EARTH", "LATITUDINAL", "LATITUDE",
    ">", 20 * rpd(), 0.0, 600.0, search)
wnsumd(north)[1] / 3600                     # hours
```

## A function of your own

[`gfuds`](@ref) searches on any scalar function of time. It needs the function and a second function that says whether the first one decreases at an epoch; [`uddc`](@ref) builds that answer from a numerical derivative over the interval that its third argument gives, here 10 seconds.

```@example gf
speed(t) = vnorm(spkezr("PROBE", t, "J2000", "NONE", "EARTH")[1][4:6])
fast = gfuds(speed, t -> uddc(speed, t, 10.0), ">", 8.0, 0.0, 600.0, search)
wncard(fast)                                # intervals with a speed above 8 km/s
```

## Several threads

Every search accepts the keyword `threads=true`, which searches the intervals of a window on separate threads when Julia runs with several threads. The window is the confinement window, and for a condition on a value (`"="`, `"<"`, `">"`) also the intervals on which the quantity only increases or only decreases, which the search finds first. A search for a local extremum (`"LOCMIN"`, `"LOCMAX"`), an occultation or a field of view on a window of one interval therefore gains nothing. A function given to `gfuds` or `gfudb` is then called from several threads at once and must be safe for that.
