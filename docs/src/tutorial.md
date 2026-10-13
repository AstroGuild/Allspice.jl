# First session

This page goes through one complete session: it loads kernels, converts a time, asks for a state and a rotation, finds the point of the Earth below a spacecraft, and searches for the times of perigee. Each call has the name that it has in CSPICE, without the `_c` suffix.

```@setup first
using Allspice
include(joinpath(@__DIR__, "..", "fixture.jl"))
kernels = DocKernels.write(mktempdir())
```

## Load kernels

A kernel is a file of data: an SPK holds trajectories, a text kernel holds constants such as leap seconds and body radii. [`furnsh`](@ref) loads a kernel of any kind, and [`ktotal`](@ref) counts the loaded kernels of one kind or of all kinds. In the block below, `kernels` is a vector of three paths.

```@example first
using Allspice

foreach(furnsh, kernels)      # leap seconds, constants of the Earth, an SPK
ktotal("ALL")
```

A meta-kernel, which is a text kernel that lists other kernels, works as in CSPICE, so that one `furnsh` call can load a whole set.

## Convert a time

All calls take an epoch as seconds past J2000 in Barycentric Dynamical Time (TDB), which SPICE calls ET. [`str2et`](@ref) converts a calendar string to that number. [`timout`](@ref) converts back with a picture, which is a string that shows the wanted format; the marker `::TDB` asks for the TDB time scale.

```@example first
et = str2et("2025 MAR 01 12:00:00 UTC")
```

```@example first
timout(et, "YYYY Mon DD HR:MN:SC.### ::TDB")
```

The input string and the output string differ by 69.185 seconds, which is the difference between TDB and UTC at that date.

## State of a body

[`spkezr`](@ref) returns the state of a target relative to an observer, as a position in kilometres and a velocity in kilometres per second, together with the one-way light time in seconds. The frame and the aberration correction are arguments: `"NONE"` gives the geometric state, and `"LT+S"` gives the state that the observer sees, corrected for light time and stellar aberration.

```@example first
state, lt = spkezr("PROBE", et, "J2000", "NONE", "EARTH")
state
```

The result is a static vector of six numbers, so the call allocates nothing. [`spkpos`](@ref) returns the position alone:

```@example first
sun, lt_sun = spkpos("SUN", et, "J2000", "LT+S", "PROBE")
lt_sun                        # seconds from the Sun to the spacecraft
```

## Rotation between frames

[`pxform`](@ref) returns the matrix that rotates a vector from one frame to another at an epoch, and [`sxform`](@ref) returns the 6 by 6 matrix that also transforms the velocity. `IAU_EARTH` is the frame that rotates with the Earth.

```@example first
r = pxform("J2000", "IAU_EARTH", et)
r * state[1:3]                # the position in the rotating frame, km
```

## The point below the spacecraft

[`subpnt`](@ref) returns the point of a body that is nearest to an observer, in the body-fixed frame; the method string selects the nearest point on the ellipsoid whose radii the constants kernel gives. The call also returns the epoch at the body and the vector from the observer to the point. [`reclat`](@ref) converts it to radius, longitude and latitude in radians, and `dpr()` is the number of degrees in a radian.

```@example first
point, epoch, vector = subpnt("NEAR POINT/ELLIPSOID", "EARTH", et, "IAU_EARTH", "NONE", "PROBE")
radius, lon, lat = reclat(point)
(lon * dpr(), lat * dpr())    # degrees
```

## Search for an event

The geometry finder calls (names that begin with `gf`) search an interval of time for the epochs at which a quantity meets a condition. The interval and the result are both a [`Window`](@ref), which is a list of disjoint intervals. The search below looks for the local minima (`"LOCMIN"`) of the distance from the Earth during one day. The two zeros are a reference value and an adjustment, which this relation does not use, and 600 seconds is the step of the search; [The step](@ref) says how to choose it.

```@example first
search = Window()
wninsd!(search, et, et + 86400.0)
perigees = gfdist("PROBE", "NONE", "EARTH", "LOCMIN", 0.0, 0.0, 600.0, search)

for i in 1:wncard(perigees)
    t, _ = wnfetd(perigees, i)
    distance = vnorm(spkpos("PROBE", t, "J2000", "NONE", "EARTH")[1])
    println(timout(t, "YYYY Mon DD HR:MN:SC ::TDB"), "   ", round(distance; digits=3), " km")
end
```

The orbit of the example has a perigee radius of 7078.137 km, which the search recovers three times in the day.

## Unload

[`kclear`](@ref) removes all kernels, and [`unload`](@ref) removes one.

```@example first
kclear()
ktotal("ALL")
```

## Errors

A call that cannot answer throws a [`SpiceError`](@ref) whose `kind` field names the reason and whose `msg` field holds the text; CSPICE would set its error flag instead. With no kernel loaded, the name of the spacecraft is unknown:

```@example first
try
    spkezr("PROBE", et, "J2000", "NONE", "EARTH")
catch err
    err.kind
end
```

The task pages in the menu continue from here, and [Differences from CSPICE](@ref) lists what a CSPICE user must know before a port.
