# Time derivatives

CSPICE returns a position and its first derivative. This package returns derivatives of any order, in two ways that give the same numbers to rounding: a direct call, and automatic differentiation through the ordinary calls.

```@setup der
using Allspice
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Direct calls

[`spkder`](@ref) returns the position of a target relative to an observer and its first `N` time derivatives, and [`pxder`](@ref) does the same for a rotation matrix. The order is a `Val`, so that the result is a tuple of static arrays and the call allocates nothing.

```@example der
using Allspice

et = str2et("2025 MAR 01 12:00:00 UTC")
position, velocity, acceleration = spkder("PROBE", et, "J2000", "EARTH", Val(2))
acceleration                                # geometric, km/s^2
```

Without a correction argument the result is geometric, and its second element is the velocity that `spkezr` returns. An aberration correction goes before the observer:

```@example der
spkder("PROBE", et, "J2000", "LT+S", "EARTH", Val(2))[3]    # as seen from the Earth, km/s^2
```

```@example der
rotation, rate = pxder("J2000", "IAU_EARTH", et, Val(1))
rate[1, 1]                                  # 1/s
```

## Dual-number epochs

The epoch of every query can be a `ForwardDiff.Dual`: `spkezr`, `spkpos`, `pxform`, `sxform`, and also `spkder`, `pxder` and the table calls. A function that calls them can therefore be differentiated as a whole, and an optimizer can differentiate through an acceleration from `spkder`. The derivative of the corrected position below includes the change of the light time.

```@example der
using ForwardDiff

apparent(t) = spkpos("PROBE", t, "IAU_EARTH", "LT+S", "EARTH")[1]
ForwardDiff.derivative(apparent, et)        # km/s in the rotating frame
```

## Where the derivatives come from

For data that are polynomials in time (SPK types 2, 3, 8, 9, 12, 13, 14, 18, 20 and the binary PCK), the derivative of any order is the exact derivative of the polynomial. For the other types, and for frames that are defined by formulas, the package applies truncated Taylor series arithmetic to the evaluation algorithm, which gives the derivatives of that algorithm without finite differences. In both cases a derivative has the smoothness of the data and no more. At the boundary between two records of an SPK, the geometric position is continuous to rounding and the higher derivatives can jump; an aberration-corrected position has steps of some millimetres there, because the target is evaluated at an earlier epoch than the observer and the two cross a boundary at different times.
