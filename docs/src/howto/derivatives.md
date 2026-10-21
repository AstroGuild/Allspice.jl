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

With ForwardDiff.jl loaded, the epoch can be a `ForwardDiff.Dual` in these calls:

- the states: `spkezr`, `spkpos`, `spkgeo`, `spkssb`, the integer-ID forms, `spkltc`, `spkaps`, `spkapo`, `spkacs`, `spkcpo`, `spkcpt`, `spkcvo` and `spkcvt`;
- the rotations `pxform`, `sxform` and `pxfrm2`;
- `spkder`, `pxder` and the table queries;
- the surface calls `subpnt`, `subslr`, `sincpt`, `ilumin`, `illumg`, `illumf`, `phaseq` and `azlcpo`;
- the time conversions `unitim` and `deltet`.

Other calls, such as `limbpt`, `tangpt`, `occult` and `et2lst`, take a `Float64` epoch. An event search takes a window of `Float64` epochs.

A function built from the calls of the list can be differentiated as a whole, and an optimizer can differentiate through an acceleration from `spkder`. The derivative of the corrected position below includes the change of the light time.

```@example der
using ForwardDiff

apparent(t) = spkpos("PROBE", t, "IAU_EARTH", "LT+S", "EARTH")[1]
ForwardDiff.derivative(apparent, et)        # km/s in the rotating frame
```

## Order 1 and the velocity of `spkezr`

Order 0 of `spkder` and of `pxder` is the result of `spkpos` and of `pxform`. For a geometric position, order 1 is the velocity that `spkezr` returns, and each higher order is a derivative of that velocity. With an aberration correction, and for some frames, the first derivative that CSPICE returns (the velocity of `spkezr`, the lower left block of `sxform`) is an approximation, so order 1 can differ from it: by about 1e-8 relative for the Moon seen from the Earth with `LT`. [Exact derivatives and the rates of CSPICE](@ref) gives the function that is differentiated for each correction and each kind of frame.
