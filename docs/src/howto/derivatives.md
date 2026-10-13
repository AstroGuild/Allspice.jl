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

## Exact derivatives and the rates of CSPICE

Order 0 of `spkder` and of `pxder` is the result of `spkpos` and of `pxform`, and every higher order is the exact derivative of that result as a function of the epoch. The first derivative that CSPICE returns (the velocity of `spkezr`, the lower left block of `sxform`) is in some cases an approximation, so order 1 can differ from it.

For a position, let `s` be -1 for reception and +1 for transmission. In the barycentric inertial frame, `spkder` differentiates

```
p(t) = r_targ(t + s ℓ(t)) - r_obs(t)
P(t) = R(t + s |pc(t)| / c) * A(p(t), s * v_obs(t))
```

where `ℓ` is the light time, `A` is the stellar aberration and `R` is the rotation from J2000 to the output frame. For the `LT` corrections `ℓ(t) = |r_targ(t) - r_obs(t)| / c`, the geometric light time, because the single iteration of CSPICE evaluates the target at that light time; for the `CN` corrections `ℓ` solves the light-time equation `ℓ(t) = |p(t)| / c`, and its derivatives satisfy that equation order by order. `A` is the identity without `+S`, and otherwise `A(p, v) = p + |p| sinφ v̂⊥ + |p| (cosφ - 1) p̂`, with `v⊥` the part of `v` perpendicular to `p` and `sinφ = |v⊥| / c`. `R` is constant for an inertial frame; for another frame its epoch is corrected for the light time to the center of the frame, and `pc` is `p` with that center in the place of the target.

The velocity of `spkezr` differs from order 1 for two reasons. CSPICE takes the rate of the light time from the converged light-time equation also for the `LT` corrections, which gives a relative difference near `|v_targ| / c` times the range rate over `c` (9e-9 for the Moon seen from the Earth). CSPICE also takes the acceleration of the observer for the `+S` corrections from a difference of velocities over ±1 s, which gives a difference at rounding level.

For a rotation, the table says what `pxder` differentiates and what `sxform` returns in its place.

| Frame | `pxder` differentiates | Derivative block of `sxform` |
|---|---|---|
| Inertial, PCK, fixed offset | the rotation | the same |
| CK (class 3) | the interpolation of the CK segment, composed with the clock conversion of `sce2c`; the stored angular velocity is not used | built from the stored angular velocity; only segments that have one are used |
| Dynamic (class 5), mean or true equator and equinox of date | the precession, nutation and obliquity series, with all terms | the nutation rates of CSPICE, which omit the drift of the coefficients |
| Dynamic, two-vector | each defining vector with its aberration correction, as in `spkder` | velocities of `spkezr`; a central difference for a velocity vector; zero for a constant vector |
| Dynamic, Euler | the polynomials of the angles | the same |
| Dynamic with `FREEZE_EPOCH` | a constant: zero relative to the base frame | zero |
| Dynamic with `ROTATION_STATE = 'INERTIAL'` | the rotating frame | zero |
| Switch (class 6) | the base that `pxform` selects at the epoch | the block of that base |

## Where the derivatives come from

For data that are polynomials in time (SPK types 2, 3, 8, 9, 12, 13, 14, 18, 20 and the binary PCK), the derivative of any order is the exact derivative of the polynomial. For the other types, and for frames that are defined by formulas, the package applies truncated Taylor series arithmetic to the evaluation algorithm, which gives the derivatives of that algorithm without finite differences. In both cases a derivative has the smoothness of the data and no more. At the boundary between two records of an SPK, the geometric position is continuous to rounding and the higher derivatives can jump; an aberration-corrected position has steps of some millimetres there, because the target is evaluated at an earlier epoch than the observer and the two cross a boundary at different times.
