# Ephemeris and rotation tables

An optimizer or a GPU kernel that needs the position of a body many times on a known interval of time does not need the kernels in its inner loop. [`ephemeris_table`](@ref) extracts the position of one body relative to another on that interval into an [`EphemerisTable`](@ref), which is one flat array of Chebyshev coefficients on intervals of equal length together with a few numbers. The evaluation then needs no kernel set, no search for a segment and no chain of bodies.

```@setup tab
using Allspice
include(joinpath(@__DIR__, "..", "..", "fixture.jl"))
foreach(furnsh, DocKernels.write(mktempdir()))
```

## Build and evaluate

```@example tab
using Allspice

t0 = str2et("2025 MAR 01 12:00:00 UTC")
t1 = t0 + 86400.0
tab = ephemeris_table("PROBE", t0, t1, "J2000", "NONE", "EARTH"; tol=1e-6)
(tab.nint, tab.ncoef, tab.maxerr)           # intervals, coefficients, measured error in km
```

[`table_position`](@ref) returns the position, [`table_state`](@ref) the position and the velocity, and [`table_jet`](@ref) the position with its derivatives to any order.

```@example tab
et = t0 + 3600.0
table_state(tab, et) - spkezr("PROBE", et, "J2000", "NONE", "EARTH")[1]
```

```@example tab
position, velocity, acceleration = table_jet(tab, et, Val(2))
acceleration                                # km/s^2
```

## How the table is built

The function has two methods. For a geometric position in an inertial frame from SPK type 2 data, which is the case of the planetary ephemerides, the position is one polynomial for each record of the kernels; the table then takes the record length and the degree of the kernels and agrees with `spkpos` to rounding. In all other cases (a rotating frame, an aberration correction, another SPK type) the function interpolates and doubles the number of intervals until the measured error is below the tolerance `tol`, and the field `maxerr` of the table records the error that it measured. The velocity of a table is always the derivative of its position, and an epoch can be a `ForwardDiff.Dual` number, also on the GPU.

A query checks that the epoch is inside the table and throws `SpiceError(:nocoverage)` if it is not. The check is a bounds check, so `@inbounds` removes it, as a GPU kernel requires; an epoch outside the table then gives the extrapolation of the nearest interval.

## On a GPU

The table is a parametric type over its array, and Adapt.jl moves it to a device. The evaluation is plain arithmetic on the coefficients, so that it compiles inside a GPU kernel.

```julia
using CUDA, Adapt

gpu = adapt(CuArray, tab)                   # `cu(tab)` gives Float32 coefficients
states = map(et -> @inbounds(table_state(gpu, et)), CuArray(epochs))
```

## Measured accuracy and speed

The table below gives measured values for 400 days of de440s.bsp (Intel CPU, one thread; NVIDIA RTX 3090 with 4,000,000 epochs; `bench/table_gpu.jl`). The error is the largest difference from the kernels in km.

| Case | Intervals | Error (km) | Table (ns per epoch) | `spkezr` (ns per epoch) |
|---|---:|---:|---:|---:|
| Moon relative to Earth, J2000, `NONE`, CPU | 101 | 2.4e-10 | 36 | 152 |
| Earth relative to Sun, ECLIPJ2000, `NONE`, CPU | 101 | 1.2e-7 | 36 | 277 |
| Moon relative to Earth, J2000, `LT+S`, CPU | 128 | 4.3e-6 | 42 | 653 |
| Earth relative to Moon, MOON_PA, `NONE`, `tol=1e-6`, CPU | 1024 | 4.2e-7 | 42 | 385 |
| Sun relative to Earth, IAU_EARTH, `NONE`, CPU | 1024 | 3.3e-3 | 42 | 341 |
| Earth relative to Sun, J2000, GPU, `Float64` | 101 | 1.2e-7 | 1.2 | |
| Earth relative to Sun, J2000, GPU, `Float32` | 101 | 16 | 0.33 | |

Two limits come from the source and not from the table. The Earth has turned about 50,000 radians since the reference epoch, so that the rounding of its rotation angle moves the Sun by some metres in the IAU_EARTH frame, and no table can be more accurate than that. An aberration-corrected position has steps of some millimetres at the record boundaries of the kernels, because the light-time correction evaluates one body on each side of a boundary. With coefficients of type `Float32` a position keeps about seven digits, which is 16 km at the distance of the Sun.

## Save a table as a kernel

A table of a geometric position can leave the session as a kernel: [`table_to_spk`](@ref)`(handle, tab, segid)` adds it to an SPK file as one segment of type 2 whose records are the records of the table, so that `spkezr` of this package and of CSPICE return the state of the table to rounding (measured 2e-16 relative for the Moon from the Earth on 200 days).

```@example tab
path = joinpath(mktempdir(), "table.bsp")
handle = spkopn(path, "from a table", 0)
table_to_spk(handle, tab, "PROBE from a table")
spkcls(handle)
spkobj(path)
```

## Rotations

[`rotation_table`](@ref)`(from, to, t0, t1)` does for `pxform` what `ephemeris_table` does for `spkpos`. It interpolates the nine elements of the rotation matrix, [`table_rotation`](@ref) returns the matrix, and [`table_rotation_jet`](@ref) returns the matrix with its time derivatives, of which the first is the derivative block of `sxform`. The matrix is not made orthogonal again, because its distance from an orthogonal matrix is of the size of the measured error of the table, which the next table gives for six days from 2022-03-08 with the default arguments.

| From | To | Intervals | Measured error of the table | Difference from CSPICE `pxform` |
|---|---|---:|---:|---:|
| `ECLIPJ2000` | `J2000` | 1 | 0 | 0 |
| `IAU_MOON` | `J2000` | 1 | 5.8e-13 | 7.5e-13 |
| `J2000` | `IAU_EARTH` | 16 | 1.5e-11 | 2.0e-11 |
| `MOON_PA` | `IAU_EARTH` | 16 | 1.6e-11 | 1.9e-11 |

The errors are Frobenius norms (the root of the sum of the squares of the nine differences) and have no unit. The two tables that reach the Earth stop at 1.5e-11 because the rotation angle of the Earth is a number of many turns, whose rounding limits the source.

```@example tab
rtab = rotation_table("J2000", "IAU_EARTH", t0, t1)
maximum(abs, table_rotation(rtab, et) - pxform("J2000", "IAU_EARTH", et))
```
