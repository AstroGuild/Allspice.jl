# Ephemeris and rotation tables

An optimizer or a GPU function that needs the position of a body many times on a known interval of time does not need the SPICE kernels in its inner loop. [`ephemeris_table`](@ref) extracts the position of one body relative to another on that interval into an [`EphemerisTable`](@ref), which is one flat array of Chebyshev coefficients on intervals of equal length. The evaluation then needs no loaded kernels, no search for a segment and no chain of bodies.

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
(tab.nint, tab.ncoef, tab.maxerr)           # intervals, coefficients per interval, measured error in km
```

The keyword `tol` is the accepted error in kilometers; without it the goal is 1e-12 times the largest distance in the table. [`table_position`](@ref) returns the position, [`table_state`](@ref) the position and the velocity, and [`table_jet`](@ref) the position with its derivatives to any order.

```@example tab
et = t0 + 3600.0
table_state(tab, et) - spkezr("PROBE", et, "J2000", "NONE", "EARTH")[1]
```

```@example tab
position, velocity, acceleration = table_jet(tab, et, Val(2))
acceleration                                # km/s^2
```

## How the table is built

The function builds the table in one of two ways, and [`ephemeris_table`](@ref) states the conditions of each. For a geometric position in an inertial frame from SPK type 2 data, which is the case of the planetary ephemerides, the position is one polynomial for each [record](@ref "Exact derivatives and the rates of CSPICE") of the kernels; the table then takes the record length and the degree of the kernels and agrees with `spkpos` to rounding. In all other cases (a rotating frame, an aberration correction, another SPK type) the function interpolates and doubles the number of intervals until the measured error is below the tolerance `tol`, and the field `maxerr` of the table records the error that it measured.

The velocity of a table is always the derivative of its position, so that it can differ from the velocity of `spkezr` with an `LT` correction and for SPK types that store the velocity apart from the position. An epoch can be a `ForwardDiff.Dual` number, also on the GPU.

A query checks that the epoch is inside the table and throws a [`SpiceError`](@ref) of kind `:spkinsuffdata` if it is not; this error is a [`SpiceKindError`](@ref), without the fields of [`SpkInsufficientDataError`](@ref). The check is a bounds check, so `@inbounds` removes it, as a GPU function requires; an epoch outside the table then gives the extrapolation of the nearest interval.

## On a GPU

The table is a parametric type over its array, and Adapt.jl moves it to a device. The evaluation is plain arithmetic on the coefficients, so that it compiles inside a GPU function. In the block below, `epochs` is a vector of `Float64` epochs. The epochs and the bounds of the table stay `Float64` also when the coefficients are `Float32`, because a `Float32` number of seconds past J2000 has a resolution of about one minute.

```julia
using CUDA, Adapt

gpu = adapt(CuArray, tab)                   # `cu(tab)` of CUDA.jl gives Float32 coefficients
states = map(et -> @inbounds(table_state(gpu, et)), CuArray(epochs))
```

## In a CasADi model

[CasADiCodegen.jl](https://github.com/EllissoideRotondo/CasADiCodegen.jl) generates Julia code from a model of CasADi, and its `callback` puts a Julia function into the model with the derivatives that CasADi needs. Because [`table_jet`](@ref) returns static vectors, a table can give the position of a body to a force model in generated code with no memory allocation per call, and the same call gives the time derivatives of the position for the Jacobian of the model. The block below is a sketch and not a complete program: the body is a satellite of the Earth that the Moon perturbs, and the times `t0` and `t1` and the gravitational parameters `GM_EARTH` and `GM_MOON` are values of the application.

```julia
using Allspice, CasADiCodegen, LinearAlgebra

tab = ephemeris_table("MOON", t0, t1, "J2000", "NONE", "EARTH")
order(k) = let v = Val(k); t -> table_jet(tab, t[1], v)[k+1]; end       # derivative k of the position
moon = callback(order(0), "moon", (1,), (3,); jacobian = (order(1), order(2)))

x, t = MX.sym("x", 6), MX.sym("t")
r, s = x[1:3], moon(t)
a = -GM_EARTH * r / norm(r)^3 + GM_MOON * ((s - r) / norm(s - r)^3 - s / norm(s)^3)
z = vertcat(x, t)
f = CasADiCodegen.Function("f", [x, t], [vertcat(x[4:6], a), jacobian(vertcat(x[4:6], a), z)])
g, g! = build_function(f; expression = Val{false})
```

The epoch `t` is in seconds past J2000 TDB, and `order(k)` returns the derivative `k` of the position in km/s^k. CasADi takes the first derivatives of the callback from `order(1)` and the second derivatives from `order(2)`.

## Measured accuracy and speed

The table below gives measured values for 400 days of de440s.bsp, from `bench/table_gpu.jl`: one thread of an Intel CPU, and an NVIDIA RTX 3090 with 4,000,000 epochs. The error is the largest difference from the kernels in km, and the tolerance is the default unless the row names one. `MOON_PA` is the principal-axis frame of the Moon from the lunar orientation kernel of NAIF.

| Case | Intervals | Error (km) | Table (ns per epoch) | `spkezr` (ns per epoch) |
|---|---:|---:|---:|---:|
| Moon relative to Earth, J2000, `NONE`, CPU | 101 | 2.4e-10 | 36 | 152 |
| Earth relative to Sun, ECLIPJ2000, `NONE`, CPU | 101 | 1.2e-7 | 36 | 277 |
| Moon relative to Earth, J2000, `LT+S`, CPU | 128 | 4.3e-6 | 42 | 653 |
| Earth relative to Moon, MOON_PA, `NONE`, `tol=1e-6`, CPU | 1024 | 4.2e-7 | 42 | 385 |
| Sun relative to Earth, IAU_EARTH, `NONE`, CPU | 1024 | 3.3e-3 | 42 | 341 |
| Earth relative to Sun, J2000, GPU, `Float64` | 101 | 1.2e-7 | 1.2 | |
| Earth relative to Sun, J2000, GPU, `Float32` | 101 | 16 | 0.33 | |

No table is more accurate than its source. The rows of the Sun in `IAU_EARTH` and of the Moon with `LT+S` reach two limits of the source, which [Known limits of accuracy](@ref) describes.

## Save a table as a kernel

[`table_to_spk`](@ref) adds a table of a geometric position to an SPK file as one segment of type 2, whose records are the records of the table. `spkezr` of this package and of CSPICE then return the state of the table to rounding; the measured difference is 2e-16 relative for the Moon from the Earth on 200 days.

```@example tab
path = joinpath(mktempdir(), "table.bsp")
handle = spkopn(path, "from a table", 0)
table_to_spk(handle, tab, "PROBE from a table")
spkcls(handle)
spkobj(path)
```

## Rotations

[`rotation_table`](@ref) does for `pxform` what `ephemeris_table` does for `spkpos`: it interpolates the nine elements of the rotation matrix. [`table_rotation`](@ref) returns the matrix, and [`table_rotation_jet`](@ref) returns the matrix with its time derivatives.

```@example tab
rtab = rotation_table("J2000", "IAU_EARTH", t0, t1)
maximum(abs, table_rotation(rtab, et) - pxform("J2000", "IAU_EARTH", et))
```

The matrix is not made orthogonal again, because its distance from an orthogonal matrix is of the size of the measured error of the table. The next table gives that error for six days from 2022-03-08 with the default arguments.

| From | To | Intervals | Measured error of the table | Difference from CSPICE `pxform` |
|---|---|---:|---:|---:|
| `ECLIPJ2000` | `J2000` | 1 | 0 | 0 |
| `IAU_MOON` | `J2000` | 1 | 5.8e-13 | 7.5e-13 |
| `J2000` | `IAU_EARTH` | 16 | 1.5e-11 | 2.0e-11 |
| `MOON_PA` | `IAU_EARTH` | 16 | 1.6e-11 | 1.9e-11 |

The errors are Frobenius norms (the root of the sum of the squares of the nine differences) and have no unit; the last column compares the table with `pxform` of CSPICE. The two tables that reach the Earth stop at 1.5e-11 because of the rounding of the rotation angle of the Earth (see [Known limits of accuracy](@ref)).
