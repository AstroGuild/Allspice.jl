# Allspice.jl

Allspice.jl is a native Julia port of the NAIF SPICE toolkit. The package keeps the call model of CSPICE, in which a program loads kernels and then queries states and orientations without any earlier description of the bodies or the frames.

The port adds two capabilities that CSPICE does not have. The first is analytical time derivatives of second and higher order, which the package computes from the kernel coefficients where a closed form exists and by differentiation of the evaluation algorithm where it does not. The second is thread safety, since all loaded kernel data belong to an explicit context and no hidden global state exists.

## Status

The package is a prototype that covers the part of the toolkit which most trajectory work needs, and the table gives the state of each area. The test suite compares every result with CSPICE through SPICE.jl, on kernels from NAIF and on kernels that the tests write with the CSPICE writers; `CONTRACT.md` records the internal interfaces.

| Area | Available | Not yet available |
|---|---|---|
| Kernel loading | SPK, CK, binary PCK, DSK, EK, text kernels (leap seconds, PCK, frames, spacecraft clock), meta-kernels | |
| SPK data types | all types of CSPICE N0067: 1, 2, 3, 5, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21; `getelm`, `evsgp4` | |
| States | `spkezr`, `spkpos`, `spkgeo`, `spkssb`, `spkltc`, `spkaps`, `spkacs`, `spkapo`, with all aberration corrections of CSPICE; constant-position and constant-velocity observers and targets (`spkcpo`, `spkcpt`, `spkcvo`, `spkcvt`) | |
| Frames | inertial, body-fixed (text and binary PCK), fixed-offset, CK-based, dynamic (all families) and switch frames; `pxform`, `sxform`; `frinfo`, `cidfrm`, `cnmfrm`, `ccifrm`, `tkfram`, `tipbod`, `tisbod`, `pckfrm`, `pckcov` | |
| Derivatives | `spkder` and `pxder` to any order; a `ForwardDiff.Dual` epoch in the public calls | |
| Orientation and clocks | CK types 1 to 6, `ckgp`, `ckgpav`, `ckcov`, `ckobj`; spacecraft clock type 1 | |
| Time | `str2et`, `utc2et`, `timout`, `et2utc`, `etcal`, `unitim`, `tparse`, `tpictr`, `timdef`, `tsetyr`, `deltet`, `et2lst`, `lspcn` | |
| File access and writers | DAF read and write calls, SPK writers for all 16 types (type 14 through `spk14b`, `spk14a`, `spk14e`), CK writers for 5 types, comment area, windows, `spkobj`, `spkcov`, `kinfo`, `kdata` | |
| Kernel pool | `gdpool`, `gcpool`, `gipool`, `bodvrd`, `bodvcd`, `bodfnd`, body and frame name calls; `pdpool`, `pcpool`, `pipool`, `dvpool`, `lmpool`, `ldpool`, `clpool`, `boddef`; `expool`, `dtpool`, `gnpool`, `szpool`, `stpool`; watchers `swpool`, `cvpool` for each context | |
| Mathematics | vectors, matrices, rotations and quaternions, coordinate conversions with their Jacobians, ellipsoid, plane and ellipse calls, conic elements and two-body propagation | |
| Shape kernels | DAS file access, DSK type 2 read and write, plate access calls, ray-surface intercept (`dskxv`, `dskxsi`, `dskx02`), `latsrf`, `srfnrm`, surface name calls (type 2 is the only DSK type of CSPICE N0067) | |
| Event kernels | EK read access (`ekopr`, `ekssum`, `ekrcec`, `ekrced`, `ekrcei`), table catalogue (`ekntab`, `ektnam`, `ekccnt`, `ekcii`), queries (`ekpsel`, `ekfind`, `eknelt`, `ekgc`, `ekgd`, `ekgi`) with a result object for each query, record and fast writers (`ekopn`, `ekops`, `ekbseg`, `ekappr`, `ekinsr`, `ekdelr`, `ekace*`, `ekuce*`, `ekifld`, `ekacl*`, `ekffld`) | |
| Geometry finder | `gfdist`, `gfsep`, `gfposc`, `gfrr`, `gfpa`, `gfsubc`, `gfsntc`, `gfilum`, `gfoclt`, `gftfov`, `gfrfov`, `gfuds`, `gfudb`, `gfevnt` for the six quantities of CSPICE, `gfocce`, `gffove`, `gfstol`, `gfrefn`, `uddf`, `uddc` | DSK shapes in `gfsubc`, `gfsntc` and `gfilum` (CSPICE has none either) |
| Geometry | `subpnt`, `subslr`, `sincpt`, `ilumin`, `illumg`, `illumf`, `phaseq`, `occult`, `limbpt`, `termpt`, `tangpt`, `getfov`, `getfvn`, `fovtrg`, `fovray`, `azlcpo`, `azlrec`, `recazl`, `dazldr`, `drdazl`, `srfrec`, on ellipsoids and on DSK data | |
| Ephemeris and rotation tables | `ephemeris_table`, `table_position`, `table_state`, `table_jet`: the position of one body relative to another on a fixed interval, as one array of coefficients for the CPU or a GPU; `table_to_spk`: the table as an SPK segment of type 2; `rotation_table`, `table_rotation`, `table_rotation_jet`: the same for the rotation between two frames | |

Three properties differ from CSPICE by design. All loaded data belong to a `SpiceContext`, whose readers take no lock, so that queries from many threads are safe while another thread loads or unloads kernels. A file that is open for access is an object that the caller owns (`h = dafopr(path)`), and the calls that use a hidden current file in CSPICE take that object as their first argument. An answer never depends on earlier calls, because the package keeps no buffers or saved search state between queries.

## Ephemeris and rotation tables

An optimizer or a GPU kernel that needs the position of a body many times on a known interval of time does not need the kernels in its inner loop. `ephemeris_table` extracts the position of one body relative to another on that interval into an `EphemerisTable`, which is one flat array of Chebyshev coefficients on intervals of equal length together with a few numbers, and `table_position`, `table_state` and `table_jet` evaluate it with no kernel set, no segment search and no chain of bodies.

```julia
tab = ephemeris_table("EARTH", t0, t1, "J2000", "NONE", "SUN")
x = table_state(tab, et)                    # position (km) and velocity (km/s)
r, v, a = table_jet(tab, et, Val(2))        # derivatives of any order

using CUDA, Adapt
gpu = adapt(CuArray, tab)                   # `cu(tab)` gives Float32 coefficients
states = map(et -> @inbounds(table_state(gpu, et)), CuArray(epochs))
```

The function has two methods. For a geometric position in an inertial frame from SPK type 2 data, which is the case of the planetary ephemerides, the position is one polynomial for each record of the kernels; the table then takes the record length and the degree of the kernels and agrees with `spkpos` to rounding. In all other cases (a rotating frame, an aberration correction, another SPK type) the function interpolates and doubles the number of intervals until the measured error is below the tolerance `tol`, and the field `maxerr` of the table records the error that it measured. The velocity of a table is always the derivative of its position, and an epoch can be a `ForwardDiff.Dual` number, also on the GPU.

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

A query checks that the epoch is inside the table and throws `SpiceError(:nocoverage)` if it is not. The check is a bounds check, so `@inbounds` removes it, as a GPU kernel requires; an epoch outside the table then gives the extrapolation of the nearest interval.

A table of a geometric position can leave the session as a kernel: `table_to_spk(handle, tab, segid)` adds it to an SPK file as one segment of type 2 whose records are the records of the table, so that `spkezr` of this package and of CSPICE return the state of the table to rounding (measured 2e-16 relative for the Moon from the Earth on 200 days).

`rotation_table(from, to, t0, t1)` does for `pxform` what `ephemeris_table` does for `spkpos`. It interpolates the nine elements of the rotation matrix, `table_rotation` returns the matrix, and `table_rotation_jet` returns the matrix with its time derivatives, of which the first is the derivative block of `sxform`. The matrix is not made orthogonal again, because its distance from an orthogonal matrix is of the size of the measured error of the table, which the next table gives for six days from 2022-03-08 with the default arguments.

| From | To | Intervals | Measured error of the table | Difference from CSPICE `pxform` |
|---|---|---:|---:|---:|
| `ECLIPJ2000` | `J2000` | 1 | 0 | 0 |
| `IAU_MOON` | `J2000` | 1 | 5.8e-13 | 7.5e-13 |
| `J2000` | `IAU_EARTH` | 16 | 1.5e-11 | 2.0e-11 |
| `MOON_PA` | `IAU_EARTH` | 16 | 1.6e-11 | 1.9e-11 |

The errors are Frobenius norms (the root of the sum of the squares of the nine differences) and have no unit. The two tables that reach the Earth stop at 1.5e-11 because the rotation angle of the Earth is a number of many turns, whose rounding limits the source.

## Relation to NAIF

This package is an independent work. NAIF does not develop it, support it or endorse it, and the CSPICE toolkit of NAIF is the reference against which the package is tested.

## Licence

The package is under the MIT licence; see `LICENSE`.
