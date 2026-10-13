# Allspice.jl

Allspice.jl reads SPICE kernels and answers the questions that the CSPICE toolkit of NAIF answers (where a body is, how a frame is oriented, when an event occurs), in Julia code with no C library behind it. The calls keep the names and the argument order of CSPICE N0067, and the results agree with CSPICE to rounding.

Use it where CSPICE is the limit of a Julia program:

| Need | What the package gives |
|---|---|
| Queries from many threads | All loaded data belong to a context object and no call keeps hidden state, so queries need no lock |
| Derivatives for an optimizer or an integrator | `spkder` and `pxder` return time derivatives of any order, and an epoch can be a `ForwardDiff.Dual` |
| Positions in an inner loop or on a GPU | `ephemeris_table` and `rotation_table` extract a trajectory or a rotation into one array of coefficients |
| No binary dependency | Readers and writers for SPK, CK, PCK, DSK and EK files are Julia code |

The package is a prototype: the interface can still change.

## Installation

Julia 1.10 or later.

```julia
using Pkg
Pkg.add("Allspice")
```

## Example

Load kernels, convert a time, and ask for a state. The kernel names are those of the generic kernels of NAIF; use the paths of your own copies.

```julia
using Allspice

furnsh("naif0012.tls")                      # leap seconds
furnsh("de440s.bsp")                        # planetary ephemeris

et = str2et("2025 MAR 01 12:00:00 UTC")     # seconds past J2000 TDB

state, lt = spkezr("MOON", et, "J2000", "LT+S", "EARTH")
state[1:3]                                  # position, km
state[4:6]                                  # velocity, km/s
lt                                          # one-way light time, s

r, v, a = spkder("MOON", et, "J2000", "EARTH", Val(2))    # acceleration, km/s^2
```

A call that cannot answer throws a `SpiceError`; results come back as return values and static arrays, and indices count from 1.

## Documentation

The manual is in `docs/`; build it with `julia --project=docs docs/make.jl` and open `docs/build/index.html`.

| Page | Content |
|---|---|
| `docs/src/tutorial.md` | A first session, from kernels to an event search |
| `docs/src/howto/` | Contexts and threads, derivatives, tables, event searches, writing kernels |
| `docs/src/explanation/cspice.md` | What behaves differently from CSPICE |
| `docs/src/explanation/coverage.md` | Which parts of the toolkit exist, and how the results were checked |
| `docs/src/reference/` | Every exported call, by subject |

`CONTRACT.md` records the internal interfaces for maintainers. Run the tests with `julia --threads=8 test/parallel.jl`; they need SPICE.jl and kernels from NAIF, whose paths `CONTRACT.md` lists.

## Relation to NAIF

This package is an independent work, which the Navigation and Ancillary Information Facility (NAIF) of the Jet Propulsion Laboratory does not develop, support or endorse. Questions and defect reports about the package go to its author through the issue tracker of the repository, and not to NAIF, because NAIF cannot support software that it did not write.

The package is derived from the CSPICE toolkit, version N0067: each routine follows the algorithm of the CSPICE routine of the same name, and the error messages keep the wording of CSPICE where a program can depend on it. The repository contains no file of the toolkit (no source, header, library or document of NAIF) and no kernel. The CSPICE toolkit itself is available from NAIF at <https://naif.jpl.nasa.gov/naif/toolkit.html>, and its rules of use are at <https://naif.jpl.nasa.gov/naif/rules.html>.

The SPICE system is the work of NAIF, and a publication that uses this package should acknowledge it with the two references that NAIF asks for:

- C. H. Acton, "Ancillary data services of NASA's Navigation and Ancillary Information Facility", *Planetary and Space Science* 44 (1), 65-70, 1996.
- C. Acton, N. Bachman, B. Semenov and E. Wright, "A look towards the future in the handling of space science mission geometry", *Planetary and Space Science* 150, 9-12, 2018.

## Licence

The code of this package is under the MIT licence; see `LICENSE`. That licence covers the work of the author and gives no right over the CSPICE toolkit, which stays under the terms of NAIF.
