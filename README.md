# Allspice.jl

[![CI](https://github.com/AstroGuild/Allspice.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/AstroGuild/Allspice.jl/actions/workflows/CI.yml)
[![Manual](https://img.shields.io/badge/docs-dev-blue.svg)](https://astroguild.github.io/Allspice.jl/dev/)

Allspice.jl is a native Julia port of CSPICE N0067 for the computation of spacecraft
and planetary geometry. It reads the data files of
[SPICE](https://naif.jpl.nasa.gov/naif/), called kernels, which NASA's Navigation and
Ancillary Information Facility (NAIF) publishes, and performs the position,
orientation, time and event calculations of CSPICE, the C implementation of SPICE,
in Julia code with no C library behind it.

Features include:

- The calls of CSPICE with their names and argument order, with Julia return values
  and exceptions, so that a program written with SPICE.jl usually needs few changes.
- Separate kernel contexts, so that several threads can query at once without a lock.
- Epochs that can be dual numbers of ForwardDiff.jl, and time derivatives of any
  order with `spkder` and `pxder`.
- Ephemeris and rotation tables for inner loops and GPU functions.
- Event searches, surface geometry, and the writing of SPK, CK, DSK and EK kernels.

Allspice.jl is meant for programs that need SPICE geometry as native Julia code. A
program that must call the CSPICE library of NAIF itself is better served by a
wrapper such as SPICE.jl.

## Installation

Allspice.jl requires Julia 1.10 or later. Install it directly from GitHub:

```julia-repl
pkg> add https://github.com/AstroGuild/Allspice.jl
```

## Quick start

Download the NAIF generic kernels
[`naif0012.tls`](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls)
and
[`de440s.bsp`](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp),
then save both files in Julia's working directory or replace the filenames below
with their paths:

```julia
using Allspice

furnsh("naif0012.tls")
furnsh("de440s.bsp")

et = str2et("2025 MAR 01 12:00:00 UTC")  # seconds past J2000 TDB
state, light_time = spkezr("MOON", et, "J2000", "LT+S", "EARTH")

position = state[1:3]  # km
velocity = state[4:6]  # km/s
light_time             # one-way light time, s
```

`spkezr` returns the state of the Moon relative to the Earth after light-time and
stellar-aberration corrections. Allspice.jl returns Julia values instead of filling
C output arguments, uses indices that start at 1, and throws [`SpiceError`](docs/src/reference/kernels.md)
when a call cannot produce a result.

## Limits

- Allspice.jl is a prototype. Its interface can still change.
- It ports the library calls of CSPICE N0067. Utilities that Julia Base covers (cells, sets, strings, error handling) are not ported; see [coverage and validation](docs/src/explanation/coverage.md).
- A result differs from CSPICE by rounding only. The measured differences are in [coverage and validation](docs/src/explanation/coverage.md#Known-limits-of-accuracy).

## Tutorials and Documentation

The manual starts with a [first session](docs/src/tutorial.md) that goes from
loading kernels to an event search. Its how-to guides then cover kernel contexts and
threads, time derivatives, ephemeris and rotation tables, event searches and the
writing of kernels; the concept pages give the
[differences from CSPICE](docs/src/explanation/cspice.md), which should be read
before CSPICE or SPICE.jl code is ported, and the coverage and validation of the
port, and the [API reference](docs/src/reference/) groups the calls by subject.

To build the manual, run from the repository root:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Then open `docs/build/index.html`.

## Relationship to NAIF

Allspice.jl is an independent project that NAIF and the Jet Propulsion Laboratory
do not develop, support, or endorse. Report defects through this repository, not to
NAIF. The repository does not redistribute CSPICE or SPICE kernels; CSPICE remains
available from [NAIF](https://naif.jpl.nasa.gov/naif/toolkit.html) under its [rules
of use](https://naif.jpl.nasa.gov/naif/rules.html).

Publications that use SPICE should cite the references requested by NAIF:

- C. H. Acton, "Ancillary data services of NASA's Navigation and Ancillary
  Information Facility", *Planetary and Space Science* 44 (1), 65–70, 1996.
- C. Acton, N. Bachman, B. Semenov, and E. Wright, "A look towards the future in
  the handling of space science mission geometry", *Planetary and Space Science*
  150, 9–12, 2018.

## Licence

Allspice.jl is available under the [MIT licence](LICENSE). This licence covers only
the work in this repository and grants no rights over CSPICE.
