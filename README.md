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
- Queries from any number of threads at once without a lock, because each query reads
  an immutable snapshot of the loaded kernels, which a later load does not change.
- Kernel contexts, which are objects that each hold a separate set of loaded kernels,
  for a program that compares two ephemerides or keeps two missions apart.
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
[`naif0012.tls`](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls),
[`de440s.bsp`](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp)
and
[`pck00011.tpc`](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc).
Save the files in Julia's working directory, or replace the filenames below with their
paths. The last file holds the orientation and the radii of the Sun, the planets and
their satellites, which body-fixed frames such as `IAU_EARTH` need.

```julia
using Allspice

furnsh("naif0012.tls")    # leap seconds
furnsh("de440s.bsp")      # planets, Moon and Sun
furnsh("pck00011.tpc")    # orientation and radii of the bodies

et = str2et("2025 MAR 01 12:00:00 UTC")  # seconds past J2000 TDB
state, light_time = spkezr("MOON", et, "J2000", "LT+S", "EARTH")

position = state[1:3]  # km
velocity = state[4:6]  # km/s
light_time             # one-way light time, s
```

`spkezr` returns the state of the Moon relative to the Earth after light-time and
stellar-aberration corrections. Allspice.jl returns Julia values instead of filling
C output arguments, uses indices that start at 1, and throws [`SpiceError`](https://astroguild.github.io/Allspice.jl/dev/reference/kernels/#Allspice.SpiceError)
when a call cannot produce a result.

The package [SPICEKernels.jl](https://github.com/JuliaAstro/EphemerisSources.jl/tree/main/lib/SPICEKernels)
downloads the generic kernels on first use. Install it with `pkg> add SPICEKernels`,
then load the same three files by name:

```julia
using SPICEKernels

furnsh(naif0012_tls())
furnsh(de440s())
furnsh(pck00011())
```

## Limits

- Allspice.jl is a prototype. Its interface can still change.
- It ports the library calls of CSPICE N0067. Utilities that Julia Base covers (cells, sets, strings, error handling) and the toolkit programs are not ported; the [differences from CSPICE](https://astroguild.github.io/Allspice.jl/dev/explanation/cspice/#What-is-absent) name the call to use in place of each.
- A result is the result that CSPICE gives in a new process, to rounding. The measured differences are in [coverage and validation](https://astroguild.github.io/Allspice.jl/dev/explanation/coverage/#Known-limits-of-accuracy).

## Tutorials and Documentation

The [manual](https://astroguild.github.io/Allspice.jl/dev/) starts with a first session
that goes from loading kernels to an event search, and continues with how-to guides,
concept pages and the API reference; before a program written for CSPICE or SPICE.jl is
ported, read the [differences from CSPICE](https://astroguild.github.io/Allspice.jl/dev/explanation/cspice/).

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

## Contributing

Report defects and ask questions in the
[issues of the repository](https://github.com/AstroGuild/Allspice.jl/issues). Run the
tests from the repository root, on 8 processes:

```sh
julia --startup-file=no --threads=8 test/parallel.jl
```

Build the manual, then open `docs/build/index.html`:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## License

Allspice.jl is available under the [MIT license](LICENSE). This license covers only
the work in this repository and grants no rights over CSPICE.
