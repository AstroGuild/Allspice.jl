# Allspice.jl

NASA's Navigation and Ancillary Information Facility (NAIF) develops
[SPICE](https://naif.jpl.nasa.gov/naif/), a system for computing spacecraft and
planetary geometry from data files called kernels; CSPICE is its C implementation.
Allspice.jl is a native Julia port of CSPICE N0067 that reads these kernels and
performs the same position, orientation, time, and event calculations without a C
library.

Choose Allspice.jl instead of a CSPICE wrapper when a native Julia implementation
matters: separate kernel contexts permit concurrent queries without a lock, epochs
can carry ForwardDiff dual numbers, and trajectory or rotation tables can run in
inner loops and on GPUs. If you need to call NAIF's CSPICE library itself, use a
wrapper such as SPICE.jl.

Calls keep CSPICE names and argument order. Results are Julia return values, and
failures throw exceptions, so programs written with SPICE.jl usually need few
changes.

Allspice.jl is a prototype, and its interface can still change. Review its
[coverage and validation](docs/src/explanation/coverage.md) before you rely on it.

## Installation

Allspice.jl requires Julia 1.10 or later. Install it directly from GitHub:

```julia-repl
pkg> add https://github.com/EllissoideRotondo/Allspice.jl
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

## Documentation

Start with the [first session](docs/src/tutorial.md) for a complete workflow from
loading kernels to searching for an event. Before porting CSPICE or SPICE.jl code,
read the [differences from CSPICE](docs/src/explanation/cspice.md); the manual also
groups the [API reference by subject](docs/src/reference/).

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
