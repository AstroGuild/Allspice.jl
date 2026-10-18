# Allspice.jl

Allspice.jl reads SPICE kernels and answers the questions that the CSPICE toolkit answers (where a body is, how a frame is oriented, when an event occurs), in Julia code with no C library behind it. The calls keep the names and the argument order of CSPICE, so that a program written for CSPICE or SPICE.jl needs few changes.

Queries are safe from many threads at once, because all loaded data belong to a context object and no call keeps hidden state. The epoch can be a dual number of ForwardDiff.jl, and `spkder` and `pxder` return time derivatives of any order, so that a state can go into an optimizer or an integrator without finite differences. A trajectory or a rotation can be extracted into a flat array of coefficients that a GPU evaluates.

## Installation

The package needs Julia 1.10 or later.

```julia
using Pkg
Pkg.add("Allspice")
```

The examples of this manual run on three small kernels that the manual writes itself: leap seconds, the constants of the Earth, and two-body orbits of the Earth around the Sun and of a spacecraft named `PROBE` around the Earth. The positions are therefore not those of the real bodies; with kernels from NAIF the same calls give real positions.

## Contents

### Tutorial

The first session loads kernels, computes states and orientations, and ends with an event search.

```@contents
Pages = ["tutorial.md"]
Depth = 2
```

### How-to guides

These guides cover several threads and kernel sets, time derivatives, tables for inner loops and GPUs, event searches, and the writing of kernels.

```@contents
Pages = [
    "howto/contexts.md",
    "howto/derivatives.md",
    "howto/tables.md",
    "howto/events.md",
    "howto/writing.md",
]
Depth = 2
```

### Concepts

These pages say where Allspice.jl differs from CSPICE, which parts of the toolkit exist, and how they were checked.

```@contents
Pages = ["explanation/cspice.md", "explanation/coverage.md"]
Depth = 2
```

### API reference

The reference groups the calls by subject.

```@contents
Pages = [
    "reference/kernels.md",
    "reference/time.md",
    "reference/states.md",
    "reference/frames.md",
    "reference/geometry.md",
    "reference/events.md",
    "reference/math.md",
    "reference/files.md",
    "reference/ek.md",
]
Depth = 1
```
