# Allspice.jl

Allspice.jl reads SPICE kernels and answers the questions that the CSPICE toolkit answers (where a body is, how a frame is oriented, when an event occurs), in Julia code with no C library behind it. The calls keep the names and the argument order of CSPICE, so that a program written for CSPICE or SPICE.jl needs few changes.

Queries are safe from many threads at once, because each query reads an immutable snapshot of the loaded kernels, which a later load does not change, and no call keeps hidden state. A context object holds one set of loaded kernels, and a program that needs two independent sets makes two contexts. With ForwardDiff.jl loaded, the epoch of the state, rotation, surface and time-conversion calls that [Dual-number epochs](@ref) lists can be a dual number, and `spkder` and `pxder` return time derivatives of any order, so that a state can go into an optimizer or an integrator without finite differences. A trajectory or a rotation can be extracted into a flat array of coefficients that a GPU evaluates.

A program that must call the CSPICE library of NAIF itself, for example to obtain results that are identical to those of another CSPICE program on every platform, is better served by the wrapper SPICE.jl. The toolkit programs of NAIF, such as `brief` and `mkspk`, are not part of this package.

## Installation

The package needs Julia 1.10 or later. Install it from GitHub:

```julia-repl
pkg> add https://github.com/Orbitology/Allspice.jl
```

The examples of this manual run on three small kernels that the manual writes itself: leap seconds; the radii and the rotation of the Sun and of the Earth, with the name of a spacecraft `PROBE`; and two-body orbits of the Earth around the Sun and of `PROBE` around the Earth. The positions are therefore not those of the real bodies; with kernels from NAIF the same calls give real positions.

## Contents

### Tutorial

The first session loads kernels, computes states and orientations, and ends with an event search.

```@contents
Pages = ["tutorial.md"]
Depth = 2
```

### How-to guides

These guides cover several threads and several sets of loaded kernels, time derivatives, tables for inner loops and GPUs, event searches, and the writing of kernels.

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

## Contributing

Report defects and ask questions in the [issues of the repository](https://github.com/Orbitology/Allspice.jl/issues). Allspice.jl is an independent project: NAIF and the Jet Propulsion Laboratory do not develop or support it, so do not send its defects to NAIF.
