# Allspice.jl

Allspice.jl reads SPICE kernels and answers the questions that the CSPICE toolkit answers (where a body is, how a frame is oriented, when an event occurs), in Julia code with no C library behind it. The calls keep the names and the argument order of CSPICE, so that a program written for CSPICE or SPICE.jl needs few changes.

Queries are safe from many threads at once, because all loaded data belong to a context object and no call keeps hidden state. The epoch can be a dual number of ForwardDiff.jl, and `spkder` and `pxder` return time derivatives of any order, so that a state can go into an optimizer or an integrator without finite differences. A trajectory or a rotation can be extracted into a flat array of coefficients that a GPU evaluates.

## Installation

The package needs Julia 1.10 or later.

```julia
using Pkg
Pkg.add("Allspice")
```

## Where to go

| Need | Page |
|---|---|
| A first session from kernels to an event search | [First session](@ref) |
| Use from several threads, or several kernel sets at once | [Contexts and threads](@ref) |
| Velocity, acceleration and higher derivatives | [Time derivatives](@ref) |
| Many positions on a known interval, in an inner loop or on a GPU | [Ephemeris and rotation tables](@ref) |
| Write an SPK, a CK or another kernel | [Writing kernels](@ref) |
| Times of closest approach, eclipses, visibility | [Event searches](@ref) |
| What behaves differently from CSPICE | [Differences from CSPICE](@ref) |
| Which parts of the toolkit exist, and how they were checked | [Coverage and validation](@ref) |
| One call in detail | the reference pages in the menu |

The examples of this manual run on three small kernels that the manual writes itself: leap seconds, the constants of the Earth, and two-body orbits of the Earth around the Sun and of a spacecraft named `PROBE` around the Earth. The positions are therefore not those of the real bodies; with kernels from NAIF the same calls give real positions.
