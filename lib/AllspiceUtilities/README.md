# AllspiceUtilities.jl

AllspiceUtilities.jl provides Julia versions of five programs of the NAIF toolkit, built
on [Allspice.jl](https://github.com/AstroGuild/Allspice.jl) with no C library behind
them: `brief`, `commnt`, `spkmerge`, `spkdiff` and `frmdiff`. Each program is a set of
functions that return structured results, such as a coverage table or the states of a
comparison, so that a Julia program can use the results directly; a printer gives a text
report in the layout of the NAIF program where one is useful. The results agree with the
NAIF programs of CSPICE N0067, which serve as the oracle of the tests, but the text of
the reports is not the same byte for byte.

The function names follow the programs: `brief` and `printbrief` summarize coverage,
`comments`, `addcomments` and `deletecomments` do the work of `commnt`, `spkmerge` takes
a list of `MergeSource`s, and `spkdiff` and `frmdiff` return comparisons that
`basicstats` and, for trajectories, `viewstats` summarize.

## Installation

AllspiceUtilities.jl requires Julia 1.11 or later. Install it from the `lib` directory of
the Allspice.jl repository:

```julia-repl
pkg> add https://github.com/AstroGuild/Allspice.jl:lib/AllspiceUtilities
```

## Quick start

```julia
using AllspiceUtilities

# Coverage of each body of an SPK file.
printbrief(brief("de440s.bsp"))

# The Moon from DE440 against DE405, at 100 epochs of their common coverage.
d = spkdiff("de440s.bsp", "de405.bsp"; body="MOON", center="EARTH BARYCENTER",
    kernels="naif0012.tls", n=100)
basicstats(d).maxpos          # largest position difference, km

# A file with the Moon from DE440 for 2000 to 2010 and from DE405 for 2010 to 2020.
spkmerge("moon.bsp", [
    MergeSource("de440s.bsp"; bodies=[301], start=0.0, stop=3.15576e8),
    MergeSource("de405.bsp"; bodies=[301], start=0.0, stop=6.31152e8),
])
addcomments("moon.bsp", ["Moon relative to the Earth-Moon barycenter."])
```

## Differences from the NAIF programs

- `brief` takes SPK and binary PCK files in one call.
- `frmdiff` keeps the last epoch before the stop time, which the NAIF program leaves out.
- `spkmerge` writes one file per call and takes epochs in seconds past J2000 TDB, so it
  needs no leapseconds kernel.
