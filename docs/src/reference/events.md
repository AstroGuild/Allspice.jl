# Event searches and windows

The geometry finder is the set of searches of CSPICE, the `gf` routines, that return the
window of epochs on which a geometric condition holds; a [`Window`](@ref) holds each
result, and the calls of the second section combine windows. Sibling calls can throw
different kinds for one condition, because each kind is the one of the CSPICE routine
that the call ports; see [`SpiceError`](@ref).

## Searches

```@docs
gfdist
gfsep
gfposc
gfrr
gfpa
gfuds
gfudb
gfevnt
gfstol
gfrefn
uddc
uddf
gfsubc
gfsntc
gfilum
gfoclt
gftfov
gfrfov
gfocce
gffove
```

## Windows

```@docs
Window
card
wncard
wninsd!
wnvald!
wnunid
wnintd
wndifd
wncomd
wnelmd
wnincd
wnfetd
wnsumd
wnreld
wnexpd!
wncond!
wnfild!
wnfltd!
wnextd!
```
