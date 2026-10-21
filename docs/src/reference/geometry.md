# Geometry

Sibling calls can throw different kinds for one condition, such as an unknown frame or
radii that are not valid, because each kind is the one of the CSPICE routine that the
call ports; see [`SpiceError`](@ref).

## Points on a body

```@docs
subpnt
subslr
sincpt
srfrec
latsrf
srfnrm
```

## Illumination

```@docs
ilumin
illumg
illumf
phaseq
trgsep
```

## Occultation, limb and terminator

```@docs
occult
limbpt
termpt
edterm
tangpt
```

## Fields of view

```@docs
fovtrg
fovray
getfov
getfvn
```

## Azimuth and elevation

```@docs
azlcpo
azlrec
recazl
dazldr
drdazl
```
