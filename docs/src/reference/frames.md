# Frames and orientation

## Rotations between frames

```@docs
pxform
pxfrm2
sxform
```

## Frame information

```@docs
frinfo
cidfrm
cnmfrm
ccifrm
tkfram
tipbod
tisbod
pckfrm
pckfrm!
pckcov
pckcov!
bltfrm
kplfrm
```

## Orientation from C-kernels

```@docs
ckgp
ckgpav
ckcov
ckcov!
ckobj
ckobj!
ckmeta
ckfrot
ckfxfm
```

## Records of C-kernel segments

These calls read one segment of an open CK file, given by its descriptor from a DAF
search ([`dafbfs`](@ref), [`daffna`](@ref), [`dafgs`](@ref)).

```@docs
cknr02
ckgr02
cknr03
ckgr03
```
