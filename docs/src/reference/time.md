# Time

A call that needs the leap-second kernel and finds none throws a kind that is the short
message of the CSPICE routine that the call ports: `:kernelvarnotfound` from
[`deltet`](@ref), `:noleapseconds` from [`str2et`](@ref) and [`utc2et`](@ref), and
`:missingtimeinfo` from [`et2utc`](@ref), [`timout`](@ref), [`unitim`](@ref) and the
spacecraft clock calls such as [`sct2e`](@ref).

## Conversion

```@docs
deltet
utc2et
str2et
et2utc
etcal
timout
tpictr
tparse
timdef
unitim
tsetyr
tparch
jul2gr
gr2jul
```

## Spacecraft clocks

```@docs
sct2e
sce2c
sce2t
sce2s
scs2e
scencd
scdecd
sctiks
scfmt
scpart
sctype
```

## Local time

```@docs
et2lst
lspcn
```

## Constants

```@docs
j1900
j1950
j2000
j2100
jyear
tyear
spd
b1900
b1950
```
