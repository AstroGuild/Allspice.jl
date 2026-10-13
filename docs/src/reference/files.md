# Files and writers

A file handle is a [`DAFHandle`](@ref) or a [`DASHandle`](@ref) object, not an integer.
The writers keep the positional arguments of CSPICE and take the context as the keyword
`ctx`; other calls take it as the optional first argument `[ctx,]`. Without it a call
uses [`default_context`](@ref). `[ctx | ks,]` means that the call also accepts a
`KernelSet`, the immutable snapshot of the kernels of a context that
`Allspice.kernelset(ctx)` returns.

An input of the wrong shape throws `SpiceError(:badarg)` in the SPK and CK writers, and
`ArgumentError` in the DAS and DSK calls that say so.

## DAF files

```@docs
DAFHandle
dafopr
dafopw
dafonw
dafcls
dafhsf
dafrfr
dafbfs
dafbbs
daffna
daffpa
dafgs
dafgn
dafps
dafus
dafgda
dafgsr
dafbna
dafada
dafena
dafws
dafrs
dafrn
dafac
dafec
dafdc
```

## SPK writers

```@docs
spkopn
spkopa
spkcls
spkpds
spkw01
spkw02
spkw03
spkw05
spkw08
spkw09
spkw12
spkw13
spkw15
spkw17
spkw18
spkw19
spkw20
spkw21
spk14b
spk14a
spk14e
SPK14Writer
spkw10
```

## Two-line element sets

[`spkw10`](@ref) writes element sets to an SPK segment.

```@docs
getelm
evsgp4
```

## CK writers

```@docs
ckopn
ckcls
ckw01
ckw02
ckw03
ckw05
ckw06
```

## DAS and DLA files

```@docs
DASHandle
dasopr
dasopw
dasonw
dasops
dascls
dasllc
dashfn
dasrfr
daslla
dashfs
dasrdd
dasrdi
dasrdc
dasadd
dasadi
dasadc
dasudd
dasudi
dasac
dasec
dasdc
DLADescr
dlaopn
dlabns
dlaens
dlabfs
dlafns
dlabbs
dlafps
```

## DSK files

```@docs
DSKDescr
DSKSegment
dskobj
dsksrf
dskgd
dskb02
dskz02
dskn02
dskp02
dskv02
dski02
dskd02
dskx02
dskxv
dskxv!
dskxsi
dskgtl
dskstl
dskopn
dskw02
dskmi2
dskrb2
dskcls
pltexp
pltnp
pltnrm
pltar
pltvol
```

## Surface names

```@docs
srfc2s
srfs2c
srfcss
srfscc
```
