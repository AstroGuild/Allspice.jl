# Files and writers

A file handle is a [`DAFHandle`](@ref) or a [`DASHandle`](@ref) object, not an integer.
The writers take the arguments of CSPICE in the same order, but leave out most counts and
flags that follow from the input: in [`ckw01`](@ref) and [`ckw03`](@ref), `avvs` is
`nothing` for a segment without angular velocity, in the place of the flag `avflag`. The
calls that continue a DAF array search or an array that is being added
([`daffna`](@ref), [`daffpa`](@ref), [`dafgs`](@ref), [`dafgn`](@ref), [`dafws`](@ref),
[`dafrs`](@ref), [`dafrn`](@ref), [`dafada`](@ref), [`dafena`](@ref)) take the handle,
which the C calls do not, because the handle holds that state.

A call that needs a context takes it as the optional first argument `[ctx,]`, a
[`SpiceContext`](@ref), and uses [`default_context`](@ref) without it. `[ctx | ks,]`
means that the call also accepts a `KernelSet`, the immutable snapshot of the kernels of
a context that `Allspice.kernelset(ctx)` returns. [`dskxv!`](@ref) takes the context
after its two output arguments.

The identifier of an SPK or CK segment (`segid`, `ident`) is the name of its DAF array,
which [`dafgn`](@ref) returns. An input of the wrong shape throws `SpiceError(:badarg)`.

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

## SPK segments

A loaded file has no handle, so [`spksfs`](@ref) returns the path of the file, which
[`spkpvn`](@ref) takes with the same context.

```@docs
spksfs
spkpvn
spkuds
spksub
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

## Binary PCK writer

```@docs
pckopn
pckw02
pckcls
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
dasudc
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

The spatial index of a type 2 segment divides a box around the plates into cubes of
equal size, the fine voxels, which form groups of `corscl^3`, the coarse voxels. The
voxel-plate list gives the plates that meet each fine voxel, and the vertex-plate list
the plates that share each vertex. The DSK Required Reading of NAIF
([https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/dsk.html](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/dsk.html))
describes the index.

```@docs
DSKDescr
DSKSegment
dskobj
dskobj!
dsksrf
dsksrf!
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
