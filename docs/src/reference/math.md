# Mathematics

The calls of this page take numbers and arrays only, and need no kernel, except
[`recpgr`](@ref), [`pgrrec`](@ref), [`drdpgr`](@ref), [`dpgrdr`](@ref) and
[`xfmsta`](@ref). These rules hold for all of them:

- Integer components are converted to `Float64`.
- A sum of products is evaluated term by term in index order, as in CSPICE, so that a
  `Float64` input gives the CSPICE result. The `*` of StaticArrays can round in another
  way.
- An input of the wrong size throws `SpiceError(:badarg)`.
- A call that cites a routine "of SPICELIB" ports a Fortran routine for which CSPICE has
  no `_c` wrapper.

## Vectors

```@docs
vadd
vsub
vscl
vminus
vdot
vcrss
vnorm
vhat
vdist
vsep
vrel
vproj
vperp
vlcom
vlcom3
vrotv
vzero
vequ
vpack
vupack
ucrss
unorm
vtmv
frame
```

## Matrices

```@docs
mxv
mtxv
mxm
mtxm
mxmt
xpose
xpose6
invstm
invert
invort
trace
isrot
ident
```

## Vectors and matrices of any dimension

The `g` calls take no dimension argument: the arguments give the dimension. Static
arguments give a static result; other arguments give a `Vector` or a `Matrix`. Two
vectors must have equal length.

```@docs
vaddg
vsubg
vsclg
vminug
vequg
vdotg
vnormg
vhatg
unormg
vdistg
vsepg
vrelg
vprojg
vlcomg
vzerog
vtmvg
mxvg
mtxvg
mxmg
mtxmg
mxmtg
xposeg
```

## Rotations and quaternions

```@docs
rotate
rotmat
rotvec
twovec
twovxf
eul2m
m2eul
eul2xf
xf2eul
axisar
raxisa
m2q
q2m
qxq
qdq2av
rav2xf
xf2rav
sharpr
```

## Derivatives of vector functions

```@docs
dvhat
dvdot
dvcrss
ducrss
dvnorm
dvsep
```

## Coordinates

```@docs
reclat
latrec
recrad
radrec
recsph
sphrec
reccyl
cylrec
cyllat
latcyl
cylsph
sphcyl
sphlat
latsph
recgeo
georec
recpgr
pgrrec
xfmsta
convrt
```

## Jacobians of coordinate conversions

Element `[i, j]` of each Jacobian is the derivative of output coordinate `i` with respect
to input coordinate `j`. Angles are in radians.

```@docs
drdlat
dlatdr
drdsph
dsphdr
drdcyl
dcyldr
drdgeo
dgeodr
drdpgr
dpgrdr
```

## Ellipsoids, planes and ellipses

An ellipsoid has its center at the origin and the semi-axes `a`, `b`, `c` along x, y, z.
The error kind for a semi-axis that is not positive follows CSPICE call by call
(`:badaxislength`, `:invalidaxislength`, `:invalidradii` or `:degeneratecase`), so read
it in the entry of the call.

```@docs
SpicePlane
SpiceEllipse
nearpt
dnearp
surfpt
surfnm
surfpv
npedln
npelpt
nplnpt
edlimb
inedpl
inelpl
inrypl
cgv2el
el2cgv
saelgv
nvc2pl
nvp2pl
psv2pl
pl2nvc
pl2nvp
pl2psv
pjelpl
ednmpt
edpnt
vprjp
vprjpi
```

## Conic orbits

```@docs
conics
oscelt
oscltx
eqncpv
prop2b
stmp03
```

## Polynomials and interpolation

A Chebyshev expansion is given by its coefficients and by `x2s = [midpoint, radius]` of
its interval, as in an SPK segment of type 2 or 3. A table for interpolation is given by
its abscissas, or by the first abscissa and the step of an equally spaced table; a
Hermite table holds the value and the derivative at each abscissa, in this order.

```@docs
chbval
chbint
chbder
chbigr
lgrint
lgrind
lgresp
hrmint
hrmesp
polyds
qderiv
rquad
diags2
```

## Hexadecimal strings

A double is written as hexadecimal digits of a mantissa in [1/16, 1) and an exponent of
16, so that the string holds every bit of the number.

```@docs
dp2hx
hx2dp
```

## Constants

```@docs
clight
halfpi
twopi
rpd
dpr
```
