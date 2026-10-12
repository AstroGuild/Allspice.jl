# NeoSpice.jl

NeoSpice.jl is a native Julia port of the NAIF SPICE toolkit, and the name is temporary. The package keeps the call model of CSPICE, in which a program loads kernels and then queries states and orientations without any earlier description of the bodies or the frames.

The port adds two capabilities that CSPICE does not have. The first is analytical time derivatives of second and higher order, which the package computes from the kernel coefficients where a closed form exists and by differentiation of the evaluation algorithm where it does not. The second is thread safety, since all loaded kernel data belong to an explicit context and no hidden global state exists.

## Status

The package is a prototype. It loads SPK kernels of types 2 and 3, binary PCK kernels of type 2 and text kernels (leap seconds, PCK, frame kernels and meta-kernels), and it answers `spkezr`, `spkpos`, `spkgeo`, `pxform` and `sxform` in inertial, body-fixed and fixed-offset frames, with the aberration corrections of CSPICE (`LT`, `LT+S`, `CN`, `CN+S` and the transmission forms). The calls `spkder` and `pxder` return time derivatives of any order. The test suite compares every result with CSPICE through SPICE.jl, and `CONTRACT.md` records the internal interfaces.

## Relation to NAIF

This package is an independent work. NAIF does not develop it, support it or endorse it, and the CSPICE toolkit of NAIF is the reference against which the package is tested.
