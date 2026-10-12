# NeoSpice.jl

NeoSpice.jl is a native Julia port of the NAIF SPICE toolkit, and the name is temporary. The package keeps the call model of CSPICE, in which a program loads kernels and then queries states and orientations without any earlier description of the bodies or the frames.

The port adds two capabilities that CSPICE does not have. The first is analytical time derivatives of second and higher order, which the package computes from the kernel coefficients where a closed form exists and by differentiation of the evaluation algorithm where it does not. The second is thread safety, since all loaded kernel data belong to an explicit context and no hidden global state exists.

## Status

The project is in the design stage and contains no code yet. The design documents are in `docs/superpowers/specs/`.

## Relation to NAIF

This package is an independent work. NAIF does not develop it, support it or endorse it, and the CSPICE toolkit of NAIF is the reference against which the package is tested.
