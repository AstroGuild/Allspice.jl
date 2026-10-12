# NeoSpice.jl

NeoSpice.jl is a native Julia port of the NAIF SPICE toolkit, and the name is temporary. The package keeps the call model of CSPICE, in which a program loads kernels and then queries states and orientations without any earlier description of the bodies or the frames.

The port adds two capabilities that CSPICE does not have. The first is analytical time derivatives of second and higher order, which the package computes from the kernel coefficients where a closed form exists and by differentiation of the evaluation algorithm where it does not. The second is thread safety, since all loaded kernel data belong to an explicit context and no hidden global state exists.

## Status

The package is a prototype that covers the part of the toolkit which most trajectory work needs, and the table gives the state of each area. The test suite compares every result with CSPICE through SPICE.jl, on kernels from NAIF and on kernels that the tests write with the CSPICE writers; `CONTRACT.md` records the internal interfaces.

| Area | Available | Not yet available |
|---|---|---|
| Kernel loading | SPK, CK, binary PCK, DSK, text kernels (leap seconds, PCK, frames, spacecraft clock), meta-kernels | EK |
| SPK data types | all types of CSPICE N0067: 1, 2, 3, 5, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21; `getelm`, `evsgp4` | |
| States | `spkezr`, `spkpos`, `spkgeo`, `spkssb`, `spkltc`, `spkaps`, `spkacs`, `spkapo`, with all aberration corrections of CSPICE | constant-position and constant-velocity observers (`spkcpo` and related calls) |
| Frames | inertial, body-fixed (text and binary PCK), fixed-offset, CK-based, dynamic (all families) and switch frames; `pxform`, `sxform`; `frinfo`, `cidfrm`, `cnmfrm`, `ccifrm`, `tkfram`, `tipbod`, `tisbod`, `pckfrm`, `pckcov` | |
| Derivatives | `spkder` and `pxder` to any order; a `ForwardDiff.Dual` epoch in the public calls | |
| Orientation and clocks | CK types 1 to 6, `ckgp`, `ckgpav`, `ckcov`, `ckobj`; spacecraft clock type 1 | |
| Time | `str2et`, `utc2et`, `timout`, `et2utc`, `etcal`, `unitim`, `tparse`, `tpictr`, `timdef`, `deltet`, `et2lst`, `lspcn` | `tsetyr` |
| File access and writers | DAF read and write calls, SPK writers for 15 types, CK writers for 5 types, comment area, windows, `spkobj`, `spkcov`, `kinfo`, `kdata` | type 14 SPK writer |
| Kernel pool | `gdpool`, `gcpool`, `gipool`, `bodvrd`, `bodvcd`, `bodfnd`, body and frame name calls; `pdpool`, `pcpool`, `pipool`, `dvpool`, `lmpool`, `ldpool`, `clpool`, `boddef`; `expool`, `dtpool`, `gnpool`, `szpool`, `stpool`; watchers `swpool`, `cvpool` for each context | |
| Mathematics | vectors, matrices, rotations and quaternions, coordinate conversions with their Jacobians, ellipsoid, plane and ellipse calls, conic elements and two-body propagation | |
| Shape kernels | DAS file access, DSK type 2 read and write, plate access calls, ray-surface intercept (`dskxv`, `dskxsi`, `dskx02`), `latsrf`, `srfnrm`, surface name calls | DSK types other than 2 |
| Geometry finder | `gfdist`, `gfsep`, `gfposc`, `gfrr`, `gfpa`, `gfuds`, `gfudb`, `gfevnt` for these quantities, `gfstol`, `gfrefn`, `uddf`, `uddc` | `gfsubc`, `gfsntc`, `gfilum`, `gfoclt`, `gftfov`, `gfrfov` |
| Geometry | | sub-points, illumination, occultation, limb and terminator points, instrument fields of view |

Three properties differ from CSPICE by design. All loaded data belong to a `SpiceContext`, whose readers take no lock, so that queries from many threads are safe while another thread loads or unloads kernels. A file that is open for access is an object that the caller owns (`h = dafopr(path)`), and the calls that use a hidden current file in CSPICE take that object as their first argument. An answer never depends on earlier calls, because the package keeps no buffers or saved search state between queries.

## Relation to NAIF

This package is an independent work. NAIF does not develop it, support it or endorse it, and the CSPICE toolkit of NAIF is the reference against which the package is tested.
