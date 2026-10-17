# Coverage and validation

## What is available

The package covers the library part of CSPICE N0067. The table lists each area with its calls; the reference pages document each call.

| Area | Available |
|---|---|
| Kernel loading | SPK, CK, binary PCK, DSK, EK, text kernels (leap seconds, PCK, frames, spacecraft clock), meta-kernels |
| SPK data types | all types of CSPICE N0067: 1, 2, 3, 5, 8, 9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21; `getelm`, `evsgp4` |
| States | `spkezr`, `spkpos`, `spkgeo`, `spkssb`, the integer-ID forms `spkez`, `spkezp` and `spkgps`, `spkltc`, `spkaps`, `spkacs`, `spkapo`, with all aberration corrections of CSPICE; constant-position and constant-velocity observers and targets (`spkcpo`, `spkcpt`, `spkcvo`, `spkcvt`) |
| Frames | inertial, body-fixed (text and binary PCK), fixed-offset, CK-based, dynamic (all families) and switch frames; `pxform`, `sxform`; `frinfo`, `cidfrm`, `cnmfrm`, `ccifrm`, `tkfram`, `tipbod`, `tisbod`, `pckfrm`, `pckcov`, `bltfrm`, `kplfrm` |
| Derivatives | `spkder` and `pxder` to any order; a `ForwardDiff.Dual` epoch in the public calls |
| Orientation and clocks | CK types 1 to 6, `ckgp`, `ckgpav`, `ckcov`, `ckobj`, `ckfrot`, `ckfxfm`; records of CK types 2 and 3 (`cknr02`, `ckgr02`, `cknr03`, `ckgr03`); spacecraft clock type 1 |
| Time | `str2et`, `utc2et`, `timout`, `et2utc`, `etcal`, `unitim`, `tparse`, `tpictr`, `timdef`, `tsetyr`, `tparch`, `deltet`, `et2lst`, `lspcn` |
| File access and writers | DAF read and write calls, SPK writers for all 16 types (type 14 through `spk14b`, `spk14a`, `spk14e`), CK writers for 5 types, the binary PCK writer of type 2 (`pckopn`, `pckw02`, `pckcls`), comment area, windows, `spkobj`, `spkcov`, `kinfo`, `kdata`; SPK segments by descriptor (`spksfs`, `spkpvn`, `spkuds`, and `spksub` for every type) |
| Kernel pool | `gdpool`, `gcpool`, `gipool`, `bodvrd`, `bodvcd`, `bodfnd`, body and frame name calls; `pdpool`, `pcpool`, `pipool`, `dvpool`, `lmpool`, `ldpool`, `clpool`, `boddef`; `expool`, `dtpool`, `gnpool`, `szpool`, `stpool`, `badkpv`; watchers `swpool`, `cvpool` for each context |
| Mathematics | vectors, matrices, rotations and quaternions, coordinate conversions with their Jacobians, ellipsoid, plane and ellipse calls, conic elements and two-body propagation, base 16 strings of doubles (`dp2hx`, `hx2dp`) |
| Shape kernels | DAS file access, DSK type 2 read and write, plate access calls, ray-surface intercept (`dskxv`, `dskxsi`, `dskx02`), `latsrf`, `srfnrm`, surface name calls (type 2 is the only DSK type of CSPICE N0067) |
| Event kernels | EK read access (`ekopr`, `ekssum`, `ekrcec`, `ekrced`, `ekrcei`), load and unload of one file (`eklef`, `ekuef`), table catalogue (`ekntab`, `ektnam`, `ekccnt`, `ekcii`), queries (`ekpsel`, `ekfind`, `eknelt`, `ekgc`, `ekgd`, `ekgi`) with a result object for each query, record and fast writers (`ekopn`, `ekops`, `ekbseg`, `ekappr`, `ekinsr`, `ekdelr`, `ekace*`, `ekuce*`, `ekifld`, `ekacl*`, `ekffld`) |
| Geometry finder | `gfdist`, `gfsep`, `gfposc`, `gfrr`, `gfpa`, `gfsubc`, `gfsntc`, `gfilum`, `gfoclt`, `gftfov`, `gfrfov`, `gfuds`, `gfudb`, `gfevnt` for the six quantities of CSPICE, `gfocce`, `gffove`, `gfstol`, `gfrefn`, `uddf`, `uddc` |
| Geometry | `subpnt`, `subslr`, `sincpt`, `ilumin`, `illumg`, `illumf`, `phaseq`, `occult`, `limbpt`, `termpt`, `tangpt`, `getfov`, `getfvn`, `fovtrg`, `fovray`, `azlcpo`, `azlrec`, `recazl`, `dazldr`, `drdazl`, `srfrec`, on ellipsoids and on DSK data |
| Ephemeris and rotation tables | `ephemeris_table`, `table_position`, `table_state`, `table_jet`: the position of one body relative to another on a fixed interval, as one array of coefficients for the CPU or a GPU; `table_to_spk`: the table as an SPK segment of type 2; `rotation_table`, `table_rotation`, `table_rotation_jet`: the same for the rotation between two frames |

## How the results were checked

CSPICE itself is the reference. The test suite calls each routine of this package and the routine of the same name in CSPICE, through SPICE.jl, on the same input, and compares the two results; it does not compare with stored numbers. The input is kernels from NAIF (planetary ephemerides, lunar orientation, leap seconds, constants) and kernels that the tests write with the CSPICE writers, so that every data type has a case.

| Kind of check | What it compares |
|---|---|
| Differential | A result of this package with the result of CSPICE for the same call |
| File | A kernel that this package wrote with the kernel that CSPICE wrote, byte for byte |
| Derivative | The first derivative with the velocity of CSPICE, and each higher derivative with automatic differentiation of the one below it |
| Error | The kind of the error, and the text of the message where a program can depend on it |
| Concurrency | Answers from many threads during changes of the kernel set with the answers of one thread |
| Allocation | A static proof (AllocCheck.jl) that 61 numeric routines cannot allocate |

| Measure | Value |
|---|---|
| Comparisons in the default run | about 80,000 |
| Comparisons in the exhaustive run (`ALLSPICE_TEST_FULL=1`) | about 100,000 |
| Line coverage of `src/` | 96.3 % |
| Julia versions tested | 1.10 and 1.13 |

## Known limits of accuracy

A result that differs from CSPICE differs by rounding. For the states of the Moon and the planets from de440s.bsp, in inertial and body-fixed frames, the default run measures the following.

| Quantity | Aberration corrections | Largest relative difference from CSPICE |
|---|---|---|
| Position, velocity and light time | all nine (`NONE` to `XCN+S`) | 7.3e-16 |

Some limits belong to the data and not to this package. The rotation angle of the Earth is a number of about 50,000 radians, whose rounding moves the Sun by some metres in the `IAU_EARTH` frame, in CSPICE and here. `spkezr` and `spkpos` differ from each other by up to 2e-11 relative with the converged corrections (`CN`, `CN+S`), in CSPICE as well as in this package, because the two routines stop their iterations by different rules.
