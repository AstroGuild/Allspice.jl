# Writing kernels

The package writes SPK, CK, DSK and EK files, and each file is equal to the file that the CSPICE writer produces from the same input. This page writes an SPK and reads it back; the other writers follow the same pattern of open, write, close.

## Write an SPK

[`spkopn`](@ref) makes a new file and returns a handle; the file must not exist. Each `spkwNN` call adds one segment of data type `NN`, and [`spkcls`](@ref) closes the file. Load the file only after the close.

The segment below is of type 5, which holds states and interpolates between them with two-body motion, so that a few states for each revolution are enough. [`conics`](@ref) computes the states from orbital elements.

```@example write
using Allspice

mu = 3.986004418e5                           # km^3/s^2
t0, t1 = 0.0, 86400.0                        # seconds past J2000 TDB

# perifocal distance (km), eccentricity, inclination, longitude of the node, argument
# of periapsis, mean anomaly at the epoch (rad), epoch, gravitational parameter
elements = [7000.0, 0.01, 0.9, 0.0, 0.0, 0.0, t0, mu]
epochs = collect(range(t0, t1; length=97))   # one state each 15 minutes
states = [conics(elements, t) for t in epochs]

path = joinpath(mktempdir(), "example.bsp")
handle = spkopn(path, "example", 0)          # path, internal file name, comment characters
# handle, body, centre, frame, coverage, segment name, mu, states, epochs
spkw05(handle, -77, 399, "J2000", t0, t1, "low orbit", mu, states, epochs)
spkcls(handle)
spkobj(path)                                 # the bodies in the file
```

The file changes on disk at each write call, as in CSPICE. The exception is the type 14 writer, which collects records in an object ([`spk14b`](@ref), [`spk14a`](@ref), [`spk14e`](@ref)) and writes the segment when it ends.

## Inspect a file

[`spkobj`](@ref) and [`spkcov`](@ref) read the bodies and the time coverage of a file that need not be loaded.

```@example write
spkcov(path, -77)
```

An SPK is a DAF (double precision array file), which is the container format of SPK, CK and binary PCK files, and the DAF calls give access to the segments themselves. The summary of an SPK segment holds 2 doubles (the coverage) and 6 integers (body, centre, frame, data type and two addresses).

```@example write
daf = dafopr(path)
dafbfs(daf)                                  # begin a forward search
while daffna(daf)                            # move to the next segment
    doubles, integers = dafus(dafgs(daf), 2, 6)
    println(dafgn(daf), ": body ", integers[1], ", centre ", integers[2], ", type ", integers[4])
end
dafcls(daf)
```

## Use the file

A body without a name is given by its ID code in a string. The epoch is a number here, so no leap seconds kernel is needed.

```@example write
furnsh(path)
state, _ = spkezr("-77", 3600.0, "J2000", "NONE", "399")
kclear()
state[1:3]                                   # km
```

## Other writers

| Kernel | Open | Write | Close |
|---|---|---|---|
| SPK, 16 data types | [`spkopn`](@ref), or [`spkopa`](@ref) to add to a file | `spkw01` to `spkw21`, [`spk14b`](@ref) | [`spkcls`](@ref) |
| SPK from two-line elements | [`spkopn`](@ref) | [`getelm`](@ref), [`spkw10`](@ref) | [`spkcls`](@ref) |
| CK, 5 data types | [`ckopn`](@ref) | `ckw01`, `ckw02`, `ckw03`, `ckw05`, `ckw06` | [`ckcls`](@ref) |
| DSK type 2 | [`dskopn`](@ref) | [`dskmi2`](@ref), [`dskw02`](@ref) | [`dskcls`](@ref) |
| EK | [`ekopn`](@ref), [`ekops`](@ref) | [`ekbseg`](@ref), [`ekappr`](@ref), `ekace*`, or whole columns with [`ekifld`](@ref) | [`ekcls`](@ref) |
| Any DAF or DAS (the container of DSK and EK files) | [`dafonw`](@ref), [`dasonw`](@ref) | `dafada`, `dasadd` and related calls | [`dafcls`](@ref), [`dascls`](@ref) |

A text kernel is plain text; write it with `write`, or put values straight into the kernel pool with [`pdpool`](@ref), [`pcpool`](@ref) and [`pipool`](@ref).
