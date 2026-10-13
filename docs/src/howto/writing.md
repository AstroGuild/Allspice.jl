# Writing kernels

The package writes the binary kernels that it reads, and the files are equal, byte for byte, to the files that the CSPICE writers produce from the same input. This page writes an SPK and reads it back; the CK, DSK and EK writers follow the same pattern of open, write, close.

## Write an SPK

[`spkopn`](@ref) makes a new file and returns a handle. Each `spkwNN` call adds one segment of data type `NN`, and [`spkcls`](@ref) closes the file. The segment below is of type 5, which holds states and interpolates between them with two-body motion; [`conics`](@ref) computes the states from orbital elements.

```@example write
using Allspice

mu = 3.986004418e5                           # km^3/s^2
first, last = 0.0, 86400.0                   # seconds past J2000 TDB
elements = [7000.0, 0.01, 0.9, 0.0, 0.0, 0.0, first, mu]
epochs = collect(range(first, last; length=97))
states = [conics(elements, t) for t in epochs]

path = joinpath(mktempdir(), "example.bsp")
handle = spkopn(path, "example", 0)
spkw05(handle, -77, 399, "J2000", first, last, "low orbit", mu, states, epochs)
spkcls(handle)
filesize(path)
```

The file changes on disk at each write call, as in CSPICE. The exception is the type 14 writer, which collects records in an object ([`spk14b`](@ref), [`spk14a`](@ref), [`spk14e`](@ref)) and writes the segment when it ends.

## Inspect a file

[`spkobj`](@ref) and [`spkcov`](@ref) read the bodies and the time coverage of a file that need not be loaded.

```@example write
(spkobj(path), spkcov(path, -77))
```

The DAF calls give access to the segments themselves. A handle is an object that the caller owns, and the search calls take it as an argument; CSPICE keeps one hidden "current" file instead.

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

```@example write
furnsh(path)
state, _ = spkezr("-77", 3600.0, "J2000", "NONE", "399")
kclear()
state[1:3]                                   # km
```

## Other writers

| Kernel | Open | Write | Close |
|---|---|---|---|
| SPK, 16 data types | [`spkopn`](@ref), [`spkopa`](@ref) | `spkw01` to `spkw21`, [`spk14b`](@ref) | [`spkcls`](@ref) |
| SPK from two-line elements | [`spkopn`](@ref) | [`getelm`](@ref), [`spkw10`](@ref) | [`spkcls`](@ref) |
| CK, 5 data types | [`ckopn`](@ref) | `ckw01`, `ckw02`, `ckw03`, `ckw05`, `ckw06` | [`ckcls`](@ref) |
| DSK type 2 | [`dskopn`](@ref) | [`dskmi2`](@ref), [`dskw02`](@ref) | [`dskcls`](@ref) |
| EK | [`ekopn`](@ref), [`ekops`](@ref) | [`ekbseg`](@ref), [`ekappr`](@ref), `ekace*`, or the fast writers | [`ekcls`](@ref) |
| Any DAF or DAS | [`dafonw`](@ref), [`dasonw`](@ref) | `dafada`, `dasadd` and related calls | [`dafcls`](@ref), [`dascls`](@ref) |

A text kernel is plain text; write it with `write`, or put values straight into the kernel pool with [`pdpool`](@ref), [`pcpool`](@ref) and [`pipool`](@ref).
