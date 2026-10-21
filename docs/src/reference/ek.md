# Event kernels

An entry is the value of one column in one row: one element in a scalar column, one or
more in an array column, or null in a column that allows nulls. A record is a row; the
names of CSPICE use both words. Segments, rows, columns
and elements count from 1 in every call of this page, where the C interface of CSPICE
counts from 0. Sibling calls can throw different kinds for one condition, such as a
wrong data type or an entry that was never written, because each kind is the one of the
CSPICE routine that the call ports; see [`SpiceError`](@ref).

## Reading a file

```@docs
EKColumnAttr
EKSegmentSummary
ekopr
ekcls
eknseg
ekssum
ekrcec
ekrced
ekrcei
```

## Writing

```@docs
EKFastLoad
ekopn
ekopw
ekops
ekbseg
ekappr
ekinsr
ekdelr
ekacec
ekaced
ekacei
ekucec
ekuced
ekucei
ekifld
ekaclc
ekacld
ekacli
ekffld
```

## Queries

```@docs
EKResult
eklef
ekuef
ekntab
ektnam
ekccnt
ekcii
ekpsel
ekfind
eknelt
ekgc
ekgd
ekgi
ekcolumn
```
