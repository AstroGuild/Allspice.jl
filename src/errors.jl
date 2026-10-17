# The error type of the package and the registry of its kinds.

"""
    Allspice.ERROR_KINDS

The kinds of [`SpiceError`](@ref) that the package can throw, as a sorted tuple of
symbols.

Use it to check a kind that a handler tests: `:spkinsuffdata in Allspice.ERROR_KINDS`. A
`SpiceError` with a kind that is not in the tuple cannot be constructed, so a comparison
with such a kind is never true.
"""
const ERROR_KINDS = (
    :arraysizemismatch, :arraytoosmall, :avaloutofrange, :axisunderflow, :badarg,
    :badarraysize, :badattribute, :badattributes, :badaxis, :badaxislength, :badaxisnumbers,
    :badboresightspec, :badboundary, :badcoarsevoxscale, :badcolumdecl, :badcolumndecl,
    :badcommentarea, :badcoordbounds, :badcoordsys, :baddaf, :baddascommentarea,
    :baddasdirectory, :baddescrtimes, :baddimensions, :baddladescriptor, :badeccentricity,
    :badendpoints, :badepochtype, :badfile, :badfinevoxelscale, :badformat, :badframe,
    :badframecount, :badgeometry, :badgm, :badindex, :badinitstate, :badkernel,
    :badlatitudebounds, :badlatituderange, :badlatusrectum, :badlimblocusmix,
    :badlongituderange, :badmeanmotion, :badmeccentricity, :badmethodsyntax, :badmsemimajor,
    :badnumber, :badorbit, :badpartnumber, :badpeccentricity, :badperiapsevalue,
    :badplatecount, :badpool, :badpriorityspec, :badquatsign, :badradius, :badradiuscount,
    :badrecordcount, :badrefvectorspec, :badsemiaxis, :badsemilatus, :badsourceradius,
    :badstoptime, :badsurfacemap, :badtermlocusmix, :badtime, :badtimebounds, :badtimecase,
    :badtimedefault, :badtimetype, :badtkframe, :badtle, :badtypeshapecombo, :badunits,
    :badvariablesize, :badvariabletype, :badvarname, :badvector, :badvertexcount,
    :badvertexindex, :badwindowsize, :barraytoosmall, :barycenterephem, :blankfilename,
    :blankfiletype, :blanknameassigned, :bodiesnotdistinct, :bodyandcentersame,
    :boresightmissing, :boundarymissing, :boundarytoobig, :boundsdisagree,
    :boundsoutoforder, :buffertoosmall, :bug, :cantfindframe, :cantgetrotationtype,
    :cellarraytoosmall, :celltoosmall, :cknonexistrec, :ckunknowndatatype, :ckwrongdatatype,
    :coarsegridoverflow,
    :columnstoredtwice, :columntoosmall, :competingframespec, :coordsysnotrec,
    :countmismatch, :counttoolarge, :coveragegap, :crossanglemissing, :dafbeggtend,
    :daffrnotfound, :dafinvalidaccess, :dafinvalidparams, :dafnegaddr, :dafnewconflict,
    :dafnoresv, :dafnosearch, :dafnosuchhandle, :dafwritefail, :dasfilereadfailed,
    :dasinvalidaccess, :dasnosuchaddress, :dasnosuchhandle, :dasnotempty, :daswritefail,
    :degeneratecase,
    :degenerateinterval, :degeneratesurface, :degreeoutofrange, :dependentvectors,
    :difflinetoolarge, :difflinetoosmall, :dividebyzero, :dskdatanotfound, :dtoutofrange,
    :ekcolnummismatch, :ekmissingcolumn, :eknosegments, :emptysegment, :emptystring,
    :endpointsmatch, :facenotfound, :filarchmismatch, :fileisnotspk, :filenotfound,
    :fileopenfailed, :filereadfailed, :fovtoowide, :framechain, :framedatanotfound,
    :framedeferror, :frameidnotfound, :frameinfonotfound, :framemissing, :framenamenotfound,
    :ftpxfererror, :gridtoolarge, :idcodenotfound, :idstringtoolong, :illegalcharacter,
    :immutablevalue, :incompatiblescale, :incompatibleunits, :indexoutofrange,
    :inputoutofbounds, :inputstoolarge, :insidebody, :insufficientangles, :intindextoosmall,
    :intlennotpos, :invalidaddress, :invalidangle, :invalidarchtype, :invalidargument,
    :invalidaxis, :invalidaxislength, :invalidconststep, :invalidcount, :invaliddatatype,
    :invaliddegree, :invaliddescrtime, :invalidelements, :invalidellipse, :invalidfiletype,
    :invalidfixref, :invalidformat, :invalidfov, :invalidframe, :invalidgeometry,
    :invalidindex, :invalidlimbtype, :invalidlocus, :invalidlonextent, :invalidmethod,
    :invalidname, :invalidnode, :invalidnumint, :invalidnumints, :invalidnumrec,
    :invalidocctype, :invalidoption, :invalidplane, :invalidradii, :invalidradius,
    :invalidrefframe, :invalidrollstep, :invalidsclkrate, :invalidsclkstring,
    :invalidsclktime, :invalidsearchstep, :invalidselection, :invalidshadow, :invalidshape,
    :invalidshapecombo, :invalidsize, :invalidstarttime, :invalidstate, :invalidstep,
    :invalidstepsize, :invalidsublist, :invalidsubtype, :invalidtarget, :invalidtermtype,
    :invalidtolerance, :invalidtype, :invalidvalue, :invalidvertex, :irfnotrec,
    :itemnotfound, :kernelvarnotfound, :lowerboundtoolow, :missingdata, :missingframevar,
    :missingtimeinfo, :missingvalue, :namesnotresolved, :noclass, :noconverg,
    :noconvergence, :nocurrentarray, :nodetoofull, :nodsksegments, :noframe,
    :noframeconnect, :noframeinfo, :nofreenodes, :nointerval, :noleapseconds,
    :noloadedfiles, :nonconicmotion, :nonemptyentry, :nonemptytree, :nonpositiveaxis,
    :nonpositivemass, :nonpositiveradius, :nonpositivescale, :nonpositivevalue,
    :nonprintablechars, :nonprintingchars, :nopartition, :nosegmentsfound, :noseparation,
    :nosuchfile, :nosuchhandle, :notanumber, :notarotation, :notcomputable, :notdisjoint,
    :notindexed, :notinpart, :notranslation, :notrecognized, :notsupported,
    :numcoeffsnotpos, :numericoverflow, :numpacketsnotpos, :numstatesnotpos,
    :objectstooclose, :orbitdecay, :outofrange, :outofroom, :partialframespec,
    :platelisttoosmall, :pointeroutofrange, :pointnotfound, :pointnotinsegment,
    :pointnotonsurface, :pointoffsurface, :pointonzaxis, :pointtoosmall, :ptrarraytoosmall,
    :querynotparsed, :radiioutoforder, :rayiszerovector, :recordnotfound, :recursiontoodeep,
    :recursiveloading, :refanglemissing, :refvectormissing, :segidtoolong, :segmentnotopen,
    :shapemissing, :shapenotsupported, :sizemismatch, :spkinsuffdata, :spkinvalidoption,
    :spktypenotsupp, :suborbital, :subpointnotfound, :syntaxerror, :tablenotloaded,
    :timesdontmatch, :timesoutoforder, :toofewpackets, :toofewplates, :toofewstates,
    :toofewvertices, :toomanyconjunctions, :toomanyhits, :toomanysurfaces, :typemismatch,
    :typenotsupported, :unallocatednode, :undefinedframe, :uninitialized,
    :uninitializedvalue, :unitsmissing, :unitsnotrec, :unknownbff, :unknownckmeta,
    :unknownfilarc, :unknownframe, :unknownframetype, :unknownmode, :unknownspktype,
    :unknownsystem, :unmatchendpts, :unnaturalrelation, :unorderedrefs, :unorderedtimes,
    :unrecognizedtype, :unresolvednames, :unresolvedtimes, :unsupportedbff,
    :unsupportedspec, :valueoutofrange, :varnametoolong, :vectortoobig, :vertexnotingrid,
    :workspacetoosmall, :wrongarchitecture, :wrongdatatype, :wrongsegment, :zeroboresight,
    :zeroboundsextent, :zerolengthcolumn, :zeroposition, :zeroquaternion, :zerostep,
    :zerovector, :zerovelocity,
)

# The check is not inlined, so that a throw site in a numeric kernel stays one call.
@noinline function _check_kind(kind::Symbol)
    kind in ERROR_KINDS || throw(ArgumentError("$(repr(kind)) is not in Allspice.ERROR_KINDS"))
    return kind
end

"""
    SpiceError <: Exception
    SpiceError(kind::Symbol, msg) -> SpiceKindError

The abstract type of the exceptions that the package throws for a call that fails.

Every subtype has the fields `kind` and `msg`. `kind` is a short symbol that names the
failure, such as `:spkinsuffdata`, `:kernelvarnotfound` or `:badarg`; test it to handle
one failure. `msg` is the text for a person. The constructor returns a
[`SpiceKindError`](@ref), which has these two fields only.

A kind is the short error message that CSPICE signals in the same call, in lower case
(`:kernelvarnotfound` for `SPICE(KERNELVARNOTFOUND)`). Two calls can therefore report
the same failure with different kinds, as in CSPICE: a file that does not exist is
`:nosuchfile` in `furnsh` and `:filenotfound` in `spkobj`. Test the kinds that the
docstring of the call lists. Where CSPICE has no message, because the failure cannot
occur in C or the call has no counterpart, the kind is a name of this package; the
main one is `:badarg`, for an argument of the wrong size or type.

A call of CSPICE that has a `found` flag returns `nothing` here when nothing is found,
unless its docstring names another result (an error, `0`, an empty string or vector, or
the flag itself).

[`Allspice.ERROR_KINDS`](@ref) lists every kind. A kind that is not in the list cannot
be constructed.

# Examples
```jldoctest
julia> try
           gdpool(SpiceContext(), "NO_SUCH_VARIABLE")
       catch err
           err isa SpiceError ? err.kind : rethrow()
       end
:kernelvarnotfound
```

See also [`SpkInsufficientDataError`](@ref).
"""
abstract type SpiceError <: Exception end

"""
    SpiceKindError <: SpiceError

A [`SpiceError`](@ref) that carries no data but its `kind` and its `msg`.

Most failures throw this type. Catch `SpiceError` and test `kind`; the name of this type
is needed only to tell it from a subtype with more fields.
"""
struct SpiceKindError <: SpiceError
    kind::Symbol
    msg::String
    SpiceKindError(kind::Symbol, msg::AbstractString) = new(_check_kind(kind), msg)
end
SpiceError(kind::Symbol, msg::AbstractString) = SpiceKindError(kind, msg)

"""
    SpkInsufficientDataError <: SpiceError

The [`SpiceError`](@ref) of a state query that the loaded SPK data cannot answer.

[`spkezr`](@ref) and every call that needs the position of a body throw it.

| Field | Meaning |
|---|---|
| `kind` | `:spkinsuffdata`, or `:noloadedfiles` if no SPK file is loaded |
| `msg` | text for a person; it lists the coverage of `body` |
| `target`, `observer` | NAIF ID codes of the two bodies of the query |
| `et` | epoch of the missing data, seconds past J2000 TDB |
| `body` | ID code of the body whose data are missing, or `nothing` if it is not known |
| `loaded` | `true` if the loaded files hold `body` at other epochs than `et` |

`loaded` tells an epoch outside the data (`true`) from a body that is in no loaded file
(`false`). With an aberration correction, `et` is the epoch of the lookup that failed,
which differs from the epoch of the call by the light time. `body` is `nothing` with no
SPK file loaded. The ephemeris tables
throw a [`SpiceKindError`](@ref) with the kind `:spkinsuffdata` for an epoch outside
the table.

# Examples
```jldoctest
julia> err = try
           spkezr(SpiceContext(), 301, 0.0, "J2000", "NONE", 399)
       catch e
           e
       end;

julia> err isa SpkInsufficientDataError
true

julia> (err.kind, err.target, err.observer, err.et, err.body)
(:noloadedfiles, 301, 399, 0.0, nothing)
```
"""
struct SpkInsufficientDataError <: SpiceError
    kind::Symbol
    msg::String
    target::Int
    observer::Int
    et::Float64
    body::Union{Nothing,Int}
    loaded::Bool
    function SpkInsufficientDataError(kind::Symbol, msg::AbstractString, target::Integer,
            observer::Integer, et::Real, body::Union{Nothing,Integer}, loaded::Bool)
        @assert kind === :spkinsuffdata || kind === :noloadedfiles
        @assert !loaded || body !== nothing
        return new(_check_kind(kind), msg, target, observer, et, body, loaded)
    end
end

Base.showerror(io::IO, e::SpiceError) = print(io, "SpiceError(", e.kind, "): ", e.msg)
Base.showerror(io::IO, e::SpkInsufficientDataError) =
    print(io, "SpkInsufficientDataError(", e.kind, "): ", e.msg)
