# The error type of the package and the registry of its kinds.

"""
    Allspice.ERROR_KINDS

The kinds of [`SpiceError`](@ref) that the package can throw, as a sorted tuple of
symbols.

Use it to check a kind that a handler tests: `:nocoverage in Allspice.ERROR_KINDS`. A
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
    :cellarraytoosmall, :celltoosmall, :ckunknowndatatype, :coarsegridoverflow,
    :columnstoredtwice, :columntoosmall, :competingframespec, :coordsysnotrec,
    :countmismatch, :counttoolarge, :coveragegap, :crossanglemissing, :dafbeggtend,
    :daffrnotfound, :dafinvalidaccess, :dafinvalidparams, :dafnegaddr, :dafnewconflict,
    :dafnoresv, :dafnosearch, :dafnosuchhandle, :dafwritefail, :dasinvalidaccess,
    :dasnosuchaddress, :dasnosuchhandle, :dasnotempty, :daswritefail, :degeneratecase,
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
    :missingvalue, :namesnotresolved, :noclass, :noconverg, :noconvergence, :nocoverage,
    :nocurrentarray, :nodetoofull, :nodsksegments, :noframe, :noframeconnect, :noframeinfo,
    :nofreenodes, :nointerval, :noloadedfiles, :nolsk, :nonconicmotion, :nonemptyentry,
    :nonemptytree, :nonpositiveaxis, :nonpositivemass, :nonpositiveradius,
    :nonpositivescale, :nonpositivevalue, :nonprintablechars, :nonprintingchars,
    :nopartition, :nopath, :nosegmentsfound, :noseparation, :nosuchfile, :nosuchhandle,
    :notanumber, :notarotation, :notcomputable, :notdisjoint, :notindexed, :notinpart,
    :notranslation, :notrecognized, :notsupported, :numcoeffsnotpos, :numericoverflow,
    :numpacketsnotpos, :numstatesnotpos, :objectstooclose, :orbitdecay, :outofrange,
    :outofroom, :partialframespec, :platelisttoosmall, :pointeroutofrange, :pointnotfound,
    :pointnotinsegment, :pointnotonsurface, :pointoffsurface, :pointonzaxis, :pointtoosmall,
    :ptrarraytoosmall, :querynotparsed, :radiioutoforder, :rayiszerovector, :recordnotfound,
    :recursiontoodeep, :recursiveloading, :refanglemissing, :refvectormissing,
    :segidtoolong, :segmentnotopen, :shapemissing, :shapenotsupported, :sizemismatch,
    :spkinvalidoption, :spktypenotsupp, :suborbital, :subpointnotfound, :syntaxerror,
    :tablenotloaded, :timesdontmatch, :timesoutoforder, :toofewpackets, :toofewplates,
    :toofewstates, :toofewvertices, :toomanyconjunctions, :toomanyhits, :toomanysurfaces,
    :typemismatch, :typenotsupported, :unallocatednode, :undefinedframe, :uninitialized,
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
    SpiceError(kind::Symbol, msg::String)

The exception that the package throws for a call that fails.

`kind` is a short symbol that names the failure, such as `:nocoverage`,
`:kernelvarnotfound` or `:badarg`; test it to handle one failure. `msg` is the text for
a person.

A kind is the short error message that CSPICE signals in the same call, in lower case
(`:kernelvarnotfound` for `SPICE(KERNELVARNOTFOUND)`). Two calls can therefore report
the same failure with different kinds, as in CSPICE: a file that does not exist is
`:nosuchfile` in `furnsh` and `:filenotfound` in `spkobj`. Test the kinds that the
docstring of the call lists. A few kinds are names of this package, where CSPICE has no
message or a less exact one: `:badarg` for an argument of the wrong size or type,
`:nocoverage` and `:nopath` for missing SPK data, `:nolsk` for missing leap seconds.

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
"""
struct SpiceError <: Exception
    kind::Symbol
    msg::String
    SpiceError(kind::Symbol, msg::AbstractString) = new(_check_kind(kind), msg)
end
Base.showerror(io::IO, e::SpiceError) = print(io, "SpiceError(", e.kind, "): ", e.msg)
