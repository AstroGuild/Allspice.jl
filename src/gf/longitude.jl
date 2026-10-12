# Geometry finder, coordinate searches: the solver for a coordinate of a vector
# (`zzgfcslv`) and the special treatment of longitude and right ascension (`zzgflong` of
# CSPICE N0067), whose values jump at a branch cut (at pi for longitude, at 0 and 2 pi
# for right ascension and for the longitudes that run from 0 to 2 pi).
#
# The jump is avoided with proxy functions. An equality is solved for the sine or the
# cosine of the coordinate, whichever changes faster at the reference value, and the
# roots of the wrong half of the circle are removed. An inequality or an extremum is
# solved region by region of the unit circle, with longitude as the function on the
# right half (where it has no jump) and right ascension on the left half.

"A coordinate with the vector of `q`, in the system `sys`, at the place `idx`; for the systems without radii."
_gf_proxy(q::GFCoordinate, sys::Int, idx::Int) =
    GFCoordinate(q.vec, sys, idx, 0.0, 0.0, 0, q.center)

"`smsgnd`: true if `x` and `y` are both positive or both negative."
@inline _gf_smsgnd(x::Float64, y::Float64) = (x > 0.0 && y > 0.0) || (x < 0.0 && y < 0.0)

"Progress message of pass `i` from the template `t` (`#` is the pass)."
_gf_pass(t::String, i::Int) = replace(t, "#" => string(i); count=1)

# `zzgfrelx` on a coordinate, on its cosine and on its sine.
_gf_relx_crd(q::GFCoordinate, relate, refval, cnfine::Window, run::_GFRun, prefix) =
    _gf_relx(et -> _gf_value(q, et), et -> _gf_decreasing(q, et), relate, refval, 0.0,
        cnfine, run, prefix, "done.")
_gf_relx_cos(q::GFCoordinate, relate, refval, cnfine::Window, run::_GFRun, prefix) =
    _gf_relx(et -> _gf_cosvalue(q, et), et -> _gf_cosdecreasing(q, et), relate, refval,
        0.0, cnfine, run, prefix, "done.")
_gf_relx_sin(q::GFCoordinate, relate, refval, cnfine::Window, run::_GFRun, prefix) =
    _gf_relx(et -> _gf_sinvalue(q, et), et -> _gf_sindecreasing(q, et), relate, refval,
        0.0, cnfine, run, prefix, "done.")

"""
    _gf_search(q::GFCoordinate, relate, refval, adjust, cnfine, run) -> Window

The coordinate solver `zzgfcslv`: the search is confined to the window on which the
vector exists, and a longitude or a right ascension goes to `_gf_long`.

Contract: `adjust >= 0` (`SpiceError(:valueoutofrange)`), a known relation
(`:notrecognized`); the result is a new window that is a subset of `cnfine`.
"""
function _gf_search(q::GFCoordinate, relate::AbstractString, refval::Real, adjust::Real,
        cnfine::Window, run::_GFRun)
    adjust < 0 && throw(SpiceError(:valueoutofrange,
        "ADJUST was $adjust; must be non-negative."))
    rel = _gf_relcode(relate)
    exist = _gf_vecexist(q.vec, cnfine, run)
    name = _GF_CRDNMS[q.sys][q.idx]
    (name == "LONGITUDE" || name == "RIGHT ASCENSION") &&
        return _gf_long(q, rel, Float64(refval), Float64(adjust), exist, run)
    return _gf_relx(et -> _gf_value(q, et), et -> _gf_decreasing(q, et), relate, refval,
        adjust, exist, run, _gf_prefix(q, rel, adjust), "done.")
end

"""
    _gf_long(q, rel, refval, adjust, cnfine, run) -> Window

The longitude and right-ascension solver `zzgflong`: window of `cnfine` on which the
angular coordinate `q` satisfies the relation with code `rel`.

The system is first mapped to one of two forms: latitudinal or geodetic longitude, in
(-pi, pi], or right ascension, in [0, 2 pi) (cylindrical and planetographic longitude
become right ascension; for a planetographic longitude that is positive west the
relation is mirrored: maximum and minimum change places, and `<`, `>` and `=` take the
reference value `2 pi - refval`). The proxy searches use the tolerance
`max(1e-7, tol / 10)`.

- `=`: roots of `sin(q) = sin(ref)` if `|sin(ref)| < sqrt(2)/2`, else of
  `cos(q) = cos(ref)`; a root is kept if the other function has the sign that it has at
  the reference value.
- `LOCMAX`, `LOCMIN`: the union of the local extrema of right ascension and of
  longitude, each a function without a jump on its half of the circle.
- `ABSMAX`, `ABSMIN`: the regions of the circle are searched from the one that holds
  the most extreme values, and the first region with an extremum gives the result. With
  `adjust > 0` the search continues as an inequality with the reference value
  `extremum -/+ adjust`, limited to the range of the coordinate.
- `<` (and `>`, as the complement): the region that holds the reference value is
  searched with a function without a jump there, and the regions that lie wholly below
  the reference value are added.

Contract: the checks of `_gf_checkadjust`. An inequality whose reference value is on
the branch cut is decided without a search: for a coordinate in [0, 2 pi) with
`cos(refval) == 1`, `<` gives the empty window and `>` gives `cnfine`; for a coordinate
in (-pi, pi] with `cos(refval) == -1`, `<` gives `cnfine` and `>` the empty window. The
result is a new window that is a subset of `cnfine`.
"""
function _gf_long(q::GFCoordinate, rel::Int, refval::Float64, adjust::Float64,
        cnfine::Window, run0::_GFRun)
    loctol = max(1.0e-7, run0.tol / 10.0)
    run = _gf_withtol(run0, loctol)
    _gf_checkadjust(rel, adjust)
    card(cnfine) == 0 && return Window()
    sys = q.sys
    if sys == _GF_CYLINDRICAL || sys == _GF_PLANETOGRAPHIC || sys == _GF_RADEC
        if cos(refval) == 1.0
            rel == _GF_LT && return Window()
            rel == _GF_GT && return copy(cnfine)
        end
    elseif cos(refval) == -1.0
        rel == _GF_LT && return copy(cnfine)
        rel == _GF_GT && return Window()
    end

    # Map the system to latitudinal (or geodetic) longitude or to right ascension.
    xrfval = refval
    uop = rel
    qn = q
    islon = true
    if sys == _GF_SPHERICAL
        qn = _gf_proxy(q, _GF_LATITUDINAL, 2)
    elseif sys == _GF_CYLINDRICAL
        qn = _gf_proxy(q, _GF_RADEC, 2)
        islon = false
    elseif sys == _GF_PLANETOGRAPHIC
        qn = _gf_proxy(q, _GF_RADEC, 2)
        islon = false
        if q.sense == -1
            # Positive west: the longitude is 2 pi minus the right ascension.
            if rel == _GF_ABSMAX
                uop = _GF_ABSMIN
            elseif rel == _GF_ABSMIN
                uop = _GF_ABSMAX
            elseif rel == _GF_LOCMAX
                uop = _GF_LOCMIN
            elseif rel == _GF_LOCMIN
                uop = _GF_LOCMAX
            elseif rel == _GF_EQ
                xrfval = twopi() - refval
            elseif rel == _GF_LT
                uop = _GF_GT
                xrfval = twopi() - refval
            else
                uop = _GF_LT
                xrfval = twopi() - refval
            end
        end
    elseif sys == _GF_RADEC
        islon = false
    end
    qlon = _gf_proxy(q, _GF_LATITUDINAL, 2)
    qra = _gf_proxy(q, _GF_RADEC, 2)
    r2ovr2 = sqrt(2.0) / 2.0

    if uop == _GF_EQ
        cv = cos(xrfval)
        sv = sin(xrfval)
        pre = ("Coordinate pass 1 of 2", "Coordinate pass 2 of 2")
        usesin = !(abs(sv) >= r2ovr2)
        if usesin
            f1 = _gf_relx_sin(qn, "=", sv, cnfine, run, pre)
        else
            f1 = _gf_relx_cos(qn, "=", cv, cnfine, run, pre)
        end
        cmpval = usesin ? cv : sv
        result = Window()
        # Remove the roots of the other half of the circle.
        for i in 1:2:card(f1)            # invariant: result holds the true roots before i
            start = f1[i]
            value = usesin ? _gf_cosvalue(qn, start) : _gf_sinvalue(qn, start)
            _gf_smsgnd(cmpval, value) && wninsd!(result, start, start)
        end
        return result
    end

    if uop == _GF_LOCMAX || uop == _GF_LOCMIN
        # The two halves of the circle, each made a little larger so that no root is
        # lost. In double precision `cos(1e-12)` is 1, so the left window is the window
        # on which the cosine is less than 1 and the right window is `cnfine`.
        left = _gf_relx_cos(qn, "<", cos(1.0e-12), cnfine, run,
            ("Coordinate pass 1 of 6", "Coordinate pass 2 of 6"))
        compl = _gf_relx_cos(qn, "<", -cos(1.0e-12), left, run,
            ("Coordinate pass 3 of 6", "Coordinate pass 4 of 6"))
        right = wndifd(cnfine, compl)
        op = _gf_relname(uop)
        res1 = _gf_relx_crd(qra, op, 0.0, left, run, ("Coordinate pass 5 of 6", ""))
        res2 = _gf_relx_crd(qlon, op, 0.0, right, run, ("Coordinate pass 6 of 6", ""))
        return wnunid(res1, res2)
    end

    # The remaining relations are ABSMAX, ABSMIN, '<' and '>'. First the windows on
    # which the coordinate is in the top, the bottom and the right half of the circle.
    isabs = uop == _GF_ABSMAX || uop == _GF_ABSMIN
    tmplat = isabs ? (adjust == 0.0 ? "Coordinate pass # of 7" : "Coordinate pass # of 7-9") :
             "Coordinate pass # of 8"
    top = _gf_relx_sin(qn, ">", sin(1.0e-11), cnfine, run,
        (_gf_pass(tmplat, 1), _gf_pass(tmplat, 2)))
    bot = _gf_relx_sin(qn, "<", -sin(1.0e-11), cnfine, run,
        (_gf_pass(tmplat, 3), _gf_pass(tmplat, 4)))
    right = _gf_relx_cos(qn, ">", 0.0, cnfine, run,
        (_gf_pass(tmplat, 5), _gf_pass(tmplat, 6)))
    total = 6
    result = Window()

    if isabs
        # The regions, from the one with the most extreme values. A region of the right
        # half is searched with longitude, a region of the left half with right
        # ascension.
        if islon
            q2 = wndifd(top, right)
            q3 = wndifd(bot, right)
            regions = uop == _GF_ABSMAX ? ((q2, false), (right, true), (q3, false)) :
                      ((q3, false), (right, true), (q2, false))
        else
            q1 = wnintd(right, top)
            q4 = wnintd(right, bot)
            left = wndifd(cnfine, right)
            regions = uop == _GF_ABSMAX ? ((q4, true), (left, false), (q1, true)) :
                      ((q1, true), (left, false), (q4, true))
        end
        op = _gf_relname(uop)
        found = false
        i = 1
        # invariant: the regions before `i` hold no extremum
        while i <= 3 && !found
            region, uselon = regions[i]
            qp = uselon ? qlon : qra
            # CSPICE passes the cosine of the coordinate as the value here, with the
            # test "is decreasing" of the coordinate itself (`zzgfcocg` with
            # `zzgfcodc`). The candidate epochs are the local extrema of the
            # coordinate, and among them the one with the most extreme cosine is
            # taken. This is kept, so that the window is the window of CSPICE.
            result = _gf_relx(et -> _gf_cosvalue(qp, et), et -> _gf_decreasing(qp, et),
                op, 0.0, 0.0, region, run, (_gf_pass(tmplat, total + 1), " "), "done.")
            card(region) > 0 && (total += 1)
            if wncard(result) > 0
                found = true
            else
                i += 1
            end
        end
        adjust == 0.0 && return result
        # CSPICE reads the first element of an empty result here; no extremum means no
        # window.
        card(result) == 0 && return result
        extval = _gf_value(qn, result[1])
        if uop == _GF_ABSMAX
            xrfval = islon ? max(extval - adjust, -Float64(pi)) : max(extval - adjust, 0.0)
            uop = _GF_GT
        else
            xrfval = islon ? min(extval + adjust, Float64(pi)) : min(extval + adjust, twopi())
            uop = _GF_LT
        end
    end

    # Inequality. '>' is solved as '<' and the result is complemented.
    flip = uop == _GF_GT
    cv = cos(xrfval)
    sv = sin(xrfval)
    locref = atan(sv, cv)
    quad = locref < -Float64(pi) / 2 ? 3 : locref < 0.0 ? 4 : locref < Float64(pi) / 2 ? 1 : 2
    pre = (_gf_pass(tmplat, total + 1), _gf_pass(tmplat, total + 2))
    if islon
        q2 = wndifd(top, right)
        q3 = wndifd(bot, right)
        # `s` is the region of the reference value and `f1` the part of the circle that
        # is wholly below the reference value.
        if quad == 1 || quad == 4
            s, uselon, f1 = right, true, q3
        elseif quad == 2
            s, uselon, f1 = q2, false, wnunid(q3, right)
        else
            s, uselon, f1 = q3, false, Window()
        end
        prxval = atan(sv, cv)
        (!uselon && prxval < 0.0) && (prxval += twopi())
        f2 = _gf_relx_crd(uselon ? qlon : qra, "<", prxval, s, run, pre)
        result = (quad == 3) ? f2 : wnunid(f1, f2)
        flip || return result
        return wndifd(wnunid(q3, wnunid(q2, right)), result)
    else
        q1 = wnintd(right, top)
        q4 = wnintd(right, bot)
        left = wndifd(cnfine, right)
        if quad == 1
            s, uselon, f1 = q1, true, Window()
        elseif quad == 2 || quad == 3
            s, uselon, f1 = left, false, q1
        else
            s, uselon, f1 = q4, true, wnunid(left, q1)
        end
        prxval = atan(sv, cv)
        (!uselon && prxval < 0.0) && (prxval += twopi())
        f2 = _gf_relx_crd(uselon ? qlon : qra, "<", prxval, s, run, pre)
        result = (quad == 1) ? f2 : wnunid(f1, f2)
        flip || return result
        return wndifd(wnunid(q4, wnunid(q1, left)), result)
    end
end
