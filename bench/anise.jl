# Minimal binding of the C interface of ANISE (through `anise_jll`), for the benchmarks.
# The header is `anise_jll.anise_h`.

using anise_jll

const libanise = anise_jll.libanise

"`AniseState` of anise.h."
struct AniseState
    position_km::NTuple{3,Cdouble}
    velocity_km_s::NTuple{3,Cdouble}
    epoch_et_s::Cdouble
    frame_ephemeris_id::Int32
    frame_orientation_id::Int32
end
AniseState() = AniseState((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0, 0, 0)

"Last error message of ANISE, or an empty string."
function anise_last_error()
    p = ccall((:anise_last_error, libanise), Cstring, ())
    return p == C_NULL ? "" : unsafe_string(p)
end

"Version string of the ANISE library."
anise_version() = unsafe_string(ccall((:anise_version, libanise), Cstring, ()))

"An ANISE almanac. Call `anise_free` when it is no longer needed."
mutable struct AniseAlmanac
    ptr::Ptr{Cvoid}
end

"New almanac with the files `paths` loaded in order. Throws if a load fails."
function anise_almanac(paths...)
    ptr = ccall((:anise_almanac_new, libanise), Ptr{Cvoid}, ())
    ptr == C_NULL && error("anise_almanac_new failed: " * anise_last_error())
    ref = Ref(ptr)
    for p in paths
        rc = ccall((:anise_almanac_load, libanise), Cint, (Ptr{Ptr{Cvoid}}, Cstring), ref, p)
        rc == 0 || error("anise_almanac_load($p) failed with code $rc: " * anise_last_error())
    end
    @assert ref[] != C_NULL
    return AniseAlmanac(ref[])
end

function anise_free(alm::AniseAlmanac)
    alm.ptr == C_NULL && return nothing
    ccall((:anise_almanac_free, libanise), Cvoid, (Ptr{Cvoid},), alm.ptr)
    alm.ptr = C_NULL
    return nothing
end

"""
    anise_spk_ezr(alm, target, et, frame, observer, out) -> Cint

State of `target` relative to `observer` in the orientation `frame` (1 is J2000) at ET
`et`, no aberration correction; the state is written to `out[]`. Returns the status
code of ANISE (0 is success). The call does not allocate.
"""
@inline function anise_spk_ezr(alm::AniseAlmanac, target::Integer, et::Float64,
        frame::Integer, observer::Integer, out::Base.RefValue{AniseState})
    return ccall((:anise_almanac_spk_ezr, libanise), Cint,
        (Ptr{Cvoid}, Int32, Cdouble, Int32, Cstring, Int32, Ptr{AniseState}),
        alm.ptr, target, et, frame, "NONE", observer, out)
end

"State as a 6-vector; throws with the ANISE message if the query fails."
function anise_state(alm::AniseAlmanac, target::Integer, et::Float64, frame::Integer,
        observer::Integer)
    out = Ref(AniseState())
    rc = anise_spk_ezr(alm, target, et, frame, observer, out)
    rc == 0 || error("anise_almanac_spk_ezr failed with code $rc: " * anise_last_error())
    s = out[]
    return SVector{6,Float64}(s.position_km..., s.velocity_km_s...)
end
