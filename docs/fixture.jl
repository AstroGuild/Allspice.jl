# Kernels for the examples of the manual. The repository holds no kernel file, so the
# examples run on three small kernels that this file writes with the package itself:
# a leap seconds kernel, a text kernel with the constants of the Earth and the name of a
# spacecraft, and an SPK with the Earth around the Sun and the spacecraft around the Earth.
# The motion is two-body motion, so the positions are not those of the real bodies.
module DocKernels

using Allspice

const MU_SUN = 1.32712440041e11      # km^3/s^2
const MU_EARTH = 3.986004418e5       # km^3/s^2

"The text of a leap seconds kernel with the leap seconds until 2017."
const LSK = """
\\begindata
DELTET/DELTA_T_A = 32.184
DELTET/K         = 1.657D-3
DELTET/EB        = 1.671D-2
DELTET/M         = ( 6.239996D0  1.99096871D-7 )
DELTET/DELTA_AT  = ( 10, @1972-JAN-1   11, @1972-JUL-1   12, @1973-JAN-1
                     13, @1974-JAN-1   14, @1975-JAN-1   15, @1976-JAN-1
                     16, @1977-JAN-1   17, @1978-JAN-1   18, @1979-JAN-1
                     19, @1980-JAN-1   20, @1981-JUL-1   21, @1982-JUL-1
                     22, @1983-JUL-1   23, @1985-JUL-1   24, @1988-JAN-1
                     25, @1990-JAN-1   26, @1991-JAN-1   27, @1992-JUL-1
                     28, @1993-JUL-1   29, @1994-JUL-1   30, @1996-JAN-1
                     31, @1997-JUL-1   32, @1999-JAN-1   33, @2006-JAN-1
                     34, @2009-JAN-1   35, @2012-JUL-1   36, @2015-JUL-1
                     37, @2017-JAN-1 )
\\begintext
"""

"The text of a kernel with the size and the rotation of the Sun and of the Earth, and the name of the spacecraft."
const CONSTANTS = """
\\begindata
BODY10_RADII     = ( 695700.0  695700.0  695700.0 )
BODY10_POLE_RA   = ( 286.13   0.0  0.0 )
BODY10_POLE_DEC  = (  63.87   0.0  0.0 )
BODY10_PM        = (  84.176  14.1844  0.0 )
BODY399_RADII    = ( 6378.1366  6378.1366  6356.7519 )
BODY399_POLE_RA  = (   0.0     -0.641   0.0 )
BODY399_POLE_DEC = (  90.0     -0.557   0.0 )
BODY399_PM       = ( 190.147  360.9856235  0.0 )
BODY399_GM       = 3.986004418D5
NAIF_BODY_NAME  += 'PROBE'
NAIF_BODY_CODE  += -999
\\begintext
"""

"""
    DocKernels.write(dir) -> Vector{String}

Write the three kernels into the directory `dir` and return their paths (leap seconds,
constants, SPK). The SPK covers the year 2025 (TDB) with two-body motion: the Earth on
a circle of one astronomical unit around the Sun in the ecliptic, and the spacecraft
`PROBE` (ID -999) on an orbit around the Earth with a perigee altitude of 700 km, an
apogee altitude of 35,786 km and an inclination of 28.5 degrees.
"""
function write(dir::AbstractString)
    lsk = joinpath(dir, "leapseconds.tls")
    pck = joinpath(dir, "constants.tpc")
    spk = joinpath(dir, "probe.bsp")
    Base.write(lsk, LSK)
    Base.write(pck, CONSTANTS)
    first = 7.889184e8                       # 2025 JAN 01 12:00 TDB
    last = first + 365 * 86400.0
    h = spkopn(spk, "manual examples", 0)
    # The Sun stays at the solar system barycentre, which aberration corrections need.
    spkw08(h, 10, 0, "J2000", first, last, "Sun at the barycentre", 1,
        [zeros(6) for _ in 1:2], first, last - first)
    earth = [1.495978707e8, 0.0, 0.0, 0.0, 100.46 * rpd(), 0.0, first, MU_SUN]
    epochs = collect(range(first, last; length=74))
    spkw05(h, 399, 10, "ECLIPJ2000", first, last, "Earth, two-body", MU_SUN,
        [conics(earth, t) for t in epochs], epochs)
    rp, ra = 6378.1366 + 700.0, 6378.1366 + 35786.0
    probe = [rp, (ra - rp) / (ra + rp), 28.5 * rpd(), 40.0 * rpd(), 270.0 * rpd(), 0.0,
        first, MU_EARTH]
    epochs = collect(range(first, last; length=4001))
    spkw05(h, -999, 399, "J2000", first, last, "PROBE, two-body", MU_EARTH,
        [conics(probe, t) for t in epochs], epochs)
    spkcls(h)
    return [lsk, pck, spk]
end

end # module
