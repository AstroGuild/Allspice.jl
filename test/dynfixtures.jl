# Frame kernels that the dynamic-frame and switch-frame groups write for themselves.
# No kernel with frames of class 5 or class 6 is on this machine, so the tests write one
# by hand that has every frame family and every kind of defining vector, and then load it
# with both engines.

module DynFix

const PATHS = Dict{Symbol,String}()
const DIR = Ref("")

"Directory of the generated kernels (made once per process)."
function dir()
    isdir(DIR[]) || (DIR[] = mktempdir())
    return DIR[]
end

"""
Frames of the kernel of `frame_kernel`, as `(name, family, center)`. `family` is the
label of the row in the table of largest differences.
"""
const FRAMES = [
    ("DT_GSE", "two-vector: position, velocity (GSE)", 399),
    ("DT_GSE_LTS", "two-vector: position LT+S, velocity LT+S", 399),
    ("DT_GSE_CN", "two-vector: position CN+S, velocity XCN+S", 399),
    ("DT_GSM", "two-vector: position, constant RA/DEC (GSM-like)", 399),
    ("DT_CONST_LAT", "two-vector: constant latitudinal, constant rectangular", 399),
    ("DT_CONST_ABCORR", "two-vector: position, constant with S", 399),
    ("DT_NADIR", "two-vector: near point, velocity", 301),
    ("DT_NADIR_LT", "two-vector: near point LT+S, position", 301),
    ("DT_VEL_PCK", "two-vector: velocity in PCK frame, position", 399),
    ("DT_EME_OF_DATE", "mean equator and equinox of date", 399),
    ("DT_ETE_OF_DATE", "true equator and equinox of date", 399),
    ("DT_ECL_OF_DATE", "mean ecliptic and equinox of date", 399),
    ("DT_EULER", "Euler", 399),
    ("DT_EULER_DEG", "Euler (degrees, base PCK frame)", 399),
    ("DT_GSE_FROZEN", "frozen two-vector", 399),
    ("DT_ETE_FROZEN", "frozen true of date", 399),
    ("DT_ETE_INERTIAL", "inertial rotation state (true of date)", 399),
    ("DT_GSE_INERTIAL", "inertial rotation state (two-vector)", 399),
    ("DT_EULER_ON_GSE", "dynamic relative to dynamic (Euler on GSE)", 399),
    ("DT_TWOVEC_IN_DYN", "dynamic relative to dynamic (vectors in dynamic frames)", 399),
    ("DT_PRODUCT", "product (three factors, two dynamic)", 399),
    ("DT_PRODUCT_ONE", "product (one factor, base PCK frame)", 399),
    ("DT_SWITCH_TIME", "switch frame, time intervals", 399),
    ("DT_SWITCH_ALL", "switch frame, no intervals", 399),
    ("DT_SWITCH_NEST", "switch frame over a switch frame", 399),
]

"ET bounds of the intervals of `DT_SWITCH_TIME` (ascending), for the tests at the edges."
const SWITCH_EDGES = [-6.0e8, -2.0e8, -1.0e8, 1.0e8, 2.5e8, 6.0e8]

"""
    frame_kernel() -> String

Write (once per process) the frame kernel of the dynamic and switch frames `FRAMES` and
return its path. Postcondition: the file exists.
"""
function frame_kernel()
    haskey(PATHS, :fk) && return PATHS[:fk]
    path = joinpath(dir(), "test_dynamic.tf")
    write(path, raw"""
    KPL/FK

    Dynamic frames (class 5) and switch frames (class 6) written by the tests.

    \begindata

       FRAME_DT_GSE                     = 1850001
       FRAME_1850001_NAME               = 'DT_GSE'
       FRAME_1850001_CLASS              = 5
       FRAME_1850001_CLASS_ID           = 1850001
       FRAME_1850001_CENTER             = 399
       FRAME_1850001_RELATIVE           = 'J2000'
       FRAME_1850001_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850001_FAMILY             = 'TWO-VECTOR'
       FRAME_1850001_PRI_AXIS           = 'X'
       FRAME_1850001_PRI_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
       FRAME_1850001_PRI_OBSERVER       = 'EARTH'
       FRAME_1850001_PRI_TARGET         = 'SUN'
       FRAME_1850001_PRI_ABCORR         = 'NONE'
       FRAME_1850001_SEC_AXIS           = 'Y'
       FRAME_1850001_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_VELOCITY'
       FRAME_1850001_SEC_OBSERVER       = 'SUN'
       FRAME_1850001_SEC_TARGET         = 'EARTH'
       FRAME_1850001_SEC_ABCORR         = 'NONE'
       FRAME_1850001_SEC_FRAME          = 'J2000'

       FRAME_DT_GSE_LTS                 = 1850002
       FRAME_1850002_NAME               = 'DT_GSE_LTS'
       FRAME_1850002_CLASS              = 5
       FRAME_1850002_CLASS_ID           = 1850002
       FRAME_1850002_CENTER             = 399
       FRAME_1850002_RELATIVE           = 'J2000'
       FRAME_1850002_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850002_FAMILY             = 'TWO-VECTOR'
       FRAME_1850002_PRI_AXIS           = 'X'
       FRAME_1850002_PRI_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
       FRAME_1850002_PRI_OBSERVER       = 'EARTH'
       FRAME_1850002_PRI_TARGET         = 'SUN'
       FRAME_1850002_PRI_ABCORR         = 'LT+S'
       FRAME_1850002_SEC_AXIS           = '-Y'
       FRAME_1850002_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_VELOCITY'
       FRAME_1850002_SEC_OBSERVER       = 'EARTH'
       FRAME_1850002_SEC_TARGET         = 'SUN'
       FRAME_1850002_SEC_ABCORR         = 'LT+S'
       FRAME_1850002_SEC_FRAME          = 'ECLIPJ2000'

       FRAME_DT_GSE_CN                  = 1850003
       FRAME_1850003_NAME               = 'DT_GSE_CN'
       FRAME_1850003_CLASS              = 5
       FRAME_1850003_CLASS_ID           = 1850003
       FRAME_1850003_CENTER             = 399
       FRAME_1850003_RELATIVE           = 'ECLIPJ2000'
       FRAME_1850003_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850003_FAMILY             = 'TWO-VECTOR'
       FRAME_1850003_PRI_AXIS           = '-Z'
       FRAME_1850003_PRI_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
       FRAME_1850003_PRI_OBSERVER       = 'EARTH'
       FRAME_1850003_PRI_TARGET         = 'MOON'
       FRAME_1850003_PRI_ABCORR         = 'CN+S'
       FRAME_1850003_SEC_AXIS           = 'X'
       FRAME_1850003_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_VELOCITY'
       FRAME_1850003_SEC_OBSERVER       = 'EARTH'
       FRAME_1850003_SEC_TARGET         = 'MARS BARYCENTER'
       FRAME_1850003_SEC_ABCORR         = 'XCN+S'
       FRAME_1850003_SEC_FRAME          = 'J2000'

       FRAME_DT_GSM                     = 1850004
       FRAME_1850004_NAME               = 'DT_GSM'
       FRAME_1850004_CLASS              = 5
       FRAME_1850004_CLASS_ID           = 1850004
       FRAME_1850004_CENTER             = 399
       FRAME_1850004_RELATIVE           = 'J2000'
       FRAME_1850004_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850004_FAMILY             = 'TWO-VECTOR'
       FRAME_1850004_PRI_AXIS           = 'X'
       FRAME_1850004_PRI_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
       FRAME_1850004_PRI_OBSERVER       = 'EARTH'
       FRAME_1850004_PRI_TARGET         = 'SUN'
       FRAME_1850004_PRI_ABCORR         = 'NONE'
       FRAME_1850004_SEC_AXIS           = 'Z'
       FRAME_1850004_SEC_VECTOR_DEF     = 'CONSTANT'
       FRAME_1850004_SEC_FRAME          = 'IAU_EARTH'
       FRAME_1850004_SEC_SPEC           = 'RA/DEC'
       FRAME_1850004_SEC_UNITS          = 'DEGREES'
       FRAME_1850004_SEC_RA             = 288.43
       FRAME_1850004_SEC_DEC            = 79.54

       FRAME_DT_CONST_LAT               = 1850005
       FRAME_1850005_NAME               = 'DT_CONST_LAT'
       FRAME_1850005_CLASS              = 5
       FRAME_1850005_CLASS_ID           = 1850005
       FRAME_1850005_CENTER             = 399
       FRAME_1850005_RELATIVE           = 'J2000'
       FRAME_1850005_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850005_FAMILY             = 'TWO-VECTOR'
       FRAME_1850005_PRI_AXIS           = 'Z'
       FRAME_1850005_PRI_VECTOR_DEF     = 'CONSTANT'
       FRAME_1850005_PRI_FRAME          = 'IAU_EARTH'
       FRAME_1850005_PRI_SPEC           = 'LATITUDINAL'
       FRAME_1850005_PRI_UNITS          = 'RADIANS'
       FRAME_1850005_PRI_LONGITUDE      = 0.7
       FRAME_1850005_PRI_LATITUDE       = -0.3
       FRAME_1850005_SEC_AXIS           = '-X'
       FRAME_1850005_SEC_VECTOR_DEF     = 'CONSTANT'
       FRAME_1850005_SEC_FRAME          = 'ECLIPB1950'
       FRAME_1850005_SEC_SPEC           = 'RECTANGULAR'
       FRAME_1850005_SEC_VECTOR         = ( 0.3, -1.2, 2.5 )
       FRAME_1850005_ANGLE_SEP_TOL      = 0.0005

       FRAME_DT_CONST_ABCORR            = 1850006
       FRAME_1850006_NAME               = 'DT_CONST_ABCORR'
       FRAME_1850006_CLASS              = 5
       FRAME_1850006_CLASS_ID           = 1850006
       FRAME_1850006_CENTER             = 399
       FRAME_1850006_RELATIVE           = 'J2000'
       FRAME_1850006_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850006_FAMILY             = 'TWO-VECTOR'
       FRAME_1850006_PRI_AXIS           = 'Y'
       FRAME_1850006_PRI_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
       FRAME_1850006_PRI_OBSERVER       = 'EARTH'
       FRAME_1850006_PRI_TARGET         = 'MOON'
       FRAME_1850006_PRI_ABCORR         = 'XLT'
       FRAME_1850006_SEC_AXIS           = 'Z'
       FRAME_1850006_SEC_VECTOR_DEF     = 'CONSTANT'
       FRAME_1850006_SEC_FRAME          = 'J2000'
       FRAME_1850006_SEC_SPEC           = 'RECTANGULAR'
       FRAME_1850006_SEC_VECTOR         = ( 0.1, 0.2, 1.0 )
       FRAME_1850006_SEC_OBSERVER       = 'EARTH'
       FRAME_1850006_SEC_ABCORR         = 'S'

       FRAME_DT_NADIR                   = 1850007
       FRAME_1850007_NAME               = 'DT_NADIR'
       FRAME_1850007_CLASS              = 5
       FRAME_1850007_CLASS_ID           = 1850007
       FRAME_1850007_CENTER             = 301
       FRAME_1850007_RELATIVE           = 'J2000'
       FRAME_1850007_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850007_FAMILY             = 'TWO-VECTOR'
       FRAME_1850007_PRI_AXIS           = 'Z'
       FRAME_1850007_PRI_VECTOR_DEF     = 'TARGET_NEAR_POINT'
       FRAME_1850007_PRI_OBSERVER       = 'MOON'
       FRAME_1850007_PRI_TARGET         = 'EARTH'
       FRAME_1850007_PRI_ABCORR         = 'NONE'
       FRAME_1850007_SEC_AXIS           = 'X'
       FRAME_1850007_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_VELOCITY'
       FRAME_1850007_SEC_OBSERVER       = 'EARTH'
       FRAME_1850007_SEC_TARGET         = 'MOON'
       FRAME_1850007_SEC_ABCORR         = 'NONE'
       FRAME_1850007_SEC_FRAME          = 'J2000'

       FRAME_DT_NADIR_LT                = 1850008
       FRAME_1850008_NAME               = 'DT_NADIR_LT'
       FRAME_1850008_CLASS              = 5
       FRAME_1850008_CLASS_ID           = 1850008
       FRAME_1850008_CENTER             = 301
       FRAME_1850008_RELATIVE           = 'J2000'
       FRAME_1850008_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850008_FAMILY             = 'TWO-VECTOR'
       FRAME_1850008_PRI_AXIS           = '-X'
       FRAME_1850008_PRI_VECTOR_DEF     = 'TARGET_NEAR_POINT'
       FRAME_1850008_PRI_OBSERVER       = 'MOON'
       FRAME_1850008_PRI_TARGET         = 'EARTH'
       FRAME_1850008_PRI_ABCORR         = 'LT+S'
       FRAME_1850008_SEC_AXIS           = 'Y'
       FRAME_1850008_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
       FRAME_1850008_SEC_OBSERVER       = 'MOON'
       FRAME_1850008_SEC_TARGET         = 'SUN'
       FRAME_1850008_SEC_ABCORR         = 'LT'

       FRAME_DT_VEL_PCK                 = 1850009
       FRAME_1850009_NAME               = 'DT_VEL_PCK'
       FRAME_1850009_CLASS              = 5
       FRAME_1850009_CLASS_ID           = 1850009
       FRAME_1850009_CENTER             = 399
       FRAME_1850009_RELATIVE           = 'J2000'
       FRAME_1850009_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850009_FAMILY             = 'TWO-VECTOR'
       FRAME_1850009_PRI_AXIS           = 'Y'
       FRAME_1850009_PRI_VECTOR_DEF     = 'OBSERVER_TARGET_VELOCITY'
       FRAME_1850009_PRI_OBSERVER       = 'EARTH'
       FRAME_1850009_PRI_TARGET         = 'MOON'
       FRAME_1850009_PRI_ABCORR         = 'LT'
       FRAME_1850009_PRI_FRAME          = 'IAU_EARTH'
       FRAME_1850009_SEC_AXIS           = 'Z'
       FRAME_1850009_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
       FRAME_1850009_SEC_OBSERVER       = 'EARTH'
       FRAME_1850009_SEC_TARGET         = 'MOON'
       FRAME_1850009_SEC_ABCORR         = 'NONE'

       FRAME_DT_EME_OF_DATE             = 1850010
       FRAME_1850010_NAME               = 'DT_EME_OF_DATE'
       FRAME_1850010_CLASS              = 5
       FRAME_1850010_CLASS_ID           = 1850010
       FRAME_1850010_CENTER             = 399
       FRAME_1850010_RELATIVE           = 'J2000'
       FRAME_1850010_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850010_FAMILY             = 'MEAN_EQUATOR_AND_EQUINOX_OF_DATE'
       FRAME_1850010_PREC_MODEL         = 'EARTH_IAU_1976'
       FRAME_1850010_ROTATION_STATE     = 'ROTATING'

       FRAME_DT_ETE_OF_DATE             = 1850011
       FRAME_1850011_NAME               = 'DT_ETE_OF_DATE'
       FRAME_1850011_CLASS              = 5
       FRAME_1850011_CLASS_ID           = 1850011
       FRAME_1850011_CENTER             = 399
       FRAME_1850011_RELATIVE           = 'J2000'
       FRAME_1850011_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850011_FAMILY             = 'TRUE_EQUATOR_AND_EQUINOX_OF_DATE'
       FRAME_1850011_PREC_MODEL         = 'EARTH_IAU_1976'
       FRAME_1850011_NUT_MODEL          = 'EARTH_IAU_1980'
       FRAME_1850011_ROTATION_STATE     = 'ROTATING'

       FRAME_DT_ECL_OF_DATE             = 1850012
       FRAME_1850012_NAME               = 'DT_ECL_OF_DATE'
       FRAME_1850012_CLASS              = 5
       FRAME_1850012_CLASS_ID           = 1850012
       FRAME_1850012_CENTER             = 399
       FRAME_1850012_RELATIVE           = 'J2000'
       FRAME_1850012_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850012_FAMILY             = 'MEAN_ECLIPTIC_AND_EQUINOX_OF_DATE'
       FRAME_1850012_PREC_MODEL         = 'EARTH_IAU_1976'
       FRAME_1850012_OBLIQ_MODEL        = 'EARTH_IAU_1980'
       FRAME_1850012_ROTATION_STATE     = 'ROTATING'

       FRAME_DT_EULER                   = 1850013
       FRAME_1850013_NAME               = 'DT_EULER'
       FRAME_1850013_CLASS              = 5
       FRAME_1850013_CLASS_ID           = 1850013
       FRAME_1850013_CENTER             = 399
       FRAME_1850013_RELATIVE           = 'J2000'
       FRAME_1850013_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850013_FAMILY             = 'EULER'
       FRAME_1850013_EPOCH              = @2010-JAN-01/12:00:00
       FRAME_1850013_AXES               = ( 3, 1, 3 )
       FRAME_1850013_UNITS              = 'RADIANS'
       FRAME_1850013_ANGLE_1_COEFFS     = ( 0.3, 7.0E-5, -2.0E-15 )
       FRAME_1850013_ANGLE_2_COEFFS     = ( -0.41 )
       FRAME_1850013_ANGLE_3_COEFFS     = ( 1.1, -3.0E-9, 4.0E-18, 1.0E-28 )

       FRAME_DT_EULER_DEG               = 1850014
       FRAME_1850014_NAME               = 'DT_EULER_DEG'
       FRAME_1850014_CLASS              = 5
       FRAME_1850014_CLASS_ID           = 1850014
       FRAME_1850014_CENTER             = 399
       FRAME_1850014_RELATIVE           = 'IAU_EARTH'
       FRAME_1850014_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850014_FAMILY             = 'EULER'
       FRAME_1850014_EPOCH              = @2000-JAN-01/12:00:00
       FRAME_1850014_AXES               = ( 1, 2, 3 )
       FRAME_1850014_UNITS              = 'DEGREES'
       FRAME_1850014_ANGLE_1_COEFFS     = ( 10.0, 1.0E-4 )
       FRAME_1850014_ANGLE_2_COEFFS     = ( -25.0, 2.0E-6, 1.0E-14 )
       FRAME_1850014_ANGLE_3_COEFFS     = ( 0.0, 4.1780745823E-3 )

       FRAME_DT_GSE_FROZEN              = 1850015
       FRAME_1850015_NAME               = 'DT_GSE_FROZEN'
       FRAME_1850015_CLASS              = 5
       FRAME_1850015_CLASS_ID           = 1850015
       FRAME_1850015_CENTER             = 399
       FRAME_1850015_RELATIVE           = 'J2000'
       FRAME_1850015_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850015_FAMILY             = 'TWO-VECTOR'
       FRAME_1850015_FREEZE_EPOCH       = @2015-MAR-03/01:02:03
       FRAME_1850015_PRI_AXIS           = 'X'
       FRAME_1850015_PRI_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
       FRAME_1850015_PRI_OBSERVER       = 'EARTH'
       FRAME_1850015_PRI_TARGET         = 'SUN'
       FRAME_1850015_PRI_ABCORR         = 'NONE'
       FRAME_1850015_SEC_AXIS           = 'Y'
       FRAME_1850015_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_VELOCITY'
       FRAME_1850015_SEC_OBSERVER       = 'SUN'
       FRAME_1850015_SEC_TARGET         = 'EARTH'
       FRAME_1850015_SEC_ABCORR         = 'NONE'
       FRAME_1850015_SEC_FRAME          = 'J2000'

       FRAME_DT_ETE_FROZEN              = 1850016
       FRAME_1850016_NAME               = 'DT_ETE_FROZEN'
       FRAME_1850016_CLASS              = 5
       FRAME_1850016_CLASS_ID           = 1850016
       FRAME_1850016_CENTER             = 399
       FRAME_1850016_RELATIVE           = 'J2000'
       FRAME_1850016_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850016_FAMILY             = 'TRUE_EQUATOR_AND_EQUINOX_OF_DATE'
       FRAME_1850016_PREC_MODEL         = 'EARTH_IAU_1976'
       FRAME_1850016_NUT_MODEL          = 'EARTH_IAU_1980'
       FRAME_1850016_FREEZE_EPOCH       = @1995-JUN-15/00:00:00

       FRAME_DT_ETE_INERTIAL            = 1850017
       FRAME_1850017_NAME               = 'DT_ETE_INERTIAL'
       FRAME_1850017_CLASS              = 5
       FRAME_1850017_CLASS_ID           = 1850017
       FRAME_1850017_CENTER             = 399
       FRAME_1850017_RELATIVE           = 'J2000'
       FRAME_1850017_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850017_FAMILY             = 'TRUE_EQUATOR_AND_EQUINOX_OF_DATE'
       FRAME_1850017_PREC_MODEL         = 'EARTH_IAU_1976'
       FRAME_1850017_NUT_MODEL          = 'EARTH_IAU_1980'
       FRAME_1850017_ROTATION_STATE     = 'INERTIAL'

       FRAME_DT_GSE_INERTIAL            = 1850018
       FRAME_1850018_NAME               = 'DT_GSE_INERTIAL'
       FRAME_1850018_CLASS              = 5
       FRAME_1850018_CLASS_ID           = 1850018
       FRAME_1850018_CENTER             = 399
       FRAME_1850018_RELATIVE           = 'J2000'
       FRAME_1850018_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850018_FAMILY             = 'TWO-VECTOR'
       FRAME_1850018_ROTATION_STATE     = 'INERTIAL'
       FRAME_1850018_PRI_AXIS           = 'X'
       FRAME_1850018_PRI_VECTOR_DEF     = 'OBSERVER_TARGET_POSITION'
       FRAME_1850018_PRI_OBSERVER       = 'EARTH'
       FRAME_1850018_PRI_TARGET         = 'SUN'
       FRAME_1850018_PRI_ABCORR         = 'NONE'
       FRAME_1850018_SEC_AXIS           = 'Y'
       FRAME_1850018_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_VELOCITY'
       FRAME_1850018_SEC_OBSERVER       = 'SUN'
       FRAME_1850018_SEC_TARGET         = 'EARTH'
       FRAME_1850018_SEC_ABCORR         = 'NONE'
       FRAME_1850018_SEC_FRAME          = 'J2000'

       FRAME_DT_EULER_ON_GSE            = 1850019
       FRAME_1850019_NAME               = 'DT_EULER_ON_GSE'
       FRAME_1850019_CLASS              = 5
       FRAME_1850019_CLASS_ID           = 1850019
       FRAME_1850019_CENTER             = 399
       FRAME_1850019_RELATIVE           = 'DT_GSE'
       FRAME_1850019_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850019_FAMILY             = 'EULER'
       FRAME_1850019_EPOCH              = @2005-JAN-01/00:00:00
       FRAME_1850019_AXES               = ( 3, 2, 1 )
       FRAME_1850019_UNITS              = 'DEGREES'
       FRAME_1850019_ANGLE_1_COEFFS     = ( 5.0, 1.0E-6 )
       FRAME_1850019_ANGLE_2_COEFFS     = ( 30.0 )
       FRAME_1850019_ANGLE_3_COEFFS     = ( -12.0, 0.0, 3.0E-16 )

       FRAME_DT_TWOVEC_IN_DYN           = 1850020
       FRAME_1850020_NAME               = 'DT_TWOVEC_IN_DYN'
       FRAME_1850020_CLASS              = 5
       FRAME_1850020_CLASS_ID           = 1850020
       FRAME_1850020_CENTER             = 399
       FRAME_1850020_RELATIVE           = 'J2000'
       FRAME_1850020_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850020_FAMILY             = 'TWO-VECTOR'
       FRAME_1850020_PRI_AXIS           = 'Z'
       FRAME_1850020_PRI_VECTOR_DEF     = 'CONSTANT'
       FRAME_1850020_PRI_FRAME          = 'DT_ETE_OF_DATE'
       FRAME_1850020_PRI_SPEC           = 'RECTANGULAR'
       FRAME_1850020_PRI_VECTOR         = ( 0, 0, 1 )
       FRAME_1850020_SEC_AXIS           = 'X'
       FRAME_1850020_SEC_VECTOR_DEF     = 'OBSERVER_TARGET_VELOCITY'
       FRAME_1850020_SEC_OBSERVER       = 'EARTH'
       FRAME_1850020_SEC_TARGET         = 'MOON'
       FRAME_1850020_SEC_ABCORR         = 'NONE'
       FRAME_1850020_SEC_FRAME          = 'DT_GSE'

       FRAME_DT_PRODUCT                 = 1850024
       FRAME_1850024_NAME               = 'DT_PRODUCT'
       FRAME_1850024_CLASS              = 5
       FRAME_1850024_CLASS_ID           = 1850024
       FRAME_1850024_CENTER             = 399
       FRAME_1850024_RELATIVE           = 'J2000'
       FRAME_1850024_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850024_FAMILY             = 'PRODUCT'
       FRAME_1850024_FROM_FRAMES        = ( 'IAU_EARTH', 'J2000', 'DT_GSE' )
       FRAME_1850024_TO_FRAMES          = ( 'ECLIPJ2000', 'IAU_MOON', 'DT_EME_OF_DATE' )

       FRAME_DT_PRODUCT_ONE             = 1850025
       FRAME_1850025_NAME               = 'DT_PRODUCT_ONE'
       FRAME_1850025_CLASS              = 5
       FRAME_1850025_CLASS_ID           = 1850025
       FRAME_1850025_CENTER             = 399
       FRAME_1850025_RELATIVE           = 'IAU_EARTH'
       FRAME_1850025_DEF_STYLE          = 'PARAMETERIZED'
       FRAME_1850025_FAMILY             = 'PRODUCT'
       FRAME_1850025_FROM_FRAMES        = ( 'IAU_EARTH' )
       FRAME_1850025_TO_FRAMES          = ( 'ECLIPB1950' )

       FRAME_DT_SWITCH_TIME             = 1850021
       FRAME_1850021_NAME               = 'DT_SWITCH_TIME'
       FRAME_1850021_CLASS              = 6
       FRAME_1850021_CLASS_ID           = 1850021
       FRAME_1850021_CENTER             = 399
       FRAME_1850021_ALIGNED_WITH       = ( 'IAU_EARTH', 'DT_GSE', 'ECLIPJ2000', 'DT_EULER' )
       FRAME_1850021_START              = ( -6.0E8, -2.0E8, -1.0E8, 1.0E8 )
       FRAME_1850021_STOP               = (  1.0E8,  2.5E8,  2.5E8, 6.0E8 )

       FRAME_DT_SWITCH_ALL              = 1850022
       FRAME_1850022_NAME               = 'DT_SWITCH_ALL'
       FRAME_1850022_CLASS              = 6
       FRAME_1850022_CLASS_ID           = 1850022
       FRAME_1850022_CENTER             = 399
       FRAME_1850022_ALIGNED_WITH       = ( 'IAU_EARTH', 'DT_ETE_OF_DATE' )

       FRAME_DT_SWITCH_NEST             = 1850023
       FRAME_1850023_NAME               = 'DT_SWITCH_NEST'
       FRAME_1850023_CLASS              = 6
       FRAME_1850023_CLASS_ID           = 1850023
       FRAME_1850023_CENTER             = 399
       FRAME_1850023_ALIGNED_WITH       = ( 'IAU_MOON', 'DT_SWITCH_TIME' )
       FRAME_1850023_START              = ( -6.0E8, -3.0E8 )
       FRAME_1850023_STOP               = (  6.0E8,  3.0E8 )

    \begintext
    """)
    PATHS[:fk] = path
    @assert isfile(path)
    return path
end

"Header lines of the frame `name` with the ID `id` and the class `class`, center 399."
head(name, id, class) = """
       FRAME_$name = $id
       FRAME_$(id)_NAME = '$name'
       FRAME_$(id)_CLASS = $class
       FRAME_$(id)_CLASS_ID = $id
       FRAME_$(id)_CENTER = 399
"""

"Lines of a two-vector frame: the primary vector is the position of the Sun from the Earth."
twovec(id, sec) = """
       FRAME_$(id)_RELATIVE = 'J2000'
       FRAME_$(id)_DEF_STYLE = 'PARAMETERIZED'
       FRAME_$(id)_FAMILY = 'TWO-VECTOR'
       FRAME_$(id)_PRI_AXIS = 'X'
       FRAME_$(id)_PRI_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
       FRAME_$(id)_PRI_OBSERVER = 'EARTH'
       FRAME_$(id)_PRI_TARGET = 'SUN'
       FRAME_$(id)_PRI_ABCORR = 'NONE'
$sec"""

"""
Frames of the kernel of `error_kernel`, as `(name, error kind)`; `:none` for a frame
that is valid. The kind is the short message of CSPICE in lower case.
"""
const ERROR_FRAMES = [
    ("DT_DEEP", :recursiontoodeep), ("DT_DEEP_BASE", :none), ("DT_PARALLEL", :degeneratecase),
    ("DT_NOFAMILY", :kernelvarnotfound), ("DT_BADFAMILY", :notsupported),
    ("DT_BADSTYLE", :notsupported), ("DT_SAMEAXIS", :undefinedframe), ("DT_NOSEC", :kernelvarnotfound),
    ("DT_BADBODY", :notranslation), ("DT_BADABCORR", :spkinvalidoption), ("DT_BADVECDEF", :notsupported),
    ("DT_NOBASE", :kernelvarnotfound), ("DT_SWITCH_GAP", :noframeconnect),
    ("DT_SWITCH_COUNT", :countmismatch),
]

"""
    error_kernel() -> String

Write (once per process) a frame kernel with definitions that are not valid, or that
need more recursion levels than CSPICE has, and return its path. It needs the frames of
`frame_kernel`.
"""
function error_kernel()
    haskey(PATHS, :err) && return PATHS[:err]
    path = joinpath(dir(), "test_dynamic_errors.tf")
    velsec(id, frame) = """
       FRAME_$(id)_SEC_AXIS = 'Y'
       FRAME_$(id)_SEC_VECTOR_DEF = 'OBSERVER_TARGET_VELOCITY'
       FRAME_$(id)_SEC_OBSERVER = 'SUN'
       FRAME_$(id)_SEC_TARGET = 'EARTH'
       FRAME_$(id)_SEC_ABCORR = 'NONE'
       FRAME_$(id)_SEC_FRAME = '$frame'
"""
    io = IOBuffer()
    println(io, "KPL/FK\n\nDefinitions that are not valid (written by the tests).\n\n\\begindata\n")
    # The velocity is given in a frame that itself has vectors in dynamic frames: level 2.
    print(io, head("DT_DEEP", 1850101, 5), twovec(1850101, velsec(1850101, "DT_TWOVEC_IN_DYN")))
    # A chain of base frames is not a recursion.
    print(io, head("DT_DEEP_BASE", 1850102, 5), """
       FRAME_1850102_RELATIVE = 'DT_TWOVEC_IN_DYN'
       FRAME_1850102_DEF_STYLE = 'PARAMETERIZED'
       FRAME_1850102_FAMILY = 'EULER'
       FRAME_1850102_EPOCH = @2000-JAN-01/12:00:00
       FRAME_1850102_AXES = ( 3, 1, 3 )
       FRAME_1850102_UNITS = 'DEGREES'
       FRAME_1850102_ANGLE_1_COEFFS = ( 1.0, 1.0E-7 )
       FRAME_1850102_ANGLE_2_COEFFS = ( 2.0 )
       FRAME_1850102_ANGLE_3_COEFFS = ( 3.0 )
""")
    print(io, head("DT_PARALLEL", 1850103, 5), twovec(1850103, """
       FRAME_1850103_SEC_AXIS = 'Y'
       FRAME_1850103_SEC_VECTOR_DEF = 'OBSERVER_TARGET_POSITION'
       FRAME_1850103_SEC_OBSERVER = 'SUN'
       FRAME_1850103_SEC_TARGET = 'EARTH'
       FRAME_1850103_SEC_ABCORR = 'NONE'
"""))
    print(io, head("DT_NOFAMILY", 1850104, 5), """
       FRAME_1850104_RELATIVE = 'J2000'
       FRAME_1850104_DEF_STYLE = 'PARAMETERIZED'
""")
    print(io, head("DT_BADFAMILY", 1850105, 5), """
       FRAME_1850105_RELATIVE = 'J2000'
       FRAME_1850105_DEF_STYLE = 'PARAMETERIZED'
       FRAME_1850105_FAMILY = 'THREE-VECTOR'
""")
    print(io, head("DT_BADSTYLE", 1850106, 5), """
       FRAME_1850106_RELATIVE = 'J2000'
       FRAME_1850106_DEF_STYLE = 'SCRIPTED'
       FRAME_1850106_FAMILY = 'EULER'
""")
    print(io, head("DT_SAMEAXIS", 1850107, 5), replace(twovec(1850107, velsec(1850107, "J2000")),
        "SEC_AXIS = 'Y'" => "SEC_AXIS = '-X'"))
    print(io, head("DT_NOSEC", 1850108, 5), twovec(1850108, ""))
    print(io, head("DT_BADBODY", 1850109, 5), replace(twovec(1850109, velsec(1850109, "J2000")),
        "PRI_TARGET = 'SUN'" => "PRI_TARGET = 'NO_SUCH_BODY'"))
    print(io, head("DT_BADABCORR", 1850110, 5), replace(twovec(1850110, velsec(1850110, "J2000")),
        "PRI_ABCORR = 'NONE'" => "PRI_ABCORR = 'LT+X'"))
    print(io, head("DT_BADVECDEF", 1850111, 5), replace(twovec(1850111, velsec(1850111, "J2000")),
        "'OBSERVER_TARGET_POSITION'" => "'OBSERVER_TARGET_DIRECTION'"))
    print(io, head("DT_NOBASE", 1850112, 5), """
       FRAME_1850112_DEF_STYLE = 'PARAMETERIZED'
       FRAME_1850112_FAMILY = 'MEAN_EQUATOR_AND_EQUINOX_OF_DATE'
       FRAME_1850112_PREC_MODEL = 'EARTH_IAU_1976'
       FRAME_1850112_ROTATION_STATE = 'ROTATING'
""")
    print(io, head("DT_SWITCH_GAP", 1850113, 6), """
       FRAME_1850113_ALIGNED_WITH = ( 'IAU_EARTH', 'DT_GSE' )
       FRAME_1850113_START = ( 1.0E9, 2.0E9 )
       FRAME_1850113_STOP = ( 1.5E9, 2.5E9 )
""")
    print(io, head("DT_SWITCH_COUNT", 1850114, 6), """
       FRAME_1850114_ALIGNED_WITH = ( 'IAU_EARTH', 'DT_GSE', 'ECLIPJ2000' )
       FRAME_1850114_START = ( -1.0E9, 2.0E9 )
       FRAME_1850114_STOP = ( 1.5E9, 2.5E9 )
""")
    println(io, "\n\\begintext")
    write(path, take!(io))
    return PATHS[:err] = path
end

"""
    switch_ck_kernel() -> String

Write (once per process) a frame kernel with switch frames over the CK frames of
`test/ckfixtures.jl` and return its path. A CK base is used only at an epoch at which
its CK has pointing, so the selection follows the coverage of the CK segments.
"""
function switch_ck_kernel()
    haskey(PATHS, :swck) && return PATHS[:swck]
    path = joinpath(dir(), "test_switch_ck.tf")
    write(path, "KPL/FK\n\n\\begindata\n" * head("DT_SWITCH_CK", 1850121, 6) * """
       FRAME_1850121_ALIGNED_WITH = ( 'IAU_EARTH', 'TESTSC_T3', 'TESTSC_T2' )
""" * head("DT_SWITCH_CK_NOAV", 1850122, 6) * """
       FRAME_1850122_ALIGNED_WITH = ( 'ECLIPJ2000', 'TESTSC_T1', 'TESTSC_T3_NOAV' )
""" * "\n\\begintext\n")
    return PATHS[:swck] = path
end

end # module
