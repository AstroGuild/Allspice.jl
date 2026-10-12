Feature: Coordinate conversions
  A user converts points between rectangular, latitudinal, right ascension and
  declination, spherical, cylindrical, geodetic and planetographic coordinates, asks
  for the Jacobian of a conversion, converts a state between coordinate systems and
  converts between units. The results agree with CSPICE to rounding level.

  Background:
    Given points with magnitudes from 1e-10 to 1e10
    And the degenerate points: the origin, the poles, the x-y plane and the negative x axis
    And CSPICE called directly as the oracle

  Scenario Outline: Conversion of a point and back
    When the user calls "<to>" for a point
    Then the coordinates agree with CSPICE to rounding level
    And "<back>" of the coordinates gives the point back

    Examples:
      | to     | back   |
      | reclat | latrec |
      | recrad | radrec |
      | recsph | sphrec |
      | reccyl | cylrec |
      | recgeo | georec |
      | recpgr | pgrrec |
      | cyllat | latcyl |
      | cylsph | sphcyl |
      | sphlat | latsph |

  Scenario: Angles of a degenerate point
    When the user converts the origin or a point on the z axis
    Then the angles that are not defined are 0, as in CSPICE
    And no error is thrown

  Scenario: Range of the longitude
    When the user calls "reclat" for a point with a negative x and a y of plus or minus zero
    Then the longitude is plus or minus pi with the sign of y
    And "recrad" gives a right ascension from 0 to 2 pi

  Scenario: Geodetic coordinates of a spheroid that is not valid
    When the user calls "recgeo" or "georec" with a radius that is not positive or a flattening of 1 or more
    Then a SpiceError with the kind "valueoutofrange" is thrown

  Scenario Outline: Sense of the planetographic longitude
    Given a context with a text PCK loaded
    When the user calls "recpgr" for "<body>"
    Then the longitude is positive "<sense>"

    Examples:
      | body  | sense |
      | EARTH | east  |
      | MOON  | east  |
      | SUN   | east  |
      | MARS  | west  |
      | VENUS | east  |

  Scenario: Longitude sense from the kernel pool
    Given the pool variable "BODY499_PGR_POSITIVE_LON" with the value "EAST"
    When the user calls "recpgr" for Mars
    Then the longitude is positive east
    And a value other than "EAST" or "WEST" gives a SpiceError with the kind "invalidoption"

  Scenario Outline: Jacobian of a conversion
    When the user calls "<call>"
    Then the matrix agrees with CSPICE to rounding level
    And the matrix agrees with the Jacobian of the conversion by ForwardDiff

    Examples:
      | call   |
      | drdlat |
      | dlatdr |
      | drdsph |
      | dsphdr |
      | drdcyl |
      | dcyldr |
      | drdgeo |
      | dgeodr |
      | drdpgr |
      | dpgrdr |

  Scenario: Jacobian on the z axis
    When the user calls "dlatdr", "dsphdr", "dcyldr", "dgeodr" or "dpgrdr" for a point on the z axis
    Then a SpiceError with the kind "pointonzaxis" is thrown

  Scenario: State between coordinate systems
    Given a state in one of the six coordinate systems
    When the user calls "xfmsta" with an input and an output system
    Then the state agrees with "xfmsta_c" of CSPICE
    And the rates are the time derivatives of the coordinates
    And a system name that is not known gives a SpiceError with the kind "coordsysnotrec"

  Scenario: Unit conversion
    When the user calls "convrt" with two units of the same class (angle, distance or time)
    Then the value is equal to the value of "convrt_c" bit for bit
    And units of different classes give a SpiceError with the kind "incompatibleunits"
    And a unit that is not known gives a SpiceError with the kind "unitsnotrec"

  Scenario: Static arguments
    When the user calls a conversion that does not need the kernel pool
    Then the result type is inferred and the call does not allocate
