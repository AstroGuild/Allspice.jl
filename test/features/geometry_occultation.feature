Feature: Occultation state
  A user asks whether one body hides another one from an observer, in part or in full.
  The answer is the code of CSPICE, from -3 to 3.

  Background:
    Given a context with the standard kernels, "de440s.bsp", a DSK file for the Moon
    And synthetic bodies whose ephemerides are written with the SPK writers
    And CSPICE with the same kernels loaded

  Scenario Outline: Codes of real events
    When the user calls "occult" at epochs about "<event>"
    Then every code is the code of "occult_c"

    Examples:
      | event                                        |
      | the solar eclipse of 2017-08-21              |
      | the annular solar eclipse of 2012-05-20      |
      | the transit of Venus of 2012-06-06           |
      | the lunar eclipse of 2019-01-21 from the Moon |

  Scenario: All nine codes
    Given two triaxial bodies, one in front of the other
    When the radii of the back body are larger, smaller, or of the same size
    Then the scans give the codes -3, -2, -1, 0, 1, 2 and 3
    And the order of the two targets changes the sign of the code

  Scenario: Point targets
    When one target has the shape "POINT"
    Then the code is 0, 3 or -3
    When the two targets have the shape "POINT"
    Then a "SpiceError" is thrown, of the kind that CSPICE gives

  Scenario: DSK shape against a point
    When the Moon has the shape "DSK/UNPRIORITIZED" and the other target is a point
    Then every code is the code of "occult_c"

  Scenario: No allocation
    When the user calls "occult" with integer IDs for ellipsoids or points
    Then the call does not allocate
