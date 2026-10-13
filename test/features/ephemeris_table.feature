Feature: Ephemeris table for a fixed interval of time
  A user who needs the position of one body relative to another many times,
  inside an optimizer or on a GPU, extracts it once for the interval of the
  problem. The table is an array of Chebyshev coefficients and a few numbers,
  and a query needs no kernel.

  Background:
    Given a context with the standard kernels and "de440s.bsp" loaded

  Scenario Outline: Geometric position in an inertial frame
    When the user calls "ephemeris_table" for <target> relative to <observer> in "<frame>" with "NONE"
    Then the table has the record length and the degree of the kernels
    And the span of the table contains the requested interval
    And "table_state" agrees with "spkezr" of CSPICE to rounding at every epoch of the span

    Examples:
      | target | observer | frame      |
      | 301    | 399      | J2000      |
      | 399    | 10       | ECLIPJ2000 |
      | 4      | 399      | J2000      |

  Scenario: Rotating frame or aberration correction
    When the user calls "ephemeris_table" with a body-fixed frame or with "LT+S"
    Then the function interpolates and doubles the number of intervals
    And the field "maxerr" of the table is the measured error in km
    And "table_position" agrees with "spkpos" of CSPICE within that error

  Scenario: Tolerance that the source cannot give
    When the user gives a tolerance below the noise of the source
    Then the call throws "SpiceError" of kind "noconvergence"

  Scenario: Derivatives
    When the user calls "table_jet" with "Val(N)" or "table_position" with a dual-number epoch
    Then every order is the derivative of the position polynomial of the table

  Scenario: Epoch outside the table
    When the user asks for an epoch before the start or after the end of the table
    Then the call throws "SpiceError" of kind "nocoverage"

  Scenario: Other array types
    When the user converts the table with "adapt" or to "Float32" coefficients
    Then the scalar fields are equal and the coefficients have the new array type
    And a table whose array is an isbits value is an isbits value

  Scenario: Export as an SPK segment
    Given an ephemeris table of a geometric position
    When the user writes it to an SPK file with "table_to_spk"
    Then "spkezr" on that file gives the state of the table to rounding, in Allspice and in CSPICE

  Scenario: Rotation table
    Given two frames and an interval of time with data
    When the user calls "rotation_table"
    Then "table_rotation" gives the matrix of "pxform" within the measured error of the table
    And "table_rotation_jet" gives its time derivatives without the kernels
