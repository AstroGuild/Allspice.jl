Feature: Load kernels, then query states and rotations
  A user loads SPICE kernels and asks for the state of one body relative to
  another. The answers agree with CSPICE for the same kernels.

  Background:
    Given a new context
    And the kernels "naif0012.tls", "pck00011.tpc" and "de440s.bsp" are loaded in that order

  Scenario: Geometric state by body names
    When the user asks for the state of "MOON" relative to "EARTH" in "J2000" with correction "NONE" at an epoch inside the coverage
    Then the result is a 6-vector in km and km/s and a one-way light time in seconds
    And each component agrees with CSPICE to within rounding

  Scenario: Integer IDs give the same answer as names
    When the user asks for the state of body 301 relative to body 399 in frame 1 at the same epoch
    Then the result is equal, bit for bit, to the result of the name form

  Scenario: State in a body-fixed frame
    When the user asks for the state of "MOON" relative to "EARTH" in "IAU_EARTH"
    Then the result agrees with CSPICE
    And the result is equal to the "sxform" matrix from "J2000" to "IAU_EARTH" times the "J2000" state

  Scenario: Rotation between two frames
    When the user asks for "pxform" from "J2000" to "IAU_MARS" at an epoch
    Then the result is a 3 by 3 rotation matrix that agrees with CSPICE
    And the matrix times its transpose is the identity

  Scenario: Rotation between two frames at two epochs
    When the user asks for "pxfrm2" from "IAU_EARTH" at one epoch to "IAU_MOON" at another epoch
    Then the result agrees with "pxfrm2_c" of CSPICE
    And with two equal epochs the result is the "pxform" matrix

  Scenario: Query without a context argument
    Given the same kernels are loaded with "furnsh(path)" and no context argument
    When the user calls "spkezr" with no context argument
    Then the answer is equal to the answer of the explicit context

  Scenario: Integer-ID forms of spkezr and spkpos, and the geometric position
    When the user calls "spkez" and "spkezp" with integer IDs, a frame and a correction
    Then the answers are equal to those of "spkezr" and "spkpos" and agree with CSPICE
    And "spkgps" gives the position and light time of "spkgeo"
