Feature: CK frames (frame class 3)
  A frame kernel can define a frame whose orientation comes from a C-kernel.
  Such a frame takes part in rotations and state transformations as other frames do.

  Background:
    Given a context with the standard kernels, an SCLK kernel, a frame kernel and CK files loaded

  Scenario: Rotation to and from a CK frame
    Given a CK frame with pointing at the epoch
    When pxform is called between that frame and an inertial, a body-fixed, a fixed-offset or a second CK frame
    Then the rotation is equal to the CSPICE rotation to rounding accuracy

  Scenario: A CK frame with a base frame that is not J2000
    Given a segment whose reference frame is IAU_EARTH or ECLIPJ2000
    When pxform is called
    Then the rotation of the base frame is included

  Scenario: State transformation
    Given a CK frame with angular velocity at the epoch
    When sxform or spkezr is called
    Then the derivative block uses the angular velocity of the segment, as CSPICE does
    And the result is equal to the CSPICE answer to rounding accuracy

  Scenario: Exact derivatives
    Given a CK frame with pointing at the epoch
    When pxder is called with an order of 1 or more
    Then the result is the derivative of the interpolated rotation
    And it agrees with a central finite difference of pxform
    And the query does not allocate memory

  Scenario: A frame without data
    Given a CK frame without pointing at the epoch, or without a clock kernel
    When a rotation is requested
    Then a SpiceError of the kind noframeconnect is thrown
    And CSPICE fails for the same request

  Scenario: State transformation without angular velocity
    Given a CK frame whose segments hold no angular velocity
    When sxform is called
    Then a SpiceError is thrown
    And pxform for the same request gives the rotation

  Scenario: Concurrent readers
    Given several tasks
    When they request rotations to CK frames at the same time
    Then every result is equal to the serial result
