Feature: Illumination angles
  A user asks for the phase, incidence and emission angles at a point of a surface, for
  the Sun or for another source, and for the phase angle of a body as a whole. The
  results agree with CSPICE to rounding level.

  Background:
    Given a context with the standard kernels, "de440s.bsp" and a DSK file for the Moon
    And CSPICE with the same kernels loaded

  Scenario Outline: Angles at a surface point
    Given a point on the surface of "<target>"
    When the user calls "<call>" with correction "<abcorr>"
    Then the epoch, the vector from the observer and the three angles agree with CSPICE
    And the call does not allocate on an ellipsoid

    Examples:
      | call   | target | abcorr |
      | ilumin | MOON   | LT+S   |
      | illumg | EARTH  | CN     |
      | illumf | VENUS  | CN+S   |
      | illumf | MOON   | XLT    |

  Scenario: Visibility and lighting flags
    When the user calls "illumf"
    Then "visibl" is true where the emission angle is below 90 degrees
    And "lit" is true where the incidence angle is below 90 degrees
    And the two flags agree with "illumf_c"

  Scenario: Normal from DSK data
    When the method is "DSK/UNPRIORITIZED/SURFACES = 2"
    Then the angles use the outward normal of the plate, as in CSPICE

  Scenario: Phase angle of a body
    When the user calls "phaseq" for a target, a source and an observer
    Then the angle agrees with "phaseq_c"
    When the correction calls for transmission
    Then a "SpiceError" is thrown, of the kind that CSPICE gives

  Scenario: Epoch as a dual number
    When the user calls "ilumin", "illumg" or "phaseq" on an ellipsoid with a "ForwardDiff.Dual" epoch
    Then the derivative parts of the angles agree with a central difference
