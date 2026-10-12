Feature: Observers and targets with a constant position or velocity
  A user has an observer or a target that is at rest, or moves with a constant
  velocity, in a frame of a body: a ground station or a landmark. The user asks for the
  state of the other object, or for its azimuth and elevation at the station.

  Background:
    Given a context with the standard kernels and "de440s.bsp" loaded
    And CSPICE with the same kernels loaded

  Scenario Outline: State relative to a fixed observer or target
    When the user calls "<call>" with the frame evaluated at "<refloc>" and any correction
    Then the state and the light time agree with CSPICE to rounding level
    And the call does not allocate and its result type is inferred

    Examples:
      | call   | refloc   |
      | spkcpo | OBSERVER |
      | spkcpt | TARGET   |
      | spkcvo | CENTER   |
      | spkcvt | OBSERVER |

  Scenario: Epoch as a dual number
    When the user calls "spkcpo" with a "ForwardDiff.Dual" epoch
    Then the derivative of the position agrees with a central difference
    And it agrees with the returned velocity to the accuracy of the CSPICE velocity

  Scenario: Azimuth and elevation at a station
    When the user calls "azlcpo" with the two sense flags
    Then the range, the angles, their rates and the light time agree with "azlcpo_c"

  Scenario Outline: Azimuth and elevation coordinates
    When the user calls "<call>"
    Then the result agrees with CSPICE for each value of the two sense flags

    Examples:
      | call   |
      | azlrec |
      | recazl |
      | dazldr |
      | drdazl |

  Scenario: Jacobian on the vertical
    When the user calls "dazldr" for a point on the Z axis
    Then a "SpiceError" of kind "pointonzaxis" is thrown
