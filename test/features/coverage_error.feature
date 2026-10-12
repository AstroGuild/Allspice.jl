Feature: Errors for queries that cannot be answered
  A query that the loaded kernels cannot answer throws a "SpiceError". The
  field "kind" tells the caller which condition occurred.

  Background:
    Given a context with "naif0012.tls", "pck00011.tpc" and "de440s.bsp" loaded

  Scenario Outline: Epoch outside the coverage
    When the user asks for the state of body 301 relative to body 399 at ET <epoch> seconds
    Then a "SpiceError" with kind "nocoverage" is thrown
    And the message is not empty

    Examples:
      | epoch    |
      | 1.0e11   |
      | -1.0e11  |
      | 4.74e9   |

  Scenario: The ends of the coverage are inside the coverage
    When the user asks for a state at the first or the last epoch of the coverage
    Then the answer agrees with CSPICE

  Scenario: Unknown body name
    When the user asks for the state of "NO_SUCH_BODY" relative to "EARTH"
    Then a "SpiceError" with kind "unknownbody" is thrown

  Scenario: Unknown frame name
    When the user asks for a state in the frame "NO_SUCH_FRAME"
    Then a "SpiceError" with kind "unknownframe" is thrown

  Scenario: Aberration correction that the prototype does not have
    When the user asks for a state with correction "LT+S"
    Then a "SpiceError" with kind "unsupported" is thrown

  Scenario: Kernel file that does not exist
    When the user loads a path that is not a file
    Then a "SpiceError" with kind "nofile" is thrown
    And the context is unchanged
