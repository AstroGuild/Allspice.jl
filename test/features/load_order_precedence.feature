Feature: Precedence by load order
  When two kernels give data for the same body, or a value for the same pool
  variable, the kernel that was loaded last has precedence, as in CSPICE.

  Scenario Outline: Two ephemerides for the same bodies
    Given a context with "<first>" loaded and then "<second>" loaded
    And CSPICE with the same kernels loaded in the same order
    When the user asks for the state of the Moon relative to the Sun at an epoch that both files cover
    Then the answer agrees with CSPICE
    And the answer is equal to the answer of a context with only "<second>" loaded

    Examples:
      | first      | second     |
      | de405.bsp  | de440s.bsp |
      | de440s.bsp | de405.bsp  |

  Scenario: The earlier kernel answers outside the coverage of the later kernel
    Given a context with "de405.bsp" loaded and then "de440s.bsp" loaded
    When the user asks for a state at an epoch that only "de405.bsp" covers
    Then the answer is equal to the answer of a context with only "de405.bsp" loaded

  Scenario: Two text kernels assign the same variable
    Given a context with "pck00010.tpc" loaded and then "pck00011.tpc" loaded
    When the user asks for "bodvrd" of "MARS", "PM"
    Then the values are those of "pck00011.tpc"
