Feature: Aberration corrections
  A user asks for the state or the position of a body as an observer sees it,
  corrected for the one-way light time and for stellar aberration, at reception
  or at transmission. The corrections are those of CSPICE, with the same
  correction strings, and the results agree with CSPICE to rounding level.

  Background:
    Given a context with the standard kernels and "de440s.bsp" loaded
    And CSPICE with the same kernels loaded

  Scenario Outline: State with a correction in any frame
    When the user calls "spkezr" for the Moon relative to the Sun in "<frame>" with correction "<abcorr>"
    Then the state agrees with "spkezr_c" of CSPICE to rounding level
    And the light time agrees with CSPICE to rounding level
    And the call does not allocate and its result type is inferred

    Examples:
      | frame         | abcorr |
      | J2000         | LT     |
      | J2000         | LT+S   |
      | ECLIPJ2000    | CN     |
      | ECLIPJ2000    | CN+S   |
      | IAU_EARTH     | XLT    |
      | IAU_MOON      | XLT+S  |
      | MOON_PA_DE440 | XCN    |
      | MOON_PA_DE440 | XCN+S  |

  Scenario: Position with a correction
    When the user calls "spkpos" with correction "CN+S"
    Then the position and the light time agree with "spkpos_c" of CSPICE
    And the position agrees with the position part of "spkezr" to rounding level

  Scenario: Frame that is not inertial
    When the user asks for a corrected state in "IAU_MOON" for the Earth seen from the Sun
    Then the frame is evaluated at the epoch corrected by the light time between the observer and the Moon
    And the velocity includes the rate of that light time

  Scenario Outline: Spelling of the correction
    When the user gives the correction as "<text>"
    Then the result is the same as for "<canonical>"

    Examples:
      | text    | canonical |
      | lt+s    | LT+S      |
      | " xcn " | XCN       |
      | C N + S | CN+S      |
      | none    | NONE      |

  Scenario: Correction given as a parsed value
    Given the value "Aberration(\"LT+S\")"
    When the user passes the value in the place of the string
    Then the result is the same as for the string "LT+S"
    And the value is plain bits with the flags "uselt", "converged", "stellar", "transmit"

  Scenario Outline: Invalid correction
    When the user asks for a state with correction "<abcorr>"
    Then a "SpiceError" with kind "badarg" is thrown

    Examples:
      | abcorr |
      | RL     |
      | XRL+S  |
      | S      |
      | XS     |
      | LT+    |
      | FOO    |

  Scenario: Body relative to itself
    When the user asks for the corrected state of the Earth relative to the Earth
    Then the state is zero and the light time is zero

  Scenario Outline: Derivatives of the corrected position
    When the user calls "spkder" with correction "<abcorr>" and "Val(3)" in "<frame>"
    Then the result is a tuple of 4 vectors of length 3
    And element 1 is the position that "spkezr" gives, to rounding level
    And each element k + 1 agrees with a central finite difference of element k
    And element 2 agrees with a central finite difference of the "spkpos" position
    And every element agrees with nested "ForwardDiff.derivative" of the "spkpos" position

    Examples:
      | frame         | abcorr |
      | J2000         | LT     |
      | J2000         | CN+S   |
      | IAU_EARTH     | XLT+S  |
      | MOON_PA_DE440 | XCN    |

  Scenario: First derivative and the CSPICE velocity
    When the user compares element 2 of "spkder" with the velocity of "spkezr"
    Then they agree to rounding level for the "CN" forms
    And they differ by a relative amount near 1e-8 for the "LT" forms, because CSPICE uses the light-time rate of the converged equation

  Scenario: Low-level calls
    When the user calls "spkssb", "spkltc", "spkaps", "spkacs", "spkapo", "ltime", "stelab" or "stlabx"
    Then the result agrees with the CSPICE routine of the same name
    And a frame that is not inertial throws a "SpiceError" with kind "badarg" where CSPICE requires an inertial frame
