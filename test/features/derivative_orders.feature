Feature: Time derivatives of any order
  A user asks for the position of a body and its time derivatives up to a
  chosen order, or for a rotation matrix and its time derivatives. The
  derivatives are analytical: they come from the same polynomials as the
  position.

  Background:
    Given a context with the standard kernels and "de440s.bsp" loaded

  Scenario Outline: Position jet of order N
    When the user calls "spkder" for the Moon relative to the Earth in "<frame>" with "Val(<N>)"
    Then the result is a tuple of <N> + 1 vectors of length 3
    And element 1 is the position that "spkezr" gives
    And element 2, if present, is the velocity that "spkezr" gives
    And each element k + 1 agrees with a central finite difference of element k

    Examples:
      | frame     | N |
      | J2000     | 0 |
      | J2000     | 1 |
      | J2000     | 4 |
      | IAU_EARTH | 3 |
      | MOON_PA   | 3 |

  Scenario: Rotation jet
    When the user calls "pxder" from "J2000" to "IAU_EARTH" with "Val(2)"
    Then element 1 is the "pxform" matrix
    And element 2 is the lower left block of the "sxform" matrix
    And element 3 agrees with a central finite difference of element 2

  Scenario: Automatic differentiation through the generic path
    Given the package "ForwardDiff" is loaded
    When the user differentiates "t -> spkpos(...)" with "ForwardDiff.derivative", nested k times
    Then the epoch that the library receives is a dual number
    And the result agrees with element k + 1 of "spkder"

  Scenario: Constant rotation
    When the user calls "pxder" between two inertial frames
    Then every element after the first is the zero matrix
