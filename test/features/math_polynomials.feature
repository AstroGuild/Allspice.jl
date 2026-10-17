Feature: Polynomials and interpolation
  A user evaluates a Chebyshev expansion with its derivatives and its integral, interpolates
  tabulated values with a Lagrange or a Hermite polynomial, and calls the small polynomial
  utilities of CSPICE. The results agree with CSPICE to rounding level and the same inputs
  are accepted and rejected.

  Background:
    Given random coefficients, abscissas and values
    And CSPICE called directly as the oracle

  Scenario: Chebyshev expansion
    When the user calls "chbval", "chbint", "chbder" or "chbigr" with coefficients, a midpoint and a radius
    Then the value, the derivatives and the integral agree with CSPICE to rounding level
    And the derivatives of "chbder" are those of a Taylor number

  Scenario: Interpolation of tabulated values
    When the user calls "lgrint", "lgrind", "lgresp", "hrmint" or "hrmesp" with a table
    Then the value and the derivative agree with CSPICE to rounding level
    And a polynomial of a degree below the number of conditions is reproduced

  Scenario Outline: Input that is not valid
    When the user calls an interpolation call with <input>
    Then a SpiceError with the kind "<kind>" is thrown, as CSPICE gives

    Examples:
      | input                                   | kind            |
      | an empty table                          | invalidsize     |
      | two equal abscissas                     | dividebyzero    |
      | an equally spaced table with step 0     | invalidstepsize |
      | a Chebyshev radius that is not positive | invalidradius   |

  Scenario: Polynomial utilities
    When the user calls "polyds", "qderiv", "rquad" or "diags2"
    Then the result agrees with CSPICE to rounding level
    And "rquad" with a zero leading and linear coefficient gives a SpiceError with the kind "degeneratecase"

  Scenario: Derivative with respect to the argument
    When the user calls an interpolation call with a ForwardDiff or Taylor argument
    Then the first derivative agrees with the derivative that the call returns

  Scenario: Static arguments
    When the user calls a polynomial call with static arrays
    Then the result type is inferred and the call does not allocate
