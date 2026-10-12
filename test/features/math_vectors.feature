Feature: Vector and matrix calls
  A user calls the vector and matrix routines of CSPICE by their CSPICE names, with
  static arrays or with plain arrays. The results agree with CSPICE, the same inputs
  are accepted and rejected, and the calls with static arguments do not allocate.

  Background:
    Given vectors with magnitudes from 1e-10 to 1e10
    And the degenerate vectors: zero, parallel, opposite, on an axis, 1e-300 and 1e300
    And CSPICE called directly as the oracle

  Scenario Outline: Arithmetic call
    When the user calls "<call>" with random and degenerate arguments
    Then the result is equal to the result of CSPICE bit for bit

    Examples:
      | call   |
      | vadd   |
      | vsub   |
      | vscl   |
      | vdot   |
      | vcrss  |
      | vnorm  |
      | vhat   |
      | unorm  |
      | vdist  |
      | vproj  |
      | vlcom3 |
      | ucrss  |
      | mxv    |
      | mtxv   |
      | mxm    |
      | mxmt   |
      | xpose6 |
      | invert |
      | twovec |
      | dvhat  |
      | dvcrss |

  Scenario Outline: Call with a trigonometric function
    When the user calls "<call>" with random and degenerate arguments
    Then the result agrees with CSPICE to a few units of rounding

    Examples:
      | call   |
      | vsep   |
      | vrotv  |
      | rotate |
      | rotmat |
      | rotvec |
      | dvsep  |

  Scenario: Zero vector
    When the user calls "vhat" or "unorm" with the zero vector
    Then the result is the zero vector and no error is thrown

  Scenario: Infinite numerator in the derivative of the separation
    Given two states of which one has a unit position with an infinite derivative
    When the user calls "dvsep"
    Then the error "numericoverflow" is thrown
    And the test does not call "dvsep_c", because CSPICE stops the process for this input

  Scenario: Singular matrix
    When the user calls "invert" with a matrix that has a determinant below 1e-16
    Then the result is the zero matrix, as in CSPICE

  Scenario: Reference frame from two vectors that are parallel
    When the user calls "twovec" with two parallel vectors
    Then a SpiceError with the kind "dependentvectors" is thrown

  Scenario: General-dimension call
    Given two vectors of the same length from 1 to 9
    When the user calls a call that ends with "g"
    Then the result agrees with CSPICE
    And vectors of different lengths give a SpiceError with the kind "badarg"

  Scenario: Angle between parallel vectors of any dimension
    When the user calls "vsepg" with two parallel vectors
    Then the angle is 0 exactly, as "vsepg_" of the Fortran-derived library gives
    And the result of "vsepg_c", which has its two branches exchanged, is not followed

  Scenario: Derivatives
    When the user calls a vector call with ForwardDiff or Taylor numbers
    Then the derivative agrees with the analytic derivative call ("dvhat", "dvdot", "dvcrss", "dvnorm", "dvsep")

  Scenario: Static arguments
    When the user calls a vector or matrix call with static arrays
    Then the result type is inferred and the call does not allocate
