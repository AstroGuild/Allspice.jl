Feature: Two-body motion
  A user propagates a state on a conic orbit, converts between a state and osculating
  elements and evaluates the Stumpff functions with the CSPICE calls. The results agree
  with CSPICE to rounding level and the same inputs are accepted and rejected.

  Background:
    Given ellipses, parabolas and hyperbolas about masses from 1e-3 to 1e11 km^3/s^2
    And the degenerate orbits: circular, equatorial, retrograde equatorial and next to a parabola
    And CSPICE called directly as the oracle

  Scenario: Propagation of a state
    When the user calls "prop2b" with a state and a time interval
    Then the state agrees with "prop2b_c" of CSPICE to rounding level
    And the energy and the angular momentum are the same as before

  Scenario: Propagation there and back
    When the user propagates a state by a time interval and then by the opposite interval
    Then the result is the first state to rounding level

  Scenario Outline: Propagation of a state that is not valid
    When the user calls "prop2b" with <input>
    Then a SpiceError with the kind "<kind>" is thrown

    Examples:
      | input                                   | kind            |
      | a mass that is not positive             | nonpositivemass |
      | a zero position                         | zeroposition    |
      | a zero velocity                         | zerovelocity    |
      | a velocity parallel to the position     | nonconicmotion  |

  Scenario: Derivative with respect to the time
    When the user calls "prop2b" or "conics" with a ForwardDiff or Taylor time
    Then the first derivative of the position is the velocity
    And the second derivative of the position is the two-body acceleration

  Scenario: State from elements
    When the user calls "conics" with eight elements and an epoch
    Then the state agrees with "conics_c" of CSPICE
    And the state at the epoch of the elements with a zero mean anomaly is at the perifocal distance

  Scenario Outline: Elements that are not valid
    When the user calls "conics" with <input>
    Then a SpiceError with the kind "<kind>" is thrown

    Examples:
      | input                                 | kind             |
      | a negative eccentricity               | badeccentricity  |
      | a perifocal distance that is not positive | badperiapsevalue |
      | a mass that is not positive           | badgm            |

  Scenario: Elements from a state and back
    When the user calls "oscelt" with a state and then "conics" with the elements at the same epoch
    Then the result is the state to rounding level
    And the elements are equal to the elements of "oscelt_c"

  Scenario: Extended elements
    When the user calls "oscltx"
    Then places 1 to 8 are the elements of "oscelt"
    And places 9 to 11 are the true anomaly, the semi-major axis and the period
    And the semi-major axis and the period are 0 where they are not defined

  Scenario: Stumpff functions
    When the user calls "stmp03" with an argument
    Then the four values agree with "stmp03_" of CSPICE to rounding level
    And an argument below the lower bound gives a SpiceError with the kind "valueoutofrange"

  Scenario: Static arguments
    When the user calls a two-body call with static arrays
    Then the result type is inferred and the call does not allocate
