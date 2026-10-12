Feature: Geometry finder, quantities from the ephemeris
  The searches for distance, range rate, angular separation and phase angle
  give the same windows as CSPICE: the same number of intervals, and end points
  that agree within twice the convergence tolerance.

  Background:
    Given a context with a leapseconds kernel, a planetary constants kernel and "de440s.bsp" loaded

  Scenario Outline: A search against CSPICE
    Given a target, an observer and an aberration correction from "NONE", "LT", "LT+S", "CN", "CN+S", "XLT", "XLT+S", "XCN", "XCN+S"
    And a random confinement window of one or of several intervals
    When "<call>" searches for each of the relations "<", "=", ">", "LOCMIN", "LOCMAX", "ABSMIN", "ABSMAX", with and without an adjust value
    Then the result has the same number of intervals as the result of CSPICE
    And each end point agrees with the end point of CSPICE within 2e-6 seconds

    Examples:
      | call   |
      | gfdist |
      | gfrr   |
      | gfsep  |
      | gfpa   |

  Scenario: An extremum whose derivative test is rounding noise
    Given a search for an extremum of the range rate, or of the angular separation with stellar aberration
    When an end point differs from the end point of CSPICE by more than 2e-6 seconds
    Then the difference is less than 1 second
    And the quantity has the same value at the two epochs, to the rounding of the states

  Scenario: An event at an end point of the confinement window
    Given a confinement window whose end points are roots of the quantity
    When the search runs for each relation
    Then the result has the same end points as the result of CSPICE

  Scenario: Angular separation of two spheres
    Given the shapes "SPHERE" for the two targets
    When "gfsep" searches for a separation below a reference value
    Then the separation is the angle between the limbs, as in CSPICE

  Scenario: Phase angle with a transmission correction
    When "gfpa" gets an aberration correction that starts with "X"
    Then it throws a SpiceError of kind "invalidoption"

  Scenario: The general call
    Given the quantity name, the parameter names and the parameter values of "gfevnt"
    When "gfevnt" searches with a step function and a refinement function
    Then the result is equal to the result of the specific call for that quantity
    And a quantity that needs surface or field of view geometry throws a SpiceError of kind "unsupported"

  Scenario: Search time
    Given a search of ten years with a step of one day
    When NeoSpice and CSPICE make the same search
    Then the two search times and their ratio are reported
