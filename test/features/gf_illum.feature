Feature: Geometry finder, illumination angles at a surface point
  "gfilum" finds the times when the phase angle, the incidence angle or the emission
  angle at a point of the surface of a target satisfies a relation. The windows are
  those of CSPICE.

  Scenario Outline: An illumination angle against CSPICE
    Given a context with the standard kernels loaded
    And a point of the surface of a target, a source of light and an observer
    When "gfilum" searches for <relation> on the angle <angle>
    Then the window has the same number of intervals as the window of "gfilum_c"
    And each end point agrees with CSPICE within twice the convergence tolerance
    And the window is a subset of the confinement window

    Examples:
      | angle     | relation |
      | PHASE     | <        |
      | INCIDENCE | =        |
      | EMISSION  | LOCMIN   |
      | PHASE     | ABSMAX   |
      | INCIDENCE | >        |

  Scenario: The value and its rate
    Given an illumination angle at a surface point
    When the value and the rate are computed at an epoch
    Then the value is the angle that "illumg_c" gives
    And the rate agrees with a central difference of the value
    And the angle crosses the reference value at each end point of an equality search

  Scenario: The source of light is the observer
    Given a search where the source of light and the observer are the same body
    When "gfilum" searches
    Then no error is thrown and the window is the window of CSPICE

  Scenario: The general call
    Given the quantity "ILLUMINATION ANGLE" of "gfevnt"
    When the parameters are those of a "gfilum" search, in any case
    Then the window is equal to the window of "gfilum"
    And fewer than 3 values of the surface point give the error "invalidcount"
    And a parameter that is not supplied gives the error "missingvalue"

  Scenario: Forms, keywords and concurrent use
    Given the same search on a kernel set, on a context and with integer codes
    When each form runs, the search runs with "threads" set to true, and many tasks search at the same time
    Then every window is equal, bit for bit, to the window of the serial search
    And the callback "bail" stops the search with an InterruptException

  Scenario Outline: Bad inputs
    Given a valid search with one input replaced by <input>
    When the search runs in NeoSpice and in CSPICE
    Then NeoSpice throws a SpiceError of kind <kind>
    And CSPICE signals the matching short error

    Examples:
      | input                              | kind              |
      | a method other than ELLIPSOID      | invalidmethod     |
      | an angle type that is not known    | notsupported      |
      | a body that is not known           | idcodenotfound    |
      | an observer equal to the target    | bodiesnotdistinct |
      | a source equal to the target       | bodiesnotdistinct |
      | a frame that is not known          | unknownframe      |
      | a frame not centered on the target | invalidframe      |
      | a bad aberration correction        | invalidoption     |
      | a transmission correction          | notsupported      |
      | a target without radii             | kernelvarnotfound |
      | a relation that is not known       | notrecognized     |
      | a negative adjust value            | valueoutofrange   |
      | a step of zero                     | invalidstep       |
