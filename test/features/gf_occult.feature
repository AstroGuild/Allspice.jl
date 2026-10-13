Feature: Geometry finder, occultations
  "gfoclt" finds the times when one body is in front of another body as seen by an
  observer. The windows are those of CSPICE.

  Scenario Outline: An occultation of two ellipsoids
    Given a context with the standard kernels loaded
    And a front body, a back body and an observer
    When "gfoclt" searches for the occultation type <type>
    Then the window has the same number of intervals as the window of "gfoclt_c"
    And each end point agrees with CSPICE within twice the convergence tolerance
    And the window is a subset of the confinement window

    Examples:
      | type    |
      | ANY     |
      | PARTIAL |
      | FULL    |
      | ANNULAR |

  Scenario: A transit of a planet across the Sun
    Given a confinement window of some days about a transit of Mercury or of Venus
    When "gfoclt" searches for each occultation type
    Then "ANY" and "ANNULAR" give one interval
    And "PARTIAL" gives two intervals, the ingress and the egress
    And "FULL" gives no interval

  Scenario: A point target
    Given one of the two bodies with the shape "POINT"
    When "gfoclt" searches for the type "ANY"
    Then the window agrees with the window of "gfoclt_c"
    And the frame of the point is not read

  Scenario: A DSK shape in front of a point
    Given a DSK shape of the Moon and a point body behind it
    When "gfoclt" searches with the shape "DSK/UNPRIORITIZED", with and without a surface list
    Then the window agrees with the window of "gfoclt_c"

  Scenario: Stellar aberration in the correction
    Given a correction with "+S"
    When "gfoclt" searches
    Then the window is equal to the window for the correction without "+S"

  Scenario: Step function and refinement of the caller
    Given a step function and a refinement function
    When "gfocce" searches with them
    Then with a constant step and "gfrefn" the window is equal to the window of "gfoclt"
    And with another refinement each end point moves by no more than twice the tolerance

  Scenario: Forms, keywords and concurrent use
    Given the same search on a kernel set, on a context and with integer codes
    When each form runs, the search runs with "threads" set to true, and many tasks search at the same time
    Then every window is equal, bit for bit, to the window of the serial search

  Scenario Outline: Bad inputs
    Given a valid search with one input replaced by <input>
    When the search runs in Allspice and in CSPICE
    Then Allspice throws a SpiceError of kind <kind>
    And CSPICE signals the matching short error

    Examples:
      | input                                 | kind              |
      | an occultation type that is not known | invalidocctype    |
      | a body that is not known              | idcodenotfound    |
      | two bodies that are equal             | bodiesnotdistinct |
      | a shape that is not known             | badmethodsyntax   |
      | two points                            | invalidshapecombo |
      | a DSK shape with an ellipsoid         | invalidshapecombo |
      | a point with a type other than ANY    | badtypeshapecombo |
      | a frame not centered on the body      | invalidframe      |
      | a bad aberration correction           | invalidoption     |
      | an ellipsoid without radii            | kernelvarnotfound |
      | a step of zero                        | invalidstep       |
