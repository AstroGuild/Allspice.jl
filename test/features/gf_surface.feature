Feature: Geometry finder, coordinates of a surface point
  "gfsubc" finds the times when a coordinate of the sub-observer point on a target
  satisfies a relation. "gfsntc" does the same for the point where a ray from the
  observer meets the target. The windows are those of CSPICE.

  Scenario Outline: A coordinate of the sub-observer point
    Given a context with the standard kernels loaded
    And a target, its body-fixed frame and an observer
    When "gfsubc" searches for <relation> on the coordinate <coordinate> of the system <system>
    Then the window has the same number of intervals as the window of "gfsubc_c"
    And each end point agrees with CSPICE within twice the convergence tolerance
    And the window is a subset of the confinement window

    Examples:
      | system         | coordinate      | relation |
      | LATITUDINAL    | LONGITUDE       | <        |
      | GEODETIC       | LATITUDE        | =        |
      | RECTANGULAR    | Z               | LOCMAX   |
      | RA/DEC         | RIGHT ASCENSION | ABSMIN   |
      | PLANETOGRAPHIC | LONGITUDE       | >        |

  Scenario: The two methods of the sub-observer point
    Given the methods "Near point: ellipsoid" and "Intercept: ellipsoid"
    When "gfsubc" searches with each method, with any case and blanks
    Then the window is the window of CSPICE for that method
    And another method string gives the error "notsupported" at the first evaluation

  Scenario: A ray that always meets the target
    Given a ray that is fixed in a dynamic frame that follows the target
    When "gfsntc" searches on a confinement window
    Then the existence window is the confinement window contracted by 1 s plus the tolerance at each end
    And the window of the search agrees with the window of "gfsntc_c"

  Scenario: A ray that meets the target for a part of the time
    Given a ray that meets the target for some days of each month only
    When "gfsntc" searches on a confinement window of two months
    Then the existence window has several intervals
    And the window of the search is a subset of the existence window
    And the window of the search agrees with the window of "gfsntc_c"

  Scenario: A ray that never meets the target
    Given a ray that points away from the target
    When "gfsntc" searches
    Then the window is empty and no error is thrown

  Scenario: A coordinate that is constant to rounding
    Given the altitude of a surface point, or the radius of a point of a sphere
    When "gfsubc" or "gfsntc" searches on that coordinate
    Then the search ends without an error
    And the window is not compared with CSPICE, because it follows rounding noise

  Scenario: The general call
    Given the quantity "COORDINATE" of "gfevnt"
    When the vector definition is "SUB-OBSERVER POINT" or "SURFACE INTERCEPT POINT"
    Then the window is equal to the window of "gfsubc" or "gfsntc" for the same inputs
    And fewer than 3 values of the ray direction give the error "invalidcount"

  Scenario: Forms and keywords
    Given the same search on a kernel set, on a context and with integer codes
    When each form runs, and the search runs with the keyword "threads" set to true
    Then every window is equal, bit for bit
    And the callback "report" receives the existence pass of "gfsntc" first
    And the callback "bail" stops the search with an InterruptException

  Scenario Outline: Bad inputs
    Given a valid search with one input replaced by <input>
    When the search runs in NeoSpice and in CSPICE
    Then NeoSpice throws a SpiceError of kind <kind>
    And CSPICE signals the matching short error

    Examples:
      | input                              | kind              |
      | a target that is not known         | idcodenotfound    |
      | an observer equal to the target    | bodiesnotdistinct |
      | a frame that is not known          | noframe           |
      | a frame not centered on the target | invalidframe      |
      | a method that is not supported     | notsupported      |
      | a bad aberration correction        | invalidoption     |
      | a coordinate system not known      | notsupported      |
      | a relation that is not known       | notrecognized     |
      | a negative adjust value            | valueoutofrange   |
      | an adjust value with "="           | invalidvalue      |
      | a step of zero                     | invalidstep       |
      | a ray direction of zero            | zerovector        |
