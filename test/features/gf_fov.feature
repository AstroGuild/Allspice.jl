Feature: Geometry finder, visibility in a field of view
  "gftfov" finds the times when a target body is in the field of view of an
  instrument. "gfrfov" does the same for a ray. The windows are those of CSPICE.

  Scenario Outline: A target body in a field of view
    Given a context with the standard kernels and an instrument kernel loaded
    And an instrument with a field of view of the shape <shape>
    When "gftfov" searches for a target with the shape <target shape>
    Then the window has the same number of intervals as the window of "gftfov_c"
    And each end point agrees with CSPICE within twice the convergence tolerance
    And the window is a subset of the confinement window

    Examples:
      | shape     | target shape |
      | CIRCLE    | ELLIPSOID    |
      | ELLIPSE   | ELLIPSOID    |
      | RECTANGLE | POINT        |
      | POLYGON   | ELLIPSOID    |
      | ELLIPSE   | POINT        |

  Scenario: A camera that turns with the Earth
    Given a camera frame that is fixed to the Earth and points near the Sun at one epoch
    When "gftfov" searches on a confinement window of one week
    Then the window has several intervals, about one for each day
    And the window agrees with the window of "gftfov_c"

  Scenario Outline: A ray in a field of view
    Given a ray with a direction that is fixed in the frame <frame>
    When "gfrfov" searches with the correction <correction>
    Then the window agrees with the window of "gfrfov_c"

    Examples:
      | frame     | correction |
      | J2000     | NONE       |
      | J2000     | S          |
      | IAU_EARTH | XS         |

  Scenario: Step function and refinement of the caller
    Given a step function and a refinement function
    When "gffove" searches with them for a body or for a ray
    Then with a constant step and "gfrefn" the window is equal to the window of "gftfov" or "gfrfov"
    And with another refinement each end point moves by no more than twice the tolerance

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
      | input                                  | kind              |
      | an instrument that is not known        | idcodenotfound    |
      | a body that has no field of view       | framemissing      |
      | a target equal to the observer         | bodiesnotdistinct |
      | a target shape that is not known       | invalidshape      |
      | the shape RAY in "gftfov"              | invalidoption     |
      | a blank frame of an ellipsoid          | invalidframe      |
      | a frame not centered on the target     | invalidframe      |
      | a light-time correction with a ray     | invalidoption     |
      | a ray direction of zero                | zerovector        |
      | a field of view 90 degrees wide        | fovtoowide        |
      | a step of zero                         | invalidstep       |
