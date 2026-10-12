Feature: Ellipsoids, planes and ellipses
  A user asks geometric questions about a triaxial ellipsoid, about planes and about
  ellipses with the CSPICE calls. Planes and ellipses are the immutable values
  "SpicePlane" and "SpiceEllipse". A call that has a "found" flag in CSPICE returns
  "nothing" when the flag is false. The results agree with CSPICE, most of them bit
  for bit.

  Background:
    Given ellipsoids of sizes from 1e-10 to 1e10 with axis ratios up to 10, among them spheres and spheroids
    And points outside, inside, on the axes, in the principal planes and at the center
    And CSPICE called directly as the oracle

  Scenario: Nearest point on an ellipsoid
    When the user calls "nearpt" for a point
    Then the near point and the altitude are equal to those of "nearpt_c"
    And the near point is on the ellipsoid
    And the altitude is negative for a point inside

  Scenario: Nearest point at the center
    When the user calls "nearpt" for the center of the ellipsoid
    Then the near point is on the shortest semi-axis, as in CSPICE

  Scenario: Derivative of the nearest point
    When the user calls "nearpt" with ForwardDiff or Taylor numbers
    Then the plain value is the one of the CSPICE iteration
    And the derivative agrees with the velocity that "dnearp" gives
    And the derivative agrees with a central difference

  Scenario: Ray and ellipsoid
    When the user calls "surfpt" with a ray that meets the ellipsoid
    Then the point is equal to the point of "surfpt_c"
    And a ray that misses gives "nothing"
    And "surfpv" gives the velocity of the point, which is the time derivative of "surfpt"

  Scenario: Semi-axis that is not positive
    When the user calls an ellipsoid call with a semi-axis that is zero or negative
    Then a SpiceError is thrown with the kind that CSPICE signals for that call

  Scenario Outline: Plane from its description and back
    When the user calls "<make>"
    Then the plane has a unit normal and a constant that is not negative
    And the plane is equal to the plane of CSPICE
    And "<back>" gives a description of the same plane

    Examples:
      | make   | back   |
      | nvc2pl | pl2nvc |
      | nvp2pl | pl2nvp |
      | psv2pl | pl2psv |

  Scenario Outline: Ray and plane
    When the user calls "inrypl" with a ray that "<relation>" the plane
    Then the count of intersection points is <count>

    Examples:
      | relation          | count |
      | crosses           | 1     |
      | points away from  | 0     |
      | lies in           | -1    |

  Scenario: Ellipse from generating vectors
    When the user calls "cgv2el" or "saelgv" with two generating vectors
    Then the semi-axes are orthogonal and the major one is not the shorter one
    And they are equal to the semi-axes of CSPICE

  Scenario: Ellipsoid and plane
    When the user calls "inedpl" with a plane that cuts the ellipsoid
    Then the points of the ellipse are on the ellipsoid and in the plane
    And a plane that misses the ellipsoid gives "nothing"

  Scenario Outline: Ellipse and plane
    When the user calls "inelpl" with a plane that meets the ellipse in <count> points
    Then the count is <count> and the points are those of CSPICE

    Examples:
      | count |
      | 0     |
      | 1     |
      | 2     |

  Scenario: Limb of an ellipsoid
    When the user calls "edlimb" with a view point outside the ellipsoid
    Then the line of sight is tangent to the ellipsoid at each limb point
    And a view point inside gives a SpiceError with the kind "degeneratecase"

  Scenario: Nearest points to a line and to an ellipse
    When the user calls "npedln", "nplnpt" or "npelpt"
    Then the near point and the distance are equal to those of CSPICE

  Scenario Outline: Terminator of a body
    Given a context with the standard kernels and "de440s.bsp" loaded
    When the user calls "edterm" for the "<type>" terminator of the Moon lit by the Sun
    Then the epoch, the observer position and the points agree with "edterm_c"
    And each point is on the ellipsoid of the Moon

    Examples:
      | type      |
      | UMBRAL    |
      | PENUMBRAL |

  Scenario: Terminator with arguments that are not valid
    When the user calls "edterm" with a terminator type, a frame, a body or a correction that is not valid
    Then a SpiceError is thrown with the kind that CSPICE signals

  Scenario: Static arguments
    When the user calls an ellipsoid, plane or ellipse call with static arrays
    Then the result type is inferred and the call does not allocate
