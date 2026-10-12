Feature: Limb, terminator and tangent points
  A user asks for points of the limb of a body as an observer sees it, for points of
  its terminator, or for the point of a ray that is nearest to the body. The results
  agree with CSPICE: the same number of points in each cutting half plane, and the same
  points to rounding level (on DSK data, to the solution tolerance).

  Background:
    Given a context with the standard kernels, "de440s.bsp" and DSK files for the Moon
    And CSPICE with the same kernels loaded

  Scenario Outline: Limb points
    When the user calls "limbpt" with method "<method>" and locus "<corloc>"
    Then the counts, the points, the epochs and the tangent vectors agree with "limbpt_c"

    Examples:
      | method                    | corloc         |
      | TANGENT/ELLIPSOID         | CENTER         |
      | GUIDED/ELLIPSOID          | CENTER         |
      | TANGENT/ELLIPSOID         | ELLIPSOID LIMB |
      | TANGENT/DSK/UNPRIORITIZED | CENTER         |
      | GUIDED/DSK/UNPRIORITIZED  | CENTER         |
      | TANGENT/DSK/UNPRIORITIZED | ELLIPSOID LIMB |

  Scenario Outline: Terminator points
    When the user calls "termpt" with method "<method>" and locus "<corloc>"
    Then the counts, the points, the epochs and the vectors agree with "termpt_c"

    Examples:
      | method                              | corloc               |
      | UMBRAL/TANGENT/ELLIPSOID            | CENTER               |
      | PENUMBRAL/GUIDED/ELLIPSOID          | CENTER               |
      | UMBRAL/TANGENT/ELLIPSOID            | ELLIPSOID TERMINATOR |
      | PENUMBRAL/TANGENT/DSK/UNPRIORITIZED | CENTER               |

  Scenario: Several points in one half plane
    Given a DSK surface of two separate parts
    When a cutting half plane passes through the two parts
    Then "limbpt" and "termpt" return three points in it, as CSPICE does

  Scenario: More points than the room
    When the points of "limbpt" do not fit in "maxn"
    Then a "SpiceError" of kind "outofroom" is thrown

  Scenario Outline: Tangent point of a ray
    When the user calls "tangpt" with locus "<corloc>"
    Then the tangent point, the altitude, the range, the surface point, the epoch and the vector agree with "tangpt_c"
    And the call does not allocate

    Examples:
      | corloc        |
      | TANGENT POINT |
      | SURFACE POINT |

  Scenario: Ray that points away from the body
    When the ray of "tangpt" points away from the target
    Then the tangent point is the observer and the range is zero, as in CSPICE
