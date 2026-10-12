Feature: Geometry finder, coordinates of a position vector
  "gfposc" searches for the times at which one coordinate of the position of a
  target, as seen by an observer in a given frame, satisfies a relation. Longitude
  and right ascension have a branch cut, so the search uses the sine and the
  cosine of the coordinate to find the events on each side of the cut.

  Background:
    Given a context with a leapseconds kernel, a planetary constants kernel and "de440s.bsp" loaded

  Scenario Outline: Each coordinate of each system
    Given the coordinate system "<system>" and each of its coordinates
    And a target, a frame, an aberration correction and a random confinement window
    When "gfposc" searches for each of the relations "<", "=", ">", "LOCMIN", "LOCMAX", "ABSMIN", "ABSMAX"
    Then the result has the same number of intervals as the result of "gfposc_c"
    And each end point agrees with the end point of CSPICE within 2e-6 seconds

    Examples:
      | system         |
      | RECTANGULAR    |
      | LATITUDINAL    |
      | RA/DEC         |
      | SPHERICAL      |
      | CYLINDRICAL    |
      | GEODETIC       |
      | PLANETOGRAPHIC |

  Scenario: A root where the coordinate is almost stationary
    Given a search whose end point differs from the end point of CSPICE by more than 2e-6 seconds
    Then the difference is less than 1 second
    And the coordinate has the same value at the two epochs, to rounding

  Scenario: An inequality checked without CSPICE
    Given the result of a search for "<" or ">"
    When the coordinate is computed at random epochs of the confinement window that are not near an end point
    Then an epoch is in the result exactly when the coordinate satisfies the relation

  Scenario Outline: A longitude across the branch cut
    Given a target whose <coordinate> crosses the branch cut at <cut> in the confinement window
    And a reference value on the cut, next to the cut, or far from it
    When "gfposc" searches for each relation
    Then the result has the same end points as the result of CSPICE

    Examples:
      | coordinate      | cut          |
      | longitude       | pi and -pi   |
      | right ascension | 0 and 2 pi   |

  Scenario: Planetographic longitude
    Given a body whose planetographic longitude is positive to the west
    When "gfposc" searches for a longitude
    Then the sense of the longitude is the one of the body, as in CSPICE

  Scenario Outline: An argument that is not valid
    When "gfposc" gets <argument>
    Then it throws a SpiceError

    Examples:
      | argument                                     |
      | a coordinate system that is not known        |
      | a coordinate that is not in the system       |
      | a frame that is not known                    |
      | a target that is the same as the observer    |
