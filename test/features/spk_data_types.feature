Feature: States from every SPK data type except two-line elements
  A user loads an SPK kernel whose segments are of a data type other than the
  Chebyshev types 2 and 3. The states agree with CSPICE, and the time derivatives
  of every order are available as for the Chebyshev types.

  Background:
    Given a new context
    And an SPK kernel with one segment of each of the types 1, 5, 8, 9, 12, 13, 14, 15, 17, 18, 19, 20 and 21, written by CSPICE

  Scenario Outline: Geometric state of a segment of one data type
    When the user asks for the state of the target of the type <type> segment relative to its center at an epoch inside the segment
    Then each component agrees with CSPICE to within rounding
    And the query does not allocate memory

    Examples:
      | type |
      | 1    |
      | 5    |
      | 8    |
      | 9    |
      | 12   |
      | 13   |
      | 14   |
      | 15   |
      | 17   |
      | 18   |
      | 19   |
      | 20   |
      | 21   |

  Scenario: Epoch on a boundary of the data
    When the user asks for the state at a segment end, at the epoch of a stored state, or at the bound between two records
    Then Allspice selects the same record or window as CSPICE
    And the state agrees with CSPICE to within rounding

  Scenario: Bound between two intervals of a type 19 segment
    Given a type 19 segment whose selection flag is "later"
    When the user asks for the state at an epoch that is the bound between two interpolation intervals
    Then the state comes from the later interval
    But with the flag "earlier" the state comes from the earlier interval

  Scenario: Higher time derivatives
    When the user asks for the derivatives to order 3 with "spkder"
    Then order 0 is the CSPICE position and order 1 is the CSPICE velocity
    And order 2 and order 3 are the first and the second derivative of the CSPICE velocity function

  Scenario: Two-line elements
    Given a segment of type 10
    When the user asks for its state
    Then Allspice throws an error of kind "unsupported"

  Scenario: Damaged segment
    Given a segment whose count at the end does not agree with its length
    When the user asks for its state
    Then Allspice throws an error of kind "baddaf"
