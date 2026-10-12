Feature: Pointing from C-kernels (CK types 1 to 6)
  A C-kernel holds the orientation of an instrument or a spacecraft structure,
  with or without angular velocity, as a function of spacecraft clock time.

  Background:
    Given a context with a leapseconds kernel, an SCLK kernel and a frame kernel loaded
    And CK files with segments of the types 1, 2, 3, 4, 5 (subtypes 0 to 3) and 6

  Scenario Outline: Pointing at a requested clock time
    Given a segment of type <type> for an instrument
    When ckgp is called with a clock time inside the coverage of the segment
    Then the C-matrix is equal to the CSPICE C-matrix to rounding accuracy
    And the output clock time is equal to the CSPICE output clock time

    Examples:
      | type |
      | 1    |
      | 2    |
      | 3    |
      | 4    |
      | 5    |
      | 6    |

  Scenario: Pointing with angular velocity
    Given a segment that holds angular velocity
    When ckgpav is called inside the coverage
    Then the C-matrix and the angular velocity are equal to the CSPICE answers

  Scenario: A segment without angular velocity
    Given the only segment that covers the request holds no angular velocity
    When ckgpav is called
    Then the found flag is false
    And ckgp for the same request finds the pointing

  Scenario: Tolerance
    Given a request in a gap of the data
    When the tolerance is smaller than the distance to the nearest pointing
    Then the found flag is false
    When the tolerance is larger than that distance
    Then the pointing at the nearest covered time is returned

  Scenario: Priority
    Given two segments for one instrument that cover the same time
    When pointing is requested
    Then the segment of the file loaded last is used
    And in one file the segment written last is used
    And the answer does not depend on the requests made before

  Scenario: Coverage
    Given a CK file
    When ckobj and ckcov are called
    Then the instruments and the coverage windows are equal to the CSPICE answers
    And the window is given in clock ticks or in ephemeris time as requested

  Scenario: An invalid request
    Given a negative tolerance
    When ckgp is called
    Then a SpiceError is thrown
