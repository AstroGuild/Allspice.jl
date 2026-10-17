Feature: Two-line elements and SPK segments of type 10
  A user has two-line element sets of Earth satellites. The user turns the two lines
  into element numbers, propagates them with SGP4 or SDP4, writes them to an SPK
  segment of type 10 and asks for states from that segment. Each result agrees with
  CSPICE.

  Background:
    Given a new context with a leap-second kernel and a planetary constants kernel
    And a text kernel that sets BODY399_J2, BODY399_J3, BODY399_J4, BODY399_KE, BODY399_QO, BODY399_SO, BODY399_ER and BODY399_AE
    And the 29 verification element sets of the SGP4 revision of 2006

  Scenario: Geophysical constants from a kernel
    When the user reads the eight constants of body 399 with bodvrd
    Then the values are equal to the values that CSPICE reads

  Scenario: Element numbers from the two lines
    When the user calls getelm with the first year 1957 and the two lines of a set
    Then the epoch and the ten element numbers are equal to the ones of CSPICE
    And the tenth element number is the epoch

  Scenario: Lines with a wrong checksum
    Given two lines whose last digit (the checksum) is wrong
    When the user calls getelm
    Then the result is the same as for the correct lines, as in CSPICE N0067

  Scenario Outline: Lines that getelm refuses
    Given two lines with <defect>
    When the user calls getelm
    Then a SpiceError of the kind badtle is thrown
    And CSPICE signals an error for the same lines

    Examples:
      | defect                                         |
      | different catalogue numbers in the two lines   |
      | a line that is shorter than 68 characters      |
      | a term that is not a number                    |
      | an inclination outside 0 to 180 degrees        |
      | a node, perigee or anomaly outside 0 to 360    |
      | a mean motion of 20 revolutions a day or more  |

  Scenario Outline: Propagation of one element set
    Given the element numbers of a set of the class <class>
    When the user calls evsgp4 at an epoch within two days of the element epoch
    Then the state in the TEME frame agrees with CSPICE to within rounding
    And the call does not allocate memory

    Examples:
      | class             |
      | low orbit         |
      | 12 h resonant     |
      | 1 day resonant    |
      | semi-synchronous  |
      | high eccentricity |
      | deep space        |

  Scenario: Element set of a decaying object
    Given the element numbers of a set with a low perigee and strong drag
    When the user calls evsgp4 at an epoch where the model has no valid orbit
    Then a SpiceError is thrown whose kind is the short message of CSPICE in lower case
    And CSPICE signals its error at the same epochs

  Scenario: State from a type 10 segment written by CSPICE
    Given an SPK file with one type 10 segment of three element sets, written by spkw10_c
    When the user asks for the state of the object relative to the Earth in J2000 or in IAU_EARTH
    Then the state agrees with CSPICE to within rounding
    And the query does not allocate memory

  Scenario: Epoch between two element sets
    When the user asks for the state at an epoch between the epochs of two element sets
    Then the state is the weighted mean of the propagations of the two sets, as in CSPICE
    But at the middle between the two epochs the later pair of sets is used, as in CSPICE

  Scenario: Epoch before the first or after the last element set
    When the user asks for the state at a segment end that is outside the element epochs
    Then the state is the propagation of the end set

  Scenario: Segment written by Allspice
    When the user writes the same element sets with spkw10 of Allspice
    Then the file has the size and the layout of the file of CSPICE
    And only stored nutation numbers differ, by rounding
    And each engine gives the same states from either file

  Scenario: Segment with more than 100 element sets
    When the user writes 250 element sets to one segment
    Then the segment has a directory of every 100th epoch
    And the states agree with CSPICE

  Scenario: Segment in the old layout
    Given a type 10 segment with packets of ten numbers and 15 meta items
    When the user asks for a state
    Then the state agrees with CSPICE

  Scenario: Time derivatives of a type 10 state
    When the user calls spkder with the order 2 or 3
    Then the orders 0 and 1 are equal to the position and the velocity of spkgeo
    And the order 2 agrees with a central difference of the velocity

  Scenario: No leap-second kernel
    Given a context with the SPK file only
    When the user asks for a state of the type 10 segment
    Then a SpiceError of the kind noleapseconds is thrown

  Scenario Outline: Segment that spkw10 refuses
    When the user calls spkw10 with <defect>
    Then a SpiceError of the kind <kind> is thrown
    And nothing is written to the file

    Examples:
      | defect                              | kind              |
      | no element set                      | numpacketsnotpos  |
      | epochs that do not increase         | unorderedrefs     |
      | a segment name of 41 characters     | segidtoolong      |
      | a tab character in the segment name | nonprintingchars  |
      | the body equal to the center        | bodyandcentersame |
      | a frame name that is not known      | invalidrefframe   |
      | the first epoch after the last one  | baddescrtimes     |

  Scenario: A subset of a type 10 segment
    Given an SPK file of type 10 that CSPICE wrote
    When spksub writes the subset of the segment for a time window
    Then the new file is equal byte for byte to the file of spksub_c
    And spkpvn gives the state of CSPICE within rounding
