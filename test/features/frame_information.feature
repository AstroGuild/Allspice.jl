Feature: Frame information, body constants and local solar time
  A user asks which frames exist and how they are defined, reads body constants by
  body ID, asks for the rotation of a body, for the coverage of a binary PCK, and for
  the local solar time on a body. The answers are those of CSPICE.

  Background:
    Given a context and CSPICE with text kernels, a binary PCK, "de440s.bsp" and frame kernels loaded

  Scenario: Attributes of a frame
    When the user asks for "frinfo" of each of 234 frame IDs
    Then center, class and class ID are equal to CSPICE, or both report "not found"

  Scenario Outline: Frame of a body
    When the user asks for "<call>" of "<body>"
    Then the frame ID and name are equal to CSPICE

    Examples:
      | call   | body  |
      | cidfrm | 399   |
      | cidfrm | 301   |
      | cidfrm | 0     |
      | cnmfrm | EARTH |
      | cnmfrm | earth |
      | cnmfrm | TITAN |

  Scenario: A frame kernel names the frame of a body
    Given a kernel that assigns "OBJECT_EARTH_FRAME" and "OBJECT_599_FRAME"
    When the user asks for "cidfrm" of 399 and of 599
    Then the frames of the kernel are returned, as in CSPICE

  Scenario: Frame from class and class ID
    When the user asks for "ccifrm" with the class and class ID of each loaded frame
    Then ID, name and center are equal to CSPICE

  Scenario: Data of a fixed-offset frame
    When the user asks for "tkfram" of a class 4 frame
    Then the matrix and the relative frame are equal to "tkfram_" of CSPICE

  Scenario: Rotation of a body
    When the user asks for "tipbod" and "tisbod" of 13 bodies in 7 inertial frames
    Then the matrices agree with CSPICE to rounding
    And the binary PCK is used for the Moon

  Scenario: Contents of a binary PCK
    When the user asks for "pckfrm" and "pckcov" of the lunar binary PCK
    Then the frame class IDs and the coverage window are equal to CSPICE

  Scenario: Body constants by ID
    When the user asks for "bodfnd" and "bodvcd" of 11 bodies and 11 items
    Then the answers are equal to CSPICE

  Scenario Outline: Local solar time
    When the user asks for "et2lst" on body <body> at a random epoch and longitude of type "<type>"
    Then hour, minute, second and both strings are equal to CSPICE

    Examples:
      | body | type           |
      | 399  | PLANETOCENTRIC |
      | 399  | PLANETOGRAPHIC |
      | 301  | PLANETOGRAPHIC |
      | 299  | PLANETOGRAPHIC |
      | 199  | PLANETOCENTRIC |
      | 10   | PLANETOCENTRIC |

  Scenario: A kernel variable sets the sense of a planetographic longitude
    Given the variable "BODY299_PGR_POSITIVE_LON" with the value "EAST" or "WEST"
    When the user asks for "et2lst" of Venus with a "PLANETOGRAPHIC" longitude
    Then the local time is equal to the local time of CSPICE
    But another value of the variable throws a "SpiceError" of the kind "invalidoption"

  Scenario: Longitude of the Sun
    When the user asks for "lspcn" of the Earth, the Moon, Venus and Mercury with five corrections
    Then the angle agrees with CSPICE to 2e-13 rad

  Scenario: Built-in frames and frames of the kernel pool
    When the user calls "bltfrm" and "kplfrm" for each frame class and for all classes
    Then the frame IDs are those of CSPICE, sorted and without duplicates
    And a class outside -1 and 1 to 6 throws a "SpiceError" of kind "badframeclass"
