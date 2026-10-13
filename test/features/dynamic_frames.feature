Feature: Dynamic frames and switch frames
  A frame kernel can define a frame whose orientation is computed from ephemerides,
  from models or from polynomials (class 5), or a frame that is aligned with one of a
  list of base frames (class 6). A user asks for rotations and states in such a frame
  with the same calls as for any other frame, and the answers agree with CSPICE.

  Background:
    Given a context with a leapseconds kernel, a text PCK and "de440s.bsp" loaded
    And a frame kernel that defines dynamic frames and switch frames
    And CSPICE with the same kernels loaded in the same order

  Scenario Outline: Rotation and state in a dynamic frame
    When the user asks for "pxform" from "J2000" to "<frame>" at a random epoch
    Then each element of the matrix agrees with CSPICE to rounding
    When the user asks for "sxform" from "J2000" to "<frame>" at the same epoch
    Then the derivative block agrees with CSPICE to the accuracy of the numerical derivative of CSPICE
    When the user asks for "spkezr" of the Moon relative to the Earth in "<frame>" with "<correction>"
    Then position and velocity agree with CSPICE

    Examples:
      | frame            | family                                | correction |
      | DT_GSE           | two-vector, position and velocity     | NONE       |
      | DT_GSM           | two-vector, constant vector           | LT+S       |
      | DT_NADIR_LT      | two-vector, near point                | LT+S       |
      | DT_EME_OF_DATE   | mean equator and equinox of date      | NONE       |
      | DT_ETE_OF_DATE   | true equator and equinox of date      | LT+S       |
      | DT_ECL_OF_DATE   | mean ecliptic and equinox of date     | NONE       |
      | DT_EULER         | Euler angles                          | CN+S       |
      | DT_PRODUCT       | product of frame transformations      | LT+S       |
      | DT_EULER_ON_GSE  | dynamic frame on a dynamic base       | XLT+S      |
      | DT_TWOVEC_IN_DYN | vectors given in dynamic frames       | NONE       |

  Scenario: A frozen frame does not move
    Given the frame "DT_GSE_FROZEN" with a "FREEZE_EPOCH"
    When the user asks for "pxform" from "J2000" to "DT_GSE_FROZEN" at 20 epochs
    Then every matrix is equal to the matrix of "DT_GSE" at the freeze epoch
    And every derivative of "pxder" is zero

  Scenario: A frame with an inertial rotation state
    Given the frame "DT_ETE_INERTIAL" with "ROTATION_STATE" equal to "INERTIAL"
    When the user asks for "sxform" from "J2000" to "DT_ETE_INERTIAL"
    Then the derivative block is zero, as in CSPICE
    And "pxform" is equal to "pxform" of the rotating frame "DT_ETE_OF_DATE"

  Scenario: Exact derivatives
    When the user asks for "pxder" from "J2000" to a dynamic frame with order 2
    Then order 0 is the rotation of "pxform"
    And order 1 agrees with a finite difference of "pxform"
    And order 2 agrees with a finite difference of order 1

  Scenario: A switch frame selects its base by time
    Given the switch frame "DT_SWITCH_TIME" with four bases and four time intervals
    When the user asks for "pxform" at an epoch inside two intervals
    Then the base that is listed last is used
    When the user asks at an epoch outside every interval
    Then the error "noframeconnect" occurs, as in CSPICE

  Scenario: A switch frame with frame IDs and time strings
    Given a switch frame whose bases are frame ID codes and whose bounds are time strings
    When the user asks for "pxform" at an epoch inside one interval, and inside two intervals
    Then the matrix is that of the base of the last interval that holds the epoch, as in CSPICE
    When the user asks at an epoch outside every interval
    Then the error "noframeconnect" occurs, as in CSPICE

  Scenario Outline: Angle units of an Euler frame
    Given an Euler frame whose angles have the unit "<unit>"
    When the user asks for "sxform" from "J2000" to that frame
    Then the rotation and the derivative block agree with CSPICE

    Examples:
      | unit        |
      | ARCMINUTES  |
      | ARCSECONDS  |
      | HOURANGLE   |
      | MINUTEANGLE |
      | SECONDANGLE |

  Scenario: A switch frame over CK frames follows the CK coverage
    Given a switch frame whose bases are a PCK frame and two CK frames
    When the user asks for "pxform" at an epoch at which only the first CK frame has pointing
    Then the first CK frame is used
    When the user asks at an epoch at which no CK frame has pointing
    Then the PCK frame is used

  Scenario: The recursion limit of CSPICE
    Given a two-vector frame whose velocity is given in a frame that has vectors in dynamic frames
    When the user asks for "pxform" to that frame
    Then the error "recursiontoodeep" occurs, as in CSPICE

  Scenario Outline: A definition that is not valid
    Given a frame kernel in which the frame "<frame>" has the defect "<defect>"
    When the user loads the kernel
    Then no error occurs
    When the user asks for "pxform" to "<frame>"
    Then the error "<error>" occurs, which is the short message of CSPICE
    And the other frames of the kernel give the same answers as before

    Examples:
      | frame           | defect                               | error             |
      | DT_NOFAMILY     | no FAMILY keyword                    | kernelvarnotfound |
      | DT_BADFAMILY    | unknown family                       | notsupported      |
      | DT_SAMEAXIS     | both vectors define the same axis    | undefinedframe    |
      | DT_PARALLEL     | the two vectors are parallel         | degeneratecase    |
      | DT_BADBODY      | unknown body name                    | notranslation     |
      | DT_BADABCORR    | unknown aberration correction        | spkinvalidoption  |
      | DT_SWITCH_COUNT | interval count differs from bases    | countmismatch     |
      | DT_BADTYPE      | a number where a string is needed    | badvariabletype   |
      | DT_BADSIZE      | two values where one is allowed      | badvariablesize   |

  Scenario: No allocation
    When the user calls "pxform", "sxform", "pxder" or "spkezr" with integer IDs for a dynamic frame
    Then the call allocates no memory
    And the return type is inferred
