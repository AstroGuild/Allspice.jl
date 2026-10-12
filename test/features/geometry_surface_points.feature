Feature: Points on the surface of a body
  A user asks for the point of a body below an observer or below the Sun, or for the
  point where a ray meets the surface. The body is a reference ellipsoid or a surface
  from DSK data. The light time is the one to the point, not to the center of the
  body, as in CSPICE, and the results agree with CSPICE to rounding level.

  Background:
    Given a context with the standard kernels, "de440s.bsp" and a DSK file for the Moon
    And CSPICE with the same kernels loaded

  Scenario Outline: Sub-observer point on an ellipsoid
    When the user calls "subpnt" with method "<method>" and correction "<abcorr>"
    Then the point, the epoch at the point and the vector from the observer agree with "subpnt_c"
    And the call does not allocate and its result type is inferred

    Examples:
      | method                | abcorr |
      | NEAR POINT/ELLIPSOID  | NONE   |
      | INTERCEPT/ELLIPSOID   | LT+S   |
      | Near point: ellipsoid | CN+S   |
      | NEAR POINT/ELLIPSOID  | XCN+S  |

  Scenario: Sub-solar point
    When the user calls "subslr" with a reception correction
    Then the results agree with "subslr_c"
    When the correction calls for transmission
    Then a "SpiceError" of kind "notsupported" is thrown

  Scenario Outline: Sub-points on DSK data
    When the user calls "subpnt" with method "<method>" for the Moon
    Then the point agrees with CSPICE for the surfaces that the method lists

    Examples:
      | method                                   |
      | NADIR/DSK/UNPRIORITIZED                  |
      | INTERCEPT/DSK/UNPRIORITIZED/SURFACES = 2 |

  Scenario: Ray that meets the surface
    Given a ray from the observer in a frame that is not the frame of the target
    When the user calls "sincpt"
    Then the intercept, its epoch and the vector from the observer agree with "sincpt_c"

  Scenario: Ray that misses the surface
    When the user calls "sincpt" with a ray that passes outside of the limb
    Then the result is "nothing", where CSPICE returns the found flag false

  Scenario: Epoch as a dual number
    When the user calls "subpnt", "subslr" or "sincpt" on an ellipsoid with a "ForwardDiff.Dual" epoch
    Then the derivative parts agree with a central difference
    When the method asks for DSK data
    Then a "SpiceError" of kind "notsupported" is thrown

  Scenario Outline: Inputs that are not valid
    When the user calls "subpnt" with <input>
    Then a "SpiceError" of kind "<kind>" is thrown, the short message of CSPICE in lower case

    Examples:
      | input                                    | kind              |
      | the correction "FOO"                     | invalidoption     |
      | a target name that is not known          | idcodenotfound    |
      | the observer equal to the target         | bodiesnotdistinct |
      | a frame that is not known                | noframe           |
      | a frame centered on another body         | invalidframe      |
      | a method without a sub-point type        | invalidsubtype    |
      | a body without radii in the kernel pool  | kernelvarnotfound |
