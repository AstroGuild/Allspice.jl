Feature: Instrument fields of view
  A user reads the field of view of an instrument from the kernel pool, and asks
  whether a body or a ray is in it at an epoch.

  Background:
    Given a context with the standard kernels, "de440s.bsp" and an instrument kernel
    And CSPICE with the same kernels loaded

  Scenario Outline: Read a field of view
    When the user calls "getfov" for an instrument with the shape "<shape>" given by "<spec>"
    Then the shape, the frame, the boresight and the boundary vectors are those of "getfov_c"

    Examples:
      | shape     | spec    |
      | CIRCLE    | CORNERS |
      | ELLIPSE   | CORNERS |
      | RECTANGLE | CORNERS |
      | POLYGON   | CORNERS |
      | CIRCLE    | ANGLES  |
      | ELLIPSE   | ANGLES  |
      | RECTANGLE | ANGLES  |

  Scenario: Read a field of view by name
    When the user calls "getfvn" with the name of the instrument
    Then the result is that of "getfov" with its ID

  Scenario: Definition that is not complete
    When a variable of the definition is absent or has a wrong size
    Then a "SpiceError" is thrown, of the kind that CSPICE gives

  Scenario Outline: Body in a field of view
    Given the instrument frame points near the target
    When the user calls "fovtrg" with the target shape "<tshape>"
    Then the flag is that of "fovtrg_c"

    Examples:
      | tshape    |
      | ELLIPSOID |
      | POINT     |

  Scenario Outline: Ray in a field of view
    When the user calls "fovray" with correction "<abcorr>"
    Then the flag is that of "fovray_c"

    Examples:
      | abcorr |
      | NONE   |
      | S      |
      | XS     |

  Scenario: Field of view that is too wide
    When the boundary vectors are 90 degrees from the axis
    Then a "SpiceError" of kind "fovtoowide" is thrown
