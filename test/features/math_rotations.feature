Feature: Rotation calls
  A user converts between rotation matrices, Euler angles, axis and angle, quaternions
  and state transformations with the CSPICE calls. The results agree with CSPICE to
  rounding level, with the same angle ranges and the same choice at the singular cases.

  Background:
    Given random rotations from Euler angles over the whole circle
    And the special rotations: the identity, the half turns, and Euler sequences with the middle angle at an end of its range
    And CSPICE called directly as the oracle

  Scenario Outline: Conversion of a rotation
    When the user calls "<call>"
    Then the result agrees with "<call>_c" of CSPICE to rounding level
    And "<inverse>" of the result gives the input back

    Examples:
      | call   | inverse |
      | eul2m  | m2eul   |
      | m2eul  | eul2m   |
      | axisar | raxisa  |
      | raxisa | axisar  |
      | m2q    | q2m     |
      | q2m    | m2q     |
      | eul2xf | xf2eul  |
      | xf2eul | eul2xf  |
      | rav2xf | xf2rav  |
      | xf2rav | rav2xf  |

  Scenario: Euler angles at the singular middle angle
    Given a rotation whose middle Euler angle is 0 or pi (a-b-a sequence) or plus or minus pi/2 (a-b-c sequence)
    When the user calls "m2eul"
    Then the first angle is 0 and the third angle carries the whole rotation, as in CSPICE

  Scenario: Unique Euler angles from a state transformation
    When the user calls "xf2eul"
    Then the flag "unique" is false only at the singular middle angle

  Scenario: Quaternion of a half turn
    When the user calls "m2q" with a half turn
    Then the scalar part is 0 and the sign of the vector part is the one of CSPICE

  Scenario: Matrix that is not a rotation
    When the user calls "m2eul", "m2q" or "raxisa" with a matrix that is not a rotation
    Then a SpiceError with the kind "notarotation" is thrown

  Scenario: Axis number outside 1 to 3
    When the user calls "eul2m" or "eul2xf" with an axis number outside 1 to 3
    Then a SpiceError with the kind "badaxisnumbers" is thrown

  Scenario: Product of quaternions and angular velocity
    When the user calls "qxq" and "qdq2av"
    Then the product agrees with CSPICE bit for bit
    And the angular velocity agrees with the derivative of the rotation by automatic differentiation

  Scenario: Sharpening a matrix
    When the user calls "sharpr" with a matrix that is almost a rotation
    Then the result is a rotation and agrees with "sharpr_" of CSPICE

  Scenario: Static arguments
    When the user calls a rotation call with static arrays
    Then the result type is inferred and the call does not allocate
