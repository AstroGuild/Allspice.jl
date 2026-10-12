Feature: Unload kernels
  A user removes a kernel that was loaded. Every later answer is the one that
  the remaining kernels give.

  Scenario: Unload restores the earlier answers
    Given a context with "de440s.bsp" loaded
    And the state of the Moon relative to the Sun is recorded at 20 epochs
    When the user loads "de405.bsp"
    Then the states at the same epochs are different
    When the user unloads "de405.bsp"
    Then the states at the same epochs are equal, bit for bit, to the recorded states

  Scenario: Unload a text kernel
    Given a context with "pck00011.tpc" and then "pck00010.tpc" loaded
    When the user unloads "pck00010.tpc"
    Then "bodvrd" of "MARS", "PM" gives the values of "pck00011.tpc"

  Scenario: Unload a meta-kernel
    Given a context in which a meta-kernel loaded four kernels
    When the user unloads the meta-kernel
    Then the four kernels are unloaded too
    And "ktotal" is equal to the count before the meta-kernel was loaded

  Scenario: Unload a kernel that is not loaded
    Given a context with "de440s.bsp" loaded
    When the user unloads a path that was never loaded
    Then no error occurs
    And all answers are unchanged

  Scenario: Clear all kernels
    Given a context with several kernels loaded
    When the user calls "kclear"
    Then "ktotal" is 0
    And a state query throws a "SpiceError"
    And built-in body names and inertial frames still resolve
