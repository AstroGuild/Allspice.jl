Feature: Errors for queries that cannot be answered
  A query that the loaded kernels cannot answer throws a "SpiceError". The
  field "kind" tells the caller which condition occurred.

  Background:
    Given a context with "naif0012.tls", "pck00011.tpc" and "de440s.bsp" loaded

  Scenario Outline: Epoch outside the coverage
    When the user asks for the state of body 301 relative to body 399 at ET <epoch> seconds
    Then a "SpiceError" with kind "nocoverage" is thrown
    And the message is not empty

    Examples:
      | epoch      |
      | 4.736e9    |
      | -4.735e9   |

  Scenario: The ends of the coverage are inside the coverage
    When the user asks for a state at the first or the last epoch of the coverage
    Then the answer agrees with CSPICE

  Scenario: Unknown body name
    When the user asks for the state of "NO_SUCH_BODY" relative to "EARTH"
    Then a "SpiceError" with kind "idcodenotfound" is thrown, as in CSPICE

  Scenario: Unknown frame name
    When the user asks for a state in the frame "NO_SUCH_FRAME"
    Then a "SpiceError" with kind "unknownframe" is thrown

  Scenario: Aberration correction that is not valid
    When the user asks for a state with correction "RL+S"
    Then a "SpiceError" with kind "invalidoption" is thrown, as in CSPICE

  Scenario: Kernel file that does not exist
    When the user loads a path that is not a file
    Then a "SpiceError" with kind "nosuchfile" is thrown, as in CSPICE
    And the context is unchanged

  Scenario: Unknown frame ID
    When the user asks for a rotation to a frame ID that has no definition
    Then a "SpiceError" with kind "unknownframe" is thrown

  Scenario: Meta-kernel that lists a meta-kernel
    When the user loads a meta-kernel whose "KERNELS_TO_LOAD" names another meta-kernel
    Then a "SpiceError" with kind "recursiveloading" is thrown, as in CSPICE
    And the number of loaded kernels is unchanged

  Scenario Outline: The kind of an error is the short message of CSPICE
    When the user makes the call "<call>" with "<fault>"
    Then the kind of the "SpiceError" is the short message that CSPICE signals, in lower case

    Examples:
      | call   | fault                                  |
      | furnsh | a path that is not a file              |
      | furnsh | a blank path                           |
      | ldpool | a path that is not a file              |
      | eklef  | a path that is not a file              |
      | ckobj  | a path that is not a file              |
      | ckobj  | a text kernel                          |
      | ckcov  | a time system that is not known        |
      | spkezr | a body name that is not known          |
      | spkezr | a correction that is not valid         |
      | spkpos | a correction that is not valid         |
      | spkapo | a correction that is not valid         |
      | bodvrd | a body name that is not known          |
      | bodvrd | an item that is not in the kernel pool |
      | bodvrd | an item that holds strings             |

  Scenario: The kinds are registered
    When the user reads "Allspice.ERROR_KINDS"
    Then it holds every kind that the package can throw, in sorted order
    And a "SpiceError" with a kind that is not in it cannot be constructed
