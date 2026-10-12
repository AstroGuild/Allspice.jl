Feature: Geometry finder, solver core
  A search finds the time window in which a scalar quantity satisfies a relation
  with a reference value. The search steps through the confinement window to
  bracket each change of state, and then refines each bracket by bisection until
  it is narrower than the convergence tolerance. The step, the tolerance and the
  callbacks are arguments of the search; there is no global state.

  Scenario Outline: A user-defined scalar function
    Given a scalar function of time and a function that tells if it is decreasing
    And a confinement window of one or of several intervals
    When "gfuds" searches for the relation "<relation>" with a step smaller than the shortest event
    Then the result has the same end points as "gfuds_c" of CSPICE, bit for bit
    And every end point of the result is in the confinement window

    Examples:
      | relation |
      | <        |
      | =        |
      | >        |
      | LOCMIN   |
      | LOCMAX   |
      | ABSMIN   |
      | ABSMAX   |

  Scenario: An absolute extremum with an adjust value
    Given the relation "ABSMAX" and an adjust value greater than zero
    When the search runs
    Then the result is the window in which the quantity is within the adjust value of its absolute maximum

  Scenario: A user-defined boolean function
    Given a boolean function of time
    When "gfudb" searches the confinement window
    Then the result is the window in which the function is true
    And it has the same end points as "gfudb_c" of CSPICE, bit for bit

  Scenario: The convergence tolerance
    Given a context with the default tolerance of 1e-6 seconds
    When "gfstol" sets a different tolerance on that context
    Then the searches on that context use the new tolerance
    And the searches on other contexts do not change
    And the keyword "tol" of one search has precedence over the tolerance of the context

  Scenario: A tolerance that is not positive
    When "gfstol" gets a value that is zero or negative
    Then it throws a SpiceError of kind "invalidtolerance"

  Scenario: Progress report and interrupt
    Given a search with the keywords "report" and "bail"
    When the search runs
    Then "report" gets one start call, update calls with epochs that do not decrease in a pass, and one end call
    And the search throws an InterruptException if "bail" returns true

  Scenario Outline: An argument that is not valid
    When a search gets <argument>
    Then it throws a SpiceError of kind "<kind>"

    Examples:
      | argument                              | kind           |
      | a step that is zero or negative       | invalidstep    |
      | a relation that is not known          | notrecognized  |
      | a negative adjust value               | valueoutofrange |
      | an adjust value with the relation "=" | invalidvalue   |

  Scenario: A search that needs the geometry of a surface or of a field of view
    When "gfsubc", "gfsntc", "gfilum", "gfoclt", "gftfov" or "gfrfov" is called
    Then it throws a SpiceError of kind "unsupported"
