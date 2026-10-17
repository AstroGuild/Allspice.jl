Feature: Windows of time intervals
  A user does coverage arithmetic with windows: sets of disjoint closed intervals.
  Every result is the window that CSPICE gives for the same input.

  Scenario: Insert intervals in any order
    Given an empty window
    When the user inserts the intervals [5, 6], [1, 2] and [2, 5] with "wninsd!"
    Then the window holds the one interval [1, 6]
    And "card" is 2 and "wncard" is 1

  Scenario: Set operations
    Given two windows with random intervals, some of them singletons or touching
    When the user calls "wnunid", "wnintd" and "wndifd"
    Then each result is equal, endpoint for endpoint, to the result of CSPICE
    And the two input windows are unchanged

  Scenario: Complement in an interval
    Given a window and an interval [left, right]
    When the user calls "wncomd"
    Then the result holds the parts of [left, right] that the window does not cover

  Scenario: Queries
    Given a window
    When the user calls "wnelmd", "wnincd", "wnfetd" and "wnsumd"
    Then each answer is equal to the answer of CSPICE

  Scenario: Invalid endpoints
    Given a window
    When the user inserts an interval whose left endpoint is greater than the right one
    Then a "SpiceError" of kind "badendpoints" is thrown
    And the window is unchanged

  Scenario: Changes in place
    Given a window with random intervals, some of them singletons or touching
    When the user calls "wnexpd!", "wncond!", "wnfild!", "wnfltd!" or "wnextd!"
    Then the window is equal, endpoint for endpoint, to the window that CSPICE changes
    And intervals that overlap after the change are merged

  Scenario: Comparison of two windows
    Given two windows
    When the user calls "wnreld" with "=", "<>", "<=", "<", ">=" or ">"
    Then the answer is equal to the answer of CSPICE
    And blanks after the operator are ignored

  Scenario Outline: Invalid side or operator
    When the user calls "<call>" with the argument "<argument>"
    Then a "SpiceError" of kind "<kind>" is thrown

    Examples:
      | call    | argument | kind              |
      | wnextd! | X        | invalidendpntspec |
      | wnreld  | " ="     | invalidoperation  |
      | wnreld  | ""       | emptystring       |
