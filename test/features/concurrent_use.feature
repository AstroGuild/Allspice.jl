Feature: Concurrent use from several threads
  Queries read an immutable snapshot of the loaded kernels and take no lock.
  Loading and unloading replace the snapshot as a whole, so a query that runs
  at the same time sees either the old set of kernels or the new one, never a
  mixture.

  Background:
    Given Julia runs with more than one thread

  Scenario: Parallel queries on one context
    Given a context with the standard kernels loaded
    When many tasks query states, rotations and derivatives at the same time
    Then every result is equal, bit for bit, to the result of the same query made serially

  Scenario: Queries while another task loads and unloads a kernel
    Given a context with "de440s.bsp" loaded
    And the reference answers with and without "de405.bsp" are computed first
    When one task repeatedly loads and unloads "de405.bsp"
    And other tasks query states on the same context
    Then every result is equal to one of the two reference answers
    And no task throws an error
    And after the writer stops, the context gives the answers without "de405.bsp"

  Scenario: Independent contexts
    Given each task creates its own context and loads its own ephemeris
    When the tasks query at the same time
    Then each task gets the answers of its own ephemeris

  Scenario: The oracle is not thread safe
    Given CSPICE keeps one process-wide state
    Then the tests call CSPICE from one task only
