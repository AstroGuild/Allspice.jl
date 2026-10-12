Feature: Geometry finder, concurrent use
  A search keeps no global state. The step, the tolerance and the callbacks are
  arguments, and the search reads one immutable snapshot of the loaded kernels.

  Scenario: The same searches on many tasks
    Given a context with the standard kernels loaded
    And the windows of a set of searches made serially
    When many tasks make the same searches at the same time
    Then every window is equal, bit for bit, to the window of the serial search

  Scenario: A search that uses tasks inside
    Given a confinement window of one or of several intervals
    When a search runs with the keyword "threads" set to true
    Then the window is equal, bit for bit, to the window of the search without the keyword

  Scenario: The tolerance belongs to the context
    Given two contexts with different tolerances set by "gfstol"
    When tasks search on the two contexts at the same time
    Then each search uses the tolerance of its own context

  Scenario: Searches while another task loads and unloads a kernel
    Given a context with the standard kernels loaded
    When one task repeatedly loads and unloads a kernel that the searches do not use
    And other tasks search on the same context
    Then every window is equal, bit for bit, to the window of the serial search
    And no task throws an error
