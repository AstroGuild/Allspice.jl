Feature: Functions that never allocate
  The numeric kernels of the package run inside loops of an optimizer or a propagator,
  where one allocation for each call is a large cost.

  Scenario: A static analysis finds no allocation in a numeric kernel
    Given a function of the list of numeric kernels and its argument types
    When AllocCheck examines every path of the compiled function that does not throw
    Then it reports no allocation and no dynamic dispatch

  Scenario: A function that allocates is reported
    Given a function that returns a new array
    When AllocCheck examines it
    Then it reports at least one allocation
