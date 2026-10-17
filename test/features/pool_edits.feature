Feature: Kernel-pool edits and watchers
  A user changes variables of the kernel pool directly, without a kernel file. Each
  edit makes a new snapshot of the context, so queries on other threads see the old
  pool or the new one, never a mixture. The rule for edits and loaded kernels is the
  rule of CSPICE.

  Background:
    Given a context and CSPICE with the same kernels loaded
    And the same sequence of operations is applied to both

  Scenario: Set, replace and delete variables
    When the user calls "pdpool", "pipool" and "pcpool" with new names
    Then "gdpool", "gipool" and "gcpool" return the values
    When the user calls "pdpool" for "BODY399_RADII"
    Then "bodvrd" of "EARTH", "RADII" returns the new values
    When the user calls "dvpool" for a variable
    Then "expool" is false for it
    And after each step the whole pool is equal to the pool of CSPICE

  Scenario: An edit stays through a later furnsh
    Given the user set "PW_KEEP" with "pdpool"
    When the user loads a text kernel that does not assign "PW_KEEP"
    Then "PW_KEEP" has its value
    And a variable that the kernel assigns has the value of the kernel

  Scenario: An edit stays through the unload of a binary kernel
    Given the user set "PW_KEEP" with "pdpool"
    When the user unloads "de440s.bsp"
    Then "PW_KEEP" has its value

  Scenario: An edit is lost with the unload of a text kernel
    Given the user set "PW_KEEP" with "pdpool" and replaced "BODY399_RADII"
    When the user unloads a text kernel
    Then "PW_KEEP" is not in the pool
    And "BODY399_RADII" has the values of the PCK that is loaded

  Scenario: kclear removes every edit
    Given the user set "PW_KEEP" with "pdpool"
    When the user calls "kclear"
    Then the pool is empty

  Scenario: clpool keeps the load list
    When the user calls "clpool"
    Then the pool is empty
    And "ktotal" is unchanged
    When the user then unloads one text kernel
    Then the variables of the other text kernels are in the pool again

  Scenario: Load variables from text lines and from a file
    When the user calls "lmpool" with assignment lines
    Then the variables are in the pool, with "+=" appended to earlier values
    When the user calls "ldpool" with a text kernel
    Then the variables of the file are in the pool
    And "ktotal" is unchanged

  Scenario: Tables that are derived from the pool follow an edit
    When the user defines a fixed-offset frame with "lmpool"
    Then "namfrm" knows the frame and "pxform" agrees with CSPICE
    When the user sets "NAIF_BODY_NAME" with "pcpool" and "NAIF_BODY_CODE" with "pipool"
    Then "bodn2c" and "bodc2n" know the new body

  Scenario: Body names of boddef
    When the user calls "boddef" with a name and a code
    Then "bodn2c" and "bodc2n" agree with CSPICE
    And the name stays after "furnsh", "unload" and "kclear"

  Scenario: Read calls
    Then "dtpool" gives the count and the type of a variable, as CSPICE
    And "gnpool" gives the names that match a template with "*" and "%", as a sorted list
    And "szpool" gives the size parameters of CSPICE
    And "stpool" joins the strings that end with a continuation marker

  Scenario: A watcher sees the changes of its variables
    Given the agent "PW_AGENT_1" watches "PW_WATCHED" with "swpool"
    Then "cvpool" for the agent is true once
    And "cvpool" for the agent is then false
    When the user changes "PW_WATCHED" with "pdpool", "dvpool", "lmpool", "ldpool" or "furnsh"
    Then "cvpool" for the agent is true once
    When the user changes another variable
    Then "cvpool" for the agent is false
    And each answer is equal to the answer of CSPICE

  Scenario: Watches belong to a context
    Given two contexts in which the agent "A" watches "V"
    When the user changes "V" in the first context
    Then "cvpool" is true in the first context and false in the second

  Scenario: Edits from many threads
    Given four tasks that set their own variables 150 times each
    And four tasks that read variables and states at the same time
    Then no reader sees a variable with mixed values
    And every state is equal, bit for bit, to the state before the edits
    And every writer finds its own update with "cvpool" after each edit
    And at the end each variable has its last value

  Scenario: Body names for messages
    When the user calls "bodc2s" with a code that has a name and with one that has none
    Then the answer is the name, or the code as a string, as CSPICE

  Scenario Outline: Checks of a kernel variable
    Given a pool variable with a number of values and a type
    When the user calls "badkpv" with a comparison, a size, a divisor and a type
    Then a variable that passes gives false
    And a variable that fails throws a "SpiceError" of the kind of CSPICE

    Examples:
      | case                   | kind             |
      | variable not in pool   | variablenotfound |
      | comparison ">="        | unknowncompare   |
      | wrong size             | badvariablesize  |
      | divisor 0              | badvariablesize  |
      | wrong type             | badvariabletype  |
      | type letter "X"        | invalidtype      |
