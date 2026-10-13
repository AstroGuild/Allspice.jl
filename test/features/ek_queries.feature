Feature: Query the loaded EK tables
  furnsh loads an EK file into a context; its segments join the tables of the context.
  ekfind runs a query of the EK query language on one snapshot of the loaded tables and
  returns a result object; eknelt, ekgc, ekgd and ekgi read that object. There is no
  global query state. CSPICE is the reference.
  Executable form: test/test_ek_query.jl, test/test_ek_threads.jl

  Background:
    Given three EK files that the CSPICE writers made from a seeded random model
    And the table EVENTS has segments in two files and GROUPS has segments of the fast writers
    And the files and a leapseconds kernel are loaded in CSPICE and in a Allspice context

  Scenario: An EK file is a loaded kernel
    When ktotal, kdata and kinfo are called with the kind "EK"
    Then they give the count, the file names and the type that CSPICE gives

  Scenario: Read the table catalogue
    When ekntab, ektnam, ekccnt and ekcii are called
    Then the tables, the columns and the column attributes are those of CSPICE
    And a table or an index that does not exist gives the error kind of CSPICE

  Scenario: Analyse a SELECT clause
    When ekpsel analyses a query
    Then the number of SELECT items, their positions, types, tables and columns are those of ekpsel_c

  Scenario Outline: Run a query
    Given a valid query with <clause>
    When ekfind runs it
    Then the row count is the one of CSPICE
    And every value that eknelt, ekgc, ekgd and ekgi return is the value of CSPICE, row by row
    Examples:
      | clause                                             |
      | each comparison operator on each data type         |
      | LIKE and NOT LIKE with both wild cards             |
      | IS NULL and IS NOT NULL                            |
      | BETWEEN and NOT BETWEEN                            |
      | AND, OR, NOT and parentheses                       |
      | a time value written as a string                   |
      | ORDER BY with several keys, ascending or descending |
      | a join of two tables on equality                   |
      | a join of two tables on an inequality              |
      | a join of three tables                             |
      | a table joined with itself through aliases         |

  Scenario: A query that is not valid
    Given a query with a syntax error, an unknown name or types that do not fit
    When ekfind runs it
    Then the error flag is set and the message is the message of CSPICE
    And the result has no rows
    And no exception is thrown

  Scenario: Generated queries agree with the model
    Given a query that a seeded generator made on one, two or three tables
    When ekfind runs it
    Then the rows are the rows that a search over the model gives with the rules of the query language

  Scenario: Results are correct where CSPICE has a defect
    Given a query that compares an indexed integer column with a number that has a fraction
    Or a query whose OR alternatives give the same row
    When ekfind runs it
    Then the rows are those of the model
    And the compatibility mode of the search gives the rows of CSPICE

  Scenario: Wrong use of the fetch calls
    Given a query result
    When a fetch call names a row, a SELECT column or an element number that is not valid
    Or the fetch call does not fit the data type of the column
    Then a SpiceError is thrown whose kind is the short message of CSPICE

  Scenario: A result outlives its files
    Given a query result
    When the EK files of the result are unloaded
    Then the result gives the same values as before

  Scenario: Unload and reload change the tables
    When an EK file is unloaded
    Then its tables and segments are not seen by later queries, as in CSPICE
    When the file is loaded again
    Then its segments are searched last, as in CSPICE

  Scenario: A segment that does not fit a loaded table
    Given a loaded table
    When furnsh loads an EK file whose segment of that table has other columns or attributes
    Then a SpiceError is thrown whose kind is the short message of CSPICE
    And the context is not changed

  Scenario: Many tasks run queries on one snapshot
    Given a context with loaded EK files
    When eight threads run queries at the same time
    Then each result is equal to the result of the serial run

  Scenario: Many tasks run queries while a file is loaded and unloaded
    Given eight tasks that run queries on one context
    And one task that loads and unloads an EK file with a segment of the table GROUPS
    When the tasks run together
    Then each result is the result with the segment or the result without it
    And both kinds of results are seen
