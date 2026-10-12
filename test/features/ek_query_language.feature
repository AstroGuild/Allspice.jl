Feature: Parse and check a query of the EK query language
  A query selects columns from loaded tables: SELECT <columns> FROM <tables> [WHERE
  <constraints>] [ORDER BY <columns>]. Before any row is read, NeoSpice turns the
  string into a checked structure in six steps (scan, parse, normalise the WHERE clause,
  resolve names, resolve times, check types and sizes), with the acceptance rules, the
  limits and the error messages of CSPICE. A query that is not valid gives an error flag
  and a message, never an exception; an exception is for the conditions that CSPICE
  signals as SPICE errors.
  Executable form: test/test_ek_parse.jl

  Background:
    Given an EK file that CSPICE wrote with the tables ALPHA, BETA and GAMMA_LONG_NAME
    And columns of the types CHR, DP, INT and TIME, scalar and array
    And the table catalogue of CSPICE (ekntab, ektnam, ekccnt, ekcii) as the catalogue of NeoSpice
    And a leap-second kernel and a spacecraft clock kernel in both engines

  Scenario: A valid query has the encoded form of CSPICE
    Given the query "SELECT A.ICOL, B.ICOL FROM ALPHA A, BETA B WHERE A.ICOL = B.ICOL ORDER BY B.XCOL DESC"
    When NeoSpice encodes the query
    Then the error flag is false
    And the tables, aliases, constraints, order-by columns and SELECT columns equal those that zzekencd_ of CSPICE gives
    And the positions of the SELECT columns in the query string equal those of ekpsel_c

  Scenario: The WHERE clause becomes a disjunction of conjunctions
    Given the query "SELECT ICOL FROM ALPHA WHERE NOT (ICOL = 1 OR DCOL BETWEEN 2 AND 3) AND CCOL LIKE 'S*'"
    When NeoSpice encodes the query
    Then NOT is moved to the relations and replaces each operator by its complement
    And BETWEEN is expanded to two conjunctions of two relations
    And the conjunctions and the constraints in them have the order of CSPICE

  Scenario Outline: Operators and literals
    Given the query "SELECT ICOL FROM ALPHA WHERE <constraint>"
    When NeoSpice encodes the query
    Then the constraint has the operator code <code> and the value type <type>

    Examples:
      | constraint            | code | type |
      | ICOL = 1              | 1    | INT  |
      | ICOL GE 1.5           | 2    | DP   |
      | ICOL ^= 1             | 6    | INT  |
      | CCOL NOT LIKE 'it''s' | 8    | CHR  |
      | ICOL IS NULL          | 9    | none |
      | ICOL != NULL          | 10   | none |
      | TCOL < '2000 JAN 1'   | 5    | TIME |

  Scenario: A time value is converted to ephemeris time
    Given the query "SELECT ICOL FROM ALPHA WHERE TCOL > '2000 JAN 1 12:00:00 TDB'"
    When NeoSpice encodes the query with the time defaults of the caller
    Then the value of the constraint has the type TIME and the double 0.0
    And a spacecraft clock string "EKQ_CRAFT SCLK 1/1000:5" gives the time that CSPICE gives

  Scenario Outline: A query that is not valid gives the message of CSPICE
    Given the query "<query>"
    When NeoSpice encodes the query
    Then the error flag is true
    And the message is "<message>"

    Examples:
      | query                                    | message                                                                                  |
      | SELECT ICOL                              | Every query must contain a FROM clause. The FROM keyword was not found.                  |
      | SELECT * FROM ALPHA                      | Table or column name expected at location 8. Actual token was: *                         |
      | SELECT ICOL FROM ALPHA WHERE ICOL = @    | Invalid character found at location 37.  <character> = '@'                               |
      | SELECT ICOL FROM ALPHA WHERE ICOL =      | Tokens were missing from comparison relation.                                            |
      | SELECT ICOL FROM ALPHA, BETA             | Column name <ICOL> is ambiguous; a qualifying table name or alias is required.           |
      | SELECT ICOL FROM NOSUCH                  | Table <NOSUCH> is not currently loaded.                                                  |
      | SELECT ICOL FROM ALPHA WHERE CCOL = 1    | Data type mismatch: column <CCOL> has data type CHARACTER; value <1> has data type INTEGER. |
      | SELECT ICOL FROM ALPHA WHERE ICOL LIKE 'A' | LIKE and NOT LIKE operators may be used only with character columns.  Column <ICOL> has type INTEGER. |

  Scenario: Random queries agree with CSPICE
    Given 2500 random queries with nested AND, OR, NOT, BETWEEN and parentheses
    And 4000 random strings of tokens of the query language
    When both engines encode each query
    Then the error flags are equal and the error messages are equal
    And for a valid query the integer cell, the character buffer and the double buffer of the encoded query are equal element by element

  Scenario: The limits of CSPICE hold
    Given a query with 11 tables, or 51 SELECT columns, or 11 order-by columns, or 501 tokens, or 101 numbers
    When NeoSpice encodes the query
    Then the error flag is true and the message names the limit
    And a WHERE clause whose normal form has more constraints than the encoded query can hold throws SpiceError(:celltoosmall), as CSPICE signals SPICE(CELLTOOSMALL)

  Scenario: The SELECT clause is analysed as ekpsel does
    Given the query "SELECT B.TCOL, DCOL, CARR, B.ICOL FROM ALPHA, BETA B"
    When NeoSpice analyses the SELECT clause
    Then it returns 4 items of the class COL with the types TIME, DP, CHR and INT
    And the table names BETA, ALPHA, ALPHA and BETA, never the alias
    And the first and last character of each item, counted from 1

  Scenario: The parser has no shared state
    Given the hand-written queries
    When several threads encode them at the same time
    Then each result equals the result of one thread
