Feature: Read an EK file through a handle that the caller owns
  An EK file holds tables in segments. A segment has typed columns (character, double,
  integer, time); an entry is one value or an array of values, or null. Allspice reads
  an EK file through a DAS handle object, with no global state. CSPICE is the reference.
  Executable form: test/test_ek_file.jl

  Background:
    Given EK files that the CSPICE writers made from a seeded random model
    And the model has every column class, null entries and entries over several pages

  Scenario: Count and summarize the segments
    Given an EK file with several segments
    When Allspice opens the file with ekopr
    Then eknseg gives the segment count that CSPICE gives
    And ekssum gives the table name, the row count and the column attributes of CSPICE

  Scenario: Read every column entry
    Given a segment made with ekappr and the ekace* calls
    And a segment made with the fast writers
    When ekrcec, ekrced and ekrcei read every entry of every row
    Then each entry is equal to the entry that CSPICE reads and to the model
    And a null entry is reported as null

  Scenario: Reject wrong input with the error of CSPICE
    Given an open EK file
    When a call names a segment, a row or a column that does not exist
    Or the column does not have the data type of the call
    Then a SpiceError is thrown whose kind is the short message of CSPICE

  Scenario: A closed handle cannot be read
    Given a handle that ekcls has closed
    When eknseg is called
    Then a SpiceError of kind dasnosuchhandle is thrown
