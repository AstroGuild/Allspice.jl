Feature: Write an EK file through a handle that the caller owns
  Allspice creates and changes EK files with the calls of CSPICE: segments, records,
  column entries, and the fast writers that store a segment column by column. A file
  that Allspice writes holds the same words as the file that CSPICE writes with the
  same calls, so each engine can read and change the files of the other. Allspice keeps
  no global state: the state of a fast load is a value that the caller holds.
  Executable form: test/test_ek_write.jl

  Background:
    Given a seeded random model with every column class and data type
    And null entries, indexed columns, and entries that continue over several data pages
    And CSPICE as the reference for files, entries, query results and errors

  Scenario: Write a new file
    Given the model of a file with several segments
    When CSPICE and Allspice write it with ekopn, ekbseg, ekappr, ekacec, ekaced and ekacei
    Then the character, double and integer address spaces of the two files are equal word for word
    And CSPICE reads every entry of the Allspice file as the model states it

  Scenario: Write a segment with the fast writers
    Given a segment whose columns are stored with ekifld, ekaclc, ekacld, ekacli and ekffld
    When the columns are scalar or array columns (classes 1 to 6)
    Or all columns have a fixed count (classes 7 to 9)
    Then the file of Allspice is equal to the file of CSPICE word for word
    And the columns can be stored in any order
    And two fast loads on two handles can run at the same time

  Scenario: Trees of depth 3
    Given a segment with more rows than a tree of depth 2 can hold
    And entries that fill more data pages than a tree of depth 2 can hold
    When both engines write the segment with the fast writers
    Then the record tree, the index trees and the data page trees have depth 3
    And the files are equal word for word

  Scenario: Edit trees of depth 3
    Given a file of Allspice with a segment that grew row by row to trees of depth 3
    When both engines delete and insert many rows
    Then the two results are equal word for word
    And both engines read every entry as the model states it

  Scenario: Trees of depth 3 that grow row by row (exhaustive run, ALLSPICE_TEST_FULL=1)
    Given the same segments, which CSPICE writes row by row in about one minute
    When both engines write the segments with ekappr, ekacec, ekaced and ekacei
    Then the files are equal word for word

  Scenario: Change an existing file
    Given a file that CSPICE wrote and the same file that Allspice wrote
    When each engine opens each file with ekopw and applies the same random script
    And the script appends, inserts and deletes rows, and replaces entries with entries of another size or with null
    And the script adds new segments with both kinds of writers
    Then the four results are equal word for word
    And the entries read back as the model after the script

  Scenario: Delete all rows and add rows again
    Given a segment with rows
    When all rows are deleted in random order and new rows are added
    Then the pages come back from the free lists in the order of CSPICE
    And the files of both engines are equal word for word

  Scenario: CSPICE queries a file of Allspice
    Given a file of Allspice with indexed columns in segments of both types
    When CSPICE loads the file and runs ekfind with equality, range and ORDER BY constraints
    Then the rows are those that the model gives, in the order of the index

  Scenario: Reject wrong input with the error of CSPICE
    Given a file open for write in both engines
    When a call has a bad table or column name, a bad column declaration, a bad segment or row number
    Or the column does not have the data type of the call, or does not allow a null entry
    Or the entry has the wrong number of elements, or the scalar entry has a value already
    Or the handle is a read handle or is closed
    Then Allspice throws a SpiceError whose kind is the short message of CSPICE
    And the handle can still be used
    And the two files are still equal word for word

  Scenario: Change a fast-loaded segment
    Given a segment that the fast writers made, with equal strings in an indexed column and array entries over several pages
    When Allspice deletes rows, inserts rows and replaces entries
    Then CSPICE reads every entry as the model states it and its queries use the indexes
    But CSPICE itself fails on these changes (ITEMNOTFOUND, or rows that it overwrites)

  Scenario: Records of a fixed-count segment cannot be changed
    Given a segment whose columns have a fixed count
    When ekappr, ekinsr or ekdelr is called for it
    Then Allspice throws SpiceError(:invalidtype) and the file is not changed
