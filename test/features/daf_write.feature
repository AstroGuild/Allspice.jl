Feature: Write SPK and CK files that CSPICE reads as its own
  Allspice writes binary kernels (DAF files) through a handle that the caller owns.
  CSPICE is the reference: a file that Allspice writes must hold the same data, the
  same descriptors, the same segment names and the same comments as the file that
  CSPICE writes from the same input.
  Executable form: test/test_dafwrite.jl

  Background:
    Given a temporary directory
    And random, valid segment data from a seeded generator

  Scenario Outline: A writer gives the data area that CSPICE gives
    Given valid input for the writer <writer>
    When Allspice writes the segment to one new file
    And CSPICE writes the same segment to a second new file
    Then each double of the data area has the same bit pattern in the two files
    And the descriptors and the segment names are the same
    And the two files are equal byte for byte
    And CSPICE computes the same result from the two files

    Examples:
      | writer |
      | spkw01 |
      | spkw02 |
      | spkw03 |
      | spkw05 |
      | spkw08 |
      | spkw09 |
      | spkw12 |
      | spkw13 |
      | spkw15 |
      | spkw17 |
      | spkw18 |
      | spkw19 |
      | spkw20 |
      | spkw21 |
      | ckw01  |
      | ckw02  |
      | ckw03  |
      | ckw05  |
      | ckw06  |

  Scenario: An invalid input is rejected and the file does not change
    Given an open file and input that breaks one precondition of a writer
    When the writer is called
    Then it throws a SpiceError whose kind is the short message of CSPICE in lower case
    And CSPICE rejects the same input
    And the closed file is the file that holds the valid segments only

  Scenario: A file holds more segments than one summary record
    Given more than 25 segments
    When Allspice writes them to one file
    Then the file is equal byte for byte to the file that CSPICE writes
    And a search forward and a search backward find all segments

  Scenario: One engine appends to the file of the other engine
    Given a file with 23 segments that CSPICE wrote
    When Allspice opens the file with spkopa and adds segments
    Then the file holds the data that one engine writes alone

  Scenario: Comments go in both directions
    Given a file with a reserved comment area of some size, also zero
    When comment lines are added with dafac before and after the segments
    Then dafec of Allspice and dafec of CSPICE give the same lines from the file
    And the segments stay readable when the comment area grows and the data move
    And dafdc removes the comment area with the same result in both engines

  Scenario: A handle is used only as its access mode permits
    Given a handle that is open for read, or a closed handle
    When a write call is made on the handle
    Then it throws SpiceError(:dafinvalidaccess) or SpiceError(:dafnosuchhandle)

  Scenario: A file of the other byte order is read but not written
    Given a copy of an SPK file in the non-native byte order
    When Allspice reads the copy
    Then summaries, names, comments and data are those of the native file
    And dafopw throws SpiceError(:unsupportedbff)

  Scenario: A written file is loaded
    Given an SPK file with a type 2 and a type 3 segment that Allspice wrote
    When the file is loaded with furnsh
    Then spkgeo gives the states that CSPICE gives from the same file
    And ktotal, kdata and kinfo report the load list as CSPICE does

  Scenario: Bodies and coverage of an SPK file
    Given an SPK file with segments for some bodies
    When spkobj and spkcov read the file
    Then the body list and each coverage window are those of CSPICE

  Scenario: Tasks write different files at the same time
    Given one handle for each task
    When the tasks write, comment and read their own files concurrently
    Then each file is equal byte for byte to the file of a serial run
