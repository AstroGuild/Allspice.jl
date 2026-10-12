Feature: Read and write DAS files through a handle that the caller owns
  A DAS file holds three address spaces (characters, doubles, integers), a comment
  area, and directory records that tell where the data records are. NeoSpice reads a
  DAS file of either byte order, segregated or not, and writes segregated files.
  CSPICE is the reference: each engine must read the files of the other.
  Executable form: test/test_das.jl

  Background:
    Given a temporary directory
    And no global state: each open file is one handle object

  Scenario: Read a file that CSPICE wrote
    Given a DSK file with three segments and five comment lines that CSPICE wrote
    When NeoSpice opens the file with dasopr
    Then dasrfr, daslla and dashfs give the values that CSPICE gives
    And dasec gives the comment lines that CSPICE gives
    And every integer and every double has the bit pattern that CSPICE reads

  Scenario: CSPICE reads a file that NeoSpice wrote
    Given a new DAS file made with dasonw
    When integers, doubles and characters are added in 120 appends of random size
    And comment lines are added with dasac
    And the file is closed with dascls
    Then CSPICE reads the same file record, the same summary and the same data
    And the special doubles NaN, -0.0, Inf and the smallest subnormal keep their bits

  Scenario: Change a file through dasopw
    Given an existing DAS file
    When values are replaced with dasudi and dasudd and a comment line is added
    Then CSPICE reads the new values and the new comment
    And the characters of the file are not changed
    When the comment area is deleted with dasdc
    Then the file has no comment records

  Scenario: Read a file that is not segregated
    Given a DAS file with 360 interleaved clusters and more than one directory record
    When NeoSpice opens the file
    Then the three address spaces hold the data in the order of the appends
    And dashfs gives the summary that CSPICE gives
    When the file is opened for write and closed
    Then the file is segregated and holds the same data

  Scenario Outline: Read a file of the other byte order
    Given the file <file> with each integer and double record byte-swapped
    When NeoSpice opens the swapped file
    Then the data are equal to the data of the native file
    And CSPICE reads the same values from the swapped file

    Examples:
      | file              |
      | segregated DSK    |
      | not segregated    |
      | generic DAS       |

  Scenario: Walk the DLA segment list
    Given a DSK file with three segments
    When the list is walked forward with dlabfs and dlafns
    And backward with dlabbs and dlafps
    Then each descriptor is equal to the descriptor that CSPICE returns
    And the backward list is the reverse of the forward list
    And a file without segments gives no descriptor

  Scenario Outline: Refuse bad input
    When <call> is made with <input>
    Then a SpiceError of kind <kind> is thrown
    And no file is created or changed

    Examples:
      | call    | input                          | kind             |
      | dasopr  | a file that does not exist     | filenotfound     |
      | dasopr  | a damaged FTP string           | ftpxfererror     |
      | dasopr  | a file cut short               | baddasdirectory  |
      | dasopr  | an unknown binary format       | unknownbff       |
      | dasrdi  | an address outside the file    | dasnosuchaddress |
      | dasadi  | a read handle                  | dasinvalidaccess |
      | dashfn  | a closed handle                | dasnosuchhandle  |
      | dasonw  | a file that exists             | fileopenfailed   |
      | dasonw  | a blank file type              | blankfiletype    |
      | dasac   | a nonprinting character        | illegalcharacter |
