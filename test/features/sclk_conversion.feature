Feature: Spacecraft clock conversion (SCLK type 1)
  A spacecraft clock kernel relates the count of an on-board clock to ephemeris
  time. The count is a string of fields with a partition number in front, or a
  number of ticks, where one tick is the smallest increment of the clock.

  Background:
    Given a context with a leapseconds kernel and a type 1 SCLK kernel loaded
    And the clock has more than one partition and more than one coefficient record

  Scenario: Ticks to ephemeris time and back
    Given an encoded clock value inside the range of the clock
    When the value is converted to ephemeris time with sct2e
    And the result is converted back with sce2c
    Then the ephemeris time is equal to the CSPICE answer
    And the round trip returns the first value to rounding accuracy

  Scenario: Discrete ticks
    Given an ephemeris time inside the range of the clock
    When the time is converted with sce2t
    Then the result is a whole number of ticks
    And the result is equal to the CSPICE answer

  Scenario: Clock strings
    Given a clock string with a partition number and fields
    When the string is encoded with scencd and decoded with scdecd
    Then the two results are equal to the CSPICE answers
    And scs2e and sce2s agree with the composition of the string and the time conversions

  Scenario: A clock with a terrestrial time parallel system
    Given an SCLK kernel that gives SCLK01_TIME_SYSTEM equal to 2
    When clock values are converted to ephemeris time
    Then the conversion goes through terrestrial time and is equal to the CSPICE answer

  Scenario: Derivative of the conversion
    Given an epoch of a Taylor or dual number type
    When the epoch is converted with sce2c
    Then the result carries the derivative of the piecewise linear clock function

  Scenario: Values that the clock cannot represent
    Given a clock value before the first partition or after the last one
    When the value is converted
    Then a SpiceError of the kind valueoutofrange is thrown
    And a string with a partition number that does not exist gives the kind badpartnumber

  Scenario: No clock kernel
    Given a context without an SCLK kernel for the spacecraft
    When a clock conversion is requested
    Then a SpiceError of the kind kernelvarnotfound is thrown

  Scenario: A clock of another type
    Given a clock whose "SCLK_DATA_TYPE" is not 1
    When a clock conversion is requested
    Then a SpiceError of the kind unsupported is thrown
