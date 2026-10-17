Feature: Conversions between time strings, ephemeris time and uniform time scales
  The time functions give the results of CSPICE N0067: an epoch is equal bit for
  bit and a string is equal character for character. A string that CSPICE rejects
  throws a "SpiceError" with kind "badtime".

  Background:
    Given a context with "naif0012.tls" loaded
    And the defaults of "timdef" are calendar "GREGORIAN", system "UTC" and no zone

  Scenario Outline: A time string becomes ephemeris time
    When the user converts the string "<string>" with "str2et"
    Then the result is equal to the result of CSPICE bit for bit

    Examples:
      | string                              |
      | 1996-12-18T12:28:28                 |
      | 1995-08T18:28:12Z                   |
      | Mon Apr 22 09:32:38 PDT 1996        |
      | Dec 18, 11:59:59.9999 A.M. 1995     |
      | 18 B.C. Jun 3, 12:29:28.291         |
      | 1997-162::12:18:28.827              |
      | 2451515.2981 (JD)                   |
      | JD 2450480.02038194 TDB             |
      | 1997 Jan 31 12:29:20.123 UTC-07:30  |
      | '93 Jan 23 12:29:47.289             |

  Scenario Outline: A string without a meaning is rejected
    When the user converts the string "<string>" with "str2et"
    Then a "SpiceError" with kind "badtime" is thrown
    And the message names the part of the string that has no meaning
    And CSPICE rejects the same string

    Examples:
      | string                        |
      | 2000 Jann 1                   |
      | 2000-13-01                    |
      | 1999-02-29                    |
      | 2000 Jan 1 12:00:00 PDT TDB   |
      | 2000 Jan 1 13:00 PM           |
      | 2000-01-01T12:00:00+01:00     |
      | JD 2451545 Jan 1              |

  Scenario: A leap second is a valid second
    When the user converts "2016 Dec 31 23:59:60.5" with "str2et"
    Then the result is 0.5 seconds before the ephemeris time of "2017 Jan 1 00:00:00"
    And "et2utc" with format "C" and 3 decimals gives "2016 DEC 31 23:59:60.500" for that epoch

  Scenario: Second 60 at another time of the year is rejected
    When the user converts "2000 Jan 1 23:59:60" with "str2et"
    Then a "SpiceError" with kind "badtime" is thrown

  Scenario: A leap second in a time zone
    When the user converts "DEC 31 2016 3:59:60 PM PST" with "str2et"
    Then the result is the ephemeris time of "2016 Dec 31 23:59:60 UTC"

  Scenario Outline: Ephemeris time becomes a string with a picture
    Given the ephemeris time of "24 Mar 2018  16:23:00 UTC"
    When the user calls "timout" with the picture "<picture>"
    Then the string is "<string>"

    Examples:
      | picture                                 | string                           |
      | MON DD,YYYY  HR:MN:SC.#### (TDB) ::TDB  | MAR 24,2018  16:24:09.1856 (TDB) |
      | MON DD,YYYY  HR:MN:SC.#### (UTC)        | MAR 24,2018  16:23:00.0000 (UTC) |
      | Wkd Mon DD HR:MN:SC PDT YYYY ::UTC-7    | Sat Mar 24 09:23:00 PDT 2018     |

  Scenario: A picture made from an example string
    When the user calls "tpictr" with "Mon Apr 22 09:32:38 PDT 1996"
    Then the picture is "Wkd Mon DD HR:MN:SC PDT YYYY ::UTC-7"
    And "timout" with that picture writes an epoch in the form of the example

  Scenario Outline: UTC strings in the five formats of "et2utc"
    When the user calls "et2utc" with format "<format>" and 0 to 14 decimals
    Then the string is equal to the string of CSPICE for epochs over several thousand years
    And for epochs in the two seconds before and after every leap second

    Examples:
      | format |
      | C      |
      | D      |
      | J      |
      | ISOC   |
      | ISOD   |

  Scenario: The ISO formats have no year before 1 A.D. or after 9999
    When the user calls "et2utc" with format "ISOC" for an epoch in the year 100 B.C.
    Then a "SpiceError" with kind "badtime" is thrown

  Scenario: Uniform time scales
    When the user converts an epoch with "unitim" between two of "TAI", "TDT", "TT", "GPS", "TDB", "ET", "JDTDT", "JDTDB", "JED"
    Then the result is equal to the result of CSPICE bit for bit
    And the call allocates no memory
    And a "ForwardDiff.Dual" epoch gives the same value with its derivative

  Scenario: UTC is not a uniform time scale
    When the user calls "unitim" with the system "UTC"
    Then a "SpiceError" with kind "badtimetype" is thrown

  Scenario: The defaults of "timdef" belong to one context
    Given a second context with the same kernel
    When the user sets the system to "TDB" with "timdef" in the first context
    Then "str2et" of "2000 Jan 1 12:00:00" gives 0.0 in the first context
    And the second context still reads the string as UTC
    And "kclear" does not change the defaults of the first context

  Scenario: A value that "timdef" does not accept
    When the user sets the zone to "UTC+13" with "timdef"
    Then a "SpiceError" with kind "badtimedefault" is thrown
    And the three defaults are unchanged

  Scenario: The defaults change while other tasks convert
    Given tasks that call "str2et" and "timout" in a loop
    When another task changes the defaults many times
    Then every conversion uses the defaults before a change or the defaults after it
    And no reader waits for a lock

  Scenario Outline: No leap-second kernel
    Given a context with no kernel loaded
    When the user calls "<call>"
    Then a "SpiceError" with kind "<kind>" is thrown, as in CSPICE

    Examples:
      | call   | kind              |
      | str2et | noleapseconds     |
      | et2utc | missingtimeinfo   |
      | timout | missingtimeinfo   |
      | unitim | missingtimeinfo   |
      | deltet | kernelvarnotfound |

  Scenario: Range check of the time fields
    Given a context where "tparch" is set to "YES"
    When the user converts "1996-02-30", "1996-13-01" or "1996-367//"
    Then "utc2et" and "tparse" reject the string, as CSPICE with the check on
    And an "@" date in a text kernel that the context loads is checked in the same way
    But another context still accepts the string
    And any value other than "YES" turns the check off
