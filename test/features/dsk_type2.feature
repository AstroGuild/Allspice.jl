Feature: DSK type 2 shape models (triangular plates) with the results of CSPICE
  A type 2 DSK segment holds vertices, plates and a voxel index. NeoSpice reads such
  segments, finds where a ray meets the surface, and writes segments. CSPICE is the
  reference: the same plate and the same point, to rounding level.
  Executable form: test/test_dsk.jl

  Background:
    Given plate models that the test generates: a tessellated ellipsoid, rough random
      surfaces, cubes, and two halves of one surface
    And a text kernel with three Mars-fixed frames, the Mars constants and surface names
    And loaded kernels that belong to an immutable snapshot, with no global state

  Scenario Outline: The writer gives the file that CSPICE gives
    Given the plate model <model> with bounds in the <system> coordinate system
    When NeoSpice makes the spatial index with dskmi2 and writes the segment with dskw02
    And CSPICE does the same with the same input
    Then each element of the two spatial indexes is equal
    And each integer and each double of the two files is equal
    And the bytes of the files differ at most in the padding of the last records
    And CSPICE reads the NeoSpice file as its own

    Examples:
      | model            | system       |
      | ellipsoid        | latitudinal  |
      | ellipsoid        | planetodetic |
      | rough surface    | latitudinal  |
      | surface halves   | latitudinal  |
      | cube             | rectangular  |

  Scenario: Access calls give the items that CSPICE gives
    Given a DSK file and a segment of it
    When dskgd, dskb02, dskz02, dskv02, dskp02, dskn02, dski02 and dskd02 are called
    Then each result is equal to the result of CSPICE, for the whole item and for
      random ranges
    And dskobj and dsksrf give the sets of CSPICE

  Scenario Outline: A ray meets the same plate at the same point
    Given nine loaded segments of one body with different surfaces and overlapping coverage
    And a surface list, an epoch and the body-fixed frame <frame>
    When dskxsi is called for a ray of the kind <kind>
    Then NeoSpice finds an intercept if and only if CSPICE finds one
    And the segment and the plate are the same
    And the points differ by at most <tolerance> km

    Examples:
      | frame     | kind                           | tolerance |
      | IAU_MARS  | from outside toward the body   | 0         |
      | IAU_MARS  | from outside, random direction | 0         |
      | IAU_MARS  | from inside                    | 0         |
      | IAU_MARS  | toward a vertex or an edge     | 0         |
      | IAU_MARS  | vertex on the surface          | 0         |
      | MARS_TILT | mixed                          | 1e-11     |
      | MARS_PCK2 | mixed                          | 1e-11     |

  Scenario: Segments are selected by coverage in time
    Given two segments of one surface with the coverage 0 to 100 s and 50 to 200 s
    When a ray is traced at an epoch
    Then the segment that covers the epoch gives the intercept
    And the later segment gives it where both cover the epoch
    And no intercept is found outside of 0 to 200 s

  Scenario: Tolerances belong to the context
    Given two contexts with the same kernels
    When dskstl sets a large plate expansion in one context
    Then the intercepts of that context are those of CSPICE with the same setting
    And the intercepts of the other context are not changed
    And dskstl refuses the keywords 5 and 6 with the kind immutablevalue

  Scenario: Surface points and normals
    Given the methods ELLIPSOID and DSK/UNPRIORITIZED with a surface list by ID or name
    When latsrf is called with longitudes and latitudes
    And srfnrm is called with the points
    Then the points and the unit normals are those of CSPICE
    And a bad method string gives the error kind of CSPICE

  Scenario: The intercept does not allocate and is type stable
    Given integer IDs for the target and the frame and static vectors for the ray
    When dskxsi or dskxv! is called after a warm-up call
    Then no memory is allocated
    And the return type is inferred

  Scenario: Many threads trace rays while another thread loads and unloads a DSK
    Given a context with one DSK file and a second DSK file
    When several tasks trace rays without a lock
    And one task loads and unloads the second file in a loop
    Then each result is the result of the set without the file or of the set with it
    And no call throws

  Scenario Outline: Refuse bad input
    When <call> is made with <input>
    Then a SpiceError of the kind that CSPICE signals is thrown

    Examples:
      | call    | input                               |
      | dskxsi  | pri = true                          |
      | dskxsi  | a frame that is not centred on the target |
      | dskxsi  | a surface list without data         |
      | dskxsi  | a zero ray direction                |
      | dskmi2  | a fine voxel scale of zero          |
      | dskmi2  | a workspace that is too small       |
      | dskw02  | bounds out of order                 |
      | dskw02  | a vertex number outside of 1:NV     |
      | dski02  | an item that does not exist         |
      | dskobj  | a file that is not a DSK            |
