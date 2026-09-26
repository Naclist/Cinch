# Recovery and next actions

Last completed milestone: M02 unified exact statistical API. M03 is TESTED but not scientifically validated against legacy uberBlast.

Current active task: finish M05 scaling prerequisites, then M06 advanced statistics.

Last verified Git commit: `f85b73c` M02 exact statistical kernels.

Files modified after M04: `cinch/association/*`; block tests, benchmark, and documentation.

Tests passed: unified 30/30; CINCH-dev2 7/7; exact small-case frozen/block outputs match.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: uberBlast/configure equivalence; O(p²) distance materialization; pair scaling; advanced modules; production inputs/resources.

Immediate next action: implement blockwise distance production or provider interface, then measure the n/p/block-size ladder.
