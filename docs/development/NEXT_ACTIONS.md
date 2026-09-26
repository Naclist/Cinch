# Recovery and next actions

Last completed milestone: M02 unified exact statistical API. M03 is TESTED but not scientifically validated against legacy uberBlast.

Current active task: finish M05 scaling prerequisites and connect M06 kernels through M07 staged workflow.

Last verified Git commit: `f85b73c` M02 exact statistical kernels.

Files modified after M05: `cinch/population/*`; `cinch/network/*`; high-order kernels; preserved dev2 sources; tests and benchmark.

Tests passed: unified 34/34; CINCH-dev2 7/7; advanced kernels have zero max difference over 100 cross-source cases.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: uberBlast/configure equivalence; O(p²) distance materialization; pair scaling; advanced workflow integration; production inputs/resources.

Immediate next action: add staged `profile`, `associate`, and `filter` CLI paths around existing validated kernels while keeping each cache boundary explicit.
