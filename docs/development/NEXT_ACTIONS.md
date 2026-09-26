# Recovery and next actions

Last completed milestone: M02 unified exact statistical API. M03 is TESTED but not scientifically validated against legacy uberBlast.

Current active task: complete M07 advanced filter/report glue and M05 streaming distance/scaling.

Last verified Git commit: `f85b73c` M02 exact statistical kernels.

Files modified after M06: staged profile/association workflows, CLI commands, controlled E2E test, and usage documentation.

Tests passed: unified 35/35; CINCH-dev2 7/7; staged profile→associate E2E and resume pass.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: uberBlast/configure equivalence; O(p²) distance materialization; pair scaling; advanced filter/report integration; production inputs/resources.

Immediate next action: add staged candidate filtering with SHC/BH/ARACNE manifests, without changing frozen `cinch filter` semantics.
