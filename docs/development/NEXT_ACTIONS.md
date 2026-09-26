# Recovery and next actions

Last completed milestone: P2 controlled population/habitat conditional association. M03 remains TESTED but not scientifically validated against legacy uberBlast.

Current active task: P2 controlled implementation is complete; next scientific work is a preregistered real metadata-rich habitat application.

The P2 implementation is committed on `integration/unified-framework`; verify local/remote SHA equality after each push.

Files modified: categorical conditional kernels, `conditional-associate` workflow/CLI, Cases A-J and audit-schema regressions, mathematical documentation, and development evidence.

Tests passed: unified 60/60 in 2.74 s plus 15/15 P2 controlled tests in 1.31 s; CINCH-dev2 historical baseline 7/7; Cases A-J, permutation-invariant support, stable audit schemas, and all-channel CLI E2E pass.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: real-study conditional exchangeability/covariates; continuous or omnibus background methods; uberBlast/configure equivalence. P0 25K and P1 full-SPN534 acceptance are user-deferred.

Immediate next action: select a real dataset with overlapping habitats inside populations, freeze its metadata/candidate family, and execute `conditional-associate`. Do not resume P0/P1 unless the user reports a P2-blocking defect.
