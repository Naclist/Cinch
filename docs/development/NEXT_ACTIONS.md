# Recovery and next actions

Last completed milestone: P2 controlled population/habitat conditional association. M03 remains TESTED but not scientifically validated against legacy uberBlast.

Current active task: complete the v0.1.0 Research Preview public-clone and GitHub Actions release gates without changing scientific methods.

The P2 implementation is committed on `integration/unified-framework`; verify local/remote SHA equality after each push.

Files modified: categorical conditional kernels, `conditional-associate` workflow/CLI, Cases A-J and audit-schema regressions, mathematical documentation, and development evidence.

Tests passed: unified 60/60 in 2.74 s plus 15/15 P2 controlled tests in 1.31 s; CINCH-dev2 historical baseline 7/7; Cases A-J, permutation-invariant support, stable audit schemas, and all-channel CLI E2E pass.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: real-study conditional exchangeability/covariates; continuous or omnibus background methods; uberBlast/configure equivalence. P0 25K and P1 full-SPN534 acceptance are user-deferred.

Immediate next action: push the release candidate, validate a fresh checkout of its immutable SHA, and verify the Python 3.10/3.12 plus package GitHub Actions jobs. Do not create a tag or GitHub Release before those checks pass.
