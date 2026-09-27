# Recovery and next actions

Last completed milestone: v0.1.0 Research Preview packaging and public-install verification. M03 remains TESTED but not scientifically validated against legacy uberBlast.

Current active task: none. The Research Preview release commit is ready for an explicit tag/GitHub Release action; no tag or GitHub Release was created during preparation.

The P2 implementation is committed on `integration/unified-framework`; verify local/remote SHA equality after each push.

Files modified: categorical conditional kernels, `conditional-associate` workflow/CLI, Cases A-J and audit-schema regressions, mathematical documentation, and development evidence.

Tests passed: unified 62/62 on Python 3.12, release-core tests on Python 3.10, 15/15 P2 controlled tests, and CINCH-dev2 historical baseline 7/7. A public HTTPS clone passed both CLI smokes, distribution build/audit, Twine checks, and isolated wheel installation. GitHub Actions run `36284486133` passed its Python 3.10, Python 3.12, and package jobs.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: real-study conditional exchangeability/covariates; continuous or omnibus background methods; uberBlast/configure equivalence. P0 25K and P1 full-SPN534 acceptance are user-deferred.

Immediate next action: if publication is authorized, create the immutable tag and GitHub Release from the reported final commit using the title `Cinch v0.1.0 — Research Preview`. Production-scale acceptance remains a separate user-deferred gate.
