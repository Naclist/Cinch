# Recovery and next actions

Last completed milestone: M01 architecture and feature audit.

Current active task: M02 unified exact statistical API.

Last verified Git commit: source baseline `98c41ef081be4face6b32cb61b557964537cb353`; integration checkpoint commit pending.

Files modified: `.gitignore`; `docs/development/*`.

Tests passed: CINCH 12/12; CINCH-dev2 7/7 under Python 3.12.14.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: uberBlast/configure availability; mapper semantics; missingness compatibility; blockwise engine; production inputs/resources.

Immediate next action: import dev2 weighted MI, directional EpiDis, and BH into differentiated `cinch.statistics` modules; add source-equivalence and frozen-v1 non-regression tests.
