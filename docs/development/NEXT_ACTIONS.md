# Recovery and next actions

Last completed milestone: M02 unified exact statistical API. M03 is TESTED but not scientifically validated against legacy uberBlast.

Current active task: close M03 documentation checkpoint, then implement M04 explicit profile/state adapters.

Last verified Git commit: `f85b73c` M02 exact statistical kernels.

Files modified: `cinch/mapping/*`; `cinch/workflow/mapping.py`; `cinch/legacy/dev2/Cinch_v8.py`; CLI, tests, benchmarks, and development records.

Tests passed: unified 24/24; CINCH-dev2 7/7; 10,000-case numerical comparison has zero difference for M02.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: uberBlast/configure availability and mapper equivalence; missingness compatibility; blockwise engine; production inputs/resources.

Immediate next action: commit M03, implement M04 explicit `present/absent/unresolved` profile validation and conversion policies, then test round trips.
