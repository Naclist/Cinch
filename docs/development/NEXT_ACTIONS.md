# Recovery and next actions

Last completed milestone: M02 unified exact statistical API. M03 is TESTED but not scientifically validated against legacy uberBlast.

Current active task: M05 bounded-memory PP/PT/TP/TT pair engine.

Last verified Git commit: `f85b73c` M02 exact statistical kernels.

Files modified after M03: `cinch/profiles/*`; `cinch/workflow/mapping.py`; profile tests and schema documentation.

Tests passed: unified 27/27; CINCH-dev2 7/7; 10,000-case numerical comparison has zero difference for M02.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: uberBlast/configure availability and mapper equivalence; blockwise engine; production inputs/resources.

Immediate next action: implement deterministic triangular block enumeration and prove exact small-case equality to frozen `score_all_pairs`.
