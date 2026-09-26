# Recovery and next actions

Last completed milestone: M02 unified exact statistical API. M03 is TESTED but not scientifically validated against legacy uberBlast.

Current active task: compiled multi-channel block scoring and broader SPN534 scaling after successful public-subset E2E.

Last verified Git commit before current SPN/mapping change: `bc933fc` staged SHC filtering and reporting.

Files modified: explicit no-hit policy, necessary-state prescreen, SPN534 benchmark assets/results, mapping/SPN audits, and regression tests.

Tests passed: unified 39/39; CINCH-dev2 7/7; real 400-locus old/new PP/PT/TP/TT outputs match exactly.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: uberBlast/configure equivalence; full-reference pair runtime; categorical population testing; historical SPN reference; private HPC inputs/resources.

Immediate next action: prototype a compiled block kernel that scores all legal channels together, and reject it unless exact output and meaningful real-data speedup both hold.
