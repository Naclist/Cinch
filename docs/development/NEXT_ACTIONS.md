# Recovery and next actions

Last completed milestone: M02 unified exact statistical API. M03 is TESTED but not scientifically validated against legacy uberBlast.

Current active task: association core scaling and higher-locus feasibility after nested 12/24/48 SPN534 sample scaling.

Last verified Git commit before the current compiled-scorer change: `5b3aa44` public SPN validation and necessary-state prescreen.

Files modified: optional Numba four-channel scorer, CLI engine selection, performance dependency, real A/B evidence, and regression tests.

Tests passed: unified 45/45; CINCH-dev2 7/7; real 12/24/48-sample distance/non-floating outputs match exactly and maximum float difference is 2.22e-15.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: uberBlast/configure equivalence; full-reference pair runtime; categorical population testing; historical SPN reference; private HPC inputs/resources.

Immediate next action: design block-level association parallelism without duplicating profiles or changing output order, then test beyond 400 loci. Do not claim full 534 or 25K pair validation from the 48-genome result.
