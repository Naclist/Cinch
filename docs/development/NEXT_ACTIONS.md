# Recovery and next actions

Last completed milestone: M02 unified exact statistical API.

Current active task: M03 mapper backend and correctness contract.

Last verified Git commit: `f2eb881` audit checkpoint; M02 commit pending.

Files modified: `cinch/statistics/*`; `tests/test_unified_statistics.py`; benchmark and development records.

Tests passed: unified 17/17; CINCH-dev2 7/7; 10,000-case numerical comparison has zero difference.

Tests failed: none after using supported Python; initial install under Python 3.9.6 was rejected by package metadata.

Unresolved issues: uberBlast/configure availability; mapper semantics; missingness compatibility; blockwise engine; production inputs/resources.

Immediate next action: define mapping hit/backend/cache contracts, move coordinate grouping outside the locus loop, and add exact/SNP/reverse/no-hit/indel/truncation controlled cases.
