# CINCH Unified master plan

## Control rules

[KNOWN | HIGH] Status vocabulary is restricted to `NOT_STARTED`, `IN_PROGRESS`, `IMPLEMENTED`, `TESTED`, `VALIDATED`, `BLOCKED`, and `DEFERRED`.

[KNOWN | HIGH] A milestone is `VALIDATED` only when its implementation, tests, scientific acceptance criteria, performance criteria, and documentation gates are all satisfied.

## Milestones

| ID | Objective | Required sources | Expected outputs | Dependencies | Acceptance and scientific validation | Performance requirement | Status | Evidence |
|---|---|---|---|---|---|---|---|---|
| M00 | Freeze sources and install the audit system | CINCH `98c41ef`; CINCH-dev2 `ec1eaa5` | ten control documents; migration manifest | none | clean source snapshots; baseline tests recorded; no source destroyed | none | VALIDATED | `TEST_EVIDENCE.md`; this branch |
| M01 | Complete architecture and feature audit | both packages, CLIs, tests, frozen contracts | `FEATURE_INVENTORY.md`; architecture comparison | M00 | every major scientific module classified; differences explicit | bottleneck hypotheses only, not invented measurements | VALIDATED | Checkpoint 1 review in `ARCHITECTURE_DECISIONS.md` |
| M02 | Establish unified import/API skeleton without changing frozen_v1 | dev2 exact kernels; frozen_v1 package | `cinch/statistics`, compatibility exports, regression tests | M01 | weighted/unweighted definitions remain separate; exact kernel tests pass | one-million-pair measured microbenchmark | VALIDATED | 17 tests; zero cross-repository numerical difference; benchmark JSON |
| M03 | Implement mapper backend contract and correct frozen mapper defects | frozen_v1 mapping; legacy `Cinch_v8.py` | backend interface; corrected internal backend; legacy adapter | M02 | required controlled mapping cases classified across backends | real threads; reusable index; resumable per genome | NOT_STARTED | none |
| M04 | Unify presence/type profiles and compatibility adapters | frozen states; legacy profiles | explicit state schema and adapters | M03 | deterministic sample/locus/type ordering; missingness round-trip | bounded profile memory | NOT_STARTED | none |
| M05 | Implement bounded-memory PP/PT/TP/TT pair engine | frozen channel semantics; dev2 kernels | blockwise pair engine; restartable block output | M04 | exhaustive small-case equality for pair universe, MI/NMI, drivers | benchmark ladder through feasible levels | NOT_STARTED | none |
| M06 | Integrate population and advanced statistics | HC/Neff; weighting; SHC; BH; ARACNE; Diff-GWES | separate named methods and workflows | M05 | formula/output regressions with documented tolerances | permutation and graph stages profiled | NOT_STARTED | none |
| M07 | Expose staged `cinch` workflow and unified schemas | all prior modules | `map`, `profile`, `associate`, `filter`, `report`, `wgs` | M03-M06 | representative end-to-end run; cache invalidation; schema compatibility | bounded stages and resume evidence | NOT_STARTED | none |
| M08 | Production acceptance at 224 genomes × 25,000 loci | production inputs or approved equivalent | benchmark report and release decision | M07 | mandatory gates in assignment satisfied | actual wall time, CPU, RSS and throughput | BLOCKED | production dataset and HPC run not yet supplied/executed |

## Predefined completion calculation

[KNOWN | HIGH] Each milestone has five equal gates: implementation, engineering tests, scientific validation, performance evidence, and documentation. Project completion is `satisfied gates / 45`; blocked or deferred gates count as unsatisfied. No subjective weighting is permitted.

[COMPUTED | HIGH] After M02, M00-M02 satisfy all fifteen gates; M03-M08 satisfy zero gates. Overall completion is therefore `15/45 = 33.3%`.
