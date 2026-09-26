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
| M03 | Implement mapper backend contract and correct frozen mapper defects | frozen_v1 mapping; legacy `Cinch_v8.py` | backend interface; corrected internal backend; preserved legacy source | M02 | controlled cases classified for frozen/indexed backends; uberBlast comparison blocked | real threads; reusable index; resumable per genome | TESTED | 24-test suite; controlled mapper benchmark; semantics below |
| M04 | Unify presence/type profiles and compatibility adapters | frozen states; legacy profiles | explicit state schema and adapters | M03 | deterministic sample/locus/type ordering; missingness round-trip | five bytes/cell raw state arrays | VALIDATED | `docs/PROFILE_SCHEMA.md`; 27-test suite |
| M05 | Implement bounded-memory PP/PT/TP/TT pair engine | frozen channel semantics; dev2 kernels | blockwise pair engine; restartable block output; sparse physical distances | M04 | exhaustive equality for pair metrics/drivers and reliable distances | preliminary 224×80 benchmark complete; scaling ladder incomplete | TESTED | exact A/B tests; 36-test suite; `ASSOCIATION_BLOCKS.md` |
| M06 | Integrate population and advanced statistics | HC/Neff; weighting; SHC; BH; ARACNE; Diff-GWES | separate named methods and workflows | M05 | zero-difference kernel migration and controlled PP workflow | kernels profiled; workflow-scale profile absent | TESTED | 37-test suite; cross-source comparison; staged SHC/BH/ARACNE test |
| M07 | Expose staged `cinch` workflow and unified schemas | all prior modules | `map`, `profile`, `associate`, `filter`, `report`, `wgs` | M03-M06 | controlled full staged workflow passes; public SPN534 E2E pending | mapping/association resume evidence; production scaling absent | IN_PROGRESS | CLI stages; 37-test suite; `STAGED_WORKFLOW.md` |
| M08 | Production acceptance at 224 genomes × 25,000 loci | production inputs or approved equivalent | benchmark report and release decision | M07 | mandatory gates in assignment satisfied | actual wall time, CPU, RSS and throughput | BLOCKED | production dataset and HPC run not yet supplied/executed |

## Predefined completion calculation

[KNOWN | HIGH] Each milestone has five equal gates: implementation, engineering tests, scientific validation, performance evidence, and documentation. Project completion is `satisfied gates / 45`; blocked or deferred gates count as unsatisfied. No subjective weighting is permitted.

[COMPUTED | HIGH] M00-M02 and M04 satisfy twenty gates. M03 and M05 each satisfy four gates. M06 satisfies implementation, engineering-test, scientific-equivalence, and documentation gates but lacks workflow-scale performance evidence. M07 satisfies implementation, engineering-test, and documentation gates but lacks public real-data scientific E2E validation and production performance evidence. M08 satisfies zero gates. Overall completion is therefore `35/45 = 77.8%`.

## M03 mapping-semantics checkpoint

| Case | frozen_v1 | unified indexed mapper | Status |
|---|---|---|---|
| exact forward hit | callable | callable | tested |
| exact reverse-complement hit | callable with reverse coordinates | callable with reverse coordinates | tested |
| SNP | callable above identity threshold | callable above identity threshold | tested |
| short indel | generally missed by fixed-length Hamming window | supported by gapped minimap2 alignment | intentionally different; tested |
| truncated locus | no hit becomes absence | detection-only/partial state | intentionally richer; tested |
| no hit | `presence=0` confident absence | `presence=-1` unresolved no hit | intentionally richer; tested |
| equal multicopy hits | non-callable | `AMBIGUOUS_MULTICOPY` | tested |
| competing reference loci | no explicit cross-locus adjudication | ambiguous/rejected competition states | tested |
| ambiguous bases | seed/window-dependent behavior | identity denominator includes query/target span | tested |

[KNOWN | HIGH] The unified mapper is not claimed to be numerically equivalent to historical uberBlast nomenclature. `configure` and `uberBlast` are unavailable, and their executable/license contract is unresolved.
