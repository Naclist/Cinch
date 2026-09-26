# CINCH Unified Development Status

## Project objective

[KNOWN | HIGH] Produce one `cinch` framework that preserves frozen_v1 PP/PT/TP/TT and recurrence logic while integrating dev2 weighted information, EpiDis, SHC permutation, ARACNE, Diff-GWES, historical reproducibility, and scalable engineering.

## Source repositories and commits

| Source | URL | Commit | Role |
|---|---|---|---|
| CINCH | `https://github.com/Naclist/Cinch` | `98c41ef081be4face6b32cb61b557964537cb353` | destination; frozen_v1 WGS/filter |
| CINCH-dev2 | `https://github.com/Naclist/Cinch-dev2` | `ec1eaa5cf3c0c46a0d4d2babd52a768d23f9f57a` | exact kernels, HPC/history source |

## Current architecture

[COMPUTED | HIGH] Production behavior remains `cinch/frozen_v1`; the integration branch currently adds governance only. No frozen implementation has been overwritten.

## Overall milestone status

| ID | Milestone | Status | Evidence | Remaining work |
|---|---|---|---|---|
| M00 | Source freeze and audit system | VALIDATED | source SHAs; baseline tests | none |
| M01 | Full architecture audit | VALIDATED | inventory, decisions, Checkpoint 1 review | keep inventory current |
| M02 | Unified exact statistical API | IN_PROGRESS | plan and source crosswalk | implement and test |
| M03 | Mapper backends and resumability | NOT_STARTED | mapping audit | all code and controlled cases |
| M04 | Unified profile/state model | NOT_STARTED | state discrepancy recorded | adapter and tests |
| M05 | Blockwise four-channel engine | NOT_STARTED | complexity audit | implementation and benchmark |
| M06 | Advanced statistics | NOT_STARTED | distinct definitions inventoried | integration/regression |
| M07 | Staged CLI/end-to-end workflow | NOT_STARTED | target commands decided | implementation/E2E |
| M08 | Production acceptance | BLOCKED | target workload recorded | data and HPC execution |

[COMPUTED | HIGH] Predefined completion is `10/45 = 22.2%`; see `MASTER_PLAN.md`.

## Scientific feature preservation

| Feature | Original implementation | Unified implementation | Validation status |
|---|---|---|---|
| PP/PT/TP/TT | `cinch/frozen_v1/statistics.py` | unchanged frozen_v1 | baseline tested only |
| HC recurrence/Neff | `filtering.py` | unchanged frozen_v1 | baseline tested only |
| Weighted MI/EpiDis/BH | dev2 `core/statistics.py` | planned under `cinch/statistics` | dev2 unit tests passed; not migrated |
| SHC permutation/ARACNE/Diff-GWES | dev2 frozen snapshots | planned distinct modules | historical evidence only; not integrated |
| HPC mapping/allele nomenclature | dev2 `Cinch_v8.py` | planned backend adapter | external dependencies unresolved |

## Confirmed performance improvements

[COMPUTED | HIGH] None have been implemented or measured on this branch.

## Unverified performance assumptions

[INFERRED | HIGH] Genome-parallel indexed alignment, integer encoding, blockwise pair enumeration, compiled contingency kernels, and sparse distance storage are candidates; none is reported as a measured improvement yet.

## Current blockers

- [KNOWN | HIGH] `uberBlast` and `configure` redistribution/installability are unresolved.
- [KNOWN | HIGH] The 224-genome/25,000-reference production dataset is not present in this repository.
- [KNOWN | HIGH] Full 312,487,500-pair execution has not been resource-estimated from a unified block engine because that engine does not yet exist.

## Known scientific discrepancies

- [KNOWN | HIGH] frozen_v1 MI is unweighted and four-channel; dev2 weighted MI is a separate binary formulation.
- [KNOWN | HIGH] HC recurrence/Neff and SHC-conditioned permutation answer different population-structure questions.
- [KNOWN | HIGH] frozen_v1 missingness distinguishes non-callable type states, whereas historical dev2 profiles often collapse non-positive values.
- [KNOWN | HIGH] EpiDis is not NMI; its squared identity is restricted to the frozen weighted empirical construction.

## Tests executed

- [COMPUTED | HIGH] CINCH baseline: 12 passed on Python 3.12.14.
- [COMPUTED | HIGH] CINCH-dev2 baseline: 7 passed on Python 3.12.14.

## Tests not executed

[KNOWN | HIGH] Controlled cross-mapper cases, blockwise statistical equivalence, end-to-end unified execution, and production-scale benchmarks remain unexecuted.

## Next three priority actions

1. [KNOWN | HIGH] Import the exact dev2 statistical primitives under differentiated unified names and add cross-repository numerical tests.
2. [KNOWN | HIGH] Define the mapper backend/cache contract and fix `by_contig` reconstruction without changing coordinates.
3. [KNOWN | HIGH] Prototype bounded block enumeration for PP while preserving the full eligible hypothesis universe.

## Current recommended usage

[KNOWN | HIGH] Use released frozen_v1 for its documented WGS/filter workflow and dev2 only for its explicitly frozen research modules. Do not use the integration branch as a unified production release.

## Release readiness

[COMPUTED | HIGH] Not ready for unified research use; mandatory M03-M08 gates remain open.
