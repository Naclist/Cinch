# CINCH Unified Development Status

## Project objective

[KNOWN | HIGH] Produce one `cinch` framework that preserves frozen_v1 PP/PT/TP/TT and recurrence logic while integrating dev2 weighted information, EpiDis, SHC permutation, ARACNE, Diff-GWES, historical reproducibility, and scalable engineering.

## Source repositories and commits

| Source | URL | Commit | Role |
|---|---|---|---|
| CINCH | `https://github.com/Naclist/Cinch` | `98c41ef081be4face6b32cb61b557964537cb353` | destination; frozen_v1 WGS/filter |
| CINCH-dev2 | `https://github.com/Naclist/Cinch-dev2` | `ec1eaa5cf3c0c46a0d4d2babd52a768d23f9f57a` | exact kernels, HPC/history source |

## Current architecture

[COMPUTED | HIGH] Production WGS behavior remains `cinch/frozen_v1`; differentiated weighted MI, EpiDis, and BH APIs exist under `cinch/statistics`, and an opt-in indexed/resumable mapper exists under `cinch/mapping`. No frozen implementation has been overwritten.

## Overall milestone status

| ID | Milestone | Status | Evidence | Remaining work |
|---|---|---|---|---|
| M00 | Source freeze and audit system | VALIDATED | source SHAs; baseline tests | none |
| M01 | Full architecture audit | VALIDATED | inventory, decisions, Checkpoint 1 review | keep inventory current |
| M02 | Unified exact statistical API | VALIDATED | 17 tests; zero dev2 diff; measured microbenchmark | extend later with categorical kernels |
| M03 | Mapper backends and resumability | TESTED | 24-test suite; mapper benchmark; preserved `Cinch_v8.py` | legacy uberBlast equivalence and representative biological validation |
| M04 | Unified profile/state model | VALIDATED | explicit schema, legacy policy adapter, frozen/V2 round trips | none |
| M05 | Blockwise four-channel engine | TESTED | exact frozen A/B, atomic resume, controlled memory benchmark | streaming distances and scaling ladder |
| M06 | Advanced statistics | TESTED | SHC/ARACNE/Diff-GWES kernels, zero-difference migration, timing | staged workflow and workflow-scale profile |
| M07 | Staged CLI/end-to-end workflow | IN_PROGRESS | map/profile/associate controlled E2E and resume | advanced filter/report and representative E2E |
| M08 | Production acceptance | BLOCKED | target workload recorded | data and HPC execution |

[COMPUTED | HIGH] Predefined completion is `33/45 = 73.3%`; see `MASTER_PLAN.md`.

## Scientific feature preservation

| Feature | Original implementation | Unified implementation | Validation status |
|---|---|---|---|
| PP/PT/TP/TT | `cinch/frozen_v1/statistics.py` | unchanged frozen_v1 | baseline tested only |
| HC recurrence/Neff | `filtering.py` | unchanged frozen_v1 | baseline tested only |
| Weighted MI/EpiDis/BH | dev2 `core/statistics.py` | `cinch/statistics` | zero max absolute difference over 10,000 random valid inputs; unified tests pass |
| SHC permutation/ARACNE/Diff-GWES | dev2 frozen snapshots | `cinch.population`, `cinch.network`, `cinch.association.high_order` | kernels have zero-difference migration evidence; workflow absent |
| Indexed mapping/allele nomenclature | frozen mapper plus dev2 `Cinch_v8.py` | mappy/minimap2 index, deterministic hash-sorted alleles, atomic per-genome cache | controlled cases pass; historical uberBlast equivalence unresolved |

## Confirmed performance improvements

[COMPUTED | HIGH] Unified weighted MI processed 1,000,000 valid mass triplets at a best observed 12,027,537 pairs/s over five runs on the local Python 3.12.14/NumPy 2.5.3 environment. This is a kernel microbenchmark, not end-to-end throughput.

[COMPUTED | HIGH] At 224 samples × 80 loci, lazy blocks reduced traced peak Python allocation from 20.27 MB to 2.99 MB with exact row equality. Runtime changed from 3.212 s to 3.171 s; no material speedup is claimed.

## Unverified performance assumptions

[COMPUTED | HIGH] On 12 synthetic genomes × 200 exact loci, the indexed mapper completed in 0.534 s with one worker and 0.468 s with four workers, a measured 1.14× speedup. Two workers were slower than one. This workload is dominated by process/index overhead and is not biological validation.

[INFERRED | HIGH] Integer encoding, blockwise pair enumeration, compiled contingency kernels, and sparse distance storage remain unmeasured candidates.

## Current blockers

- [KNOWN | HIGH] `uberBlast` and `configure` redistribution/installability are unresolved.
- [KNOWN | HIGH] The 224-genome/25,000-reference production dataset is not present in this repository.
- [KNOWN | HIGH] Full 312,487,500-pair execution has not been resource-estimated; the block engine exists, but distance materialization and per-pair Python scoring remain scaling blockers.

## Known scientific discrepancies

- [KNOWN | HIGH] frozen_v1 MI is unweighted and four-channel; dev2 weighted MI is a separate binary formulation.
- [KNOWN | HIGH] HC recurrence/Neff and SHC-conditioned permutation answer different population-structure questions.
- [KNOWN | HIGH] frozen_v1 missingness distinguishes non-callable type states, whereas historical dev2 profiles often collapse non-positive values.
- [KNOWN | HIGH] EpiDis is not NMI; its squared identity is restricted to the frozen weighted empirical construction.

## Tests executed

- [COMPUTED | HIGH] Unified branch: 35 passed on Python 3.12.14.
- [COMPUTED | HIGH] CINCH-dev2 baseline: 7 passed on Python 3.12.14.

## Tests not executed

[KNOWN | HIGH] Legacy uberBlast cross-mapper equivalence, blockwise statistical equivalence, end-to-end unified execution, and production-scale benchmarks remain unexecuted.

## Next three priority actions

1. [KNOWN | HIGH] Replace the O(p²) distance DataFrame/lookup with blockwise coordinate-distance production.
2. [KNOWN | HIGH] Run the pair-engine n/p/block-size scaling ladder and profile the contingency loop.
3. [KNOWN | HIGH] Complete staged SHC/BH/ARACNE filtering and report-only commands.

## Current recommended usage

[KNOWN | HIGH] Use released frozen_v1 for its documented WGS/filter workflow and dev2 only for its explicitly frozen research modules. Do not use the integration branch as a unified production release.

## Release readiness

[COMPUTED | HIGH] Not ready for unified research use; the M03 and M05 performance/equivalence gaps and mandatory M06-M08 gates remain open.
