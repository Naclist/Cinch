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
| M03 | Mapper backends and resumability | TESTED | public SPN12 genomes; resume; 25K synthetic reference; audit | legacy uberBlast equivalence and broader biological validation |
| M04 | Unified profile/state model | VALIDATED | explicit schema, legacy policy adapter, frozen/V2 round trips | none |
| M05 | Blockwise four-channel engine | TESTED | exact A/B, sparse distance, real 50–400-locus ladder, 2.32× prescreen gain | n/core scaling and full 25K pair run |
| M06 | Advanced statistics | TESTED | SHC/ARACNE/Diff-GWES kernels plus PP SHC/BH/ARACNE workflow | workflow-scale profile |
| M07 | Staged CLI/end-to-end workflow | VALIDATED | 12 public SPN534 assemblies, real FASTA E2E, stage timings and reports | none within M07; production acceptance remains M08 |
| M08 | Production acceptance | BLOCKED | target workload recorded | data and HPC execution |

[COMPUTED | HIGH] Predefined completion is `37/45 = 82.2%`; see `MASTER_PLAN.md`.

## Scientific feature preservation

| Feature | Original implementation | Unified implementation | Validation status |
|---|---|---|---|
| PP/PT/TP/TT | `cinch/frozen_v1/statistics.py` | unchanged frozen_v1 | baseline tested only |
| HC recurrence/Neff | `filtering.py` | unchanged frozen_v1 | baseline tested only |
| Weighted MI/EpiDis/BH | dev2 `core/statistics.py` | `cinch/statistics` | zero max absolute difference over 10,000 random valid inputs; unified tests pass |
| SHC permutation/ARACNE/Diff-GWES | dev2 frozen snapshots | kernels plus staged PP SHC/BH/ARACNE filter | PP workflow tested; Diff-GWES remains kernel-level |
| Indexed mapping/allele nomenclature | frozen mapper plus dev2 `Cinch_v8.py` | mappy/minimap2 index, deterministic hash-sorted alleles, atomic per-genome cache | controlled cases pass; historical uberBlast equivalence unresolved |

## Confirmed performance improvements

[COMPUTED | HIGH] Unified weighted MI processed 1,000,000 valid mass triplets at a best observed 12,027,537 pairs/s over five runs on the local Python 3.12.14/NumPy 2.5.3 environment. This is a kernel microbenchmark, not end-to-end throughput.

[COMPUTED | HIGH] At 224 samples × 80 loci, lazy blocks reduced traced peak Python allocation from 20.27 MB to 2.99 MB with exact row equality. Runtime changed from 3.212 s to 3.171 s; no material speedup is claimed.

[COMPUTED | HIGH] On 12 public SPN534 assemblies, the scientifically safe global-state prescreen reduced 400-locus association wall time from 13.80 s to 5.95 s (2.32×) with exact PP/PT/TP/TT DataFrame equality. At 200 loci it reduced 4.23 s to 2.41 s (1.76×).

[COMPUTED | HIGH] A synthetic 25,000-locus × 1-genome mapping stress test completed in 17.86 s external wall time at 302.4 MB maximum RSS. This is engineering evidence only, not biological validation or full pairwise acceptance.

## Unverified performance assumptions

[COMPUTED | HIGH] On 12 synthetic genomes × 200 exact loci, the indexed mapper completed in 0.534 s with one worker and 0.468 s with four workers, a measured 1.14× speedup. Two workers were slower than one. This workload is dominated by process/index overhead and is not biological validation.

[COMPUTED | HIGH] Blockwise enumeration and sparse distance storage are implemented and measured. Compiled all-channel block scoring remains an unimplemented candidate; a local dense-code contingency rewrite was rejected after only 1.9% real-data improvement.

## Current blockers

- [KNOWN | HIGH] `uberBlast` and `configure` redistribution/installability are unresolved.
- [KNOWN | HIGH] The 224-genome/25,000-reference production dataset is not present in this repository.
- [KNOWN | HIGH] Full 312,487,500-pair execution has not been resource-estimated; sparse distance and block output exist, but per-pair/channel Python scoring remains a runtime blocker.

## Known scientific discrepancies

- [KNOWN | HIGH] frozen_v1 MI is unweighted and four-channel; dev2 weighted MI is a separate binary formulation.
- [KNOWN | HIGH] HC recurrence/Neff and SHC-conditioned permutation answer different population-structure questions.
- [KNOWN | HIGH] frozen_v1 missingness distinguishes non-callable type states, whereas historical dev2 profiles often collapse non-positive values.
- [KNOWN | HIGH] EpiDis is not NMI; its squared identity is restricted to the frozen weighted empirical construction.

## Tests executed

- [COMPUTED | HIGH] Unified branch: 39 passed on Python 3.12.14.
- [COMPUTED | HIGH] CINCH-dev2 baseline: 7 passed on Python 3.12.14.

## Tests not executed

[KNOWN | HIGH] Legacy uberBlast equivalence, 534-genome/full-reference reproduction, 25K full-pair execution, and private Bacillus HPC acceptance remain unexecuted. A public 12-genome SPN534-subset E2E passes with a documented alternative reference.

## Next three priority actions

1. [KNOWN | HIGH] Design compiled multi-channel block scoring to remove remaining per-pair Python dispatch, retaining exact outputs.
2. [KNOWN | HIGH] Extend real SPN534 subsets in sample count and reference size within local resources.
3. [KNOWN | HIGH] Execute full 534-genome/public-reference or historical-reference validation when the required reference/profile becomes available.

## Current recommended usage

[KNOWN | HIGH] Use released frozen_v1 for its documented WGS/filter workflow and dev2 only for its explicitly frozen research modules. Do not use the integration branch as a unified production release.

## Release readiness

[COMPUTED | HIGH] The staged public-subset workflow is validated for research-preview use, but the repository is not production-ready because M03, M05, M06 performance, and M08 acceptance gates remain open.
