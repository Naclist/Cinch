# Unified integration changelog

## 2026-09-26 — Checkpoint 1

- [COMPUTED | HIGH] Cloned CINCH at `98c41ef081be4face6b32cb61b557964537cb353` into an independent worktree.
- [COMPUTED | HIGH] Created branch `integration/unified-framework` without rewriting source history.
- [COMPUTED | HIGH] Frozen CINCH-dev2 source at `ec1eaa5cf3c0c46a0d4d2babd52a768d23f9f57a`.
- [COMPUTED | HIGH] Added the permanent development audit system and migration manifest.
- [COMPUTED | HIGH] Executed both baseline test suites; no production behavior changed.

## 2026-09-26 — M02 exact statistical API

- [COMPUTED | HIGH] Added explicitly named weighted MI, directional EpiDis, and BH modules without altering frozen_v1 unweighted statistics.
- [COMPUTED | HIGH] Added independent scalar/hand regression tests and cross-repository numerical comparison.
- [COMPUTED | HIGH] Recorded a one-million-pair weighted-MI kernel microbenchmark.

## 2026-09-26 — M03 indexed mapping

- [COMPUTED | HIGH] Added process-parallel mappy/minimap2 mapping, reusable indexes, atomic per-genome caches, deterministic allele nomenclature, and explicit callability states.
- [COMPUTED | HIGH] Corrected frozen coordinate grouping for zero-hit and multi-locus genomes.
- [COMPUTED | HIGH] Preserved the historical dev2 `Cinch_v8.py`; legacy uberBlast equivalence remains blocked.
- [COMPUTED | HIGH] Added exact/SNP/reverse/indel/truncation/multicopy/competition tests and a controlled 1/2/4-worker benchmark.

## 2026-09-26 — M04 state schema

- [COMPUTED | HIGH] Added validated `present/absent/unresolved` and nominal-type arrays.
- [COMPUTED | HIGH] Added an explicit legacy non-call policy and frozen/V2 atomic round trips.

## 2026-09-26 — M05 bounded association

- [COMPUTED | HIGH] Added deterministic triangular blocks, atomic per-channel Parquet output, provenance digests, and true block resume.
- [COMPUTED | HIGH] Proved exact small-case equality to frozen PP/PT/TP/TT outputs and reduced traced result-accumulation peak by 85.2% in the controlled benchmark.
- [COMPUTED | HIGH] Replaced dense missing physical-distance rows with sparse observed-pair storage and explicit reliable/insufficient/unobserved states.

## 2026-09-26 — M06 advanced kernels

- [COMPUTED | HIGH] Integrated separately named SHC conditional MI/permutation, ARACNE, and Diff-GWES modifier kernels.
- [COMPUTED | HIGH] Recorded zero maximum difference across 100 seeded cross-source cases and isolated SHC/ARACNE timings.
- [COMPUTED | HIGH] Preserved the three corresponding historical dev2 scripts.

## 2026-09-26 — M07 staged workflow preview

- [COMPUTED | HIGH] Added `cinch map`, `cinch profile`, and `cinch associate` stages without changing frozen `wgs`/`filter` behavior.
- [COMPUTED | HIGH] Added a controlled profile-to-association CLI test with manifest validation and resume.
- [KNOWN | HIGH] Staged advanced filtering/reporting and production acceptance remain incomplete.

## 2026-09-26 — Staged PP filtering and reporting

- [COMPUTED | HIGH] Added deterministic PP-only SHC permutation, BH correction, ARACNE pruning, hashed provenance, and atomic outputs.
- [COMPUTED | HIGH] Added a report-only stage with tabular summary and SHC diagnostic figures.
- [KNOWN | HIGH] Binary SHC permutation was not generalized to categorical PT/TP/TT states without a validated statistical definition.

## 2026-09-26 — Public SPN534 real-genome validation

- [COMPUTED | HIGH] Downloaded 12 versioned SPN534 RefSeq assemblies spanning tree SHC 1–6 and folds A/B; all genome FASTA and GFF3 checksums match historical provenance.
- [COMPUTED | HIGH] Added explicit mapper no-hit policy so PP absence assumptions cannot be introduced silently.
- [COMPUTED | HIGH] Completed map→profile output→associate→SHC/BH/ARACNE→report on real FASTA using a documented alternative NCBI CDS reference.
- [COMPUTED | HIGH] Added a 50/100/200/400-locus real-data ladder and a 25K-locus synthetic mapper stress test.
- [COMPUTED | HIGH] Added necessary global-state prescreening; real 400-locus association improved 2.32× with exact all-channel outputs.

## 2026-09-27 — Optional compiled four-channel scorer

- [COMPUTED | HIGH] Added an opt-in Numba engine that computes all legal PP/PT/TP/TT contingency metrics per block while preserving the default Python path.
- [COMPUTED | HIGH] Added randomized multiallelic/missingness equivalence tests and retained the kernel only after the full 43-test suite passed.
- [COMPUTED | HIGH] An initial 12-genome/400-locus A/B showed a meaningful scorer gain, then the larger 48-genome test exposed a zero-margin driver defect; the defect was fixed and regression-tested before acceptance.
- [COMPUTED | HIGH] Added compiled sparse order/bp aggregation after profiling identified it as the dominant remaining association-stage cost.
- [COMPUTED | HIGH] On nested 12/24/48-genome real subsets, compiled medians improved 2.71×/2.86×/3.24×; distance/non-floating outputs were exact and maximum floating difference was 2.22e-15.
- [COMPUTED | HIGH] Completed the 48-genome map→associate→SHC/BH/ARACNE→report chain; association core scaling and 25K all-pair acceptance remain open.
