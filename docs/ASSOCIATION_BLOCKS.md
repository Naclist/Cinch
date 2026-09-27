# Bounded PP/PT/TP/TT association engine

## Preserved scientific contract

[KNOWN | HIGH] The block engine calls frozen_v1 `pool_rare_states`, `channel_vectors`, and `contingency` directly. It preserves triangular pair order, channel orientation, eligibility masks, minimum counts, MI/NMI definitions, enriched/depleted driver tie handling, and distance columns.

[COMPUTED | HIGH] A seeded 40-sample × 6-locus exhaustive regression compared every output cell with frozen `score_all_pairs` using exact equality. All PP/PT/TP/TT frames and the rare-state audit matched.

## Bounded execution and recovery

[KNOWN | HIGH] `iter_pair_blocks` stores at most `block_pairs` integer pairs. `write_association_blocks` writes one atomic Parquet file per block and channel. A SHA-256 digest covers state arrays, labels, distance values/schema, and scoring configuration. Complete compatible blocks are skipped on rerun; input drift is rejected.

[KNOWN | HIGH] An interrupted four-file block is recomputed and atomically replaced. Existing blocks without a provenance manifest are rejected rather than trusted.

## Controlled benchmark

[COMPUTED | HIGH] At 224 samples × 80 loci (3,160 pairs; 12,640 eligible channel rows), frozen accumulation took 3.212 s and peaked at 20,266,160 traced Python bytes. Block consumption took 3.171 s and peaked at 2,991,529 traced Python bytes. The block/frozen memory ratio was 0.1476 and runtime ratio was 0.9873.

[KNOWN | HIGH] This benchmark excludes Parquet I/O and native allocations invisible to `tracemalloc`. It proves bounded Python result accumulation, not 25K-locus readiness.

## Sparse physical-distance source

[KNOWN | HIGH] Staged association uses `SparseOrderDistance`, which stores only locus pairs observed on at least one same-sample contig. Unobserved pairs are generated lazily as `NO_SAME_CONTIG_OBSERVATION`; observed pairs below the reliability count remain distinct as `INSUFFICIENT_SAME_CONTIG_OBSERVATIONS`.

[COMPUTED | HIGH] Controlled tests show exact order/bp equality with frozen_v1 for reliable pairs. The full pair hypothesis universe is unchanged; sparse storage changes representation only.

[KNOWN | HIGH] Per-pair/channel contingency construction remains a Python loop, so the 25,000-locus runtime gate is still open.
