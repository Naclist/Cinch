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

## Remaining scaling defect

[KNOWN | HIGH] The current API still receives a materialized all-pair distance DataFrame and builds a keyed lookup, both O(p²). A streaming/block distance provider must replace this before the full 25,000-locus target can pass its memory gate.
