# Indexed mapping backend audit

[COMPUTED | HIGH] The active staged mapper uses mappy/minimap2 2.31. Each worker maps all reference loci against one genome index with `n_threads=1`; `ProcessPoolExecutor` limits concurrent genome workers to `min(--threads, genomes)`, preventing nested mapper-thread oversubscription.

[COMPUTED | HIGH] Cache identity includes reference SHA256, genome SHA256, every mapping configuration field, and backend version. Per-genome gzip JSON writes are atomic. The 12-genome cold run used four workers; the warm rerun reported all cache hits and reduced external wall time from 7.99 s to 2.56 s for the 2,241-CDS reference.

[COMPUTED | HIGH] Allele IDs are assigned only after worker completion, by per-locus sorting on sequence SHA256 then sequence. Serial and parallel controlled runs produce identical matrices and IDs.

[KNOWN | HIGH] `no_hit_policy=unresolved` is conservative and is the default. `no_hit_policy=absence` explicitly enables gene-content PP analysis for assemblies assumed complete enough; only zero raw hits become absence. Partial/low-confidence, competing-locus, and ambiguity states remain unresolved or present-but-type-uncallable.

[KNOWN | HIGH] Historical uberBlast equivalence remains untested because its executable/configure redistribution contract is unavailable. The indexed backend is a scientifically documented replacement candidate, not a byte/nomenclature-compatible reproduction.
