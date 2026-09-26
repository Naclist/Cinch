# Public SPN534 validation record

## Input audit

[COMPUTED | HIGH] The repository contains 534 versioned RefSeq accessions, a 534-row sample/SHC/fold mapping, a 4,272-file source checksum manifest, compact frozen result tables, a tree, and publication diagnostics. It does not contain raw assemblies, the historical reference-allele FASTA, the complete allele profile, raw coordinates, or the multi-million-row pair cache.

[COMPUTED | HIGH] Twelve assemblies were selected deterministically: folds A and B from each tree SHC 1–6. NCBI Datasets 18.37.0 downloaded genome FASTA, GFF3, CDS, and sequence reports. All 12 genome FASTA and all 12 GFF3 files exactly match the historical committed SHA256 values; no mismatch occurred.

## Reference boundary

[KNOWN | HIGH] The historical reference-allele FASTA could not be recovered from the public CINCH snapshot. The real-genome run therefore uses deterministic CDS subsets from the NCBI annotation of `GCF_900037015.1`. This is an alternative public reference and cannot reproduce frozen SPN534 WGS calls exactly.

## Successful real-genome workflow

[COMPUTED | HIGH] With 12 real assemblies and 200 alternative-reference CDS loci, explicit `no_hit_policy=absence` produced PP=120, PT=1,393, TP=1,415, and TT=16,048 eligible rows from 19,900 unordered pairs. SHC filtering tested 120 PP candidates; four met the small-subset permutation eligibility rule, and none passed BH at 0.05 with 99 permutations. Reporting completed and produced TSV, Markdown, PNG, PDF, and SVG artifacts.

[KNOWN | HIGH] `no_hit_policy=absence` is explicit because PP requires distinguishable absence. The default remains `unresolved`; partial or low-confidence alignments remain unresolved under both policies. Absence mode assumes assembly completeness and adequate mapper sensitivity.

## Frozen comparison boundary

[COMPUTED | HIGH] The frozen SPN534 result has 2,902,845/1,180,977/2,109,466/67,853 eligible PP/PT/TP/TT rows and 223 final records under its historical full reference, 534 genomes, HC69, and order threshold 100. The new 12-genome alternative-reference counts are not directly comparable and are not claimed as reproduction.

[KNOWN | HIGH] Frozen PP/PT/TP/TT formulas, state drivers, order/bp aggregation, HC recurrence, and Neff remain regression-tested at kernel or frozen-asset level. Exact full SPN534 result reproduction requires the omitted historical reference/profile/pair inputs and identical configuration.

## Larger public subset and optimization evidence

[COMPUTED | HIGH] A deterministic nested 12/24/48-assembly ladder selects one/two/four samples from every tree-SHC × fold stratum. For the 48-set, all 48 genome FASTA and all 48 GFF3 files match the historical SHA256 manifest. Eight-worker mapping against the same 400-CDS alternative reference completed in 7.31 s wall time, used 41.05 CPU-seconds, and peaked at 373.0 MB RSS.

[COMPUTED | HIGH] At 48 samples, the compiled A/B retained 1,016/7,226/10,155/76,038 PP/PT/TP/TT rows. The 66,621-row physical-distance frame, state drivers, identifiers, eligibility, and all other non-floating fields were exact; maximum floating difference was `2.22e-15`. Three-run median time fell from 8.74 to 2.70 s and median RSS from 458.2 to 386.3 MB.

[COMPUTED | HIGH] The 48-genome advanced stage processed 1,016 PP candidates with 99 permutations: 371 were permutation-eligible, 79 passed BH at 0.05, and ARACNE retained 50 direct edges. Filter and report took 3.96 and 1.41 s. These counts use the alternative reference and reduced permutation count, so they are not frozen biological results.

## Four-review record

- [COMPUTED | HIGH] Scientific: absence semantics are explicit; low-confidence calls are not converted to absence; no frozen reproduction claim is made.
- [COMPUTED | HIGH] Engineering: map, profile serialization, associate, advanced-filter, and report complete on 48 real FASTA; downloaded genome/GFF provenance matches; 45 tests pass after compiled kernel/distance integration.
- [COMPUTED | HIGH] Integration: the real run exercises indexed mapping, deterministic alleles, sparse distances, all four channels, SHC/BH/ARACNE, and reporting.
- [COMPUTED | HIGH] Scope: pairwise runtime was the largest measured bottleneck and was reduced locally; broader sample/core scaling, full SPN534, 25K all-pair, and private Bacillus HPC acceptance remain separate open gates.
