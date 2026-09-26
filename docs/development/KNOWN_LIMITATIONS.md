# Known limitations

- [KNOWN | HIGH] The frozen internal mapper requires fixed-length windows and exact seeds; indels and distributed seed mutations can be missed.
- [KNOWN | HIGH] The frozen mapper remains serial; the opt-in indexed mapper parallelizes by genome.
- [KNOWN | HIGH] The indexed mapper assigns allele IDs by deterministic sequence-hash ordering; these IDs are not guaranteed to match historical first-observation or uberBlast nomenclature IDs.
- [KNOWN | HIGH] frozen_v1 records no-hit as confident absence even when mapper sensitivity, truncation, or ambiguity could be responsible.
- [KNOWN | HIGH] Unified indexed mapping defaults no-hit to unresolved. Explicit absence mode enables PP but assumes adequate assembly completeness and mapper sensitivity; it does not convert partial/low-confidence hits to absence.
- [KNOWN | HIGH] frozen_v1 order/bp distance materializes every locus pair; staged association uses sparse observed-pair storage.
- [KNOWN | HIGH] frozen pair scoring accumulates all rows; the unified block engine bounds output and uses sparse observed distances, but still executes Python per pair/channel.
- [KNOWN | HIGH] `configure` and `uberBlast` are not redistributed and their license/install contracts remain unresolved.
- [KNOWN | HIGH] Staged map/profile/associate, sparse distance, PP SHC/BH/ARACNE filtering, and reporting are implemented; categorical PT/TP/TT population tests, Diff-GWES workflow glue, and full-scale benchmark are not implemented.
- [COMPUTED | HIGH] Small synthetic mapping scaled poorly: 2 workers were 0.95× and 4 workers 1.14× relative to one worker; larger representative workloads remain unmeasured.
- [KNOWN | HIGH] Existing historical validation is not an unseen biological validation cohort.
- [KNOWN | HIGH] The public repository omits the historical SPN534 reference-allele FASTA and full pair cache; the public-subset E2E uses a documented NCBI CDS alternative and is not exact frozen reproduction.
