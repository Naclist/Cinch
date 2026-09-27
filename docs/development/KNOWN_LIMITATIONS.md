# Known limitations

- [KNOWN | HIGH] The frozen internal mapper requires fixed-length windows and exact seeds; indels and distributed seed mutations can be missed.
- [KNOWN | HIGH] The frozen mapper remains serial; the opt-in indexed mapper parallelizes by genome.
- [KNOWN | HIGH] The indexed mapper assigns allele IDs by deterministic sequence-hash ordering; these IDs are not guaranteed to match historical first-observation or uberBlast nomenclature IDs.
- [KNOWN | HIGH] frozen_v1 records no-hit as confident absence even when mapper sensitivity, truncation, or ambiguity could be responsible.
- [KNOWN | HIGH] Unified indexed mapping defaults no-hit to unresolved. Explicit absence mode enables PP but assumes adequate assembly completeness and mapper sensitivity; it does not convert partial/low-confidence hits to absence.
- [KNOWN | HIGH] frozen_v1 order/bp distance materializes every locus pair; staged association uses sparse observed-pair storage.
- [KNOWN | HIGH] frozen pair scoring accumulates all rows; the unified block engine bounds output and uses sparse observed distances. Its optional Numba scorer and distance aggregator have been measured only to 48 genomes × 400 loci; association core scaling and 25K all-pair execution remain untested.
- [KNOWN | HIGH] `configure` and `uberBlast` are not redistributed and their license/install contracts remain unresolved.
- [KNOWN | HIGH] Staged map/profile/associate, sparse distance, PP SHC/BH/ARACNE filtering, reporting, and categorical PP/PT/TP/TT population/background tests are implemented. P2 supports one explicit categorical two-background contrast per run; continuous metadata and multi-background omnibus tests are not implemented.
- [COMPUTED | HIGH] Small synthetic mapping scaled poorly: 2 workers were 0.95× and 4 workers 1.14× relative to one worker; larger representative workloads remain unmeasured.
- [KNOWN | HIGH] Existing historical validation is not an unseen biological validation cohort.
- [KNOWN | HIGH] The public repository omits the historical SPN534 reference-allele FASTA and full pair cache; the public-subset E2E uses a documented NCBI CDS alternative and is not exact frozen reproduction.
- [KNOWN | HIGH] P2 q-values control the declared testing family only. A candidate file selected using the same association statistic can introduce selection bias; CINCH records that family but does not claim broader post-selection FDR control.
- [KNOWN | HIGH] Background-label permutation assumes conditional exchangeability within the selected population column. CINCH diagnoses overlap but cannot prove that unmeasured within-population confounding is absent.
- [KNOWN | HIGH] P0 25K-locus and P1 full-SPN534 acceptance are user-deferred, not validated.
