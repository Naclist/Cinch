# Known limitations

- [KNOWN | HIGH] The frozen internal mapper requires fixed-length windows and exact seeds; indels and distributed seed mutations can be missed.
- [KNOWN | HIGH] The frozen mapper remains serial; the opt-in indexed mapper parallelizes by genome.
- [KNOWN | HIGH] The indexed mapper assigns allele IDs by deterministic sequence-hash ordering; these IDs are not guaranteed to match historical first-observation or uberBlast nomenclature IDs.
- [KNOWN | HIGH] frozen_v1 records no-hit as confident absence even when mapper sensitivity, truncation, or ambiguity could be responsible.
- [KNOWN | HIGH] order/bp distance currently materializes every locus pair and cannot scale to 25K loci.
- [KNOWN | HIGH] pair scoring currently accumulates Python dictionaries/Pandas rows for all eligible pairs.
- [KNOWN | HIGH] `configure` and `uberBlast` are not redistributed and their license/install contracts remain unresolved.
- [KNOWN | HIGH] `cinch map`, its resumable per-genome cache, and explicit profile adapters are implemented; the block engine, advanced statistics, full staged CLI, and full-scale benchmark are not implemented.
- [COMPUTED | HIGH] Small synthetic mapping scaled poorly: 2 workers were 0.95× and 4 workers 1.14× relative to one worker; larger representative workloads remain unmeasured.
- [KNOWN | HIGH] Existing historical validation is not an unseen biological validation cohort.
