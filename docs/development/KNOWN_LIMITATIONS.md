# Known limitations

- [KNOWN | HIGH] The frozen internal mapper requires fixed-length windows and exact seeds; indels and distributed seed mutations can be missed.
- [KNOWN | HIGH] `--threads` does not currently control mapping computation.
- [KNOWN | HIGH] Allele IDs are assigned on first observation and would become nondeterministic if naive parallel completion order were introduced.
- [KNOWN | HIGH] frozen_v1 records no-hit as confident absence even when mapper sensitivity, truncation, or ambiguity could be responsible.
- [KNOWN | HIGH] order/bp distance currently materializes every locus pair and cannot scale to 25K loci.
- [KNOWN | HIGH] pair scoring currently accumulates Python dictionaries/Pandas rows for all eligible pairs.
- [KNOWN | HIGH] `configure` and `uberBlast` are not redistributed and their license/install contracts remain unresolved.
- [KNOWN | HIGH] The unified stage CLI, resumable cache, block engine, advanced statistics, and full-scale benchmark are not implemented.
- [KNOWN | HIGH] Existing historical validation is not an unseen biological validation cohort.
