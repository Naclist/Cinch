# SPN534 real-data development snapshot

This directory contains compact outputs from the 534-genome *Streptococcus
pneumoniae* development case. It intentionally excludes the multi-gigabyte NCBI
assembly mirror and the multi-million-row all-pair Parquet cache.

- `results/FINAL_CANDIDATES.tsv`: 223 frozen candidate records.
- `results/FILTER_COUNTS.tsv`: attrition by channel and stage.
- `results/VALIDATION.tsv`: invariant/control checks.
- `figures/`: publication-ready PNG/PDF/SVG diagnostics.
- `plot_data/`: exact compact tables used by the filter diagnostics.
- `provenance/`: accession mapping and byte-level SHA256 manifests.

These outputs are development evidence, not unseen external validation and not
proof of biological epistasis.
