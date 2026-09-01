# SPN534 real-data development snapshot

This directory contains compact outputs from the 534-genome *Streptococcus
pneumoniae* development case. It intentionally excludes the multi-gigabyte NCBI
assembly mirror and the multi-million-row all-pair Parquet cache.

- `results/FINAL_CANDIDATES.tsv`: 223 frozen candidate records.
- `results/FILTER_COUNTS.tsv`: attrition by channel and stage.
- `results/VALIDATION.tsv`: invariant/control checks.
- `figures/`: original frozen PNG/PDF/SVG diagnostics.
- `plot_data/`: exact compact plotting tables, including PP-versus-TT and
  distance-bin summaries.
- `diagnostics/`: exact legacy distance definitions, distance availability,
  Coinfinder concordance and the V-ATPase linkage control.
- `tree_state/`: published-tree-aligned HC69 and P/T state tracks for the frozen
  four PP plus one TT examples.
- `network/`: the complete 167-node/223-edge frozen graph as TSV, GraphML, GEXF
  and an interactive Cytoscape.js page.
- `provenance/`: accession mapping and byte-level SHA256 manifests.

These outputs are development evidence, not unseen external validation and not
proof of biological epistasis. The checked-in 223-edge table is pre-ARACNE; the
interactive network must not be described as an ARACNE-pruned significance set.
