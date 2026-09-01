# Tiny demo

```bash
cinch wgs -r examples/tiny_wgs/input/ref.cds.fasta -p tiny -t 4 \
  --phenotype examples/tiny_wgs/input/phenotype.tsv \
  examples/tiny_wgs/input/genomes/*.fasta
cinch filter --wgs_results tiny.cinch --hc 2 --order-threshold 2
```

The 48-genome/20-locus demo includes PP association and dissociation, P↔T and
TT dependence, a PP=0/TT-high core pair, two ambiguous type calls, a local pair,
a one-population decoy and cross-population recurrence. The frozen manual run
produces 33 oriented candidates (PP 12, PT 18, TP 0, TT 3), including
`pp_co_A--pp_co_B`, `pt_presence--pt_type`, and `tt_core_A--tt_core_B`; local and
population-confined controls fail.
