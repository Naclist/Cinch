# Tiny WGS demo

This is a fully synthetic and redistributable 48-genome example.

Run from the repository root:

```bash
cinch wgs -r examples/tiny_wgs/input/ref.cds.fasta -p tiny -t 4 \
  --phenotype examples/tiny_wgs/input/phenotype.tsv \
  examples/tiny_wgs/input/genomes/*.fasta
cinch filter --wgs_results tiny.cinch --hc 2 --order-threshold 2
```

`input/generate_demo.py` is the source generator. `expected/` contains compact
reference outputs and the QC figure; it is not used as hidden input.

