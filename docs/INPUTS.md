# Inputs

## `cinch wgs`

- `-r/--reference`: FASTA of complete reference CDS; record IDs are locus IDs.
- positional `genomes`: one or more genome FASTA assemblies; basename is sample ID.
- `-p/--prefix`: output prefix (`PREFIX.cinch`).
- `-t/--threads`: requested thread count, recorded in provenance.
- `--phenotype`: optional table whose first column is sample ID. v1 validates and freezes it only.
- `--min-identity`: full-length nucleotide mapping identity, default 0.90.
- `--min-informative`: pairwise-complete sample minimum, default 20.
- `--min-state-count`: marginal state minimum and rare-pooling cutoff, default 3.
- `--min-order-observations`: same-contig observations required for distance, default 5.

## `cinch filter`

- `--wgs_results`: directory containing `00_manifest/WGS_MANIFEST.yaml` and all paths it names.
- `--hc`: dataset-specific allelic-distance threshold.
- `--order-threshold`: dataset-specific long-range gene-order threshold.
- `--minimum-populations`: same-driver HC recurrence count, frozen default 3.

