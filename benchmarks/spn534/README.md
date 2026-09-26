# Public SPN534 real-genome benchmark

This benchmark uses a deterministic 12-assembly subset of the 534 versioned
RefSeq accessions under `examples/spn534/provenance`. It samples folds A and B
from each of tree SHC 1–6. Raw NCBI packages and generated results are excluded
from Git; manifests, commands, checksums, and measured summaries are retained.

The historical SPN534 reference-allele FASTA is not present in the public CINCH
repository. A reference reconstructed from public NCBI annotation is therefore
an engineering/real-genome validation reference, not an exact reconstruction of
the frozen SPN534 WGS mapping.

`select_public_subset.py` creates nested deterministic 12/24/48-assembly ladders
by selecting one/two/four lexicographically ordered samples from every
`tree_shc × fold` stratum in the versioned public manifest. The selection is for
engineering scaling coverage; it is not a random biological cohort.

Recreate the 48-assembly selection and raw package with:

```bash
python benchmarks/spn534/select_public_subset.py \
  --manifest examples/spn534/provenance/TABLE_EXT01_SAMPLE_MANIFEST.tsv \
  --per-stratum 4 \
  --manifest-output benchmarks/spn534/subset48_manifest.tsv \
  --accessions-output benchmarks/spn534/accessions48.txt

datasets download genome accession \
  --inputfile benchmarks/spn534/accessions48.txt \
  --include genome,gff3 \
  --filename benchmarks/spn534/downloads/spn48.zip
```

The executed package SHA256 was
`fd232af7420e62705c3426607932d88a08cc76dcbc1d62e0fe67ab0b65d21ee1`.
All 48 genome FASTA and 48 GFF3 hashes matched
`examples/spn534/provenance/ASSEMBLY_FILE_SHA256.tsv`. Timings, output counts,
environment versions, and equivalence tolerances are frozen in
`benchmarks/results/spn534_public_subset_numba_20260927.json`.
