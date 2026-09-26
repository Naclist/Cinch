# Public SPN534 real-genome benchmark

This benchmark uses a deterministic 12-assembly subset of the 534 versioned
RefSeq accessions under `examples/spn534/provenance`. It samples folds A and B
from each of tree SHC 1–6. Raw NCBI packages and generated results are excluded
from Git; manifests, commands, checksums, and measured summaries are retained.

The historical SPN534 reference-allele FASTA is not present in the public CINCH
repository. A reference reconstructed from public NCBI annotation is therefore
an engineering/real-genome validation reference, not an exact reconstruction of
the frozen SPN534 WGS mapping.
