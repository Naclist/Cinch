# SPN534 provenance

## Published source

The deposited presence/absence matrix, core tree and Coinfinder outputs are from:

F. J. Whelan, M. Rusilowicz and J. O. McInerney (2020), *Coinfinder: detecting
significant associations and dissociations in pangenomes*, Microbial Genomics
6:e000338, DOI: [10.1099/mgen.0.000338](https://doi.org/10.1099/mgen.0.000338).

Repository: <https://github.com/fwhelan/coinfinder-manuscript>

## Independent assembly provenance

`provenance/PUBLICATION_TO_EXISTING_ASSEMBLY.tsv` maps each publication sample to
an exact versioned NCBI assembly record and records paired accessions, BioSample,
study, assembly status and annotation release.

`provenance/ASSEMBLY_FILE_SHA256.tsv` contains the relative source path, byte
count and SHA256 for every FASTA, GFF3 and GBFF used. The checksums preserve the
exact byte identity without redistributing the large source files.

`provenance/PUBLICATION_SOURCE_SHA256.tsv` hashes the deposited publication
inputs. `provenance/TABLE_EXT01_SAMPLE_MANIFEST.tsv` records the analysis sample
IDs and folds used in the earlier external-validation audit.

The accession manifest is sufficient to retrieve source records, but NCBI may
later update a current annotation. Reproduction should verify the committed
SHA256 and treat a checksum mismatch as a provenance failure, not silently use
new bytes.

