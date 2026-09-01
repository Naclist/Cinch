# Cinch WGS v1

`cinch wgs` accepts a reference CDS FASTA plus genome assemblies and writes one
complete WGS-result directory. It validates and hashes every input, calls
presence and nominal complete-CDS types, reconstructs within-contig gene order,
builds an independent allele-profile distance basis, scores every eligible
unordered locus pair in PP/PT/TP/TT and extracts enriched and depleted state
drivers.

Example:

```bash
cinch wgs -r ref.cds.fasta -p study -t 8 --phenotype phenotype.tsv genomes/*.fasta
```

The internal v1 mapper is a deterministic full-CDS seed-and-extend nucleotide
mapper. A unique full-length hit at or above `--min-identity` is present and
type-callable; tied distinct multicopy hits are present but type-noncallable;
no hit is confidently absent. Complete CDS byte identity defines the nominal
type. Every callable call is traced to sample, contig, start, end, strand and
CDS SHA256. Threads are recorded but do not change deterministic results.

WGS stops before HC filtering, Neff filtering, ARACNE or network construction.

