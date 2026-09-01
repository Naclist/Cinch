# Development history and frozen decisions

## Unified state model

The original information-theory core emphasized binary gene content. The reboot
made P and nominal CDS type T explicit, producing PP/PT/TP/TT without SNPs. The
first checkpoint required a structure-free toy to recover core-core TT dependence
that P/A cannot represent.

## Population development

Inverse-neighbour weighting reduced effective information and harmed power on the
development toy. SHC reduced phylogenetic decoys but not local linkage. Same-data
CMI screening followed by permutation/BH raised selection-induced-bias concerns.
These results prevented those modules from becoming one opaque corrected score.

The current design uses HC69 as cross-background recurrence evidence, not a p/q
generator. Profile distance remains QC context, not a rule that more distant means
more credible.

## Linkage and hierarchy

Different contigs were never assigned an artificial huge distance. Verified GFF
coordinates restored bp and order diagnostics. Order was retained as the primary
linkage axis; SPN534's 100-gene threshold is frozen and dataset-specific.

HC4 means a four-allele single-linkage threshold, not four clusters. HC69 and
HC218 were inspected; HC69 was frozen for SPN534 recurrence. Neff was then added
for all-pair background-support dispersion against order-matched expectations.

## Current freeze and future boundary

The snapshot is PP 4, PT 58, TP 160, TT 1 (223 total). Future runs must retain
edges recurrent in exactly two populations (`>=2`); this freeze is unchanged.
Development data may guide choices, but validation inputs must be preregistered
and used once. If validation results change the method, that dataset becomes
development data and a new unseen validation is required.
