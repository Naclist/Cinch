# Order distance

Within each sample and deposited contig, mapped features are ordered by genomic
start. For a locus pair with duplicate placements, the minimum absolute order
difference is used for that sample-contig. The reported pair distance is the
arithmetic mean across sample-contig observations when at least five exist.
Absent loci, different contigs and insufficient support are NA. Different
contigs are never assigned an artificial long distance.

SPN534's 100-gene threshold is dataset-specific and manually frozen. No exact
deterministic legacy common elbow selector was recoverable, so automatic order
selection remains AMBIGUOUS rather than silently inventing a replacement.

