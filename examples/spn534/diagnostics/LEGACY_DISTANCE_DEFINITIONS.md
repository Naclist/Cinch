# SPN534 original distance diagnostic definitions

## MIraw
Direct byte-preserving read of `track_B_retrospective/results/all_pair_descriptive_scores.parquet::raw_mi`.
No association statistic was recomputed.

## phylo_dist
Exact historical `Leca/code/Cinch_v8.py::compute_vectorized_phylo_dist` definition:
`(X.T @ D_tree @ X)[A,B] / (n_joint[A,B] * (n_joint[A,B]-1))`.
`D_tree` is the published-tree patristic-distance matrix in Newick branch-length units; X is the
frozen binary publication P/A matrix. The numerator spans all A-carriers x B-carriers while the
denominator uses joint carriers, exactly as implemented historically. It is undefined for
`n_joint < 2`. This is not the recent cg/wg profile carrier-distance or geometry distance.

## order_dist
For every exact `(sample, contig)` co-observation, take the minimum absolute difference between
the full-annotation `gene_order` values across duplicate placements. Report the arithmetic mean
across sample-contig observations only when at least 5 observations
exist. Different contigs and insufficient support are NA. Units are annotated genes.

## bp_dist
For every exact `(sample, contig)` co-observation, interval distance is
`max(0, max(start_A,start_B) - min(end_A,end_B) - 1)`. Duplicate placements use the minimum.
Report the arithmetic mean across sample-contig observations only when at least
5 observations exist. Different contigs are never assigned a distance.
Replicon identity is unavailable in this draft-assembly dataset. Units are bp.
