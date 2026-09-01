# HC selection

The frozen profile distance follows the audited development implementation. For
two samples, jointly callable mismatches and overlap are counted. A 3% missing
allowance sets `L=max(callable_i,callable_j)-0.03*G`; when `L>overlap`, both the
effective mismatch and overlap are increased by `L-overlap`. The integer
distance is `round(G * effective_mismatch/effective_overlap)`.

HC assignments cut a SciPy single-linkage hierarchy at distance x. Consequently
HCx is an allele-distance connectivity level and never means x clusters. Scans
report cluster sizes, singletons, adjacent NMI/ARI, silhouette and within/between
distance. The assistive selector seeks the smallest point of a unique longest
stable adjacent-partition plateau; ties are AMBIGUOUS and absence is UNRESOLVED.
Manual `--hc` always overrides it.

