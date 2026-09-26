# SHC, ARACNE, and Diff-GWES kernels

## Integrated definitions

[KNOWN | HIGH] `cinch.population.shc` preserves dev2 sample-weighted binary MI, SHC-conditional MI, informative-cluster classification, and within-SHC permutation with the `(1 + exceedances)/(B + 1)` p-value.

[KNOWN | HIGH] `cinch.network.aracne` searches triangles in the observed adjacency graph. It deletes a uniquely weakest edge and preserves exact/near ties under the historical `np.isclose(rtol=1e-12, atol=1e-12)` rule.

[KNOWN | HIGH] `cinch.association.high_order` preserves the Diff-GWES modifier statistic: within each modifier state it mass-averages MI over informative SHCs, converts nats to `EpiDis=sqrt(MI/ln(2))`, and reports `delta=EpiDis(C=1)-EpiDis(C=0)`.

## Migration evidence

[COMPUTED | HIGH] Across 100 seeded random binary datasets, unified weighted MI, conditional MI, informative-SHC classification, and high-order statistics had maximum absolute difference `0.0` from `Cinch-dev2@ec1eaa5`. A controlled triangle had exact ARACNE frame equality, and stable seeds matched exactly.

[COMPUTED | HIGH] The unified test suite contains hand-calculable MI/CMI cases, deterministic permutation, ARACNE unique-minimum/tie cases, EpiDis identity checks, positive modifier delta, and within-cluster permutation invariants.

## Controlled timing

[COMPUTED | HIGH] For 1,000 samples, 20 SHCs, and 199 permutations, one edge required 0.07095 s (2,805 permutations/s). For a 500-node, 1,000-edge sparse graph, ARACNE required 0.01542 s (64,831 input edges/s).

[KNOWN | HIGH] These are isolated single-run kernel timings. They exclude block I/O, multiprocessing, many-edge permutation scheduling, candidate generation, and full high-order searches.

## Not yet integrated

[KNOWN | HIGH] No production CLI stage yet connects blockwise pair results to SHC permutation, BH correction, ARACNE, or discovery/validation Diff-GWES. The kernels are tested, but workflow-level checkpoint M06 remains open.
