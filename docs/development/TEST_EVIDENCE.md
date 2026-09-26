# Test evidence

| Date | Source/commit | Environment | Command | Result | Scope |
|---|---|---|---|---|---|
| 2026-09-26 | CINCH `98c41ef` | Python 3.12.14 local venv | `pytest -q` | 12 passed in 5.88 s | core semantics and frozen release assets |
| 2026-09-26 | CINCH-dev2 `ec1eaa5` | Python 3.12.14 dev2 venv | `pytest -q tests` | 7 passed | weighted MI/EpiDis/BH and frozen hashes |
| 2026-09-26 | CINCH `98c41ef` | system Python 3.9.6 | editable install | failed before tests | environment below declared `>=3.10`; not a code failure |
| 2026-09-26 | unified M02 | Python 3.12.14, NumPy 2.5.3 | `pytest -q` | 17 passed in 1.31 s | frozen tests plus weighted MI/EpiDis/BH |
| 2026-09-26 | unified M02 vs dev2 `ec1eaa5` | 10,000 seeded valid inputs | cross-repository comparison | max abs diff 0 for all three APIs | numerical migration equivalence |
| 2026-09-26 | unified M02 | 1,000,000 pairs × 5 | `benchmark_weighted_mi.py` | best 0.083142542 s; 12,027,537 pairs/s | kernel only |
| 2026-09-26 | unified M03 | Python 3.12.14, mappy 2.31 | `pytest -q` | 24 passed in 1.54 s | frozen/statistics plus controlled indexed mapping and resume workflow |
| 2026-09-26 | unified M03 | 12 synthetic genomes × 200 exact loci | `benchmark_mapper.py --genomes 12 --loci 200 --length 300 --threads 1 2 4` | 0.534/0.560/0.468 s; 4-worker speedup 1.14× | computational benchmark only |

## Evidence limitations

[KNOWN | HIGH] Existing tests validate process-based mapping, per-genome cache reuse, indels, truncations, competing loci, strand/coordinate handling, and deterministic nomenclature. They do not validate legacy uberBlast installation/equivalence, biological sensitivity/specificity, 25K-locus memory, blockwise numerical equivalence, or end-to-end integration.
