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
| 2026-09-26 | unified M04 | Python 3.12.14 | `pytest -q` | 27 passed in 1.55 s | state invariants, explicit legacy missing policy, V2/frozen round trips |
| 2026-09-26 | unified M05 | Python 3.12.14 | `pytest -q` | 30 passed in 1.67 s | exact frozen/block equality, triangular universe, resume/input drift |
| 2026-09-26 | unified M05 | 224 samples × 80 loci | `benchmark_pair_blocks.py` | 3.212→3.171 s; traced peak 20.27→2.99 MB; row counts equal | computational benchmark; Parquet I/O excluded |
| 2026-09-26 | unified M06 vs dev2 `ec1eaa5` | 100 seeded cases | direct cross-source comparison | max abs diff 0; ARACNE and stable seeds exact | SHC/MI/high-order migration |
| 2026-09-26 | unified M06 | Python 3.12.14 | `pytest -q` | 34 passed in 1.61 s | hand cases and advanced-kernel invariants |
| 2026-09-26 | unified M06 | controlled synthetic kernels | `benchmark_advanced_statistics.py` | SHC 2,805 permutations/s; ARACNE 64,831 edges/s | isolated single-run timing |
| 2026-09-26 | unified M07 | Python 3.12.14 | `pytest -q` | 36 passed in 1.63 s | profile conversion → sparse distance → association CLI, manifest and resume |
| 2026-09-26 | unified M06/M07 | Python 3.12.14 | `pytest -q` | 37 passed in 3.38 s | staged PP SHC permutation, BH, ARACNE and report artifacts |
| 2026-09-26 | public SPN534 subset | 12 versioned assemblies | NCBI Datasets 18.37.0 + SHA256 audit | 12/12 genome FASTA and 12/12 GFF3 exact historical matches | raw-input provenance |
| 2026-09-26 | public SPN534 subset | 12 genomes × 200 alternative CDS | staged CLI | map→profile output→associate→advanced-filter→report completed | real FASTA E2E; not frozen reproduction |
| 2026-09-26 | public SPN534 subset | 12 genomes × 400 CDS | old/new association A/B | PP/PT/TP/TT frames exact; 13.80→5.95 s | necessary-state prescreen equivalence and performance |
| 2026-09-26 | synthetic stress | 1 genome × 25,000 loci | `benchmark_mapper.py` | 17.86 s external wall; 302.4 MB RSS; 25,000 callable | engineering only |
| 2026-09-26 | unified current | Python 3.12.14 | `pytest -q` | 39 passed in 2.05 s | full local suite after SPN/no-hit changes |
| 2026-09-27 | optional compiled association | Python 3.12.14, NumPy 2.3.5, Numba 0.63.1 | randomized and real SPN A/B | all non-floats exact; max float difference 6.66e-16 | three random multiallelic/missingness cases plus 79,800 real pairs |
| 2026-09-27 | public SPN534 ladder | nested 12/24/48 genomes × 400 alternative CDS; three repeats/engine | staged `associate` | 2.71×/2.86×/3.24× medians | exact distance/non-float outputs; max float difference 2.22e-15 |
| 2026-09-27 | public SPN534 subset | 48 genomes × 400 alternative CDS | staged CLI | map→associate→advanced-filter→report completed | 48/48 FASTA and GFF3 historical SHA256 matches; alternative reference |
| 2026-09-27 | unified current | Python 3.12.14, NumPy 2.3.5, Numba 0.63.1 | `pytest -q` | 45 passed in 2.45 s | full suite including pairwise-absent-state and compiled-distance regressions |
| 2026-09-27 | P2 conditional association | deterministic controlled Cases A-J plus support/audit regressions | `pytest -q tests/test_conditional_association.py` | 15 passed in 1.31 s | population confounding, habitat coexistence/exclusion/null, overlap, TT/PT/TP, missingness, imbalance, non-identifiability CLI output, permutation-invariant support, stable audit schemas, legacy Diff-GWES |
| 2026-09-27 | unified current with P2 | Python 3.12.14, optional Numba installed | `pytest -q` | 60 passed in 2.74 s | complete repository regression after P2 integration |
| 2026-09-27 | v0.1.0 release candidate | Python 3.12.14 | `pytest -q` | 61 passed in 2.42 s | package-version regression added; no scientific method changed |
| 2026-09-27 | installed-command smoke | local editable and isolated wheel installs | `scripts/release_smoke.py` | PASS twice | tiny WGS/filter counts 36/32/31/21 and 33 final; conditional TT 1 row/2 tests |
| 2026-09-27 | v0.1.0 distributions | isolated build and wheel venv | `python -m build`; artifact verifier; `twine check`; `pip check` | PASS | package modules, selected docs, license notices, console entry point and SPN534 exclusion verified |
| 2026-09-27 | public GitHub clone at `6166c96` | fresh Python 3.12.14 environments | public HTTPS clone; `pytest -q`; both CLI smokes; build; artifact/Twine checks; wheel install | 62 passed; all release checks PASS | no local source checkout or private input used |
| 2026-09-27 | GitHub Actions `ae0a319` | hosted Ubuntu, Python 3.10/3.12 | workflow run `36284486133` | all three jobs PASS | regression matrix plus package build, artifact audit, isolated wheel install and CLI smoke |

## Evidence limitations

[KNOWN | HIGH] Tests and executed benchmarks validate process-based mapping, resume, mapping edge cases, deterministic nomenclature, exact/tolerance-bounded blockwise equivalence, staged filtering/reporting, public SPN-subset E2E, and a 25K-locus mapper stress case. They do not validate legacy uberBlast equivalence, full 534-genome biological sensitivity/specificity, sample/core scaling, or 25K-locus all-pair execution.
