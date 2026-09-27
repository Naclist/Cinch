# Architecture decisions

## ADR-001: Destination and history

[INFERRED | HIGH] `Naclist/Cinch` is the destination because it owns the public `cinch` command and frozen_v1 release contract. Dev2 is imported by provenance; neither source history is rewritten.

## ADR-002: Frozen code boundary

[KNOWN | HIGH] `cinch/frozen_v1` remains intact until replacements pass explicit scientific regressions. Unified modules live beside it; the `wgs` compatibility command may delegate only after equivalence evidence exists.

## ADR-003: Statistically distinct names

[KNOWN | HIGH] Unweighted four-channel MI, weighted binary MI, NMI, EpiDis, HC recurrence, SHC permutation, BH, and ARACNE remain separately named operations.

## ADR-004: Mapping backends

[INFERRED | HIGH] A small backend interface is justified because the internal mapper and legacy uberBlast mapper differ materially in sensitivity and deployability. The interface will return explicit hit/callability records; it will not normalize away discrepancies.

## ADR-005: Pairwise storage

[KNOWN | HIGH] The target 25,000 loci imply 312,487,500 unordered pairs. Unified pair computation must be blockwise, restartable, and stream outputs; dense all-pair DataFrames are rejected for production scale.

## Architecture comparison

| Component | CINCH frozen_v1 | CINCH-dev2 | Proposed unified implementation | Reason |
|---|---|---|---|---|
| WGS mapping | serial seed/Hamming full length | uberBlast + genome thread pool | backend contract; maintainable indexed backend plus legacy adapter | preserves semantics while enabling real parallelism |
| Allele calling | first-seen per-locus IDs; ambiguous multicopy uncallable | nomenclature by sequence hashes/mapping | deterministic global sequence ordering with compatibility mode | worker order must not change IDs |
| Profiles | P/T arrays | pan/core allele profiles | explicit P/T schema and historical adapter | prevents missing=absence collapse |
| Population | HierCC/HC recurrence/Neff | weights and SHC | separate population modules | methods are not interchangeable |
| Binary association | PP and mixed channels | weighted binary MI | both | distinct estimators |
| Multiallelic association | PT/TP/TT contingency | Numba NMI/G-test | compiled categorical kernel preserving channel masks | scale without changing hypotheses |
| Physical distance | dense order/bp pair table | order/bp diagnostic scripts | sparse observed distances and explicit reason codes | unknown is not far |
| Phylogenetic correction | HC recurrence | SHC permutation/tree mixing | both, explicitly selected | different inferential targets |
| Permutation/FDR | absent in WGS v1 | SHC permutation/BH | unified explicit stage | formal testing remains separate from recurrence |
| Network | presentation graph | ARACNE/community tools | post-test graph stage | topology is not significance |
| High order | absent | Diff-GWES | optional held-out high-order stage | preserve without contaminating pairwise semantics |
| Reporting | strong frozen reports/site | scientific result archive | schema-driven reports plus preserved historical assets | publication and audit continuity |
| Validation | 12 release/core tests | 7 exact/freeze tests | cross-repository golden/equivalence suite | imports alone are insufficient |
| Scalability | serial and dense pair objects | compiled kernels but monolithic scripts | block engine, cache, actual worker controls | production target requires bounded memory |

## Checkpoint 1 four-review record

- [COMPUTED | HIGH] Scientific review: no definition or pair universe has been changed; discrepancies are recorded in `SCIENTIFIC_EQUIVALENCE.md`.
- [COMPUTED | HIGH] Engineering review: both baselines install under Python 3.12 and pass their existing tests; frozen_v1 threading and dense-pair limits are confirmed.
- [COMPUTED | HIGH] Integration review: every required capability maps to a destination or preserved legacy area; no source implementation was deleted.
- [COMPUTED | HIGH] Scope review: the next critical path is exact statistical API integration, followed by mapping correctness and scalability; cosmetic work is deferred.

## Checkpoints 2–6 four-review record

- [COMPUTED | HIGH] Scientific review: block outputs are exactly equal to frozen PP/PT/TP/TT small cases; advanced kernels have zero cross-source difference; indexed mapping differences are explicit and not called equivalent to uberBlast.
- [COMPUTED | HIGH] Engineering review: 36 tests pass; mapper and pair/advanced microbenchmarks are recorded; atomic cache/output recovery is exercised.
- [COMPUTED | HIGH] Integration review: map/profile/associate share the V2 state schema; frozen wgs/filter remain intact; advanced workflow glue is still absent.
- [COMPUTED | HIGH] Scope review: production pair throughput, advanced filter/report glue, and 224×25K acceptance dominate remaining work; cosmetic repository work remains deferred.

## Staged advanced-filter review

- [COMPUTED | HIGH] Scientific review: SHC permutation is restricted to binary PP edges; BH is applied only across eligible tested edges; ARACNE occurs after significance and remains graph pruning.
- [COMPUTED | HIGH] Engineering review: deterministic edge-derived RNG seeds, aligned sample checks, input hashes, atomic outputs, and report artifacts pass the 37-test suite.
- [COMPUTED | HIGH] Integration review: `advanced-filter` consumes association blocks/profile/SHC/weights and `report` consumes immutable filter outputs; frozen `filter` is unchanged.
- [COMPUTED | HIGH] Scope review: public SPN534 real-genome E2E now outranks further local workflow embellishment.

## Public SPN534 and scaling review

- [COMPUTED | HIGH] Scientific review: real FASTA provenance is exact for the selected assemblies; the alternative CDS reference and absence policy are explicit; frozen full-result reproduction is not claimed.
- [COMPUTED | HIGH] Engineering review: all five staged operations complete; 39 tests pass; real 400-locus optimized outputs are exactly equal to the pre-prescreen outputs.
- [COMPUTED | HIGH] Integration review: SPN sample IDs align to committed SHC metadata, and mapper state semantics feed every downstream channel and PP-only SHC workflow without adapters outside the repository.
- [COMPUTED | HIGH] Scope review: per-pair Python channel dispatch remains the largest measured bottleneck; the rejected dense-code contingency experiment improved only 1.9% and is not retained.

## Optional compiled scorer review

- [COMPUTED | HIGH] Scientific review: randomized and real SPN A/B comparisons preserve pair/channel eligibility, state drivers, distance annotations, and all non-floating fields; maximum floating difference is 6.66e-16.
- [COMPUTED | HIGH] Engineering review: 43 tests pass; the dependency is optional; engine identity enters the resume digest; three-run timing and RSS are recorded.
- [COMPUTED | HIGH] Integration review: `cinch associate --engine numba` writes the same schema and consumes the same V2 profile/distance inputs; the default remains `python` and frozen commands are untouched.
- [COMPUTED | HIGH] Scope review: the later 12/24/48 ladder and compiled distance aggregation raise the measured local gain to 2.71×/2.86×/3.24×, but do not close core, high-cardinality, 534-genome, or 25K all-pair gates; completion remains 37/45.

## SPN48 and compiled-distance review

- [COMPUTED | HIGH] Scientific review: 48-set raw provenance is exact; full distance frames and non-floating association outputs match Python; the pairwise-absent global-state driver bug is covered by a regression test.
- [COMPUTED | HIGH] Engineering review: the 48-genome map→associate→filter→report chain completes; 45 tests pass; three repeats per association engine are recorded.
- [COMPUTED | HIGH] Integration review: the same staged schemas and SHC metadata scale from 12 to 48 genomes without an adapter or output-schema change.
- [COMPUTED | HIGH] Scope review: the alternative reference, 99 permutations, single-core association, 400 loci, and absent historical profile prevent a frozen/full-scale claim; completion remains 37/45.

## ADR-006: P2 standardized background estimand

[KNOWN | HIGH] Habitat/background contrasts use common-support population weights proportional to the smaller background-specific sample count in each population with adequate samples in both backgrounds. State variation is diagnostic rather than an eligibility selector, keeping the tested support invariant under background-label permutation. Empirical habitat-specific population mixtures are reported diagnostically but are not compared as the differential estimand.

## ADR-007: P2 randomization nulls

[KNOWN | HIGH] Population-conditioned association permutes channel-valid `Y` states within population and uses a one-sided MI statistic. Differential association permutes categorical background labels within population, fixes observed common support/weights, and uses a two-sided absolute delta-MI statistic. Alleles are never coerced to binary for PT/TP/TT.

## P2 four-review record

- [COMPUTED | HIGH] Scientific review: global, population-conditioned, standardized background and frozen Diff-GWES estimands remain distinct; Cases A-J pass and perfect confounding is not identifiable.
- [COMPUTED | HIGH] Engineering review: exact metadata QC, deterministic seeds, atomic outputs, complete tested/excluded audit and per-family BH are implemented; 15 P2 tests pass.
- [COMPUTED | HIGH] Integration review: `conditional-associate` consumes PROFILE_V2 and candidate association tables and emits all eight required artifacts across PP/PT/TP/TT.
- [KNOWN | HIGH] Scope review: categorical two-background inference is validated on controlled data; continuous metadata, multi-background omnibus tests and real habitat biology remain outside this P2 release.
