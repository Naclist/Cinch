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
