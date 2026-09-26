# Scientific equivalence register

| Method | Original formula/implementation | Planned unified implementation | Expected equivalence | Test data | Tolerance | Observed difference | Status/explanation |
|---|---|---|---|---|---|---|---|
| frozen unweighted MI/NMI | `frozen_v1.statistics.contingency`; empirical counts, nats | preserved function then compiled equivalent | exact within floating summation | tiny categorical matrices + frozen examples | `1e-12` | not yet measured cross-kernel | NOT_STARTED |
| four-channel eligibility | `channel_vectors` + `score_all_pairs` | block engine masks | identical eligible pair/channel universe | exhaustive small P/T matrices | exact identities | not measured | NOT_STARTED |
| weighted binary MI/NMI | dev2 `core.statistics.weighted_binary_information` | same implementation under `cinch.statistics` | exact | dev2 unit fixtures | `1e-12` | dev2 baseline passes | TESTED at source only |
| EpiDis | dev2 `directional_epidis` | same implementation | exact; `EpiDis²=MI_bits` only under frozen construction | dev2 perfect association/independence fixtures | `1e-12` | source test passes | TESTED at source only |
| BH | dev2 `bh_adjust` | same stable-sort implementation | exact | known p vector/ties | exact or `1e-15` | source test passes | TESTED at source only |
| HierCC distance | frozen_v1 `hiercc_distances` | preserved | exact | current core fixture | exact | baseline passes | TESTED |
| HC recurrence/Neff | frozen_v1 filtering | preserved | exact | release assets and hand counts | exact | baseline passes | TESTED |
| SHC conditional permutation | dev2 snapshot | extracted module with same exchangeability and RNG | exact indices/p-values with frozen seed | frozen summary + synthetic strata | exact p, `1e-12` CMI | not measured | NOT_STARTED |
| ARACNE | dev2 snapshot | adjacency-triangle implementation | identical deleted edge set including ties | small hand graph + frozen edges | exact | not measured | NOT_STARTED |
| Diff-GWES | dev2 snapshot | optional high-order module | identical candidate universe and Δ statistic | frozen toy truth/negative and true summary | exact counts; `1e-12` statistic | not measured | NOT_STARTED |

## Non-equivalences that must remain explicit

- [KNOWN | HIGH] Weighted MI is not equivalent to unweighted MI.
- [KNOWN | HIGH] EpiDis is not equivalent to NMI.
- [KNOWN | HIGH] HC recurrence/Neff is not a p-value or multiple-testing correction.
- [KNOWN | HIGH] SHC-conditioned permutation does not replace HierCC recurrence.
- [KNOWN | HIGH] ARACNE edge deletion is not a significance test.
- [KNOWN | HIGH] Historical non-positive dev2 profile values cannot be silently reinterpreted as the richer frozen_v1 missingness states.
