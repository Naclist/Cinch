# P2 population- and habitat-dependent association status

## Priority boundary

[KNOWN | HIGH] P0 25K-locus acceptance and P1 full 534-genome SPN534 acceptance are deferred to the user. Their tests and documentation remain present, but they are not counted against P2 scientific completion and are not marked validated.

## P2 acceptance table

| Gate | Status | Evidence | Remaining limitation |
|---|---|---|---|
| Mathematical specification | VALIDATED | `docs/CONDITIONAL_ASSOCIATION.md` | conditional exchangeability remains a study-design assumption |
| Population-conditioned MI | VALIDATED | exact controlled population-only case | unweighted empirical population mass only |
| Habitat-specific MI | VALIDATED | common-support standardized controlled cases | categorical backgrounds only |
| Differential MI | VALIDATED | coexistence, exclusion and null cases | two-background contrast per run |
| PP inference | VALIDATED | all-channel CLI E2E | biological validation requires user dataset |
| PT inference | VALIDATED | separate conditioned channel test | same limitation |
| TP inference | VALIDATED | separate conditioned channel test | same limitation |
| TT inference | VALIDATED | three-state categorical case | high-cardinality performance not assessed under P2 |
| Permutation validity | VALIDATED | fixed channel mask; Y-within-population and background-within-population tests | exchangeability must be justified by the study |
| Population-habitat overlap | VALIDATED | perfect-confounding and unbalanced-design tests | no model-based extrapolation beyond common support |
| Multiple testing | VALIDATED | complete tested/excluded audit and per-family BH | q-values apply only to the declared candidate family |
| Metadata integration | VALIDATED | exact ID QC, duplicate/missing/extra reporting | selected population/background values cannot be missing |
| CLI integration | VALIDATED | required eight-file output E2E | no continuous-background method |
| Frozen Diff-GWES compatibility | VALIDATED | exact result equality to preserved snapshot kernel | frozen EpiDis contrast remains a separate mode/API |
| End-to-end validation | VALIDATED | 15 P2 tests; four-channel CLI case | controlled method validation, not biological discovery |
| Scientific documentation | VALIDATED | estimand, null, overlap, FDR and interpretation boundaries documented | study-specific covariates still require design review |

[COMPUTED | HIGH] P2 controlled-method completion is `16/16`. This does not make the whole repository production-accepted; P0/P1 remain user-deferred and no real habitat dataset has yet supplied a biological validation claim.

## Mandatory checkpoint reviews

| Checkpoint | Scientific review | Engineering review | Integration review | Scope review |
|---|---|---|---|---|
| 1 Mathematical specification | global, empirical-CMI and common-support standardized delta are distinct | formulas are executable kernels | frozen Diff-GWES boundary retained | categorical two-background extension only |
| 2 Metadata/schema | exact sample alignment and selected-column completeness required | QC is written before inference failure | consumes PROFILE_V2 sample order | optional metadata columns remain optional |
| 3 Conditional MI | nominal categories and four channel masks preserved | controlled scalar/population cases pass | same profile states as staged association | unweighted estimator is not weighted MI/EpiDis |
| 4 Differential association | same population weights used for both backgrounds | coexistence, exclusion, null and multiallelic cases pass | candidate pair adapter accepts global association tables | no causal interpretation |
| 5 Permutation | nulls and sidedness are explicit; missingness mask fixed | deterministic per-hypothesis seeds | applies to PP/PT/TP/TT without binary coercion | validity depends on conditional exchangeability |
| 6 Multiple testing | BH families are declared and complete | tested and excluded hypotheses are both serialized | no ARACNE/significance conflation | external candidate selection limits the q-value family |
| 7 CLI | parameters expose population, contrast, channels and support | atomic outputs and failed-QC manifest implemented | `cinch conditional-associate` uses existing profile schema | continuous covariates rejected by absence of a continuous mode |
| 8 End-to-end | Cases A-J pass, including non-identifiability and legacy equality | required eight outputs generated; full suite passes | production CLI path is exercised | controlled validation only; P0/P1 remain deferred |
