# Population- and habitat-conditioned association

## Scope and interpretation

[KNOWN | HIGH] This workflow tests context-dependent statistical association. It does not by itself establish molecular epistasis, adaptive epistasis, a fitness interaction, or causality.

[KNOWN | HIGH] The four state channels retain their existing definitions: PP uses callable presence at both loci; PT conditions the right locus on presence and uses its nominal allele type; TP applies the mirrored conditioning; TT requires both loci present with callable nominal types. Missing, unresolved, and non-callable states are excluded, not encoded as biological categories.

## Estimands

For channel-specific categorical variables `X` and `Y`, population `C`, and categorical background `H`:

```text
global_MI = I(X;Y)
population_conditioned_MI = sum_c P(c) I(X;Y | C=c)
```

[KNOWN | HIGH] The population-conditioned quantity uses the empirical population mass among the channel's eligible samples.

For reference background `h0` and comparison background `h1`, let `S` contain populations in which both backgrounds meet the configured sample minimum. Define common-support weights:

```text
m_c = min(n[c,h0], n[c,h1])
pi_c = m_c / sum_{k in S} m_k
I_std(h) = sum_{c in S} pi_c I(X;Y | C=c,H=h)
delta_MI = I_std(h1) - I_std(h0)
```

[KNOWN | HIGH] Both backgrounds therefore use the same population distribution. A raw difference based on incompatible empirical population mixtures is not used as evidence of habitat-dependent association.

## Identifiability

[KNOWN | HIGH] A differential hypothesis is `NOT_IDENTIFIABLE` when either background lacks the configured total support, too few populations have adequate samples in both backgrounds, or no positive common support remains. State variation is reported separately but does not select the shared populations because that selection would change under background-label permutation. A monomorphic margin has a defined MI of zero. No effect estimate or `p=1` substitute is emitted for an unidentifiable hypothesis.

Reported diagnostics include background counts, population counts, shared sample-supported populations, per-cell state variation, common sample support, effective population count `1/sum(pi_c^2)`, and observed joint-state-cell counts.

## Null hypotheses and permutations

[KNOWN | HIGH] The population association test holds `X`, population labels, and channel eligibility fixed and permutes `Y` within population. Its one-sided statistic is population-conditioned MI because MI is non-negative.

[KNOWN | HIGH] The differential test holds `X`, `Y`, population labels, channel eligibility, shared-population set, and standardization weights fixed. It permutes background labels within population, preserving the observed background counts in every population. The two-sided statistic is `abs(delta_MI)`. This is a conditional-exchangeability test; it is not valid when unmodelled within-population variables make habitat labels non-exchangeable.

Both tests use:

```text
p = (1 + count(null statistic >= observed statistic)) / (B + 1)
```

[KNOWN | HIGH] Degenerate permutation tables are assigned their mathematically defined MI, including zero for a monomorphic margin. Missing states are never permuted into callable states because permutations occur only after the channel-specific eligible mask is fixed.

## Multiple testing

[KNOWN | HIGH] BH correction is applied separately to declared families: population-conditioned tests are grouped by channel and population column; differential tests are grouped by channel, population column, background column, and ordered background contrast. Every tested or excluded hypothesis is written to an audit table. Eligibility filtering is separated from significance.

[KNOWN | HIGH] When `--pairs` supplies a screened candidate set, q-values describe that declared candidate family only. They are not represented as FDR control over every possible locus pair. Full-pair mode and candidate-input mode use the same output audit schema.

## Relationship to other CINCH quantities

- [KNOWN | HIGH] Unweighted MI/NMI describe association strength for the current nominal states.
- [KNOWN | HIGH] Weighted MI and EpiDis are separate binary estimators and are not substituted into this categorical workflow.
- [KNOWN | HIGH] HC recurrence and Neff describe recurrence breadth/evenness across HierCC backgrounds; they are not p/q-values.
- [KNOWN | HIGH] SHC-conditioned association tests residual dependence after population stratification.
- [KNOWN | HIGH] Habitat-conditioned and differential MI ask whether association strength changes across an explicitly selected categorical background.
- [KNOWN | HIGH] ARACNE is post-test graph redundancy pruning and does not create significance.

[KNOWN | HIGH] The array-level standardized-background kernel accepts any nominal categorical background, including a binary genetic state or environmental category. The `conditional-associate` CLI is the metadata-table adapter; frozen binary genetic-background Diff-GWES remains the separate legacy-compatible API below.

[KNOWN | HIGH] Globally lineage-associated pairs are allowed as conditional-analysis candidates. They are not discarded merely because global MI follows population structure, and they are not called significant unless the corresponding conditional test supports that claim.

## Frozen Diff-GWES boundary

[KNOWN | HIGH] Frozen Diff-GWES remains a separate legacy-compatible statistic:

```text
delta_frozen = sqrt(I(A;B|SHC,C=1)/ln(2)) - sqrt(I(A;B|SHC,C=0)/ln(2))
```

It uses a binary genetic background, held-out discovery/validation selection, and within-SHC background permutation. The new habitat `delta_MI` neither renames nor replaces this EpiDis contrast.
