# Mathematical specification

For the legal pairwise-complete sample subset, let `n_ab` be a contingency
count, `n` its total, `p_ab=n_ab/n`, and `p_a`, `p_b` the marginals.

`MIraw = sum_{a,b:n_ab>0} p_ab ln[p_ab/(p_a p_b)]` nats.

`H(A)=-sum_a p_a ln p_a`, and `NMI=2 MIraw/[H(A)+H(B)]`. NMI is undefined
when both entropy terms sum to zero. Zero observed cells contribute zero to MI.

A channel is eligible only when pairwise-complete `N>=20`, both marginals have
at least two observed states and every marginal state has at least the configured
minimum count after per-locus rare-state pooling. Missing calls are removed from
the legal subset and never recoded as absence.

The driver table uses `E_ab=n_a n_b/n`, Pearson residual
`R_ab=(n_ab-E_ab)/sqrt(E_ab)` and enrichment `n_ab/E_ab`. The maximum residual
is the enriched driver and the minimum residual the depleted driver.

High-information envelopes use fixed-width `log10(order_distance+1)` bins with
`min(200,max(4,floor(N/1000)))` bins per oriented channel. SPN534 therefore uses
the audited 200-bin implementation. Q95, Q99 and Q99.9 are empirical within-bin
quantiles; the primary pass is Q99.9.

For HC driver counts `n_h`, `p_h=n_h/sum_h n_h` and
`Neff=1/sum_h p_h^2`. Its primary background is the arithmetic mean Neff in the
same fixed-width log-order bin, using every finite eligible pair and pooling PT
and TP into P↔T. The strict rule is `Neff > mean`, not `>=`; Q95 is sensitivity only.

