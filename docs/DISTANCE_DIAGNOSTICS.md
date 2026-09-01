# Distance diagnostics

Exact definitions are in
[`LEGACY_DISTANCE_DEFINITIONS.md`](../examples/spn534/diagnostics/LEGACY_DISTANCE_DEFINITIONS.md);
numeric results are in
[`TABLE_DISTANCE_COMPARISON.tsv`](../examples/spn534/diagnostics/TABLE_DISTANCE_COMPARISON.tsv).

| distance | valid pairs | fraction of 3,955,078 |
|---|---:|---:|
| legacy phylogenetic | 3,624,696 | 91.65% |
| gene order | 990,661 | 25.05% |
| bp | 990,661 | 25.05% |

Order and bp use the same exact same-contig coordinate observations and minimum
five-observation rule. There are no order-only or bp-only pairs. Different-contig
pairs have no numerical bp/order distance and remain NA.

Descriptive Spearman correlations with MIraw are 0.0879 (phylogenetic), -0.0733
(order) and -0.0782 (bp). With millions of pairs, p-values are not biologically
useful; effect size and the full distance envelope are the diagnostics.

Among jointly mapped pairs, order and bp have rho=0.9973. Order is retained not
for greater availability but because locus units are less sensitive to variable
intergenic length and transfer more cleanly across fragmented assemblies.
