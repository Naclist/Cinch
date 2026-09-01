# Development results and interpretation

## Synthetic demonstration

The 48-genome synthetic example checks the representation rather than claiming
realistic biological calibration. It includes:

- a PP co-occurrence pair;
- a P-to-T dependency;
- two always-present loci with a TT dependency and PP MI near zero;
- a short-range linked pair intended to fail the order filter;
- a population-confounded pair intended to fail cross-background recurrence.

The expected result directory is a compact audit artifact. Regenerating the run
from `examples/tiny_wgs/input` is the primary executable demonstration.

## SPN534 development snapshot

The raw eligible-pair universe contains:

| Channel | Eligible pairs | High-information | Long-range | HC recurrent | Final |
|---|---:|---:|---:|---:|---:|
| PP | 2,902,845 | 646 | 309 | 50 | 4 |
| PT | 1,180,977 | 606 | 303 | 71 | 58 |
| TP | 2,109,466 | 853 | 448 | 160 | 160 |
| TT | 67,853 | 76 | 42 | 1 | 1 |

The frozen development choice uses HC69 and an order distance of 100 genes. The
final set contains 223 oriented records and 223 unique unordered locus pairs.
The asymmetry in PT versus TP arises from orientation-specific callability and
state entropy; plots may merge them as P↔T, but the calculation retains the
direction.

The retained set is strongly enriched for P↔T candidates. This is evidence that
the unified state representation exposes cross-scale information that a PP-only
workflow cannot represent. It is not evidence that every retained pair is a
functional interaction.

The TT example `gyrB`--`aguA` has MIraw 2.1602104671 nats, NMI 0.7982993337,
mean order distance 120.2857 genes, support in three HC69 backgrounds and Neff
2.5789473. This is a useful illustration of allele-type dependence after the
frozen recurrence and order rules, not a causal claim.

`rsuA_2`--`glyS`, which was attractive in exploratory plots, is not in the final
frozen set. Retaining negative outcomes is part of the development record.

## What the figures show

- `ORDER_ELBOW.png`: channel-specific order-distance background scan.
- `HC_STABILITY.png`: hierarchy stability across allele-difference thresholds.
- `NEFF_DIAGNOSTIC.png`: recurrence concentration versus the order-conditioned background.
- `FILTER_ATTRITION.png`: counts remaining after each frozen filter.
- `FINAL_CANDIDATE_LANDSCAPE.png`: final candidates by channel and score context.
- `NMI_VS_ORDER.png` and `NEFF_VS_ORDER.png`: score/background geometry with final candidates highlighted.

SPN534 was examined repeatedly while Cinch was developed. It must therefore be
reported as a real-data development case with orthogonal published comparison,
not as a one-shot external validation.

