# Figure index

Result figures are available as PNG, PDF and SVG under `docs/assets/figures`.
The mathematical definitions in the README and scientific narrative are the
authoritative method overview.

| figure | question | evidence |
|---|---|---|
| FIG06 | How do mean MIraw curves and phylo/order/bp availability compare? | fixed-width-bin lead view + 3,955,078-pair audit |
| FIG07 | What are channel-specific order envelopes? | frozen summaries |
| FIG08 | Why HC69? | frozen HC scan |
| FIG09 | How is Neff derived? | exact formula/examples |
| FIG10 | Single-background versus recurrent TT? | frozen states/drivers |
| FIG11 | How does Neff vary with order? | frozen all-pair summary |
| FIG12 | Does Cinch recover Coinfinder structure? | 985-pair audit |
| FIG13 | What does TT add beyond PP? | 67,853 pairs |
| FIG14 | Where do frozen PP/TT states lie on the tree? | 534 isolates |
| FIG15 | How do filters narrow each channel? | frozen counts |
| FIG16 | What is the full candidate graph? | frozen 223 edges |
| FIG17 | How is the graph partitioned? | deterministic communities |

The build script asserts that `FINAL_CANDIDATES.tsv` has identical SHA256 before
and after presentation generation.
