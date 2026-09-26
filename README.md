# Cinch

> **Research preview.** Cinch reports candidate bacterial genome dependencies.
> A retained edge is not, by itself, molecular epistasis, a fitness interaction,
> or causality.

## The problem Cinch addresses

A bacterial whole-genome collection may contain $M=10^5\text{–}10^6$ variable SNP
sites. Exhaustive pairwise screening grows as

$$
N_{\mathrm{SNP\ pairs}}=\binom{M}{2}=\frac{M(M-1)}{2}=O(M^2),
$$

creating an enormous, sparse and difficult-to-interpret search space. At the same
time, bacterial genomics already provides cgMLST/wgMLST profiles that group sequence
variation into biologically named loci. Cinch uses that locus-level abstraction:

$$
\text{SNP-scale variation}
\xrightarrow{\text{cg/wgMLST locus abstraction}}
\text{nominal locus states}
\xrightarrow{\text{information dependence}}
\text{candidate dependency graph}.
$$

This is a candidate-discovery framework: it compresses the search from nucleotide
coordinates to interpretable loci before evaluating dependence and graph structure.

## One locus, two biological state variables

For isolate $i$ and locus $g$, Cinch records

$$
P_{ig}\in\{-1,0,1\},\qquad
T_{ig}\in\{a_1,a_2,\ldots\}\cup\{-1\},
$$

where $P$ is non-callable/absent/present and $T$ is a nominal complete-CDS
type when present. Type IDs are categories, not ordered numbers, SNP distances or
branch lengths. For loci $A,B$, the four tested channels are

$$
\begin{aligned}
PP&: I(P_A;P_B),\\
PT&: I(P_A;T_B\mid P_B=1),\\
TP&: I(T_A;P_B\mid P_A=1),\\
TT&: I(T_A;T_B\mid P_A=P_B=1).
\end{aligned}
$$

PT and TP remain separate computational orientations, but may be displayed jointly
as P↔T. The TT channel is what allows two nearly universal genes to have
$MI_{PP}\approx0$ but strong allele-type dependence.

## Association and its state driver

For state variables $X,Y$, Cinch reports mutual information in nats:

$$
MI_{\rm raw}(X;Y)=\sum_{x,y}p(x,y)
\ln\frac{p(x,y)}{p(x)p(y)},
\qquad
NMI=\frac{2MI_{\rm raw}}{H(X)+H(Y)}.
$$

MI/NMI quantify dependence magnitude, not sign or direction. The enriched state
cell is described separately through the observed/expected contingency residual:

$$
R_{xy}=\frac{O_{xy}-E_{xy}}{\sqrt{E_{xy}}},
\qquad
E_{xy}=\frac{O_{x\cdot}O_{\cdot y}}{N}.
$$

## Local-linkage evidence

For loci reliably co-located on the same contig, the frozen order distance is

$$
d_{\rm order}(A,B)=
\operatorname{mean}_{i\in\mathcal C_{AB}}
\left[
\min_{u\in A_i,v\in B_i}|\operatorname{order}(u)-\operatorname{order}(v)|
\right],
$$

requiring at least five same-contig observations. Different-contig pairs remain
NA; they are never assigned an artificial large distance. Each channel has its
own empirical conditional envelope $Q_{99.9}(MI\mid d_{\rm order})$. The frozen
SPN534 filter requires both an exceptional channel-specific MI and
$d_{\rm order}\ge100$ genes.

![Fine-resolution mean MIraw versus distance](docs/assets/figures/FIG06_DISTANCE_MEAN_CURVES.png)

**Take-home:** local order/bp dependence decays toward a low baseline, whereas the
legacy phylogenetic curve has a distinct far-distance tail. The plotted point area
is proportional to $\sqrt N$, so tail bins do not masquerade as dense evidence.

Across all 3,955,078 SPN534 locus pairs, legacy phylogenetic distance is available
for 3,624,696; order and bp are available for the same 990,661 pairs. When both are
defined, $\rho_{\rm Spearman}(d_{order},d_{bp})=0.9973$. Order is retained not for
greater coverage, but because locus units normalize genome-scale length, are less
sensitive to variable intergenic expansion, and transfer more directly across
assemblies/species. The full density and missingness audit is in
[`FIG06_DISTANCE_COMPARISON`](docs/assets/figures/FIG06_DISTANCE_COMPARISON.png).

## Population recurrence and Neff

HC69 is a single-linkage allele-difference threshold, not “69 clusters” and not a
model producing adjusted p/q-values. Let $n_h$ be support for the same enriched
state-pair driver in HC69 block $h$:

$$
K_{HC}=\sum_h\mathbf 1(n_h>0),\qquad
p_h=\frac{n_h}{\sum_j n_j},\qquad
N_{eff}=\frac{1}{\sum_h p_h^2}.
$$

$K_{HC}$ measures recurrence breadth; $N_{eff}$ measures its evenness. Three
blocks with counts $(98,1,1)$ and $(34,33,33)$ have the same $K_{HC}=3$, but
very different $N_{eff}$. The **current frozen SPN534 rule** is

$$
K_{HC}\ge3
\quad\land\quad
N_{eff}>E[N_{eff}\mid \operatorname{bin}(d_{order}),\text{channel}].
$$

The frozen threshold is therefore stated explicitly wherever this SPN534 result
is shown; later design decisions are kept in the development-history document.

## One edge from state table to final graph

The retained TT example follows the complete frozen chain:

$$
\boxed{
gyrB\!-!aguA:
P\approx1
\Rightarrow MI_{PP}\approx0
\Rightarrow MI_{TT}=2.1602
\Rightarrow NMI=0.7983
\Rightarrow type_{46}\leftrightarrow type_{51}
\Rightarrow d_{order}=120.29
\Rightarrow (n_h)=(1,3,3)
\Rightarrow N_{eff}=2.5789
\Rightarrow \text{retained}
}
$$

The negative development contrast is

$$
rsuA\_2\!-!glyS:
MI_{TT}\ \text{high}
\Rightarrow (n_h)=(5)
\Rightarrow K_{HC}=1,\ N_{eff}=1
\Rightarrow \text{rejected}.
$$

![TT background contrast](docs/assets/figures/FIG10_TT_BACKGROUND_CONTRAST.png)

**Take-home:** high TT information is insufficient when the driver is only a
population identity tag; recurrent support changes the interpretation.

## Published comparator and information beyond P/A

Among 985 high-MI pairs in the Coinfinder audit, 982 (99.70%) are direct published
association/dissociation edges and 963 (97.77%) lie in the same published component.
This is concordance with a published dependency method, not causal ground truth or
an unseen accuracy estimate. The V-ATPase control further shows that a dependency
may be detected correctly and then interpreted as local linkage.

![Coinfinder concordance](docs/assets/figures/FIG12_COINFINDER_CONCORDANCE.png)

**Take-home:** Cinch recovers the established P/A dependency structure while the
P↔T/TT channels define its additional information domain.

![PP versus TT](docs/assets/figures/FIG13_PP_VS_TT.png)

**Take-home:** the $MI_{PP}\approx0,MI_{TT}\gg0$ region contains dependencies that
a presence/absence-only analysis cannot represent.

## Frozen SPN534 result

SPN534 is a 534-genome *Streptococcus pneumoniae* **development case**, not an
unseen validation cohort. Frozen settings are HC69, $K_{HC}\ge3$, and a 100-gene
long-range threshold.

| stage | PP | PT | TP | TT |
|---|---:|---:|---:|---:|
| eligible | 2,902,845 | 1,180,977 | 2,109,466 | 67,853 |
| high information | 646 | 606 | 853 | 76 |
| order filtered | 309 | 303 | 448 | 42 |
| HC69 recurrent | 50 | 71 | 160 | 1 |
| Neff filtered | **4** | **58** | **160** | **1** |

The result contains 223 oriented records and 223 unique unordered locus pairs.
The published tree, HC69 strip and the exact P/T states for the four PP plus one TT
examples use the same 534-isolate order.

![Tree-aligned states](docs/assets/figures/FIG14_TREE_HC_STATE.png)

**Take-home:** retained state drivers are visible in their phylogenetic-background
context rather than being presented as context-free scores.

## Candidate dependency network

The frozen table forms a 167-locus, 223-edge **candidate dependency network**.
It is pre-ARACNE: the frozen configuration has `ARACNE: false` and
`network_construction: false`; graph construction here is a presentation of the
network-ready table, not another statistical filter. Communities are descriptive
graph partitions, not automatically pathways.

![Candidate dependency network](docs/assets/figures/FIG16_FINAL_NETWORK.png)

**Take-home:** most retained candidates occupy one large component, with a small
set of disconnected dependency components. Use the
**[interactive network](examples/spn534/network/CINCH_INTERACTIVE.html)** for
channel filters, ego views, product search and exact edge metadata.

## Install and run the redistributable demo

The actual WGS input contract is a reference CDS FASTA plus genome FASTA files:

```bash
git clone https://github.com/Naclist/Cinch.git
cd Cinch
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
python -m pip install -e ".[test]"
pytest -q

cinch wgs \
  -r examples/tiny_wgs/input/ref.cds.fasta \
  -p tiny -t 4 \
  --phenotype examples/tiny_wgs/input/phenotype.tsv \
  examples/tiny_wgs/input/genomes/*.fasta

cinch filter --wgs_results tiny.cinch --hc 2 --order-threshold 2
```

The toy's `--hc 2` is a toy-specific HC level, not the SPN534 recurrence-count
criterion. SPN534 provenance, versioned assembly mappings and SHA256 manifests are
under [`examples/spn534/provenance`](examples/spn534/provenance).

## Unified staged integration preview

The integration branch also exposes explicit restart boundaries:

```bash
cinch map -r reference.cds.fasta -o mapped -t 4 genomes/*.fasta
cinch profile legacy.tsv -o converted --missing-policy unresolved
cinch associate \
  --profile mapped/profiles/PROFILE_V2.npz \
  --coordinates mapped/mapping/COORDINATES.tsv \
  -o associated --block-pairs 100000

# optional compiled scorer
python -m pip install -e ".[performance]"
cinch associate \
  --profile mapped/profiles/PROFILE_V2.npz \
  --coordinates mapped/mapping/COORDINATES.tsv \
  -o associated-numba --block-pairs 100000 --engine numba
```

`map` caches each indexed genome independently. `profile` refuses to guess whether
legacy non-calls mean absence or unresolved. `associate` preserves frozen
PP/PT/TP/TT semantics and writes atomic resumable blocks. `--engine numba` is an
optional tolerance-equivalent compiled scorer; the default remains Python. This
staged path is not yet a production replacement: categorical advanced testing and
production-scale pair throughput are still incomplete. See [`docs/STAGED_WORKFLOW.md`](docs/STAGED_WORKFLOW.md)
and [`docs/development/COMPLETION_README.md`](docs/development/COMPLETION_README.md).

## Repository map and boundaries

| path | purpose |
|---|---|
| [`cinch/frozen_v1`](cinch/frozen_v1) | state construction, MI/NMI, hierarchy, order distance, recurrence and Neff |
| [`cinch/mapping`](cinch/mapping) | indexed genome mapper, deterministic alleles and per-genome cache |
| [`cinch/profiles`](cinch/profiles) | explicit presence/type schema and legacy adapters |
| [`cinch/association`](cinch/association) | bounded four-channel blocks and high-order kernels |
| [`cinch/population`](cinch/population) | SHC-conditioned MI and permutation kernels |
| [`cinch/network`](cinch/network) | ARACNE triangle pruning |
| [`docs/SCIENTIFIC_NARRATIVE.md`](docs/SCIENTIFIC_NARRATIVE.md) | full mathematical/scientific walkthrough |
| [`docs/FIGURE_INDEX.md`](docs/FIGURE_INDEX.md) | result figures and exact evidence source |
| [`docs/DEVELOPMENT_HISTORY.md`](docs/DEVELOPMENT_HISTORY.md) | decisions, failures and frozen/future boundaries |
| [`examples/tiny_wgs`](examples/tiny_wgs) | fully redistributable demo |
| [`examples/spn534`](examples/spn534) | frozen tables, diagnostics, tree states, network and checksums |
| [`site`](site) | GitHub Pages source |
| [`FROZEN_CINCH_WGS_V1.yaml`](FROZEN_CINCH_WGS_V1.yaml) | machine-readable freeze contract |

- MI/NMI measures dependence, not causality.
- HC recurrence is not an adjusted significance test.
- Order distance reduces obvious local linkage; it does not prove function.
- ARACNE, if used later, is post-filter redundancy pruning and cannot define significance.
- SPN534 informed development and must not be presented as unseen validation.

Comparator inputs originate from
[fwhelan/coinfinder-manuscript](https://github.com/fwhelan/coinfinder-manuscript)
and Whelan, Rusilowicz & McInerney (2020),
[doi:10.1099/mgen.0.000338](https://doi.org/10.1099/mgen.0.000338).

No redistribution license has yet been selected; public source visibility is not
a grant of reuse rights.
