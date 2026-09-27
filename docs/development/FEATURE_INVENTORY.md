# Scientific feature inventory

[KNOWN | HIGH] This inventory covers significant production modules and frozen scientific modules; plotting-only helpers are grouped with their parent analysis rather than falsely presented as independent methods.

| ID | Scientific purpose | Origin and original symbol | Input contract | Output contract | Statistical definition | Unified destination | Migration | Validation |
|---|---|---|---|---|---|---|---|---|
| MAP-A1 | Internal WGS locus detection | CINCH `mapping.find_hits` | reference sequence, contigs, identity | best full-length exact-window hits | five exact seeds then Hamming identity | `cinch/mapping/backends/internal.py` | preserve then adapt | baseline tests partial |
| MAP-A2 | Genome mapping/state trace | CINCH `mapping.map_genomes` | reference FASTA, genome FASTAs | P/T arrays, trace, coordinates, QC | best hits; unique single hit callable | `cinch/workflow/mapping.py` | adapt | not cross-validated |
| MAP-B1 | Indexed HPC mapping | dev2 `Cinch_v8.nomenclature` | genome, reference DB, thresholds | mapping and allele files | uberBlast-dependent | `cinch/mapping/backends/legacy_uberblast.py` | preserve/adapter | BLOCKED dependency |
| PROF-B1 | Core/pan profiles | dev2 `gene_conservative`, `generate_profile` | mapped alleles | genome×locus profiles | prevalence and allele nomenclature | `cinch/profiles` | planned | historical only |
| STATE-A1 | Presence/type matrices | CINCH WGS | mapping hits | P∈{-1,0,1}; nominal T or -1 | callable/present state contract | `cinch/profiles` | preserve | baseline tested |
| ASSOC-A1 | PP/PT/TP/TT vectors | CINCH `channel_vectors` | P/T matrices and pair | pairwise-complete vectors | frozen four-channel definitions | `cinch/association/channels.py` | preserve | baseline tested |
| ASSOC-A2 | Unweighted contingency, MI/NMI, drivers | CINCH `contingency` | categorical vectors | MIraw/NMI/residual drivers | empirical unweighted probabilities, nats | `cinch/statistics` | preserve distinct | baseline tested |
| ASSOC-A3 | Exhaustive four-channel scan | CINCH `score_all_pairs` | P/T/loci/dense distances | channel tables | eligibility after rare pooling | `cinch/association/pairwise.py` | replace after equivalence | scalability unvalidated |
| STAT-B1 | Weighted binary MI/NMI | dev2 `weighted_binary_information` | weighted joint/marginal masses | MI, NMI, entropies | no pseudocount; nats | `cinch/statistics/weighted_mi.py` | adapted complete | zero dev2 diff; scalar reference test |
| STAT-B2 | EpiDis | dev2 `directional_epidis` | weighted binary masses | directional EpiDis | frozen weighted JSD construction | `cinch/statistics/epidis.py` | adapted complete | zero dev2 diff; identity test |
| STAT-B3 | BH adjustment | dev2 `bh_adjust` | valid p vector | q vector | stable-sort BH | `cinch/statistics/multiple_testing.py` | adapted complete | zero dev2 diff; hand test |
| STAT-B4 | Multiallelic NMI/G-test | dev2 `cgmlst_allelic_scan` kernels | integer allele matrix | pair NMI/test statistics | categorical contingency | `cinch/statistics/multiallelic_mi.py` | planned extraction | historical only |
| WEIGHT-B1 | Sample similarity/reweighting | dev2 `profile_wg_cg_similarity` | allele profile | similarity and inverse-neighbour weights | thresholded profile similarity | `cinch/population/weighting.py` | planned | historical only |
| POP-A1 | HierCC-like distance | CINCH `hiercc_distances` | allele profile | distance matrices | frozen missingness-adjusted distance | `cinch/population/hierarchy.py` | preserve | baseline tested |
| POP-A2 | HC selection | CINCH `scan_hc`, `select_hc_level` | distance matrix | assignments and three-state decision | single linkage/stability diagnostics | `cinch/population/hierarchy.py` | preserve | baseline tested |
| POP-A3 | HC recurrence/Neff | CINCH `driver_recurrence`, `neff_from_counts` | edge drivers, HC blocks | recurrence distribution, Neff | inverse Simpson across carrying blocks | `cinch/population/recurrence.py` | preserve | baseline tested |
| POP-B1 | SHC conditional MI/permutation | dev2 `hiercc_shc_phylo_correct_mi` | pair vectors, weights, SHC | CMI, empirical p/q | within-SHC exchangeability | `cinch/population/shc.py` | planned extraction | historical only |
| LINK-A1 | Order/bp distance | CINCH `pair_order_distance` | mapped coordinates | dense pair DataFrame | mean within-sample same-contig minima | `cinch/distance/SparseOrderDistance` | adapted to sparse observed-pair storage | exact controlled equality |
| LINK-B1 | Phylogenetic mixing | dev2 allele diversity workflow | profile, tree | per-locus diversity/mixing | Simpson diversity and same-allele MPD | `cinch/linkage/phylogenetic_distance.py` | planned | historical results frozen |
| NET-B1 | ARACNE pruning | dev2 `spydrpick_locus_wg_nmi.aracne_prune` | weighted edge graph | pruned graph | DPI/topology post-processing | `cinch/network/aracne.py` | planned | historical only |
| HIGH-B1 | Diff-GWES | dev2 `run_diff_gwes` | seed pairs, modifier, SHC | Δ conditional association, p/q | held-out/background comparison | `cinch/association/high_order.py` | planned | toy and frozen summaries only |
| FILT-A1 | Information envelope | CINCH `envelope`, `attach_high_information` | pair tables/order | conditional thresholds/pass | channel/order conditional quantiles | `cinch/association/eligibility.py` | preserve | release assets tested |
| FILT-A2 | Neff background filter | CINCH `neff_background`, `attach_neff_filter` | recurrent edges | final pass fields | order/channel-conditioned Neff | `cinch/population/recurrence.py` | preserve | release assets tested |
| REPORT-A1 | Frozen WGS/filter reporting | CINCH `wgs`, `filtering`, site assets | stage outputs | tables/figures/manifests | presentation only | `cinch/reporting` | adapt | release tests |
| HIST-B1 | Historical exploratory modules | dev2 `scripts/snapshots/**` | heterogeneous frozen inputs | historical result files | module-specific | `cinch/legacy/dev2` | preserve, not productionize wholesale | hash verified only |

## Dependency and path audit

- [COMPUTED | HIGH] CINCH production package has no machine-specific absolute path.
- [COMPUTED | HIGH] dev2 historical runner contains heterogeneous frozen scripts and optional non-PyPI `configure`/`uberBlast` imports.
- [COMPUTED | HIGH] frozen_v1 `threads_requested` is provenance only; `map_genomes` is serial.
- [COMPUTED | HIGH] frozen_v1 constructs all locus pairs both in `pair_order_distance` and `score_all_pairs`.
