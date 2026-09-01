# Cinch WGS v1 / Filter v1 freeze report

## 1. Frozen commands

```bash
cinch wgs -r ref.cds.fasta -p RUN -t THREADS [--phenotype phenotype.tsv] genomes/*.fasta
cinch filter --wgs_results RUN.cinch/ --hc HC_LEVEL --order-threshold ORDER_DISTANCE
```

## 2. WGS inputs

Required: reference CDS FASTA, genome FASTAs and prefix. Optional phenotype is
validated and frozen only. Mapping identity, minimum informative N, marginal
state support and same-contig order support are exposed CLI parameters. Exact
formats and defaults are in `docs/INPUTS.md` and CLI help.

## 3. WGS outputs

The complete file contract, purpose and row semantics are enumerated in
`docs/OUTPUTS.md`. Massive states and pair records use Parquet/NPZ; compact QC,
configuration and summaries use TSV/YAML/Markdown. Every input is byte-hashed.

## 4. Filter input

`--wgs_results` expects `00_manifest/RUN_MANIFEST.yaml`, the four `04_pairs`
Parquets, state/profile basis and all relative paths named by the manifest. Raw
users never provide internal tables manually.

## 5. Filter outputs

Order and HC scans, Q99.9 candidates, order-filtered rows, HC driver
distributions, all-pair Neff backgrounds, final candidates, network-ready
edges/nodes, plots, manifests, hashes, validation and logs are emitted below
`RUN.cinch/filter/`. No ARACNE or graph is emitted.

## 6. Real-time metrics

WGS prints input/sample/locus counts, mapped presence and type-callable counts,
profile loci and exact-profile unique N, per-channel eligible counts and MI
ranges, driver count, finite order-distance count, background-bin count and
runtime. Filter prints manifest dimensions, order bins and selection status, HC
levels and selection status, per-channel Q99.9 counts, order attrition,
HC-recurrence count, Neff-final count and network-ready total.

Tiny example excerpt:

```text
[WGS-08] PP: eligible=36; MI range=0.0956026..0.693147
[WGS-08] PT: eligible=32; MI range=0..0.693147
[WGS-08] TP: eligible=31; MI range=0..0.693147
[WGS-08] TT: eligible=21; MI range=0..1.38629
[FILTER-07] order >= 2: 55 -> 49
[FILTER-08] same-driver recurrence in >= 3 HC blocks: 33
[FILTER-09] Neff above channel/order-conditioned all-pair mean: final=33
```

## 7. User-exposed parameters

Universal definitions: natural-log MI, symmetric NMI, channel legality,
missingness semantics, arithmetic-mean order aggregation, single-linkage HC,
driver-count Neff and strict above-mean Neff rule. Dataset-specific parameters:
mapping identity, eligibility supports, HC level and order threshold. SPN534's
HC69 and 100 genes are not universal defaults.

## 8. HC implementation

HC uses the audited 3%-missingness-adjusted integer allele distance and SciPy
single linkage. HCx is threshold x, not x clusters and not a complete-linkage
diameter constraint. Scans report adjacent NMI/ARI, silhouette, cluster sizes,
singletons and within/between distance. Auto selection is assistive; manual HC
is the publication path.

## 9. Order threshold

The distance is the arithmetic mean of within-sample/contig minimum absolute
gene-order differences, with at least five observations. Different contigs are
NA. Q95/Q99/Q99.9 are fixed-log-width bin statistics. No recoverable
deterministic common legacy elbow algorithm existed, so automatic common order
selection is explicitly AMBIGUOUS; manual values override it.

## 10. Phenotype interface

Sample IDs and columns are validated, exact bytes and SHA256 are preserved. No
phenotype GWAS, edge selection or network overlay is implemented in v1.

## 11. Tiny demo

`demo/tiny_wgs/run_wgs.sh`, `run_filter_auto.sh` and
`run_filter_manual.sh` provide exact commands. The raw 48-genome demo completes
in seconds. The manual run yields PP 12, PT 18, TP 0, TT 3, total 33; intended
PP, PT and PP=0/TT-high controls are recovered while local and single-population
controls fail.

## 12. SPN534 demo

`demo/SPN534/run_filter.sh` prepares the self-contained frozen WGS adapter and
runs HC69/order-100 filtering. Exact final counts are PP 4, PT 58, TP 160,
P↔T 218, TT 1, total 223 oriented and 223 unique undirected pairs. The four PP
smoke tests and `gyrB--aguA` survive; `rsuA_2--glyS` fails. The Q99.9 values
recomputed from all eligible pairs agree with frozen values within 1e-16 and
finite-order pass membership is identical.

## 13. Reproducibility

`FREEZE_SHA256.tsv`, both manifests and input SHA tables record exact bytes.
Validated runtime: Python 3.12.13, NumPy 2.3.5, pandas 3.0.1, SciPy 1.18.1,
scikit-learn 1.9.0, Matplotlib 3.11.1 and Biopython 1.88. The complete repository
test run is 52 passed.

## 14. Status

FROZEN / REPRODUCIBLE

