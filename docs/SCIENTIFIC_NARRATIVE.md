# Scientific narrative

## The search-space problem

For `M` variable nucleotide sites, exhaustive SNP-pair screening requires
`M(M-1)/2 = O(M^2)` tests. With `M` in the hundreds of thousands or millions,
the space is computationally large and fragmented across hard-to-interpret
nucleotide coordinates. Cinch uses cg/wgMLST as a gene/locus-level abstraction:

```text
SNP-scale variation -> nominal locus states -> information dependence
                    -> candidate dependency graph
```

Presence/absence is then recognized as only one locus resolution. It is useful
for accessory genes but cannot resolve two nearly universal loci. Cinch therefore
makes P and nominal CDS type T explicit and applies one information-theory core
across PP, PT, TP and TT.

## What Cinch measures

MIraw is mutual information of the observed nominal state table in nats. NMI
rescales by available marginal entropy. Neither value contains a sign or causal
direction. Enriched/depleted state cells are retained as driver descriptors.

## Three distinct evidential axes

1. MIraw/NMI measures state dependence.
2. Gene-order distance and its channel-specific envelope address local linkage.
3. HC69 recurrence asks whether the same driver repeats across backgrounds;
   Neff measures how evenly that support is distributed.

The accepted SPN534 configuration uses HC69 as a recurrence filter, requires the
same driver in at least three HC69 blocks, and uses 100 genes as a manually frozen
long-range threshold. HC69 is a 69-allele single-linkage threshold; it does not
mean 69 groups.

## Frozen result

The filter retains 223 oriented records: 4 PP, 58 PT, 160 TP and 1 TT. PT/TP are
kept distinct in the table because tested state orientation differs, while plots
may group both as P↔T. The graph contains 167 loci and 223 unordered pairs.

The retained TT edge `gyrB–aguA` has MIraw 2.1602 nats, NMI 0.7983, mean order
distance 120.29 genes, three supporting HC69 backgrounds and Neff 2.5789. Its
dominant state driver is `type_46 ↔ type_51`. These facts motivate a hypothesis;
they do not establish direct molecular interaction.

## What the real-data case does and does not validate

SPN534 provides published pangenome inputs, a published comparator, versioned
assembly/GFF provenance, a local-linkage control and a population-background audit.
It also informed method choices and is therefore development data, not an unseen
validation cohort. Coinfinder concordance does not make published edges ground
truth. V-ATPase dependence can be real statistically while best interpreted as
local physical linkage.

## Frozen boundary

The checked-in 223-edge table stops before ARACNE. ARACNE may be applied only as
post-filter redundancy pruning, with pre/post networks preserved. It must never
participate in significance or error-control definitions.
