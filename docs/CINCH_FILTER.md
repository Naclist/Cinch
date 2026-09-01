# Cinch Filter v1

`cinch filter` consumes one complete WGS result directory:

```bash
cinch filter --wgs_results study.cinch/ --hc 69 --order-threshold 100
```

The publication-reproduction mode should provide both dataset-specific values.
Manual values always override automatic diagnostics and are recorded as
`USER_DEFINED`. HCx means a single-linkage allele-distance threshold, not x
clusters.

The frozen sequence is:

1. oriented-channel conditional Q99.9 MIraw envelope;
2. `order_distance >= threshold`;
3. the same enriched state-pair driver in at least three selected HC blocks;
4. `Neff` above the PP/P↔T/TT order-bin mean computed from all eligible pairs.

Q95 Neff is exported only as sensitivity information. The stage stops at
`NETWORK_READY_EDGES.tsv` and `NETWORK_READY_NODES.tsv`. ARACNE and graph
construction belong to a future `cinch network` command and are not run here.

Automatic order selection is currently `AMBIGUOUS`: no deterministic legacy
common-threshold selector could be recovered. Diagnostics are written and exit
code 3 requests `--order-threshold`. HC auto-selection is assistive and has
PASS/AMBIGUOUS/UNRESOLVED outcomes; publication runs should freeze `--hc`.

