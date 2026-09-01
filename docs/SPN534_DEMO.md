# SPN534 compact development demo

The public repository contains the compact frozen filter outputs under
`examples/spn534`. Multi-million-row all-pair caches and the NCBI assembly
mirror are intentionally omitted; the exact assembly accessions, versions and
source-file SHA256 remain available under `examples/spn534/provenance`.

With HC69 and order distance 100, filter v1 returns PP 4, PT 58, TP 160, P↔T
218, TT 1, total 223 oriented records and 223 unique undirected locus pairs.
The four PP smoke-test edges and `gyrB--aguA` survive; `rsuA_2--glyS` does not.
For `gyrB--aguA`, MIraw=2.1602104671 nats, NMI=0.7982993337,
order=120.2857 genes, driver counts `7:1;16:3;27:3`, and Neff=2.5789473684.

SPN534 is a development dataset, not an unseen validation cohort. The candidate
table is not a list of experimentally validated epistatic interactions.
