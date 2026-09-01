# Outputs

## WGS directory

- `00_manifest/RUN_MANIFEST.yaml`, `WGS_CONFIG.yaml`, `COMMAND.txt`, `INPUT_SHA256.tsv`, `SOFTWARE_VERSIONS.tsv`: API, parameters and provenance.
- `01_mapping/LOCUS_MAPPING.parquet`, `MAPPING_QC.tsv`: source-level call trace and assembly summaries.
- `02_profiles/PRESENCE.parquet`, `TYPES.parquet`, `WGMLST_PROFILE.parquet`, `LOCUS_METADATA.tsv`, `UNIFIED_STATES.npz`: biological profiles and state object.
- `03_state_qc/SAMPLE_QC.tsv`, `LOCUS_QC.tsv`, `MISSINGNESS.tsv`, `TYPE_CALLABILITY.tsv`: state quality audit.
- `04_pairs/PP.parquet`, `PT.parquet`, `TP.parquet`, `TT.parquet`: all eligible oriented pair records.
- `05_drivers/STATE_DRIVERS.parquet`: enriched and depleted contingency cells.
- `06_order/ORDER_DISTANCE.parquet`, `ORDER_BACKGROUND.tsv`: pair distance and raw-MI envelopes.
- `07_population_basis/PROFILE_DISTANCE.npz`, `HC_SCAN_INPUT.tsv`: independent nominal allele profile and audited distances.
- `08_phenotype/*`: optional validation and frozen copy.
- `09_figures/WGS_QC_STATE_AND_CHANNEL_COUNTS.*`: WGS QC.
- `10_report/WGS_SUMMARY.md`, `WGS_COUNTS.tsv`, `VALIDATION.tsv`: summary and gates.
- `RUN.log`: complete console transcript.

## Filter directory

- `00_manifest/FILTER_MANIFEST.yaml`, `FILTER_CONFIG.yaml`, `COMMAND.txt`, `FILTER_INPUT_SHA256.tsv`.
- `01_order_scan/ORDER_SCAN.tsv`, per-channel backgrounds, elbow summary and plot.
- `02_hc_scan/HC_LEVEL_SCAN.tsv`, stable ranges, assignments and HC figures.
- `03_high_information/HIGH_INFORMATION.tsv/.parquet`.
- `04_order_filtered/ORDER_FILTERED.tsv`.
- `05_hc_recurrent/HC_RECURRENT.tsv`, `HC_DRIVER_DISTRIBUTIONS.tsv`.
- `06_neff/NEFF_ORDER_BACKGROUND.tsv`, `NEFF_FILTERED.tsv`, diagnostic plot.
- `07_final/FINAL_CANDIDATES.tsv`, `NETWORK_READY_EDGES.tsv`, `NETWORK_READY_NODES.tsv`.
- `08_figures/FILTER_ATTRITION`, `NMI_VS_ORDER`, `NEFF_VS_ORDER`, `NMI_VS_NEFF`, `FINAL_CANDIDATE_LANDSCAPE` in PNG/PDF/SVG.
- `09_report/FILTER_SUMMARY.md`, `FILTER_COUNTS.tsv`, `VALIDATION.tsv`.
- `FILTER.log`: complete console transcript.
