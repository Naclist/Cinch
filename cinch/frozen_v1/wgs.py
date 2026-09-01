from __future__ import annotations

import shutil
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import VERSION
from .hierarchy import hiercc_distances
from .mapping import map_genomes, pair_order_distance
from .statistics import CHANNELS, score_all_pairs
from .utils import RunLog, hash_rows, sha256, software_versions, write_json, write_tsv, write_yaml


def save_figure(fig, stem: Path):
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(stem.with_suffix(f".{suffix}"), dpi=280 if suffix == "png" else None,
                    bbox_inches="tight")
    plt.close(fig)


def envelope(frame: pd.DataFrame, channel: str, target_pairs_per_bin: int = 1000) -> pd.DataFrame:
    valid = frame.MIraw.notna() & frame.order_distance.notna() & frame.order_distance.ge(0)
    x = frame.loc[valid, "order_distance"].to_numpy(float)
    y = frame.loc[valid, "MIraw"].to_numpy(float)
    if len(x) == 0:
        return pd.DataFrame()
    z = np.log10(x + 1.0)
    n_bins = min(200, max(4, int(len(x) // target_pairs_per_bin)))
    edges = np.linspace(float(z.min()), float(z.max()) if z.max() > z.min() else float(z.min() + 1), n_bins + 1)
    ids = np.clip(np.searchsorted(edges, z, side="right") - 1, 0, n_bins - 1)
    work = pd.DataFrame({"bin": ids, "distance": x, "MIraw": y})
    grouped = work.groupby("bin", sort=True)
    central = grouped.agg(N=("MIraw", "size"), distance_median=("distance", "median"),
                          median_MIraw=("MIraw", "median"), mean_MIraw=("MIraw", "mean"))
    quantiles = grouped.MIraw.quantile([.25, .75, .95, .99, .999]).unstack()
    quantiles.columns = ["Q25", "Q75", "Q95", "Q99", "Q99_9"]
    base = pd.DataFrame({"channel": channel, "bin": np.arange(n_bins),
                         "log10_lower": edges[:-1], "log10_upper": edges[1:]})
    return base.merge(central.join(quantiles), left_on="bin", right_index=True, how="left")


def validate_phenotype(path: Path, samples: np.ndarray, out: Path) -> dict:
    frame = pd.read_csv(path, sep=None, engine="python", dtype=str)
    if frame.empty or frame.shape[1] < 2:
        raise ValueError("phenotype table must contain sample ID plus at least one phenotype")
    sample_column = frame.columns[0]
    duplicate = frame[sample_column].duplicated().sum()
    known = frame[sample_column].isin(samples)
    validation = pd.DataFrame([{"check": "unique_sample_ids", "pass": duplicate == 0, "detail": int(duplicate)},
                               {"check": "all_samples_known", "pass": bool(known.all()),
                                "detail": int((~known).sum())}])
    write_tsv(validation, out / "08_phenotype" / "PHENOTYPE_VALIDATION.tsv")
    if not validation["pass"].all():
        raise ValueError("phenotype validation failed")
    shutil.copy2(path, out / "08_phenotype" / "PHENOTYPE_FROZEN.tsv")
    return {"supplied": True, "sha256": sha256(path), "sample_column": str(sample_column),
            "phenotype_columns": list(map(str, frame.columns[1:]))}


def run_wgs(reference: Path, genomes: list[Path], prefix: str, threads: int,
            output: Path | None = None, phenotype: Path | None = None,
            minimum_identity: float = .90, minimum_informative: int = 20,
            minimum_state_count: int = 3, minimum_same_contig_observations: int = 5) -> Path:
    started = time.perf_counter()
    reference, genomes = reference.resolve(), [x.resolve() for x in genomes]
    out = (output or Path(f"{prefix}.cinch")).resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"output directory is not empty: {out}")
    for name in ("00_manifest", "01_mapping", "02_profiles", "03_state_qc", "02_loci", "03_states", "04_pairs", "05_drivers",
                 "06_order", "07_population_basis", "08_phenotype", "09_figures", "10_report"):
        (out / name).mkdir(parents=True, exist_ok=True)
    log = RunLog(out / "RUN.log")
    try:
        log.stage("WGS-01", f"validated reference={reference.name}; genomes={len(genomes)}; threads={threads}")
        for path in [reference, *genomes]:
            if not path.is_file(): raise FileNotFoundError(path)
        config = {"workflow": "CINCH_WGS_V1", "version": VERSION, "prefix": prefix,
                  "mapping": {"method": "internal full-length seed-and-extend nucleotide mapper",
                              "minimum_identity": minimum_identity, "complete_CDS_required_for_type": True},
                  "eligibility": {"minimum_informative_N": minimum_informative,
                                  "minimum_marginal_state_count": minimum_state_count,
                                  "rare_state_policy": "pool per locus before scoring"},
                  "order_distance": {"aggregation": "arithmetic mean of within-sample/contig minimum order differences",
                                     "minimum_same_contig_observations": minimum_same_contig_observations},
                  "threads_requested": threads}
        write_yaml(config, out / "00_manifest" / "WGS_CONFIG.yaml")
        log.stage("WGS-02", "indexed reference CDS records")
        mapped = map_genomes(reference, genomes, minimum_identity)
        samples, loci, presence, types, trace, coordinates, annotations, genome_qc = mapped
        log.stage("WGS-03", f"mapped {len(samples)} genomes to {len(loci)} loci")
        write_tsv(genome_qc, out / "01_mapping" / "GENOME_MAPPING_QC.tsv")
        write_tsv(genome_qc, out / "01_mapping" / "MAPPING_QC.tsv")
        trace.to_parquet(out / "01_mapping" / "CDS_TYPE_TRACE.parquet", index=False)
        trace.to_parquet(out / "01_mapping" / "LOCUS_MAPPING.parquet", index=False)
        write_tsv(annotations, out / "02_loci" / "LOCUS_METADATA.tsv")
        write_tsv(annotations, out / "02_profiles" / "LOCUS_METADATA.tsv")
        log.stage("WGS-04", f"presence calls={(presence == 1).sum():,}; type-callable={(types >= 0).sum():,}; missing-type={((presence == 1) & (types < 0)).sum():,}")
        np.savez_compressed(out / "03_states" / "UNIFIED_STATES.npz", samples=samples, loci=loci,
                            presence=presence, type_state=types)
        np.savez_compressed(out / "02_profiles" / "UNIFIED_STATES.npz", samples=samples, loci=loci,
                            presence=presence, type_state=types)
        presence_frame = pd.DataFrame(presence, index=samples, columns=loci).rename_axis("sample_id").reset_index()
        type_frame = pd.DataFrame(types, index=samples, columns=loci).rename_axis("sample_id").reset_index()
        presence_frame.to_parquet(out / "03_states" / "PRESENCE_MATRIX.parquet", index=False)
        type_frame.to_parquet(out / "03_states" / "TYPE_MATRIX.parquet", index=False)
        presence_frame.to_parquet(out / "02_profiles" / "PRESENCE.parquet", index=False)
        type_frame.to_parquet(out / "02_profiles" / "TYPES.parquet", index=False)
        type_frame.to_parquet(out / "02_profiles" / "WGMLST_PROFILE.parquet", index=False)
        coordinates.to_parquet(out / "03_states" / "COORDINATES.parquet", index=False)
        sample_qc = pd.DataFrame({"sample_id": samples, "presence_callable_fraction": (presence >= 0).mean(1),
                                  "presence_fraction": (presence == 1).mean(1), "type_callability_fraction": (types >= 0).mean(1)})
        locus_qc = pd.DataFrame({"locus_id": loci, "presence_callable_fraction": (presence >= 0).mean(0),
                                 "prevalence": (presence == 1).mean(0), "type_callability_fraction": (types >= 0).mean(0),
                                 "n_nominal_types": [len(np.unique(types[types[:,j]>=0,j])) for j in range(len(loci))]})
        write_tsv(sample_qc, out / "03_state_qc" / "SAMPLE_QC.tsv")
        write_tsv(locus_qc, out / "03_state_qc" / "LOCUS_QC.tsv")
        write_tsv(sample_qc[["sample_id","presence_callable_fraction","type_callability_fraction"]], out / "03_state_qc" / "MISSINGNESS.tsv")
        write_tsv(locus_qc[["locus_id","type_callability_fraction","n_nominal_types"]], out / "03_state_qc" / "TYPE_CALLABILITY.tsv")
        log.stage("WGS-05", f"reconstructed {len(coordinates):,} callable locus coordinates")
        explicitly_background = annotations.annotation.astype(str).str.contains("[BACKGROUND]", regex=False).to_numpy()
        base_profile_mask = explicitly_background if explicitly_background.any() else np.ones(len(loci), bool)
        profile_mask = base_profile_mask & ((types >= 0).mean(0) >= .80) & (np.array([len(np.unique(types[types[:, j] >= 0, j])) for j in range(len(loci))]) >= 2)
        if profile_mask.sum() < 2:
            profile_mask = (types >= 0).sum(0) > 0
        profile = types[:, profile_mask]
        distance = hiercc_distances(profile)
        np.savez_compressed(out / "07_population_basis" / "PROFILE_DISTANCE.npz", samples=samples,
                            loci=loci[profile_mask], profile=profile, **distance)
        write_tsv(pd.DataFrame({"locus_id": loci, "profile_included": profile_mask,
                                "type_callability": (types >= 0).mean(0)}),
                  out / "07_population_basis" / "HC_SCAN_INPUT.tsv")
        log.stage("WGS-06", f"built locus profile with {profile_mask.sum()} loci; exact-profile unique N={len(np.unique(profile, axis=0))}")
        distances = pair_order_distance(coordinates, len(loci), minimum_same_contig_observations)
        distances.insert(2, "locus_A", loci[distances.i]); distances.insert(3, "locus_B", loci[distances.j])
        distances.to_parquet(out / "06_order" / "ORDER_DISTANCE.parquet", index=False)
        log.stage("WGS-07", "evaluating PP/PT/TP/TT channel eligibility")
        pair_tables, rare_audit = score_all_pairs(presence, types, loci, distances,
                                                  minimum_informative, minimum_state_count)
        write_tsv(rare_audit, out / "02_loci" / "RARE_STATE_POOLING.tsv")
        counts, envelopes, drivers = [], [], []
        for channel in CHANNELS:
            table = pair_tables[channel]
            table.to_parquet(out / "04_pairs" / f"{channel}.parquet", index=False)
            counts.append({"channel": channel, "eligible_pairs": len(table)})
            envelopes.append(envelope(table, channel))
            if len(table):
                drivers.append(table[[c for c in table.columns if "driver" in c or c in ("edge_id", "locus_A", "locus_B", "channel")]])
            log.stage("WGS-08", f"{channel}: eligible={len(table):,}; MI range={table.MIraw.min() if len(table) else float('nan'):.6g}..{table.MIraw.max() if len(table) else float('nan'):.6g}")
        pd.concat(drivers, ignore_index=True).to_parquet(out / "05_drivers" / "STATE_DRIVERS.parquet", index=False)
        background = pd.concat(envelopes, ignore_index=True)
        write_tsv(background, out / "06_order" / "ORDER_BACKGROUND.tsv")
        log.stage("WGS-09", f"saved enriched/depleted drivers for {sum(len(x) for x in pair_tables.values()):,} eligible records")
        log.stage("WGS-10", f"order distances finite={distances.order_distance.notna().sum():,}/{len(distances):,}; different contigs remain NA")
        log.stage("WGS-11", f"generated {len(background):,} channel/order background bins")
        phenotype_info = {"supplied": False}
        if phenotype:
            phenotype_info = validate_phenotype(phenotype.resolve(), samples, out)
        input_hashes = hash_rows([reference, *genomes] + ([phenotype.resolve()] if phenotype else []))
        write_tsv(input_hashes, out / "00_manifest" / "INPUT_SHA256.tsv")
        command = f"cinch wgs -r {reference} -p {prefix} -t {threads} " + " ".join(map(str, genomes))
        (out / "00_manifest" / "COMMAND.txt").write_text(command + "\n", encoding="utf-8")
        versions = software_versions()
        write_tsv(pd.DataFrame([{"software": key, "version": value} for key, value in versions.items()]),
                  out / "00_manifest" / "SOFTWARE_VERSIONS.tsv")
        manifest = {"schema": "CINCH_WGS_RESULT_V1", "version": VERSION, "prefix": prefix,
                    "samples": len(samples), "loci": len(loci), "pair_storage": "04_pairs/{channel}.parquet",
                    "profile_distance": "07_population_basis/PROFILE_DISTANCE.npz",
                    "states": "02_profiles/UNIFIED_STATES.npz", "coordinates": "03_states/COORDINATES.parquet",
                    "order_background": "06_order/ORDER_BACKGROUND.tsv", "phenotype": phenotype_info,
                    "config_sha256": sha256(out / "00_manifest" / "WGS_CONFIG.yaml"),
                    "software": versions, "scientific_boundaries": {"HC_filter": False,
                    "Neff_filter": False, "ARACNE": False, "population_correction": False}}
        write_yaml(manifest, out / "00_manifest" / "WGS_MANIFEST.yaml")
        write_yaml(manifest, out / "00_manifest" / "RUN_MANIFEST.yaml")
        write_tsv(pd.DataFrame(counts), out / "10_report" / "WGS_COUNTS.tsv")
        validation = pd.DataFrame([
            {"check": "missing_not_absence", "pass": not np.any((presence == 0) & (types >= 0)), "detail": "P=0 always T=NA"},
            {"check": "typed_only_when_present", "pass": not np.any((types >= 0) & (presence != 1)), "detail": "T callable implies P=1"},
            {"check": "all_critical_inputs_hashed", "pass": len(input_hashes) == len(genomes) + 1 + int(phenotype is not None), "detail": len(input_hashes)},
        ])
        write_tsv(validation, out / "10_report" / "VALIDATION.tsv")
        summary = f"# Cinch WGS v1 summary\n\nSamples: {len(samples)}  \nLoci: {len(loci)}  \nEligible oriented records: {sum(x['eligible_pairs'] for x in counts)}  \nRuntime seconds: {time.perf_counter()-started:.3f}\n\nNo HC, Neff, ARACNE, permutation or population correction was applied.\n"
        (out / "10_report" / "WGS_SUMMARY.md").write_text(summary, encoding="utf-8")
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), constrained_layout=True)
        axes[0].hist((presence == 1).mean(0), bins=15); axes[0].set(xlabel="Locus prevalence", ylabel="Loci")
        axes[1].hist((types >= 0).mean(0), bins=15); axes[1].set(xlabel="Type callability", ylabel="Loci")
        axes[2].bar([x["channel"] for x in counts], [x["eligible_pairs"] for x in counts]); axes[2].set(ylabel="Eligible pairs")
        save_figure(fig, out / "09_figures" / "WGS_QC_STATE_AND_CHANNEL_COUNTS")
        log.stage("WGS-12", f"manifest/config/hashes/QC complete; runtime={time.perf_counter()-started:.2f}s")
        return out
    finally:
        log.close()
