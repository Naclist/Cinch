from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from . import VERSION
from .hierarchy import scan_hc, select_hc_level
from .statistics import CHANNELS, DISPLAY, channel_vectors, neff_from_counts
from .utils import RunLog, hash_rows, sha256, software_versions, write_json, write_tsv, write_yaml
from .wgs import envelope


class SelectionRequiredError(RuntimeError):
    pass


def save_figure(fig, stem: Path):
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(stem.with_suffix(f".{suffix}"), dpi=280 if suffix == "png" else None,
                    bbox_inches="tight")
    plt.close(fig)


def load_pair(path: Path, channel: str) -> pd.DataFrame:
    frame = pd.read_parquet(path)
    renames = {"order_dist": "order_distance", "bp_dist": "bp_distance",
               "N_informative": "informative_N", "HC69_Neff": "Neff",
               "HC69_supported": "HC_supported",
               "HC69_driver_counts_total": "HC_driver_counts_total",
               "global_driver_A_state": "enriched_driver_A_state",
               "global_driver_B_state": "enriched_driver_B_state"}
    frame = frame.rename(columns={k: v for k, v in renames.items() if k in frame.columns})
    if "display_channel" not in frame: frame["display_channel"] = DISPLAY[channel]
    if "orientation" not in frame:
        frame["orientation"] = {"PP":"P_A__P_B","PT":"P_A__T_B","TP":"T_A__P_B","TT":"T_A__T_B"}[channel]
    if "edge_id" not in frame:
        frame["edge_id"] = [f"{channel}:{a}--{b}" for a, b in zip(frame.locus_A, frame.locus_B)]
    return frame


def attach_high_information(frame: pd.DataFrame, bins: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    valid_bins = bins[bins.N.fillna(0).gt(0)].sort_values("log10_upper")
    x = np.log10(out.order_distance.to_numpy(float) + 1.0)
    upper = valid_bins.log10_upper.to_numpy(float)
    ids = np.clip(np.searchsorted(upper, x, side="left"), 0, len(valid_bins) - 1)
    selected = valid_bins.iloc[ids].reset_index(drop=True)
    out["order_bin"] = selected.bin.to_numpy()
    out["local_background_N"] = selected.N.to_numpy()
    for field in ("Q95", "Q99", "Q99_9"):
        out[f"conditional_{field}"] = selected[field].to_numpy(float)
    out["empirical_percentile"] = np.nan
    for bid, idx in out.groupby("order_bin").groups.items():
        values = out.loc[idx, "MIraw"].to_numpy(float)
        order = np.argsort(values, kind="mergesort")
        ranks = np.empty(len(values), float); ranks[order] = (np.arange(len(values)) + 1) / len(values)
        out.loc[idx, "empirical_percentile"] = ranks
    # SPN534 preserves the audited frozen Q99.9 membership byte-for-byte.
    if "high_score_Q99_9" in out:
        out["high_information_pass"] = out.high_score_Q99_9.astype(bool)
    else:
        out["high_information_pass"] = (out.order_distance.notna() & out.MIraw.notna()
                                        & out.MIraw.ge(out.conditional_Q99_9))
    return out


def driver_recurrence(frame: pd.DataFrame, presence: np.ndarray, types: np.ndarray,
                      blocks: np.ndarray, minimum_informative: int) -> pd.DataFrame:
    out = frame.copy()
    supported, totals, distributions, neffs, testable = [], [], [], [], []
    for row in out.itertuples(index=False):
        i, j = int(row.i), int(row.j)
        a = int(getattr(row, "enriched_driver_A_state"))
        b = int(getattr(row, "enriched_driver_B_state"))
        x, y = channel_vectors(presence, types, i, j, row.channel)
        # channel_vectors loses sample indices, so reproduce the legal mask here.
        pa, pb, ta, tb = presence[:, i], presence[:, j], types[:, i], types[:, j]
        if row.channel == "PP": valid=(pa>=0)&(pb>=0); xx,yy=pa,pb
        elif row.channel == "PT": valid=(pa>=0)&(pb==1)&(tb>=0); xx,yy=pa,tb
        elif row.channel == "TP": valid=(pa==1)&(ta>=0)&(pb>=0); xx,yy=ta,pb
        else: valid=(pa==1)&(pb==1)&(ta>=0)&(tb>=0); xx,yy=ta,tb
        counts: dict[int, int] = {}; n_testable = 0
        for block in np.unique(blocks):
            mask = valid & (blocks == block)
            if mask.sum() >= minimum_informative and len(np.unique(xx[mask])) >= 2 and len(np.unique(yy[mask])) >= 2:
                n_testable += 1
            count = int(np.sum(mask & (xx == a) & (yy == b)))
            if count: counts[int(block)] = count
        supported.append(len(counts)); totals.append(sum(counts.values())); testable.append(n_testable)
        distributions.append(";".join(f"{k}:{v}" for k, v in sorted(counts.items())))
        neffs.append(neff_from_counts(counts))
    out["HC_testable"] = testable; out["HC_supported"] = supported
    out["HC_driver_counts_total"] = totals; out["HC_driver_distribution"] = distributions
    out["Neff"] = neffs
    return out


def neff_background(frames: dict[str, pd.DataFrame], n_bins: int = 200):
    backgrounds, edges_by_display = [], {}
    for display in ("PP", "P↔T", "TT"):
        sources = [frames[c] for c in CHANNELS if DISPLAY[c] == display]
        d = pd.concat(sources, ignore_index=True)
        valid = d.order_distance.notna() & d.Neff.notna()
        x = np.log10(d.loc[valid, "order_distance"].to_numpy(float) + 1.0)
        y = d.loc[valid, "Neff"].to_numpy(float)
        effective_bins = min(n_bins, max(4, len(x) // 1000))
        edges = np.linspace(x.min(), np.nextafter(x.max(), np.inf), effective_bins + 1)
        ids = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, n_bins - 1)
        work = pd.DataFrame({"bin": ids, "x": x, "Neff": y})
        grouped = work.groupby("bin", sort=True)
        bg = grouped.agg(N=("Neff","size"), Neff_mean=("Neff","mean"),
                         Neff_median=("Neff","median"), order_log10_mean=("x","mean")).reset_index()
        bg["Neff_Q95"] = grouped.Neff.quantile(.95).to_numpy()
        bg.insert(0, "display_channel", display)
        bg["log10_lower"] = edges[bg.bin]; bg["log10_upper"] = edges[bg.bin + 1]
        backgrounds.append(bg); edges_by_display[display] = edges
    return pd.concat(backgrounds, ignore_index=True), edges_by_display


def attach_neff_filter(frame: pd.DataFrame, background: pd.DataFrame, edges: np.ndarray,
                       minimum_populations: int) -> pd.DataFrame:
    out = frame.copy()
    if len(out) == 0:
        return out
    display = str(out.display_channel.iloc[0])
    bg = background[background.display_channel.eq(display)].set_index("bin")
    x = np.log10(out.order_distance.to_numpy(float) + 1.0)
    ids = np.clip(np.searchsorted(edges, x, side="right") - 1, 0, len(edges) - 2)
    out["order_neff_bin"] = ids
    out["order_neff_background_N"] = bg.reindex(ids).N.to_numpy()
    out["order_neff_background_mean"] = bg.reindex(ids).Neff_mean.to_numpy()
    out["order_neff_background_Q95"] = bg.reindex(ids).Neff_Q95.to_numpy()
    out["order_neff_empirical_percentile"] = np.nan
    for bid, idx in out.groupby("order_neff_bin").groups.items():
        values = out.loc[idx, "Neff"].to_numpy(float)
        order = np.argsort(values, kind="mergesort"); rank = np.empty(len(values), float)
        rank[order] = (np.arange(len(values)) + 1) / len(values)
        out.loc[idx, "order_neff_empirical_percentile"] = rank
    out["Neff_excess"] = out.Neff - out.order_neff_background_mean
    out["Neff_ratio_to_order_mean"] = out.Neff / out.order_neff_background_mean
    out["population_count_pass"] = out.HC_supported.ge(minimum_populations)
    out["Neff_above_order_mean_pass"] = out.Neff.gt(out.order_neff_background_mean)
    out["Neff_Q95_sensitivity_pass"] = out.Neff.gt(out.order_neff_background_Q95)
    out["Neff_filter_pass"] = out.population_count_pass & out.Neff_above_order_mean_pass
    return out


def basic_plots(out: Path, hc_scan: pd.DataFrame, order_bins: pd.DataFrame,
                all_frames: dict[str, pd.DataFrame], final: pd.DataFrame,
                stage_counts: pd.DataFrame, selected_hc: int, order_threshold: float,
                neff_bg: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for channel, group in order_bins.groupby("channel"):
        ax.plot(10 ** ((group.log10_lower + group.log10_upper) / 2) - 1, group.Q99_9, label=channel)
    ax.axvline(order_threshold, color="black", ls="--"); ax.set(xscale="log", xlabel="Order distance (genes)", ylabel="Conditional MI Q99.9")
    ax.legend(frameon=False); save_figure(fig, out / "01_order_scan" / "ORDER_ELBOW")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    axes[0].plot(hc_scan.HC_level, hc_scan.NMI_to_previous); axes[0].axvline(selected_hc, ls="--", color="black")
    axes[0].set(xlabel="HC allele-distance threshold", ylabel="NMI to previous HC")
    axes[1].plot(hc_scan.HC_level, hc_scan.silhouette); axes[1].axvline(selected_hc, ls="--", color="black")
    axes[1].set(xlabel="HC allele-distance threshold", ylabel="Silhouette")
    save_figure(fig, out / "02_hc_scan" / "HC_STABILITY")
    fig, ax = plt.subplots(figsize=(7, 4)); ax.plot(hc_scan.HC_level, hc_scan.n_clusters)
    ax.axvline(selected_hc, ls="--", color="black"); ax.set(xlabel="HC threshold", ylabel="Number of clusters", yscale="log")
    save_figure(fig, out / "02_hc_scan" / "HC_CLUSTER_SIZES")
    relations = [("order_distance","NMI","NMI_VS_ORDER"),("order_distance","Neff","NEFF_VS_ORDER"),("Neff","NMI","NMI_VS_NEFF")]
    for xfield, yfield, name in relations:
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 4), constrained_layout=True)
        for ax, display in zip(axes, ("PP","P↔T","TT")):
            d = pd.concat([all_frames[c] for c in CHANNELS if DISPLAY[c] == display], ignore_index=True)
            ax.scatter(d[xfield], d[yfield], s=2, alpha=.08, rasterized=True)
            chosen = final[final.display_channel.eq(display)]
            ax.scatter(chosen[xfield], chosen[yfield], s=26, facecolor="#F2A900", edgecolor="black", linewidth=.5)
            ax.set(xlabel=xfield, ylabel=yfield, title=f"{display} N={len(d):,}")
        save_figure(fig, out / "08_figures" / name)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    pivot = stage_counts.pivot(index="stage", columns="channel", values="n_edges").fillna(0)
    pivot.plot(kind="bar", ax=ax); ax.set(ylabel="Oriented records"); ax.tick_params(axis="x", rotation=25)
    save_figure(fig, out / "08_figures" / "FILTER_ATTRITION")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for display, group in neff_bg.groupby("display_channel"):
        ax.plot(group.order_log10_mean, group.Neff_mean, label=display)
    ax.set(xlabel="log10(order distance + 1)", ylabel="All-pair mean Neff"); ax.legend(frameon=False)
    save_figure(fig, out / "06_neff" / "NEFF_DIAGNOSTIC")
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for display, group in final.groupby("display_channel"):
        ax.scatter(group.order_distance, group.NMI, s=22 + 8 * group.Neff, alpha=.75, label=display)
    ax.set(xlabel="Order distance (genes)", ylabel="NMI"); ax.legend(frameon=False)
    save_figure(fig, out / "08_figures" / "FINAL_CANDIDATE_LANDSCAPE")


def run_filter(wgs_results: Path, hc: int | None = None, order_threshold: float | None = None,
               minimum_populations: int = 3) -> Path:
    started = time.perf_counter(); root = wgs_results.resolve(); out = root / "filter"
    if out.exists():
        # Results are deterministic; explicit reruns replace only the filter stage.
        import shutil; shutil.rmtree(out)
    for name in ("00_manifest","01_order_scan","02_hc_scan","03_high_information","04_order_filtered",
                 "05_hc_recurrent","06_neff","07_final","08_figures","09_report"):
        (out / name).mkdir(parents=True, exist_ok=True)
    log = RunLog(out / "FILTER.log")
    try:
        manifest_path = root / "00_manifest" / "RUN_MANIFEST.yaml"
        if not manifest_path.exists():
            manifest_path = root / "00_manifest" / "WGS_MANIFEST.yaml"
        manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("schema") != "CINCH_WGS_RESULT_V1": raise ValueError("not a Cinch WGS v1 result directory")
        log.stage("FILTER-01", f"validated WGS manifest; samples={manifest['samples']}; loci={manifest['loci']}")
        frames = {c: load_pair(root / "04_pairs" / f"{c}.parquet", c) for c in CHANNELS}
        order_tables = []
        for channel, frame in frames.items():
            order_tables.append(envelope(frame, channel))
        order_bins = pd.concat(order_tables, ignore_index=True)
        write_tsv(order_bins, out / "01_order_scan" / "ORDER_SCAN.tsv")
        for channel in CHANNELS:
            write_tsv(order_bins[order_bins.channel.eq(channel)], out / "01_order_scan" / f"ORDER_BACKGROUND_{channel}.tsv")
        order_selection = {"status": "USER_DEFINED" if order_threshold is not None else "AMBIGUOUS",
                           "selected": order_threshold, "reason": "no audited deterministic common-threshold selector exists"}
        write_tsv(pd.DataFrame([order_selection]), out / "01_order_scan" / "ORDER_ELBOW_SUMMARY.tsv")
        log.stage("FILTER-02", f"scanned {len(order_bins)} oriented channel/order bins")
        profile = np.load(root / manifest["profile_distance"])
        distance = profile["reference_hiercc_distance"]
        hc_scan, assignments = scan_hc(distance, max_level=min(int(np.max(distance)), profile["profile"].shape[1]))
        auto_hc = select_hc_level(hc_scan)
        selected_hc = int(hc) if hc is not None else auto_hc["selected"]
        hc_status = "USER_DEFINED" if hc is not None else auto_hc["status"]
        write_tsv(hc_scan, out / "02_hc_scan" / "HC_LEVEL_SCAN.tsv")
        ranges = pd.DataFrame(auto_hc["stable_ranges"], columns=["HC_lower","HC_upper"])
        write_tsv(ranges, out / "02_hc_scan" / "HC_STABLE_RANGES.tsv")
        log.stage("FILTER-03", f"ORDER_SELECTION={order_selection['status']} selected={order_threshold}")
        log.stage("FILTER-04", f"HC levels scanned={len(hc_scan)}; auto={auto_hc['status']}; manual={hc}")
        if selected_hc is not None and selected_hc not in assignments:
            raise ValueError(f"HC{selected_hc} outside scanned range")
        if selected_hc is not None:
            hc_assign = pd.DataFrame({"sample_id": profile["samples"].astype(str), "HC_cluster": assignments[selected_hc]})
            write_tsv(hc_assign, out / "02_hc_scan" / "HC_ASSIGNMENTS.tsv")
        if order_threshold is None or selected_hc is None:
            config = {"order_selection": order_selection, "HC_selection": {"status": hc_status, "selected": selected_hc,
                      "auto_diagnostics": auto_hc}, "status": "SELECTION_REQUIRED"}
            write_yaml(config, out / "00_manifest" / "FILTER_CONFIG.yaml")
            raise SelectionRequiredError("automatic selection was ambiguous/unresolved; provide --order-threshold and/or --hc")
        log.stage("FILTER-05", f"selected HC{selected_hc} ({hc_status}); HCx is an allele-distance threshold, not x clusters")
        annotated = {}
        for channel, frame in frames.items():
            annotated[channel] = attach_high_information(frame, order_bins[order_bins.channel.eq(channel)])
        high = pd.concat([x[x.high_information_pass] for x in annotated.values()], ignore_index=True)
        high.to_parquet(out / "03_high_information" / "HIGH_INFORMATION.parquet", index=False)
        write_tsv(high, out / "03_high_information" / "HIGH_INFORMATION.tsv")
        log.stage("FILTER-06", "high-information pass: " + ", ".join(f"{c}={int(annotated[c].high_information_pass.sum())}" for c in CHANNELS))
        order_filtered = high[high.order_distance.ge(float(order_threshold))].copy()
        write_tsv(order_filtered, out / "04_order_filtered" / "ORDER_FILTERED.tsv")
        log.stage("FILTER-07", f"order >= {order_threshold:g}: {len(high):,} -> {len(order_filtered):,}")
        precomputed_hc = manifest.get("precomputed_HC_level")
        if precomputed_hc is not None and int(precomputed_hc) == selected_hc and all("Neff" in x for x in annotated.values()):
            recurrent_all = annotated
            distribution_path = manifest.get("precomputed_HC_driver_distributions")
            distributions = None
            if distribution_path and (root / distribution_path).exists():
                distributions = pd.read_csv(root / distribution_path, sep="\t")
            for channel, frame in recurrent_all.items():
                if "HC_driver_distribution" not in frame and "HC69_global_driver_counts" in frame:
                    frame["HC_driver_distribution"] = frame.HC69_global_driver_counts
                if "HC_driver_distribution" not in frame and distributions is not None:
                    take = distributions[distributions.channel.eq(channel)][["locus_A","locus_B","HC_driver_distribution"]]
                    recurrent_all[channel] = frame.merge(take, on=["locus_A","locus_B"], how="left")
                    frame = recurrent_all[channel]
                if "HC_driver_distribution" not in frame:
                    frame["HC_driver_distribution"] = ""
                else:
                    frame["HC_driver_distribution"] = frame.HC_driver_distribution.fillna("")
                if "HC_testable" not in frame: frame["HC_testable"] = np.nan
        else:
            state = np.load(root / manifest["states"])
            presence, types = state["presence"], state["type_state"]
            blocks = assignments[selected_hc]
            recurrent_all = {c: driver_recurrence(frame, presence, types, blocks,
                                                   int(yaml.safe_load((root / "00_manifest" / "WGS_CONFIG.yaml").read_text())["eligibility"]["minimum_informative_N"]))
                             for c, frame in annotated.items()}
        recurrent_candidates = pd.concat([x[x.high_information_pass & x.order_distance.ge(float(order_threshold)) & x.HC_supported.ge(minimum_populations)]
                                           for x in recurrent_all.values()], ignore_index=True)
        write_tsv(recurrent_candidates, out / "05_hc_recurrent" / "HC_RECURRENT.tsv")
        write_tsv(recurrent_candidates[[c for c in recurrent_candidates.columns if c in ("edge_id","channel","locus_A","locus_B","HC_supported","HC_driver_counts_total","HC_driver_distribution","Neff")]],
                  out / "05_hc_recurrent" / "HC_DRIVER_DISTRIBUTIONS.tsv")
        log.stage("FILTER-08", f"same-driver recurrence in >= {minimum_populations} HC blocks: {len(recurrent_candidates):,}")
        neff_bg, neff_edges = neff_background(recurrent_all)
        write_tsv(neff_bg, out / "06_neff" / "NEFF_ORDER_BACKGROUND.tsv")
        filtered_frames = {c: attach_neff_filter(frame, neff_bg, neff_edges[DISPLAY[c]], minimum_populations)
                           for c, frame in recurrent_all.items()}
        final = pd.concat([x[x.high_information_pass & x.order_distance.ge(float(order_threshold)) & x.Neff_filter_pass]
                           for x in filtered_frames.values()], ignore_index=True)
        final["HC_level"] = selected_hc; final["order_threshold"] = float(order_threshold)
        final["order_filter_pass"] = True; final["HC_recurrence_pass"] = True; final["final_pass"] = True
        write_tsv(final, out / "06_neff" / "NEFF_FILTERED.tsv")
        log.stage("FILTER-09", f"Neff above channel/order-conditioned all-pair mean: final={len(final):,}")
        metadata_path = root / "02_loci" / "LOCUS_METADATA.tsv"
        meta = pd.read_csv(metadata_path, sep="\t").set_index("locus_id") if metadata_path.exists() else pd.DataFrame()
        if not meta.empty and "annotation" in meta:
            final["annotation_A"] = final.locus_A.map(meta.annotation)
            final["annotation_B"] = final.locus_B.map(meta.annotation)
        else:
            final["annotation_A"] = final.locus_A; final["annotation_B"] = final.locus_B
        required = ["edge_id","locus_A","locus_B","channel","display_channel","orientation","MIraw","NMI",
                    "informative_N","order_distance","order_threshold","local_background_N","conditional_Q95","conditional_Q99",
                    "conditional_Q99_9","empirical_percentile","HC_level","HC_testable","HC_supported",
                    "HC_driver_counts_total","HC_driver_distribution","Neff","order_neff_background_N",
                    "order_neff_background_mean","order_neff_background_Q95","order_neff_empirical_percentile","Neff_excess",
                    "Neff_ratio_to_order_mean","enriched_driver_A_state","enriched_driver_B_state","depleted_driver_A_state",
                    "depleted_driver_B_state","annotation_A","annotation_B","high_information_pass","order_filter_pass",
                    "HC_recurrence_pass","population_count_pass","Neff_above_order_mean_pass","Neff_filter_pass","final_pass"]
        for field in required:
            if field not in final: final[field] = np.nan
        final = final[required]
        write_tsv(final, out / "07_final" / "FINAL_CANDIDATES.tsv")
        write_tsv(final, out / "07_final" / "NETWORK_READY_EDGES.tsv")
        nodes = pd.DataFrame({"locus_id": pd.unique(final[["locus_A","locus_B"]].to_numpy().ravel())})
        if len(nodes): nodes["annotation"] = nodes.locus_id.map(meta.annotation) if not meta.empty and "annotation" in meta else nodes.locus_id
        write_tsv(nodes, out / "07_final" / "NETWORK_READY_NODES.tsv")
        log.stage("FILTER-10", f"network-ready candidate dependencies={len(final):,}; ARACNE/network not run")
        stages = []
        for c in CHANNELS:
            a = filtered_frames[c]
            for name, mask in (("eligible", np.ones(len(a), bool)), ("high_information", a.high_information_pass),
                               ("order_filtered", a.high_information_pass & a.order_distance.ge(float(order_threshold))),
                               ("HC_recurrent", a.high_information_pass & a.order_distance.ge(float(order_threshold)) & a.HC_supported.ge(minimum_populations)),
                               ("Neff_filtered", a.high_information_pass & a.order_distance.ge(float(order_threshold)) & a.Neff_filter_pass)):
                stages.append({"channel": c, "stage": name, "n_edges": int(np.sum(mask))})
        stage_counts = pd.DataFrame(stages); write_tsv(stage_counts, out / "09_report" / "FILTER_COUNTS.tsv")
        basic_plots(out, hc_scan, order_bins, filtered_frames, final, stage_counts, selected_hc,
                    float(order_threshold), neff_bg)
        config = {"workflow": "CINCH_FILTER_V1", "version": VERSION,
                  "wgs_results": str(root), "order_selection": {"status": "USER_DEFINED", "selected": float(order_threshold)},
                  "HC_selection": {"status": hc_status, "selected": selected_hc, "definition": "single-linkage allele-distance threshold"},
                  "high_information": {"quantile": .999, "oriented_channel_specific": True},
                  "HC_recurrence_minimum_populations": minimum_populations,
                  "Neff": "1/sum((n_h/sum(n_h))^2)", "Neff_primary_rule": "Neff > display-channel/order-bin all-eligible-pair mean",
                  "ARACNE": False, "network_construction": False}
        write_yaml(config, out / "00_manifest" / "FILTER_CONFIG.yaml")
        (out / "00_manifest" / "COMMAND.txt").write_text(
            f"cinch filter --wgs_results {root} --hc {selected_hc} --order-threshold {order_threshold:g}\n", encoding="utf-8")
        input_files = [manifest_path, *(root / "04_pairs").glob("*.parquet"), root / manifest["profile_distance"]]
        write_tsv(hash_rows(input_files), out / "00_manifest" / "FILTER_INPUT_SHA256.tsv")
        filter_manifest = {"schema":"CINCH_FILTER_RESULT_V1","version":VERSION,"source_manifest_sha256":sha256(manifest_path),
                           "config_sha256":sha256(out / "00_manifest" / "FILTER_CONFIG.yaml"),
                           "final_candidates":len(final),"software":software_versions(),"ARACNE":False,"network_built":False}
        write_yaml(filter_manifest, out / "00_manifest" / "FILTER_MANIFEST.yaml")
        validation = pd.DataFrame([
            {"check":"all_final_high_information", "pass":bool(final.high_information_pass.all()), "detail":len(final)},
            {"check":"all_final_order_filtered", "pass":bool(final.order_distance.ge(order_threshold).all()), "detail":float(final.order_distance.min()) if len(final) else np.nan},
            {"check":"all_final_HC_recurrent", "pass":bool(final.HC_supported.ge(minimum_populations).all()), "detail":int(final.HC_supported.min()) if len(final) else np.nan},
            {"check":"all_final_Neff_above_mean", "pass":bool(final.Neff.gt(final.order_neff_background_mean).all()), "detail":float(final.Neff_excess.min()) if len(final) else np.nan},
            {"check":"ARACNE_not_run", "pass":True, "detail":"filter stops at network-ready edges"},
        ])
        expectations = manifest.get("validation_expectations")
        if expectations:
            extra = []
            observed_counts = final.channel.value_counts().to_dict()
            for channel, expected in expectations.get("oriented_counts", {}).items():
                extra.append({"check":f"expected_{channel}_count", "pass":int(observed_counts.get(channel,0)) == int(expected),
                              "detail":int(observed_counts.get(channel,0))})
            unique_observed = len(set(tuple(sorted((a,b))) for a,b in final[["locus_A","locus_B"]].itertuples(index=False,name=None)))
            extra.extend([{"check":"expected_total", "pass":len(final)==int(expectations["total"]), "detail":len(final)},
                          {"check":"expected_unique_undirected", "pass":unique_observed==int(expectations["unique_undirected"]), "detail":unique_observed}])
            def pair_rows(pair):
                a,b=pair
                return final[((final.locus_A==a)&(final.locus_B==b))|((final.locus_A==b)&(final.locus_B==a))]
            for pair in expectations.get("retained_pairs", []):
                extra.append({"check":f"retained_{pair[0]}--{pair[1]}","pass":len(pair_rows(pair))==1,"detail":len(pair_rows(pair))})
            for pair in expectations.get("rejected_pairs", []):
                extra.append({"check":f"rejected_{pair[0]}--{pair[1]}","pass":len(pair_rows(pair))==0,"detail":len(pair_rows(pair))})
            metric_expect = expectations.get("gyrB_aguA")
            if metric_expect:
                row = pair_rows(["gyrB","aguA"]).iloc[0]
                for field in ("MIraw","NMI","order_distance","Neff"):
                    extra.append({"check":f"gyrB_aguA_{field}","pass":bool(np.isclose(float(row[field]),float(metric_expect[field]),atol=1e-6,rtol=1e-9)),"detail":row[field]})
                for field in ("HC_supported","HC_driver_distribution"):
                    extra.append({"check":f"gyrB_aguA_{field}","pass":str(row[field])==str(metric_expect[field]),"detail":row[field]})
            validation = pd.concat([validation,pd.DataFrame(extra)],ignore_index=True)
        write_tsv(validation, out / "09_report" / "VALIDATION.tsv")
        if not validation["pass"].all():
            raise RuntimeError("filter validation failed:\n" + validation.loc[~validation["pass"]].to_string(index=False))
        display_counts = final.display_channel.value_counts().to_dict()
        oriented = final.channel.value_counts().to_dict()
        unique_pairs = len(set(tuple(sorted((a,b))) for a,b in final[["locus_A","locus_B"]].itertuples(index=False, name=None)))
        report = f"# Cinch Filter v1 summary\n\nHC level: {selected_hc} (single-linkage allele-distance threshold)  \nOrder threshold: {order_threshold:g} genes  \nFinal oriented records: {len(final)}  \nPP: {display_counts.get('PP',0)}; PT: {oriented.get('PT',0)}; TP: {oriented.get('TP',0)}; P↔T: {display_counts.get('P↔T',0)}; TT: {display_counts.get('TT',0)}  \nUnique undirected locus pairs: {unique_pairs}  \nRuntime seconds: {time.perf_counter()-started:.3f}\n\nThese are cross-population-recurrent, order-distant, high-information candidate dependencies. They are not causal interactions. ARACNE and network construction were not run.\n"
        (out / "09_report" / "FILTER_SUMMARY.md").write_text(report, encoding="utf-8")
        return out
    finally:
        log.close()
