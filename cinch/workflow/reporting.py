"""Report-only stage for SHC/BH/ARACNE outputs."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def run_report(filter_directory: Path, output: Path) -> Path:
    filter_directory, output = Path(filter_directory).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((filter_directory / "MANIFEST.json").read_text(encoding="utf-8"))
    all_edges = pd.read_parquet(filter_directory / "SHC_ALL_PP_CANDIDATES.parquet")
    significant = pd.read_parquet(filter_directory / "SHC_SIGNIFICANT_ARACNE.parquet")
    direct = pd.read_parquet(filter_directory / "SHC_DIRECT_EDGES.parquet")
    summary = pd.DataFrame([{
        "candidate_edges": len(all_edges),
        "permutation_eligible": int(all_edges.eligible_for_permutation.sum()) if len(all_edges) else 0,
        "bh_significant": int(all_edges.passes_shc_bh.sum()) if len(all_edges) else 0,
        "aracne_direct": len(direct),
        "aracne_indirect": len(significant) - len(direct),
    }])
    summary.to_csv(output / "SUMMARY.tsv", sep="\t", index=False)
    report = (
        "# CINCH staged SHC report\n\n"
        f"Scientific scope: {manifest['scientific_scope']}  \n"
        f"Candidate PP edges: {len(all_edges)}  \n"
        f"Permutation eligible: {summary.at[0, 'permutation_eligible']}  \n"
        f"BH-significant: {summary.at[0, 'bh_significant']}  \n"
        f"ARACNE direct: {len(direct)}  \n\n"
        "SHC permutation p/q-values are distinct from frozen HC recurrence and Neff filtering. "
        "ARACNE is graph pruning, not a significance test.\n"
    )
    (output / "REPORT.md").write_text(report, encoding="utf-8")
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    eligible = all_edges[all_edges.eligible_for_permutation] if len(all_edges) else all_edges
    if len(eligible):
        axes[0].hist(eligible.phylo_perm_p, bins=min(20, max(5, len(eligible))), color="#5477C4")
        valid = eligible.order_distance.notna() if "order_distance" in eligible else pd.Series(False, index=eligible.index)
        axes[1].scatter(eligible.loc[valid, "order_distance"], eligible.loc[valid, "conditional_mi"], s=14, alpha=.7)
    axes[0].set(xlabel="SHC permutation p", ylabel="PP edges")
    axes[1].set(xlabel="Mean gene-order distance", ylabel="Conditional MI (nats)")
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(output / f"SHC_FILTER_DIAGNOSTICS.{suffix}", dpi=240 if suffix == "png" else None, bbox_inches="tight")
    plt.close(fig)
    return output
