"""Staged SHC permutation, BH correction, and ARACNE filtering for PP edges."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from cinch.network import aracne_prune
from cinch.population import informative_clusters, phylo_permutation
from cinch.profiles import load_profile
from cinch.statistics import bh_adjust


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp.parquet")
    frame.to_parquet(temporary, index=False)
    os.replace(temporary, path)


def _atomic_json(payload: dict, path: Path) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _aligned_column(path: Path, samples: np.ndarray, value_column: str, dtype) -> np.ndarray:
    table = pd.read_csv(path, sep=None, engine="python", dtype={"sample_id": str})
    if "sample_id" not in table or value_column not in table:
        raise ValueError(f"{path} must contain sample_id and {value_column}")
    if table.sample_id.duplicated().any():
        raise ValueError(f"duplicate sample_id in {path}")
    indexed = table.set_index("sample_id")
    missing = sorted(set(samples) - set(indexed.index))
    extra = sorted(set(indexed.index) - set(samples))
    if missing or extra:
        raise ValueError(f"sample mismatch in {path}: missing={missing[:5]}, extra={extra[:5]}")
    return indexed.loc[samples, value_column].to_numpy(dtype=dtype)


def _edge_seed(seed: int, edge_id: str) -> int:
    return int.from_bytes(hashlib.sha256(f"{seed}|{edge_id}".encode()).digest()[:8], "little")


def run_advanced_filter(
    association_directory: Path,
    profile_path: Path,
    shc_path: Path,
    output: Path,
    *,
    weights_path: Path | None = None,
    permutations: int = 999,
    minimum_cluster_size: int = 4,
    minimum_informative_clusters: int = 2,
    alpha: float = 0.05,
    seed: int = 20260926,
) -> Path:
    """Test PP edges only under the preserved binary SHC procedure."""

    if permutations < 1 or minimum_cluster_size < 2 or minimum_informative_clusters < 1:
        raise ValueError("invalid permutation or SHC eligibility configuration")
    if not 0 < alpha <= 1:
        raise ValueError("alpha must be in (0, 1]")
    association_directory = Path(association_directory).resolve()
    profile_path, shc_path, output = Path(profile_path).resolve(), Path(shc_path).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    block_paths = sorted((association_directory / "blocks" / "PP").glob("block_*.parquet"))
    if not block_paths:
        raise FileNotFoundError("no PP association blocks found")
    frames = [pd.read_parquet(path) for path in block_paths]
    candidates = pd.concat([frame for frame in frames if len(frame)], ignore_index=True) if any(len(frame) for frame in frames) else pd.DataFrame()
    required = {"edge_id", "i", "j", "locus_A", "locus_B", "channel"}
    if candidates.empty:
        candidates = pd.DataFrame(columns=sorted(required))
    if not required.issubset(candidates.columns):
        raise ValueError(f"PP blocks are missing columns: {sorted(required - set(candidates.columns))}")
    if not candidates.channel.eq("PP").all():
        raise ValueError("advanced SHC filter currently accepts PP edges only")
    profile = load_profile(profile_path)
    clusters = _aligned_column(shc_path, profile.samples, "shc", object)
    if weights_path is None:
        weights = np.ones(len(profile.samples), dtype=float)
    else:
        weights_path = Path(weights_path).resolve()
        weights = _aligned_column(weights_path, profile.samples, "weight", float)
        if np.any(~np.isfinite(weights)) or np.any(weights < 0) or weights.sum() <= 0:
            raise ValueError("weights must be finite, non-negative, and have positive total")
    records = []
    for row in candidates.itertuples(index=False):
        left, right = int(row.i), int(row.j)
        valid = (profile.presence[:, left] >= 0) & (profile.presence[:, right] >= 0)
        x = profile.presence[valid, left].astype(bool)
        y = profile.presence[valid, right].astype(bool)
        local_clusters, local_weights = clusters[valid], weights[valid]
        informative, positive = informative_clusters(x, y, local_clusters, minimum_cluster_size)
        record = row._asdict()
        record.update({
            "shc_callable_N": int(valid.sum()),
            "informative_shc": len(informative),
            "positive_shc": positive,
            "eligible_for_permutation": len(informative) >= minimum_informative_clusters and positive >= minimum_informative_clusters,
            "conditional_mi": np.nan,
            "conditional_mi_null_mean": np.nan,
            "conditional_mi_zscore": np.nan,
            "phylo_perm_p": np.nan,
        })
        if record["eligible_for_permutation"]:
            observed, null_mean, zscore, pvalue = phylo_permutation(
                x, y, local_weights, local_clusters, informative, permutations,
                np.random.default_rng(_edge_seed(seed, str(row.edge_id))),
            )
            record.update({
                "conditional_mi": observed,
                "conditional_mi_null_mean": null_mean,
                "conditional_mi_zscore": zscore,
                "phylo_perm_p": pvalue,
            })
        records.append(record)
    results = pd.DataFrame(records)
    defaults = {
        "eligible_for_permutation": pd.Series(dtype=bool),
        "conditional_mi_zscore": pd.Series(dtype=float),
        "phylo_perm_p": pd.Series(dtype=float),
    }
    for column, empty in defaults.items():
        if column not in results:
            results[column] = empty
    results["phylo_perm_q"] = np.nan
    eligible = results.eligible_for_permutation.astype(bool) if len(results) else pd.Series(dtype=bool)
    if eligible.any():
        results.loc[eligible, "phylo_perm_q"] = bh_adjust(results.loc[eligible, "phylo_perm_p"].to_numpy(float))
    results["passes_shc_bh"] = eligible & (results.phylo_perm_q < alpha) & (results.conditional_mi_zscore > 0)
    significant = results[results.passes_shc_bh].copy()
    if len(significant):
        significant["gene1"], significant["gene2"] = significant.locus_A, significant.locus_B
        significant = aracne_prune(significant, "conditional_mi")
    else:
        significant["aracne_direct"] = pd.Series(dtype=bool)
    direct = significant[significant.aracne_direct].copy() if len(significant) else significant.copy()
    _atomic_parquet(results, output / "SHC_ALL_PP_CANDIDATES.parquet")
    _atomic_parquet(significant, output / "SHC_SIGNIFICANT_ARACNE.parquet")
    _atomic_parquet(direct, output / "SHC_DIRECT_EDGES.parquet")
    provenance = {
        "schema": "CINCH_SHC_FILTER_V1",
        "scientific_scope": "PP binary presence edges only",
        "association_manifest_sha256": _sha256(association_directory / "MANIFEST.json"),
        "profile_sha256": _sha256(profile_path),
        "shc_sha256": _sha256(shc_path),
        "weights_sha256": _sha256(weights_path) if weights_path else None,
        "permutations": permutations,
        "minimum_cluster_size": minimum_cluster_size,
        "minimum_informative_clusters": minimum_informative_clusters,
        "alpha": alpha,
        "seed": seed,
        "edge_seed_rule": "sha256(seed|edge_id)",
        "candidate_edges": len(results),
        "permutation_eligible": int(eligible.sum()),
        "bh_significant": int(results.passes_shc_bh.sum()) if len(results) else 0,
        "aracne_direct": len(direct),
    }
    _atomic_json(provenance, output / "MANIFEST.json")
    return output
