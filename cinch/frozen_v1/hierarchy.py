from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score


def canonical_labels(labels: np.ndarray) -> np.ndarray:
    mapping: dict[int, int] = {}
    out = np.empty(len(labels), np.int32)
    for i, value in enumerate(labels):
        if int(value) not in mapping:
            mapping[int(value)] = len(mapping) + 1
        out[i] = mapping[int(value)]
    return out


def hiercc_distances(profile: np.ndarray) -> dict[str, np.ndarray]:
    """Frozen missingness-adjusted integer allele distance from HC4 development."""
    n, n_loci = profile.shape
    x = np.where(profile >= 0, profile + 1, 0).astype(np.int32)
    dual = np.zeros((n, n), np.int32)
    p_dist = np.zeros((n, n), np.int32)
    raw = np.zeros((n, n), np.int32)
    overlap = np.zeros((n, n), np.int32)
    callable_n = (x > 0).sum(1)
    for i in range(n):
        overlap[i, i] = callable_n[i]
        for j in range(i):
            valid = (x[i] > 0) & (x[j] > 0)
            joint = int(valid.sum())
            mismatches = int(((x[i] != x[j]) & valid).sum())
            ad, al = mismatches + 0.0001, joint + 0.0001
            lower = max(callable_n[i], callable_n[j]) - 0.03 * n_loci
            if lower > al:
                ad += lower - al; al = lower
            d = int(ad / al * n_loci + 0.5)
            p = int(-np.log(1.0 - (mismatches + 0.5) / (joint + 1.0)) * n_loci * 100.0 + 0.5) if joint else 0
            dual[i, j] = dual[j, i] = d
            p_dist[i, j] = p_dist[j, i] = p
            raw[i, j] = raw[j, i] = mismatches
            overlap[i, j] = overlap[j, i] = joint
    return {"reference_hiercc_distance": dual, "reference_hcceval_distance": p_dist,
            "raw_joint_callable_mismatches": raw, "callable_overlap": overlap}


def scan_hc(distance: np.ndarray, max_level: int | None = None) -> tuple[pd.DataFrame, dict[int, np.ndarray]]:
    if len(distance) < 2:
        raise ValueError("HC scan needs at least two samples")
    hierarchy = linkage(squareform(distance, checks=False), method="single")
    max_level = int(np.nanmax(distance)) if max_level is None else int(max_level)
    assignments: dict[int, np.ndarray] = {}
    rows, previous = [], None
    for level in range(max_level + 1):
        labels = canonical_labels(fcluster(hierarchy, t=level, criterion="distance"))
        assignments[level] = labels
        sizes = np.bincount(labels)[1:]
        n_clusters = len(sizes)
        sil = float(silhouette_score(distance.astype(float), labels, metric="precomputed")) if 2 <= n_clusters < len(labels) else 0.0
        tri = np.triu_indices(len(labels), 1)
        same = labels[tri[0]] == labels[tri[1]]
        values = distance[tri].astype(float)
        rows.append({"HC_level": level, "n_clusters": n_clusters,
                     "largest_cluster_size": int(sizes.max()), "median_cluster_size": float(np.median(sizes)),
                     "singleton_count": int(np.sum(sizes == 1)),
                     "singleton_fraction": float(np.sum(sizes == 1) / len(labels)),
                     "NMI_to_previous": np.nan if previous is None else float(normalized_mutual_info_score(previous, labels)),
                     "ARI_to_previous": np.nan if previous is None else float(adjusted_rand_score(previous, labels)),
                     "silhouette": sil,
                     "mean_within_distance": float(np.mean(values[same])) if same.any() else np.nan,
                     "mean_between_distance": float(np.mean(values[~same])) if (~same).any() else np.nan})
        previous = labels
    return pd.DataFrame(rows), assignments


def select_hc_level(scan: pd.DataFrame) -> dict:
    """Assistive stable-plateau selector with explicit three-state outcome."""
    usable = scan[(scan.n_clusters >= 2) & (scan.singleton_fraction <= 0.25)].copy()
    if usable.empty:
        return {"status": "UNRESOLVED", "selected": None, "stable_ranges": []}
    stable = usable.NMI_to_previous.fillna(0).ge(0.995) & usable.ARI_to_previous.fillna(0).ge(0.995)
    runs, current = [], []
    for idx, ok in zip(usable.index, stable):
        if ok and (not current or idx == current[-1] + 1):
            current.append(idx)
        else:
            if len(current) >= 2: runs.append(current)
            current = [idx] if ok else []
    if len(current) >= 2: runs.append(current)
    ranges = [(int(scan.loc[r[0], "HC_level"]), int(scan.loc[r[-1], "HC_level"])) for r in runs]
    if not ranges:
        return {"status": "UNRESOLVED", "selected": None, "stable_ranges": []}
    lengths = np.asarray([b - a for a, b in ranges])
    best = np.flatnonzero(lengths == lengths.max())
    if len(best) != 1:
        return {"status": "AMBIGUOUS", "selected": None, "stable_ranges": ranges}
    chosen = ranges[int(best[0])]
    return {"status": "PASS", "selected": chosen[0], "stable_ranges": ranges}

