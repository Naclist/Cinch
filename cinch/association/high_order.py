"""Frozen SHC-conditioned background/modifier statistic for Diff-GWES."""

from __future__ import annotations

import hashlib

import numpy as np

from cinch.population.shc import weighted_mi_samples


LN2 = np.log(2.0)


def make_cluster_indices(shc: np.ndarray, minimum_size: int = 2) -> list[np.ndarray]:
    shc = np.asarray(shc)
    return [
        indices for cluster in np.unique(shc)
        if len(indices := np.flatnonzero(shc == cluster)) >= minimum_size
    ]


def background_stat(
    x: np.ndarray,
    y: np.ndarray,
    c: np.ndarray,
    weights: np.ndarray,
    shc: np.ndarray,
    min_group: int = 25,
    min_shc_cell: int = 3,
    min_informative: int = 2,
    cluster_indices: list[np.ndarray] | None = None,
) -> dict:
    x, y, c = np.asarray(x, dtype=bool), np.asarray(y, dtype=bool), np.asarray(c, dtype=bool)
    weights, shc = np.asarray(weights, dtype=float), np.asarray(shc)
    if not (x.ndim == 1 and x.shape == y.shape == c.shape == weights.shape == shc.shape):
        raise ValueError("x, y, c, weights, and shc must be aligned one-dimensional arrays")
    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("weights must be finite and non-negative")
    result = {"n0": int((~c).sum()), "n1": int(c.sum())}
    groups = cluster_indices if cluster_indices is not None else make_cluster_indices(shc, 1)
    for state in (0, 1):
        numerator, denominator, informative = 0.0, 0.0, 0
        for indices in groups:
            selected = indices[c[indices] == bool(state)]
            size = len(selected)
            if size < min_shc_cell:
                continue
            sx, sy = int(x[selected].sum()), int(y[selected].sum())
            if not (0 < sx < size and 0 < sy < size):
                continue
            mass = float(weights[selected].sum())
            mi = weighted_mi_samples(x[selected], y[selected], weights[selected])
            if np.isfinite(mi) and mass > 0:
                numerator += mass * mi
                denominator += mass
                informative += 1
        conditional = numerator / denominator if denominator > 0 else np.nan
        result[f"informative_shc_{state}"] = informative
        result[f"conditional_mi_{state}"] = conditional
        result[f"epidis_{state}"] = np.sqrt(max(conditional, 0.0) / LN2) if np.isfinite(conditional) else np.nan
    result["eligible"] = bool(
        result["n0"] >= min_group and result["n1"] >= min_group
        and result["informative_shc_0"] >= min_informative
        and result["informative_shc_1"] >= min_informative
    )
    result["delta"] = result["epidis_1"] - result["epidis_0"] if result["eligible"] else np.nan
    return result


def stable_seed(text: str) -> int:
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def permuted_background(
    c: np.ndarray, cluster_indices: list[np.ndarray], rng: np.random.Generator,
) -> np.ndarray:
    c = np.asarray(c)
    output = c.copy()
    for indices in cluster_indices:
        output[indices] = c[indices][rng.permutation(len(indices))]
    return output
