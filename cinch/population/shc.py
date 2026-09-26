"""Frozen SHC-conditioned binary MI and within-cluster permutation.

Source: Naclist/Cinch-dev2@ec1eaa5, hiercc_shc_phylo_correct_mi.py.
"""

from __future__ import annotations

import numpy as np


def _validated_vectors(x, y, weights, clusters=None):
    x = np.asarray(x, dtype=bool)
    y = np.asarray(y, dtype=bool)
    weights = np.asarray(weights, dtype=float)
    if x.ndim != 1 or not (x.shape == y.shape == weights.shape):
        raise ValueError("x, y, and weights must be aligned one-dimensional arrays")
    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("weights must be finite and non-negative")
    if clusters is None:
        return x, y, weights
    clusters = np.asarray(clusters)
    if clusters.ndim != 1 or clusters.shape != x.shape:
        raise ValueError("clusters must align with x")
    return x, y, weights, clusters


def weighted_mi_samples(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    x, y, weights = _validated_vectors(x, y, weights)
    total = weights.sum()
    if total <= 0:
        return 0.0
    joint = np.array([
        weights[(~x) & (~y)].sum(), weights[(~x) & y].sum(),
        weights[x & (~y)].sum(), weights[x & y].sum(),
    ], dtype=float).reshape(2, 2) / total
    px, py = joint.sum(axis=1), joint.sum(axis=0)
    expected = px[:, None] * py[None, :]
    valid = (joint > 0) & (expected > 0)
    return float(np.sum(joint[valid] * np.log(joint[valid] / expected[valid])))


def conditional_mi(
    x: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    clusters: np.ndarray,
    eligible_clusters: list[int] | None = None,
) -> float:
    x, y, weights, clusters = _validated_vectors(x, y, weights, clusters)
    total = weights.sum()
    if total <= 0:
        return 0.0
    value = 0.0
    values = np.unique(clusters) if eligible_clusters is None else eligible_clusters
    for cluster in values:
        mask = clusters == cluster
        cluster_weight = weights[mask].sum()
        if cluster_weight > 0:
            value += cluster_weight / total * weighted_mi_samples(x[mask], y[mask], weights[mask])
    return float(value)


def informative_clusters(
    x: np.ndarray, y: np.ndarray, clusters: np.ndarray, minimum_size: int,
) -> tuple[list[int], int]:
    x, y, _, clusters = _validated_vectors(x, y, np.ones(len(x)), clusters)
    informative, positive = [], 0
    for cluster in np.unique(clusters):
        mask = clusters == cluster
        size = int(mask.sum())
        if size < minimum_size:
            continue
        sx, sy = int(x[mask].sum()), int(y[mask].sum())
        if not (0 < sx < size and 0 < sy < size):
            continue
        informative.append(int(cluster))
        n11 = int(np.sum(x[mask] & y[mask]))
        n10 = int(np.sum(x[mask] & ~y[mask]))
        n01 = int(np.sum(~x[mask] & y[mask]))
        n00 = int(np.sum(~x[mask] & ~y[mask]))
        if n11 * n00 > n10 * n01:
            positive += 1
    return informative, positive


def phylo_permutation(
    x: np.ndarray,
    y: np.ndarray,
    weights: np.ndarray,
    clusters: np.ndarray,
    eligible_clusters: list[int],
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float, float, float]:
    if permutations < 1:
        raise ValueError("permutations must be >= 1")
    x, y, weights, clusters = _validated_vectors(x, y, weights, clusters)
    observed = conditional_mi(x, y, weights, clusters, eligible_clusters)
    null = np.empty(permutations, dtype=float)
    eligible_mask = np.isin(clusters, eligible_clusters)
    for permutation in range(permutations):
        shuffled = y.copy()
        for cluster in eligible_clusters:
            indices = np.flatnonzero(clusters == cluster)
            shuffled[indices] = rng.permutation(shuffled[indices])
        null[permutation] = conditional_mi(
            x[eligible_mask], shuffled[eligible_mask], weights[eligible_mask],
            clusters[eligible_mask], eligible_clusters,
        )
    pvalue = (1.0 + np.sum(null >= observed)) / (permutations + 1.0)
    standard_deviation = float(null.std(ddof=1))
    zscore = (observed - float(null.mean())) / standard_deviation if standard_deviation > 0 else np.inf
    return observed, float(null.mean()), zscore, float(pvalue)
