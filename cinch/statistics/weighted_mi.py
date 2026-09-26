"""Weighted binary mutual-information primitives from CINCH-dev2.

Source: Naclist/Cinch-dev2@ec1eaa5, ``core/statistics.py``.
The implementation preserves the frozen weighted-MI-v2 definition: inputs are
weighted masses, logarithms are natural, and no pseudocount is applied.
"""

from __future__ import annotations

import numpy as np


def weighted_binary_information(
    joint_11: np.ndarray,
    marginal_x1: np.ndarray,
    marginal_y1: np.ndarray,
    total_weight: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return weighted MI in nats, symmetric NMI, H(X), and H(Y)."""

    c11 = np.asarray(joint_11, dtype=np.float64)
    m1 = np.asarray(marginal_x1, dtype=np.float64)
    m2 = np.asarray(marginal_y1, dtype=np.float64)
    if total_weight <= 0 or not np.isfinite(total_weight):
        raise ValueError("total_weight must be finite and positive")
    if not (c11.shape == m1.shape == m2.shape):
        raise ValueError("joint and marginal arrays must have identical shapes")
    counts = np.column_stack(
        [total_weight - m1 - m2 + c11, m2 - c11, m1 - c11, c11]
    )
    counts = np.clip(counts, 0.0, total_weight)
    probabilities = counts / total_weight
    px = np.column_stack(
        [probabilities[:, 0] + probabilities[:, 1], probabilities[:, 2] + probabilities[:, 3]]
    )
    py = np.column_stack(
        [probabilities[:, 0] + probabilities[:, 2], probabilities[:, 1] + probabilities[:, 3]]
    )
    expected = np.column_stack(
        [px[:, 0] * py[:, 0], px[:, 0] * py[:, 1], px[:, 1] * py[:, 0], px[:, 1] * py[:, 1]]
    )
    terms = np.zeros_like(probabilities)
    valid = (probabilities > 0.0) & (expected > 0.0)
    terms[valid] = probabilities[valid] * np.log(probabilities[valid] / expected[valid])
    mi = terms.sum(axis=1)

    hx_terms = np.zeros_like(px)
    hy_terms = np.zeros_like(py)
    px_valid = px > 0.0
    py_valid = py > 0.0
    hx_terms[px_valid] = -px[px_valid] * np.log(px[px_valid])
    hy_terms[py_valid] = -py[py_valid] * np.log(py[py_valid])
    hx = hx_terms.sum(axis=1)
    hy = hy_terms.sum(axis=1)
    nmi = np.divide(2.0 * mi, hx + hy, out=np.zeros_like(mi), where=(hx + hy) > 0.0)
    return mi, nmi, hx, hy
