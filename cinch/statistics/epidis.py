"""Frozen directional EpiDis construction from CINCH-dev2.

Source: Naclist/Cinch-dev2@ec1eaa5, ``core/statistics.py``.
"""

from __future__ import annotations

import numpy as np


EPSILON = 1e-27


def _xlog2_ratio(x: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    out = np.zeros_like(x)
    valid = (x > 0.0) & (denominator > 0.0)
    out[valid] = x[valid] * np.log2(x[valid] / denominator[valid])
    return out


def directional_epidis(
    joint_11: np.ndarray,
    conditioning_1: np.ndarray,
    target_1: np.ndarray,
    epsilon: float = EPSILON,
) -> np.ndarray:
    """Evaluate the frozen EpiDis equations in one conditioning direction."""

    joint_11 = np.asarray(joint_11, dtype=np.float64)
    alpha = np.asarray(conditioning_1, dtype=np.float64)
    target_1 = np.asarray(target_1, dtype=np.float64)
    if not (joint_11.shape == alpha.shape == target_1.shape):
        raise ValueError("joint and marginal arrays must have identical shapes")
    if epsilon < 0 or not np.isfinite(epsilon):
        raise ValueError("epsilon must be finite and non-negative")
    beta = 1.0 - alpha
    valid = (alpha > 0.0) & (beta > 0.0)
    result = np.full(alpha.shape, np.nan, dtype=np.float64)
    if not np.any(valid):
        return result

    a = alpha[valid]
    b = beta[valid]
    c11 = joint_11[valid]
    t1 = target_1[valid]
    p1 = np.clip(c11 / a, 0.0, 1.0)
    q1 = np.clip((t1 - c11) / b, 0.0, 1.0)
    p = np.column_stack((p1 + epsilon, 1.0 - p1 + epsilon))
    q = np.column_stack((q1 + epsilon, 1.0 - q1 + epsilon))
    mixture = a[:, None] * p + b[:, None] * q
    jsd = a * _xlog2_ratio(p, mixture).sum(axis=1) + b * _xlog2_ratio(q, mixture).sum(axis=1)
    result[valid] = np.sqrt(np.maximum(jsd, 0.0))
    return result
