"""Multiple-testing procedures preserved from CINCH-dev2."""

from __future__ import annotations

import numpy as np


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    """Return Benjamini-Hochberg adjusted values using stable tie ordering."""

    pvalues = np.asarray(pvalues, dtype=np.float64)
    if pvalues.ndim != 1:
        raise ValueError("pvalues must be one-dimensional")
    if np.any(~np.isfinite(pvalues)) or np.any((pvalues < 0.0) | (pvalues > 1.0)):
        raise ValueError("pvalues must be finite values in [0, 1]")
    if pvalues.size == 0:
        return pvalues.copy()
    order = np.argsort(pvalues, kind="mergesort")
    ranked = pvalues[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output
