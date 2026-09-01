from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd


CHANNELS = ("PP", "PT", "TP", "TT")
DISPLAY = {"PP": "PP", "PT": "P↔T", "TP": "P↔T", "TT": "TT"}
ORIENTATION = {
    "PP": "P_A__P_B", "PT": "P_A__T_B", "TP": "T_A__P_B", "TT": "T_A__T_B"
}


def pool_rare_states(types: np.ndarray, minimum_count: int) -> tuple[np.ndarray, list[dict]]:
    pooled = np.full(types.shape, -1, np.int32)
    audit = []
    for j in range(types.shape[1]):
        observed = types[:, j]
        values, counts = np.unique(observed[observed >= 0], return_counts=True)
        common = values[counts >= minimum_count]
        rare = values[counts < minimum_count]
        mapping = {int(v): k for k, v in enumerate(common)}
        rare_code = len(common) if len(rare) else -1
        for old, new in mapping.items():
            pooled[observed == old, j] = new
        if len(rare):
            pooled[np.isin(observed, rare), j] = rare_code
        audit.append({"locus_index": j, "original_states": len(values),
                      "pooled_states": len(common) + int(len(rare) > 0),
                      "rare_observations": int(counts[np.isin(values, rare)].sum()) if len(rare) else 0,
                      "rare_code": rare_code})
    return pooled, audit


def channel_vectors(presence: np.ndarray, types: np.ndarray, i: int, j: int, channel: str):
    pa, pb, ta, tb = presence[:, i], presence[:, j], types[:, i], types[:, j]
    if channel == "PP":
        valid = (pa >= 0) & (pb >= 0); x, y = pa, pb
    elif channel == "PT":
        valid = (pa >= 0) & (pb == 1) & (tb >= 0); x, y = pa, tb
    elif channel == "TP":
        valid = (pa == 1) & (ta >= 0) & (pb >= 0); x, y = ta, pb
    elif channel == "TT":
        valid = (pa == 1) & (pb == 1) & (ta >= 0) & (tb >= 0); x, y = ta, tb
    else:
        raise ValueError(channel)
    return x[valid].astype(np.int32), y[valid].astype(np.int32)


def contingency(x: np.ndarray, y: np.ndarray) -> dict | None:
    if len(x) == 0:
        return None
    ux, xi = np.unique(x, return_inverse=True)
    uy, yi = np.unique(y, return_inverse=True)
    table = np.zeros((len(ux), len(uy)), np.int64)
    np.add.at(table, (xi, yi), 1)
    n = int(table.sum())
    px, py = table.sum(1), table.sum(0)
    expected = px[:, None] * py[None, :] / n
    nz = table > 0
    mi = float(np.sum((table[nz] / n) * np.log(table[nz] / expected[nz])))
    qx, qy = px / n, py / n
    hx = float(-np.sum(qx[qx > 0] * np.log(qx[qx > 0])))
    hy = float(-np.sum(qy[qy > 0] * np.log(qy[qy > 0])))
    residual = np.divide(table - expected, np.sqrt(expected), out=np.zeros_like(expected), where=expected > 0)
    enrichment = np.divide(table, expected, out=np.full_like(expected, np.nan), where=expected > 0)
    ea, eb = np.unravel_index(np.argmax(residual), residual.shape)
    da, db = np.unravel_index(np.argmin(residual), residual.shape)
    denom = hx + hy
    return {
        "MIraw": max(0.0, mi), "NMI": 2.0 * max(0.0, mi) / denom if denom > 0 else np.nan,
        "informative_N": n, "entropy_A": hx, "entropy_B": hy,
        "n_states_A": len(ux), "n_states_B": len(uy),
        "min_state_count": int(min(px.min(), py.min())),
        "enriched_driver_A_state": int(ux[ea]), "enriched_driver_B_state": int(uy[eb]),
        "enriched_observed": int(table[ea, eb]), "enriched_expected": float(expected[ea, eb]),
        "enriched_residual": float(residual[ea, eb]), "enriched_ratio": float(enrichment[ea, eb]),
        "depleted_driver_A_state": int(ux[da]), "depleted_driver_B_state": int(uy[db]),
        "depleted_observed": int(table[da, db]), "depleted_expected": float(expected[da, db]),
        "depleted_residual": float(residual[da, db]), "depleted_ratio": float(enrichment[da, db]),
    }


def score_all_pairs(presence: np.ndarray, raw_types: np.ndarray, loci: np.ndarray,
                    distances: pd.DataFrame, minimum_informative: int,
                    minimum_state_count: int) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    types, rare_audit = pool_rare_states(raw_types, minimum_state_count)
    distance_index = distances.set_index(["i", "j"])
    outputs: dict[str, list[dict]] = {c: [] for c in CHANNELS}
    for i in range(len(loci) - 1):
        for j in range(i + 1, len(loci)):
            distance = distance_index.loc[(i, j)].to_dict()
            for channel in CHANNELS:
                x, y = channel_vectors(presence, types, i, j, channel)
                metrics = contingency(x, y)
                if metrics is None:
                    continue
                eligible = (metrics["informative_N"] >= minimum_informative and
                            metrics["n_states_A"] >= 2 and metrics["n_states_B"] >= 2 and
                            metrics["min_state_count"] >= minimum_state_count)
                if not eligible:
                    continue
                outputs[channel].append({
                    "edge_id": f"{channel}:{loci[i]}--{loci[j]}", "i": i, "j": j,
                    "locus_A": str(loci[i]), "locus_B": str(loci[j]), "channel": channel,
                    "display_channel": DISPLAY[channel], "orientation": ORIENTATION[channel],
                    **metrics, **distance,
                })
    return {c: pd.DataFrame(rows) for c, rows in outputs.items()}, pd.DataFrame(rare_audit)


def neff_from_counts(counts: dict[int, int]) -> float:
    values = np.asarray([v for v in counts.values() if v > 0], float)
    if values.size == 0:
        return np.nan
    p = values / values.sum()
    return float(1.0 / np.sum(p * p))


def parse_counts(text: str) -> dict[int, int]:
    if not isinstance(text, str) or not text.strip():
        return {}
    return {int(piece.split(":", 1)[0]): int(piece.split(":", 1)[1])
            for piece in text.split(";") if ":" in piece}

