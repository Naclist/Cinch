"""Sparse same-contig order and base-pair distance aggregation."""

from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd

from .numba_order import aggregate_sparse_distances, available as numba_available


class SparseOrderDistance:
    """Store only locus pairs observed on the same sample contig."""

    def __init__(self, coordinates: pd.DataFrame, n_loci: int, minimum_observations: int,
                 engine: str = "python"):
        required = {"sample_id", "contig", "locus_index", "start", "end", "gene_order"}
        if not required.issubset(coordinates.columns):
            raise ValueError(f"coordinates are missing columns: {sorted(required - set(coordinates.columns))}")
        if minimum_observations < 1:
            raise ValueError("minimum_observations must be >= 1")
        if engine not in {"python", "numba"}:
            raise ValueError("engine must be 'python' or 'numba'")
        if engine == "numba" and not numba_available():
            raise RuntimeError("Numba engine requires `pip install -e '.[performance]'`")
        self.n_loci = int(n_loci)
        self.minimum_observations = int(minimum_observations)
        locus_indices = coordinates["locus_index"].to_numpy(dtype=np.int64)
        if ((locus_indices < 0) | (locus_indices >= self.n_loci)).any():
            raise ValueError("coordinate locus_index is outside profile locus range")
        self._observed = (
            self._aggregate_numba(coordinates) if engine == "numba"
            else self._aggregate_python(coordinates)
        )

    def _aggregate_python(self, coordinates: pd.DataFrame):
        accumulated: dict[tuple[int, int], list[tuple[float, float]]] = {}
        for (_, _), group in coordinates.groupby(["sample_id", "contig"], sort=False):
            records = group[["locus_index", "start", "end", "gene_order"]].to_numpy()
            local: dict[tuple[int, int], list[tuple[float, float]]] = {}
            for a in range(len(records) - 1):
                for b in range(a + 1, len(records)):
                    left, right = int(records[a, 0]), int(records[b, 0])
                    if left == right:
                        continue
                    key = (min(left, right), max(left, right))
                    bp = max(0.0, max(records[a, 1], records[b, 1]) - min(records[a, 2], records[b, 2]) - 1)
                    order = abs(records[a, 3] - records[b, 3])
                    local.setdefault(key, []).append((float(order), float(bp)))
            for key, values in local.items():
                accumulated.setdefault(key, []).append((min(x[0] for x in values), min(x[1] for x in values)))
        return {
            key: (
                float(np.mean([value[0] for value in values])),
                float(np.mean([value[1] for value in values])),
                len(values),
            )
            for key, values in accumulated.items()
        }

    def _aggregate_numba(self, coordinates: pd.DataFrame):
        groups, _ = pd.factorize(
            pd.MultiIndex.from_frame(coordinates[["sample_id", "contig"]]), sort=False,
        )
        order = np.argsort(groups, kind="stable")
        keys, mean_order, mean_bp, observations = aggregate_sparse_distances(
            np.asarray(groups[order], dtype=np.int64),
            coordinates["locus_index"].to_numpy(dtype=np.int64)[order],
            coordinates["start"].to_numpy(dtype=np.float64)[order],
            coordinates["end"].to_numpy(dtype=np.float64)[order],
            coordinates["gene_order"].to_numpy(dtype=np.float64)[order],
            self.n_loci,
        )
        return {
            (int(key // self.n_loci), int(key % self.n_loci)): (
                float(order_distance), float(bp_distance), int(count),
            )
            for key, order_distance, bp_distance, count
            in zip(keys, mean_order, mean_bp, observations)
        }

    def get_pair(self, left: int, right: int) -> dict:
        key = (left, right) if left < right else (right, left)
        value = self._observed.get(key)
        if value is None:
            return {
                "order_distance": np.nan, "bp_distance": np.nan,
                "n_same_contig_observations": 0,
                "physical_distance_status": "NO_SAME_CONTIG_OBSERVATION",
            }
        order, bp, count = value
        if count < self.minimum_observations:
            return {
                "order_distance": np.nan, "bp_distance": np.nan,
                "n_same_contig_observations": count,
                "physical_distance_status": "INSUFFICIENT_SAME_CONTIG_OBSERVATIONS",
            }
        return {
            "order_distance": order, "bp_distance": bp,
            "n_same_contig_observations": count,
            "physical_distance_status": "RELIABLE_SAME_CONTIG_DISTANCE",
        }

    def digest(self) -> str:
        payload = {
            "n_loci": self.n_loci,
            "minimum_observations": self.minimum_observations,
            "observed": [[*key, *value] for key, value in sorted(self._observed.items())],
        }
        return hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()

    def observed_frame(self, loci: np.ndarray | None = None) -> pd.DataFrame:
        rows = []
        for (left, right), _ in sorted(self._observed.items()):
            row = {"i": left, "j": right, **self.get_pair(left, right)}
            if loci is not None:
                row.update({"locus_A": str(loci[left]), "locus_B": str(loci[right])})
            rows.append(row)
        return pd.DataFrame(rows)
