"""Optional compiled sparse same-contig distance aggregation."""

from __future__ import annotations

import numpy as np

try:
    from numba import njit, types
    from numba.typed import Dict
except ImportError:  # pragma: no cover - optional dependency
    njit = None
    types = None
    Dict = None


def available() -> bool:
    return njit is not None


if njit is not None:
    @njit(cache=True)
    def aggregate_sparse_distances(group_ids, loci, starts, ends, orders, n_loci):
        order_sums = Dict.empty(key_type=types.int64, value_type=types.float64)
        bp_sums = Dict.empty(key_type=types.int64, value_type=types.float64)
        counts = Dict.empty(key_type=types.int64, value_type=types.int64)
        group_start = 0
        while group_start < group_ids.size:
            group_end = group_start + 1
            while group_end < group_ids.size and group_ids[group_end] == group_ids[group_start]:
                group_end += 1
            local_order = Dict.empty(key_type=types.int64, value_type=types.float64)
            local_bp = Dict.empty(key_type=types.int64, value_type=types.float64)
            for a in range(group_start, group_end - 1):
                for b in range(a + 1, group_end):
                    left, right = loci[a], loci[b]
                    if left == right:
                        continue
                    if left > right:
                        left, right = right, left
                    key = left * n_loci + right
                    order_distance = abs(orders[a] - orders[b])
                    bp_distance = max(0.0, max(starts[a], starts[b]) - min(ends[a], ends[b]) - 1.0)
                    if key not in local_order or order_distance < local_order[key]:
                        local_order[key] = order_distance
                    if key not in local_bp or bp_distance < local_bp[key]:
                        local_bp[key] = bp_distance
            for key in local_order:
                order_sums[key] = order_sums.get(key, 0.0) + local_order[key]
                bp_sums[key] = bp_sums.get(key, 0.0) + local_bp[key]
                counts[key] = counts.get(key, 0) + 1
            group_start = group_end
        size = len(counts)
        keys = np.empty(size, dtype=np.int64)
        mean_order = np.empty(size, dtype=np.float64)
        mean_bp = np.empty(size, dtype=np.float64)
        observations = np.empty(size, dtype=np.int64)
        index = 0
        for key in counts:
            count = counts[key]
            keys[index] = key
            mean_order[index] = order_sums[key] / count
            mean_bp[index] = bp_sums[key] / count
            observations[index] = count
            index += 1
        return keys, mean_order, mean_bp, observations
else:
    def aggregate_sparse_distances(*args, **kwargs):  # pragma: no cover
        raise RuntimeError("compiled order distance requires the optional performance dependency")
