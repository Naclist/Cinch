"""Optional compiled four-channel block scorer."""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - optional dependency
    njit = None


def available() -> bool:
    return njit is not None


if njit is not None:
    @njit(cache=True)
    def score_block_metrics(presence, types, pairs, presence_variable, type_variable,
                            minimum_informative, minimum_state_count, type_cardinality):
        # metric columns follow _METRIC_NAMES below; status is column 0.
        output = np.full((pairs.shape[0], 4, 21), np.nan, dtype=np.float64)
        for pair_index in range(pairs.shape[0]):
            left, right = pairs[pair_index, 0], pairs[pair_index, 1]
            for channel in range(4):
                if channel == 0:
                    if not (presence_variable[left] and presence_variable[right]):
                        continue
                    kx, ky = 2, 2
                elif channel == 1:
                    if not (presence_variable[left] and type_variable[right]):
                        continue
                    kx, ky = 2, type_cardinality[right]
                elif channel == 2:
                    if not (type_variable[left] and presence_variable[right]):
                        continue
                    kx, ky = type_cardinality[left], 2
                else:
                    if not (type_variable[left] and type_variable[right]):
                        continue
                    kx, ky = type_cardinality[left], type_cardinality[right]
                table = np.zeros((kx, ky), dtype=np.int64)
                for sample in range(presence.shape[0]):
                    pa, pb = presence[sample, left], presence[sample, right]
                    ta, tb = types[sample, left], types[sample, right]
                    valid, x, y = False, 0, 0
                    if channel == 0:
                        valid, x, y = pa >= 0 and pb >= 0, pa, pb
                    elif channel == 1:
                        valid, x, y = pa >= 0 and pb == 1 and tb >= 0, pa, tb
                    elif channel == 2:
                        valid, x, y = pa == 1 and ta >= 0 and pb >= 0, ta, pb
                    else:
                        valid, x, y = pa == 1 and pb == 1 and ta >= 0 and tb >= 0, ta, tb
                    if valid:
                        table[x, y] += 1
                n = table.sum()
                if n == 0:
                    continue
                px, py = table.sum(axis=1), table.sum(axis=0)
                nxa, nya, min_count = 0, 0, n + 1
                for value in px:
                    if value > 0:
                        nxa += 1
                        if value < min_count: min_count = value
                for value in py:
                    if value > 0:
                        nya += 1
                        if value < min_count: min_count = value
                if n < minimum_informative or nxa < 2 or nya < 2 or min_count < minimum_state_count:
                    continue
                hx, hy, mi = 0.0, 0.0, 0.0
                for value in px:
                    if value > 0:
                        probability = value / n
                        hx -= probability * np.log(probability)
                for value in py:
                    if value > 0:
                        probability = value / n
                        hy -= probability * np.log(probability)
                max_residual, min_residual = -np.inf, np.inf
                ea = eb = da = db = 0
                enriched_expected = depleted_expected = 0.0
                for a in range(kx):
                    if px[a] == 0:
                        continue
                    for b in range(ky):
                        if py[b] == 0:
                            continue
                        expected = px[a] * py[b] / n
                        if table[a, b] > 0:
                            probability = table[a, b] / n
                            mi += probability * np.log(table[a, b] / expected)
                        residual = (table[a, b] - expected) / np.sqrt(expected)
                        if residual > max_residual:
                            max_residual, ea, eb, enriched_expected = residual, a, b, expected
                        if residual < min_residual:
                            min_residual, da, db, depleted_expected = residual, a, b, expected
                if mi < 0: mi = 0.0
                denom = hx + hy
                output[pair_index, channel, 0] = 1.0
                output[pair_index, channel, 1] = mi
                output[pair_index, channel, 2] = 2.0 * mi / denom if denom > 0 else np.nan
                output[pair_index, channel, 3] = n
                output[pair_index, channel, 4] = hx
                output[pair_index, channel, 5] = hy
                output[pair_index, channel, 6] = nxa
                output[pair_index, channel, 7] = nya
                output[pair_index, channel, 8] = min_count
                output[pair_index, channel, 9] = ea
                output[pair_index, channel, 10] = eb
                output[pair_index, channel, 11] = table[ea, eb]
                output[pair_index, channel, 12] = enriched_expected
                output[pair_index, channel, 13] = max_residual
                output[pair_index, channel, 14] = table[ea, eb] / enriched_expected
                output[pair_index, channel, 15] = da
                output[pair_index, channel, 16] = db
                output[pair_index, channel, 17] = table[da, db]
                output[pair_index, channel, 18] = depleted_expected
                output[pair_index, channel, 19] = min_residual
                output[pair_index, channel, 20] = table[da, db] / depleted_expected
        return output
else:
    def score_block_metrics(*args, **kwargs):  # pragma: no cover
        raise RuntimeError("Numba engine requires the optional performance dependency")


METRIC_NAMES = (
    "MIraw", "NMI", "informative_N", "entropy_A", "entropy_B", "n_states_A",
    "n_states_B", "min_state_count", "enriched_driver_A_state",
    "enriched_driver_B_state", "enriched_observed", "enriched_expected",
    "enriched_residual", "enriched_ratio", "depleted_driver_A_state",
    "depleted_driver_B_state", "depleted_observed", "depleted_expected",
    "depleted_residual", "depleted_ratio",
)
