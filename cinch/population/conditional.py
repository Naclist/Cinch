"""Categorical population- and background-conditioned association kernels."""

from __future__ import annotations

import hashlib

import numpy as np


def _vectors(x, y, population=None, background=None):
    x, y = np.asarray(x), np.asarray(y)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError("x and y must be aligned one-dimensional arrays")
    output = [x, y]
    for name, values in (("population", population), ("background", background)):
        if values is not None:
            values = np.asarray(values)
            if values.ndim != 1 or values.shape != x.shape:
                raise ValueError(f"{name} must align with x")
            output.append(values)
    return tuple(output)


def categorical_mi(x: np.ndarray, y: np.ndarray) -> float:
    """Unweighted mutual information in nats for nominal categorical states."""

    x, y = _vectors(x, y)
    if x.size == 0:
        return np.nan
    _, xi = np.unique(x, return_inverse=True)
    _, yi = np.unique(y, return_inverse=True)
    table = np.zeros((int(xi.max()) + 1, int(yi.max()) + 1), dtype=np.int64)
    np.add.at(table, (xi, yi), 1)
    total = float(table.sum())
    joint = table / total
    expected = joint.sum(axis=1)[:, None] * joint.sum(axis=0)[None, :]
    valid = (joint > 0) & (expected > 0)
    return float(np.sum(joint[valid] * np.log(joint[valid] / expected[valid])))


def categorical_summary(x: np.ndarray, y: np.ndarray) -> dict:
    """MI, entropy, and observed/expected state drivers for nominal states."""

    x, y = _vectors(x, y)
    if x.size == 0:
        raise ValueError("categorical summary requires at least one sample")
    ux, xi = np.unique(x, return_inverse=True)
    uy, yi = np.unique(y, return_inverse=True)
    table = np.zeros((len(ux), len(uy)), dtype=np.int64)
    np.add.at(table, (xi, yi), 1)
    total = float(table.sum())
    joint = table / total
    px, py = joint.sum(axis=1), joint.sum(axis=0)
    expected = px[:, None] * py[None, :]
    valid = (joint > 0) & (expected > 0)
    mi = float(np.sum(joint[valid] * np.log(joint[valid] / expected[valid])))
    hx = float(-np.sum(px[px > 0] * np.log(px[px > 0])))
    hy = float(-np.sum(py[py > 0] * np.log(py[py > 0])))
    residual = np.divide(joint - expected, np.sqrt(expected / total),
                         out=np.zeros_like(expected), where=expected > 0)
    enriched = np.unravel_index(np.argmax(residual), residual.shape)
    depleted = np.unravel_index(np.argmin(residual), residual.shape)
    return {
        "mi": max(0.0, mi),
        "nmi": 2.0 * max(0.0, mi) / (hx + hy) if hx + hy > 0 else np.nan,
        "n_states_x": len(ux), "n_states_y": len(uy),
        "joint_state_combinations": int((table > 0).sum()),
        "enriched_x_state": int(ux[enriched[0]]),
        "enriched_y_state": int(uy[enriched[1]]),
        "depleted_x_state": int(ux[depleted[0]]),
        "depleted_y_state": int(uy[depleted[1]]),
    }


def population_conditioned_mi(x: np.ndarray, y: np.ndarray, population: np.ndarray) -> dict:
    """Empirically population-weighted conditional MI and support diagnostics."""

    x, y, population = _vectors(x, y, population)
    total = len(x)
    value, informative = 0.0, 0
    details = []
    for label in np.unique(population):
        mask = population == label
        count = int(mask.sum())
        mi = categorical_mi(x[mask], y[mask])
        x_states, y_states = np.unique(x[mask]).size, np.unique(y[mask]).size
        is_informative = x_states >= 2 and y_states >= 2
        informative += int(is_informative)
        value += count / total * mi
        details.append({"population": str(label), "n": count, "mi": mi,
                        "x_states": x_states, "y_states": y_states,
                        "informative": is_informative})
    return {"mi": float(value), "informative_populations": informative, "details": details}


def _aligned_joint(x, y, mask, ux, uy):
    x_index = {int(value): index for index, value in enumerate(ux)}
    y_index = {int(value): index for index, value in enumerate(uy)}
    table = np.zeros((len(ux), len(uy)), dtype=float)
    for left, right in zip(x[mask], y[mask]):
        table[x_index[int(left)], y_index[int(right)]] += 1.0
    return table / table.sum()


def _distribution_drivers(distribution: np.ndarray, ux: np.ndarray, uy: np.ndarray) -> dict:
    px, py = distribution.sum(axis=1), distribution.sum(axis=0)
    row_mask, column_mask = px > 0, py > 0
    observed = distribution[np.ix_(row_mask, column_mask)]
    local_x, local_y = ux[row_mask], uy[column_mask]
    expected = observed.sum(axis=1)[:, None] * observed.sum(axis=0)[None, :]
    residual = np.divide(observed - expected, np.sqrt(expected),
                         out=np.zeros_like(expected), where=expected > 0)
    enriched = np.unravel_index(np.argmax(residual), residual.shape)
    depleted = np.unravel_index(np.argmin(residual), residual.shape)
    return {
        "enriched_x_state": int(local_x[enriched[0]]),
        "enriched_y_state": int(local_y[enriched[1]]),
        "depleted_x_state": int(local_x[depleted[0]]),
        "depleted_y_state": int(local_y[depleted[1]]),
    }


def standardized_background_mi(
    x: np.ndarray,
    y: np.ndarray,
    population: np.ndarray,
    background: np.ndarray,
    reference: str,
    comparison: str,
    *,
    minimum_cell: int,
    minimum_shared_populations: int,
    minimum_background_samples: int,
) -> dict:
    """Common-support standardized background MI contrast."""

    x, y, population, background = _vectors(x, y, population, background)
    reference, comparison = str(reference), str(comparison)
    background = background.astype(str)
    reference_n = int(np.sum(background == reference))
    comparison_n = int(np.sum(background == comparison))
    shared, diagnostics = [], []
    for label in np.unique(population):
        pop = population == label
        masks = {value: pop & (background == value) for value in (reference, comparison)}
        counts = {value: int(mask.sum()) for value, mask in masks.items()}
        states = {
            value: (np.unique(x[mask]).size, np.unique(y[mask]).size)
            for value, mask in masks.items()
        }
        # Eligibility must not depend on the observed allocation of X/Y states
        # to the background whose labels are permuted.  Monomorphic margins
        # have a defined MI of zero and remain part of the target population.
        usable = all(counts[value] >= minimum_cell for value in (reference, comparison))
        diagnostics.append({
            "population": str(label), "reference_n": counts[reference],
            "comparison_n": counts[comparison], "reference_x_states": states[reference][0],
            "reference_y_states": states[reference][1], "comparison_x_states": states[comparison][0],
            "comparison_y_states": states[comparison][1],
            "shared_by_sample_support": usable,
            "state_variable_in_both_backgrounds": all(
                states[value][0] >= 2 and states[value][1] >= 2
                for value in (reference, comparison)
            ),
        })
        if usable:
            shared.append((label, min(counts[reference], counts[comparison])))
    result = {
        "status": "TESTED", "reason": "", "reference_n": reference_n,
        "comparison_n": comparison_n, "shared_population_count": len(shared),
        "effective_sample_support": int(sum(value for _, value in shared)),
        "diagnostics": diagnostics,
    }
    if reference_n < minimum_background_samples or comparison_n < minimum_background_samples:
        result.update(status="NOT_IDENTIFIABLE", reason="INSUFFICIENT_BACKGROUND_SAMPLES")
        return result
    if len(shared) < minimum_shared_populations:
        result.update(status="NOT_IDENTIFIABLE", reason="INSUFFICIENT_SHARED_POPULATIONS")
        return result
    masses = np.asarray([value for _, value in shared], dtype=float)
    weights = masses / masses.sum()
    labels = [label for label, _ in shared]
    ux, uy = np.unique(x), np.unique(y)
    mi = {reference: 0.0, comparison: 0.0}
    distributions = {
        reference: np.zeros((len(ux), len(uy)), dtype=float),
        comparison: np.zeros((len(ux), len(uy)), dtype=float),
    }
    combinations = {reference: 0, comparison: 0}
    for label, weight in zip(labels, weights):
        for value in (reference, comparison):
            mask = (population == label) & (background == value)
            mi[value] += weight * categorical_mi(x[mask], y[mask])
            distributions[value] += weight * _aligned_joint(x, y, mask, ux, uy)
            combinations[value] += int(len(set(zip(x[mask].tolist(), y[mask].tolist()))))
    result.update({
        "reference_mi": float(mi[reference]), "comparison_mi": float(mi[comparison]),
        "delta_mi": float(mi[comparison] - mi[reference]),
        "shared_populations": labels, "population_weights": weights,
        "effective_population_count": float(1.0 / np.sum(weights * weights)),
        "reference_joint_combinations": combinations[reference],
        "comparison_joint_combinations": combinations[comparison],
        "reference_drivers": _distribution_drivers(distributions[reference], ux, uy),
        "comparison_drivers": _distribution_drivers(distributions[comparison], ux, uy),
    })
    return result


def fixed_standardized_delta(
    x, y, population, background, reference, comparison, shared_populations, population_weights,
) -> float:
    """Compute delta on a fixed observed support/standardization target."""

    x, y, population, background = _vectors(x, y, population, background)
    background = background.astype(str)
    values = {str(reference): 0.0, str(comparison): 0.0}
    for label, weight in zip(shared_populations, population_weights):
        for value in values:
            mask = (population == label) & (background == value)
            values[value] += float(weight) * categorical_mi(x[mask], y[mask])
    return float(values[str(comparison)] - values[str(reference)])


def population_permutation_test(x, y, population, permutations: int, rng) -> dict:
    """Test residual association by permuting Y within population."""

    if permutations < 1:
        raise ValueError("permutations must be >= 1")
    x, y, population = _vectors(x, y, population)
    observed = population_conditioned_mi(x, y, population)["mi"]
    null = np.empty(permutations, dtype=float)
    groups = [np.flatnonzero(population == label) for label in np.unique(population)]
    for index in range(permutations):
        shuffled = y.copy()
        for indices in groups:
            shuffled[indices] = rng.permutation(y[indices])
        null[index] = population_conditioned_mi(x, shuffled, population)["mi"]
    pvalue = (1.0 + np.sum(null >= observed)) / (permutations + 1.0)
    return {"observed": observed, "null_mean": float(null.mean()), "p": float(pvalue)}


def differential_permutation_test(
    x, y, population, background, reference, comparison, shared_populations,
    population_weights, permutations: int, rng,
) -> dict:
    """Two-sided background-label randomization within population."""

    if permutations < 1:
        raise ValueError("permutations must be >= 1")
    x, y, population, background = _vectors(x, y, population, background)
    observed = fixed_standardized_delta(
        x, y, population, background, reference, comparison,
        shared_populations, population_weights,
    )
    null = np.empty(permutations, dtype=float)
    groups = [np.flatnonzero(population == label) for label in shared_populations]
    for index in range(permutations):
        shuffled = background.copy()
        for indices in groups:
            shuffled[indices] = rng.permutation(background[indices])
        null[index] = fixed_standardized_delta(
            x, y, population, shuffled, reference, comparison,
            shared_populations, population_weights,
        )
    pvalue = (1.0 + np.sum(np.abs(null) >= abs(observed))) / (permutations + 1.0)
    return {"observed": observed, "null_mean": float(null.mean()), "p": float(pvalue)}


def stable_hypothesis_seed(seed: int, text: str) -> int:
    return int.from_bytes(hashlib.sha256(f"{seed}|{text}".encode()).digest()[:8], "little")
