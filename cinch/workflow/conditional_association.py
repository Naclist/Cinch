"""Population- and categorical-background-conditioned four-channel workflow."""

from __future__ import annotations

import hashlib
import itertools
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from cinch.frozen_v1.statistics import CHANNELS, pool_rare_states
from cinch.population import (
    categorical_summary,
    differential_permutation_test,
    population_conditioned_mi,
    population_permutation_test,
    stable_hypothesis_seed,
    standardized_background_mi,
)
from cinch.profiles import StateProfile, load_profile
from cinch.statistics import bh_adjust


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_text(text: str, path: Path) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def _atomic_tsv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    frame.to_csv(temporary, sep="\t", index=False)
    os.replace(temporary, path)


def _atomic_json(payload: dict, path: Path) -> None:
    _atomic_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", path)


def _metadata_qc(metadata_path: Path, samples: np.ndarray, population_column: str,
                 background_column: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = pd.read_csv(metadata_path, sep=None, engine="python", dtype=str)
    required = {"sample_id", population_column, background_column}
    missing_columns = sorted(required - set(metadata.columns))
    if missing_columns:
        qc = pd.DataFrame([{
            "check": "required_columns", "status": "FAIL", "count": len(missing_columns),
            "details": ",".join(missing_columns),
        }])
        return metadata, qc
    duplicated = sorted(metadata.loc[metadata.sample_id.duplicated(keep=False), "sample_id"].dropna().unique())
    profile_ids, metadata_ids = set(samples.astype(str)), set(metadata.sample_id.dropna().astype(str))
    missing = sorted(profile_ids - metadata_ids)
    extra = sorted(metadata_ids - profile_ids)
    missing_population = metadata[population_column].isna() | metadata[population_column].str.strip().eq("")
    missing_background = metadata[background_column].isna() | metadata[background_column].str.strip().eq("")
    rows = [
        {"check": "required_columns", "status": "PASS", "count": 0, "details": ""},
        {"check": "duplicate_sample_ids", "status": "FAIL" if duplicated else "PASS",
         "count": len(duplicated), "details": ",".join(duplicated[:20])},
        {"check": "profile_samples_missing_metadata", "status": "FAIL" if missing else "PASS",
         "count": len(missing), "details": ",".join(missing[:20])},
        {"check": "metadata_samples_not_in_profile", "status": "FAIL" if extra else "PASS",
         "count": len(extra), "details": ",".join(extra[:20])},
        {"check": "missing_population_values", "status": "FAIL" if missing_population.any() else "PASS",
         "count": int(missing_population.sum()),
         "details": ",".join(metadata.loc[missing_population, "sample_id"].dropna().head(20))},
        {"check": "missing_background_values", "status": "FAIL" if missing_background.any() else "PASS",
         "count": int(missing_background.sum()),
         "details": ",".join(metadata.loc[missing_background, "sample_id"].dropna().head(20))},
    ]
    return metadata, pd.DataFrame(rows)


def _channel_states(profile: StateProfile, pooled_types: np.ndarray, left: int, right: int,
                    channel: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    pa, pb = profile.presence[:, left], profile.presence[:, right]
    ta, tb = pooled_types[:, left], pooled_types[:, right]
    if channel == "PP":
        valid, x, y = (pa >= 0) & (pb >= 0), pa, pb
    elif channel == "PT":
        valid, x, y = (pa >= 0) & (pb == 1) & (tb >= 0), pa, tb
    elif channel == "TP":
        valid, x, y = (pa == 1) & (ta >= 0) & (pb >= 0), ta, pb
    elif channel == "TT":
        valid, x, y = (pa == 1) & (pb == 1) & (ta >= 0) & (tb >= 0), ta, tb
    else:
        raise ValueError(f"unknown channel: {channel}")
    return x[valid].astype(np.int32), y[valid].astype(np.int32), valid


def _candidate_hypotheses(profile: StateProfile, channels: tuple[str, ...], pairs_path: Path | None):
    locus_index = {name: index for index, name in enumerate(profile.loci)}
    if pairs_path is None:
        return [
            (left, right, channel)
            for left, right in itertools.combinations(range(len(profile.loci)), 2)
            for channel in channels
        ]
    table = pd.read_csv(pairs_path, sep=None, engine="python")
    required_names = {"locus_A", "locus_B"}
    required_indices = {"i", "j"}
    if not (required_names.issubset(table.columns) or required_indices.issubset(table.columns)):
        raise ValueError("pairs must contain locus_A/locus_B or i/j")
    output = set()
    for row in table.itertuples(index=False):
        if required_indices.issubset(table.columns):
            left, right = int(row.i), int(row.j)
        else:
            if str(row.locus_A) not in locus_index or str(row.locus_B) not in locus_index:
                raise ValueError(f"unknown candidate locus: {row.locus_A}, {row.locus_B}")
            left, right = locus_index[str(row.locus_A)], locus_index[str(row.locus_B)]
        if left == right or not (0 <= left < len(profile.loci) and 0 <= right < len(profile.loci)):
            raise ValueError(f"invalid candidate pair: {(left, right)}")
        left, right = min(left, right), max(left, right)
        row_channels = (str(row.channel),) if "channel" in table.columns else channels
        for channel in row_channels:
            if channel not in CHANNELS:
                raise ValueError(f"unknown candidate channel: {channel}")
            if channel in channels:
                output.add((left, right, channel))
    return sorted(output)


def _overlap_table(population: np.ndarray, background: np.ndarray, reference: str,
                   comparison: str, minimum_cell: int) -> tuple[pd.DataFrame, dict]:
    rows = []
    for label in np.unique(population):
        pop = population == label
        reference_n = int(np.sum(pop & (background == reference)))
        comparison_n = int(np.sum(pop & (background == comparison)))
        shared = reference_n >= minimum_cell and comparison_n >= minimum_cell
        rows.append({
            "population": str(label), "reference_background": reference,
            "comparison_background": comparison, "reference_N": reference_n,
            "comparison_N": comparison_n, "common_support_N": min(reference_n, comparison_n),
            "shared_by_sample_support": shared,
        })
    table = pd.DataFrame(rows)
    shared = table[table.shared_by_sample_support]
    masses = shared.common_support_N.to_numpy(float)
    summary = {
        "reference_N": int(np.sum(background == reference)),
        "comparison_N": int(np.sum(background == comparison)),
        "population_count": len(table), "shared_population_count": len(shared),
        "effective_sample_support": int(masses.sum()),
        "effective_population_count": float((masses.sum() ** 2) / np.sum(masses ** 2)) if len(masses) else 0.0,
    }
    return table, summary


def _minimum_marginal_count(x: np.ndarray, y: np.ndarray) -> int:
    return int(min(np.unique(x, return_counts=True)[1].min(), np.unique(y, return_counts=True)[1].min()))


def run_conditional_association(
    profile_path: Path,
    metadata_path: Path,
    output: Path,
    *,
    population_column: str,
    background_column: str,
    reference_background: str,
    comparison_background: str,
    channels: tuple[str, ...] = CHANNELS,
    pairs_path: Path | None = None,
    permutations: int = 999,
    seed: int = 42,
    minimum_eligible_samples: int = 20,
    minimum_population_cell: int = 4,
    minimum_shared_populations: int = 2,
    minimum_background_samples: int = 10,
    minimum_state_count: int = 1,
) -> Path:
    """Execute audited population and categorical-background association tests."""

    channels = tuple(dict.fromkeys(channels))
    if not channels or any(channel not in CHANNELS for channel in channels):
        raise ValueError(f"channels must be selected from {CHANNELS}")
    if reference_background == comparison_background:
        raise ValueError("reference and comparison backgrounds must differ")
    if min(permutations, minimum_eligible_samples, minimum_population_cell,
           minimum_shared_populations, minimum_background_samples, minimum_state_count) < 1:
        raise ValueError("permutation and support parameters must be positive")
    profile_path, metadata_path, output = map(Path, (profile_path, metadata_path, output))
    profile_path, metadata_path, output = profile_path.resolve(), metadata_path.resolve(), output.resolve()
    pairs_path = Path(pairs_path).resolve() if pairs_path else None
    output.mkdir(parents=True, exist_ok=True)
    profile = load_profile(profile_path)
    metadata, metadata_qc = _metadata_qc(
        metadata_path, profile.samples, population_column, background_column,
    )
    _atomic_tsv(metadata_qc, output / "METADATA_QC.tsv")
    if (metadata_qc.status == "FAIL").any():
        _atomic_json({
            "schema": "CINCH_CONDITIONAL_ASSOCIATION_V1", "status": "FAILED_METADATA_QC",
            "profile_sha256": _sha256(profile_path), "metadata_sha256": _sha256(metadata_path),
        }, output / "MANIFEST.json")
        raise ValueError("metadata failed exact sample/schema QC; inspect METADATA_QC.tsv")
    aligned = metadata.set_index("sample_id").loc[profile.samples]
    population = aligned[population_column].astype(str).to_numpy()
    background = aligned[background_column].astype(str).to_numpy()
    observed_backgrounds = set(background)
    if reference_background not in observed_backgrounds or comparison_background not in observed_backgrounds:
        raise ValueError("reference and comparison backgrounds must both occur in metadata")
    overlap, overlap_summary = _overlap_table(
        population, background, reference_background, comparison_background,
        minimum_population_cell,
    )
    _atomic_tsv(overlap, output / "POPULATION_BACKGROUND_OVERLAP.tsv")
    pooled_types, _ = pool_rare_states(profile.types, minimum_state_count)
    hypotheses = _candidate_hypotheses(profile, channels, pairs_path)
    records, tested, excluded = [], [], []
    for left, right, channel in hypotheses:
        x, y, valid = _channel_states(profile, pooled_types, left, right, channel)
        local_population = population[valid]
        local_background = background[valid]
        hypothesis_id = f"{channel}:{profile.loci[left]}--{profile.loci[right]}"
        base = {
            "hypothesis_id": hypothesis_id, "edge_id": hypothesis_id,
            "i": left, "j": right, "locus_A": str(profile.loci[left]),
            "locus_B": str(profile.loci[right]), "channel": channel,
            "population_column": population_column, "background_variable": background_column,
            "reference_background": reference_background,
            "comparison_background": comparison_background, "informative_N": int(valid.sum()),
            "permutations": permutations,
        }
        record = dict(base)
        record.update({
            "global_MI": np.nan, "global_NMI": np.nan,
            "population_conditioned_MI": np.nan, "population_raw_p": np.nan,
            "population_adjusted_q": np.nan, "population_test_status": "NOT_TESTED",
            "population_exclusion_reason": "", "background_0_conditioned_MI": np.nan,
            "background_1_conditioned_MI": np.nan, "delta_MI": np.nan,
            "delta_direction": "NOT_ESTIMATED", "shared_population_count": 0,
            "effective_sample_support": 0, "effective_population_count": np.nan,
            "raw_p": np.nan, "adjusted_q": np.nan, "test_status": "NOT_TESTED",
            "exclusion_reason": "",
        })
        base_reason = ""
        if len(x) < minimum_eligible_samples:
            base_reason = "INSUFFICIENT_ELIGIBLE_SAMPLES"
        elif np.unique(x).size < 2 or np.unique(y).size < 2:
            base_reason = "INSUFFICIENT_GLOBAL_STATE_VARIATION"
        elif _minimum_marginal_count(x, y) < minimum_state_count:
            base_reason = "INSUFFICIENT_GLOBAL_STATE_SUPPORT"
        if base_reason:
            record.update(population_exclusion_reason=base_reason, exclusion_reason=base_reason,
                          population_test_status="NOT_IDENTIFIABLE", test_status="NOT_IDENTIFIABLE")
            for test_type in ("population_conditioned", "differential"):
                excluded.append(base | {"test_type": test_type, "test_status": "NOT_IDENTIFIABLE",
                                        "exclusion_reason": base_reason})
            records.append(record)
            continue
        global_summary = categorical_summary(x, y)
        record.update({
            "global_MI": global_summary["mi"], "global_NMI": global_summary["nmi"],
            "global_enriched_A_state": global_summary["enriched_x_state"],
            "global_enriched_B_state": global_summary["enriched_y_state"],
            "global_depleted_A_state": global_summary["depleted_x_state"],
            "global_depleted_B_state": global_summary["depleted_y_state"],
            "global_joint_state_combinations": global_summary["joint_state_combinations"],
        })
        population_result = population_conditioned_mi(x, y, local_population)
        record["population_conditioned_MI"] = population_result["mi"]
        record["population_informative_count"] = population_result["informative_populations"]
        record["eligible_population_support_json"] = json.dumps(
            population_result["details"], sort_keys=True, separators=(",", ":"),
        )
        eligible_population_details = [
            item for item in population_result["details"]
            if item["n"] >= minimum_population_cell and item["informative"]
        ]
        if len(eligible_population_details) >= minimum_shared_populations:
            population_test = population_permutation_test(
                x, y, local_population, permutations,
                np.random.default_rng(stable_hypothesis_seed(seed, f"population|{hypothesis_id}")),
            )
            record.update(population_test_status="TESTED", population_raw_p=population_test["p"],
                          population_null_mean=population_test["null_mean"])
            tested.append(base | {
                "test_type": "population_conditioned", "metadata_variable": population_column,
                "contrast": "within_population_Y_permutation", "test_status": "TESTED",
                "raw_p": population_test["p"],
                "testing_family": f"population|{population_column}|{channel}",
            })
        else:
            reason = "INSUFFICIENT_INFORMATIVE_POPULATIONS"
            record.update(population_test_status="NOT_IDENTIFIABLE", population_exclusion_reason=reason)
            excluded.append(base | {"test_type": "population_conditioned",
                                    "metadata_variable": population_column,
                                    "contrast": "within_population_Y_permutation",
                                    "test_status": "NOT_IDENTIFIABLE", "exclusion_reason": reason})
        contrast_mask = np.isin(local_background, [reference_background, comparison_background])
        differential = standardized_background_mi(
            x[contrast_mask], y[contrast_mask], local_population[contrast_mask],
            local_background[contrast_mask], reference_background, comparison_background,
            minimum_cell=minimum_population_cell,
            minimum_shared_populations=minimum_shared_populations,
            minimum_background_samples=minimum_background_samples,
        )
        record.update({
            "reference_background_N": differential["reference_n"],
            "comparison_background_N": differential["comparison_n"],
            "reference_population_count": int(np.unique(
                local_population[contrast_mask][local_background[contrast_mask] == reference_background]
            ).size),
            "comparison_population_count": int(np.unique(
                local_population[contrast_mask][local_background[contrast_mask] == comparison_background]
            ).size),
            "shared_population_count": differential["shared_population_count"],
            "effective_sample_support": differential["effective_sample_support"],
            "population_background_support_json": json.dumps(
                differential["diagnostics"], sort_keys=True, separators=(",", ":"),
            ),
            "test_status": differential["status"], "exclusion_reason": differential["reason"],
        })
        if differential["status"] == "TESTED":
            differential_test = differential_permutation_test(
                x[contrast_mask], y[contrast_mask], local_population[contrast_mask],
                local_background[contrast_mask], reference_background, comparison_background,
                differential["shared_populations"], differential["population_weights"],
                permutations,
                np.random.default_rng(stable_hypothesis_seed(seed, f"differential|{hypothesis_id}")),
            )
            delta = differential["delta_mi"]
            record.update({
                "background_0_conditioned_MI": differential["reference_mi"],
                "background_1_conditioned_MI": differential["comparison_mi"],
                "delta_MI": delta,
                "delta_direction": "COMPARISON_STRONGER" if delta > 0 else "REFERENCE_STRONGER" if delta < 0 else "EQUAL",
                "effective_population_count": differential["effective_population_count"],
                "reference_joint_state_combinations": differential["reference_joint_combinations"],
                "comparison_joint_state_combinations": differential["comparison_joint_combinations"],
                "reference_enriched_A_state": differential["reference_drivers"]["enriched_x_state"],
                "reference_enriched_B_state": differential["reference_drivers"]["enriched_y_state"],
                "reference_depleted_A_state": differential["reference_drivers"]["depleted_x_state"],
                "reference_depleted_B_state": differential["reference_drivers"]["depleted_y_state"],
                "comparison_enriched_A_state": differential["comparison_drivers"]["enriched_x_state"],
                "comparison_enriched_B_state": differential["comparison_drivers"]["enriched_y_state"],
                "comparison_depleted_A_state": differential["comparison_drivers"]["depleted_x_state"],
                "comparison_depleted_B_state": differential["comparison_drivers"]["depleted_y_state"],
                "raw_p": differential_test["p"], "differential_null_mean": differential_test["null_mean"],
            })
            tested.append(base | {
                "test_type": "differential", "metadata_variable": background_column,
                "contrast": f"{comparison_background}-{reference_background}",
                "test_status": "TESTED", "raw_p": differential_test["p"],
                "testing_family": f"differential|{population_column}|{background_column}|{comparison_background}-{reference_background}|{channel}",
            })
        else:
            excluded.append(base | {
                "test_type": "differential", "metadata_variable": background_column,
                "contrast": f"{comparison_background}-{reference_background}",
                "test_status": "NOT_IDENTIFIABLE", "exclusion_reason": differential["reason"],
            })
        records.append(record)
    result_columns = [
        "hypothesis_id", "edge_id", "i", "j", "locus_A", "locus_B", "channel",
        "population_column", "background_variable", "reference_background",
        "comparison_background", "informative_N", "permutations", "global_MI", "global_NMI",
        "global_enriched_A_state", "global_enriched_B_state", "global_depleted_A_state",
        "global_depleted_B_state", "global_joint_state_combinations",
        "population_conditioned_MI", "population_informative_count", "population_raw_p",
        "population_adjusted_q", "population_null_mean", "population_test_status",
        "population_exclusion_reason", "eligible_population_support_json",
        "background_0_conditioned_MI", "background_1_conditioned_MI", "delta_MI",
        "delta_direction", "reference_background_N", "comparison_background_N",
        "reference_population_count", "comparison_population_count", "shared_population_count",
        "effective_sample_support", "effective_population_count",
        "reference_joint_state_combinations", "comparison_joint_state_combinations",
        "reference_enriched_A_state", "reference_enriched_B_state",
        "reference_depleted_A_state", "reference_depleted_B_state",
        "comparison_enriched_A_state", "comparison_enriched_B_state",
        "comparison_depleted_A_state", "comparison_depleted_B_state",
        "population_background_support_json", "raw_p", "adjusted_q",
        "differential_null_mean", "test_status", "exclusion_reason",
    ]
    results = pd.DataFrame(records)
    result_extras = sorted(set(results.columns) - set(result_columns))
    results = results.reindex(columns=result_columns + result_extras)
    tested_columns = [
        "hypothesis_id", "edge_id", "i", "j", "locus_A", "locus_B", "channel",
        "population_column", "background_variable", "reference_background",
        "comparison_background", "informative_N", "permutations", "test_type",
        "metadata_variable", "contrast", "test_status", "raw_p", "adjusted_q",
        "testing_family",
    ]
    excluded_columns = [
        "hypothesis_id", "edge_id", "i", "j", "locus_A", "locus_B", "channel",
        "population_column", "background_variable", "reference_background",
        "comparison_background", "informative_N", "permutations", "test_type",
        "metadata_variable", "contrast", "test_status", "exclusion_reason",
    ]
    tested_frame = pd.DataFrame(tested)
    if len(tested_frame):
        tested_frame["adjusted_q"] = np.nan
        for _, indices in tested_frame.groupby("testing_family", sort=False).groups.items():
            tested_frame.loc[indices, "adjusted_q"] = bh_adjust(
                tested_frame.loc[indices, "raw_p"].to_numpy(float),
            )
        q_lookup = tested_frame.set_index(["hypothesis_id", "test_type"])["adjusted_q"]
        for index, row in results.iterrows():
            population_key = (row.hypothesis_id, "population_conditioned")
            differential_key = (row.hypothesis_id, "differential")
            if population_key in q_lookup:
                results.at[index, "population_adjusted_q"] = q_lookup[population_key]
            if differential_key in q_lookup:
                results.at[index, "adjusted_q"] = q_lookup[differential_key]
    tested_extras = sorted(set(tested_frame.columns) - set(tested_columns))
    tested_frame = tested_frame.reindex(columns=tested_columns + tested_extras)
    excluded_frame = pd.DataFrame(excluded)
    excluded_extras = sorted(set(excluded_frame.columns) - set(excluded_columns))
    excluded_frame = excluded_frame.reindex(columns=excluded_columns + excluded_extras)
    differential_frame = results[results.test_status == "TESTED"].copy()
    _atomic_tsv(results, output / "CONDITIONAL_ASSOCIATIONS.tsv")
    _atomic_tsv(differential_frame, output / "DIFFERENTIAL_ASSOCIATIONS.tsv")
    _atomic_tsv(tested_frame, output / "TESTED_HYPOTHESES.tsv")
    _atomic_tsv(excluded_frame, output / "EXCLUDED_HYPOTHESES.tsv")
    summary = (
        "# CINCH conditional association summary\n\n"
        f"Requested hypotheses: {len(hypotheses)}  \n"
        f"Population tests: {(tested_frame.test_type == 'population_conditioned').sum() if len(tested_frame) else 0}  \n"
        f"Differential tests: {(tested_frame.test_type == 'differential').sum() if len(tested_frame) else 0}  \n"
        f"Non-identifiable/excluded tests: {len(excluded_frame)}  \n\n"
        "The habitat estimand uses common-support population standardization. The differential null "
        "permutes background labels within population and uses a two-sided absolute delta-MI statistic. "
        "Results are context-dependent statistical associations, not proof of epistasis or causality.\n"
    )
    _atomic_text(summary, output / "SUMMARY.md")
    manifest = {
        "schema": "CINCH_CONDITIONAL_ASSOCIATION_V1", "status": "COMPLETE",
        "profile_sha256": _sha256(profile_path), "metadata_sha256": _sha256(metadata_path),
        "pairs_sha256": _sha256(pairs_path) if pairs_path else None,
        "population_column": population_column, "background_column": background_column,
        "background_type": "categorical",
        "reference_background": reference_background, "comparison_background": comparison_background,
        "channels": list(channels), "permutations": permutations, "seed": seed,
        "minimum_eligible_samples": minimum_eligible_samples,
        "minimum_population_cell": minimum_population_cell,
        "minimum_shared_populations": minimum_shared_populations,
        "minimum_background_samples": minimum_background_samples,
        "minimum_state_count": minimum_state_count,
        "background_estimand": "common-support standardized MI; pi_c proportional to min(n_c_reference,n_c_comparison)",
        "population_null": "permute Y within population; one-sided MI statistic",
        "differential_null": "permute background labels within population; two-sided absolute delta_MI statistic",
        "pvalue_convention": "(1 + exceedances) / (B + 1)",
        "requested_hypotheses": len(hypotheses), "tested_hypotheses": len(tested_frame),
        "excluded_hypotheses": len(excluded_frame), "metadata_overlap": overlap_summary,
    }
    _atomic_json(manifest, output / "MANIFEST.json")
    return output
