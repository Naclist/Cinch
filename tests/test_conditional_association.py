import json

import numpy as np
import pandas as pd
import pytest

from cinch.association.high_order import background_stat as unified_legacy_stat
from cinch.cli import main
from cinch.legacy.dev2.shc_diff_gwes.run_diff_gwes import background_stat as frozen_legacy_stat
from cinch.population import (
    differential_permutation_test,
    population_conditioned_mi,
    population_permutation_test,
    standardized_background_mi,
)
from cinch.profiles import StateProfile, save_profile
from cinch.workflow.conditional_association import _channel_states


def _binary_independent(repeats=10):
    cells = np.tile(np.array([[0, 0], [0, 1], [1, 0], [1, 1]]), (repeats, 1))
    return cells[:, 0], cells[:, 1]


def _overlapping_binary(effect="coexistence", populations=3, per_cell=40):
    pop, habitat, x, y = [], [], [], []
    independent_x, independent_y = _binary_independent(per_cell // 4)
    for cluster in range(populations):
        for background in ("soil", "hospital"):
            local_x = np.tile([0, 1], per_cell // 2)
            if background == "soil" or effect == "null":
                local_x, local_y = independent_x.copy(), independent_y.copy()
            elif effect == "coexistence":
                local_y = local_x.copy()
            elif effect == "exclusion":
                local_y = 1 - local_x
            else:
                raise ValueError(effect)
            pop.extend([cluster] * per_cell)
            habitat.extend([background] * per_cell)
            x.extend(local_x)
            y.extend(local_y)
    return np.asarray(x), np.asarray(y), np.asarray(pop), np.asarray(habitat)


def test_case_a_population_only_global_signal_disappears_after_conditioning():
    grid = np.arange(100)
    x0, y0 = (grid // 10 >= 9).astype(int), (grid % 10 >= 9).astype(int)
    x1, y1 = (grid // 10 >= 1).astype(int), (grid % 10 >= 1).astype(int)
    x, y = np.r_[x0, x1], np.r_[y0, y1]
    population = np.repeat([0, 1], 100)
    from cinch.population import categorical_mi
    assert categorical_mi(x, y) > 0.2
    conditioned = population_conditioned_mi(x, y, population)
    assert conditioned["mi"] == pytest.approx(0.0, abs=1e-15)
    result = population_permutation_test(x, y, population, 99, np.random.default_rng(1001))
    assert result["p"] > 0.5


def test_case_b_habitat_specific_coexistence_is_recovered():
    x, y, population, habitat = _overlapping_binary("coexistence")
    result = standardized_background_mi(
        x, y, population, habitat, "soil", "hospital",
        minimum_cell=20, minimum_shared_populations=2, minimum_background_samples=40,
    )
    assert result["status"] == "TESTED"
    assert result["reference_mi"] == pytest.approx(0.0, abs=1e-15)
    assert result["comparison_mi"] == pytest.approx(np.log(2))
    test = differential_permutation_test(
        x, y, population, habitat, "soil", "hospital",
        result["shared_populations"], result["population_weights"], 199,
        np.random.default_rng(2002),
    )
    assert test["observed"] > 0
    assert test["p"] <= 0.01


def test_case_c_habitat_specific_mutual_exclusion_uses_depleted_driver_not_negative_mi():
    x, y, population, habitat = _overlapping_binary("exclusion")
    result = standardized_background_mi(
        x, y, population, habitat, "soil", "hospital",
        minimum_cell=20, minimum_shared_populations=2, minimum_background_samples=40,
    )
    assert result["delta_mi"] > 0
    assert result["comparison_mi"] >= 0
    depleted = result["comparison_drivers"]
    assert depleted["depleted_x_state"] == depleted["depleted_y_state"]


def test_background_driver_never_selects_a_state_absent_from_that_background():
    x, y, population, habitat = _overlapping_binary("coexistence")
    # State 2 occurs only in hospital; soil drivers must remain in its observed {0,1} support.
    x[(habitat == "hospital") & (x == 1)] = 2
    result = standardized_background_mi(
        x, y, population, habitat, "soil", "hospital",
        minimum_cell=20, minimum_shared_populations=2, minimum_background_samples=40,
    )
    assert result["reference_drivers"]["enriched_x_state"] in {0, 1}
    assert result["reference_drivers"]["depleted_x_state"] in {0, 1}


def test_case_d_equal_dependence_has_zero_differential_effect():
    x, y, population, habitat = _overlapping_binary("null")
    # Replace both backgrounds with the same perfect dependence.
    for cluster in np.unique(population):
        for value in ("soil", "hospital"):
            mask = (population == cluster) & (habitat == value)
            x[mask] = np.tile([0, 1], mask.sum() // 2)
            y[mask] = x[mask]
    result = standardized_background_mi(
        x, y, population, habitat, "soil", "hospital",
        minimum_cell=20, minimum_shared_populations=2, minimum_background_samples=40,
    )
    assert result["delta_mi"] == pytest.approx(0.0)
    test = differential_permutation_test(
        x, y, population, habitat, "soil", "hospital",
        result["shared_populations"], result["population_weights"], 99,
        np.random.default_rng(3003),
    )
    assert test["p"] > 0.5


def test_case_e_perfect_habitat_lineage_confounding_is_not_identifiable():
    x = np.tile([0, 1], 40)
    y = x.copy()
    population = np.repeat([0, 1], 40)
    habitat = np.repeat(["soil", "hospital"], 40)
    result = standardized_background_mi(
        x, y, population, habitat, "soil", "hospital",
        minimum_cell=10, minimum_shared_populations=1, minimum_background_samples=20,
    )
    assert result["status"] == "NOT_IDENTIFIABLE"
    assert result["reason"] == "INSUFFICIENT_SHARED_POPULATIONS"
    assert "delta_mi" not in result


def test_case_f_multiallelic_tt_habitat_dependence():
    x, y, population, habitat = [], [], [], []
    independent = np.tile(np.array([(a, b) for a in range(3) for b in range(3)]), (4, 1))
    dependent = np.tile(np.array([[0, 0], [1, 1], [2, 2]]), (12, 1))
    for cluster in range(3):
        for value, cells in (("soil", independent), ("hospital", dependent)):
            population.extend([cluster] * len(cells)); habitat.extend([value] * len(cells))
            x.extend(cells[:, 0]); y.extend(cells[:, 1])
    result = standardized_background_mi(
        np.asarray(x), np.asarray(y), np.asarray(population), np.asarray(habitat),
        "soil", "hospital", minimum_cell=20, minimum_shared_populations=2,
        minimum_background_samples=40,
    )
    assert result["status"] == "TESTED"
    assert result["reference_mi"] == pytest.approx(0.0, abs=1e-15)
    assert result["comparison_mi"] == pytest.approx(np.log(3))


def test_monomorphic_background_margin_remains_estimable_on_sample_common_support():
    population = np.repeat([0, 1], 40)
    habitat = np.tile(np.repeat(["soil", "hospital"], 20), 2)
    x = np.zeros(80, dtype=int)
    y = np.zeros(80, dtype=int)
    hospital = habitat == "hospital"
    x[hospital] = np.tile([0, 1], 20)
    y[hospital] = x[hospital]
    result = standardized_background_mi(
        x, y, population, habitat, "soil", "hospital",
        minimum_cell=10, minimum_shared_populations=2, minimum_background_samples=20,
    )
    assert result["status"] == "TESTED"
    assert result["shared_population_count"] == 2
    assert result["reference_mi"] == pytest.approx(0.0)
    assert result["comparison_mi"] == pytest.approx(np.log(2))
    assert all(not item["state_variable_in_both_backgrounds"] for item in result["diagnostics"])


def test_case_i_unbalanced_habitats_use_common_support_weights():
    x, y, population, habitat = [], [], [], []
    for cluster, sizes in ((0, (80, 20)), (1, (20, 80))):
        for value, size in zip(("soil", "hospital"), sizes):
            local_x = np.tile([0, 1], size // 2)
            local_y = np.resize(np.array([0, 0, 1, 1]), size) if value == "soil" else local_x
            population.extend([cluster] * size); habitat.extend([value] * size)
            x.extend(local_x); y.extend(local_y)
    result = standardized_background_mi(
        np.asarray(x), np.asarray(y), np.asarray(population), np.asarray(habitat),
        "soil", "hospital", minimum_cell=10, minimum_shared_populations=2,
        minimum_background_samples=40,
    )
    assert result["effective_sample_support"] == 40
    assert result["population_weights"].tolist() == pytest.approx([0.5, 0.5])
    assert result["delta_mi"] > 0
    test = differential_permutation_test(
        np.asarray(x), np.asarray(y), np.asarray(population), np.asarray(habitat),
        "soil", "hospital", result["shared_populations"], result["population_weights"],
        99, np.random.default_rng(6006),
    )
    assert test["p"] <= 0.05
    repeated = differential_permutation_test(
        np.asarray(x), np.asarray(y), np.asarray(population), np.asarray(habitat),
        "soil", "hospital", result["shared_populations"], result["population_weights"],
        99, np.random.default_rng(6006),
    )
    assert repeated == test


def _all_channel_profile():
    samples, populations, habitats, presence_rows, type_rows = [], [], [], [], []
    loci = np.array(["PPA", "PPB", "P", "T", "TA", "Q", "T1", "T2"])
    binary_independent = np.tile(np.array([[0, 0], [0, 1], [1, 0], [1, 1]]), (9, 1))
    ternary_independent = np.tile(np.array([(a, b) for a in range(3) for b in range(3)]), (4, 1))
    for cluster in range(4):
        for habitat in ("soil", "hospital"):
            binary = binary_independent if habitat == "soil" else np.tile([[0, 0], [1, 1]], (18, 1))
            ternary = ternary_independent if habitat == "soil" else np.tile([[0, 0], [1, 1], [2, 2]], (12, 1))
            for index in range(36):
                bx, by = binary[index]; tx, ty = ternary[index]
                p = np.ones(8, dtype=np.int8)
                t = np.zeros(8, dtype=np.int32)
                p[0], p[1], p[2], p[5] = bx, by, bx, by
                t[p == 0] = -1
                t[3], t[4], t[6], t[7] = by, bx, tx, ty
                samples.append(f"s{len(samples)}"); populations.append(f"c{cluster}"); habitats.append(habitat)
                presence_rows.append(p); type_rows.append(t)
    # Explicit unresolved/multicopy-ambiguous rows must be excluded from every channel.
    for cluster in range(4):
        for habitat in ("soil", "hospital"):
            samples.append(f"amb{cluster}_{habitat}"); populations.append(f"c{cluster}"); habitats.append(habitat)
            presence_rows.append(np.full(8, -1, dtype=np.int8)); type_rows.append(np.full(8, -1, dtype=np.int32))
    profile = StateProfile(np.asarray(samples), loci, np.asarray(presence_rows), np.asarray(type_rows))
    metadata = pd.DataFrame({"sample_id": samples, "tree_shc": populations, "habitat": habitats})
    pairs = pd.DataFrame({
        "locus_A": ["PPA", "P", "TA", "T1"],
        "locus_B": ["PPB", "T", "Q", "T2"],
        "channel": ["PP", "PT", "TP", "TT"],
    })
    return profile, metadata, pairs


def test_cases_g_h_and_cli_end_to_end_all_channels(tmp_path):
    profile, metadata, pairs = _all_channel_profile()
    profile_path, metadata_path, pairs_path = tmp_path / "profile.npz", tmp_path / "metadata.tsv", tmp_path / "pairs.tsv"
    save_profile(profile, profile_path)
    metadata.to_csv(metadata_path, sep="\t", index=False)
    pairs.to_csv(pairs_path, sep="\t", index=False)
    output = tmp_path / "conditional"
    assert main([
        "conditional-associate", "--profile", str(profile_path), "--metadata", str(metadata_path),
        "--population-column", "tree_shc", "--background-column", "habitat",
        "--reference-background", "soil", "--comparison-background", "hospital",
        "--channels", "PP", "PT", "TP", "TT", "--pairs", str(pairs_path),
        "--permutations", "99", "--seed", "4004", "--min-eligible-samples", "100",
        "--min-population-cell", "20", "--min-shared-populations", "3",
        "--min-background-samples", "100", "-o", str(output),
    ]) == 0
    required = {
        "MANIFEST.json", "METADATA_QC.tsv", "POPULATION_BACKGROUND_OVERLAP.tsv",
        "CONDITIONAL_ASSOCIATIONS.tsv", "DIFFERENTIAL_ASSOCIATIONS.tsv",
        "TESTED_HYPOTHESES.tsv", "EXCLUDED_HYPOTHESES.tsv", "SUMMARY.md",
    }
    assert required == {path.name for path in output.iterdir()}
    results = pd.read_csv(output / "CONDITIONAL_ASSOCIATIONS.tsv", sep="\t")
    assert {
        "locus_A", "locus_B", "channel", "background_variable",
        "reference_background", "comparison_background", "global_MI",
        "population_conditioned_MI", "background_0_conditioned_MI",
        "background_1_conditioned_MI", "delta_MI", "delta_direction",
        "informative_N", "shared_population_count", "raw_p", "adjusted_q",
        "test_status",
    }.issubset(results.columns)
    assert set(results.channel) == {"PP", "PT", "TP", "TT"}
    assert results.test_status.eq("TESTED").all()
    assert results.population_test_status.eq("TESTED").all()
    assert results.informative_N.eq(288).all()
    assert results.delta_MI.gt(0.5).all()
    assert results.adjusted_q.le(0.05).all()
    tested = pd.read_csv(output / "TESTED_HYPOTHESES.tsv", sep="\t")
    assert len(tested) == 8
    assert tested.adjusted_q.notna().all()
    assert pd.read_csv(output / "EXCLUDED_HYPOTHESES.tsv", sep="\t").empty
    manifest = json.loads((output / "MANIFEST.json").read_text())
    assert manifest["status"] == "COMPLETE"
    assert manifest["tested_hypotheses"] == 8


def test_case_h_callable_presence_with_ambiguous_type_is_not_an_allele_category():
    profile = StateProfile(
        np.array(["s0", "s1"]), np.array(["A", "B"]),
        np.ones((2, 2), dtype=np.int8), np.full((2, 2), -1, dtype=np.int32),
    )
    pooled = profile.types.copy()
    assert _channel_states(profile, pooled, 0, 1, "PP")[2].tolist() == [True, True]
    assert not _channel_states(profile, pooled, 0, 1, "PT")[2].any()
    assert not _channel_states(profile, pooled, 0, 1, "TP")[2].any()
    assert not _channel_states(profile, pooled, 0, 1, "TT")[2].any()


def test_metadata_mismatch_is_reported_before_inference(tmp_path):
    profile, metadata, pairs = _all_channel_profile()
    profile_path, metadata_path = tmp_path / "profile.npz", tmp_path / "metadata.tsv"
    save_profile(profile, profile_path)
    metadata.iloc[:-1].to_csv(metadata_path, sep="\t", index=False)
    output = tmp_path / "failed"
    assert main([
        "conditional-associate", "--profile", str(profile_path), "--metadata", str(metadata_path),
        "--population-column", "tree_shc", "--background-column", "habitat",
        "--reference-background", "soil", "--comparison-background", "hospital",
        "--permutations", "9", "-o", str(output),
    ]) == 2
    qc = pd.read_csv(output / "METADATA_QC.tsv", sep="\t")
    assert qc.loc[qc.check == "profile_samples_missing_metadata", "count"].iloc[0] == 1
    assert json.loads((output / "MANIFEST.json").read_text())["status"] == "FAILED_METADATA_QC"


def test_case_e_cli_writes_not_identifiable_without_fake_effect_or_pvalue(tmp_path):
    x = np.tile([0, 1], 40).astype(np.int8)
    presence = np.column_stack([x, x])
    types = np.where(presence == 1, 0, -1).astype(np.int32)
    samples = np.array([f"s{i}" for i in range(80)])
    profile = StateProfile(samples, np.array(["A", "B"]), presence, types)
    profile_path = tmp_path / "profile.npz"
    save_profile(profile, profile_path)
    metadata = pd.DataFrame({
        "sample_id": samples, "tree_shc": np.repeat(["c0", "c1"], 40),
        "habitat": np.repeat(["soil", "hospital"], 40),
    })
    metadata_path = tmp_path / "metadata.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)
    output = tmp_path / "conditional"
    assert main([
        "conditional-associate", "--profile", str(profile_path), "--metadata", str(metadata_path),
        "--population-column", "tree_shc", "--background-column", "habitat",
        "--reference-background", "soil", "--comparison-background", "hospital",
        "--channels", "PP", "--permutations", "9", "--min-eligible-samples", "20",
        "--min-population-cell", "10", "--min-shared-populations", "1",
        "--min-background-samples", "20", "-o", str(output),
    ]) == 0
    row = pd.read_csv(output / "CONDITIONAL_ASSOCIATIONS.tsv", sep="\t").iloc[0]
    assert row.test_status == "NOT_IDENTIFIABLE"
    assert row.exclusion_reason == "INSUFFICIENT_SHARED_POPULATIONS"
    assert pd.isna(row.delta_MI) and pd.isna(row.raw_p) and pd.isna(row.adjusted_q)
    excluded = pd.read_csv(output / "EXCLUDED_HYPOTHESES.tsv", sep="\t")
    assert ((excluded.test_type == "differential") &
            (excluded.test_status == "NOT_IDENTIFIABLE")).any()


def test_all_excluded_run_keeps_stable_audit_schemas(tmp_path):
    samples = np.array([f"s{i}" for i in range(40)])
    profile = StateProfile(
        samples, np.array(["A", "B"]), np.ones((40, 2), dtype=np.int8),
        np.zeros((40, 2), dtype=np.int32),
    )
    profile_path = tmp_path / "profile.npz"
    save_profile(profile, profile_path)
    metadata = pd.DataFrame({
        "sample_id": samples, "tree_shc": np.repeat(["c0", "c1"], 20),
        "habitat": np.tile(np.repeat(["soil", "hospital"], 10), 2),
    })
    metadata_path = tmp_path / "metadata.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)
    output = tmp_path / "conditional"
    assert main([
        "conditional-associate", "--profile", str(profile_path), "--metadata", str(metadata_path),
        "--population-column", "tree_shc", "--background-column", "habitat",
        "--reference-background", "soil", "--comparison-background", "hospital",
        "--channels", "PP", "--permutations", "9", "--min-eligible-samples", "20",
        "--min-population-cell", "5", "--min-shared-populations", "2",
        "--min-background-samples", "10", "-o", str(output),
    ]) == 0
    tested = pd.read_csv(output / "TESTED_HYPOTHESES.tsv", sep="\t")
    excluded = pd.read_csv(output / "EXCLUDED_HYPOTHESES.tsv", sep="\t")
    assert tested.empty
    assert {"channel", "metadata_variable", "contrast", "raw_p", "adjusted_q",
            "testing_family"}.issubset(tested.columns)
    assert len(excluded) == 2
    assert {"channel", "metadata_variable", "contrast", "exclusion_reason"}.issubset(excluded.columns)


def test_case_j_frozen_diff_gwes_statistic_is_exactly_preserved():
    rng = np.random.default_rng(5005)
    x = rng.integers(0, 2, 120).astype(bool)
    y = rng.integers(0, 2, 120).astype(bool)
    c = rng.integers(0, 2, 120).astype(bool)
    weights = rng.random(120) + 0.1
    shc = np.repeat(np.arange(6), 20)
    expected = frozen_legacy_stat(x, y, c, weights, shc, min_group=20, min_shc_cell=3, min_informative=2)
    observed = unified_legacy_stat(x, y, c, weights, shc, min_group=20, min_shc_cell=3, min_informative=2)
    assert observed == expected
