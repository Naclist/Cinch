import json

import numpy as np
import pandas as pd

from cinch.profiles import StateProfile, save_profile
from cinch.workflow.advanced_filter import run_advanced_filter
from cinch.workflow.reporting import run_report


def test_staged_shc_bh_aracne_and_report(tmp_path):
    samples = np.array([f"s{i}" for i in range(16)])
    loci = np.array(["A", "B", "C"])
    pattern = np.tile([0, 0, 1, 1], 4)
    presence = np.column_stack([pattern, pattern, np.roll(pattern, 1)]).astype(np.int8)
    types = np.full(presence.shape, -1, dtype=np.int32)
    types[presence == 1] = 0
    profile_path = save_profile(StateProfile(samples, loci, presence, types), tmp_path / "profile.npz")
    association = tmp_path / "association"
    (association / "blocks" / "PP").mkdir(parents=True)
    (association / "MANIFEST.json").write_text(json.dumps({"schema": "test"}) + "\n")
    pd.DataFrame([
        {"edge_id": "PP:A--B", "i": 0, "j": 1, "locus_A": "A", "locus_B": "B", "channel": "PP", "order_distance": 10.0},
        {"edge_id": "PP:A--C", "i": 0, "j": 2, "locus_A": "A", "locus_B": "C", "channel": "PP", "order_distance": 20.0},
    ]).to_parquet(association / "blocks" / "PP" / "block_00000000.parquet", index=False)
    shc = tmp_path / "shc.tsv"
    pd.DataFrame({"sample_id": samples, "shc": np.repeat([0, 1], 8)}).to_csv(shc, sep="\t", index=False)
    filtered = run_advanced_filter(
        association, profile_path, shc, tmp_path / "filtered",
        permutations=39, minimum_cluster_size=4,
        minimum_informative_clusters=2, alpha=1.0, seed=7,
    )
    result = pd.read_parquet(filtered / "SHC_ALL_PP_CANDIDATES.parquet")
    assert result.eligible_for_permutation.tolist() == [True, False]
    assert result.loc[result.eligible_for_permutation, "phylo_perm_q"].notna().all()
    assert (filtered / "SHC_SIGNIFICANT_ARACNE.parquet").is_file()
    manifest = json.loads((filtered / "MANIFEST.json").read_text())
    assert manifest["scientific_scope"] == "PP binary presence edges only"
    report = run_report(filtered, tmp_path / "report")
    assert (report / "REPORT.md").is_file()
    assert (report / "SUMMARY.tsv").is_file()
    assert (report / "SHC_FILTER_DIAGNOSTICS.png").is_file()


def test_advanced_filter_accepts_empty_pp_candidate_blocks(tmp_path):
    samples = np.array(["s1", "s2"])
    profile_path = save_profile(
        StateProfile(samples, np.array(["A"]), np.ones((2, 1), np.int8), np.zeros((2, 1), np.int32)),
        tmp_path / "profile.npz",
    )
    association = tmp_path / "association"
    (association / "blocks" / "PP").mkdir(parents=True)
    (association / "MANIFEST.json").write_text('{"schema":"test"}\n')
    pd.DataFrame(columns=["edge_id", "i", "j", "locus_A", "locus_B", "channel"]).to_parquet(
        association / "blocks" / "PP" / "block_00000000.parquet", index=False,
    )
    shc = tmp_path / "shc.tsv"
    pd.DataFrame({"sample_id": samples, "shc": [0, 0]}).to_csv(shc, sep="\t", index=False)
    output = run_advanced_filter(
        association, profile_path, shc, tmp_path / "filtered",
        permutations=9, minimum_cluster_size=2, minimum_informative_clusters=1,
    )
    result = pd.read_parquet(output / "SHC_ALL_PP_CANDIDATES.parquet")
    assert result.empty
    assert {"eligible_for_permutation", "phylo_perm_q", "passes_shc_bh"}.issubset(result.columns)
