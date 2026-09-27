import json

import numpy as np
import pandas as pd

from cinch.cli import main
from cinch.profiles import load_profile


def test_profile_and_associate_cli_stages(tmp_path, capsys):
    source = tmp_path / "legacy.tsv"
    pd.DataFrame({
        "sample": [f"s{i}" for i in range(8)],
        "A": [1, 1, 2, 2, 1, 1, 2, 2],
        "B": [1, 1, 2, 2, 1, 1, 2, 2],
        "C": [0, 1, 0, 1, 0, 1, 0, 1],
    }).to_csv(source, sep="\t", index=False)
    profile_dir = tmp_path / "profile"
    assert main(["profile", str(source), "-o", str(profile_dir), "--missing-policy", "absence"]) == 0
    profile = load_profile(profile_dir / "PROFILE_V2.npz")
    assert profile.presence[:, 2].tolist() == [0, 1, 0, 1, 0, 1, 0, 1]

    coordinates = pd.DataFrame([
        {"sample_id": sample, "contig": "c", "locus_index": locus,
         "start": 1 + locus * 100, "end": 50 + locus * 100, "gene_order": locus + 1}
        for sample in profile.samples for locus in range(len(profile.loci))
        if profile.presence[list(profile.samples).index(sample), locus] == 1
    ])
    coordinate_path = tmp_path / "coordinates.tsv"
    coordinates.to_csv(coordinate_path, sep="\t", index=False)
    association_dir = tmp_path / "associate"
    assert main([
        "associate", "--profile", str(profile_dir / "PROFILE_V2.npz"),
        "--coordinates", str(coordinate_path), "-o", str(association_dir),
        "--min-informative", "4", "--min-state-count", "1",
        "--min-order-observations", "2", "--block-pairs", "2",
    ]) == 0
    manifest = json.loads((association_dir / "association" / "MANIFEST.json").read_text())
    assert manifest["status"] == "COMPLETE"
    assert manifest["pair_universe"] == 3
    # A second run resumes both triangular blocks.
    assert main([
        "associate", "--profile", str(profile_dir / "PROFILE_V2.npz"),
        "--coordinates", str(coordinate_path), "-o", str(association_dir),
        "--min-informative", "4", "--min-state-count", "1",
        "--min-order-observations", "2", "--block-pairs", "2",
    ]) == 0
    resumed = json.loads((association_dir / "association" / "MANIFEST.json").read_text())
    assert resumed["completed_blocks_this_run"] == 0
    assert resumed["resumed_blocks_this_run"] == 2
