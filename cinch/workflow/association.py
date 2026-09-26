"""Staged distance and four-channel association workflow."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

from cinch.association import write_association_blocks
from cinch.frozen_v1.mapping import pair_order_distance
from cinch.profiles import load_profile


def _read_coordinates(path: Path) -> pd.DataFrame:
    if path.suffix.casefold() in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t")


def run_association(
    profile_path: Path,
    coordinates_path: Path,
    output: Path,
    *,
    minimum_informative: int,
    minimum_state_count: int,
    minimum_order_observations: int,
    block_pairs: int,
) -> Path:
    profile_path, coordinates_path = Path(profile_path).resolve(), Path(coordinates_path).resolve()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    profile = load_profile(profile_path)
    coordinates = _read_coordinates(coordinates_path)
    distances = pair_order_distance(coordinates, len(profile.loci), minimum_order_observations)
    distances.insert(2, "locus_A", profile.loci[distances.i.to_numpy(int)])
    distances.insert(3, "locus_B", profile.loci[distances.j.to_numpy(int)])
    distance_path = output / "DISTANCES.parquet"
    distances.to_parquet(distance_path, index=False)
    manifest = write_association_blocks(
        profile, distances, output / "association",
        minimum_informative=minimum_informative,
        minimum_state_count=minimum_state_count,
        block_pairs=block_pairs,
    )
    stage = {
        "schema": "CINCH_ASSOCIATION_STAGE_V1",
        "profile": str(profile_path),
        "coordinates": str(coordinates_path),
        "distance_output": distance_path.name,
        "association_output": "association",
        "association_manifest": manifest,
    }
    temporary = output / f".MANIFEST.json.{os.getpid()}.tmp"
    temporary.write_text(json.dumps(stage, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, output / "MANIFEST.json")
    return output
