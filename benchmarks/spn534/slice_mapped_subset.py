#!/usr/bin/env python3
"""Create a nested sample subset from one real-genome mapping result."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from cinch.profiles import StateProfile, load_profile, save_profile


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping-result", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    profile = load_profile(args.mapping_result / "profiles" / "PROFILE_V2.npz")
    manifest = pd.read_csv(args.manifest, sep="\t")
    requested = set(manifest["sample_id"].astype(str) + "_genomic")
    mask = [sample in requested for sample in profile.samples]
    observed = set(profile.samples[mask])
    if observed != requested:
        raise ValueError(f"sample mismatch: missing={sorted(requested - observed)}")
    metadata = dict(profile.metadata)
    metadata.update({"benchmark_subset": True, "sample_manifest": str(args.manifest.resolve())})
    subset = StateProfile(
        profile.samples[mask], profile.loci, profile.presence[mask], profile.types[mask], metadata,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    save_profile(subset, args.output / "PROFILE_V2.npz")
    coordinates = pd.read_csv(args.mapping_result / "mapping" / "COORDINATES.tsv", sep="\t")
    coordinates = coordinates[coordinates["sample_id"].isin(requested)].copy()
    if set(coordinates["sample_id"]) != requested:
        raise ValueError("coordinates do not cover every requested sample")
    coordinates.to_csv(args.output / "COORDINATES.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
