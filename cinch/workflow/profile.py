"""Explicit legacy-profile conversion stage."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

from cinch.profiles import from_legacy_alleles, save_profile


def run_profile_conversion(
    source: Path, output: Path, *, missing_policy: str, separator: str | None = None,
) -> Path:
    source, output = Path(source).resolve(), Path(output).resolve()
    frame = pd.read_csv(source, sep=separator, engine="python" if separator is None else "c")
    if frame.shape[1] < 2:
        raise ValueError("profile table requires sample ID plus at least one locus")
    if frame.iloc[:, 0].duplicated().any():
        raise ValueError("sample identifiers must be unique")
    profile = from_legacy_alleles(
        frame.iloc[:, 1:].to_numpy(object), frame.iloc[:, 0].astype(str).to_numpy(),
        frame.columns[1:].astype(str).to_numpy(), missing_policy=missing_policy,
    )
    output.mkdir(parents=True, exist_ok=True)
    target = save_profile(profile, output / "PROFILE_V2.npz")
    manifest = {
        "schema": "CINCH_PROFILE_CONVERSION_V1",
        "source": str(source),
        "missing_policy": missing_policy,
        "samples": len(profile.samples),
        "loci": len(profile.loci),
        "output": target.name,
    }
    temporary = output / f".MANIFEST.json.{os.getpid()}.tmp"
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, output / "MANIFEST.json")
    return output
