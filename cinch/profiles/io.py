"""Atomic serialization for validated state profiles."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path

import numpy as np

from .schema import StateProfile


def save_profile(profile: StateProfile, path: Path) -> Path:
    """Save one validated profile without partially replacing an older result."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".npz", dir=path.parent)
    os.close(descriptor)
    temporary_path = Path(temporary)
    try:
        np.savez_compressed(
            temporary_path,
            samples=profile.samples,
            loci=profile.loci,
            presence=profile.presence,
            types=profile.types,
            metadata_json=np.asarray(json.dumps(profile.metadata, sort_keys=True)),
        )
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)
    return path


def load_profile(path: Path) -> StateProfile:
    """Load PROFILE_V2 or frozen_v1 UNIFIED_STATES with validation."""

    with np.load(Path(path), allow_pickle=False) as archive:
        names = set(archive.files)
        required = {"samples", "loci", "presence"}
        if not required.issubset(names):
            raise ValueError(f"profile is missing arrays: {sorted(required - names)}")
        type_key = "types" if "types" in names else "type_state" if "type_state" in names else None
        if type_key is None:
            raise ValueError("profile is missing types/type_state array")
        metadata = {}
        if "metadata_json" in names:
            try:
                metadata = json.loads(str(archive["metadata_json"].item()))
            except (ValueError, TypeError, json.JSONDecodeError) as error:
                raise ValueError("invalid profile metadata_json") from error
        metadata.setdefault("schema", "CINCH_FROZEN_V1_STATE" if type_key == "type_state" else "CINCH_STATE_PROFILE_V2")
        return StateProfile(
            samples=archive["samples"],
            loci=archive["loci"],
            presence=archive["presence"],
            types=archive[type_key],
            metadata=metadata,
        )
