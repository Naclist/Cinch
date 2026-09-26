"""Restartable PP/PT/TP/TT scoring with frozen_v1 scientific semantics."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd

from cinch.frozen_v1.statistics import (
    CHANNELS,
    DISPLAY,
    ORIENTATION,
    channel_vectors,
    contingency,
    pool_rare_states,
)
from cinch.profiles import StateProfile

from .blocks import PairBlock, iter_pair_blocks


def _distance_lookup(distances: pd.DataFrame) -> dict[tuple[int, int], dict]:
    required = {"i", "j"}
    if not required.issubset(distances.columns):
        raise ValueError("distances must contain i and j columns")
    lookup: dict[tuple[int, int], dict] = {}
    for row in distances.to_dict("records"):
        left, right = int(row.pop("i")), int(row.pop("j"))
        if left == right:
            raise ValueError("self-pairs are not valid distances")
        key = (left, right) if left < right else (right, left)
        if key in lookup:
            raise ValueError(f"duplicate distance row for pair {key}")
        lookup[key] = row
    return lookup


def _score_block(
    block: PairBlock,
    presence: np.ndarray,
    pooled_types: np.ndarray,
    loci: np.ndarray,
    distances: dict[tuple[int, int], dict],
    minimum_informative: int,
    minimum_state_count: int,
) -> dict[str, pd.DataFrame]:
    outputs: dict[str, list[dict]] = {channel: [] for channel in CHANNELS}
    for left, right in block.pairs:
        try:
            distance = distances[(left, right)]
        except KeyError as error:
            raise ValueError(f"missing distance row for pair {(left, right)}") from error
        for channel in CHANNELS:
            x, y = channel_vectors(presence, pooled_types, left, right, channel)
            metrics = contingency(x, y)
            if metrics is None:
                continue
            eligible = (
                metrics["informative_N"] >= minimum_informative
                and metrics["n_states_A"] >= 2
                and metrics["n_states_B"] >= 2
                and metrics["min_state_count"] >= minimum_state_count
            )
            if not eligible:
                continue
            outputs[channel].append({
                "edge_id": f"{channel}:{loci[left]}--{loci[right]}",
                "i": left,
                "j": right,
                "locus_A": str(loci[left]),
                "locus_B": str(loci[right]),
                "channel": channel,
                "display_channel": DISPLAY[channel],
                "orientation": ORIENTATION[channel],
                **metrics,
                **distance,
            })
    return {channel: pd.DataFrame(rows) for channel, rows in outputs.items()}


def iter_scored_blocks(
    profile: StateProfile,
    distances: pd.DataFrame,
    *,
    minimum_informative: int,
    minimum_state_count: int,
    block_pairs: int,
) -> tuple[Iterator[tuple[PairBlock, dict[str, pd.DataFrame]]], pd.DataFrame]:
    """Return a lazy block iterator and the once-computed rare-state audit."""

    pooled_types, rare_audit = pool_rare_states(profile.types, minimum_state_count)
    lookup = _distance_lookup(distances)

    def generate() -> Iterator[tuple[PairBlock, dict[str, pd.DataFrame]]]:
        for block in iter_pair_blocks(len(profile.loci), block_pairs):
            yield block, _score_block(
                block, profile.presence, pooled_types, profile.loci, lookup,
                minimum_informative, minimum_state_count,
            )

    return generate(), pd.DataFrame(rare_audit)


def _input_digest(profile: StateProfile, distances: pd.DataFrame, configuration: dict) -> str:
    digest = hashlib.sha256()
    for array in (profile.samples, profile.loci, profile.presence, profile.types):
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode())
        digest.update(str(contiguous.shape).encode())
        digest.update(contiguous.tobytes())
    canonical = distances.sort_values(["i", "j"], kind="mergesort").reset_index(drop=True)
    digest.update(json.dumps(canonical.columns.tolist()).encode())
    digest.update(json.dumps([str(dtype) for dtype in canonical.dtypes]).encode())
    digest.update(pd.util.hash_pandas_object(canonical, index=True).to_numpy().tobytes())
    digest.update(json.dumps(configuration, sort_keys=True).encode())
    return digest.hexdigest()


def _write_parquet_atomic(frame: pd.DataFrame, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".parquet", dir=target.parent)
    os.close(descriptor)
    temporary_path = Path(temporary)
    try:
        frame.to_parquet(temporary_path, index=False)
        os.replace(temporary_path, target)
    finally:
        temporary_path.unlink(missing_ok=True)


def _write_json_atomic(payload: dict, target: Path) -> None:
    temporary = target.with_name(f".{target.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, target)


def write_association_blocks(
    profile: StateProfile,
    distances: pd.DataFrame,
    output_directory: Path,
    *,
    minimum_informative: int,
    minimum_state_count: int,
    block_pairs: int = 100_000,
) -> dict:
    """Write atomic channel blocks and skip only complete compatible blocks."""

    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    configuration = {
        "minimum_informative": minimum_informative,
        "minimum_state_count": minimum_state_count,
        "block_pairs": block_pairs,
    }
    digest = _input_digest(profile, distances, configuration)
    manifest_path = output_directory / "MANIFEST.json"
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if existing.get("input_digest") != digest:
            raise ValueError("existing association output was created from different inputs or configuration")
    else:
        if any((output_directory / "blocks").glob("*/*.parquet")):
            raise ValueError("association blocks exist without a provenance manifest")
        _write_json_atomic({
            "schema": "CINCH_ASSOCIATION_BLOCKS_V1",
            "status": "RUNNING",
            "input_digest": digest,
            "configuration": configuration,
        }, manifest_path)
    pooled_types, rare_rows = pool_rare_states(profile.types, minimum_state_count)
    rare_audit = pd.DataFrame(rare_rows)
    lookup = _distance_lookup(distances)
    _write_parquet_atomic(rare_audit, output_directory / "RARE_STATE_AUDIT.parquet")
    completed, resumed, row_counts = 0, 0, {channel: 0 for channel in CHANNELS}
    for block in iter_pair_blocks(len(profile.loci), block_pairs):
        targets = {
            channel: output_directory / "blocks" / channel / f"block_{block.block_id:08d}.parquet"
            for channel in CHANNELS
        }
        if all(target.exists() for target in targets.values()):
            resumed += 1
            for channel, target in targets.items():
                row_counts[channel] += len(pd.read_parquet(target, columns=[]))
            continue
        frames = _score_block(
            block, profile.presence, pooled_types, profile.loci, lookup,
            minimum_informative, minimum_state_count,
        )
        for channel, frame in frames.items():
            _write_parquet_atomic(frame, targets[channel])
            row_counts[channel] += len(frame)
        completed += 1
    manifest = {
        "schema": "CINCH_ASSOCIATION_BLOCKS_V1",
        "status": "COMPLETE",
        "input_digest": digest,
        "samples": len(profile.samples),
        "loci": len(profile.loci),
        "pair_universe": len(profile.loci) * (len(profile.loci) - 1) // 2,
        "configuration": configuration,
        "completed_blocks_this_run": completed,
        "resumed_blocks_this_run": resumed,
        "eligible_rows": row_counts,
    }
    _write_json_atomic(manifest, manifest_path)
    return manifest
