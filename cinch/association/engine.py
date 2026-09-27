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
from .numba_kernel import METRIC_NAMES, available as numba_available, score_block_metrics


class _FrameDistanceLookup:
    def __init__(self, distances: pd.DataFrame):
        self.values = _distance_lookup(distances)

    def get_pair(self, left: int, right: int) -> dict:
        try:
            return self.values[(left, right)]
        except KeyError as error:
            raise ValueError(f"missing distance row for pair {(left, right)}") from error


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


def _distance_source(distances):
    if hasattr(distances, "get_pair"):
        return distances
    if isinstance(distances, pd.DataFrame):
        return _FrameDistanceLookup(distances)
    raise TypeError("distances must be a DataFrame or expose get_pair(i, j)")


def _score_block(
    block: PairBlock,
    presence: np.ndarray,
    pooled_types: np.ndarray,
    loci: np.ndarray,
    distances,
    presence_variable: np.ndarray,
    type_variable: np.ndarray,
    minimum_informative: int,
    minimum_state_count: int,
) -> dict[str, pd.DataFrame]:
    outputs: dict[str, list[dict]] = {channel: [] for channel in CHANNELS}
    for left, right in block.pairs:
        distance = None
        for channel in CHANNELS:
            if channel == "PP" and not (presence_variable[left] and presence_variable[right]):
                continue
            if channel == "PT" and not (presence_variable[left] and type_variable[right]):
                continue
            if channel == "TP" and not (type_variable[left] and presence_variable[right]):
                continue
            if channel == "TT" and not (type_variable[left] and type_variable[right]):
                continue
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
            if distance is None:
                distance = distances.get_pair(left, right)
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


def _score_block_numba(
    block: PairBlock, presence: np.ndarray, pooled_types: np.ndarray, loci: np.ndarray,
    distances, presence_variable: np.ndarray, type_variable: np.ndarray,
    minimum_informative: int, minimum_state_count: int, type_cardinality: np.ndarray,
) -> dict[str, pd.DataFrame]:
    if not numba_available():
        raise RuntimeError("Numba engine requires `pip install -e '.[performance]'`")
    pairs = np.asarray(block.pairs, dtype=np.int32)
    scored = score_block_metrics(
        presence, pooled_types, pairs, presence_variable, type_variable,
        minimum_informative, minimum_state_count, type_cardinality,
    )
    integer_metrics = {
        "informative_N", "n_states_A", "n_states_B", "min_state_count",
        "enriched_driver_A_state", "enriched_driver_B_state", "enriched_observed",
        "depleted_driver_A_state", "depleted_driver_B_state", "depleted_observed",
    }
    outputs: dict[str, list[dict]] = {channel: [] for channel in CHANNELS}
    for pair_index, (left, right) in enumerate(block.pairs):
        distance = None
        for channel_index, channel in enumerate(CHANNELS):
            if scored[pair_index, channel_index, 0] != 1.0:
                continue
            if distance is None:
                distance = distances.get_pair(left, right)
            metrics = {
                name: int(scored[pair_index, channel_index, index + 1]) if name in integer_metrics
                else float(scored[pair_index, channel_index, index + 1])
                for index, name in enumerate(METRIC_NAMES)
            }
            outputs[channel].append({
                "edge_id": f"{channel}:{loci[left]}--{loci[right]}", "i": left, "j": right,
                "locus_A": str(loci[left]), "locus_B": str(loci[right]), "channel": channel,
                "display_channel": DISPLAY[channel], "orientation": ORIENTATION[channel],
                **metrics, **distance,
            })
    return {channel: pd.DataFrame(rows) for channel, rows in outputs.items()}


def iter_scored_blocks(
    profile: StateProfile,
    distances: pd.DataFrame,
    *,
    minimum_informative: int,
    minimum_state_count: int,
    block_pairs: int,
    engine: str = "python",
) -> tuple[Iterator[tuple[PairBlock, dict[str, pd.DataFrame]]], pd.DataFrame]:
    """Return a lazy block iterator and the once-computed rare-state audit."""

    if engine not in {"python", "numba"}:
        raise ValueError("engine must be 'python' or 'numba'")

    pooled_types, rare_audit = pool_rare_states(profile.types, minimum_state_count)
    presence_variable = np.array([
        len(np.unique(column[column >= 0])) >= 2 for column in profile.presence.T
    ], dtype=bool)
    type_variable = np.array([
        len(np.unique(column[column >= 0])) >= 2 for column in pooled_types.T
    ], dtype=bool)
    type_cardinality = np.array([
        int(column[column >= 0].max()) + 1 if np.any(column >= 0) else 1
        for column in pooled_types.T
    ], dtype=np.int32)
    lookup = _distance_source(distances)

    def generate() -> Iterator[tuple[PairBlock, dict[str, pd.DataFrame]]]:
        for block in iter_pair_blocks(len(profile.loci), block_pairs):
            arguments = (
                block, profile.presence, pooled_types, profile.loci, lookup,
                presence_variable, type_variable, minimum_informative, minimum_state_count,
            )
            frames = (
                _score_block_numba(*arguments, type_cardinality)
                if engine == "numba" else _score_block(*arguments)
            )
            yield block, frames

    return generate(), pd.DataFrame(rare_audit)


def _input_digest(profile: StateProfile, distances, configuration: dict) -> str:
    digest = hashlib.sha256()
    for array in (profile.samples, profile.loci, profile.presence, profile.types):
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode())
        digest.update(str(contiguous.shape).encode())
        digest.update(contiguous.tobytes())
    if hasattr(distances, "digest"):
        digest.update(str(distances.digest()).encode())
    elif isinstance(distances, pd.DataFrame):
        canonical = distances.sort_values(["i", "j"], kind="mergesort").reset_index(drop=True)
        digest.update(json.dumps(canonical.columns.tolist()).encode())
        digest.update(json.dumps([str(dtype) for dtype in canonical.dtypes]).encode())
        digest.update(pd.util.hash_pandas_object(canonical, index=True).to_numpy().tobytes())
    else:
        raise TypeError("distance source requires a stable digest")
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
    engine: str = "python",
) -> dict:
    """Write atomic channel blocks and skip only complete compatible blocks."""

    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    if engine not in {"python", "numba"}:
        raise ValueError("engine must be 'python' or 'numba'")
    configuration = {
        "minimum_informative": minimum_informative,
        "minimum_state_count": minimum_state_count,
        "block_pairs": block_pairs,
        "engine": engine,
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
    presence_variable = np.array([
        len(np.unique(column[column >= 0])) >= 2 for column in profile.presence.T
    ], dtype=bool)
    type_variable = np.array([
        len(np.unique(column[column >= 0])) >= 2 for column in pooled_types.T
    ], dtype=bool)
    type_cardinality = np.array([
        int(column[column >= 0].max()) + 1 if np.any(column >= 0) else 1
        for column in pooled_types.T
    ], dtype=np.int32)
    rare_audit = pd.DataFrame(rare_rows)
    lookup = _distance_source(distances)
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
        arguments = (
            block, profile.presence, pooled_types, profile.loci, lookup,
            presence_variable, type_variable, minimum_informative, minimum_state_count,
        )
        frames = (
            _score_block_numba(*arguments, type_cardinality)
            if engine == "numba" else _score_block(*arguments)
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
