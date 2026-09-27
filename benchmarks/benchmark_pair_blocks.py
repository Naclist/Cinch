#!/usr/bin/env python3
"""Compare frozen all-result accumulation with bounded block consumption."""

from __future__ import annotations

import argparse
import itertools
import json
import time
import tracemalloc

import numpy as np
import pandas as pd

from cinch.association import iter_scored_blocks
from cinch.frozen_v1.statistics import CHANNELS, score_all_pairs
from cinch.profiles import StateProfile


def measured(function):
    tracemalloc.start()
    started = time.perf_counter()
    value = function()
    elapsed = time.perf_counter() - started
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return value, elapsed, peak


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=224)
    parser.add_argument("--loci", type=int, default=80)
    parser.add_argument("--block-pairs", type=int, default=256)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    presence = rng.choice(np.array([-1, 0, 1], dtype=np.int8), (args.samples, args.loci), p=[0.03, 0.17, 0.80])
    types = rng.integers(0, 8, size=presence.shape, dtype=np.int32)
    types[presence != 1] = -1
    profile = StateProfile(
        np.array([f"s{i}" for i in range(args.samples)]),
        np.array([f"L{i}" for i in range(args.loci)]), presence, types,
    )
    distances = pd.DataFrame([
        {"i": i, "j": j, "order_distance": float(j - i), "bp_distance": float((j - i) * 1000)}
        for i, j in itertools.combinations(range(args.loci), 2)
    ])

    def frozen():
        result, _ = score_all_pairs(presence, types, profile.loci, distances, 20, 2)
        return {channel: len(result[channel]) for channel in CHANNELS}

    def blocked():
        iterator, _ = iter_scored_blocks(
            profile, distances, minimum_informative=20,
            minimum_state_count=2, block_pairs=args.block_pairs,
        )
        counts = {channel: 0 for channel in CHANNELS}
        for _, frames in iterator:
            for channel in CHANNELS:
                counts[channel] += len(frames[channel])
        return counts

    frozen_counts, frozen_seconds, frozen_peak = measured(frozen)
    block_counts, block_seconds, block_peak = measured(blocked)
    print(json.dumps({
        "scope": "controlled computational benchmark only; excludes block I/O",
        "samples": args.samples,
        "loci": args.loci,
        "pairs": args.loci * (args.loci - 1) // 2,
        "block_pairs": args.block_pairs,
        "frozen": {"seconds": frozen_seconds, "tracemalloc_peak_bytes": frozen_peak, "rows": frozen_counts},
        "blockwise": {"seconds": block_seconds, "tracemalloc_peak_bytes": block_peak, "rows": block_counts},
        "row_counts_equal": frozen_counts == block_counts,
        "peak_python_memory_ratio": block_peak / frozen_peak,
        "runtime_ratio": block_seconds / frozen_seconds,
    }, indent=2))


if __name__ == "__main__":
    main()
