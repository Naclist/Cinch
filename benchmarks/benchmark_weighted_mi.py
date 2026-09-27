#!/usr/bin/env python3
"""Reproducible microbenchmark for the unified weighted binary MI kernel."""

from __future__ import annotations

import argparse
import json
import platform
import time

import numpy as np

from cinch.statistics import weighted_binary_information


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=int, default=1_000_000)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    mx = rng.uniform(0.05, 0.95, args.pairs)
    my = rng.uniform(0.05, 0.95, args.pairs)
    low = np.maximum(0.0, mx + my - 1.0)
    high = np.minimum(mx, my)
    c11 = rng.uniform(low, high)
    weighted_binary_information(c11[:10], mx[:10], my[:10])
    times = []
    for _ in range(args.repeats):
        started = time.perf_counter()
        result = weighted_binary_information(c11, mx, my)
        times.append(time.perf_counter() - started)
    best = min(times)
    print(json.dumps({
        "pairs": args.pairs,
        "repeats": args.repeats,
        "seed": args.seed,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "seconds": times,
        "best_seconds": best,
        "pairs_per_second_best": args.pairs / best,
        "finite_mi": int(np.isfinite(result[0]).sum()),
    }, indent=2))


if __name__ == "__main__":
    main()
