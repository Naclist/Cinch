#!/usr/bin/env python3
"""Controlled timing for SHC permutation and sparse ARACNE stages."""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd

from cinch.network import aracne_prune
from cinch.population import phylo_permutation


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--clusters", type=int, default=20)
    parser.add_argument("--permutations", type=int, default=199)
    parser.add_argument("--nodes", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    clusters = np.arange(args.samples) % args.clusters
    x = rng.integers(0, 2, args.samples).astype(bool)
    y = np.logical_xor(x, rng.random(args.samples) < 0.2)
    weights = rng.random(args.samples) + 0.01
    eligible = list(range(args.clusters))
    started = time.perf_counter()
    statistic = phylo_permutation(
        x, y, weights, clusters, eligible, args.permutations, np.random.default_rng(args.seed)
    )
    permutation_seconds = time.perf_counter() - started

    # A ring with next-nearest chords has O(nodes) edges and O(nodes) local triangles.
    edge_map = {}
    for node in range(args.nodes):
        for offset in (1, 2):
            edge = tuple(sorted((node, (node + offset) % args.nodes)))
            edge_map[edge] = float(rng.random())
    edges = pd.DataFrame([
        {"gene1": f"g{left}", "gene2": f"g{right}", "score": score}
        for (left, right), score in sorted(edge_map.items())
    ])
    started = time.perf_counter()
    pruned = aracne_prune(edges, "score")
    aracne_seconds = time.perf_counter() - started
    print(json.dumps({
        "scope": "controlled computational benchmark only",
        "shc_permutation": {
            "samples": args.samples, "clusters": args.clusters,
            "permutations": args.permutations, "seconds": permutation_seconds,
            "permutations_per_second": args.permutations / permutation_seconds,
            "observed": statistic[0], "pvalue": statistic[3],
        },
        "aracne": {
            "nodes": args.nodes, "edges": len(edges), "seconds": aracne_seconds,
            "edges_per_second": len(edges) / aracne_seconds,
            "direct_edges": int(pruned.aracne_direct.sum()),
        },
    }, indent=2))


if __name__ == "__main__":
    main()
