"""Triangle-based ARACNE pruning with frozen tie semantics."""

from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd


def aracne_prune(edges: pd.DataFrame, weight_col: str) -> pd.DataFrame:
    required = {"gene1", "gene2", weight_col}
    if not required.issubset(edges.columns):
        raise ValueError(f"edges are missing columns: {sorted(required - set(edges.columns))}")
    edge_weight, neighbours = {}, {}
    for row in edges[["gene1", "gene2", weight_col]].itertuples(index=False):
        if row.gene1 == row.gene2:
            raise ValueError("ARACNE does not accept self-edges")
        key = tuple(sorted((row.gene1, row.gene2)))
        if key in edge_weight:
            raise ValueError(f"duplicate undirected edge: {key}")
        edge_weight[key] = float(getattr(row, weight_col))
        neighbours.setdefault(row.gene1, set()).add(row.gene2)
        neighbours.setdefault(row.gene2, set()).add(row.gene1)
    indirect = set()
    for center, adjacent in neighbours.items():
        for left, right in combinations(sorted(adjacent), 2):
            closing = tuple(sorted((left, right)))
            if closing not in edge_weight:
                continue
            triangle = [
                (tuple(sorted((center, left))), edge_weight[tuple(sorted((center, left)))]),
                (tuple(sorted((center, right))), edge_weight[tuple(sorted((center, right)))]),
                (closing, edge_weight[closing]),
            ]
            values = np.array([value for _, value in triangle])
            minimum = values.min()
            if np.sum(np.isclose(values, minimum, rtol=1e-12, atol=1e-12)) == 1:
                indirect.add(triangle[int(np.argmin(values))][0])
    output = edges.copy()
    output["aracne_direct"] = [
        tuple(sorted((left, right))) not in indirect
        for left, right in zip(output["gene1"], output["gene2"])
    ]
    return output
