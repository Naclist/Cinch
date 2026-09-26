#!/usr/bin/env python3
"""Select deterministic SPN534 subsets stratified by tree SHC and fold."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--per-stratum", required=True, type=int)
    parser.add_argument("--manifest-output", required=True, type=Path)
    parser.add_argument("--accessions-output", required=True, type=Path)
    args = parser.parse_args()
    if args.per_stratum < 1:
        raise ValueError("per-stratum must be positive")
    source = pd.read_csv(args.manifest, sep="\t")
    required = {"sample_id", "tree_shc", "fold", "assembly_accession", "publication_matrix_column"}
    if not required.issubset(source.columns):
        raise ValueError(f"manifest is missing columns: {sorted(required - set(source.columns))}")
    if "in_genome_list" in source.columns:
        source = source[source["in_genome_list"].astype(str).str.casefold() == "true"]
    source = source.sort_values(["tree_shc", "fold", "sample_id"], kind="mergesort")
    selected = source.groupby(["tree_shc", "fold"], sort=True).head(args.per_stratum)
    sizes = selected.groupby(["tree_shc", "fold"]).size()
    if len(sizes) != 12 or not (sizes == args.per_stratum).all():
        raise ValueError("every tree_shc x fold stratum must supply the requested sample count")
    selected = selected[[
        "sample_id", "tree_shc", "fold", "assembly_accession", "publication_matrix_column",
    ]]
    if selected["assembly_accession"].duplicated().any():
        raise ValueError("selected assembly accessions are not unique")
    args.manifest_output.parent.mkdir(parents=True, exist_ok=True)
    args.accessions_output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.manifest_output, sep="\t", index=False)
    args.accessions_output.write_text(
        "\n".join(selected["assembly_accession"].astype(str)) + "\n", encoding="utf-8",
    )


if __name__ == "__main__":
    main()
