#!/usr/bin/env python3
"""Create a deterministic nonhistorical CDS reference and SHC table for SPN534."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cds", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--loci", type=int, default=200)
    parser.add_argument("--reference-output", required=True, type=Path)
    parser.add_argument("--shc-output", required=True, type=Path)
    args = parser.parse_args()
    records = [
        record for record in SeqIO.parse(args.cds, "fasta")
        if 150 <= len(record.seq) <= 3000 and set(str(record.seq).upper()) <= set("ACGT")
    ]
    records = sorted(records, key=lambda record: record.id)[: args.loci]
    if len(records) < args.loci:
        raise ValueError(f"only {len(records)} eligible CDS records")
    args.reference_output.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, args.reference_output, "fasta")
    manifest = pd.read_csv(args.manifest, sep="\t")
    shc = pd.DataFrame({
        "sample_id": manifest.sample_id.astype(str) + "_genomic",
        "shc": manifest.tree_shc,
    }).sort_values("sample_id")
    args.shc_output.parent.mkdir(parents=True, exist_ok=True)
    shc.to_csv(args.shc_output, sep="\t", index=False)


if __name__ == "__main__":
    main()
