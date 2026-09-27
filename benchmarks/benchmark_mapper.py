#!/usr/bin/env python3
"""Controlled computational scaling benchmark for indexed genome mapping."""

from __future__ import annotations

import argparse
import json
import random
import tempfile
import time
from pathlib import Path

from cinch.mapping import MappingConfig, map_genomes_indexed


def sequences(count: int, length: int, seed: int) -> list[str]:
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(length)) for _ in range(count)]


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(f">{name}\n{value}\n" for name, value in records), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--genomes", type=int, default=12)
    parser.add_argument("--loci", type=int, default=200)
    parser.add_argument("--length", type=int, default=300)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument("--seed", type=int, default=20260926)
    args = parser.parse_args()
    loci = sequences(args.loci, args.length, args.seed)
    config = MappingConfig(kmer_size=9, minimizer_window=3, minimum_chain_score=10, minimum_dp_score=10)
    rows = []
    with tempfile.TemporaryDirectory(prefix="cinch-mapper-benchmark-") as temporary:
        root = Path(temporary)
        reference = root / "reference.fasta"
        write_fasta(reference, [(f"locus_{index:05d}", value) for index, value in enumerate(loci)])
        genomes = []
        for sample in range(args.genomes):
            genome = root / f"genome_{sample:03d}.fasta"
            # Distinct contigs prevent separator sequence from becoming alignment evidence.
            write_fasta(genome, [(f"contig_{index:05d}", value) for index, value in enumerate(loci)])
            genomes.append(genome)
        for workers in args.threads:
            started = time.perf_counter()
            mapped = map_genomes_indexed(
                reference, genomes, root / f"run_{workers}", threads=workers, config=config
            )
            elapsed = time.perf_counter() - started
            rows.append({
                "workers": workers,
                "seconds": elapsed,
                "genome_locus_queries": args.genomes * args.loci,
                "queries_per_second": args.genomes * args.loci / elapsed,
                "callable": int((mapped["types"] >= 0).sum()),
                "cache_hits": int(mapped["cache_hits"].sum()),
            })
    baseline = rows[0]["seconds"]
    for row in rows:
        row["speedup_vs_first"] = baseline / row["seconds"]
    print(json.dumps({
        "scope": "controlled computational benchmark only; not biological validation",
        "genomes": args.genomes,
        "loci": args.loci,
        "length": args.length,
        "seed": args.seed,
        "runs": rows,
    }, indent=2))


if __name__ == "__main__":
    main()
