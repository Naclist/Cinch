#!/usr/bin/env python3
"""Verify local source files against a Cinch relative-path SHA256 manifest."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()
    table = pd.read_csv(args.manifest, sep="\t")
    required = {"relative_path", "bytes", "sha256"}
    if not required.issubset(table.columns):
        parser.error(f"manifest needs columns: {sorted(required)}")

    failures = []
    for row in table.itertuples(index=False):
        path = args.root / Path(str(row.relative_path))
        if not path.is_file():
            failures.append((str(row.relative_path), "MISSING"))
            continue
        if path.stat().st_size != int(row.bytes):
            failures.append((str(row.relative_path), "SIZE_MISMATCH"))
            continue
        if sha256(path).lower() != str(row.sha256).lower():
            failures.append((str(row.relative_path), "SHA256_MISMATCH"))

    for path, reason in failures:
        print(f"{reason}\t{path}")
    print(f"verified={len(table) - len(failures)} failed={len(failures)} total={len(table)}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
