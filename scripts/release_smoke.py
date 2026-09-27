#!/usr/bin/env python3
"""Run installed-command smoke tests without importing the source checkout."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import tempfile
from pathlib import Path


def _run(command: list[str], cwd: Path) -> None:
    subprocess.run(command, cwd=cwd, check=True)


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _resolve_executable(value: str) -> str:
    """Resolve explicit paths while leaving command names to PATH lookup."""
    candidate = Path(value).expanduser()
    if candidate.is_absolute() or candidate.parent != Path("."):
        resolved = candidate.resolve(strict=True)
        if not resolved.is_file():
            raise SystemExit(f"cinch executable is not a file: {resolved}")
        return str(resolved)
    resolved = shutil.which(value)
    if resolved is None:
        raise SystemExit(f"cinch executable not found on PATH: {value}")
    return resolved


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--cinch", default="cinch", help="installed cinch executable")
    args = parser.parse_args()
    repository = args.repository.resolve()
    executable = _resolve_executable(args.cinch)
    tiny = repository / "examples" / "tiny_wgs" / "input"
    conditional = repository / "examples" / "conditional_smoke"
    with tempfile.TemporaryDirectory(prefix="cinch-release-smoke-") as temporary:
        work = Path(temporary)
        wgs = work / "tiny.cinch"
        _run([
            executable, "wgs", "-r", str(tiny / "ref.cds.fasta"), "-p", "tiny",
            "-t", "2", "-o", str(wgs), "--phenotype", str(tiny / "phenotype.tsv"),
            *map(str, sorted((tiny / "genomes").glob("*.fasta"))),
        ], work)
        counts = {row["channel"]: int(row["eligible_pairs"])
                  for row in _rows(wgs / "10_report" / "WGS_COUNTS.tsv")}
        if counts != {"PP": 36, "PT": 32, "TP": 31, "TT": 21}:
            raise SystemExit(f"tiny WGS count mismatch: {counts}")
        _run([
            executable, "filter", "--wgs_results", str(wgs), "--hc", "2",
            "--order-threshold", "2",
        ], work)
        final = _rows(wgs / "filter" / "07_final" / "FINAL_CANDIDATES.tsv")
        if len(final) != 33:
            raise SystemExit(f"tiny filter candidate mismatch: {len(final)}")
        expected = {
            "PP:pp_co_A--pp_co_B", "PT:pt_presence--pt_type",
            "TT:tt_core_A--tt_core_B",
        }
        if not expected.issubset({row["edge_id"] for row in final}):
            raise SystemExit("tiny filter omitted an expected control edge")

        profile = work / "conditional-profile"
        _run([
            executable, "profile", str(conditional / "profile.tsv"), "-o", str(profile),
            "--missing-policy", "absence",
        ], work)
        result = work / "conditional-results"
        _run([
            executable, "conditional-associate", "--profile", str(profile / "PROFILE_V2.npz"),
            "--metadata", str(conditional / "metadata.tsv"), "--population-column", "population",
            "--background-column", "habitat", "--reference-background", "soil",
            "--comparison-background", "hospital", "--channels", "TT", "--pairs",
            str(conditional / "pairs.tsv"), "--permutations", "9", "--seed", "42",
            "--min-eligible-samples", "16", "--min-population-cell", "4",
            "--min-shared-populations", "2", "--min-background-samples", "8",
            "-o", str(result),
        ], work)
        manifest = json.loads((result / "MANIFEST.json").read_text(encoding="utf-8"))
        rows = _rows(result / "DIFFERENTIAL_ASSOCIATIONS.tsv")
        if manifest.get("status") != "COMPLETE" or manifest.get("tested_hypotheses") != 2:
            raise SystemExit(f"conditional manifest mismatch: {manifest}")
        if len(rows) != 1 or rows[0]["channel"] != "TT" or rows[0]["test_status"] != "TESTED":
            raise SystemExit(f"conditional result mismatch: {rows}")
        print(json.dumps({
            "tiny_wgs_counts": counts, "tiny_final_candidates": len(final),
            "conditional_tested_hypotheses": manifest["tested_hypotheses"],
            "conditional_differential_rows": len(rows),
        }, sort_keys=True))


if __name__ == "__main__":
    main()
