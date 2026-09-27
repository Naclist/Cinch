#!/usr/bin/env python3
"""Fail when CINCH release archives omit required files or include excluded data."""

from __future__ import annotations

import argparse
import tarfile
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VERSION = "0.1.0"


def _source_names(path: Path) -> set[str]:
    with tarfile.open(path, "r:gz") as archive:
        names = set()
        for member in archive.getmembers():
            parts = Path(member.name).parts
            if len(parts) > 1:
                names.add(Path(*parts[1:]).as_posix())
        return names


def _wheel_names(path: Path) -> set[str]:
    with zipfile.ZipFile(path) as archive:
        return set(archive.namelist())


def _required_modules() -> set[str]:
    return {
        path.relative_to(ROOT).as_posix()
        for path in (ROOT / "cinch").rglob("*.py")
    }


def _require(names: set[str], required: set[str], label: str) -> None:
    missing = sorted(required - names)
    if missing:
        raise SystemExit(f"{label} missing required files: {missing}")


def verify_sdist(path: Path) -> None:
    names = _source_names(path)
    required = _required_modules() | {
        "LICENSE", "THIRD_PARTY_NOTICES.md", "RELEASE_NOTES.md", "README.md",
        "pyproject.toml", "FROZEN_CINCH_WGS_V1.yaml", "docs/CONDITIONAL_ASSOCIATION.md",
        "docs/development/RELEASE_LICENSE_AUDIT.md", "examples/tiny_wgs/README.md",
        "examples/conditional_smoke/README.md", "examples/conditional_smoke/profile.tsv",
        "examples/conditional_smoke/metadata.tsv", "examples/conditional_smoke/pairs.tsv",
    }
    _require(names, required, "sdist")
    forbidden = sorted(name for name in names if name.startswith("examples/spn534/"))
    if forbidden:
        raise SystemExit(f"sdist contains excluded SPN534 evidence: {forbidden[:5]}")


def verify_wheel(path: Path) -> None:
    names = _wheel_names(path)
    _require(names, _required_modules(), "wheel")
    if any(name.startswith("examples/") for name in names):
        raise SystemExit("wheel unexpectedly contains example/research data")
    license_names = {name for name in names if ".dist-info/licenses/" in name}
    if not any(name.endswith("/LICENSE") for name in license_names):
        raise SystemExit("wheel does not contain LICENSE")
    if not any(name.endswith("/THIRD_PARTY_NOTICES.md") for name in license_names):
        raise SystemExit("wheel does not contain THIRD_PARTY_NOTICES.md")
    wheel_docs = {
        "README.md", "RELEASE_NOTES.md", "FROZEN_CINCH_WGS_V1.yaml",
        "CONDITIONAL_ASSOCIATION.md", "STAGED_WORKFLOW.md", "KNOWN_LIMITATIONS.md",
        "SCIENTIFIC_EQUIVALENCE.md",
    }
    installed_docs = {
        Path(name).name for name in names
        if ".data/data/share/doc/cinch-wgs/" in name
    }
    if installed_docs != wheel_docs:
        raise SystemExit(
            f"wheel documentation mismatch: missing={sorted(wheel_docs - installed_docs)} "
            f"extra={sorted(installed_docs - wheel_docs)}"
        )
    metadata_name = next((name for name in names if name.endswith(".dist-info/METADATA")), None)
    entry_name = next((name for name in names if name.endswith(".dist-info/entry_points.txt")), None)
    if metadata_name is None or entry_name is None:
        raise SystemExit("wheel lacks METADATA or console entry points")
    with zipfile.ZipFile(path) as archive:
        metadata = archive.read(metadata_name).decode("utf-8")
        entries = archive.read(entry_name).decode("utf-8")
    if f"Version: {VERSION}" not in metadata or "License-Expression: MIT" not in metadata:
        raise SystemExit("wheel metadata version/license mismatch")
    if "cinch = cinch.cli:main" not in entries:
        raise SystemExit("wheel does not expose the cinch command")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", type=Path)
    args = parser.parse_args()
    sdists = [path for path in args.artifacts if path.name.endswith(".tar.gz")]
    wheels = [path for path in args.artifacts if path.suffix == ".whl"]
    if len(sdists) != 1 or len(wheels) != 1:
        raise SystemExit("expected exactly one sdist and one wheel")
    verify_sdist(sdists[0])
    verify_wheel(wheels[0])
    print(f"verified sdist={sdists[0].name} wheel={wheels[0].name}")


if __name__ == "__main__":
    main()
