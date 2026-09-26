"""Content-addressed, atomic per-genome mapping cache."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import tempfile
from pathlib import Path

from cinch.frozen_v1.mapping import file_hash

from .models import GenomeMappingResult, MappingConfig


def cache_key(reference_sha256: str, genome_sha256: str, config: MappingConfig, backend_version: str) -> str:
    payload = json.dumps(
        {
            "reference_sha256": reference_sha256,
            "genome_sha256": genome_sha256,
            "config": config.to_dict(),
            "backend_version": backend_version,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def genome_identity(path: Path) -> tuple[str, str]:
    return path.name.rsplit(".", 1)[0], file_hash(path)


def cache_path(directory: Path, sample_id: str, key: str) -> Path:
    safe = "".join(char if char.isalnum() or char in "._-" else "_" for char in sample_id)
    return directory / f"{safe}.{key}.json.gz"


def read_cache(path: Path) -> GenomeMappingResult:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return GenomeMappingResult.from_dict(json.load(handle))


def write_cache_atomic(path: Path, result: GenomeMappingResult) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(descriptor)
    temporary_path = Path(temporary)
    try:
        with gzip.open(temporary_path, "wt", encoding="utf-8") as handle:
            json.dump(result.to_dict(), handle, sort_keys=True, separators=(",", ":"))
        os.replace(temporary_path, path)
    finally:
        temporary_path.unlink(missing_ok=True)
