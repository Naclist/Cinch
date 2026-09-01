from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_rows(paths: Iterable[Path], root: Path | None = None) -> pd.DataFrame:
    rows = []
    for path in sorted({Path(p).resolve() for p in paths}):
        label = str(path.relative_to(root.resolve())) if root and path.is_relative_to(root.resolve()) else str(path)
        rows.append({"path": label, "bytes": path.stat().st_size, "sha256": sha256(path)})
    return pd.DataFrame(rows)


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def write_yaml(value: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(value, sort_keys=False, allow_unicode=True), encoding="utf-8")


def write_json(value: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, default=json_default) + "\n", encoding="utf-8")


def json_default(value):
    if isinstance(value, (np.integer, np.floating, np.bool_)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def software_versions() -> dict[str, str]:
    import Bio
    import matplotlib
    import scipy
    import sklearn

    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "biopython": Bio.__version__,
    }


class RunLog:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = path.open("w", encoding="utf-8")

    def emit(self, message: str = "") -> None:
        print(message, flush=True)
        self.handle.write(message + "\n")
        self.handle.flush()

    def stage(self, stage: str, message: str) -> None:
        self.emit(f"[{stage}] {message}")

    def close(self) -> None:
        self.handle.close()

