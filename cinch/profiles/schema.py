"""Scientific state contract for CINCH presence and nominal allele profiles."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

import numpy as np


class LegacyMissingPolicy(str, Enum):
    """Required interpretation of non-positive values in historical matrices."""

    ABSENT = "absence"
    UNRESOLVED = "unresolved"


@dataclass(frozen=True)
class StateProfile:
    """Two-layer biological state with missingness distinct from absence.

    Presence values are -1 unresolved, 0 confidently absent, and 1 present.
    Type values are -1 non-callable or a non-negative nominal allele code.
    """

    samples: np.ndarray
    loci: np.ndarray
    presence: np.ndarray
    types: np.ndarray
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        samples = np.asarray(self.samples, dtype=str)
        loci = np.asarray(self.loci, dtype=str)
        presence = np.asarray(self.presence)
        types = np.asarray(self.types)
        if samples.ndim != 1 or loci.ndim != 1:
            raise ValueError("samples and loci must be one-dimensional")
        if presence.ndim != 2 or types.ndim != 2:
            raise ValueError("presence and types must be two-dimensional")
        expected = (samples.size, loci.size)
        if presence.shape != expected or types.shape != expected:
            raise ValueError(f"profile arrays must have shape {expected}")
        if len(set(samples.tolist())) != samples.size:
            raise ValueError("sample identifiers must be unique")
        if len(set(loci.tolist())) != loci.size:
            raise ValueError("locus identifiers must be unique")
        if presence.dtype.kind not in "biu" or types.dtype.kind not in "biu":
            raise TypeError("presence and type arrays must contain integers")
        if not np.isin(presence, (-1, 0, 1)).all():
            raise ValueError("presence values must be in {-1, 0, 1}")
        if (types < -1).any():
            raise ValueError("type values must be -1 or non-negative nominal codes")
        if ((types >= 0) & (presence != 1)).any():
            raise ValueError("a callable type requires presence=1")
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "loci", loci)
        object.__setattr__(self, "presence", presence.astype(np.int8, copy=False))
        object.__setattr__(self, "types", types.astype(np.int32, copy=False))
        object.__setattr__(self, "metadata", dict(self.metadata))


def _legacy_value_key(value: object) -> tuple[str, str]:
    return type(value).__name__, str(value)


def _legacy_sort_key(value: object) -> tuple[int, float, str, str]:
    try:
        return 0, float(value), type(value).__name__, str(value)
    except (TypeError, ValueError):
        return 1, 0.0, type(value).__name__, str(value)


def _legacy_is_present(value: object) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    if not text or text.casefold() in {"na", "nan", "none", "null", "."}:
        return False
    try:
        return bool(np.isfinite(float(value)) and float(value) > 0)
    except (TypeError, ValueError):
        return True


def from_legacy_alleles(
    alleles: np.ndarray,
    samples: np.ndarray,
    loci: np.ndarray,
    *,
    missing_policy: LegacyMissingPolicy | str,
) -> StateProfile:
    """Convert a legacy genome-by-locus allele matrix without guessing missingness.

    Positive values are observed nominal alleles. Zero, negative, NaN, and None
    are non-calls and require an explicit absent-versus-unresolved policy.
    """

    try:
        policy = LegacyMissingPolicy(missing_policy)
    except ValueError as error:
        raise ValueError("missing_policy must be 'absence' or 'unresolved'") from error
    matrix = np.asarray(alleles, dtype=object)
    expected = (len(samples), len(loci))
    if matrix.shape != expected:
        raise ValueError(f"legacy allele matrix must have shape {expected}")
    present = np.zeros(matrix.shape, dtype=bool)
    for index, value in np.ndenumerate(matrix):
        present[index] = _legacy_is_present(value)
    fill = 0 if policy is LegacyMissingPolicy.ABSENT else -1
    presence = np.full(matrix.shape, fill, dtype=np.int8)
    presence[present] = 1
    types = np.full(matrix.shape, -1, dtype=np.int32)
    codebooks: dict[str, list[str]] = {}
    for column, locus in enumerate(np.asarray(loci, dtype=str)):
        values = sorted(set(matrix[present[:, column], column].tolist()), key=_legacy_sort_key)
        mapping = {(_legacy_value_key(value)): code for code, value in enumerate(values)}
        for row in np.flatnonzero(present[:, column]):
            types[row, column] = mapping[_legacy_value_key(matrix[row, column])]
        codebooks[str(locus)] = [str(value) for value in values]
    return StateProfile(
        samples=np.asarray(samples, dtype=str),
        loci=np.asarray(loci, dtype=str),
        presence=presence,
        types=types,
        metadata={
            "schema": "CINCH_STATE_PROFILE_V2",
            "source": "legacy_allele_matrix",
            "legacy_missing_policy": policy.value,
            "legacy_allele_codebooks": codebooks,
        },
    )
