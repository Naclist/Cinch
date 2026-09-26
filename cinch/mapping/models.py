"""Stable data contracts for mapping backends."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field


@dataclass(frozen=True)
class MappingConfig:
    minimum_identity: float = 0.90
    detection_coverage: float = 0.60
    callable_coverage: float = 0.95
    minimum_mapq: int = 0
    best_n: int = 20
    kmer_size: int = 11
    minimizer_window: int = 5
    minimum_chain_score: int = 20
    minimum_dp_score: int = 20
    competition_score_ratio: float = 0.99
    backend: str = "mappy"
    schema_version: str = "CINCH_MAPPING_V2"

    def __post_init__(self) -> None:
        for name in ("minimum_identity", "detection_coverage", "callable_coverage"):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be in [0, 1]")
        if self.callable_coverage < self.detection_coverage:
            raise ValueError("callable_coverage must be >= detection_coverage")
        if not 0.0 < self.competition_score_ratio <= 1.0:
            raise ValueError("competition_score_ratio must be in (0, 1]")
        if self.best_n < 1 or self.kmer_size < 7 or self.minimizer_window < 1:
            raise ValueError("invalid mapper search parameters")

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class MappingHit:
    locus_id: str
    contig: str
    start: int
    end: int
    strand: str
    query_start: int
    query_end: int
    query_length: int
    matches: int
    block_length: int
    mapq: int
    edit_distance: int
    cigar: str
    sequence: str
    is_primary: bool

    @property
    def identity(self) -> float:
        query_span = self.query_end - self.query_start
        target_span = self.end - self.start + 1
        denominator = max(self.block_length, query_span, target_span)
        return self.matches / denominator if denominator else 0.0

    @property
    def coverage(self) -> float:
        return (self.query_end - self.query_start) / self.query_length if self.query_length else 0.0

    def to_dict(self) -> dict:
        value = asdict(self)
        value["identity"] = self.identity
        value["coverage"] = self.coverage
        return value

    @classmethod
    def from_dict(cls, value: dict) -> "MappingHit":
        fields = cls.__dataclass_fields__
        return cls(**{key: value[key] for key in fields})


@dataclass(frozen=True)
class LocusCall:
    locus_id: str
    presence_state: int
    type_callable: bool
    callability_reason: str
    allele_sequence: str | None
    hits: tuple[MappingHit, ...] = field(default_factory=tuple)

    def to_dict(self) -> dict:
        return {
            "locus_id": self.locus_id,
            "presence_state": self.presence_state,
            "type_callable": self.type_callable,
            "callability_reason": self.callability_reason,
            "allele_sequence": self.allele_sequence,
            "hits": [hit.to_dict() for hit in self.hits],
        }

    @classmethod
    def from_dict(cls, value: dict) -> "LocusCall":
        return cls(
            locus_id=value["locus_id"],
            presence_state=int(value["presence_state"]),
            type_callable=bool(value["type_callable"]),
            callability_reason=value["callability_reason"],
            allele_sequence=value.get("allele_sequence"),
            hits=tuple(MappingHit.from_dict(hit) for hit in value.get("hits", [])),
        )


@dataclass(frozen=True)
class GenomeMappingResult:
    sample_id: str
    genome_path: str
    genome_sha256: str
    reference_sha256: str
    backend: str
    backend_version: str
    config: dict
    calls: tuple[LocusCall, ...]

    def to_dict(self) -> dict:
        return {
            "schema": "CINCH_GENOME_MAPPING_RESULT_V2",
            "sample_id": self.sample_id,
            "genome_path": self.genome_path,
            "genome_sha256": self.genome_sha256,
            "reference_sha256": self.reference_sha256,
            "backend": self.backend,
            "backend_version": self.backend_version,
            "config": self.config,
            "calls": [call.to_dict() for call in self.calls],
        }

    @classmethod
    def from_dict(cls, value: dict) -> "GenomeMappingResult":
        if value.get("schema") != "CINCH_GENOME_MAPPING_RESULT_V2":
            raise ValueError("incompatible mapping result schema")
        return cls(
            sample_id=value["sample_id"],
            genome_path=value["genome_path"],
            genome_sha256=value["genome_sha256"],
            reference_sha256=value["reference_sha256"],
            backend=value["backend"],
            backend_version=value["backend_version"],
            config=value["config"],
            calls=tuple(LocusCall.from_dict(call) for call in value["calls"]),
        )
