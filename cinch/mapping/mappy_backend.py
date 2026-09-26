"""Indexed nucleotide mapping backend powered by minimap2's mappy binding."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from dataclasses import replace

from .cache import cache_key, cache_path, genome_identity, read_cache, write_cache_atomic
from .models import GenomeMappingResult, LocusCall, MappingConfig, MappingHit


def _mappy():
    try:
        import mappy
    except ImportError as error:  # pragma: no cover - exercised by installation checks
        raise RuntimeError("indexed mapping requires the optional 'mappy' dependency") from error
    return mappy


def read_reference(path: Path) -> tuple[tuple[str, str], ...]:
    records = tuple((str(record.id), str(record.seq).upper().replace("U", "T")) for record in SeqIO.parse(path, "fasta"))
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    identifiers = [identifier for identifier, _ in records]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("reference CDS identifiers must be unique")
    if any(not sequence for _, sequence in records):
        raise ValueError("reference CDS sequences must be non-empty")
    return records


def _reference_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _build_or_load_index(genome: Path, genome_sha256: str, index_directory: Path, config: MappingConfig):
    mappy = _mappy()
    index_directory.mkdir(parents=True, exist_ok=True)
    target = index_directory / f"{genome_sha256}.mmi"
    if not target.exists():
        temporary = index_directory / f".{genome_sha256}.{os.getpid()}.mmi.tmp"
        aligner = mappy.Aligner(
            fn_idx_in=str(genome), preset="asm5", k=config.kmer_size, w=config.minimizer_window,
            n_threads=1, fn_idx_out=str(temporary), best_n=config.best_n,
            min_chain_score=config.minimum_chain_score, min_dp_score=config.minimum_dp_score,
        )
        if not aligner:
            raise RuntimeError(f"failed to index genome {genome}")
        os.replace(temporary, target)
    aligner = mappy.Aligner(
        fn_idx_in=str(target), preset="asm5", n_threads=1, best_n=config.best_n,
        min_chain_score=config.minimum_chain_score, min_dp_score=config.minimum_dp_score,
    )
    if not aligner:
        raise RuntimeError(f"failed to load genome index {target}")
    return aligner


def _oriented_sequence(aligner, hit) -> str:
    sequence = aligner.seq(hit.ctg, hit.r_st, hit.r_en)
    if sequence is None:
        return ""
    sequence = sequence.upper()
    return sequence if hit.strand >= 0 else str(Seq(sequence).reverse_complement())


def _convert_hit(locus_id: str, query_length: int, aligner, hit) -> MappingHit:
    return MappingHit(
        locus_id=locus_id,
        contig=str(hit.ctg),
        start=int(hit.r_st) + 1,
        end=int(hit.r_en),
        strand="+" if hit.strand >= 0 else "-",
        query_start=int(hit.q_st),
        query_end=int(hit.q_en),
        query_length=query_length,
        matches=int(hit.mlen),
        block_length=int(hit.blen),
        mapq=int(hit.mapq),
        edit_distance=int(hit.NM),
        cigar=str(hit.cigar_str),
        sequence=_oriented_sequence(aligner, hit),
        is_primary=bool(hit.is_primary),
    )


def _call_locus(locus_id: str, query: str, aligner, config: MappingConfig) -> LocusCall:
    raw = tuple(_convert_hit(locus_id, len(query), aligner, hit) for hit in aligner.map(query, cs=True))
    accepted = tuple(
        hit for hit in raw
        if hit.identity >= config.minimum_identity
        and hit.coverage >= config.detection_coverage
        and hit.mapq >= config.minimum_mapq
    )
    if not accepted:
        if raw:
            return LocusCall(locus_id, -1, False, "PARTIAL_OR_LOW_CONFIDENCE", None, raw)
        if config.no_hit_policy == "absence":
            return LocusCall(locus_id, 0, False, "CONFIDENT_ABSENCE_BY_POLICY", None, raw)
        return LocusCall(locus_id, -1, False, "UNRESOLVED_NO_HIT", None, raw)
    best_matches = max(hit.matches for hit in accepted)
    best = tuple(hit for hit in accepted if hit.matches == best_matches)
    full = tuple(hit for hit in best if hit.coverage >= config.callable_coverage)
    if len(full) == 1:
        return LocusCall(locus_id, 1, True, "CALLABLE", full[0].sequence, accepted)
    reason = "AMBIGUOUS_MULTICOPY" if len(best) > 1 else "PARTIAL_DETECTION"
    return LocusCall(locus_id, 1, False, reason, None, accepted)


def _overlap_fraction(left: MappingHit, right: MappingHit) -> float:
    if left.contig != right.contig:
        return 0.0
    overlap = max(0, min(left.end, right.end) - max(left.start, right.start) + 1)
    shorter = min(left.end - left.start + 1, right.end - right.start + 1)
    return overlap / shorter if shorter else 0.0


def _resolve_competing_loci(calls: tuple[LocusCall, ...], ratio: float) -> tuple[LocusCall, ...]:
    """Classify distinct loci whose best hits occupy the same genomic region."""

    ambiguous: set[int] = set()
    rejected: set[int] = set()
    best_hits = [call.hits[0] if call.hits else None for call in calls]
    for left_index in range(len(calls) - 1):
        left = best_hits[left_index]
        if left is None:
            continue
        for right_index in range(left_index + 1, len(calls)):
            right = best_hits[right_index]
            if right is None or _overlap_fraction(left, right) < 0.5:
                continue
            left_score, right_score = left.matches, right.matches
            weaker, stronger = sorted((left_score, right_score))
            if stronger == 0 or weaker / stronger >= ratio:
                ambiguous.update((left_index, right_index))
            elif left_score > right_score:
                rejected.add(right_index)
            else:
                rejected.add(left_index)
    resolved = []
    for index, call in enumerate(calls):
        if index in ambiguous:
            resolved.append(replace(call, presence_state=1, type_callable=False,
                                    callability_reason="AMBIGUOUS_COMPETING_LOCI", allele_sequence=None))
        elif index in rejected:
            resolved.append(replace(call, presence_state=-1, type_callable=False,
                                    callability_reason="COMPETING_LOCUS_REJECTED", allele_sequence=None))
        else:
            resolved.append(call)
    return tuple(resolved)


def map_one_genome(
    reference_path: Path,
    genome_path: Path,
    cache_directory: Path,
    index_directory: Path,
    config: MappingConfig,
) -> tuple[GenomeMappingResult, bool]:
    mappy = _mappy()
    reference_path, genome_path = reference_path.resolve(), genome_path.resolve()
    sample_id, genome_sha256 = genome_identity(genome_path)
    reference_sha256 = _reference_sha256(reference_path)
    version = str(mappy.__version__)
    key = cache_key(reference_sha256, genome_sha256, config, version)
    saved = cache_path(cache_directory, sample_id, key)
    if saved.exists():
        return read_cache(saved), True
    references = read_reference(reference_path)
    aligner = _build_or_load_index(genome_path, genome_sha256, index_directory, config)
    calls = tuple(_call_locus(identifier, sequence, aligner, config) for identifier, sequence in references)
    calls = _resolve_competing_loci(calls, config.competition_score_ratio)
    result = GenomeMappingResult(
        sample_id=sample_id,
        genome_path=str(genome_path),
        genome_sha256=genome_sha256,
        reference_sha256=reference_sha256,
        backend="mappy",
        backend_version=version,
        config=config.to_dict(),
        calls=calls,
    )
    write_cache_atomic(saved, result)
    return result, False
