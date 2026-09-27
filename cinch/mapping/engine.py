"""Parallel orchestration and deterministic profile assembly."""

from __future__ import annotations

import hashlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from .mappy_backend import map_one_genome, read_reference
from .models import GenomeMappingResult, MappingConfig


def _worker(arguments):
    index, reference, genome, cache, indexes, config = arguments
    result, cache_hit = map_one_genome(Path(reference), Path(genome), Path(cache), Path(indexes), MappingConfig(**config))
    return index, result, cache_hit


def _deterministic_profiles(
    results: list[GenomeMappingResult], loci: tuple[str, ...]
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame, pd.DataFrame]:
    presence = np.full((len(results), len(loci)), -1, dtype=np.int8)
    types = np.full((len(results), len(loci)), -1, dtype=np.int32)
    allele_sets: list[set[str]] = [set() for _ in loci]
    for sample_index, result in enumerate(results):
        if tuple(call.locus_id for call in result.calls) != loci:
            raise ValueError(f"locus order mismatch for {result.sample_id}")
        for locus_index, call in enumerate(result.calls):
            presence[sample_index, locus_index] = call.presence_state
            if call.type_callable and call.allele_sequence is not None:
                allele_sets[locus_index].add(call.allele_sequence)
    allele_maps = [
        {sequence: code for code, sequence in enumerate(sorted(values, key=lambda value: (hashlib.sha256(value.encode()).hexdigest(), value)))}
        for values in allele_sets
    ]
    trace, coordinates = [], []
    for sample_index, result in enumerate(results):
        for locus_index, call in enumerate(result.calls):
            if call.type_callable and call.allele_sequence is not None:
                types[sample_index, locus_index] = allele_maps[locus_index][call.allele_sequence]
            trace.append({
                "sample_id": result.sample_id,
                "locus_id": call.locus_id,
                "presence_state": call.presence_state,
                "type_callable": call.type_callable,
                "callability_reason": call.callability_reason,
                "type_ID": int(types[sample_index, locus_index]) if types[sample_index, locus_index] >= 0 else "",
                "n_hits": len(call.hits),
            })
            for hit in call.hits if call.presence_state == 1 else ():
                coordinates.append({
                    "sample_id": result.sample_id,
                    "locus_id": call.locus_id,
                    "locus_index": locus_index,
                    "contig": hit.contig,
                    "start": hit.start,
                    "end": hit.end,
                    "strand": hit.strand,
                    "identity": hit.identity,
                    "coverage": hit.coverage,
                    "mapq": hit.mapq,
                    "type_callable": call.type_callable,
                    "callability_reason": call.callability_reason,
                })
    coordinate_frame = pd.DataFrame(coordinates)
    if not coordinate_frame.empty:
        coordinate_frame = coordinate_frame.sort_values(
            ["sample_id", "contig", "start", "end", "locus_index"], kind="mergesort"
        ).reset_index(drop=True)
        coordinate_frame["gene_order"] = coordinate_frame.groupby(["sample_id", "contig"]).cumcount() + 1
    return presence, types, pd.DataFrame(trace), coordinate_frame


def map_genomes_indexed(
    reference_path: Path,
    genome_paths: list[Path],
    output_directory: Path,
    threads: int = 1,
    config: MappingConfig | None = None,
) -> dict:
    """Map genomes with real process-level parallelism and resumable per-genome caches."""

    if threads < 1:
        raise ValueError("threads must be >= 1")
    if not genome_paths:
        raise ValueError("at least one genome is required")
    config = config or MappingConfig()
    reference_path = reference_path.resolve()
    genome_paths = [path.resolve() for path in genome_paths]
    sample_ids = [path.name.rsplit(".", 1)[0] for path in genome_paths]
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("genome filenames produce duplicate sample identifiers")
    output_directory.mkdir(parents=True, exist_ok=True)
    cache = output_directory / "cache" / "genomes"
    indexes = output_directory / "cache" / "indexes"
    tasks = [
        (index, str(reference_path), str(genome), str(cache), str(indexes), config.to_dict())
        for index, genome in enumerate(genome_paths)
    ]
    completed: dict[int, tuple[GenomeMappingResult, bool]] = {}
    if threads == 1:
        for task in tasks:
            index, result, cache_hit = _worker(task)
            completed[index] = (result, cache_hit)
    else:
        with ProcessPoolExecutor(max_workers=min(threads, len(tasks))) as executor:
            futures = {executor.submit(_worker, task): task[0] for task in tasks}
            for future in as_completed(futures):
                index, result, cache_hit = future.result()
                completed[index] = (result, cache_hit)
    results = [completed[index][0] for index in range(len(tasks))]
    loci = tuple(identifier for identifier, _ in read_reference(reference_path))
    presence, types, trace, coordinates = _deterministic_profiles(results, loci)
    return {
        "samples": np.asarray([result.sample_id for result in results], dtype=str),
        "loci": np.asarray(loci, dtype=str),
        "presence": presence,
        "types": types,
        "trace": trace,
        "coordinates": coordinates,
        "results": results,
        "cache_hits": np.asarray([completed[index][1] for index in range(len(tasks))], dtype=bool),
    }
