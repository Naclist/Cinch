"""User-facing indexed mapping stage."""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd

from cinch.frozen_v1.mapping import file_hash
from cinch.frozen_v1.utils import software_versions, write_json, write_tsv
from cinch.mapping import MappingConfig, map_genomes_indexed
from cinch.profiles import StateProfile, save_profile


def run_mapping(
    reference: Path,
    genomes: list[Path],
    output: Path,
    threads: int,
    config: MappingConfig,
) -> Path:
    """Execute mapping independently and emit reusable profiles plus provenance."""

    started = time.perf_counter()
    output = output.resolve()
    for directory in (output, output / "profiles", output / "mapping", output / "manifest"):
        directory.mkdir(parents=True, exist_ok=True)
    mapped = map_genomes_indexed(reference, genomes, output, threads=threads, config=config)
    save_profile(
        StateProfile(
            samples=mapped["samples"], loci=mapped["loci"],
            presence=mapped["presence"], types=mapped["types"],
            metadata={
                "schema": "CINCH_STATE_PROFILE_V2",
                "source": "indexed_mapping",
                "presence_states": {"-1": "unresolved", "0": "absent", "1": "present"},
                "type_missing": -1,
            },
        ),
        output / "profiles" / "PROFILE_V2.npz",
    )
    write_tsv(mapped["trace"], output / "mapping" / "CALL_TRACE.tsv")
    write_tsv(mapped["coordinates"], output / "mapping" / "COORDINATES.tsv")
    genome_rows = []
    for index, result in enumerate(mapped["results"]):
        genome_rows.append({
            "sample_id": result.sample_id,
            "path": result.genome_path,
            "sha256": result.genome_sha256,
            "cache_hit": bool(mapped["cache_hits"][index]),
            "present_loci": sum(call.presence_state == 1 for call in result.calls),
            "callable_loci": sum(call.type_callable for call in result.calls),
            "unresolved_loci": sum(call.presence_state < 0 for call in result.calls),
        })
    write_tsv(pd.DataFrame(genome_rows), output / "mapping" / "GENOME_MAPPING_QC.tsv")
    manifest = {
        "schema": "CINCH_MAPPING_RUN_V2",
        "reference": {"path": str(reference.resolve()), "sha256": file_hash(reference)},
        "genomes": genome_rows,
        "samples": len(mapped["samples"]),
        "loci": len(mapped["loci"]),
        "threads": threads,
        "config": config.to_dict(),
        "backend": mapped["results"][0].backend,
        "backend_version": mapped["results"][0].backend_version,
        "software": software_versions(),
        "elapsed_seconds": time.perf_counter() - started,
        "profile": "profiles/PROFILE_V2.npz",
        "call_trace": "mapping/CALL_TRACE.tsv",
        "coordinates": "mapping/COORDINATES.tsv",
    }
    write_json(manifest, output / "manifest" / "MAPPING_MANIFEST.json")
    return output
