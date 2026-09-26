from __future__ import annotations

import random
import json
from pathlib import Path

import numpy as np
from Bio.Seq import Seq

from cinch.frozen_v1.mapping import map_genomes
from cinch.mapping import MappingConfig, map_genomes_indexed
from cinch.mapping.mappy_backend import map_one_genome
from cinch.workflow.mapping import run_mapping


def sequence(seed: int, length: int = 600) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(f">{name}\n{value}\n" for name, value in records), encoding="utf-8")


def config() -> MappingConfig:
    return MappingConfig(
        minimum_identity=0.88,
        detection_coverage=0.60,
        callable_coverage=0.90,
        minimum_mapq=0,
        kmer_size=9,
        minimizer_window=3,
        minimum_chain_score=10,
        minimum_dp_score=10,
    )


def test_frozen_mapper_handles_a_zero_hit_genome(tmp_path: Path):
    ref = tmp_path / "ref.fasta"
    genome = tmp_path / "none.fasta"
    write_fasta(ref, [("locus", "A" * 80 + "C" * 80)])
    write_fasta(genome, [("contig", "G" * 500)])
    result = map_genomes(ref, [genome], 0.90)
    assert result[2].tolist() == [[0]]
    assert result[5].empty


def test_indexed_mapper_exact_snp_reverse_indel_and_no_hit(tmp_path: Path):
    base = sequence(1)
    snp_reference = sequence(5)
    snp = list(snp_reference)
    for index in (31, 143, 257, 399, 511):
        snp[index] = {"A": "C", "C": "G", "G": "T", "T": "A"}[snp[index]]
    snp = "".join(snp)
    reverse = str(Seq(sequence(2)).reverse_complement())
    indel_reference = sequence(3)
    indel = indel_reference[:260] + "GTA" + indel_reference[260:]
    absent = sequence(4)
    ref = tmp_path / "ref.fasta"
    genome = tmp_path / "sample.fasta"
    write_fasta(ref, [("exact", base), ("snp", snp_reference), ("reverse", sequence(2)),
                      ("indel", indel_reference), ("absent", absent)])
    write_fasta(genome, [("c1", "T" * 50 + base + "A" * 50),
                         ("c_snp", "G" * 50 + snp + "C" * 50),
                         ("c2", "C" * 50 + reverse + "G" * 50),
                         ("c3", "A" * 50 + indel + "T" * 50)])
    result, cache_hit = map_one_genome(ref, genome, tmp_path / "cache", tmp_path / "indexes", config())
    calls = {call.locus_id: call for call in result.calls}
    assert not cache_hit
    assert calls["exact"].type_callable
    assert calls["snp"].type_callable
    assert calls["reverse"].type_callable
    assert calls["reverse"].hits[0].strand == "-"
    assert calls["indel"].type_callable
    assert "I" in calls["indel"].hits[0].cigar or "D" in calls["indel"].hits[0].cigar
    assert calls["absent"].presence_state == -1
    assert calls["absent"].callability_reason == "UNRESOLVED_NO_HIT"


def test_indexed_mapper_marks_multicopy_and_partial_states(tmp_path: Path):
    full = sequence(10)
    partial = sequence(11)
    ref = tmp_path / "ref.fasta"
    genome = tmp_path / "sample.fasta"
    write_fasta(ref, [("duplicated", full), ("partial", partial)])
    write_fasta(genome, [
        ("copy1", full),
        ("copy2", full),
        ("fragment", partial[:450]),
    ])
    result, _ = map_one_genome(ref, genome, tmp_path / "cache", tmp_path / "indexes", config())
    calls = {call.locus_id: call for call in result.calls}
    assert calls["duplicated"].presence_state == 1
    assert not calls["duplicated"].type_callable
    assert calls["duplicated"].callability_reason == "AMBIGUOUS_MULTICOPY"
    assert calls["partial"].presence_state == 1
    assert not calls["partial"].type_callable
    assert calls["partial"].callability_reason == "PARTIAL_DETECTION"


def test_parallel_profile_is_deterministic_and_cache_is_reused(tmp_path: Path):
    first = sequence(20)
    second = list(first)
    second[200] = {"A": "C", "C": "G", "G": "T", "T": "A"}[second[200]]
    second = "".join(second)
    ref = tmp_path / "ref.fasta"
    g1, g2 = tmp_path / "g1.fasta", tmp_path / "g2.fasta"
    write_fasta(ref, [("locus", first)])
    write_fasta(g1, [("contig", first)])
    write_fasta(g2, [("contig", second)])
    serial = map_genomes_indexed(ref, [g1, g2], tmp_path / "serial", threads=1, config=config())
    parallel = map_genomes_indexed(ref, [g1, g2], tmp_path / "parallel", threads=2, config=config())
    np.testing.assert_array_equal(serial["presence"], parallel["presence"])
    np.testing.assert_array_equal(serial["types"], parallel["types"])
    assert serial["samples"].tolist() == ["g1", "g2"]
    assert sorted(serial["types"][:, 0].tolist()) == [0, 1]
    resumed = map_genomes_indexed(ref, [g1, g2], tmp_path / "serial", threads=2, config=config())
    assert resumed["cache_hits"].tolist() == [True, True]
    np.testing.assert_array_equal(serial["types"], resumed["types"])


def test_contig_boundary_is_not_called_absent(tmp_path: Path):
    locus = sequence(30)
    ref = tmp_path / "ref.fasta"
    genome = tmp_path / "split.fasta"
    write_fasta(ref, [("boundary", locus)])
    write_fasta(genome, [("left", locus[:300]), ("right", locus[300:])])
    result, _ = map_one_genome(ref, genome, tmp_path / "cache", tmp_path / "indexes", config())
    call = result.calls[0]
    assert call.presence_state == -1
    assert not call.type_callable


def test_competing_loci_and_ambiguous_bases_are_explicit(tmp_path: Path):
    target = sequence(40)
    close = list(target)
    for index in (50, 150, 250, 350, 450):
        close[index] = {"A": "C", "C": "G", "G": "T", "T": "A"}[close[index]]
    ambiguous = sequence(41)
    ambiguous_target = ambiguous[:200] + "N" * 200 + ambiguous[400:]
    ref = tmp_path / "ref.fasta"
    genome = tmp_path / "competing.fasta"
    write_fasta(ref, [("target_a", target), ("target_b", "".join(close)), ("ambiguous", ambiguous)])
    write_fasta(genome, [("shared", target), ("with_ns", ambiguous_target)])
    result, _ = map_one_genome(ref, genome, tmp_path / "cache", tmp_path / "indexes", config())
    calls = {call.locus_id: call for call in result.calls}
    assert calls["target_a"].callability_reason == "AMBIGUOUS_COMPETING_LOCI"
    assert calls["target_b"].callability_reason == "AMBIGUOUS_COMPETING_LOCI"
    assert not calls["target_a"].type_callable and not calls["target_b"].type_callable
    assert calls["ambiguous"].presence_state == -1
    assert calls["ambiguous"].callability_reason in {"PARTIAL_OR_LOW_CONFIDENCE", "UNRESOLVED_NO_HIT"}


def test_mapping_stage_writes_reusable_profile_and_manifest(tmp_path: Path):
    locus = sequence(50)
    ref = tmp_path / "ref.fasta"
    genome = tmp_path / "genome.fasta"
    output = tmp_path / "mapped"
    write_fasta(ref, [("locus", locus)])
    write_fasta(genome, [("contig", locus)])
    assert run_mapping(ref, [genome], output, threads=1, config=config()) == output.resolve()
    with np.load(output / "profiles" / "PROFILE_V2.npz") as profile:
        assert profile["samples"].tolist() == ["genome"]
        assert profile["loci"].tolist() == ["locus"]
        assert profile["presence"].tolist() == [[1]]
        assert profile["types"].tolist() == [[0]]
    manifest = json.loads((output / "manifest" / "MAPPING_MANIFEST.json").read_text())
    assert manifest["schema"] == "CINCH_MAPPING_RUN_V2"
    assert manifest["threads"] == 1
    assert manifest["genomes"][0]["cache_hit"] is False
    run_mapping(ref, [genome], output, threads=1, config=config())
    resumed = json.loads((output / "manifest" / "MAPPING_MANIFEST.json").read_text())
    assert resumed["genomes"][0]["cache_hit"] is True
