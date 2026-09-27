from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


@dataclass
class Hit:
    contig: str
    start: int
    end: int
    strand: str
    identity: float
    sequence: str


def fasta_records(path: Path) -> list[tuple[str, str, str]]:
    records = []
    for record in SeqIO.parse(path, "fasta"):
        sequence = str(record.seq).upper().replace("U", "T")
        if not sequence:
            raise ValueError(f"empty FASTA record {record.id} in {path}")
        records.append((str(record.id), sequence, str(record.description)))
    if not records:
        raise ValueError(f"no FASTA records in {path}")
    return records


def _candidate_starts(reference: str, target: str, seed_length: int = 13) -> set[int]:
    length = len(reference)
    k = min(seed_length, max(7, length // 8))
    offsets = sorted(set([0, max(0, length // 4 - k // 2), max(0, length // 2 - k // 2),
                          max(0, 3 * length // 4 - k // 2), max(0, length - k)]))
    starts: set[int] = set()
    for offset in offsets:
        seed = reference[offset:offset + k]
        position = target.find(seed)
        while position >= 0:
            start = position - offset
            if 0 <= start <= len(target) - length:
                starts.add(start)
            position = target.find(seed, position + 1)
    return starts


def find_hits(reference: str, contigs: list[tuple[str, str, str]], minimum_identity: float) -> list[Hit]:
    candidates: list[Hit] = []
    for contig, target, _ in contigs:
        for strand, query in (("+", reference), ("-", str(Seq(reference).reverse_complement()))):
            for start in _candidate_starts(query, target):
                observed = target[start:start + len(query)]
                callable_bases = sum(base in "ACGT" for base in observed)
                if callable_bases != len(query):
                    continue
                identity = sum(a == b for a, b in zip(query, observed)) / len(query)
                if identity >= minimum_identity:
                    biological = observed if strand == "+" else str(Seq(observed).reverse_complement())
                    candidates.append(Hit(contig, start + 1, start + len(query), strand, identity, biological))
    if not candidates:
        return []
    best = max(hit.identity for hit in candidates)
    return [hit for hit in candidates if abs(hit.identity - best) <= 1e-12]


def map_genomes(reference_path: Path, genome_paths: list[Path], minimum_identity: float):
    references = fasta_records(reference_path)
    locus_ids = np.asarray([x[0] for x in references], str)
    if len(set(locus_ids)) != len(locus_ids):
        raise ValueError("reference CDS identifiers must be unique")
    sample_ids = np.asarray([path.name.rsplit(".", 1)[0] for path in genome_paths], str)
    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("genome filenames produce duplicate sample identifiers")
    presence = np.zeros((len(genome_paths), len(references)), np.int8)
    types = np.full_like(presence, -1, dtype=np.int32)
    trace_rows, coordinate_rows = [], []
    per_locus_sequences: list[dict[str, int]] = [dict() for _ in references]
    genome_hashes = []
    for s, path in enumerate(genome_paths):
        contigs = fasta_records(path)
        genome_hashes.append({"sample_id": sample_ids[s], "path": str(path.resolve()),
                              "sha256": file_hash(path), "contigs": len(contigs),
                              "total_bp": sum(len(x[1]) for x in contigs)})
        sample_hits: list[tuple[int, Hit]] = []
        for j, (locus, reference, description) in enumerate(references):
            hits = find_hits(reference, contigs, minimum_identity)
            if not hits:
                trace_rows.append({"sample_id": sample_ids[s], "locus_id": locus, "presence_state": 0,
                                   "type_callable": False, "callability_reason": "CONFIDENT_ABSENCE",
                                   "contig": "", "start": np.nan, "end": np.nan, "strand": "",
                                   "CDS_SHA256": "", "type_ID": ""})
                continue
            presence[s, j] = 1
            sequences = sorted(set(hit.sequence for hit in hits))
            type_callable = len(sequences) == 1 and len(hits) == 1
            if type_callable:
                sequence = sequences[0]
                if sequence not in per_locus_sequences[j]:
                    per_locus_sequences[j][sequence] = len(per_locus_sequences[j])
                types[s, j] = per_locus_sequences[j][sequence]
            reason = "CALLABLE" if type_callable else "AMBIGUOUS_MULTICOPY"
            for hit in hits:
                sample_hits.append((j, hit))
                trace_rows.append({"sample_id": sample_ids[s], "locus_id": locus, "presence_state": 1,
                                   "type_callable": type_callable, "callability_reason": reason,
                                   "contig": hit.contig, "start": hit.start, "end": hit.end, "strand": hit.strand,
                                   "CDS_SHA256": hashlib.sha256(hit.sequence.encode()).hexdigest(),
                                   "type_ID": f"type_{types[s, j]}" if type_callable else "",
                                   "mapping_identity": hit.identity})
        by_contig: dict[str, list[tuple[int, Hit]]] = {}
        for item in sample_hits:
            by_contig.setdefault(item[1].contig, []).append(item)
        for contig, items in by_contig.items():
            for gene_order, (j, hit) in enumerate(sorted(items, key=lambda x: (x[1].start, x[1].end, x[0])), 1):
                coordinate_rows.append({"sample_id": sample_ids[s], "locus_id": locus_ids[j], "locus_index": j,
                                        "contig": contig, "start": hit.start, "end": hit.end,
                                        "strand": hit.strand, "gene_order": gene_order})
    annotations = pd.DataFrame({"locus_id": locus_ids,
                                "annotation": [description.partition(" ")[2] or locus for locus, _, description in references],
                                "reference_length": [len(sequence) for _, sequence, _ in references],
                                "reference_sha256": [hashlib.sha256(sequence.encode()).hexdigest() for _, sequence, _ in references]})
    return sample_ids, locus_ids, presence, types, pd.DataFrame(trace_rows), pd.DataFrame(coordinate_rows), annotations, pd.DataFrame(genome_hashes)


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def pair_order_distance(coordinates: pd.DataFrame, n_loci: int, minimum_observations: int) -> pd.DataFrame:
    accum: dict[tuple[int, int], list[tuple[float, float]]] = {}
    for (_, _), group in coordinates.groupby(["sample_id", "contig"], sort=False):
        records = group[["locus_index", "start", "end", "gene_order"]].to_numpy()
        local: dict[tuple[int, int], list[tuple[float, float]]] = {}
        for a in range(len(records) - 1):
            for b in range(a + 1, len(records)):
                ia, ib = int(records[a, 0]), int(records[b, 0])
                if ia == ib:
                    continue
                key = (min(ia, ib), max(ia, ib))
                bp = max(0.0, max(records[a, 1], records[b, 1]) - min(records[a, 2], records[b, 2]) - 1)
                order = abs(records[a, 3] - records[b, 3])
                local.setdefault(key, []).append((float(order), float(bp)))
        for key, values in local.items():
            accum.setdefault(key, []).append((min(x[0] for x in values), min(x[1] for x in values)))
    rows = []
    for i in range(n_loci - 1):
        for j in range(i + 1, n_loci):
            values = accum.get((i, j), [])
            # One mapping per locus is the v1 callable case. If duplicates occur,
            # mapping is ambiguous and therefore excluded upstream.
            reliable = len(values) >= minimum_observations
            rows.append({"i": i, "j": j,
                         "order_distance": float(np.mean([x[0] for x in values])) if reliable else np.nan,
                         "bp_distance": float(np.mean([x[1] for x in values])) if reliable else np.nan,
                         "n_same_contig_observations": len(values)})
    return pd.DataFrame(rows)
