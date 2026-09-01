#!/usr/bin/env python3
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parent
N = 48
SEED = 314159
LOCI = [
    *(f"bg_{i}" for i in range(1, 7)),
    "pp_co_A", "local_A", "local_B", "pp_ex_A", "pt_presence",
    "neutral_1", "pp_co_B", "pt_type", "pp_ex_B", "tt_core_A", "neutral_2", "tt_core_B",
    "pop_A", "pop_B",
]


def reference_sequences():
    rng = np.random.default_rng(SEED)
    sequences = {}
    for locus in LOCI:
        while True:
            sequence = "".join(rng.choice(list("ACGT"), 180))
            if all(sequence not in value and value not in sequence for value in sequences.values()):
                sequences[locus] = sequence; break
    return sequences


def allele(reference: str, state: int) -> str:
    if state == 0: return reference
    sequence = list(reference)
    for offset in (23 + state, 83 + 2 * state, 137 + state):
        old = sequence[offset % len(sequence)]
        sequence[offset % len(sequence)] = "ACGT"[("ACGT".index(old) + (state % 3) + 1) % 4]
    return "".join(sequence)


def main():
    genomes = ROOT / "genomes"; genomes.mkdir(parents=True, exist_ok=True)
    refs = reference_sequences()
    with (ROOT / "ref.cds.fasta").open("w", encoding="ascii") as handle:
        for locus, sequence in refs.items():
            background = " [BACKGROUND]" if locus.startswith("bg_") else ""
            handle.write(f">{locus} synthetic {locus} product{background}\n{sequence}\n")
    phenotype = ["sample_id\tbackground\tphenotype"]
    truth = []
    for s in range(N):
        block, within = s // 12, s % 12
        bit = within % 2
        states = {locus: 0 for locus in LOCI}
        present = {locus: True for locus in LOCI}
        for locus in LOCI:
            if locus.startswith("bg_"): states[locus] = block + 1
        # One within-block background variant preserves a clear HC plateau.
        if within == 11: states[f"bg_{block % 6 + 1}"] = block + 5
        present["pp_co_A"] = present["pp_co_B"] = bool(bit)
        present["pp_ex_A"] = not bool(bit); present["pp_ex_B"] = bool(bit)
        present["local_A"] = present["local_B"] = bool(bit)
        present["pt_presence"] = bool(bit); states["pt_type"] = bit + 1
        states["tt_core_A"] = within % 4 + 1; states["tt_core_B"] = within % 4 + 1
        present["pop_A"] = present["pop_B"] = block == 0 and bit == 1
        states["neutral_1"] = (s * 7) % 3; states["neutral_2"] = (s * 5 + 1) % 3
        fragments = []
        for locus in LOCI:
            if present[locus]: fragments.append(allele(refs[locus], states[locus]))
        # Deliberate wgMLST non-call: two equally good, distinct full CDS copies.
        if s in (3, 15):
            fragments.append(allele(refs["bg_6"], block + 2))
        genome = ("N" * 35).join(fragments)
        sample = f"tiny_{s+1:03d}"
        (genomes / f"{sample}.fasta").write_text(f">chromosome\n{genome}\n", encoding="ascii")
        phenotype.append(f"{sample}\tB{block+1}\t{int(bit and block in (1,3))}")
    (ROOT / "phenotype.tsv").write_text("\n".join(phenotype) + "\n", encoding="utf-8")
    for channel, a, b, expected in [
        ("PP","pp_co_A","pp_co_B","association"), ("PP","pp_ex_A","pp_ex_B","dissociation"),
        ("PT","pt_presence","pt_type","dependency"), ("TT","tt_core_A","tt_core_B","PP_zero_TT_high"),
        ("PP","local_A","local_B","fails_order"), ("PP","pop_A","pop_B","population_confined")]:
        truth.append({"channel":channel,"locus_A":a,"locus_B":b,"expected":expected})
    expected = {"seed": SEED, "samples": N, "loci": len(LOCI), "manual_HC": 2,
                "manual_order_threshold": 2, "truth": truth,
                "expected_final_oriented_counts": {"PP": 12, "PT": 18, "TP": 0, "TT": 3, "total": 33},
                "expected_retained": ["PP:pp_co_A--pp_co_B", "PT:pt_presence--pt_type", "TT:tt_core_A--tt_core_B"],
                "expected_rejected": ["PP:local_A--local_B", "PP:pop_A--pop_B"]}
    (ROOT / "EXPECTED_RESULTS.yaml").write_text(yaml.safe_dump(expected, sort_keys=False), encoding="utf-8")


if __name__ == "__main__": main()
