from hashlib import sha256
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPN = ROOT / "examples" / "spn534"


def test_frozen_candidate_invariants():
    path = SPN / "results" / "FINAL_CANDIDATES.tsv"
    assert sha256(path.read_bytes()).hexdigest() == "07bb301b1d98f9df17e0e14cc5e8b4a780e6ad507062eed21fd2365df943dfe5"
    df = pd.read_csv(path, sep="\t")
    assert len(df) == 223
    assert df.channel.value_counts().to_dict() == {"TP": 160, "PT": 58, "PP": 4, "TT": 1}
    assert df[["locus_A", "locus_B"]].apply(lambda r: tuple(sorted(r)), axis=1).nunique() == 223


def test_network_is_exact_frozen_pair_universe():
    final = pd.read_csv(SPN / "results" / "FINAL_CANDIDATES.tsv", sep="\t")
    edges = pd.read_csv(SPN / "network" / "NETWORK_EDGES.tsv", sep="\t")
    nodes = pd.read_csv(SPN / "network" / "NETWORK_NODES.tsv", sep="\t")
    frozen = {tuple(sorted(x)) for x in final[["locus_A", "locus_B"]].itertuples(index=False, name=None)}
    graph = {tuple(sorted(x)) for x in edges[["locus_A", "locus_B"]].itertuples(index=False, name=None)}
    assert graph == frozen
    assert len(edges) == 223
    assert len(nodes) == 167


def test_methods_paper_assets_are_complete():
    figs = ROOT / "docs" / "assets" / "figures"
    for number in range(1, 18):
        stem = f"FIG{number:02d}_"
        matching = [p for p in figs.iterdir() if p.name.startswith(stem)]
        assert {p.suffix for p in matching} == {".png", ".pdf", ".svg"}
    assert (SPN / "network" / "CINCH_INTERACTIVE.html").stat().st_size > 50_000
    assert (ROOT / "site" / "index.html").exists()


def test_order_and_bp_availability_is_explicit():
    table = pd.read_csv(SPN / "diagnostics" / "TABLE_DISTANCE_COMPARISON.tsv", sep="\t")
    order = table.set_index("distance").loc["order"]
    bp = table.set_index("distance").loc["bp"]
    assert int(order.available_pairs) == int(bp.available_pairs) == 990_661
    assert int(order.order_only_available) == int(order.bp_only_available) == 0
    assert float(order.order_bp_spearman_rho_when_both) > 0.99
