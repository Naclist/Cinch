from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cinch import __version__
from cinch.cli import build_parser
from cinch.frozen_v1.hierarchy import scan_hc, select_hc_level
from cinch.frozen_v1.mapping import pair_order_distance
from cinch.frozen_v1.statistics import channel_vectors, contingency, neff_from_counts


def test_missing_is_not_absence_and_tt_is_pairwise_complete():
    presence = np.array([[1, 1], [1, 1], [1, 1], [-1, 1]], np.int8)
    types = np.array([[0, 0], [1, 1], [-1, 0], [-1, 1]], np.int16)
    pp_x, _ = channel_vectors(presence, types, 0, 1, "PP")
    tt_x, tt_y = channel_vectors(presence, types, 0, 1, "TT")
    assert len(pp_x) == 3
    assert len(tt_x) == 2
    assert -1 not in pp_x and -1 not in tt_x and -1 not in tt_y


def test_always_present_pair_can_be_pp_zero_and_tt_high():
    n = 24
    presence = np.ones((n, 2), np.int8)
    types = np.column_stack([np.arange(n) % 4, np.arange(n) % 4]).astype(np.int16)
    pp = contingency(*channel_vectors(presence, types, 0, 1, "PP"))
    tt = contingency(*channel_vectors(presence, types, 0, 1, "TT"))
    assert pp["MIraw"] == pytest.approx(0.0)
    assert tt["MIraw"] == pytest.approx(np.log(4))


def test_association_and_mutual_exclusion_have_dependence():
    x = np.tile([0, 1], 20)
    associated = contingency(x, x)
    excluded = contingency(x, 1 - x)
    assert associated["MIraw"] == pytest.approx(np.log(2))
    assert excluded["MIraw"] == pytest.approx(np.log(2))
    assert associated["enriched_driver_A_state"] == associated["enriched_driver_B_state"]
    assert excluded["enriched_driver_A_state"] != excluded["enriched_driver_B_state"]


def test_order_distance_and_different_contig_na():
    coords = pd.DataFrame([
        {"sample_id": "s1", "contig": "c1", "locus_index": 0, "start": 1, "end": 10, "gene_order": 1},
        {"sample_id": "s1", "contig": "c1", "locus_index": 1, "start": 21, "end": 30, "gene_order": 4},
        {"sample_id": "s2", "contig": "c1", "locus_index": 0, "start": 1, "end": 10, "gene_order": 2},
        {"sample_id": "s2", "contig": "c1", "locus_index": 1, "start": 31, "end": 40, "gene_order": 6},
        {"sample_id": "s1", "contig": "c2", "locus_index": 2, "start": 1, "end": 10, "gene_order": 1},
    ])
    result = pair_order_distance(coords, 3, 2).set_index(["i", "j"])
    assert result.loc[(0, 1), "order_distance"] == pytest.approx(3.5)
    assert np.isnan(result.loc[(0, 2), "order_distance"])


def test_hcx_is_a_distance_threshold_not_number_of_clusters():
    distance = np.array([[0, 1, 6, 6], [1, 0, 6, 6], [6, 6, 0, 1], [6, 6, 1, 0]])
    _, assignments = scan_hc(distance, 6)
    assert len(np.unique(assignments[1])) == 2


def test_hc_selection_has_pass_ambiguous_and_unresolved_states():
    def frame(pattern):
        n = len(pattern)
        return pd.DataFrame({"HC_level": range(n), "n_clusters": [4] * n,
                             "singleton_fraction": [0.0] * n,
                             "NMI_to_previous": pattern, "ARI_to_previous": pattern})
    assert select_hc_level(frame([0, 0, 1, 1, 0]))["status"] == "PASS"
    assert select_hc_level(frame([0, 0, 1, 1, 0, 1, 1]))["status"] == "AMBIGUOUS"
    assert select_hc_level(frame([0, 0, 0, 0]))["status"] == "UNRESOLVED"


def test_neff_exact():
    assert neff_from_counts({7: 1, 16: 3, 27: 3}) == pytest.approx(2.5789473684210527)


def test_cli_exposes_wgs_and_filter():
    parser = build_parser()
    assert parser.parse_args(["map", "g.fasta", "-r", "ref.fasta", "-o", "mapped"]).command == "map"
    assert parser.parse_args(["profile", "x.tsv", "-o", "p", "--missing-policy", "unresolved"]).command == "profile"
    assert parser.parse_args(["associate", "--profile", "p.npz", "--coordinates", "c.tsv", "-o", "a"]).command == "associate"
    assert parser.parse_args(["advanced-filter", "--association", "a", "--profile", "p.npz", "--shc", "s.tsv", "-o", "f"]).command == "advanced-filter"
    assert parser.parse_args(["report", "--filter-results", "f", "-o", "r"]).command == "report"
    assert parser.parse_args(["filter", "--wgs_results", "x", "--hc", "69",
                              "--order-threshold", "100"]).command == "filter"


def test_public_package_version_is_distinct_from_frozen_workflow_version():
    from cinch.frozen_v1 import VERSION as frozen_workflow_version
    assert __version__ == "0.1.0"
    assert frozen_workflow_version == "1.0.0"


def test_release_smoke_resolves_command_name_from_path(monkeypatch):
    from scripts.release_smoke import _resolve_executable

    monkeypatch.setattr("scripts.release_smoke.shutil.which", lambda value: "/tmp/bin/cinch")
    assert _resolve_executable("cinch") == "/tmp/bin/cinch"
