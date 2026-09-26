import numpy as np
import pandas as pd
import pytest

from cinch.distance import SparseOrderDistance
from cinch.frozen_v1.mapping import pair_order_distance


def test_sparse_distance_matches_frozen_and_preserves_missing_categories():
    coordinates = pd.DataFrame([
        {"sample_id": "s1", "contig": "c", "locus_index": 0, "start": 1, "end": 10, "gene_order": 1},
        {"sample_id": "s1", "contig": "c", "locus_index": 1, "start": 31, "end": 40, "gene_order": 3},
        {"sample_id": "s2", "contig": "c", "locus_index": 0, "start": 2, "end": 11, "gene_order": 2},
        {"sample_id": "s2", "contig": "c", "locus_index": 1, "start": 42, "end": 51, "gene_order": 5},
        {"sample_id": "s1", "contig": "c", "locus_index": 2, "start": 61, "end": 70, "gene_order": 6},
    ])
    sparse = SparseOrderDistance(coordinates, 4, minimum_observations=2)
    frozen = pair_order_distance(coordinates, 4, 2).set_index(["i", "j"])
    reliable = sparse.get_pair(0, 1)
    assert reliable["order_distance"] == frozen.loc[(0, 1), "order_distance"]
    assert reliable["bp_distance"] == frozen.loc[(0, 1), "bp_distance"]
    assert reliable["physical_distance_status"] == "RELIABLE_SAME_CONTIG_DISTANCE"
    assert sparse.get_pair(0, 2)["physical_distance_status"] == "INSUFFICIENT_SAME_CONTIG_OBSERVATIONS"
    absent = sparse.get_pair(0, 3)
    assert absent["physical_distance_status"] == "NO_SAME_CONTIG_OBSERVATION"
    assert np.isnan(absent["order_distance"])
    assert len(sparse.observed_frame()) == 3


def test_compiled_sparse_distance_matches_python_with_multicopy_coordinates():
    pytest.importorskip("numba")
    coordinates = pd.DataFrame([
        {"sample_id": "s1", "contig": "a", "locus_index": 0, "start": 1, "end": 10, "gene_order": 1},
        {"sample_id": "s1", "contig": "a", "locus_index": 0, "start": 101, "end": 110, "gene_order": 5},
        {"sample_id": "s1", "contig": "a", "locus_index": 1, "start": 31, "end": 40, "gene_order": 3},
        {"sample_id": "s1", "contig": "a", "locus_index": 2, "start": 61, "end": 70, "gene_order": 4},
        {"sample_id": "s2", "contig": "b", "locus_index": 0, "start": 5, "end": 14, "gene_order": 2},
        {"sample_id": "s2", "contig": "b", "locus_index": 1, "start": 45, "end": 54, "gene_order": 6},
    ])
    python = SparseOrderDistance(coordinates, 4, minimum_observations=2, engine="python")
    compiled = SparseOrderDistance(coordinates, 4, minimum_observations=2, engine="numba")
    pd.testing.assert_frame_equal(
        python.observed_frame(), compiled.observed_frame(), check_exact=True,
    )
    assert python.digest() == compiled.digest()
