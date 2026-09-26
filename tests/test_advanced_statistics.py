import numpy as np
import pandas as pd
import pytest

from cinch.association.high_order import background_stat, make_cluster_indices, permuted_background, stable_seed
from cinch.network import aracne_prune
from cinch.population import conditional_mi, informative_clusters, phylo_permutation, weighted_mi_samples


def test_weighted_and_conditional_mi_hand_calculation():
    x = np.array([0, 0, 1, 1], dtype=bool)
    y = x.copy()
    weights = np.ones(4)
    assert weighted_mi_samples(x, y, weights) == pytest.approx(np.log(2))
    clusters = np.array([0, 0, 1, 1])
    assert conditional_mi(x, y, weights, clusters) == 0.0


def test_informative_clusters_and_permutation_are_seed_reproducible():
    x = np.tile(np.array([0, 0, 1, 1], dtype=bool), 2)
    y = x.copy()
    clusters = np.repeat([10, 20], 4)
    weights = np.ones(8)
    assert informative_clusters(x, y, clusters, 4) == ([10, 20], 2)
    one = phylo_permutation(x, y, weights, clusters, [10, 20], 39, np.random.default_rng(7))
    two = phylo_permutation(x, y, weights, clusters, [10, 20], 39, np.random.default_rng(7))
    assert one == two
    assert one[0] == pytest.approx(np.log(2))
    assert 0 < one[3] <= 1


def test_aracne_unique_minimum_deleted_but_ties_preserved():
    edges = pd.DataFrame({"gene1": ["a", "a", "b"], "gene2": ["b", "c", "c"], "score": [0.1, 0.8, 0.9]})
    assert aracne_prune(edges, "score").aracne_direct.tolist() == [False, True, True]
    tied = edges.assign(score=[0.1, 0.1, 0.9])
    assert aracne_prune(tied, "score").aracne_direct.all()


def test_high_order_delta_epidis_identity_and_cluster_permutation():
    # Two SHCs each contain both C states; A/B are dependent only when C=1.
    c_block = np.repeat([0, 1], 8)
    c = np.tile(c_block, 2).astype(bool)
    x0 = np.tile([0, 1], 4)
    y0 = np.tile([0, 0, 1, 1], 2)
    x1 = np.tile([0, 0, 1, 1], 2)
    y1 = x1.copy()
    x = np.tile(np.r_[x0, x1], 2).astype(bool)
    y = np.tile(np.r_[y0, y1], 2).astype(bool)
    shc = np.repeat([0, 1], 16)
    weights = np.ones(32)
    groups = make_cluster_indices(shc, 2)
    result = background_stat(x, y, c, weights, shc, min_group=8, min_shc_cell=4,
                             min_informative=2, cluster_indices=groups)
    assert result["eligible"]
    assert result["epidis_0"] ** 2 == pytest.approx(result["conditional_mi_0"] / np.log(2))
    assert result["epidis_1"] ** 2 == pytest.approx(result["conditional_mi_1"] / np.log(2))
    assert result["delta"] > 0
    assert stable_seed("block-7") == stable_seed("block-7")
    shuffled = permuted_background(c, groups, np.random.default_rng(2))
    for indices in groups:
        assert shuffled[indices].sum() == c[indices].sum()
