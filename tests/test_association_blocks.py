import itertools

import numpy as np
import pandas as pd
import pytest

from cinch.association import iter_pair_blocks, iter_scored_blocks, write_association_blocks
from cinch.frozen_v1.statistics import CHANNELS, score_all_pairs
from cinch.profiles import StateProfile


def fixture_profile():
    rng = np.random.default_rng(90210)
    presence = rng.integers(0, 2, size=(40, 6), dtype=np.int8)
    presence[rng.random(presence.shape) < 0.08] = -1
    types = rng.integers(0, 4, size=presence.shape, dtype=np.int32)
    types[presence != 1] = -1
    samples = np.array([f"s{i}" for i in range(len(presence))])
    loci = np.array([f"L{i}" for i in range(presence.shape[1])])
    rows = [
        {"i": i, "j": j, "order_distance": float(i + j), "bp_distance": float((j - i) * 100)}
        for i, j in itertools.combinations(range(len(loci)), 2)
    ]
    return StateProfile(samples, loci, presence, types), pd.DataFrame(rows)


def test_pair_blocks_cover_triangular_universe_once():
    blocks = list(iter_pair_blocks(8, 5))
    pairs = [pair for block in blocks for pair in block.pairs]
    assert len(blocks) == 6
    assert pairs == list(itertools.combinations(range(8), 2))
    assert all(len(block.pairs) <= 5 for block in blocks)


def test_block_engine_is_exactly_equal_to_frozen_small_case():
    profile, distances = fixture_profile()
    old, old_audit = score_all_pairs(
        profile.presence, profile.types, profile.loci, distances,
        minimum_informative=8, minimum_state_count=2,
    )
    iterator, new_audit = iter_scored_blocks(
        profile, distances, minimum_informative=8, minimum_state_count=2, block_pairs=4,
    )
    parts = {channel: [] for channel in CHANNELS}
    for _, frames in iterator:
        for channel in CHANNELS:
            if len(frames[channel]):
                parts[channel].append(frames[channel])
    for channel in CHANNELS:
        new = pd.concat(parts[channel], ignore_index=True) if parts[channel] else pd.DataFrame()
        pd.testing.assert_frame_equal(new, old[channel], check_exact=True)
    pd.testing.assert_frame_equal(new_audit, old_audit, check_exact=True)


def test_block_outputs_resume_and_reject_input_drift(tmp_path):
    profile, distances = fixture_profile()
    first = write_association_blocks(
        profile, distances, tmp_path, minimum_informative=8,
        minimum_state_count=2, block_pairs=4,
    )
    second = write_association_blocks(
        profile, distances, tmp_path, minimum_informative=8,
        minimum_state_count=2, block_pairs=4,
    )
    assert first["completed_blocks_this_run"] == 4
    assert second["completed_blocks_this_run"] == 0
    assert second["resumed_blocks_this_run"] == 4
    drifted = distances.copy()
    drifted.loc[0, "order_distance"] += 1
    with pytest.raises(ValueError, match="different inputs"):
        write_association_blocks(
            profile, drifted, tmp_path, minimum_informative=8,
            minimum_state_count=2, block_pairs=4,
        )
