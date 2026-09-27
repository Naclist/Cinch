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


def test_numba_block_engine_matches_python_with_float_tolerance():
    pytest.importorskip("numba")
    profile, distances = fixture_profile()
    python_blocks, _ = iter_scored_blocks(
        profile, distances, minimum_informative=8, minimum_state_count=2,
        block_pairs=4, engine="python",
    )
    numba_blocks, _ = iter_scored_blocks(
        profile, distances, minimum_informative=8, minimum_state_count=2,
        block_pairs=4, engine="numba",
    )
    for (python_block, python_frames), (numba_block, numba_frames) in zip(python_blocks, numba_blocks):
        assert python_block == numba_block
        for channel in CHANNELS:
            pd.testing.assert_frame_equal(
                python_frames[channel], numba_frames[channel],
                check_exact=False, rtol=1e-12, atol=1e-12,
            )


@pytest.mark.parametrize("seed,minimum_state_count", [(7, 1), (91, 2), (1301, 3)])
def test_numba_block_engine_matches_random_multiallelic_profiles(seed, minimum_state_count):
    pytest.importorskip("numba")
    rng = np.random.default_rng(seed)
    presence = rng.integers(0, 2, size=(36, 8), dtype=np.int8)
    presence[rng.random(presence.shape) < 0.12] = -1
    types = np.full(presence.shape, -1, dtype=np.int32)
    for locus in range(types.shape[1]):
        called = presence[:, locus] == 1
        types[called, locus] = rng.integers(0, locus + 2, size=int(called.sum()))
    profile = StateProfile(
        np.array([f"s{i}" for i in range(len(presence))]),
        np.array([f"L{i}" for i in range(presence.shape[1])]),
        presence,
        types,
    )
    distances = pd.DataFrame([
        {"i": i, "j": j, "order_distance": float(i + j), "bp_distance": float(j - i)}
        for i, j in itertools.combinations(range(types.shape[1]), 2)
    ])
    python_blocks, _ = iter_scored_blocks(
        profile, distances, minimum_informative=6,
        minimum_state_count=minimum_state_count, block_pairs=7, engine="python",
    )
    numba_blocks, _ = iter_scored_blocks(
        profile, distances, minimum_informative=6,
        minimum_state_count=minimum_state_count, block_pairs=7, engine="numba",
    )
    for (_, python_frames), (_, numba_frames) in zip(python_blocks, numba_blocks):
        for channel in CHANNELS:
            pd.testing.assert_frame_equal(
                python_frames[channel], numba_frames[channel],
                check_exact=False, rtol=1e-12, atol=1e-12,
            )


def test_numba_driver_ignores_pairwise_absent_global_type_states():
    pytest.importorskip("numba")
    presence = np.array([
        [1, 0], [1, 0], [1, 0], [1, 0],
        [1, 1], [1, 1], [1, 1], [1, 1],
    ], dtype=np.int8)
    types = np.array([
        [0, -1], [0, -1], [0, -1], [0, -1],
        [1, 0], [1, 1], [2, 0], [2, 1],
    ], dtype=np.int32)
    profile = StateProfile(
        np.array([f"s{i}" for i in range(8)]), np.array(["left", "right"]),
        presence, types,
    )
    distances = pd.DataFrame([{"i": 0, "j": 1, "order_distance": 1.0, "bp_distance": 100.0}])
    python_blocks, _ = iter_scored_blocks(
        profile, distances, minimum_informative=4, minimum_state_count=1,
        block_pairs=1, engine="python",
    )
    numba_blocks, _ = iter_scored_blocks(
        profile, distances, minimum_informative=4, minimum_state_count=1,
        block_pairs=1, engine="numba",
    )
    python_frames = next(python_blocks)[1]
    numba_frames = next(numba_blocks)[1]
    pd.testing.assert_frame_equal(
        python_frames["TT"], numba_frames["TT"],
        check_exact=False, rtol=1e-12, atol=1e-12,
    )
    assert int(numba_frames["TT"].iloc[0].enriched_driver_A_state) == 1


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
