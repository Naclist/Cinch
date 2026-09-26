import numpy as np
import pytest

from cinch.profiles import LegacyMissingPolicy, StateProfile, from_legacy_alleles, load_profile, save_profile


def test_state_profile_rejects_missingness_collapse_and_shape_errors():
    with pytest.raises(ValueError, match="callable type requires"):
        StateProfile(np.array(["s1"]), np.array(["a"]), np.array([[0]]), np.array([[0]]))
    with pytest.raises(ValueError, match="shape"):
        StateProfile(np.array(["s1"]), np.array(["a"]), np.array([[1, 1]]), np.array([[0]]))
    with pytest.raises(ValueError, match=r"\{-1, 0, 1\}"):
        StateProfile(np.array(["s1"]), np.array(["a"]), np.array([[2]]), np.array([[-1]]))


def test_legacy_adapter_requires_and_records_missing_policy():
    matrix = np.array([[12, 0], [7, -1], [12, np.nan]], dtype=object)
    samples, loci = np.array(["b", "a", "c"]), np.array(["L1", "L2"])
    absent = from_legacy_alleles(matrix, samples, loci, missing_policy=LegacyMissingPolicy.ABSENT)
    unresolved = from_legacy_alleles(matrix, samples, loci, missing_policy="unresolved")
    assert absent.presence.tolist() == [[1, 0], [1, 0], [1, 0]]
    assert unresolved.presence.tolist() == [[1, -1], [1, -1], [1, -1]]
    # Codes are deterministic by source value, not row encounter order.
    assert absent.types[:, 0].tolist() == [1, 0, 1]
    assert absent.metadata["legacy_missing_policy"] == "absence"
    with pytest.raises(ValueError, match="missing_policy"):
        from_legacy_alleles(matrix, samples, loci, missing_policy="guess")

    strings = from_legacy_alleles(
        np.array([["allele_A"], ["NA"], ["0"]], dtype=object), samples, np.array(["L1"]),
        missing_policy="unresolved",
    )
    assert strings.presence[:, 0].tolist() == [1, -1, -1]


def test_profile_round_trip_and_frozen_key_compatibility(tmp_path):
    profile = StateProfile(
        np.array(["s2", "s1"]), np.array(["b", "a"]),
        np.array([[1, 0], [-1, 1]], dtype=np.int8),
        np.array([[3, -1], [-1, 2]], dtype=np.int32),
        {"schema": "CINCH_STATE_PROFILE_V2", "note": "roundtrip"},
    )
    path = save_profile(profile, tmp_path / "profile.npz")
    loaded = load_profile(path)
    np.testing.assert_array_equal(loaded.samples, profile.samples)
    np.testing.assert_array_equal(loaded.loci, profile.loci)
    np.testing.assert_array_equal(loaded.presence, profile.presence)
    np.testing.assert_array_equal(loaded.types, profile.types)
    assert loaded.metadata == profile.metadata

    frozen = tmp_path / "frozen.npz"
    np.savez_compressed(frozen, samples=profile.samples, loci=profile.loci,
                        presence=profile.presence, type_state=profile.types)
    compatible = load_profile(frozen)
    np.testing.assert_array_equal(compatible.types, profile.types)
    assert compatible.metadata["schema"] == "CINCH_FROZEN_V1_STATE"
