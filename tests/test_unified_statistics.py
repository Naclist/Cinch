from __future__ import annotations

import math

import numpy as np
import pytest

from cinch.frozen_v1.statistics import contingency
from cinch.statistics import bh_adjust, directional_epidis, weighted_binary_information


def scalar_binary_mi(c11: float, mx: float, my: float) -> float:
    cells = np.array([1.0 - mx - my + c11, my - c11, mx - c11, c11])
    px = np.array([cells[0] + cells[1], cells[2] + cells[3]])
    py = np.array([cells[0] + cells[2], cells[1] + cells[3]])
    expected = np.array([px[0] * py[0], px[0] * py[1], px[1] * py[0], px[1] * py[1]])
    return float(sum(p * math.log(p / e) for p, e in zip(cells, expected) if p > 0))


def test_weighted_binary_matches_independent_scalar_reference():
    c11 = np.array([0.25, 0.50, 0.00, 0.18])
    mx = np.array([0.50, 0.50, 0.50, 0.30])
    my = np.array([0.50, 0.50, 0.50, 0.60])
    mi, nmi, hx, hy = weighted_binary_information(c11, mx, my)
    expected = np.array([scalar_binary_mi(a, b, c) for a, b, c in zip(c11, mx, my)])
    np.testing.assert_allclose(mi, expected, atol=1e-15, rtol=0)
    assert nmi[0] == pytest.approx(0.0, abs=1e-15)
    assert nmi[1] == pytest.approx(1.0, abs=1e-15)
    assert np.all(np.isfinite(hx)) and np.all(np.isfinite(hy))


def test_unweighted_frozen_mi_is_not_silently_replaced():
    x = np.array([0, 0, 1, 1], dtype=np.int32)
    y = np.array([0, 1, 0, 1], dtype=np.int32)
    frozen = contingency(x, y)
    weighted, *_ = weighted_binary_information(np.array([0.25]), np.array([0.5]), np.array([0.5]))
    assert frozen is not None
    assert frozen["MIraw"] == pytest.approx(weighted[0], abs=1e-15)
    assert contingency is not weighted_binary_information


def test_epidis_squared_identity_is_limited_to_frozen_construction():
    c11 = np.array([0.25, 0.50, 0.08, 0.18])
    mx = np.array([0.50, 0.50, 0.20, 0.30])
    my = np.array([0.50, 0.50, 0.40, 0.60])
    mi, *_ = weighted_binary_information(c11, mx, my)
    epidis = directional_epidis(c11, mx, my)
    np.testing.assert_allclose(epidis**2, mi / np.log(2.0), atol=1e-14, rtol=0)


def test_epidis_rejects_shape_mismatch_and_marks_constant_conditioner():
    with pytest.raises(ValueError):
        directional_epidis(np.array([0.1]), np.array([0.2, 0.3]), np.array([0.4]))
    result = directional_epidis(np.array([0.0, 1.0]), np.array([0.0, 1.0]), np.array([0.2, 1.0]))
    assert np.isnan(result).all()


def test_bh_matches_hand_calculation_and_validates_domain():
    p = np.array([0.01, 0.04, 0.03, 0.002, 0.04])
    np.testing.assert_allclose(bh_adjust(p), [0.025, 0.04, 0.04, 0.01, 0.04], atol=0, rtol=0)
    with pytest.raises(ValueError):
        bh_adjust(np.array([np.nan]))
    with pytest.raises(ValueError):
        bh_adjust(np.array([-0.1]))
