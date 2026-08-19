"""Lift-curve fitting, checked against synthetic data with known answers."""

import math

import pytest

from polar import TWO_PI_PER_DEG, lift_slope


def _points(slope, intercept, alphas):
    return [{"alpha": a, "cl": slope * a + intercept, "cd": 0.01}
            for a in alphas]


def test_thin_aerofoil_constant_is_two_pi_per_radian():
    """2*pi per radian expressed per degree."""
    assert TWO_PI_PER_DEG == pytest.approx(2 * math.pi * math.pi / 180.0)
    assert TWO_PI_PER_DEG == pytest.approx(0.10966, abs=1e-5)
    # Same thing the long way round, as a guard against a units slip.
    assert TWO_PI_PER_DEG == pytest.approx(2 * math.pi / (180.0 / math.pi))


@pytest.mark.parametrize("slope", [0.09, 0.105, TWO_PI_PER_DEG])
def test_recovers_a_known_slope(slope):
    fit = lift_slope(_points(slope, 0.0, [-4, -2, 0, 2, 4, 6, 8]))
    assert fit["slope_per_deg"] == pytest.approx(slope, rel=1e-12)


def test_symmetric_section_gives_zero_lift_at_zero_incidence():
    fit = lift_slope(_points(0.105, 0.0, [-4, -2, 0, 2, 4]))
    assert fit["cl_at_zero"] == pytest.approx(0.0, abs=1e-12)
    assert fit["alpha_zero_lift"] == pytest.approx(0.0, abs=1e-12)


def test_detects_a_shifted_zero_lift_angle():
    """A cambered section, or a broken symmetry, must show up here."""
    fit = lift_slope(_points(0.105, 0.21, [-4, -2, 0, 2, 4]))
    assert fit["alpha_zero_lift"] == pytest.approx(-2.0, rel=1e-9)


def test_ratio_against_thin_aerofoil_theory():
    fit = lift_slope(_points(TWO_PI_PER_DEG * 0.92, 0.0, [-4, 0, 4, 8]))
    assert fit["vs_thin_aerofoil"] == pytest.approx(0.92, rel=1e-12)


def test_only_the_linear_range_is_fitted():
    """Points past the requested range must not drag the slope down.

    The fit exists to measure the linear region; including stalled points would
    report a slope that describes neither the linear part nor the stall.
    """
    pts = _points(0.105, 0.0, [-4, -2, 0, 2, 4, 6, 8])
    pts += [{"alpha": 14.0, "cl": 0.9, "cd": 0.03},
            {"alpha": 16.0, "cl": 0.8, "cd": 0.05}]
    fit = lift_slope(pts, -4.0, 8.0)
    assert fit["n_points"] == 7
    assert fit["slope_per_deg"] == pytest.approx(0.105, rel=1e-12)


def test_too_few_points_returns_nothing():
    assert lift_slope(_points(0.1, 0.0, [0.0, 2.0])) is None
