"""Grid convergence maths, checked against sequences with known answers."""

import math

import pytest

from analyse import gci


@pytest.mark.parametrize("p_true", [1.0, 1.5, 2.0, 2.5])
def test_recovers_a_known_order_and_limit(p_true):
    """Build f(h) = f_exact + C h^p exactly, then see if GCI inverts it.

    This is the only honest way to test an extrapolation: give it a sequence
    whose limit you already know.
    """
    f_exact, coeff, r = 0.00800, 0.5, 1.5
    h1 = 0.01
    f1 = f_exact + coeff * h1 ** p_true
    f2 = f_exact + coeff * (h1 * r) ** p_true
    f3 = f_exact + coeff * (h1 * r * r) ** p_true

    res = gci(f1, f2, f3, r)
    assert res["monotone"]
    assert res["p"] == pytest.approx(p_true, rel=1e-9)
    assert res["extrapolated"] == pytest.approx(f_exact, rel=1e-9)


def test_gci_shrinks_as_the_grid_refines():
    """A finer sequence must report a smaller uncertainty band."""
    f_exact, coeff, r, p = 0.008, 0.5, 1.5, 2.0
    coarse = [f_exact + coeff * (0.02 * r ** i) ** p for i in range(3)]
    fine = [f_exact + coeff * (0.005 * r ** i) ** p for i in range(3)]
    assert gci(*fine, r)["gci_fine"] < gci(*coarse, r)["gci_fine"]


def test_oscillatory_sequence_is_reported_not_extrapolated():
    """Richardson does not apply to a non-monotone sequence.

    Returning a number anyway would be the failure mode worth guarding: it
    looks like a converged answer and is not one.
    """
    res = gci(0.0100, 0.0090, 0.0095, 1.5)
    assert res["monotone"] is False
    assert res["extrapolated"] is None
    assert res["p"] is None


def test_identical_values_give_no_result():
    assert gci(0.008, 0.008, 0.008, 1.5) is None


def test_extrapolation_lies_beyond_the_finest_grid():
    """The limit must sit on the far side of the finest value, not between."""
    f_exact, coeff, r, p = 0.008, 0.5, 1.5, 2.0
    f = [f_exact + coeff * (0.01 * r ** i) ** p for i in range(3)]
    res = gci(*f, r)
    assert res["extrapolated"] < f[0] < f[1] < f[2]
