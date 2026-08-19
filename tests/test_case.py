"""Flow conditions and the dictionaries derived from them."""

import math

import pytest

from src.case import Case


def test_reynolds_sets_viscosity():
    c = Case(reynolds=6.0e6, velocity=1.0, chord=1.0)
    assert c.nu == pytest.approx(1.0 / 6.0e6, rel=1e-12)
    assert c.velocity * c.chord / c.nu == pytest.approx(6.0e6, rel=1e-12)


@pytest.mark.parametrize("alpha", [0.0, 4.0, -6.0, 12.0])
def test_lift_and_drag_directions_are_orthonormal(alpha):
    """A polar computed against non-orthogonal reference axes is meaningless."""
    c = Case(alpha_deg=alpha)
    d, l = c.flow_direction, c.lift_direction
    assert d[0] * l[0] + d[1] * l[1] == pytest.approx(0.0, abs=1e-15)
    assert math.hypot(*d) == pytest.approx(1.0, rel=1e-15)
    assert math.hypot(*l) == pytest.approx(1.0, rel=1e-15)


def test_zero_incidence_aligns_with_the_x_axis():
    c = Case(alpha_deg=0.0)
    assert c.flow_direction == pytest.approx((1.0, 0.0), abs=1e-15)
    assert c.lift_direction == pytest.approx((0.0, 1.0), abs=1e-15)


def test_positive_alpha_tilts_the_flow_upward():
    c = Case(alpha_deg=10.0)
    assert c.flow_direction[1] > 0
    assert c.lift_direction[0] < 0


def test_freestream_eddy_viscosity_stays_below_the_molecular_one():
    """Menter's guidance for external aerodynamics.

    A freestream nut comparable to nu seeds turbulence into the boundary layer
    and the transition location becomes an artefact of the inlet condition.
    """
    c = Case()
    assert c.nut_inf < c.nu
    assert c.k_inf == pytest.approx(c.nut_inf * c.omega_inf, rel=1e-12)


def test_case_writes_every_required_file(tmp_path):
    written = Case().write(tmp_path)
    names = {p.relative_to(tmp_path).as_posix() for p in written}
    assert {"0/U", "0/p", "0/k", "0/omega", "0/nut",
            "constant/transportProperties", "constant/turbulenceProperties",
            "system/fvSchemes", "system/fvSolution",
            "system/controlDict"} <= names


def test_gradient_scheme_is_unlimited(tmp_path):
    """cellLimited broke symmetry: Cl = 0.049 at zero incidence, Cd 69% high.

    Pinned in a test because 'add a limiter' is the standard reflex when a run
    misbehaves, and doing it here silently reintroduces both errors.
    """
    Case().write(tmp_path)
    schemes = (tmp_path / "system" / "fvSchemes").read_text()
    assert "gradSchemes     { default Gauss linear; }" in schemes
    assert "cellLimited" not in schemes.split("gradSchemes")[1][:60]


def test_force_coefficients_use_the_flow_axes(tmp_path):
    Case(alpha_deg=8.0).write(tmp_path)
    control = (tmp_path / "system" / "controlDict").read_text()
    c = Case(alpha_deg=8.0)
    assert f"{c.lift_direction[0]:.10g}" in control
    assert f"{c.flow_direction[0]:.10g}" in control
    assert "patches         (aerofoil);" in control
