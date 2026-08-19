"""The C-grid generator.

These run without Docker: they check the dictionary the generator emits, which
is where the mesh's correctness is decided. Whether OpenFOAM then likes it is
checkMesh's job, and that is recorded in the README.
"""

import pytest

from src.airfoil import NACA4
from src.cgrid import CGrid

G = CGrid()


def test_every_block_is_anticlockwise():
    """blockMesh rejects a clockwise base face as inside-out.

    to_dict raises on a clockwise block, so generating at all is the assertion.
    """
    assert G.to_dict().count("hex") == 6


def test_surface_segments_are_exact_mirrors():
    """A symmetric section must produce a mesh symmetric about y = 0.

    Splitting the surface by comparing x against 0.5*chord silently broke this:
    cosine spacing puts the mid-chord point at 0.49999999999999994, so both
    `<= 0.5` and `< 0.5` accept it and the block corner reappears as an
    interior polyLine point — on one surface and not its mirror.
    """
    upper, lower = G.surface_points(True), G.surface_points(False)
    iu, il = G._mid_index(upper), G._mid_index(lower)
    assert iu == il
    fwd_u, fwd_l = list(reversed(upper[1:iu])), lower[1:il]
    aft_u, aft_l = list(reversed(upper[iu + 1:-1])), lower[il + 1:-1]
    assert len(fwd_u) == len(fwd_l) and len(aft_u) == len(aft_l)
    for (xu, zu), (xl, zl) in zip(fwd_u, reversed(fwd_l)):
        assert xu == pytest.approx(xl, abs=1e-15)
        assert zu == pytest.approx(-zl, abs=1e-15)


def test_mid_index_is_not_at_an_endpoint():
    upper = G.surface_points(True)
    assert 0 < G._mid_index(upper) < len(upper) - 1


def test_first_cell_height_matches_the_geometric_series():
    """Sum of the graded radial cells must equal the block's radial length."""
    n, ratio = G.n_normal, G.expansion_normal ** (1.0 / (G.n_normal - 1))
    total = G.first_cell_height() * (ratio ** n - 1.0) / (ratio - 1.0)
    assert total == pytest.approx(G.radius, rel=1e-9)


@pytest.mark.parametrize("target", [0.5, 1.0, 2.0])
def test_expansion_solves_for_the_requested_y_plus(target):
    """The solver and the estimator must agree, or one of them is wrong."""
    e = G.expansion_for_y_plus(target, 6.0e6)
    tuned = CGrid(expansion_normal=e)
    assert tuned.y_plus_estimate(6.0e6) == pytest.approx(target, rel=1e-3)


def test_finer_y_plus_needs_stronger_grading():
    assert (G.expansion_for_y_plus(0.5, 6.0e6)
            > G.expansion_for_y_plus(2.0, 6.0e6))


def test_dictionary_declares_every_patch():
    text = G.to_dict()
    for patch in ("aerofoil", "inlet", "outlet", "topAndBottom", "frontAndBack"):
        assert patch in text
    assert "type empty;" in text          # 2D
    assert "type wall;" in text           # the aerofoil


def test_rejects_a_far_field_that_is_too_close():
    with pytest.raises(ValueError, match="too close"):
        CGrid(far_radius=1.5)


def test_rejects_an_outlet_inside_the_arc():
    with pytest.raises(ValueError, match="wake_length must exceed"):
        CGrid(far_radius=12.0, wake_length=8.0)


def test_rejects_a_useless_resolution():
    with pytest.raises(ValueError, match="too coarse"):
        CGrid(n_normal=4)


def test_cambered_section_is_not_symmetric():
    """Guards the mirror test above from passing vacuously."""
    g = CGrid(section=NACA4.from_code("4412"))
    upper, lower = g.surface_points(True), g.surface_points(False)
    assert any(abs(zu + zl) > 1e-6 for (_, zu), (_, zl) in zip(upper, lower))
