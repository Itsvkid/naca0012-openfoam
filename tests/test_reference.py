"""The reference-data loader and comparison.

Synthetic data throughout: these test the machinery, not any particular
experiment.
"""

import json

import pytest

from src.reference import Reference, compare, load

GOOD = {
    "source": "Abbott & Von Doenhoff, Theory of Wing Sections, Dover 1959",
    "figure": "Appendix IV, NACA 0012",
    "reynolds": 6.0e6,
    "surface_condition": "smooth",
    "digitised_by": "test",
    "points": [{"alpha": -4.0, "cl": -0.44, "cd": 0.0068},
               {"alpha": 0.0, "cl": 0.0, "cd": 0.0060},
               {"alpha": 4.0, "cl": 0.44, "cd": 0.0068},
               {"alpha": 8.0, "cl": 0.86, "cd": 0.0090}],
}


def _write(tmp_path, **overrides):
    data = {**GOOD, **overrides}
    p = tmp_path / "ref.json"
    p.write_text(json.dumps(data))
    return p


def test_loads_a_complete_dataset(tmp_path):
    ref = load(_write(tmp_path))
    assert ref.reynolds == 6.0e6
    assert len(ref.points) == 4


@pytest.mark.parametrize("field", ["source", "figure", "surface_condition",
                                   "digitised_by"])
def test_refuses_missing_provenance(tmp_path, field):
    """Data of unknown origin produces a plot that looks like validation."""
    with pytest.raises(ValueError, match="missing provenance"):
        load(_write(tmp_path, **{field: "TODO"}))


def test_refuses_an_empty_template(tmp_path):
    with pytest.raises(ValueError, match="template"):
        load(_write(tmp_path, points=[]))


def test_refuses_a_vague_surface_condition(tmp_path):
    """Smooth and roughened differ by ~40% in minimum drag."""
    with pytest.raises(ValueError, match="surface_condition"):
        load(_write(tmp_path, surface_condition="fairly clean"))


def test_interpolates_between_measured_points(tmp_path):
    ref = load(_write(tmp_path))
    assert ref.interpolate(2.0, "cl") == pytest.approx(0.22, rel=1e-9)


def test_never_extrapolates(tmp_path):
    """Reading past the end of a digitised curve invents data."""
    ref = load(_write(tmp_path))
    assert ref.interpolate(14.0, "cl") is None
    assert ref.interpolate(-10.0, "cl") is None


def test_comparison_reports_drag_in_counts(tmp_path):
    """A percentage on Cd 0.008 turns 2 counts into an apparent 25% failure."""
    ref = load(_write(tmp_path))
    rows = compare([{"alpha": 0.0, "cl": 0.0, "cd": 0.0062}], ref)
    assert rows[0]["cd_counts"] == pytest.approx(2.0, rel=1e-9)
    assert rows[0]["cl_delta"] == pytest.approx(0.0, abs=1e-12)


def test_points_outside_the_measured_range_are_flagged(tmp_path):
    ref = load(_write(tmp_path))
    rows = compare([{"alpha": 12.0, "cl": 1.17, "cd": 0.029}], ref)
    assert rows[0]["in_range"] is False
    assert rows[0]["cl_delta"] is None


def test_the_shipped_template_is_refused():
    """The template must stay unusable until someone fills it in."""
    with pytest.raises(ValueError):
        load("validation/reference/naca0012-template.json")
