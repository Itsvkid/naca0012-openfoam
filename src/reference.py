"""Loading published experimental data for comparison.

The loader refuses a dataset that does not say where it came from. That is not
bureaucracy: a polar compared against numbers of unknown provenance looks like
validation and is not one, and by the time anyone asks "which Reynolds number,
which surface condition, which figure?" the plot is already in a report.

Required provenance
-------------------
source              full citation, enough to find the figure again
figure              the specific figure or table the points were read from
reynolds            the experiment's Reynolds number, not the computation's
surface_condition   "smooth" or "standard roughness" — they differ by more than
                    the discrepancy most CFD comparisons are trying to resolve
digitised_by        who read the points off, and how it was checked

Optional but worth carrying
---------------------------
source_sha256       checksum of the exact document the points came from, so a
                    reader can confirm they are looking at the same edition
                    rather than a different scan with different pagination
l_over_d            the source's own lift-to-drag column, if it prints one —
                    the loader then verifies Cl/Cd reproduces it

A NACA 0012 at Re 6e6 has a drag coefficient near 0.0060 smooth and near 0.0085
with standard roughness. Comparing a computation against the wrong one moves the
answer by 40%, which is larger than every numerical effect in this project put
together.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

REQUIRED = ("source", "figure", "reynolds", "surface_condition", "digitised_by")
PLACEHOLDERS = {"", "TODO", "FIXME", "TBD", "unknown", "?"}


@dataclass(frozen=True)
class Reference:
    """Digitised experimental points with the provenance to justify them."""

    source: str
    figure: str
    reynolds: float
    surface_condition: str
    digitised_by: str
    points: list          # [{"alpha": deg, "cl": ..., "cd": ...}, ...]
    notes: str = ""

    @property
    def alphas(self) -> list:
        return [p["alpha"] for p in self.points]

    def interpolate(self, alpha: float, key: str) -> float | None:
        """Linear interpolation onto `alpha`. None outside the measured range.

        Deliberately does not extrapolate. Reading past the end of a digitised
        curve invents data, and it does so most eagerly near stall where the
        curve is least linear and the temptation is greatest.
        """
        pts = sorted((p for p in self.points if key in p),
                     key=lambda p: p["alpha"])
        if len(pts) < 2 or not (pts[0]["alpha"] <= alpha <= pts[-1]["alpha"]):
            return None
        for a, b in zip(pts, pts[1:]):
            if a["alpha"] <= alpha <= b["alpha"]:
                span = b["alpha"] - a["alpha"]
                if span == 0:
                    return a[key]
                t = (alpha - a["alpha"]) / span
                return a[key] + t * (b[key] - a[key])
        return None


def load(path) -> Reference:
    """Load a reference dataset, refusing anything without real provenance."""
    path = Path(path)
    raw = json.loads(path.read_text())

    missing = [f for f in REQUIRED
               if str(raw.get(f, "")).strip() in PLACEHOLDERS or f not in raw]
    if missing:
        raise ValueError(
            f"{path.name} is missing provenance for {missing}. A comparison "
            f"against data of unknown origin is not a validation. Fill these in "
            f"from the source you actually read the points off."
        )

    points = raw.get("points") or []
    if len(points) < 2:
        raise ValueError(
            f"{path.name} has {len(points)} points. This is the template — "
            f"digitise the figure named in it and fill in the points."
        )
    for p in points:
        if "alpha" not in p or not ({"cl", "cd"} & set(p)):
            raise ValueError(f"{path.name}: every point needs alpha and cl or cd")

    # If the source printed its own lift-to-drag column, use it. Cl/Cd must
    # reproduce it, and a transposed or misread digit in either coefficient
    # breaks the ratio even though both numbers still look plausible on their
    # own. This is the cheapest check there is on a digitised table, and it
    # catches the error that matters most.
    for p in points:
        if {"cl", "cd", "l_over_d"} <= set(p) and p["cd"]:
            if abs(p["cl"] / p["cd"] - p["l_over_d"]) > 0.1:
                raise ValueError(
                    f"{path.name}: at alpha {p['alpha']}, Cl/Cd = "
                    f"{p['cl'] / p['cd']:.2f} but the source prints "
                    f"{p['l_over_d']:.2f}. One of the three was misread."
                )

    if str(raw["surface_condition"]).lower() not in ("smooth", "standard roughness"):
        raise ValueError(
            f"surface_condition {raw['surface_condition']!r} is not one of "
            f"'smooth' or 'standard roughness'. They differ by roughly 40% in "
            f"minimum drag, so the distinction cannot be left vague."
        )

    return Reference(
        source=raw["source"], figure=raw["figure"],
        reynolds=float(raw["reynolds"]),
        surface_condition=raw["surface_condition"],
        digitised_by=raw["digitised_by"], points=points,
        notes=raw.get("notes", ""))


def compare(computed: list, reference: Reference) -> list:
    """Match computed points to interpolated reference values.

    Cl is compared as an absolute difference and Cd in drag counts (1e-4),
    because a percentage on a drag coefficient of 0.008 makes a 2-count
    difference look like a 25% failure when 2 counts is close to the digitising
    error of reading a printed figure.
    """
    rows = []
    for c in sorted(computed, key=lambda p: p["alpha"]):
        a = c["alpha"]
        cl_ref, cd_ref = reference.interpolate(a, "cl"), reference.interpolate(a, "cd")
        rows.append({
            "alpha": a,
            "cl": c["cl"], "cl_ref": cl_ref,
            "cl_delta": None if cl_ref is None else c["cl"] - cl_ref,
            "cd": c["cd"], "cd_ref": cd_ref,
            "cd_counts": None if cd_ref is None else (c["cd"] - cd_ref) * 1e4,
            "in_range": cl_ref is not None or cd_ref is not None,
        })
    return rows
