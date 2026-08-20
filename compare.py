"""Compare the computed polar against digitised experimental data.

Usage:
    python compare.py validation/reference/<your-digitised-file>.json

Refuses to run against the template. Fill in a reference file first — see
src/reference.py for what provenance is required and why.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from polar import THEMES, _axes  # noqa: E402
from src.reference import compare, load  # noqa: E402

HERE = Path(__file__).parent
SWEEP = HERE / "validation" / "alpha-sweep.json"


def figure(rows, ref, theme="light", suffix=""):
    t = THEMES[theme]
    have = [r for r in rows if r["cl_ref"] is not None]
    fig, ax = _axes(t, f"Lift curve vs experiment — {ref.surface_condition}, "
                       f"Re {ref.reynolds:.0e}",
                    "Angle of attack α (deg)", "Lift coefficient $C_l$")
    ax.plot([r["alpha"] for r in rows], [r["cl"] for r in rows], "o-",
            color=t["series"][0], linewidth=2.0, markersize=7, zorder=4,
            label="Computed (k-ω SST)")
    if have:
        ax.plot([r["alpha"] for r in have], [r["cl_ref"] for r in have], "s--",
                color=t["series"][1], linewidth=1.8, markersize=6, zorder=3,
                label=f"{ref.source[:40]}")
    ax.legend(frameon=False, labelcolor=t["ink_muted"], fontsize=8.5,
              loc="upper left")
    p = HERE / "validation" / f"vs-experiment{suffix}.png"
    fig.savefig(p, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)
    return p


def main() -> None:
    if len(sys.argv) < 2:
        print(__doc__)
        raise SystemExit(2)

    ref = load(sys.argv[1])
    computed = json.loads(SWEEP.read_text())["points"]
    rows = compare(computed, ref)

    print(f"source   {ref.source}")
    print(f"figure   {ref.figure}")
    print(f"surface  {ref.surface_condition}   Re {ref.reynolds:.3g}")
    print(f"digitised by {ref.digitised_by}\n")

    print(f"{'alpha':>7} {'Cl':>9} {'Cl exp':>9} {'ΔCl':>8} "
          f"{'Cd':>9} {'Cd exp':>9} {'Δ counts':>9}")
    for r in rows:
        def fmt(v, width, prec):
            return f"{v:>+{width}.{prec}f}" if v is not None else " " * width
        print(f"{r['alpha']:>+7.1f} {r['cl']:>+9.4f} "
              f"{fmt(r['cl_ref'], 9, 4)} {fmt(r['cl_delta'], 8, 4)} "
              f"{r['cd']:>9.5f} {fmt(r['cd_ref'], 9, 5)} "
              f"{fmt(r['cd_counts'], 9, 1)}")

    inside = [r for r in rows if r["in_range"]]
    outside = [r for r in rows if not r["in_range"]]
    if outside:
        angles = ", ".join(f"{r['alpha']:+.0f}deg" for r in outside)
        print(f"\n{len(outside)} computed angles lie outside the measured "
              f"range and are not compared: {angles}")
    if inside:
        cl_d = [abs(r["cl_delta"]) for r in inside if r["cl_delta"] is not None]
        cd_c = [abs(r["cd_counts"]) for r in inside if r["cd_counts"] is not None]
        if cl_d:
            print(f"mean |ΔCl|        {sum(cl_d)/len(cl_d):.4f}")
        if cd_c:
            print(f"mean |ΔCd|        {sum(cd_c)/len(cd_c):.1f} counts")

    for theme, suffix in (("light", ""), ("dark", "-dark")):
        print("wrote", figure(rows, ref, theme, suffix).name)


if __name__ == "__main__":
    main()
