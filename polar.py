"""Polars from the alpha sweep, and the checks that make them evidence.

The reference used here is **thin aerofoil theory**, dCl/dalpha = 2*pi per
radian. That is exact analysis, not remembered data, so it can be plotted
honestly. Experimental points from Abbott & Von Doenhoff are deliberately not
hard-coded: reconstructing a published curve from memory produces a comparison
that looks authoritative and is not checkable. Digitise the source and add it
here before presenting this as a validation against experiment.

What the theory does and does not bound:
  * the *slope* in the linear region — thin aerofoil theory is a good check,
    and a real 12% section sits a few percent below 2*pi
  * the zero-lift angle — zero for a symmetric section, exactly
  * stall — not predicted at all, and steady RANS does not resolve it either
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).parent
RESULTS = HERE / "validation" / "alpha-sweep.json"
TWO_PI_PER_DEG = 2.0 * math.pi * math.pi / 180.0     # 0.1097 per degree
# -4..+4 only. Widening to +8 pulls the slope from 0.1071 to 0.1057 and drags
# the fitted zero-lift angle from +0.0001 to +0.0119 degrees. Symmetry requires
# that angle to be exactly zero, so its drift is not a physical result — it is
# the fit reporting that it has swallowed points which are no longer on a
# straight line. The zero-lift angle is therefore the diagnostic for choosing
# the range, not a free output.
LINEAR_RANGE = (-4.0, 4.0)

THEMES = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink_muted": "#52514e",
              "grid": "#e4e3df", "series": ("#2a78d6", "#eb6834", "#1baf7a")},
    "dark": {"surface": "#101316", "ink": "#edeef0", "ink_muted": "#788087",
             "grid": "#23282e", "series": ("#3987e5", "#d95926", "#199e70")},
}


def lift_slope(points, lo=LINEAR_RANGE[0], hi=LINEAR_RANGE[1]):
    """Least-squares dCl/dalpha over the linear range, with the intercept.

    The intercept is the zero-lift angle check: a symmetric section must give
    zero, and a value that drifts from it is a symmetry error rather than an
    aerodynamic result.
    """
    pts = [(p["alpha"], p["cl"]) for p in points if lo <= p["alpha"] <= hi]
    if len(pts) < 3:
        return None
    n = len(pts)
    sx = sum(a for a, _ in pts)
    sy = sum(c for _, c in pts)
    sxx = sum(a * a for a, _ in pts)
    sxy = sum(a * c for a, c in pts)
    slope = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    intercept = (sy - slope * sx) / n
    return {"slope_per_deg": slope, "cl_at_zero": intercept,
            "alpha_zero_lift": -intercept / slope if slope else None,
            "vs_thin_aerofoil": slope / TWO_PI_PER_DEG, "n_points": n}


def _axes(t, title, xlabel, ylabel, size=(7.2, 4.6)):
    fig, ax = plt.subplots(figsize=size, dpi=160)
    fig.patch.set_facecolor(t["surface"])
    ax.set_facecolor(t["surface"])
    ax.set_title(title, color=t["ink"], fontsize=12, fontweight="600",
                 loc="left", pad=14)
    ax.set_xlabel(xlabel, color=t["ink_muted"], fontsize=10)
    ax.set_ylabel(ylabel, color=t["ink_muted"], fontsize=10)
    ax.grid(True, color=t["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["grid"])
    ax.tick_params(colors=t["ink_muted"], labelsize=9, length=0)
    return fig, ax


def figures(data, theme="light", suffix=""):
    t = THEMES[theme]
    pts = sorted(data["points"], key=lambda p: p["alpha"])
    a = [p["alpha"] for p in pts]
    cl = [p["cl"] for p in pts]
    cd = [p["cd"] for p in pts]
    fit = data.get("lift_slope")
    out = []

    # Cl vs alpha, against thin aerofoil theory
    fig, ax = _axes(t, "Lift curve — NACA 0012, Re 6×10⁶",
                    "Angle of attack α (deg)", "Lift coefficient $C_l$")
    ax.axhline(0, color=t["grid"], linewidth=1.0, zorder=1)
    ax.axvline(0, color=t["grid"], linewidth=1.0, zorder=1)
    theory = [TWO_PI_PER_DEG * x for x in a]
    ax.plot(a, theory, color=t["series"][1], linewidth=1.6,
            linestyle=(0, (5, 3)), zorder=3,
            label=f"Thin aerofoil theory  2π  ({TWO_PI_PER_DEG:.4f}/deg)")
    ax.plot(a, cl, "o-", color=t["series"][0], linewidth=2.0, markersize=7,
            zorder=4, label="Computed (k-ω SST)")
    if fit:
        ax.plot([], [], " ",
                label=f"Fitted slope {fit['slope_per_deg']:.4f}/deg  "
                      f"({fit['vs_thin_aerofoil']*100:.0f}% of 2π)")
    ax.legend(frameon=False, labelcolor=t["ink_muted"], fontsize=8.5,
              loc="upper left")
    p = HERE / "validation" / f"lift-curve{suffix}.png"
    fig.savefig(p, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)
    out.append(p)

    # Drag polar
    fig, ax = _axes(t, "Drag polar — NACA 0012, Re 6×10⁶",
                    "Drag coefficient $C_d$", "Lift coefficient $C_l$")
    ax.plot(cd, cl, "o-", color=t["series"][0], linewidth=2.0, markersize=7,
            zorder=4)
    for x, y, ang in zip(cd, cl, a):
        ax.annotate(f"{ang:+.0f}°", xy=(x, y), xytext=(7, -3),
                    textcoords="offset points", color=t["ink_muted"],
                    fontsize=8)
    ax.set_xlim(left=0)
    p = HERE / "validation" / f"drag-polar{suffix}.png"
    fig.savefig(p, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)
    out.append(p)

    # Efficiency
    fig, ax = _axes(t, "Aerodynamic efficiency — NACA 0012, Re 6×10⁶",
                    "Angle of attack α (deg)", "$C_l / C_d$")
    ax.plot(a, [c / d for c, d in zip(cl, cd)], "o-", color=t["series"][2],
            linewidth=2.0, markersize=7, zorder=4)
    best = max(zip(a, cl, cd), key=lambda r: r[1] / r[2])
    ax.annotate(f"best L/D {best[1]/best[2]:.1f} at {best[0]:+.0f}°",
                xy=(best[0], best[1] / best[2]), xytext=(8, -14),
                textcoords="offset points", color=t["ink"], fontsize=9)
    p = HERE / "validation" / f"efficiency{suffix}.png"
    fig.savefig(p, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)
    out.append(p)
    return out


def main() -> None:
    data = json.loads(RESULTS.read_text())
    pts = sorted(data["points"], key=lambda p: p["alpha"])

    print(f"{'alpha':>7} {'Cl':>9} {'Cd':>9} {'Cl/Cd':>8} {'y+max':>7} "
          f"{'Cl drift':>9}")
    for p in pts:
        print(f"{p['alpha']:>+7.1f} {p['cl']:>+9.4f} {p['cd']:>9.5f} "
              f"{p['cl']/p['cd']:>8.1f} {p['y_plus_max']:>7.2f} "
              f"{p['cl_drift']:>9.1e}")

    print("\nfit-range sensitivity (zero-lift angle must be 0 by symmetry):")
    for lo, hi in ((-4.0, 4.0), (-4.0, 6.0), (-4.0, 8.0), (-4.0, 10.0)):
        f = lift_slope(pts, lo, hi)
        if f:
            print(f"  {lo:+.0f}..{hi:+.0f}  slope {f['slope_per_deg']:.5f}/deg  "
                  f"{f['vs_thin_aerofoil']*100:5.1f}% of 2π   "
                  f"α₀ {f['alpha_zero_lift']:+.4f}°")

    strict = lift_slope(pts, *LINEAR_RANGE)
    if strict:
        print("\ndeparture from the linear fit:")
        for q in pts:
            if q["alpha"] > LINEAR_RANGE[1]:
                lin = strict["slope_per_deg"] * q["alpha"] + strict["cl_at_zero"]
                print(f"  α={q['alpha']:+5.1f}  linear {lin:.4f}  "
                      f"computed {q['cl']:.4f}  {100*(q['cl']-lin)/lin:+.1f}%")

    fit = lift_slope(pts)
    if fit:
        data["lift_slope"] = fit
        print(f"\nlift-curve slope over {LINEAR_RANGE[0]:+.0f}..{LINEAR_RANGE[1]:+.0f}°"
              f" ({fit['n_points']} points)")
        print(f"  computed          {fit['slope_per_deg']:.5f} /deg")
        print(f"  thin aerofoil 2π  {TWO_PI_PER_DEG:.5f} /deg")
        print(f"  ratio             {fit['vs_thin_aerofoil']*100:.1f}%")
        print(f"  zero-lift angle   {fit['alpha_zero_lift']:+.4f}°  "
              f"(symmetric section requires 0)")
        RESULTS.write_text(json.dumps(data, indent=2))

    for theme, suffix in (("light", ""), ("dark", "-dark")):
        for p in figures(data, theme, suffix):
            print("wrote", p.name)


if __name__ == "__main__":
    main()
