"""Grid convergence analysis and figures for the mesh independence study.

Implements Roache's Grid Convergence Index. The observed order of convergence
is computed from the results rather than assumed to be the scheme's formal
order — the two differ whenever the meshes are not yet in the asymptotic range,
and a study that assumes p = 2 hides exactly that.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).parent

# Abbott & Von Doenhoff, NACA 0012, Re 6e6, alpha 0, smooth surface.
# Carried in code, not transcribed from the source — verify before publishing.
PUBLISHED_CD = 0.0080

THEMES = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink_muted": "#52514e",
              "grid": "#e4e3df", "series": ("#2a78d6", "#eb6834", "#1baf7a")},
    "dark": {"surface": "#101316", "ink": "#edeef0", "ink_muted": "#788087",
             "grid": "#23282e", "series": ("#3987e5", "#d95926", "#199e70")},
}


def gci(f1, f2, f3, r):
    """Roache GCI from three grids, finest first.

    Returns observed order p, the Richardson-extrapolated value, and the GCI on
    the finest grid as a fraction.
    """
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0:
        return None
    ratio = e32 / e21
    if ratio <= 0:
        # Oscillatory convergence: the sequence is not monotone, so Richardson
        # extrapolation does not apply and saying so is the honest result.
        return {"p": None, "extrapolated": None, "gci_fine": None,
                "monotone": False}
    p = math.log(abs(ratio)) / math.log(r)
    f_ext = f1 + (f1 - f2) / (r ** p - 1.0)
    gci_fine = 1.25 * abs(e21 / f1) / (r ** p - 1.0)
    return {"p": p, "extrapolated": f_ext, "gci_fine": gci_fine,
            "monotone": True}


def figure(data, theme="light", suffix=""):
    t = THEMES[theme]
    levels = data["levels"]
    h = [1.0 / math.sqrt(l["cells"]) for l in levels]      # representative size
    cd = [l["cd"] for l in levels]

    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=160)
    fig.patch.set_facecolor(t["surface"])
    ax.set_facecolor(t["surface"])
    ax.set_title("Grid convergence — NACA 0012, Re 6×10⁶, α = 0°",
                 color=t["ink"], fontsize=12, fontweight="600", loc="left",
                 pad=14)
    ax.set_xlabel("Representative cell size  h ∝ 1/√N", color=t["ink_muted"],
                  fontsize=10)
    ax.set_ylabel("Drag coefficient  $C_d$", color=t["ink_muted"], fontsize=10)
    ax.grid(True, color=t["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(t["grid"])
    ax.tick_params(colors=t["ink_muted"], labelsize=9, length=0)

    ax.plot(h, cd, "o-", color=t["series"][0], linewidth=2.0, markersize=8,
            zorder=4, label="Computed")
    for hi, ci, l in zip(h, cd, levels):
        ax.annotate(f"{l['label']}  {l['cells']:,}", xy=(hi, ci),
                    xytext=(6, 8), textcoords="offset points",
                    color=t["ink"], fontsize=8)

    res = data.get("gci")
    if res and res.get("extrapolated") is not None:
        ax.axhline(res["extrapolated"], color=t["series"][2], linewidth=1.4,
                   linestyle=(0, (5, 3)), zorder=3,
                   label=f"Richardson h→0: {res['extrapolated']:.5f}")
        ax.plot([0], [res["extrapolated"]], "o", color=t["series"][2],
                markersize=8, zorder=5)

    ax.axhline(PUBLISHED_CD, color=t["series"][1], linewidth=1.4,
               linestyle=(0, (2, 2)), zorder=3,
               label=f"Published: {PUBLISHED_CD:.4f}")

    # Room on the right for the coarsest point's label, and the legend moved
    # off the two horizontal reference lines it was sitting on.
    ax.set_xlim(0, max(h) * 1.18)
    ax.legend(frameon=False, labelcolor=t["ink_muted"], fontsize=9,
              loc="upper left")

    out = HERE / "validation" / f"grid-convergence{suffix}.png"
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, facecolor=t["surface"], bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    data = json.loads((HERE / "validation" / "mesh-study.json").read_text())
    levels = data["levels"]
    r = data["ratio"]

    print(f"{'':4} {'cells':>8} {'first cell':>12} {'y+ max':>8} "
          f"{'Cd':>10} {'Cl':>11} {'drift':>9}")
    for l in levels:
        print(f"{l['label']:<4} {l['cells']:>8,} {l['first_cell']:>12.3e} "
              f"{l['y_plus_max']:>8.3f} {l['cd']:>10.6f} {l['cl']:>+11.2e} "
              f"{l['cd_spread_last500']:>9.1e}")

    if len(levels) >= 3:
        f1, f2, f3 = (levels[-1]["cd"], levels[-2]["cd"], levels[-3]["cd"])
        res = gci(f1, f2, f3, r)
        data["gci"] = res
        print()
        if res and res["monotone"]:
            print(f"observed order of convergence p = {res['p']:.2f}")
            print(f"Richardson extrapolation to h=0: Cd = {res['extrapolated']:.6f}")
            print(f"GCI on the finest grid          = {res['gci_fine']*100:.2f}%")
            print(f"published Cd                    = {PUBLISHED_CD:.4f}")
            err = (res["extrapolated"] - PUBLISHED_CD) / PUBLISHED_CD
            print(f"extrapolated vs published       = {err*100:+.1f}%")
        else:
            print("convergence is not monotone — Richardson extrapolation "
                  "does not apply to this sequence")

    (HERE / "validation" / "mesh-study.json").write_text(json.dumps(data, indent=2))
    for theme, suffix in (("light", ""), ("dark", "-dark")):
        print("wrote", figure(data, theme, suffix))


if __name__ == "__main__":
    main()
