"""Mesh independence study.

Refines every direction by a constant ratio, including the near-wall spacing,
so the sequence is *systematic* and Richardson extrapolation applies. Holding
the wall spacing fixed while refining only tangentially is the common shortcut
and it invalidates the extrapolation, because the meshes then differ in more
than one parameter.

Usage:  /opt/anaconda3/envs/pyocc_env/bin/python study.py
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from src.case import Case          # noqa: E402
from src.cgrid import CGrid        # noqa: E402

HERE = Path(__file__).parent
REYNOLDS = 6.0e6
RATIO = 1.5          # linear refinement ratio between consecutive levels
ITERATIONS = 3000

# (label, n_surface, n_normal, n_wake, y+ target). Counts and wall spacing both
# scale by RATIO, so h ~ 1/n throughout.
LEVELS = [
    ("L1", 54, 40, 40, 2.25),
    ("L2", 80, 60, 60, 1.50),
    ("L3", 120, 90, 90, 1.00),
    ("L4", 180, 135, 135, 0.667),
]


def foam(command: str) -> str:
    result = subprocess.run([str(HERE / "foam.sh"), command],
                            capture_output=True, text=True)
    return result.stdout + result.stderr


def run_level(label, ns, nn, nw, y_target) -> dict:
    base = CGrid(n_surface=ns, n_normal=nn, n_wake=nw)
    grid = CGrid(n_surface=ns, n_normal=nn, n_wake=nw,
                 expansion_normal=base.expansion_for_y_plus(y_target, REYNOLDS))
    grid.write(HERE / "case" / "system" / "blockMeshDict")
    Case(alpha_deg=0.0, reynolds=REYNOLDS, iterations=ITERATIONS).write(HERE / "case")

    print(f"  {label}: meshing…", flush=True)
    out = foam("rm -rf postProcessing [1-9]* log.* && blockMesh 2>&1 | tail -2")
    cells = 2 * ns * nn + 2 * nn * nw

    print(f"  {label}: solving {cells} cells…", flush=True)
    foam(f"simpleFoam > log.simpleFoam 2>&1")

    log = foam("cat log.simpleFoam | tail -400")
    yplus = re.findall(r"y\+ : min = ([\d.e+-]+), max = ([\d.e+-]+), "
                       r"average = ([\d.e+-]+)", log)
    coeffs = foam("cat postProcessing/forceCoeffs/0/coefficient.dat "
                  "| grep -v '^#'").strip().splitlines()
    rows = [r.split() for r in coeffs if r.strip()]
    last = rows[-1]
    tail = [float(r[1]) for r in rows[-500:]]
    spread = (max(tail) - min(tail)) / abs(sum(tail) / len(tail))

    return {
        "label": label, "n_surface": ns, "n_normal": nn, "n_wake": nw,
        "cells": cells,
        "first_cell": grid.first_cell_height(),
        "cd": float(last[1]), "cl": float(last[4]),
        "y_plus_max": float(yplus[-1][1]) if yplus else None,
        "y_plus_avg": float(yplus[-1][2]) if yplus else None,
        "cd_spread_last500": spread,
    }


RESULTS = HERE / "validation" / "mesh-study.json"


def load() -> dict:
    if RESULTS.exists():
        return json.loads(RESULTS.read_text())
    return {"reynolds": REYNOLDS, "ratio": RATIO, "levels": []}


def save(data: dict) -> None:
    RESULTS.parent.mkdir(exist_ok=True)
    RESULTS.write_text(json.dumps(data, indent=2))


def main() -> None:
    """Runs the levels not already recorded, saving after each one.

    Saving incrementally rather than at the end because the finest level takes
    far longer than the rest, and an interruption there previously discarded
    every completed result with it.
    """
    data = load()
    done = {l["label"] for l in data["levels"]}

    for level in LEVELS:
        label = level[0]
        if label in done:
            print(f"  {label}: already recorded, skipping", flush=True)
            continue
        data["levels"].append(run_level(*level))
        data["levels"].sort(key=lambda l: l["cells"])
        save(data)
        r = data["levels"][-1]
        print(f"  {label}: {r['cells']:>6} cells  Cd {r['cd']:.6f}  "
              f"Cl {r['cl']:+.2e}  y+max {r['y_plus_max']}  "
              f"drift {r['cd_spread_last500']:.1e}  (saved)", flush=True)

    print(f"{len(data['levels'])} levels recorded in {RESULTS.name}")


if __name__ == "__main__":
    main()
