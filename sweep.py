"""Angle-of-attack sweep on a fixed mesh.

Incidence is applied by rotating the freestream vector, not by re-meshing, so
every point in the polar sits on the identical grid. Re-meshing per alpha would
fold mesh differences into the polar and there would be no way afterwards to
tell those apart from aerodynamics.

The C-grid's wake line is aligned with x, so at high incidence the wake leaves
at an angle to it and is resolved progressively less well. That is a real limit
of this mesh and it is reported rather than hidden.

Saves after every angle and skips angles already recorded, so an interruption
costs one run.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from src.case import Case          # noqa: E402
from src.cgrid import CGrid        # noqa: E402

HERE = Path(__file__).parent
RESULTS = HERE / "validation" / "alpha-sweep.json"
REYNOLDS = 6.0e6
ITERATIONS = 2000

# L3 from the mesh study: 37 800 cells, y+ ~ 0.95, GCI 6.4% on Cd.
MESH = dict(n_surface=120, n_normal=90, n_wake=90)
Y_PLUS_TARGET = 1.0

ALPHAS = [-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0]


def foam(command: str, check: bool = False) -> str:
    """Run a command in the container and return stdout only.

    stdout and stderr are deliberately kept apart. Merging them meant a shell
    error — `grep: ...: No such file or directory` — was parsed as a data row,
    and the resulting float() failure pointed at the parser rather than at the
    solve that had actually failed.
    """
    r = subprocess.run([str(HERE / "foam.sh"), command],
                       capture_output=True, text=True)
    if check and r.returncode != 0:
        raise RuntimeError(f"{command!r} failed:\n{r.stderr.strip()[-600:]}")
    return r.stdout


def build_mesh() -> int:
    base = CGrid(**MESH)
    grid = CGrid(**MESH,
                 expansion_normal=base.expansion_for_y_plus(Y_PLUS_TARGET, REYNOLDS))
    grid.write(HERE / "case" / "system" / "blockMeshDict")
    foam("rm -rf postProcessing [1-9]* log.* && blockMesh > /dev/null 2>&1")
    return 2 * MESH["n_surface"] * MESH["n_normal"] + 2 * MESH["n_normal"] * MESH["n_wake"]


def run_alpha(alpha: float) -> dict:
    Case(alpha_deg=alpha, reynolds=REYNOLDS, iterations=ITERATIONS).write(HERE / "case")
    started = time.time()
    # The mesh is untouched: only the 0/ fields and controlDict change, so
    # blockMesh is not re-run.
    foam("rm -rf postProcessing [1-9]* log.simpleFoam; "
         "simpleFoam > log.simpleFoam 2>&1", check=False)
    elapsed = time.time() - started

    # The solver reports failure inside its log, not through an exit code that
    # survives the redirect, so the log is the thing to check.
    if "FOAM FATAL" in foam("cat log.simpleFoam"):
        raise RuntimeError(
            f"alpha {alpha}: simpleFoam failed\n"
            + foam("grep -A6 'FOAM FATAL' log.simpleFoam | head -20"))

    log = foam("tail -300 log.simpleFoam")
    yplus = re.findall(r"y\+ : min = ([\d.e+-]+), max = ([\d.e+-]+), "
                       r"average = ([\d.e+-]+)", log)
    rows = [r.split() for r in
            foam("grep -v '^#' postProcessing/forceCoeffs/0/coefficient.dat")
            .strip().splitlines() if r.strip()]
    if not rows:
        raise RuntimeError(
            f"alpha {alpha}: no coefficients written\n"
            + foam("tail -15 log.simpleFoam"))

    cd_tail = [float(r[1]) for r in rows[-400:]]
    cl_tail = [float(r[4]) for r in rows[-400:]]
    return {
        "alpha": alpha,
        "cl": float(rows[-1][4]),
        "cd": float(rows[-1][1]),
        "cm": float(rows[-1][7]),
        "y_plus_max": float(yplus[-1][1]) if yplus else None,
        "cl_drift": (max(cl_tail) - min(cl_tail)),
        "cd_drift": (max(cd_tail) - min(cd_tail)),
        "iterations": int(rows[-1][0]),
        "seconds": round(elapsed, 1),
    }


def main() -> None:
    data = (json.loads(RESULTS.read_text()) if RESULTS.exists()
            else {"reynolds": REYNOLDS, "mesh": MESH, "points": []})
    done = {p["alpha"] for p in data["points"]}

    cells = build_mesh()
    print(f"mesh: {cells} cells (rebuilt once; incidence rotates the freestream)",
          flush=True)

    for alpha in ALPHAS:
        if alpha in done:
            print(f"  a={alpha:+5.1f}  already recorded", flush=True)
            continue
        p = run_alpha(alpha)
        data["points"].append(p)
        data["points"].sort(key=lambda q: q["alpha"])
        RESULTS.parent.mkdir(exist_ok=True)
        RESULTS.write_text(json.dumps(data, indent=2))
        print(f"  a={alpha:+5.1f}  Cl {p['cl']:+.4f}  Cd {p['cd']:.5f}  "
              f"y+max {p['y_plus_max']:.2f}  drift Cl {p['cl_drift']:.1e}  "
              f"{p['seconds']:.0f}s", flush=True)

    print(f"{len(data['points'])} angles recorded")


if __name__ == "__main__":
    main()
