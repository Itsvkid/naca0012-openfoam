# 05 — OpenFOAM 2D airfoil simulation

Steady-state incompressible flow over a 2D airfoil using `simpleFoam`.

**Status:** Not started
**Tool:** OpenFOAM — <https://www.openfoam.com> (free). On macOS run it in
Docker; native builds are more trouble than they are worth.

## Objective

The benchmark comparison is the deliverable, not the contour plot. Computed Cl
and Cd curves plotted against published wind tunnel data — Abbott & Von
Doenhoff for NACA sections, or NASA Langley's NACA 0012 validation case.

Anyone can produce a colourful pressure field. Showing yours matches measured
data — or explaining exactly where and why it does not — is what makes this
worth a reviewer's time.

## Core concepts

- OpenFOAM case hierarchy and dictionary configuration
- Structured C-grid generation with `blockMesh`
- Wall-bounded turbulence, y+ ≈ 1, k-omega SST
- Convergence assessment and force coefficient extraction

## Case structure

```
case/
  0/          U, p, k, omega, nut          initial and boundary fields
  constant/   turbulenceProperties, polyMesh
  system/     blockMeshDict, fvSchemes, fvSolution, controlDict
```

## Workflow

1. **Mesh** — C-grid around NACA 0012 via `blockMeshDict`. Grade toward the
   wall so first-cell height gives y+ ≤ 1 at your Reynolds number. Compute that
   height beforehand; do not guess it.
2. **Fields** — set `U`, `p`, `k`, `omega`, `nut` in `0/`. Inlet turbulence
   quantities should reflect the wind tunnel you are comparing against.
3. **Solver** — bounded second-order upwind for divergence in `fvSchemes`;
   convergence criteria in `fvSolution`; add the `forces` function object to
   `controlDict` so Cl and Cd are logged every iteration.
4. **Run:**
   ```bash
   blockMesh && checkMesh && simpleFoam > log.simpleFoam 2>&1 &
   ```
   `checkMesh` must pass before `simpleFoam`. A mesh with non-orthogonality
   above ~70 will produce numbers, and they will be wrong.
5. **Post-process** — `paraFoam`. Pressure contours, boundary layer profiles,
   trailing-edge flow.

## Validation

- [ ] y+ verified from the solution, not just estimated beforehand
- [ ] Residuals converged, and force coefficients flat — residuals alone are
      not convergence
- [ ] Mesh independence across three refinements
- [ ] Cl vs alpha plotted against published data
- [ ] Cd vs alpha plotted against published data, with the discrepancy
      discussed honestly

Steady RANS will not match experiment past stall. Say that in the write-up
instead of stopping the sweep where the agreement stops looking good.

## Deliverables

- [ ] Complete, runnable case in `case/`
- [ ] Validation plots in `validation/`
- [ ] Figures for the site and LinkedIn
- [ ] Its own GitHub repository

## Repository hygiene

Commit `0/`, `constant/` (minus `polyMesh`) and `system/`. Never commit time
directories, `polyMesh/` or VTK output — the root `.gitignore` excludes them.
The case dictionaries are the work; everything else regenerates.

## Log

| Date | What was done |
|---|---|
| | |
