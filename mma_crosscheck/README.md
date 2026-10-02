# Neko-TOP MMA vs topopt_in_petsc MMA.cc — numerical cross-check

Runs the **unmodified** Neko-TOP CPU MMA (`sources/mma/mma.f90` +
`sources/mma/bcknd/cpu/mma_cpu.f90`, `dip` subsolver) and the **unmodified**
`MMA.cc` from topopt_in_petsc on the same problem, and compares every iterate.

Checked against neko-top `4a3d4f9` and topopt_in_petsc `73307a2`.

## Problem and settings

The 1D cantilever beam from `tests/regression/mma`: beam weight objective, tip
deflection plus 10 stress constraints, n = 221184, m = 11, x0 = 0.5. The
parameters are those of `1d_beam_dip.case`: asymptotes 0.5/1.2/0.7, c = 1000,
a = 0, no move limit. These equal the defaults of `MMA.cc`.

`neko_stubs.f90` replaces the few Neko modules that MMA uses (vector, matrix,
logger, json, ...). `beam_driver.f90` calls `mma_t` exactly the way
`mma_optimizer.f90` does.

**Harness validation:** the harness reproduces Neko-TOP's own stored
`reference_data/optimization_data_1d_beam_dip.csv` to ≤ 1.2e-10 (absolute) in
every column over all 20 iterations.

## Variants (`aage_variants.patch`, applied to `mma_cpu.f90` only)

- `-DAAGE_REG`: p0/q0/pij/qij exactly as `MMA::GenSub` with
  `constraintModification = false`. That means 0.5e-6/(U-L) on the objective
  and no regularisation on the constraints.
- `-DAAGE_TERM`: subsolver stopping rules exactly as `MMA::SolveDIP`:
  - inner loop runs while err > 0.9·ε;
  - err is carried over between ε levels;
  - outer loop stops at epsimin.

## Results (20 iterations)

| Comparison with MMA.cc | max rel. Δf0 | max abs. Δg | max abs. Δx |
|---|---|---|---|
| upstream Neko-TOP | 1.1e-1 | 4.1e-1 | 3.1e-1 |
| + Aage p/q coefficients | 5.5e-8 | 1.7e-7 | 4.8e-3 |
| + Aage constraint p/q only | 4.4e-2 | 2.0e-1 | 1.8e-1 |
| + Aage stopping rules only | 1.1e-1 | 4.1e-1 | 3.1e-1 |
| + coefficients and stopping rules | 1.2e-14 | 2.6e-13 | 5.1e-12 |
| same, move limit 0.2 (vs `SetOuterMovelimit`) | 1.8e-14 | 2.5e-13 | 5.1e-12 |

Full per-iteration tables are in `results_20it.txt`. That file also compares
the `MMA.cc` KKT residual with Neko-TOP's dip "KKT" measure.

## Run

```sh
export PKG_CONFIG_PATH=$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig   # PETSc.pc
./run.sh /path/to/neko-top /path/to/topopt_in_petsc 20          # optional 4th arg: move limit
```

Each run takes about 2 min on one core.

## Not covered

- the device backends;
- `pdip`;
- a ≠ 0, so the z-term difference is not exercised;
- constraint scaling in `mma_optimizer.f90`.
