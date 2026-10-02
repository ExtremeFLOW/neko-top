# Neko-TOP MMA vs topopt_in_petsc MMA.cc: numerical cross-check

This harness runs the Neko-TOP CPU MMA (`sources/mma/mma.f90` and
`sources/mma/bcknd/cpu/mma_cpu.f90`, `dip` subsolver) and the **unmodified**
`MMA.cc` from topopt_in_petsc on the same problem. It then compares every
iterate and the KKT residual.

The comparison was checked against topopt_in_petsc `73307a2`.

## Problem and settings

The problem is the 1D cantilever beam from `tests/regression/mma`:
- objective: beam weight;
- constraints: tip deflection and 10 stress constraints;
- n = 221184, m = 11, x0 = 0.5.

The parameters are those of `1d_beam_dip.case`, which equal the defaults of
`MMA.cc`:
- asymptotes 0.5/1.2/0.7;
- c = 1000, a = 0;
- no move limit.

The harness also runs three variants:
- a move limit of 0.2, applied through `SetOuterMovelimit` in the reference;
- a = 1, which exercises the z term;
- 3 MPI ranks on the Neko-TOP side.

What each file does:
- `neko_stubs.f90` replaces the few Neko modules that MMA uses (vector,
  matrix, logger, json, ...).
- `beam_driver.f90` block-partitions the design over the MPI ranks. It calls
  `mma_t` the way `mma_optimizer.f90` does.
- `beam_ref.cc` drives `MMA.cc` and calls `KKTresidual` with the true bounds.
- `petsc_shim/petsc.h` is a minimal stand-in for the PETSc Vec API. It is
  used when PETSc is not visible to pkg-config.

## Run

```sh
./run.sh /path/to/neko-top /path/to/topopt_in_petsc 20
```

Each run takes about 1–2 minutes on one core. `compare.py` prints one
summary line per comparison. Drop `-q` from it to get per-iteration tables.

## Results (20 iterations, working tree after the alignment)

| Comparison with MMA.cc | max rel. Δf0 | max abs. Δg | max abs. Δx | max rel. ΔKKTmax |
|---|---|---|---|---|
| a = 0 | 1.2e-14 | 2.6e-13 | 5.1e-12 | 2.8e-9 |
| a = 0, 3 ranks | 4.1e-14 | 2.6e-13 | 5.1e-12 | 9.9e-9 |
| move limit 0.2 | 1.8e-14 | 2.5e-13 | 5.1e-12 | 2.8e-9 |
| a = 1 | 1.2e-14 | 1.8e-12 | 4.0e-12 | 3.6e-9 |

The KKT residual converges to about m times the last barrier level
(1.3e-5 here). Near that level it comes from cancellation between O(1) terms,
which is why its relative agreement is about 1e-9 rather than 1e-14.

The `pdip` path keeps Svanberg's approximation, and its output is bit-identical
to that before the alignment.

`CONTEXT.md` and `results_20it.txt` are the handoff from the earlier review.
They describe the state before the alignment.

## Not covered

- the device backends;
- the objective normalisation `fscale = 10/f0(x0)` done by topopt_in_petsc's
  `main.cc`, which is not part of `MMA.cc` (both sides here use the unscaled
  objective).
