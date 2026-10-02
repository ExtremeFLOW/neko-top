The CPU `dip` path now reproduces `MMA.cc`/`MMA.h` to round-off. I tested a = 0, a = 1, move limit 0.2, and 1, 2 and 3 MPI ranks; f0 agrees to ≤ 4e-14 and the KKT residual to ≤ 5e-8 relative. `pdip` is bit-identical to before. Nothing is committed.

**What I changed, on top of your `2a04bf5f progress claude` commit:**
- **`mma_cpu.f90`:**
  - Aage's p/q coefficients are used only when `subsolver = "dip"`. `pdip` goes back to Svanberg's original form, as you chose.
  - The `dip` KKT check is now a port of `MMA::KKTresidual`, using the true bounds. Before, it never tested stationarity.
  - The z term now matches `MMA.cc`: 0.05z², z = max(0, 10(λᵀa − a0)), and −10aaᵀ in the Hessian, using `a0` where `MMA.cc` hard-codes 1.
  - If the barrier loop never runs, the subsolver now returns the input design, as `MMA.cc` does. It used to return an uninitialised x.
  - Removed the residual computation the staged stopping-rule change had made dead, plus a few comment fixes.
- **`mma.f90`:** `z`/`zeta` are now initialised to 0; the first KKT call used to read them uninitialised. Init now rejects epsimin ≤ 0, since the staged change dropped the old 1e-12 floor that guarded against it.
- **`mma_optimizer.f90`:** the constraint scaling used in the update is now also applied before the KKT check. It was missing whenever `scale` ≠ 1 or `auto_scale` was on.
- **Docs:**
  - `MMA.md` updated: p/q form per subsolver, the KKT definition and its floor of about m × the last barrier level, and the z weight.
  - New `MMA.md` section on how Neko-TOP relates to topopt_in_petsc, covering the asymptote defaults you kept and the objective normalisation note.
  - `c` default corrected to 1000 in `configuration.md`.
- **Regression data:** I regenerated the `dip` reference CSV in CI's configuration (2 ranks, `OMP_NUM_THREADS=1`, `NEKO_GS_COMM=MPI`), and `check.py` passes. The real optimizer pipeline reproduces `MMA.cc` end to end.

**Rebasing onto develop:** the `dip` CSV will conflict with #531's. Keep this branch's version, which was generated with the same configuration. For `pdip`, take develop's; the new code reproduces it exactly.

**Deliberately left different from `MMA.cc`, all documented in `MMA.md`:**
- the asymptote defaults (0.2/1.05/0.65 vs 0.5/1.2/0.7), as you chose;
- the starting multipliers, max(1, c/2) instead of c/2, which only differ for c < 2;
- a 1e-5 floor on the bound span for the initial asymptotes, so variables with xmin = xmax don't produce NaN;
- the active-set test at machine precision, and pivoted LAPACK instead of the reference's unpivoted LU. Both only change results at round-off.

**Not done:**
- **Device backend, out of scope and untested:** it still has the old p/q form, stopping rules, z term and KKT, so it no longer matches the new CPU reference. It also has two bugs:
  - z is computed with `device_glsum`, a global sum over values every rank already holds, so with a ≠ 0 it comes out P times too large on P ranks;
  - `zerom` is never zeroed.
- **Pre-existing problems I found but didn't fix:**
  - a rank with zero design variables reads unallocated arrays and will likely crash;
  - `auto_scale` divides by the first constraint value, which blows up as that value approaches 0;
  - `n_global` is a 32-bit integer;
  - the DGESV interface breaks single-precision builds;
  - `check.py` stops at the first failing file, so `pdip` isn't checked when `dip` fails.

**Side effects:**
- `tests/regression/mma/reg_mma_bin` was rebuilt in place from a separate Release build in the scratchpad; your `build/` is untouched.
- Your `run.sh` edit is unchanged.
- In `mma_crosscheck/` I only made `run.sh` work against the current code without PETSc installed, and deleted its 237 MB of output.