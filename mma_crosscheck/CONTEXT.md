# Context: Neko-TOP MMA vs. topopt_in_petsc `MMA.cc` (review handoff)

This file hands off a review done in a claude.ai session, so the work can
continue in Claude Code without being redone. The findings are verified at the
commits listed below. Re-check line numbers if the repositories have moved on.

## Task

Check whether Neko-TOP's MMA implementation follows the reference `MMA.cc` in
https://github.com/topopt/topopt_in_petsc. That reference is Svanberg's MMA
with the dual interior-point subsolver from Aage & Lazarov (2013).

## Status

**Done:**
- line-by-line static review of the CPU backend, the device backend and the
  optimizer driver;
- a numerical cross-check of the CPU `dip` path against the unmodified
  reference.

**Not done:**
- device numerics (no GPU was available);
- `pdip` (it has no counterpart in the reference);
- runs with a ≠ 0;
- a runtime test of the constraint-scaling bug;
- any fixes in neko-top.

## Versions checked

- neko-top `4a3d4f9` (2026-09-29)
- topopt_in_petsc `73307a2` (2025-06-12)
- Neko `6a661f1` (2026-10-01), used only for the semantics of
  `src/registries/scratch_registry.f90`

## Bottom line

On the beam regression problem, run with the reference's parameters, the CPU
`dip` path becomes **numerically identical** to `MMA.cc` once two things are
switched to Aage's versions:

1. the p/q coefficients in `mma_gensub_cpu`;
2. the stopping rules in `mma_subsolve_dip_cpu`.

With both switched, the two codes agree to ≤ 2e-14 relative in f0 and ≤ 5.1e-12
in x over 20 iterations. The p/q coefficients are the material difference:
upstream, f0 differs by 11% after iteration 1.

Separately, the `dip` "KKT residual" that the optimizer uses as its
convergence test is not a KKT residual.

## Code map

### Reference (topopt_in_petsc)

- **`MMA.cc` 195–291:** constructors. asyminit 0.5, asymdec 0.7, asyminc 1.2;
  class defaults a = 0, c = 1000, d = 0.
- **`MMA.cc` 386–405:** `SetOuterMovelimit`. The driver calls it before every
  `Update` (`main.cc` 81), with movlim 0.2 (`TopOpt.cc` 128).
- **`MMA.cc` 428–496:** `KKTresidual`.
- **`MMA.cc` 522–649:** `GenSub`.
  - asymptotes: 532–591
  - α/β: 595–596
  - p0/q0: 599–600
  - pij/qij: 601–613
  - b: 615–632
- **`MMA.cc` 651–981:** the dual solver.
  - `SolveDIP`: 651–688
  - `XYZofLAMBDA`: 690–740
  - `DualGrad`: 742–777
  - `DualHess`: 779–880
  - `DualLineSearch`: 882–900
  - `DualResidual`: 902–946
  - unpivoted LU: 948–981
- **`main.cc` 45:** `// mma->SetAsymptotes(0.2, 0.65, 1.05);`, commented out.
- **`TopOpt.cc` 391–397:** a = 0, c = 1000, d = 0.

### Neko-TOP

- **`sources/mma/mma.f90`:** defaults at 147–160; epsimin default
  1e-9·√(m + n_global) at 281 and 456.
- **`sources/mma/bcknd/cpu/mma_cpu.f90`:**
  - `mma_dip_KKT_cpu`: 193–231
  - `mma_gensub_cpu`: 237–374
  - `mma_subsolve_pdip_cpu`: 378–770
  - `mma_subsolve_dip_cpu`: 774–1089
- **`sources/mma/bcknd/device/mma_device.f90`:**
  - `mma_dip_KKT_device`: 112–137
  - gensub driver: 225–299
  - `mma_subsolve_dip_device`: 770–1056
- **Device kernels:** `device/cuda/mma_kernel.h` and `device/hip/mma_kernel.h`.
  The p/q coefficients are in `mma_sub3_kernel`, at line 375 in both files.
- **`sources/optimizer/mma_optimizer.f90`:** tolerance default 1e-3 at
  134–135, scaling at 328–336, update at 340, KKT at 369, convergence test at
  372.
- **`documentation/pages/theory/MMA.md`:** parameter table. Line 211 says
  c = 1000.
- **`tests/regression/mma/`:** the 1D beam case. `check.py` compares only
  against its own stored CSVs, not against the reference.

## What matches `MMA.cc` (dip path)

The numerical cross-check confirms all of these.

- **Asymptote initialisation for `iter < 3`:** `mma_cpu.f90` 268–271 vs
  `MMA.cc` 532–537. The optimizer's iteration counter starts at 1, so this
  equals the reference's internal k.
- **Asymptote update and clamps** to [0.01, 10]·span (`RobustAsymptotesType = 0`):
  273–293 vs 555–573. The clamp order for `upp` differs but gives the same
  result.
- **Internal move limit** (254–261), equivalent to `SetOuterMovelimit`.
  - The move-limited span feeds the asymptotes and α/β exactly as in the
    reference.
- **α = max(xmin, 0.9L + 0.1x), β = min(xmax, 0.9U + 0.1x):** 311–312.
- **bᵢ:** 356–371.
- **x(λ) closed form with clamping, and y = max(0, λ − c):** 883–901.
  - Both codes treat d as 1 and ignore the user's d.
- **Dual gradient and residual:** 904–937.
- **Dual Hessian:** 943–1002.
  - PQᵀ·diag(df2)·PQ with active variables removed;
  - −1 where λ > c (Neko-TOP tests y > 0);
  - −μ/λ on the diagonal;
  - diagonal shift max(−1e-4·trace/m, 1e-7).
- **Newton step, μ update and line search** with the 1.005/1.01 rule:
  1004–1029. DGESV replaces the reference's unpivoted LU.
- **Solver schedule:** ε reduced by ×0.1, at most 100 inner iterations,
  tolerance 1e-9·√(m + n).
- **xold update order:** 1079–1081 vs 512–513.

## Deviations

### D1 — p/q coefficients (material)

**Neko-TOP** (`mma_cpu.f90` 323–349, and the same in `mma_sub3_kernel`):

```
p0 = (U-x)^2 * (1.001*max(df,0) + 0.001*max(-df,0) + 1e-5/max(xmax-xmin, 1e-5))
q0 = (x-L)^2 * (0.001*max(df,0) + 1.001*max(-df,0) + 1e-5/max(xmax-xmin, 1e-5))
```

pij and qij use the same form for every constraint. This is Svanberg's style,
with a 10⁻⁵/(xmax − xmin) term as in his 2007 MMA/GCMMA notes.

**Reference** (`MMA.cc` 593–613, default `constraintModification = false`):

```
p0  = (U-x)^2 * (max(df,0) + 0.001*|df| + 0.5e-6/(U-L))
q0  = (x-L)^2 * (max(-df,0) + 0.001*|df| + 0.5e-6/(U-L))
pij = (U-x)^2 * max(dg,0);  qij = (x-L)^2 * max(-dg,0)
```

The constraints get no regularisation in the reference.

The gradient terms are identical, since max(df,0) + 0.001·|df| equals
1.001·max(df,0) + 0.001·max(−df,0). The differences are:
- the absolute 1e-5 term;
- that term also being added to every constraint.

**Why it matters on the beam problem:**
- The objective gradient is 6.35e-5 per variable, so the 1e-5 term is about 16%
  of it.
- The term is added for all 221,184 variables in each of the 11 constraints.
  Each of the 10 stress constraints actually depends on a single element.
- In Neko-TOP, xmax − xmin is the move-limited span, so a move limit makes the
  term larger still (1e-5/0.4 with a 0.2 move limit).

### D2 — z term (only matters if a ≠ 0)

**Neko-TOP:** z = max(0, λᵀa − a0) (888, 1043), which corresponds to a +½z²
term. The Hessian gets −aaᵀ (975–981).

**Reference:** z = max(0, 10·(λᵀa − 1)) (`MMA.cc` 713), which corresponds to
+0.05z², with a0 hard-coded to 1. The Hessian gets −10·aaᵀ (849–855).

Both codes enable the Hessian term when λᵀa > 0, although z > 0 only when
λᵀa > a0. This is a quirk of the reference that Neko-TOP copied. Neko-TOP is
otherwise internally consistent. With a = 0 (the default in both) nothing
differs.

### D3 — Defaults

**Asymptote settings:**
- Neko-TOP (`mma.f90` 154–156): asyinit 0.2, asyincr 1.05, asydecr 0.65.
- Reference: 0.5, 1.2, 0.7.
- The 0.2/0.65/1.05 set appears only in the commented-out line `main.cc` 45.

**c:**
- Neko-TOP: `c_default = 100` (`mma.f90` 149).
- Reference: 1000.
- Neko-TOP's own docs (`MMA.md` 211) also say 1000, so code and docs disagree.

**Initial λ:** Neko-TOP uses max(1, c/2) (`mma_cpu.f90` 848); the reference uses
c/2 (`MMA.cc` 655). These differ only when c < 2.

**Regression case:** `tests/regression/mma/cases/1d_beam_dip.case` uses the
reference's values: 0.5/1.2/0.7, c = 1000, move_limit −1 (off), a = 0, d = 1.

### D4 — Subsolver stopping rules (immaterial on the test problem)

**Neko-TOP:**
- The outer loop runs while ε > max(0.9·epsimin, 1e-12) (856).
- The residual is recomputed with the new ε at each outer level (903–917).
- The inner loop exits when the residual is below ε (924).

**Reference** (`MMA.cc` 658–686):
- The outer loop runs while ε > epsimin.
- The inner loop runs while err > 0.9·ε.
- err is carried over stale from the previous ε level.

Neko-TOP's recomputation is arguably more correct.

### D5 — `dip` KKT residual is not a KKT residual (important)

`mma_dip_KKT_cpu` (193–231) and `mma_dip_KKT_device` (112–137) compute only:
- relambda = f(x_new) − a·z − y + μ
- remu = λ·μ

The arguments x, df0dx and dfdx are never used.

The subsolve ends with h̃(x_new) − a·z − y + μ ≈ 0, where h̃ is the MMA
approximation. So:
- relambda ≈ f(x_new) − h̃(x_new), which is approximation error;
- remu ≈ the final barrier ε.

There is no stationarity test.

The reference `KKTresidual` (428–496) does include stationarity:
- it computes ri = df0 + Σλ·dg;
- it estimates bound multipliers at x ≤ xmin + 1e-5 and x ≥ xmax − 1e-5;
- it adds bound complementarity and λᵀ(a·z + y − f).

`mma_optimizer.f90` 372 declares convergence when `residumax < tolerance`
(default 1e-3), and `dip` is the default subsolver.

**Evidence** (stored regression CSV, reproduced by the harness):
- KKTmax is pinned at exactly 1.000e-6 (the final ε) from iteration 11 onward.
- So the case's tolerance of 1e-6 is never met.
- With 1e-3 the run would stop at iteration 8 (KKTmax 2.0e-4), while the
  maximum design change is still 3.7e-2.

### D6 — Constraint scaling is not applied before the KKT check

`mma_optimizer.f90` 328–336 scales the constraint values and sensitivities
before `update` (340). After re-evaluation, the `KKT` call (369) receives
unscaled values, while λ, μ, y and z come from the scaled subproblem.

This is wrong whenever `mma.scale` ≠ 1 or `auto_scale` is true. It was found by
reading the code only and has not been exercised.

### D7 — Device: `zerom` is never zeroed (latent)

`mma_device.f90` 798 requests `zerom` with `clear = .false.` and never fills
it. It is then used as the zero in y = max(λ − c, zerom), at 852 and 993.

Neko's `request_vector` reuses an existing slot of the same size without
zeroing it. Today `zerom` happens to be zero, because:
- the ninth m-sized slot is never written elsewhere;
- m ≠ n_local.

This is fragile; it would break if, for example, n_local == m.

Fix: `device_cfill(zerom%x_d, 0.0_rp, this%m)`, or request with
`clear = .true.`, or use a scalar max.

### Device backend in general (read, not run)

These mirror the CPU code:
- the gensub kernels (`mma_sub2/3/4_kernel`);
- x(λ), Ljjxinv and the Hessian updates;
- `maxval`, which uses |·|;
- `maxval2`, used for the step length.

The linear solve uses cuSOLVER/hipSOLVER.

### Minor

- **x_diff floor at initialisation:** the 1e-5 floor on x_diff is also applied
  to the initial asymptotes; the reference applies it only from iteration 3.
  Negligible.
- **Active-set test:** Neko-TOP uses `x − α < NEKO_EPS` on CPU and
  `|x − α| ≤ 1e-16` on device; the reference tests the unclamped x(λ) < α.
  Equivalent in practice.
- **Misleading comment:** "this n is global" at `mma_cpu.f90` 604 and 962. It
  is the local n; the following Allreduce makes the result correct.
- **λ ≥ 0 clamp:** the reference clamps λ explicitly. Neko-TOP relies on the
  line search, which keeps λ > 0.

## Numerical cross-check

### Package

`mma_crosscheck.zip`, produced in the same session. It contains:
- `README.md`, `run.sh`;
- `neko_stubs.f90`, `mma_device_stub.f90`, `beam_driver.f90`;
- `beam_ref.cc`;
- `aage_variants.patch`;
- `compare.py`;
- `results_20it.txt`.

### Problem

The 1D cantilever beam from `tests/regression/mma`:
- weight objective, with tip-deflection and 10 stress constraints;
- n = 221,184 (8192 elements × 27 points from `box.nmsh`), m = 11;
- x0 = 0.5, bounds [0, 1];
- 20 iterations.

The parameters are those of `1d_beam_dip.case`, which are also `MMA.cc`'s
defaults:
- asymptotes 0.5/1.2/0.7;
- a0 = 1, a = 0, c = 1000, d = 1;
- max_iter 100;
- epsimin default (4.70e-7);
- no move limit.

### How each side is built

**Neko-TOP:**
- The unmodified `mma.f90` and `mma_cpu.f90` are compiled against stub Neko
  modules (vector, matrix, logger, json, comm, …).
- `beam_driver.f90` makes the same calls as `mma_optimizer.f90`:
  `init_from_components`, `update`, then `kkt`.

**Reference:**
- The unmodified `MMA.cc`/`MMA.h` are built with PETSc 3.19 and `beam_ref.cc`.
- `KKTresidual` is called with the true bounds.
- `SetOuterMovelimit` is used when a move limit is given.

**Validation:** the harness reproduces
`reference_data/optimization_data_1d_beam_dip.csv` to ≤ 1.2e-10 absolute in
every column over 20 iterations.

**Variants**, via `aage_variants.patch` applied to `mma_cpu.f90`:
- `-DAAGE_REG`: Aage's p/q coefficients.
- `-DAAGE_TERM`: Aage's stopping rules.

### Results (20 iterations, compared with `MMA.cc`)

| Neko-TOP variant | max rel Δf0 | max abs Δg | max abs Δx | rel Δf0 at it 20 |
|---|---|---|---|---|
| upstream | 1.1e-1 | 4.1e-1 | 3.1e-1 | 2.9e-8 |
| Aage p/q coefficients only | 5.5e-8 | 1.7e-7 | 4.8e-3 | 4.3e-8 |
| Aage constraint p/q only | 4.4e-2 | 2.0e-1 | 1.8e-1 | 2.3e-7 |
| Aage stopping rules only | 1.1e-1 | 4.1e-1 | 3.1e-1 | 4.8e-9 |
| coefficients + stopping rules | 1.2e-14 | 2.6e-13 | 5.1e-12 | 1.0e-14 |
| same, move limit 0.2 (vs `SetOuterMovelimit`) | 1.8e-14 | 2.5e-13 | 5.1e-12 | 4.8e-15 |

- **Upstream:** f0 after iteration 1 is 8.869 vs 7.979. Both runs converge to
  f0 ≈ 8.1914263.
- **Convergence measures:** each value below is taken along its own code's
  trajectory.

  | it | `MMA.cc` KKTresidual (inf) | Neko-TOP dip KKTmax | Neko-TOP max Δx |
  |---|---|---|---|
  | 1 | 3.0e-1 | 1.9e-1 | 3.9e-1 |
  | 7 | 1.2e-4 | 2.1e-3 | 6.4e-2 |
  | 8 | 1.8e-5 | 2.0e-4 | 3.7e-2 |
  | 11 | 1.3e-5 | 1.000e-6 | 4.8e-3 |
  | 20 | 1.3e-5 | 1.000e-6 | 2.1e-6 |

  The reference value plateaus at about 1.3e-5; this was not investigated.

### Environment used

- Ubuntu 24.04 container, 1 core.
- Packages:
  `apt-get install gfortran libopenmpi-dev openmpi-bin liblapack-dev libblas-dev petsc-dev`
  (PETSc 3.19; the install takes several minutes).
- `export PKG_CONFIG_PATH=/usr/lib/petscdir/petsc3.19/x86_64-linux-gnu-real/lib/pkgconfig`
- Running as root needs `OMPI_ALLOW_RUN_AS_ROOT=1` and
  `OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1`; `run.sh` sets both.
- Each 20-iteration run takes about 2 minutes. Builds take seconds.

### Rerun

```sh
unzip mma_crosscheck.zip && cd mma_crosscheck
./run.sh /path/to/neko-top /path/to/topopt_in_petsc 20        # optional 4th arg: move limit
```

`run.sh` builds and runs three cases: upstream Neko-TOP, the `aage` variant
(`-DAAGE_REG -DAAGE_TERM`) and the reference. It then prints both comparisons.
A 3-iteration smoke test of the script, from clean checkouts, reproduced the
numbers above.

## Suggested next steps (not started)

1. **p/q coefficients:** adopt Aage's, or make the choice a JSON option.
   - Change `mma_gensub_cpu` and `mma_sub3_kernel` in both the CUDA and HIP
     kernels.
   - Rerun the harness; the matching variant should stay at round-off
     agreement with `MMA.cc`.
2. **`dip` KKT check:** replace `mma_dip_KKT_cpu`/`_device` with a port of
   `MMA::KKTresidual`.
   - Use the true box bounds `this%xmin`/`this%xmax`, not the move-limited
     ones.
   - Alternatively, document what the current measure is and revisit the
     default tolerance.
3. **Scaling:** in `mma_optimizer_step`, apply the constraint scaling before
   the `KKT` call as well.
4. **Device `zerom`:** fix it in `mma_subsolve_dip_device`.
5. **Defaults and docs:** align the code and docs on c (100 vs 1000), and state
   that the asymptote defaults differ from `MMA.cc`.
6. **Optional, z term** (needed only if a ≠ 0 should behave like the
   reference): apply the factor 10, use a0, and switch the Hessian term on when
   λᵀa > a0.
7. **Optional, CI:** turn the harness into a cross-check test. The current
   regression test only compares against its own stored output.
8. **Still unverified:**
   - device numerics vs CPU on a GPU;
   - `pdip` vs Svanberg's `subsolv`;
   - behaviour with a ≠ 0.
