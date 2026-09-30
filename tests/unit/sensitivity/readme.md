# Finite difference Sensitivity Test

This test verifies the sensitivity analysis of a design variable using finite
difference methods. It checks that the computed sensitivities match the
expected values within a specified tolerance.

The test is tagged as a `unit` test, which means it is mandatory for the
build to pass our CI/CD pipeline.

## Test Overview

The test is done through a few common files:

- `prepare.sh`: Script designed to construct a mesh for us to work on.
- `sensitivity.f90`: The shared finite-difference sensitivity checker module
  (MPI-safe reduction handling + the tolerance assertion). This same file is
  also used by `tests/regression/sensitivity/` — it is the single source of
  truth for both, see that directory's `CMakeLists.txt`.
- `problem_tester.f90`: The generic driver. It auto-detects whether the case
  under test defines an objective or a constraint and drives
  `compute_sensitivity` accordingly — it does not need to change when a new
  case is added.
- `*.case`: One Neko case file per test, each exercising a single
  objective/constraint type. The `fd_test_tolerance` key under
  `optimization` (optional, JSON) overrides the default assertion tolerance
  — linear/state-independent functionals (e.g. the volume constraint) can use
  a tight round-off tolerance; PDE-coupled objectives need a looser one
  matched to their discretisation and noise floor. Determine it empirically
  per case, not by guessing, from the FD check's error column `e`,
  `(FD - adjoint) / |adjoint|` at the probed dof:

  1. The model is `e = C + A m^p` with `m = |eps|` (the perturbation
     magnitude; the harness may use negative perturbations). For two
     points `a` (larger `m`) and `b`, Richardson extrapolation gives
     `A = (e_a - e_b) / (m_a^p - m_b^p)` and `C = e_b - A m_b^p`.
  2. Order the sweep by decreasing `m` and set `d_k = e_k - e_(k+1)`; the
     order of the ratio `d_k / d_(k+1)` is the `p` solving
     `d_k / d_(k+1) = (m_k^p - m_(k+1)^p) / (m_(k+1)^p - m_(k+2)^p)`
     (the sweeps are not geometric; ratio `k` uses points `k`, `k+1`
     and `k+2`). Take the longest (the first, if tied) run of at
     least two consecutive ratios whose order lies within 0.1 of 1;
     `p` is the median order over the run (for an even count, the
     mean of the two middle values). A run's points are every point
     used by one of its ratios; `C` and `A` come from the two
     smallest of them, i.e. the last two points used by the run's
     last ratio. Omit the Richardson term if there is no run.
  3. Set `tol = max(2|e(m_min)|, 2(|C| + |A| m_min^p), tol_min)`,
     rounded up to one significant figure, where `m_min` is the
     smallest perturbation magnitude probed. `tol_min` is the
     reachable-tolerance floor: `tol_min = sqrt(40 N |A1|)`, where `A1`
     is the `m` coefficient of a least-squares fit
     `e = c0 + A1 m + B m^2` over the four smallest points, and `N` is
     the median of `|e - fit| m` over the six smallest points (every
     point, if the sweep has fewer), `fit` being the same model
     refitted to those points.
  4. Read these from the last `perturbation,F,dFdx,error` block of
     `FD_check_<case>.csv` (the file is appended per run); the
     regression suite (`tests/regression/sensitivity`) applies the
     same rule to its own `FD_check_<case>.csv` or its committed
     `reference_data/ref_FD_check_<case>.csv`.

  Re-derive the tolerance whenever the case's mesh resolution, rank
  count, Neko version, linear-solver tolerance, steady-state tolerance,
  polynomial order or perturbation sweep changes.
- `CMakeLists.txt`: Defines the build process and, via the `test_list`
  variable, registers one CTest per case file.

Currently, we have added tests for the following components:

- `volume_constraint_t` (`volume.case`, `volume_filtered.case`)
- `viscous_dissipation_objective_t` (`viscous_dissipation.case`) — isolates
  the always-on `augmented_lagrangian_objective_t` state-coupling term, since
  this objective's own `update_sensitivity` is empty.
- `brinkman_dissipation_objective_t` (`brinkman_dissipation.case`) — isolates
  the direct-partial-derivative sensitivity path.

`scalar_mixing_objective_t` is not yet covered here (its only existing case,
`tests/regression/sensitivity/cases/passive_scalar.case`, needed a Neko-core
scalar-scheme fix first — see `known-bugs-backlog.md` #9). The heavier,
higher-order, higher-Reynolds-number (and more tightly converged) versions
of these and other cases (`dissipation`, `dissipation_weights`, unsteady
variants) live in `tests/regression/sensitivity/` instead — that suite is
opt-in (`NEKO_TOP_RUN_SENSITIVITY_REGRESSION=1`) and not part of the
default/PR-blocking test budget, since it's too slow to gate every PR.

## Adding New Tests

Reuse the existing generic driver — you almost never need new Fortran code:

1. Add a new `.case` file exercising the objective/constraint you want to
   cover. For a linear, state-independent functional (e.g. the volume
   constraint), keep `time.end_time`/`time.timestep` short — a handful of
   steps, not a physically converged run (see `volume.case`) — so the test
   stays fast. For a PDE-coupled objective (e.g. `viscous_dissipation.case`,
   `brinkman_dissipation.case`), the finite-difference check needs the
   fluid to have reached the same steady state the objective is evaluated
   at, or the comparison is ill-posed; those two cases run to `end_time =
   5.0` because their `steady` simulation component only freezes around
   `t = 2.25` on the default mesh. Tune `end_time` by running the case and
   confirming the fluid has reached steady state; derive `fd_test_tolerance`
   from the per-perturbation `error` output using the rule above, not by
   guessing.
2. Add the case file name to the `test_list` variable in `CMakeLists.txt`:
   ```cmake
   set(test_list
       "volume.case"
       "volume_filtered.case"
       "viscous_dissipation.case"
       "brinkman_dissipation.case"
       "new_case.case"  # Add your new case here
   )
   ```
3. Only write new Fortran (in `problem_tester.f90`/`sensitivity.f90`) if the
   generic driver genuinely can't express what you need — it currently
   handles any single objective or single constraint automatically.
