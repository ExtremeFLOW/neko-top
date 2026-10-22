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
  matched to the bias their discretisation and unconverged steady state
  leave (determine this empirically per case, don't guess).
- `CMakeLists.txt`: Defines the build process and, via the `test_list`
  variable, registers one CTest per case file.

Currently, we have added tests for the following components:

- `volume_constraint_t` (`volume.case`, `volume_filtered.case`)
- `viscous_dissipation_objective_t` (`viscous_dissipation.case`) — isolates
  the always-on `augmented_lagrangian_objective_t` state-coupling term, since
  this objective's own `update_sensitivity` is empty.
- `brinkman_dissipation_objective_t` (`brinkman_dissipation.case`) — isolates
  the direct-partial-derivative sensitivity path.
- `viscous_dissipation_variable_dt.case` — the same objective with
  `case.time.variable_timestep`, a regression gate for backlog #55. It fails
  if the driver ignores the `case.time` settings (the run then marches at a
  placeholder `dt = 1.0`), or if one run's last `dt` leaks into the next
  run's time-step history. Its time steps do not depend on the design, by
  construction: the CFL stays far below `target_cfl`, so every run ramps
  from `timestep` by the controller's largest factor, 1.2 per step, until
  `max_timestep` caps it. The finite differences therefore compare runs with
  identical steps. The flow starts from the inflow velocity rather than from
  rest, because Neko's controller divides by the CFL of the first step, and
  a zero CFL stops a Debug build on a floating-point trap (backlog #103).
  The finite-difference estimate converges at first order to a limit about
  0.6% above the adjoint sensitivity. That limit is a bias, not a floor:
  the forward run has not reached a steady state at `end_time`, while the
  steady adjoint assumes it has. Starting from rest doubles it to 1.3%; the
  fixed-timestep `viscous_dissipation.case` has 0.6% too.
- `viscous_dissipation_variable_dt_weighted.case` — the same case with the
  objective weighted by 1e4. The adjoint equations are linear, so its
  relative finite-difference error must equal that of the unweighted case.
  It fails if the adjoint sizes its time steps by the adjoint velocity, whose
  amplitude follows the weight, instead of by the primal velocity that
  convects it.
- `adjoint_joblimit_variable_dt.case`, run by its own driver
  `adjoint_joblimit_tester.f90`. A steady adjoint takes as many steps as the
  forward run did, with time steps of its own, so its clock can stop short
  of `end_time`. It must still count as complete and write no job-limit
  checkpoint. The forward run starts almost at rest with a large `timestep`,
  so its first step is much longer than any adjoint step, which follows the
  CFL of the developed flow. The driver also checks that the adjoint's clock
  did stop short, so that the check cannot pass vacuously.
- The guards against a variable timestep in an unsteady problem
  (`unsteady_variable_timestep.case`) and against a variable timestep
  without the `timestep` that caps its first step
  (`variable_timestep_without_timestep.case`). They are not in `test_list`:
  each passes only if the run stops with its error.

`scalar_mixing_objective_t` is not yet covered here (its only existing case,
`tests/regression/sensitivity/cases/passive_scalar.case`, needed a Neko-core
scalar-scheme fix first — see `known-bugs-backlog.md` #9). The heavier,
realistic (~1000-timestep) versions of these and other cases
(`dissipation`, `dissipation_weights`, unsteady variants) live in
`tests/regression/sensitivity/` instead — that suite is opt-in
(`NEKO_TOP_RUN_SENSITIVITY_REGRESSION=1`) and not part of the default/PR-blocking
test budget, since it's too slow to gate every PR.

## Adding New Tests

Reuse the existing generic driver — you almost never need new Fortran code:

1. Add a new `.case` file exercising the objective/constraint you want to
   cover. Keep `time.end_time`/`time.timestep` short (a handful of steps, not
   a physically converged run — see `volume.case` and
   `viscous_dissipation.case` for the pattern) so the test stays fast; tune
   this and `fd_test_tolerance` empirically by running the case and reading
   the per-perturbation `error` output, not by guessing.
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
