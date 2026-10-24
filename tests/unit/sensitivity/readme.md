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
  matched to their discretisation/steady-state floor (determine this
  empirically per case, don't guess).
- `CMakeLists.txt`: Defines the build process and, via the `test_list`
  variable, registers one CTest per case file.

Currently, we have added tests for the following components:

- `volume_constraint_t` (`volume.case`, `volume_filtered.case`)
- `viscous_dissipation_objective_t` (`viscous_dissipation.case`) — isolates
  the always-on `augmented_lagrangian_objective_t` state-coupling term, since
  this objective's own `update_sensitivity` is empty.
- `brinkman_dissipation_objective_t` (`brinkman_dissipation.case`) — isolates
  the direct-partial-derivative sensitivity path.
- `brinkman_dissipation_objective_t`, unsteady
  (`unsteady_brinkman_dissipation.case`) — 20 time steps at polynomial order
  3, with a PDE filter and the objective's default over-integration. It checks
  that the reported value is integrated with the quadrature its adjoint
  forcing and sensitivity differentiate. With the value on the GLL points and
  its derivatives on the GL points, the error at the smallest perturbation was
  `+1.775e-2`; with one quadrature it is `+1.05e-5`, which is first-order
  truncation (`+5.25e-5` at `5e-4`). The tolerance of `1e-4` lies between the
  two.
- `brinkman_dissipation_objective_t`, unsteady and masked
  (`unsteady_brinkman_dissipation_masked.case`) — the same case with the
  objective masked to the design domain. The PDE filter's transpose carries
  any sensitivity reported outside the mask onto the design variables next to
  the domain boundary, and the probed variable lies `0.015` inside it. With
  the sensitivity unmasked the error was `-6.24e-2`; masked it is `+1.13e-5`.
- `brinkman_dissipation_dealias_mismatch_rejected`
  (`brinkman_dissipation_dealias_mismatch.case`) — the same case with
  `dealias_sensitivity` false and `dealias_forcing` left at its default. A
  mixed setting has no value its sensitivity is the gradient of, so the
  objective must reject it; the test passes only if the run stops with that
  error.

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
   A case the driver must refuse is registered with
   `add_rejected_case(<name> "<expected error>")` instead, which passes only
   if the run stops with that error.
3. Only write new Fortran (in `problem_tester.f90`/`sensitivity.f90`) if the
   generic driver genuinely can't express what you need — it currently
   handles any single objective or single constraint automatically.
