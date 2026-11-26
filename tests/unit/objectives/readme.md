# Objective time-window tests

These tests check that an objective's reported value depends only on the time
window it is accumulated over, and not on how long the simulation happens to
run. They replace the former `examples/time_test`, which covered the same
matrix but could only be inspected by hand. Two steady cases share the mesh
and the inflow, and check what the `steady` simulation component reports
about an unconverged scalar.

The tests are tagged as `unit`, so they are mandatory for CI to pass. Together
they take about 6 s.

## Common files

- `prepare.sh`: builds the box mesh (`genmeshbox`), a short channel with all
  six boundary faces exposed as separate zones.
- `objectives_user.f90`: the user-defined velocity and scalar inflow, used
  by every case here.
- `time_window_tester.f90`: the driver of the time-window tests. It runs the
  case once per entry in `optimization.time_window_test.end_times` and
  compares the objective values across those runs and, optionally, against
  each other within a run. It only takes unsteady problems.
- `steady_tester.f90`: the driver of the steady cases. It runs the case once
  as a steady problem and checks nothing itself; `check_steady_log.cmake`
  reads the log it leaves.
- `step_count.pf`: the pFUnit test for `end_time_step_count` alone, which
  checks where a run stops rather than what its objectives report.

Only the objectives declared in the case file are checked. `problem_t` appends
an internal augmented-Lagrangian objective of its own, which the driver skips.

The drivers call `problem_t%compute` and never `compute_sensitivity`, so no
adjoint is solved. That is what keeps these affordable as unit tests. The one
exception is `end_time_step_count`, below, which solves a single adjoint of a
ten-step run with its own driver.

## The scalar setup, restored from the example

`examples/time_test` drove a paraboloid velocity profile into the duct
carrying two species split either side of a smooth interface, let them leave
through the outlet with every other face insulated and no-slip, and measured
mixing over a zone at the far end. The point was the first iteration of an
optimization: a clear split entering, still a clear split leaving, the design
having done nothing yet.

The conversion to unit tests had no user module to carry either profile, so it
substituted a uniform velocity plug, a zero-flux condition on **all six**
faces for the scalar, a `point_zone` initial blob, and an unmasked objective.
That is a different problem — a conserved blob being stirred, whose only route
to a steady state is slow diffusion, measured over the whole domain.

All four cases here now run the original setup, held in
`objectives_user.f90`. Each driver registers it on `sim%neko_case%user` before
`sim%init`, which works because `user_intf_init` only substitutes its own
defaults for pointers still null. No `makeneko` is involved.

- **Velocity:** `user_velocity` on the inlet, `36 y(y-1) z(z-1)`. The shape is
  matched to the duct: it vanishes on all four walls, so it agrees with the
  no-slip condition applied there, where a uniform plug is discontinuous at
  the inlet edges. The factor 36 makes the mean one, so the duct carries unit
  flow rate — the same as the plug it replaces. Its peak is 2.25, but the CFL
  only moves from 0.29 to 0.36, since the plug had the steeper near-wall
  gradients.
- **Scalar:** `user_dirichlet` on the inlet, zero-flux Neumann on the other
  five zones, and the same profile as the initial condition so no run opens
  with a transient it then has to sit through.
- **Objective:** `scalar_mixing` masked to `outlet_region`, the equivalent of
  the example's `objective_domain`.

The split uses a logistic profile in `z`, as the example did, but with a
steepness of 20 rather than 200. The example ran 32x8x8 at polynomial order 5;
these tests run 8x4x4 at order 3, thirteen points across the split, where 200
would be a step function in all but name and would ring.

For that profile a completely unmixed scalar gives a `scalar_mixing` objective
of `0.1000`, and a uniformly mixed one gives `0`. The time-window cases
report about `0.099` at the outlet — the split fully intact, since the flow
has barely developed.

## `time_window_run_length`

The regression guard. Every objective is given the closed window
`[0.02, 0.04]`, which lies inside both runs, so each objective measures the
same interval whether the run stops at `t = 0.05` or continues to `t = 0.1`.
The values must therefore be identical.

This is the property that was broken before: accumulation was normalised by
the *simulation's* window rather than the objective's own, so doubling the run
halved every windowed objective. Against that code all three objectives fail
here with a relative difference of exactly `0.5` — re-checked after the
scalar setup was restored, by putting the old normalisation back into
`functional_accumulate_value` and rerunning. The property is independent of
what the objectives are worth, so the restored physics neither strengthens
nor weakens this guard.

The tolerance is `1e-9`. Two independent runs of the same case agree to about
`1e-11` — the residual is the iterative solvers, not the windowing — so this
leaves ample headroom while staying far below the factor of two a
normalisation regression would produce.

## `time_window_equivalence`

The semantic check, covering the four ways a window can be written. The
objectives are read as consecutive pairs, each pair selecting the same samples
by different means:

| pair | objective            | first form        | second form                                     |
| ---- | -------------------- | ----------------- | ----------------------------------------------- |
| 1    | viscous dissipation  | `end_time` only   | the same window written out in full             |
| 2    | Brinkman dissipation | no window at all  | an explicit `start_time` of zero                |
| 3    | scalar mixing        | `start_time` only | a closed window running past the end of the run |

Pair 3 also exercises a window clipped by the end of the run: its `end_time`
of `0.06` is deliberately beyond the run's own `end_time` of `0.05`.

Note that this test cannot catch a normalisation error on its own — every
objective in a single run shares the same divisor, so a global rescaling
cancels out of a within-run comparison. It passes against the pre-fix code.
`time_window_run_length` is what guards that; this one pins down what each
window form means.

## `end_time_step_count`

The run-length guard. It runs `time_window_run_length.case`, whose
`end_time = 0.05` lies on a step boundary for its `dt = 0.005`, once through
each branch of `problem_t%compute`, and requires each run to take exactly ten
steps and stop at `t = 0.05`. Its pFUnit test, `step_count.pf`, checks nothing
but where the runs stop, and refuses a case whose
`end_time` is not on a step boundary, since only there does the check bite. It
costs 4 s.

This is the property that was broken before. The time accumulated after ten
steps is `0.049999999999999996`, and the unsteady path's loop compared it
against `end_time` with no tolerance, so it took an eleventh step. The steady
path already stopped through `time_state_t%is_done`, which allows for that
round-off, so the two paths disagreed by a step. Against that code the test
fails with

```
The unsteady run took 11 steps to t =  0.550000000000000E-01 instead of 10 to t =  0.500000000000000E-01
```

while the steady run passes. Both paths are reported before the test fails,
so a regression in either is named.

The driver then solves the adjoint of the unsteady run once and requires that
it writes no job-limit checkpoint. With `output_at_end` off, which this case
has, the adjoint's finalisation writes one whenever it judges the run to have
stopped short of `end_time`. The adjoint replays the forward run's step count,
so its clock ends where the forward run's did, a few ulps short of `end_time`
here, and it needs the same tolerance to see the run as complete. Comparing
there without one, as it did before, made every such run write a spurious
`joblimit00000.chkp` once the forward loop was fixed. With that comparison put
back, the step counts pass and the test fails with

```
The adjoint wrote joblimit00000.chkp, so it judged the run to have stopped short of end_time
```

The file name is built the way the adjoint builds it, so a change to Neko's
checkpoint naming cannot leave the test looking for the wrong file, and any
copy left by an earlier failure is deleted before the adjoint runs.

## `steady_scalar_unconverged` and `steady_scalar_unconverged_simtime`

The `steady` freeze trap, run on purpose.
`steady` freezes only the fluid, and with `scalar_coupled` off it freezes it on
the fluid's residuals alone, so a scalar objective can be evaluated on a scalar
that is still evolving. Nothing used to say so.

Both are steady problems (`"unsteady": false`), run once each by
`steady_tester`. The `steady` component belongs to steady problems only; an
unsteady problem that declares it stops at setup, which
`tests/unit/sensitivity/unsteady_steady_component.case` guards. The cases
declare no objectives, since the warnings are what is under test.

The settings are chosen for the fluid to converge quickly while the scalar
does not: `Re = 5`, so the fluid residual reaches the `tol` of `1e-6` by
`t = 1.29`; `Pe = 100`, so the scalar has no diffusion to settle it quickly;
`timestep = 0.01`, since `0.02` makes the scalar diverge at that `Pe`;
`f_max = 50`, keeping `chi * dt` at `0.5`; and the scalar solver at `1e-12`.

`steady_scalar_unconverged` runs to `end_time = 1.495`, under `tsteps` control
at every step. The fluid is frozen at `t = 1.29`, when the scalar residual is
still `1.5e-2`, and the run ends at `t = 1.5` with the scalar still far from
converged. `steady_scalar_unconverged_steady_log` then reads the log and
requires exactly one warning at the freeze and exactly one at the end of the
run, carrying the time of the component's last check, here the final step at
`t = 1.5`, and no report that the scalar converged. An exact count catches a
warning repeated every step, and the time catches an end-of-run warning
issued early.

`steady_scalar_unconverged_simtime` is the same trap under time-based compute
control: `compute_control` is `simulationtime` with a `compute_value` of
`0.1`, and `end_time` is moved off the step grid to `2.495`. The run ends at
`t = 2.5`, but the last check is at `t = 2.4`, since the scheduled `2.5` lies
past `end_time` and the controller skips it. The end-of-run report must
therefore come at `t = 2.4`. The fluid is frozen at `t = 1.6`, later than
under `tsteps`, because each increment then spans ten steps but is still
divided by one `dt`.

The expected counts are one of each warning and no converged report. Until
these cases were made steady problems, each was run twice, through the steady
and the unsteady path, and expected two of each. With those counts the checks
now fail (`FROZEN: found 1, expected 2`), as does the `simtime` check with the
end time moved to `2.5` (`AT_END: found 0, expected 1`).

All three checks were seen to fail when they were written. Against the code
before the warnings were added, the two `tsteps` checks report none of them.
With the end-of-run check moved to the first step after the freeze, the
end-of-run warning carries `t = 1.3` and the `steady_scalar_unconverged` check
fails, while the case itself still passes. The first version of the end-of-run
check predicted the next computation from the compute interval, which only
matches a `tsteps` schedule. Under `simulationtime` it saw a computation still
to come at `t = 2.4`, there was none, and no end-of-run report was issued at
all, so the `simtime` check failed. The check now asks the compute controller
itself.

The check is a separate test, `check_steady_log.cmake` run as a fixture of
the case's own test, because the warnings only exist in the log. The driver
has no handle on the simulation component that issues them.

## Choosing window boundaries in new cases

Prefer a closed `end_time` strictly inside the run rather than exactly equal
to the run's own `end_time`. Objectives hand their window to the adjoint
forcing source terms, and `source_term_t`'s gate compares against the
accumulated simulation time with no tolerance, so a boundary landing exactly
on the final step can silently lose that step's forcing.

Keep every window boundary away from a step time, by roughly half a timestep.
The accumulation gate has a `1e-6 * dt` tolerance, so a boundary sitting on a
step is decided by round-off in the accumulated time rather than by the case
file. The simulation's own `end_time` no longer needs that care: a run stops
on a step-boundary `end_time` rather than one step past it, on both paths —
see `end_time_step_count` above.

To select the final step alone, set `start_time` to the simulation's
`end_time` and leave the objective's `end_time` unset. The run always overshoots
a non-boundary `end_time` by exactly one step, so that window holds one sample.
