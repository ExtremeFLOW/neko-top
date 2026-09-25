# Objective time-window tests

These tests check that an objective's reported value depends only on the time
window it is accumulated over, and not on how long the simulation happens to
run, nor on which of the two code paths — steady or unsteady — evaluated it.
They replace the former `examples/time_test`, which covered the same matrix but
could only be inspected by hand.

The tests are tagged as `unit`, so they are mandatory for CI to pass. Together
they take about 44 s, almost all of it in `steady_unsteady_converged`, which
has to run a flow to a genuine steady state.

## Common files

- `prepare.sh`: builds the box mesh (`genmeshbox`), a short channel with all
  six boundary faces exposed as separate zones.
- `objectives_user.f90`: the user-defined scalar inflow, used by
  `steady_unsteady_converged` and ignored by the cases whose boundary
  conditions do not ask for it.
- `time_window_tester.f90`: the driver. It runs the case once per entry in
  `optimization.time_window_test.end_times` and compares the objective values
  across those runs and, optionally, against each other within a run.

With `compare_steady` set it also runs the case with `simulation_t%unsteady`
cleared — the flag `problem_t%compute` branches on — and requires that result
to agree with the unsteady ones. `require_steady_state` additionally asserts
that every run reached a steady state, which is what makes a window average
comparable to a single converged field at all.

Only the objectives declared in the case file are checked. `problem_t` appends
an internal augmented-Lagrangian objective of its own, which the driver skips.

The driver calls `problem_t%compute` and never `compute_sensitivity`, so no
adjoint is solved. That is what keeps these affordable as unit tests.

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
`objectives_user.f90`. The driver registers it on `sim%neko_case%user` before
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
of `0.1000`, and a uniformly mixed one gives `0`. The three short cases all
report about `0.099` at the outlet — the split fully intact, since the flow
has barely developed — and `steady_unsteady_converged`, which runs to a steady
state, reports `0.0648`, or 65% of the split still there.

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

| pair | objective | first form | second form |
|------|-----------|------------|-------------|
| 1 | viscous dissipation | `end_time` only | the same window written out in full |
| 2 | Brinkman dissipation | no window at all | an explicit `start_time` of zero |
| 3 | scalar mixing | `start_time` only | a closed window running past the end of the run |

Pair 3 also exercises a window clipped by the end of the run: its `end_time`
of `0.06` is deliberately beyond the run's own `end_time` of `0.05`.

Note that this test cannot catch a normalisation error on its own — every
objective in a single run shares the same divisor, so a global rescaling
cancels out of a within-run comparison. It passes against the pre-fix code.
`time_window_run_length` is what guards that; this one pins down what each
window form means.

## `steady_unsteady_converged`

The steady and the unsteady approach are two ways of putting a number on the
same problem, and once that problem has reached a steady state they must give
the same number. This test is that statement: the flow is run to convergence,
the steady path evaluates each objective on the converged field, the unsteady
path averages the same objective over a window lying inside the converged
tail, and the two are required to agree.

Nothing checked this before. The steady path was reachable from
`tests/unit/sensitivity` but never compared against the unsteady one, so the
two were free to disagree about which field they evaluate, or when.

### The split has to survive to the outlet

This is the property that decides `Pe`, and it is worth stating as a number
rather than an intention. For the `k = 20` profile, a scalar that reached the
outlet completely unmixed would give a `scalar_mixing` objective of `0.1000`;
one mixed to uniformity would give `0`. Measured over the whole domain:

| `Pe` | objective | fraction of the split preserved |
|------|-----------|---------------------------------|
| `1` | `0.00869` | 9% |
| `100` | `0.07725` | 77% |

At `Pe = 1` the split is gone by the outlet, which is the opposite of the
first-iteration state the case is supposed to represent. `Pe = 100`, the value
the sibling cases already use, keeps it. Restricted to the outlet fifth of the
domain the figure is 65%, so the two streams are still clearly separated where
it matters. (The `Pe = 1` figure was measured against the plug inflow, before
the paraboloid was restored; the paraboloid shears a little harder, so every
number here is slightly lower than it was.)

This is the only case here carrying `scalar_mixing` twice, once unmasked and
once masked to `outlet_region`. The unmasked one is what the `Pe` table above
is measured on; the masked one is what the example actually optimized.

### Why the rest of the case looks the way it does

The boundary conditions and `Pe` are the example's. The rest is this test's,
each number measured rather than inherited.

| setting | value | why |
|---------|-------|-----|
| `Re` | `5.0` | At the sibling cases' `Re = 50` the fluid residual is still `9e-6` at `t = 3.0`, 39 s of run time short of the `1e-6` it needs. At `Re = 5` it reaches that by `t ≈ 1.27`. |
| `timestep` | `0.01` | Not free. At `Pe = 100` there is no diffusion to damp grid-scale advection, and `dt = 0.02` — a perfectly comfortable CFL of 0.58 for the fluid — makes the *scalar* diverge, reaching a residual of `6e36`. It also saves less than it looks, since larger steps need more solver iterations: 600 steps at `0.02` cost 28 s against 33 s for 1000 at `0.01`. |
| `f_max` | `50.0` | Halved alongside the timestep, so the Brinkman penalty keeps the `chi * dt` of `0.5` the other cases run at. |
| `end_time` | `19.995` | 2000 steps, ending at `t = 20.0`, half a step off a boundary. See below — this is set by how long the scalar takes to stop drifting, not by the fluid. |
| scalar `absolute_tolerance` | `1e-12` | The solver's own noise, not the physics, was the floor. See below. |
| `scalar_coupled` | `false` | So the fluid freezes as soon as *it* has converged, at `t ≈ 1.27`, rather than waiting 18 more time units for the scalar. See below. |

`require_steady_state` makes the driver assert that each run really did
converge, by checking that `steady_simcomp` froze the fluid.

### Why the run is 20 time units long, and why that is still cheap

The fluid settles almost immediately. The scalar does not: at `Pe = 100` its
approach to steady state is advective flushing past the no-slip walls, which
e-folds about every 0.74 time units, and it starts five decades away. Measured
directly, as the disagreement between the two paths for the same objective
windowed at three places in one run:

| window | relative difference |
|--------|---------------------|
| `[8, 9]` | `1.7e-5` |
| `[11, 12]` | `4.2e-7` |
| `[14, 15]` | `4.8e-9` |

Clean exponential decay, which is why the window sits at `[19, 20]`.

The length is affordable because `scalar_coupled` is off. The fluid is
genuinely converged at `t ≈ 1.27`, `steady_simcomp` freezes it there, and the
remaining 1870 steps skip the fluid solve entirely — around 10 ms a step
against 33 ms unfrozen. Two runs of 2000 steps cost 35 s. Leaving
`scalar_coupled` on would keep re-solving an already-converged fluid for
another 18 time units at four times the cost, for the same answer.

That is the one place this case departs from the advice in item 24 of the
workspace bug backlog, and deliberately: it does not use the freeze to
certify that the scalar has settled. The comparison itself does that, far more
strictly than a residual threshold would.

### The scalar solver tolerance was a floor, not the physics

Worth its own note, because it looked like a convergence problem and was not.
With the scalar solver at the `1e-9` this case inherited, the two paths would
not agree better than about `7.5e-10` however long the run went — lengthening
the window moved it far less than the residual had dropped over the same
stretch. That plateau is per-step solver noise. The steady path reads one
field; the unsteady path averages a hundred of them, so the noise partly
averages out of one side and not the other, and no amount of extra convergence
closes the gap. Tightening the scalar solver to `1e-12` removed it.

### What it measures

The fluid converges and freezes at `t ≈ 1.27`; the window is `[19, 20]`, 101
samples.

| objective | steady | unsteady average | relative difference |
|-----------|--------|------------------|---------------------|
| viscous dissipation | `3.29060805369292` | `3.29060805369292` | `0` |
| Brinkman dissipation | `0.575829630750006` | `0.575829630750006` | `0` |
| scalar mixing, whole domain | `0.0772496512776399` | `0.0772496512779822` | `4.4e-12` |
| scalar mixing, outlet | `0.0647557680999154` | `0.0647557681012532` | `2.1e-11` |

The two velocity-based objectives are bit-identical because
`steady_simcomp` freezes the fluid on convergence, so the velocity field is
literally constant across the window and the average of a constant is that
constant. It also means the test checks its own margin — had the freeze landed
inside the window, these two would stop agreeing.

The scalar is converged but never frozen (freezing it is an open `@todo` in
`steady_simcomp`), so it keeps relaxing through the window. The two scalar
figures are what is left of that, against a `1e-9` tolerance, the same one the
other tests here use — 226x and 48x of margin.

The test costs 35 s, by far the most expensive in this directory, which is the
price of a scalar that both reaches a steady state and still looks like the
problem it is meant to represent.

## `steady_unsteady_final_step`

The same steady-versus-unsteady comparison with convergence taken out of it.
The case is `time_window_run_length` unchanged apart from a `steady`
simulation component and a window holding the final timestep alone; a
one-sample average is that sample, so the two paths must agree exactly
whatever the flow is doing. It costs 3 s.

This exists to localise a failure of `steady_unsteady_converged`. If both
fail, the two paths disagree about which field they evaluate. If only the
converged one fails, the paths are fine and the run is not reaching a steady
state. An unconverged run is also the stricter comparison of the two: a frozen
field gives the same objective at every step near the end, so a path sampling
the wrong step would slip through, whereas a field still in motion catches it.

### Ordering, and why it matters

Both tests run the steady path first. Every run leaves its values in the same
objectives, so a steady path that quietly stopped evaluating them would,
running second, still be holding the unsteady run's numbers and would agree
with it. Running first it reports the objectives' initial zero instead.
Checked by deleting the `update_objectives` call from `problem_compute`'s
steady branch: the test fails with a relative difference of exactly `1.0`.
With the other order it passed.

### Step boundaries

`steady_unsteady_final_step` runs to `end_time = 0.0475` against a `dt` of
`0.005`, deliberately off a step boundary. The time loop stops at the first
step reaching `end_time`, so the run takes ten steps and finishes at
`t = 0.05`, and every objective's `start_time` of `0.0475` admits that step
and no other. Both ends of the accepted interval sit half a timestep from a
sample, so no amount of round-off in the accumulated time can change which
steps are counted.

An `end_time` sitting *on* a step boundary is what to avoid, and not for a
subtle reason: with `end_time = 0.05` and `dt = 0.005` the accumulated time
after ten steps lands just below `0.05`, the loop takes an eleventh step, and
the run ends at `t = 0.055`. The window then holds two samples while the
steady path still evaluates one, and the objectives disagree by 0.3%, 10% and
0.2% — a real failure with a thoroughly misleading cause. The underlying
overshoot is Neko's, not this test's.

Measured agreement between the two paths, with the window right, is `4e-13`,
`4e-12` and `6e-14` relative, the same floor two identical unsteady runs of
this case reach.

## Choosing window boundaries in new cases

Prefer a closed `end_time` strictly inside the run rather than exactly equal
to the run's own `end_time`. Objectives hand their window to the adjoint
forcing source terms, and `source_term_t`'s gate compares against the
accumulated simulation time with no tolerance, so a boundary landing exactly
on the final step can silently lose that step's forcing.

Keep every window boundary away from a step time, by roughly half a timestep.
The accumulation gate has a `1e-6 * dt` tolerance, so a boundary sitting on a
step is decided by round-off in the accumulated time rather than by the case
file. The same applies to the simulation's own `end_time`, which decides how
many steps the run takes — see `steady_unsteady_equivalence` above.

To select the final step alone, set `start_time` to the simulation's
`end_time` and leave the objective's `end_time` unset. The run always overshoots
a non-boundary `end_time` by exactly one step, so that window holds one sample.
