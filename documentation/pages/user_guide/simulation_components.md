# Simulation components {#simulation_components}

\tableofcontents

This page documents the `neko-top` simulation component currently provided.

## Steady simulation component

The `steady` simulation component monitors the change between consecutive
forward fields and freezes the forward fluid solve when the change is below a
user-defined tolerance.

In case files, enable it under `case.simulation_components`:

```json
"simulation_components": [
    {
        "type": "steady",
        "tol": 1e-7,
        "log_frequency": 50,
        "scalar_coupled": false,
        "compute_control": "tsteps",
        "compute_value": 1
    }
]
```

### Convergence measure

For each field, `steady` computes an energy-like norm of the increment
\f$\Delta q = q^n - q^{n-1}\f$,

\f[
\lVert \Delta q \rVert =
\frac{\sqrt{\int_\Omega (\Delta q)^2 \, d\Omega}}{\Delta t\,|\Omega|}.
\f]

The component evaluates this for fluid fields (`u`, `v`, `w`, `p`) and, when
present, one scalar field.

### Stopping behavior

- Fluid convergence is reached when `max(u, v, w, p)` is below `tol`.
- Scalar convergence is checked separately when a scalar is present.
- If `scalar_coupled = true`, both fluid and scalar must satisfy `tol` before
  fluid freezing is applied.
- Once converged, the component sets fluid `freeze = true`.

\note Scalar freezing is not yet applied by this component. The current
implementation only freezes the fluid solver.

### The scalar after the fluid is frozen

Since the scalar is not frozen, it keeps evolving in the frozen flow, and
`steady` keeps computing its increment after freezing the fluid.

With `scalar_coupled = false` the scalar takes no part in the decision to
freeze. The fluid is often frozen long before the scalar converges, more so
the higher `Pe`, so any scalar objective is evaluated on a scalar that may
still be evolving at the end of the run. `steady` therefore issues a warning:

- when it freezes the fluid while the scalar increment is above `tol`, and
- again on its last computation of that run, if the scalar increment is still
  above `tol`. If it has fallen below `tol` by then, a message says so
  instead.

The remedy is `scalar_coupled = true`, with an `end_time` long enough for both
fields to converge. With `scalar_coupled = true` neither warning is issued,
since the fluid is then only frozen once the scalar has converged too.

### Parameters

| Parameter         | Type         | Default   | Description                                                                       |
| ----------------- | ------------ | --------- | --------------------------------------------------------------------------------- |
| `type`            | string       | -         | Must be `"steady"`.                                                               |
| `tol`             | real         | `1e-6`    | Convergence threshold for the field-increment norm.                               |
| `log_frequency`   | integer      | `50`      | Write cadence (in timesteps) for convergence diagnostics.                         |
| `scalar_coupled`  | bool         | `false`   | If `true`, require both fluid and scalar convergence before freezing fluid.       |
| `compute_control` | string       | inherited | Standard Neko simulation-component scheduling control. Typical value: `"tsteps"`. |
| `compute_value`   | integer/real | inherited | Scheduling interval/value associated with `compute_control`. Typical value: `1`.  |

### Output

The component writes convergence data to `steady_state_data.csv` with header:

`iter,time,u,v,w,p,scalar`

where each field column is the normed increment at that logging point.
Rows continue after the fluid is frozen, so the scalar column follows the
scalar to the end of the run. In those rows `u`, `v`, `w` and `p` are exactly
zero, the increment of a frozen field.

### Limitations

- The component currently supports at most one scalar field.
- For truly time-dependent optimization workflows, prefer the unsteady
  configuration described in \subpage unsteady_simulation.