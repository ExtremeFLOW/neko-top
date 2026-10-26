# Case file configuration {#configuration}

In general case files are used as the main input in order to configure a
simulation and optimization problem. Most settings pertaining to the simulation
are defined by Neko, so please refer to the Neko documentation for more
information on these settings.

## Case file structure {#configuration-case-file}

The overall structure of the case file is the following
```json
{
    "version": "1.0",
    "case": { ... },
    "optimization": { 
        "design": { ... },
        "solver": { ... },
        "objectives": [ ... ],
        "constraints": [ ... ]
    }    
}
```

## Optimization Solver settings {#configuration-solver}

Neko-TOP has two optimizers, selected with `optimization.solver.type`:

- `"mma"` — the Method of Moving Asymptotes, which handles both constrained
  and unconstrained problems.
- `"steepest_descent"` — the steepest descent method with move limits (SDM),
  which handles **unconstrained** problems only.

Both read a handful of settings directly from `optimization.solver`, and keep
their own options in a nested object named after the type.

### Common solver settings {#configuration-solver-common}

| Json entry            | Description                                                             | Admissible Values                 | Default value      |
| ---------------------- | ------------------------------------------------------------------------ | ------------------------------------ | --------------------- |
| `type`                | Selects the optimizer.                                                  | `"mma"` or `"steepest_descent"`   | *(required)*        |
| `max_iterations`      | Maximum number of outer optimization iterations.                       | Integer, `> 0`                     | `100`               |
| `enable_output`       | Write the convergence log (`optimization_data.csv`) and the simulation output at every iteration. | `.true.` or `.false.` | `.true.`            |
| `stop_design_change`  | Stop once the maximum design change drops strictly below this value.   | Real, `> 0`                        | machine epsilon     |
| `max_runtime`         | Maximum wall-clock runtime, as a duration string (e.g. `1-12:00:00` for one day and twelve hours). Empty disables the limit. | String                | `""` (disabled)     |

Restart and checkpoint keys (`restart_file`, `checkpoint.*`) are likewise read
the same way by every optimizer; see
\ref configuration-checkpointing "Checkpointing and restart" below.

### The Method of Moving Asymptotes (MMA) {#configuration-solver-mma}

In a case file, MMA is configured in `optimization.solver`, with MMA-specific
settings under `optimization.solver.mma`.

```json
{
  "optimization": {
    "solver": {
      "type": "mma",
      "mma": {
        // MMA options go here
      }
    }
  }
}
```

| Json entry   | Description                                             | Admissible Values     | Default value                                      |
| ------------ | ------------------------------------------------------- | --------------------- | -------------------------------------------------- |
| `max_iter`   | Maximum number of iterations in the MMA subsolver.      | Integer, `> 0`        | `100`                                              |
| `epsimin`    | Minimum barrier parameter used by MMA.                  | Real, `> 0`           | `1.0e-9 * sqrt(m + n_global)`                      |
| `asyinit`    | Initial asymptote distance factor.                      | Real, `> 0`           | `0.2`                                              |
| `asyincr`    | Asymptote increase factor.                              | Real, `> 0`           | `1.05`                                             |
| `asydecr`    | Asymptote decrease factor.                              | Real, `> 0`           | `0.65`                                             |
| `move_limit` | Per-iteration move cap around current design variables. | Real                  | `0.2`                                              |
| `backend`    | MMA execution backend.                                  | `"cpu"` or `"device"` | `"cpu"` (or `"device"` when `NEKO_BCKND_DEVICE=1`) |
| `subsolver`  | Interior-point variant used in MMA subproblems.         | `"dip"` or `"pdip"`   | `"dip"`                                            |
| `xmin`       | Global lower bound applied to all design variables.     | Real                  | `0.0`                                              |
| `xmax`       | Global upper bound applied to all design variables.     | Real                  | `1.0`                                              |
| `a0`         | MMA scalar coefficient \f$a_0\f$.                       | Real                  | `1.0`                                              |
| `a`          | MMA vector coefficient \f$a_i\f$ (applied uniformly).   | Real                  | `0.0`                                              |
| `c`          | MMA vector coefficient \f$c_i\f$ (applied uniformly).   | Real                  | `1000.0`                                           |
| `d`          | MMA vector coefficient \f$d_i\f$ (applied uniformly).   | Real                  | `0.0`                                              |
| `scale`      | Constraint scaling target value.                        | Real, `> 0`           | `1.0`                                              |
| `auto_scale` | Enable adaptive scaling of constraints each iteration.  | `.true.` or `.false.` | `.false.`                                          |

#### MMA configuration options
- `subsolver` selects the interior-point variant used inside MMA (`dip`/`pdip`).
- `backend` selects whether MMA update/KKT operations run on CPU or device
  routines.

#### MMA subproblem and asymptote options
- `max_iter` and `epsimin` define limits for the MMA subsolver iterations and
  barrier parameter.
- `asyinit`, `asyincr`, and `asydecr` control moving-asymptote updates.

#### Design bounds and MMA coefficients
- `a0`, `a`, `c`, and `d` are coefficients in the standard MMA subproblem
  formulation.
- `xmin` and `xmax` set global lower/upper bounds for all design variables.
- `move_limit` restricts each update to a local interval around the current
  iterate. The effective interval used by MMA is `max(xmin, x - move_limit)` to
  `min(xmax, x + move_limit)`. If `move_limit <= 0`, only global bounds are
  used.

#### Scaling options
- `scale` sets the constraint scaling magnitude.
- `auto_scale` enables adaptive scaling from the current constraint value at
  each iteration.

### Steepest descent with move limits {#configuration-solver-steepest-descent}

The steepest descent method with move limits (SDM) is Appendix D of da Silva,
Aage, Beck & Sigmund, "Three-dimensional manufacturing tolerant topology
optimization with hundreds of millions of local stress constraints", IJNME 122
(2021) 548-578, doi:10.1002/nme.6548 (Eqs. D22-D24). It is configured in
`optimization.solver`, with its own settings under
`optimization.solver.steepest_descent`.

```json
{
  "optimization": {
    "solver": {
      "type": "steepest_descent",
      "max_iterations": 200,
      "stop_design_change": 0.01,
      "steepest_descent": {
        "xmin": 0.0,
        "xmax": 1.0,
        "step_length": 1.0,
        "move_limit": {
          "initial": 0.1,
          "min": 0.001,
          "max": 0.1,
          "decrease": 0.7,
          "increase": 1.1
        }
      }
    }
  }
}
```

| Json entry               | Description                                                                    | Admissible Values | Default value |
| -------------------------- | --------------------------------------------------------------------------------- | -------------------- | ---------------- |
| `xmin`                   | Lower bound of the design variables.                                           | Real               | `0.0`          |
| `xmax`                   | Upper bound of the design variables. Must be larger than `xmin`.               | Real               | `1.0`          |
| `step_length`            | Step length \f$\Psi\f$ applied to the normalised direction.                    | Real, `> 0`        | `1.0`          |
| `move_limit.initial`     | Initial per-variable move limit. Must not be smaller than `move_limit.min`.    | Real               | `0.1`          |
| `move_limit.min`         | Smallest allowed move limit.                                                   | Real, `> 0`        | `0.001`        |
| `move_limit.max`         | Largest allowed move limit. Must not be smaller than `move_limit.initial`.     | Real               | `0.1`          |
| `move_limit.decrease`    | Factor applied to a variable's move limit when it oscillates.                 | Real, `(0, 1]`     | `0.7`          |
| `move_limit.increase`    | Factor applied to a variable's move limit when it moves monotonically.        | Real, `>= 1`       | `1.1`          |

Violating any of the relations in the "Admissible Values" column above stops
the run with an error naming the offending key.

#### The update {#configuration-solver-steepest-descent-update}

Let \f$\rho^{(b)}\f$ be the design at iteration \f$b\f$ and
\f$G = \nabla f(\rho^{(b)})\f$ its gradient (the directional derivative of the
objective for a Brinkman-masked design, the raw objective sensitivity
otherwise). Each iteration:

1. **Direction.** The steepest-descent direction \f$S = -G\f$ is projected
   onto the box: a component pointing out of an active bound is zeroed,
   \f[
   S_e = 0 \quad\text{if}\quad
   (\rho_e^{(b)} \ge x^{\max} \text{ and } S_e > 0)
   \quad\text{or}\quad
   (\rho_e^{(b)} \le x^{\min} \text{ and } S_e < 0).
   \f]
2. **Normalisation.** \f$S\f$ is normalised by its own global maximum absolute
   value,
   \f[
   S_{\max} = \max_e |S_e|,
   \f]
   taken over all MPI ranks.
3. **Move limits.** From the third iteration on, each per-variable move limit
   \f$\delta_e\f$ is adapted from the sign of the last two design changes,
   \f[
   \delta_e \leftarrow
   \begin{cases}
   \texttt{move\_limit.decrease} \cdot \delta_e
   & (\rho_e^{(b)} - \rho_e^{(b-1)})(\rho_e^{(b-1)} - \rho_e^{(b-2)}) < 0,\\[4pt]
   \texttt{move\_limit.increase} \cdot \delta_e
   & (\rho_e^{(b)} - \rho_e^{(b-1)})(\rho_e^{(b-1)} - \rho_e^{(b-2)}) > 0,\\[4pt]
   \delta_e & \text{otherwise},
   \end{cases}
   \f]
   clamped to \f$[\texttt{move\_limit.min}, \texttt{move\_limit.max}]\f$
   afterwards. A sign change (oscillation) shrinks the limit;
   monotonic progress grows it.
4. **Update.** The new design is the move-limit window clamped to the box,
   \f[
   \rho^{(b+1)} = \min\!\Bigl(x^{\max}, \max\!\Bigl(x^{\min},
   \min\bigl(\rho^{(b)} + \delta,\,
   \max\bigl(\rho^{(b)} - \delta,\,
   \rho^{(b)} + \Psi \tfrac{S}{S_{\max}}\bigr)\bigr)\Bigr)\Bigr),
   \f]
   with \f$\Psi = \f$ `step_length`. This is exactly Eq. D24 whenever
   \f$\rho^{(b)}\f$ lies in \f$[x^{\min}, x^{\max}]\f$ (see
   \ref configuration-solver-steepest-descent-bounds below for when it does
   not).

#### Only unconstrained problems {#configuration-solver-steepest-descent-unconstrained}

SDM handles box bounds only. If the case file defines any
`optimization.constraints`, initialisation stops with an error pointing to the
`"mma"` optimizer instead.

#### Convergence {#configuration-solver-steepest-descent-convergence}

There is no tolerance-based convergence criterion. A run ends when one of:

1. the step reports "converged" — only when the bound-projected gradient
   \f$S\f$ is exactly zero, i.e. \f$\rho^{(b)}\f$ is a KKT point of the
   box-constrained problem;
2. the maximum design change drops strictly below `stop_design_change`;
3. `max_iterations` or `max_runtime` is reached.

\note The paper's practical stopping rule,
\f$\lVert\Delta\rho\rVert_\infty \le 0.01\f$, corresponds to
`"stop_design_change": 0.01`. The base check is a strict `<`.

\warning Near a solution, SDM tends to settle into a period-2 cycle in which
each variable's move limit is repeatedly shrunk by `move_limit.decrease`,
until it is clamped to `move_limit.min`; the resulting step then sits at, or
just above, `move_limit.min`. If `stop_design_change` is not larger than
`move_limit.min`, criterion 2 above can never fire, and the run continues
until `max_iterations`/`max_runtime` instead.

#### Bounds and the Brinkman mask {#configuration-solver-steepest-descent-bounds}

For a Brinkman design with an optimization-domain mask, points outside that
domain are held at `0` by the design's forward mapping, regardless of `xmin`.
With `xmin > 0` those points therefore sit below `xmin` permanently, and —
because both the current and the updated design are masked to the same value
there — they do not contribute to the reported maximum design change.

An initial design that is infeasible (outside `[xmin, xmax]`) is projected
into the box on the very first step: the box clamp in the update above is
applied *after* the move-limit window, so that first move is not bounded by
`move_limit.initial` and can be larger than `move_limit.max`.

#### The CSV column {#configuration-solver-steepest-descent-csv}

The optimizer adds a `max_abs_gradient` column to `optimization_data.csv`:
the value logged for iteration \f$k\f$ is \f$S_{\max}\f$ of the
bound-projected direction computed at the *start* of iteration \f$k\f$, i.e.
evaluated at the design produced by iteration \f$k - 1\f$. Row 0, logged
before any step is taken, is `0`.

#### Checkpointing and restart {#configuration-solver-steepest-descent-checkpoint}

The design history \f$\rho^{(b-1)}\f$, \f$\rho^{(b-2)}\f$ and the per-variable
move limits \f$\delta\f$ are saved to, and restored from, the checkpoint under
`/SteepestDescent` (see \ref configuration-checkpointing below for the
general checkpoint settings). A restart must use the same number of MPI ranks
as the run that produced the checkpoint.

### Checkpointing and restart {#configuration-checkpointing}

Neko-TOP supports checkpointing and restart of optimization runs. This is
particularly useful for long-running optimization problems, as it allows the
user to save the state of the optimization at regular intervals and restart from
the last checkpoint in case of a failure or if the optimization needs to be
paused. Additionally, a checkpoint can be emitted when a given maximum runtime
is reached, allowing the user to effectively set a maximum runtime for the
optimization.

When enabled, the checkpointing system will emit checkpoint files at regular
`interval`s during the optimization. Each checkpoint will be named according to
the following pattern: `path/base_iter.extension`, where `path` is the path
specified by the user for checkpoint files, `base` is the base name specified by
the user for checkpoint files, `iter` is the iteration number at which the
checkpoint was emitted, and `extension` is the file extension corresponding to
the chosen checkpoint format. A special case is the runtime-based checkpointing,
which will be named as `path/optimizer_rt_checkpoint.extension`.

In order to restart an optimization from a checkpoint, the user can specify the
checkpoint file to restart from in the case file under
`optimization.solver.restart_file`. When this parameter is set, the
optimization will be restarted from the specified checkpoint file. Please note
that the checkpoint file must be compatible with the current case file and
Neko-TOP version in order for the restart to work correctly. If a runtime-based
checkpoint exists and the user has not specified a checkpoint file to restart
from, the optimization will automatically restart from the runtime-based
checkpoint.

On a restart, any \ref continuation "continuation" parameters are first
restored to the values they held at the checkpointed iteration, before the
problem is re-evaluated at the restored design — so the first post-restart
gradient matches the one from the original run.

Parameters related to checkpointing can be set under the `optimization.solver`
section of the case file.

| Parameter             | Type    | Default                | Description                                                                                                                                                                                               |
| --------------------- | ------- | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `restart_file`        | string  | `""`                   | Filename of a checkpoint file to restart from. If this parameter is set, the optimization will be restarted from the specified checkpoint file.                                                           |
| `checkpoint.interval` | integer | -1                     | Interval (in iterations) at which checkpoints are emitted. If set to a positive value, a checkpoint will be emitted every `checkpoint.interval` iterations. If set to -1, no checkpoints will be emitted. |
| `checkpoint.path`     | string  | `"./checkpoints/"`     | Path where checkpoint files will be saved.                                                                                                                                                                 |
| `checkpoint.base`     | string  | `"optimizer_checkpoint"` | Base name for checkpoint files.                                                                                                                                                                          |
| `checkpoint.format`   | string  | `"hdf5"`               | Format for checkpoint files. Accepted spellings are `"h5"`, `"hdf5"`, `"hf5"` and `"hdf"`, all of which write an HDF5 file.                                                                                |

## Additional reading {#configuration-additional}

The individual components are described in greater detail in the linked
sections.

- \subpage design_types
- \subpage mapping_cascade
- \subpage objectives_and_constraints
- \subpage simulation_components
- \subpage unsteady_simulation
