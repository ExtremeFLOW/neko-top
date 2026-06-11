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

### MMA Options overview

The Method of Moving Asymptotes (MMA) is the optimization algorithm currently
used in `neko-top`. In a case file, MMA is configured in `optimization.solver`,
with MMA-specific settings under `optimization.solver.mma`.

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
| `subsolver`  | Interior-point variant used in MMA subproblems.         | `"dip"` or `"dpip"`   | `"dip"`                                            |
| `xmin`       | Global lower bound applied to all design variables.     | Real                  | `0.0`                                              |
| `xmax`       | Global upper bound applied to all design variables.     | Real                  | `1.0`                                              |
| `a0`         | MMA scalar coefficient \f$a_0\f$.                       | Real                  | `1.0`                                              |
| `a`          | MMA vector coefficient \f$a_i\f$ (applied uniformly).   | Real                  | `0.0`                                              |
| `c`          | MMA vector coefficient \f$c_i\f$ (applied uniformly).   | Real                  | `100.0`                                            |
| `d`          | MMA vector coefficient \f$d_i\f$ (applied uniformly).   | Real                  | `0.0`                                              |
| `scale`      | Constraint scaling target value.                        | Real, `> 0`           | `1.0`                                              |
| `auto_scale` | Enable adaptive scaling of constraints each iteration.  | `.true.` or `.false.` | `.false.`                                          |

#### MMA configuration options
- `subsolver` selects the interior-point variant used inside MMA (`dip`/`dpip`).
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
`optimization.solver.checkpoint.file`. When this parameter is set, the
optimization will be restarted from the specified checkpoint file. Please note
that the checkpoint file must be compatible with the current case file and
Neko-TOP version in order for the restart to work correctly. If a runtime-based
checkpoint exists and the user has not specified a checkpoint file to restart
from, the optimization will automatically restart from the runtime-based
checkpoint.

Parameters related to checkpointing can be set under the `optimization.solver`
section of the case file.

| Parameter             | Type    | Default                | Description                                                                                                                                                                                               |
| --------------------- | ------- | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `checkpoint.file`     | string  | -                      | Filename of a checkpoint file to restart from. If this parameter is set, the optimization will be restarted from the specified checkpoint file.                                                           |
| `checkpoint.interval` | integer | -1                     | Interval (in iterations) at which checkpoints are emitted. If set to a positive value, a checkpoint will be emitted every `checkpoint.interval` iterations. If set to -1, no checkpoints will be emitted. |
| `checkpoint.path`     | string  | -                      | Path where checkpoint files will be saved. If not set, checkpoint files will be saved in the current working directory.                                                                                   |
| `checkpoint.base`     | string  | "optimizer_checkpoint" | Base name for checkpoint files.                                                                                                                                                                           |
| `checkpoint.format`   | string  | "h5"                   | Format for checkpoint files. Supported formats are "h5" (HDF5)                                                                                                                                            |

## Additional reading {#configuration-additional}

The individual components are described in greater detail in the linked
sections.

- \subpage design_types
- \subpage mapping_cascade
- \subpage objectives_and_constraints
- \subpage simulation_components
- \subpage unsteady_simulation
