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
