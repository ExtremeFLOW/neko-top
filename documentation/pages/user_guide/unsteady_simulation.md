# Unsteady simulations {#unsteady_simulation}

Unsteady simulations are simulations which cannot be simply progressed to a
steady state solution which will be used to estimate the entire adjoint
solution. In the case of an unsteady problem, the adjoint simulation require the
entire forward state to be restored to each timestep in reverse order. This
document outline the possible approaches implemented in Neko-TOP.

## Checkpoint restore system

A straight forward approach is to simply store the entire forward state at each
timestep, and then restore it when needed by the adjoint simulation. However,
storing all states in memory is impossible for large simulations, while it is
prohibitively slow to read and write the state to disc for every timestep.
Therefore we support a checkpoint based system.

During the forward simulation, we determine if the state should be placed in
memory or written to disc. The disc files, are restart capable checkpoints
directly from Neko and are used to restart the forward simulation. The memory
checkpoints are a RAM based copy of the required state variables, and not a
complete memory of the forward state.

Current implementation supports the following algorithms:

- Linear: Store a regular interval of restart capable checkpoints to disc. 

By default, we store pressure, velocity, and scalar fields in the checkpoint,
but users can specify extra fields to be included in the checkpoint.

### Parameters

| Parameter          | Type           | Default        | Description                                                               |
| ------------------ | -------------- | -------------- | ------------------------------------------------------------------------- |
| `enabled`          | `bool`         | `false`        | Whether checkpointing is enabled.                                         |
| `algorithm`        | `string`       | `"linear"`     | The checkpointing algorithm to use. Currently only "linear" is supported. |
| `n_memory`         | `int`          | `10`           | The number of checkpoints to store in memory.                             |
| `filename`         | `string`       | `"checkpoint"` | The base name for checkpoint files.                                       |
| `format`           | `string`       | `"chkp"`       | The file format for checkpoint files. Currently "chkp" is supported.      |
| `keep_checkpoints` | `bool`         | `false`        | Whether to keep checkpoint files after the simulation.                    |
| `extra_fields`     | `string array` | `[]`           | A list of extra fields to include in the checkpoint.                      |

### Linear algorithm

The forward simulation will be restarted from the disc checkpoint when the
adjoint simulation requires the state. The forward simulation will be progressed
to the next checkpoint, and the intermediate states will be stored in memory.
