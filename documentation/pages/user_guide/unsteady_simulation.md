# Unsteady simulations {#unsteady_simulation}

Unsteady simulations are simulations which cannot be simply progressed to a
steady state solution which will be used to estimate the entire adjoint
solution. In the case of an unsteady problem, the adjoint simulation requires
the forward state to be restored at each timestep in reverse order. Neko-TOP
handles this through the `state_recovery_t` abstraction. A state recovery
implementation is responsible for storing enough forward information and then
restoring or reconstructing the forward state when the adjoint run asks for it.

At the user level, unsteady problems declare a `state_recovery` block at the
top level of the case file:

```json
"unsteady": true,
"state_recovery": {
  "type": "checkpoint"
}
```

The available implementation is `checkpoint`, which restarts the forward
problem from stored checkpoints.

## Checkpoint state recovery

A straight forward approach is to simply store the entire forward state at each
timestep, and then restore it when needed by the adjoint simulation. However,
storing all states in memory is impossible for large simulations, while it is
prohibitively slow to read and write the state to disc for every timestep.
Therefore we support a checkpoint based system.

During the forward simulation, we determine if the state should be placed in
memory or written to disc. The disc files are restart capable checkpoints
directly from Neko and are used to restart the forward simulation. The memory
checkpoints are a RAM based copy of the required state variables, and not a
complete memory of the forward state.

The current checkpoint implementation supports the following algorithms:

- Linear: Store a regular interval of restart capable checkpoints to disc.

By default, Neko-TOP stores pressure, velocity, and scalar fields in the
checkpoint, but users can specify extra fields to be included in the
checkpoint.

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
adjoint simulation requires the state. The forward simulation will be
progressed to the next checkpoint, and the intermediate states will be stored
in memory.

## POD state recovery

The POD state recovery option replaces forward-state checkpointing with an
in-situ reduced-order workflow. During the forward simulation, Neko-TOP streams
snapshots to a Python driver. The Python side uses
[`pySEMTools`](https://github.com/ExtremeFLOW/pySEMTools) to build a POD basis
and store the time coefficients. During the adjoint run, the reduced model is
used to reconstruct the forward state when it is needed.

This reduces the amount of stored state substantially, but it introduces extra
software requirements compared with checkpoint recovery and it changes the
runtime model from a Neko-only MPI job to a coupled Neko + Python MPI launch.

At the case-file level, POD recovery looks like:

```json
"unsteady": true,
"state_recovery": {
  "type": "pod",
  "batch_size": 10,
  "n_modes": 5,
  "i_stream": 10,
  "dtype": "double",
  "output_reconstruction": true,
  "write_modes": true
}
```

### Parameters

| Parameter                         | Type      | Default                | Description |
| --------------------------------- | --------- | ---------------------- | ----------- |
| `type`                            | `string`  | `"checkpoint"`         | Selects the recovery implementation. Use `"pod"` for POD state recovery. |
| `batch_size`                      | `int`     | required               | Number of streamed snapshots accumulated by the Python POD driver before each POD update. |
| `n_modes`                         | `int`     | required               | Number of POD modes kept for the reduced basis. |
| `i_stream`                        | `int`     | required               | Snapshot stride in timesteps. Every `i_stream` steps is sent to the Python driver. |
| `dtype`                           | `string`  | `"double"`             | Floating-point precision used in the POD driver. Supported values are `"single"` and `"double"`. |
| `write_modes`                     | `bool`    | `false`                | Whether to write the final POD modes from the Neko side. |
| `output_reconstruction`           | `bool`    | `false`                | Whether to write reconstructed fields during replay. |
| `output_precision`                | `string`  | `"sp"`                 | Precision used when writing POD modes. |
| `output_format`                   | `string`  | `"fld"`                | Output format used when writing POD modes. |
| `output_file_name`                | `string`  | `"POD_modes"`          | Base name for POD mode output. |
| `reconstruction_output_precision` | `string`  | inherited from output  | Overrides reconstruction output precision when set. Otherwise reconstruction inherits the POD output precision, or falls back to `case.output_precision` when POD output precision was not set explicitly. |
| `reconstruction_output_format`    | `string`  | inherited from output  | Format used for reconstructed-field output. |
| `reconstruction_output_file_name` | `string`  | `"pod_reconstruction"` | Base name for reconstructed-field output. |
| `debug`                           | `bool`    | `false`                | Enables verbose control-stream logging for the coupled Neko/Python run. |

When `output_reconstruction = true`, the reconstruction output cadence is taken
from `case.fluid.output_control` and `case.fluid.output_value`.

## POD technical setup

### ADIOS2 requirement

POD state recovery depends on ADIOS2. The Fortran side uses the Neko-TOP
streaming layer, and the Python side imports ADIOS2-backed
[`pySEMTools`](https://github.com/ExtremeFLOW/pySEMTools) streaming support.
The POD build requires ADIOS2; CMake reports a configuration error when it is
unavailable.

The build and runtime helpers try to locate ADIOS2 from:

- `ADIOS2_PATH`
- `ADIOS2_DIR`
- `${repo}/external/adios2`
- `adios2-config` on `PATH`

If ADIOS2 is not already installed, `scripts/dependencies.sh` can build it.
The relevant configuration is:

- `NEKO_WITH_ADIOS2=true`
- `ADIOS2_DIR=/path/to/adios2` if you want to point at an existing install
- optional ADIOS2 build options such as `ADIOS2_ENABLE_FORTRAN=ON`,
  `ADIOS2_ENABLE_PYTHON=ON`, and `ADIOS2_ENABLE_SST=ON`

The runtime helper in `scripts/mpmd_run_helpers.sh` also extends
`PYTHONPATH` and `LD_LIBRARY_PATH` so that `adios2.bindings` can be imported
on the Python side.

### Python requirements

The POD Python driver is `scripts/python/pod_state_recover.py`. The runtime
dependencies are:

- `numpy`
- `mpi4py`
- `adios2.bindings`
- [`pySEMTools`](https://github.com/ExtremeFLOW/pySEMTools)

More specifically, the driver imports:

- `pysemtools.datatypes.coef.Coef`
- `pysemtools.datatypes.msh.Mesh`
- `pysemtools.io.adios2.stream.DataStreamer`
- `pysemtools.io.utils.get_fld_from_ndarray`
- `pysemtools.rom.io_help.IoHelp`
- `pysemtools.rom.pod.POD`

These modules come from the `pySEMTools` project linked above.

If you are doing POD post-processing or plotting, some example scripts also use
packages such as `matplotlib` and `h5py`, but those are not required for the
runtime POD driver itself.

### MPI and Python environment

Because the POD workflow uses `mpi4py` and a shared MPI launch with both the
Python driver and the Neko executable, the Python environment must be
consistent with the MPI stack used to run Neko. This is the part most likely to
fail if packages are mixed from different environments.

The recommended workflow is:

1. Build Neko-TOP with ADIOS2 enabled.
2. Create or activate a conda environment that contains `mpi4py`, `numpy`,
   and the Python dependencies required by `pySEMTools`.
3. Use the example `prepare.sh` scripts to build the case in the same software
   environment you intend to run.
4. Launch the POD case through the repository-level `run.sh` script from the
   Neko-TOP root, for example `./run.sh POD_rugby_ball`, rather than manually
   piecing together the MPI command.

The POD helper scripts print the active `python3`, `mpirun`, `CONDA_PREFIX`,
`PYTHONPATH`, and `LD_LIBRARY_PATH` before launching. If `mpi4py` or ADIOS2
import errors occur, or if the Python side appears to hang at startup, the
first thing to check is whether the Python environment and the MPI launcher
come from the same stack.
