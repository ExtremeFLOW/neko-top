# Documentation for the Neko TOP Repository.

[Examples](./examples.md)

## Long-start

The setup script is an automated setup of Neko, Neko-TOP and other required
dependencies. The script relies in a number of environment variables, which can
be used to modify the behaviour of the system and allow the user to specify
custom install locations for the given dependencies.

| Variable           | Description                                                         | Default               |
| ------------------ | ------------------------------------------------------------------- | --------------------- |
| `NEKO_DIR`         | Location of the Neko library.                                       | external/neko         |
| `JSON_FORTRAN_DIR` | JSON-Fortran library, required dependency of Neko.                  | external/json-fortran |
| `NEK5000_DIR`      | Nek5000, primarily used for meshing and for GSLib.                  | external/Nek5000      |
| `PFUNIT_DIR`       | Unit testing library used in Neko.                                  | external/pFUnit       |
| `CUDA_DIR`         | Location of the CUDA library folders, needed for Nvidia GPU support | -                     |

These can be defined either on the command line by the user or in a `prepare.sh`
file which is loaded by the setup script if it exists in the root of Neko-TOP.
The prepare script provide a convenient way to use module systems such as
`spack` or similar to activate environments and such before compilation.

An example of a `prepare.sh` file is shown below:

```bash
#!/bin/bash
module load cuda/10.1

export CUDA_DIR=$CUDA_HOME
export NEKO_DIR=$HOME/neko
```
