# Newton-Krylov Stability Example

This example uses `neko-top`, `neko`, and
[`LightKrylov`](https://github.com/nekStab/LightKrylov) to compute a steady
baseflow with Newton-Krylov iteration and then estimate leading eigenpairs of the
linearized time-stepper.

The current integration targets the `pr-neko_fix_sk` LightKrylov branch from
[`nekStab/LightKrylov#256`](https://github.com/nekStab/LightKrylov/pull/256).

An additional dependency is [`fpm`](https://github.com/fortran-lang/fpm), which is required to build LightKrylov.

## Running

For cuda builds, from the repository root:

```bash
CUDA_DIR=/path/to/cuda CUDA_ARCH=80 HDF5_DIR=/path/to/hdf5 ./setup.sh -d CUDA
./run.sh newton
```

`setup.sh` installs the required external dependencies into `external/` when they
are missing. This includes LightKrylov, which is installed under
`external/lightkrylov/install` and cached across normal rebuilds.

The mesh can be regenerated from `newton.geo` with `gmsh` if needed.

## Case Options

The LightKrylov controls for this example live in `newton.case`:

```json
"lightkrylov": {
    "newton": {
        "relative_tolerance": 1e-3
    },
    "eigensolver": {
        "number_of_eigenvalues": 15,
        "tolerance": 1e-14
    }
}
```

The driver reads these values with defaults matching the case file above. The
Newton tolerance is passed as `rtol` to `newton`, and the eigensolver tolerance is
passed to `eigs`.

## State Vector

The LightKrylov state is velocity-only:

- `u`
- `v`
- `w`

Pressure is still part of the underlying Neko flow solve, but it is not stored in
the Krylov vector and is not included in vector algebra, Krylov subspace storage,
or eigenvector output.

The state vector implements LightKrylov's `init_like(self, mold)` API. This lets
LightKrylov allocate polymorphic vectors from a mold while `neko-top` handles the
Neko-specific field initialization through the coefficient/dofmap structure.

## Output

After Newton converges, the example writes the baseflow through the normal Neko
forward output path.

Eigenvectors are written by `state_vector_t%write`. These files contain the
velocity block only. The writer explicitly enables velocity output and skips
pressure and temperature in the underlying `.fld` writer.

The eigensolver also writes its convergence information to LightKrylov's standard
text output, which can be inspected or plotted with `plot.py`.