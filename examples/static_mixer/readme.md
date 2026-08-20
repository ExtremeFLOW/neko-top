# Static mixer example {#static-mixers}

This example simulates flow through a static mixer - a duct containing a fixed
solid insert that passively mixes two fluid streams. The geometry is the design
by Casper Schousboe Andreasen. The example demonstrates the use of the Brinkman
immersed boundary method in [Neko](https://github.com/ExtremeFLOW/neko) to
represent the solid mixer structure on a structured background mesh.

## Physics

The simulation solves the incompressible Navier-Stokes equations coupled with a
passive scalar transport equation (labelled `temperature`) using the pN/pN
spectral-element scheme with third-order time integration. The key
non-dimensional parameters are the Reynolds number **Re** and the Peclet number
**Pe** (set equal to Re by default).

### Domain and mesh

| Parameter        | Value                                           |
| ---------------- | ----------------------------------------------- |
| Domain           | Rectangular duct, 4 x 1 x 1                     |
| Mesh resolution  | `Nx x Ny x Nz` elements, where `Ny = Nz = Nx/4` |
| Polynomial order | 7 (GLL nodes per element)                       |

### Boundary and initial conditions

**Velocity**

| Boundary                                | Condition                                                                                                                                                                                                                      |
| --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Inlet (`x = 0`, zone 1)                 | Smooth parabolic profile - smooth-step blending in `y` and `z`, `u_max = 1`, `v = w = 0`. The boundary layer thickness is `max(1.0 / sqrt(Re), element_size)`, ensuring the blending region is at least one mesh element wide. |
| Outlet (`x = 4`, zone 2)                | Pressure outflow (outflow+dong)                                                                                                                                                                                                |
| Walls (`y = 0/1`, `z = 0/1`, zones 3-6) | No-slip                                                                                                                                                                                                                        |

**Temperature (passive scalar)**

| Boundary                             | Condition                                                                                        |
| ------------------------------------ | ------------------------------------------------------------------------------------------------ |
| Inlet (`x = 0`, zone 1)              | Dirichlet: smooth-step transition from `T = 0` (bottom, `z < 0.45`) to `T = 1` (top, `z > 0.55`) |
| All remaining boundaries (zones 2-6) | Zero-flux Neumann                                                                                |

**Initial conditions**

- **Velocity**: `u = 1` everywhere outside the solid (the Brinkman indicator is
  subtracted so that velocity is zero inside the mixer structure), `v = w = 0`.
- **Temperature**: uniform `T = 0`.

## Brinkman immersed boundary method

The solid mixer geometry is represented using a Brinkman penalization source
term rather than body-conforming mesh boundaries. A scalar indicator field
\f$ \chi(\mathbf{x}) \in [0, 1] \f$ is constructed from the signed distance to
the mixer surface (provided as an STL file), where \f$ \chi = 1 \f$ marks solid
and \f$ \chi = 0 \f$ marks fluid.

The indicator is mapped to a local penalization coefficient \f$ k(\chi) \f$
using the RAMP (Rational Approximation of Material Properties) scheme:

\f[
  k(\chi) = k_0 + (k_1 - k_0)\, \chi\, \frac{q + 1}{q + \chi}
\f]

where \f$ k_0 \f$ = `limits[0]` (minimum, fluid region), \f$ k_1 \f$ =
`limits[1]` (maximum, solid region), and \f$ q \f$ = `penalty` is the RAMP
convexity parameter. The Brinkman source term then adds a volumetric resistance
force to the Navier-Stokes momentum equation:

\f[
  \mathbf{f} = -k(\chi)\, \mathbf{u}
\f]

When \f$ k \to \infty \f$ inside the solid region the velocity is driven to
zero, recovering a no-slip condition on the embedded surface. In practice
\f$ k_1 \f$ is set to `1e6 / Re` following the rule-of-thumb that the Brinkman
layer thickness scales as \f$ 1/\sqrt{k_1} \f$.

### Implementation in Neko

The Brinkman term is added to the fluid solver through the `source_terms` array
in the case file:

```json
"source_terms": [
    {
        "type": "brinkman",
        "limits": [0.0, <eta_max>],
        "penalty": 3.0,
        "output": { "enable": true, "format": "vtkhdf" },
        "filter": {
            "type": "PDE",
            "radius": <r>
        },
        "objects": [
            {
                "type": "boundary_mesh",
                "name": "<path/to/mixer.stl>",
                "cache": true,
                "cache_file": "<path/to/cache_>",
                "distance_transform": {
                    "type": "smooth_step",
                    "value": 0.01
                }
            }
        ]
    }
]
```

Key parameters:

| Parameter                  | Description                                                                                              |
| -------------------------- | -------------------------------------------------------------------------------------------------------- |
| `limits`                   | `[k_0, k_1]` - RAMP penalization range; `k_1 = 1e6 / Re`                                                 |
| `penalty`                  | RAMP convexity parameter \f$ q \f$ controlling the sharpness of the solid-fluid transition (default 3.0) |
| `filter.radius`            | PDE filter length scale; set to `4.0 / Nx` so it spans one mesh element in the x direction               |
| `distance_transform.value` | Half-width of the smooth-step interface layer                                                            |

At startup, Neko reads the STL, computes the signed distance field, applies the
smooth-step transform to obtain \f$ \chi \f$, smooths it with a PDE filter, and
caches the result. The indicator field is subsequently available in user Fortran
code via:

```fortran
brinkman_indicator => neko_registry%get_field("brinkman_indicator")
```

This allows the initial condition and any other user routines to mask off the
fluid region inside the solid, as shown in `petsc/mixer.f90`.

@note The STL file for the mixer geometry (`Test_mixer_095_smooth.stl`) is not
distributed with this repository. Please open an issue or contact the
ExtremeFLOW team directly to request a copy.

## Generating the case files

The script `create_examples.sh` generates case files and (optionally) meshes for
all combinations of resolution and flow parameters. It must be run from the
**project root**:

```bash
./examples/static_mixer/create_examples.sh        # generate case files only
./examples/static_mixer/create_examples.sh mesh   # also generate meshes and seed caches
```

The parameter space iterated by the script is:

| Variable | Values                 | Description                                     |
| -------- | ---------------------- | ----------------------------------------------- |
| `N`      | `64`, `128`            | Number of elements in the longest (x) direction |
| `Re`     | `1000`, `1500`, `3000` | Reynolds number                                 |
| `Pe`     | mirrors `Re`           | Peclet number                                   |

For each combination the script:

1. Computes `Nx = N`, `Ny = Nz = N/4` and derives the PDE filter radius as
   `r = 4.0 / Nx` (one element width in x).
2. Sets the Brinkman upper limit to `k_1 = 1e6 / Re`.
3. Copies `petsc/case.template` and substitutes all derived values (mesh path,
   `Re`, `Pe`, filter radius, Brinkman limits, cache path).
4. Writes a LUMI-G job script to
   `scripts/jobscripts/LUMI-G/static_mixer/petsc/`.
5. When the `mesh` argument is given, calls `mesh.sh` to generate the
   background mesh and runs a short warm-up simulation to populate the
   Brinkman distance cache. If the cache is not found at the end of this step
   the generated case and job files are removed automatically.

Generated case files are written to `examples/static_mixer/petsc/` with names
of the form `<N>_re_<Re>_pe_<Pe>.case`.

## Running the example

Use the top-level `run.sh` from the project root. The script locates the case
file under `examples/`, sets up the environment and launches Neko. A Neko
installation must be available; set the `NEKO_DIR` environment variable if it is
not in the default location (`external/neko`).

```bash
# Run a single case
./run.sh static_mixer/petsc/32_re_1000_pe_1000.case

# Run all static mixer cases
./run.sh static_mixer/petsc

# Run with a specific number of MPI ranks
./run.sh --procs 16 static_mixer/petsc/128_re_1000_pe_1000.case

# Submit to a cluster (e.g. LUMI)
./run.sh --submit LUMI-G static_mixer/petsc/128_re_1000_pe_1000.case
```

Useful `run.sh` flags:

| Flag                              | Description                                                    |
| --------------------------------- | -------------------------------------------------------------- |
| `-p N` / `--procs N`              | Number of MPI ranks                                            |
| `-r` / `--re-run`                 | Re-run even if results already exist                           |
| `-c` / `--clean`                  | Remove logs from previous runs                                 |
| `-d` / `--delete`                 | Delete result directories before running                       |
| `--dry-run`                       | Print the commands that would be executed without running them |
| `-s CLUSTER` / `--submit CLUSTER` | Submit job scripts to a cluster queue                          |

## Expected output and visualisation

During the simulation Neko writes two types of VTKHDF output files (by default
into a `fields/` subdirectory next to the case file):

| File                    | Contents                                                                    |
| ----------------------- | --------------------------------------------------------------------------- |
| `fields/field_0.vtkhdf` | Time-series of velocity (`Velocity`) and temperature (`temperature`) fields |
| `brinkman_0.vtkhdf`     | Static Brinkman indicator field \f$ \chi \f$ (written once at startup)      |
| Checkpoint files        | Restart snapshots written every `checkpoint_value` simulation time units    |

The `visualise.py` script post-processes the output with ParaView (`pvbatch`).
It reads both VTKHDF files, clips the domain at the mid-plane, shows cross-
sections at the inlet (`x ~= 0`) and outlet (`x = 4`), and renders two separate
views - one coloured by velocity magnitude, one by temperature. A PNG frame is
saved for each timestep.

```bash
# Single-process post-processing
pvbatch examples/static_mixer/visualise.py <path/to/fields/parent>

# Parallel post-processing (use MPICH-compatible launcher)
mpiexec.mpich -n 8 pvbatch examples/static_mixer/visualise.py <path/to/fields/parent>
```

Available options:

| Option                      | Default         | Description                                                       |
| --------------------------- | --------------- | ----------------------------------------------------------------- |
| `--stride N`                | 1               | Process every Nth timestep                                        |
| `--brinkman-clip VALUE`     | 0.5             | Threshold on \f$ \chi \f$ used to surface the solid geometry      |
| `--text-color black\|white` | black           | Color of annotations; each variant is saved in its own sub-folder |
| `--output-dir DIR`          | `visualisation` | Output directory (relative to `input_dir`)                        |
| `--overwrite`               | false           | Re-render frames that already exist                               |

Output PNG files are written to `<input_dir>/visualisation/<text-color>/` with
names `Velocity_XXXXXX.png` and `Temperature_XXXXXX.png` where `XXXXXX` is the
zero-padded timestep index. A well-mixed case at steady state should show a
nearly uniform temperature field at the outlet cross-section, confirming that
the mixer structure has homogenised the two inlet streams.
