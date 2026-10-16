# Cylinder stability examples

This directory ports the cylinder examples from
`https://github.com/nekStab/neklab/tree/main/examples/cylinder` to the
Neko/Neko-TOP/LightKrylov workflow used in this repository.

The shared Nek mesh is stored in `_shared/1cyl.re2`. Run

```bash
./prepare.sh
```

to generate `_shared/1cyl.nmsh` before running the JSON case files. The
converter may need the same HDF5/JSON-Fortran runtime libraries used to build
Neko.

The LightKrylov drivers write visualization-ready `.fld` files under their
local `viz/` directories at algorithmic milestones, not through case-file
output controls. Fixed-point Newton writes the converged baseflow; stability
runs write baseflows and modes; the periodic-orbit driver writes a sampled
movie over one converged period; the resolvent driver writes baseflow, forcing,
and response real/imaginary parts.

Implemented here:

- `dns/Re40` and `dns/Re180`: Neko JSON translations of the upstream DNS runs.
- `Re40/baseflow`, `Re40/LNS`, `Re40/adjoint`, and `Re40/transient_growth`:
  steady-baseflow stability examples.
- `Re60/baseflow`, `Re60/LNS`, and `Re60/adjoint`: artificially stabilized
  steady-baseflow examples.
- `Re180/periodic_orbit`, `Re180/LNS_floquete`, and
  `Re180/adjoint_floquete`: periodic-orbit Floquet examples.
- `resolvent`: complex shifted resolvent solve using the time-stepper
  approximation to the continuous linearized generator.
