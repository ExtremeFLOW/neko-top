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
- `newton/Re40_fixed_point`: fixed-point Newton solve from the upstream
  `BF.fld` initial guess.
- `stability/direct` and `stability/adjoint`: fixed-point eigenspectrum runs
  using the local `sources/stability` LightKrylov/Neko propagator layer.
- `newton/Re180_periodic_orbit`: extended-state Newton solve with the period as
  an unknown.
- `resolvent`: complex shifted resolvent solve using the time-stepper
  approximation to the continuous linearized generator.
