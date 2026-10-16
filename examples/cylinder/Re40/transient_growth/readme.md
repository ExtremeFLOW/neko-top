# Cylinder transient growth

This mirrors neklab's transient-growth example: compute the leading singular
triplets of the finite-time linearized propagator about a fixed baseflow.

The baseflow is not bundled here. Generate it first with the Re40 Newton fixed
point example, then run transient growth:

```sh
cd ../baseflow
./prepare.sh
./neko 1cyl_newton.case

cd ../transient_growth
./prepare.sh
./neko 1cyl_transient_growth.case
```

The transient-growth case reads
`../baseflow/viz/baseflow_newton_Re400.f00000`.

Outputs:

- `singular_spectrum.dat`: singular values and residuals.
- `svds_output.txt`: LightKrylov intermediate convergence history.
- `viz/baseflow_transient_growth*.fld`: fixed baseflow.
- `viz/prt*.fld`: optimal initial perturbations.
- `viz/rsp*.fld`: optimal responses at the time horizon.
