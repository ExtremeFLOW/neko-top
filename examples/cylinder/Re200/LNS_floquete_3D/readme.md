# Cylinder Re200 3D Floquet analysis

This example computes Floquet multipliers for 3D perturbations of the Re200
cylinder periodic orbit.

The source slab field is `baseflow0.f00000`, which is repeated in the spanwise
direction by `extrude_baseflow.py` to create `baseflow3d0.f00000`. The generated
3D field uses `u=u(x,y)`, `v=v(x,y)`, `w=0`, and spans length 4 with 4 elements.
The case interpolates this field onto the generated 3D mesh when initializing.
The driver reads the integration time from this field. If the field time is not
positive, it falls back to `case.lightkrylov.floquet.period` and then
`case.time.end_time`.

The driver uses the transient-baseflow LNS propagator over one full period,
writes the Floquet spectrum to `floquet_spectrum.csv`, and writes mode movies
to `viz/`.

Run `prepare.sh` to regenerate `newton.nmsh` and `baseflow3d0.f00000`.
