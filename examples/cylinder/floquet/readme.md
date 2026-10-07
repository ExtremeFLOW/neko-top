# Cylinder Re180 Floquet analysis

This example computes Floquet multipliers for the Re180 cylinder periodic orbit.

The case starts from `periodic_orbit_Re1800.f00000`, which is the phase-zero
periodic orbit field produced by the Newton periodic-orbit example. The driver
uses the transient-baseflow LNS propagator over one full period, writes the
Floquet spectrum to `floquet_spectrum.csv`, and writes mode movies to `viz/`.

Run `prepare.sh` if you need to regenerate `newton.nmsh`.
