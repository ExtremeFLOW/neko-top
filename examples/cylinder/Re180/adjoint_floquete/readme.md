# Cylinder Re180 adjoint Floquet analysis

This example computes adjoint Floquet modes for the Re180 cylinder periodic orbit.

The case starts from `baseflow0.f00000`, which is the phase-zero periodic orbit
field produced by the Newton periodic-orbit example. The Newton driver writes
that file with its field time set to the converged period `T`, so the Floquet
driver reads the integration time from the baseflow file itself. If the field
time is not positive, the driver falls back to `case.lightkrylov.floquet.period`
and then `case.time.end_time`.

The driver uses the checkpointed transient-baseflow adjoint propagator over one
full period, writes the spectrum to `adjoint_floquet_spectrum.csv`, and writes
adjoint mode movies to `viz/`.

Run `prepare.sh` if you need to regenerate `newton.nmsh`.
