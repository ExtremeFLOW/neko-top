# Resolvent

This is the Neko/Neko-TOP translation of neklab's `cylinder/resolvent` setup.
The included `BF_1cyl0.f00001` baseflow is copied from upstream.

The driver uses `complex_state_vector_t` and solves
`(i omega I - L) x = f` with LightKrylov GMRES. The continuous-time generator is
approximated from the existing time-stepper as `L v ~= (Phi_tau v - v) / tau`.
