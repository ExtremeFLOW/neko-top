# Re180 periodic orbit

This is the Neko/Neko-TOP translation of neklab's
`cylinder/newton/Re180_periodic_orbit` setup.

The driver uses `ext_state_vector_t`, containing the flow state and the unknown
period, and `periodic_orbit_system_t`, whose residual is `Phi_T(q) - q` plus one
phase condition.

After Newton converges, the driver advances the converged orbit over one period
and writes `viz/periodic_orbit_Re180_movie*.fld` frames. The number of frames is
controlled by `case.lightkrylov.periodic_orbit.movie_frames`.

The initial Newton guess is whatever `case.fluid.initial_condition` defines.
The driver first marches that state for
`case.lightkrylov.periodic_orbit.transient_time`, then uses the final field as
the Newton initial condition. Set `transient_time` to `0.0` to start Newton
directly from the case initial condition.

The Newton iterates are written to `viz/outer_periodic_orbit_Re180*.fld`, using
the orbit state at `t=0`, so the convergence of the initial condition can be
inspected directly.
