# Re60 adjoint

Runnable translation of the artificially stabilized Re60 adjoint example using the
local Neko-TOP/LightKrylov propagator. The driver sets
`adjoint_case%if_adjoint = true`, keeps the loaded state as the fixed baseflow,
and then calls `eigs` on the adjoint tangent propagator over `tau = 1`.
