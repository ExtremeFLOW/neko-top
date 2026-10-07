# Direct fixed-point stability

Runnable translation of neklab's `cylinder/stability/direct` example using the
local Neko-TOP/LightKrylov propagator. The driver sets
`adjoint_case%if_adjoint = false`, keeps the loaded state as the fixed
baseflow, and then calls `eigs` on the tangent propagator over `tau = 1`.
