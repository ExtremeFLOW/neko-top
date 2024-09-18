# Simple setup to compile the current example.
set(DRIVER_TYPE "user")
set(EXTRA_SOURCES adjoint.f90 user.f90 adjoint_scalar_scheme.f90 adjoint_scalar_pnpn.f90)

build_example(${DRIVER_TYPE})
