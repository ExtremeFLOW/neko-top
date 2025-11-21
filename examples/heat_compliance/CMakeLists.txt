# Simple setup to compile the current example.
set(DRIVER_TYPE "custom")
set(DRIVER driver.f90)
set(EXTRA_SOURCES heat_compliance.f90)

build_example(${DRIVER_TYPE})