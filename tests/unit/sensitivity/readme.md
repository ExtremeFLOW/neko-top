# Finite difference Sensitivity Test

This test verifies the sensitivity analysis of a design variable using finite
difference methods. It checks that the computed sensitivities match the
expected values within a specified tolerance.

The test is tagged as a `unit` test, which means it is mandatory for the
build to pass our CI/CD pipeline.

## Test Overview

The test is done through a few common files:

- `prepare.sh`: Script designed to construct a mesh for us to work on.
- `sensitivity.f90`: Contain the finite difference sensitivity subroutine.
- `sensitivity.case`: The Neko case file to control the simulation.
- `CMakeLists.txt`: The CMake file defining the build process and adding the
  sensitivity test as a CTest.

Currently, we have added tests for the following components:

- `volume_constraint_t`