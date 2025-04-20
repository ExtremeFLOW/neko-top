# `neko-top` and `LightKrylov`
This is a step in the direction of demonstrating how the modulararity of both 
`neko-top` and [`LightKrylov`](https://github.com/nekStab/LightKrylov) can be 
utilized to perform linear stability analysis, attempting to replicate some of the
functionality of either the
[KTH_Framework](https://github.com/KTH-Nek5000/KTH_Framework)
or
[nekStab](https://github.com/nekStab/nekStab) projects with the addition of
GPU acceleration.

This example simply compiles and runs the `ginzburg_landau` example from the
`LightKrylov` repository using the CMake structure of `neko-top`.
