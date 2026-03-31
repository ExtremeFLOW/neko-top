# Code structure {#code-structure}
\tableofcontents

This page provides an overview of the code structure of the Neko-TOP library. 
The structure is designed to be modular and extensible, with a clear separation
of concerns between the different components of the library. While supporting
the control of the optimization loop. The library is designed to be
flexible and allow for easy extension and modification.


## Components
`neko-top` distinguishes between four key objects:
- \subpage design. The design space over which one optimizes.
- \subpage problem. The definition of the optimization problem being solved.
- \subpage optimizer. The optimization algorithm used to solve the optimization
  problem.
- \subpage simulation. Specifically for problems involving fluid mechanics, a
  simulation can allow for interfaces with `neko` to perform forward simulation
  and additional `neko-top` libraries to perform adjoint sensitivity analysis.
- \subpage continuation. The continuation strategy for gradually changing
  different parameters in the run time using a continuation scheduler.
