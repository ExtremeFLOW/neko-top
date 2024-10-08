# Code structure
\tableofcontents

This page provides an overview of the code structure of the Neko-TOP library. 
The structure is designed to be modular and extensible, with a clear separation
of concerns between the different components of the library. While supporting
the control flow of the optimization loop, the library is designed to be
flexible and allow for easy extension and modification.

## Overall structure of the control flow

- Initialization
- Optimization loop
  - Solve Fluid
    - Map design to IBM
    - Run Neko
    - Update objective
    - Update constraints
  - Update the design
    - Compute adjoint solution
    - Compute the sensitivity
    - Run MMA optimizer.
  - Check convergence

## List of tasks

- [ ] Build dummy simulation components for the optimization loop.
  - [ ] Immersed boundary method.
  - [ ] Optimization objective and constraints.
  - [ ] Checkpoints for use in adjoint solver.
- [ ] Setup adjoint solver.
- [ ] Tie the adjoint into the sensitivity computation.
- [ ] Setup MMA optimizer.

## Immersed Boundary Method
