# CEEC LHC-3 Example

This example holds the case files etc for the CEEC LHC-3 benchmark problem.
The problem is a series of mixer problems, with varying settings.

## Things to implement

Notebook to ensure we can monitor the run and its progress.

- Statistics of design change, average change, max change. (Will probably come
  with the stagnation check.)
- Statistics for sensitivities, norm of the sensitivity, maximum sensitivity.
  Unconstrained: KKT is norm and max of the gradient, so should be golden.
- Objective function history, weighted history plot of total and individuals,
  and the unweighted history plot of individuals.
- MND: measure of non-discreteness.
