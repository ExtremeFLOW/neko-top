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

- Read the case file in question and plot the interpolations used (RAMP, BP,
  etc).
- Continued parameters should be plotted as a function of design iterations.
- We need to read through both the current run but also older runs (run_XX)
  since the current might not be the only run.


Forslag til ny notebook eller job-monitor python script baseret på egne
erfaringer med at jeg tror koden gør A, mens den i virkeligheden gør B. Dette
kan (for det mestes) kommes til livs, med et lille script man kan køre efter
jobstart der tager 1) input filen / setup filen og 2) stdout som input. Scriptet
tjekker så de to mod hinanden og plotter f.eks. interpolationer, continuation
schemes, historikker on-the-fly. Når et job netop er startet er det et super
godt sanity check, og når jobbet kører er det samme script super til at tjekke
at alt går som det skal og alle continuations har den ønskede effekt.

Jarvis er sygt god til at lave den slags
Niels
