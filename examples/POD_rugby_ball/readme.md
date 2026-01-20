# Unsteady rugby ball {#unsteady-rugby-ball}

This example is an unsteady extension of the classic "rugby ball" example by
[Borrvall & Petersson 2003](https://doi.org/10.1002/fld.426).

The objective is to minimize the time integral of the dissipation

\f[
\mathcal{F} = \int_0^T \frac{1}{|\Omega_\text{obj}|}\int \frac{1}{2} 
\left[
\nabla \mathbf{u} \cdot \left(\nabla \mathbf{u} + (\nabla \mathbf{u})^T \right)
+ \chi \mathbf{u} \cdot \mathbf{u} 
\right] d\Omega, dt
\f]

where the velocity boundary condition progressively increase from 0 to 1 over
the course of 0.5 time units.

It is also subject to a volume constraint

\f[
\mathcal{C} = \frac{1}{|\Omega_\text{opt}|}\int_{\Omega_\text{opt}} \rho d\Omega.
\f]


More information regarding objectives and constraints can be found in
[Objectives and constraints](@ref objectives_and_constraints).
