# POD Rugby Ball {#pod-rugby-ball}

This example is the POD state-recovery variant of
[Unsteady rugby ball](@ref unsteady-rugby-ball). It uses the same setup and
optimization problem, but replaces checkpoint-based state recovery with POD.

The objective is to minimize the time integral of the dissipation

\f[
\mathcal{F} = \int_0^T \frac{1}{|\Omega_\text{obj}|}\int \frac{1}{2}
\left[
\nabla \mathbf{u} \cdot \left(\nabla \mathbf{u} + (\nabla \mathbf{u})^T \right)
+ \chi \mathbf{u} \cdot \mathbf{u}
\right] d\Omega, dt
\f]
downstream while subject to a volume constraint

\f[
\mathcal{C} =
\frac{1}{|\Omega_\text{opt}|}\int_{\Omega_\text{opt}} \rho d\Omega.
\f]

More information regarding objectives and constraints can be found in
[Objectives and constraints](@ref objectives_and_constraints).
