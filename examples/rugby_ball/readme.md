# Rugby ball {#rugby-ball}

This example follows the classic "rugby ball" example by
[Borrvall & Petersson 2003](https://doi.org/10.1002/fld.426).

The objective is to minimize the dissipation

\f[
\mathcal{F} = \frac{1}{|\Omega_\text{obj}|}\int \frac{1}{2} 
\left[
\nabla \mathbf{u} \cdot \left(\nabla \mathbf{u} + (\nabla \mathbf{u})^T \right)
+ \chi \mathbf{u} \cdot \mathbf{u} 
\right] d\Omega,
\f]

subject to a volume constraint

\f[
\mathcal{C} = \frac{1}{|\Omega_\text{opt}|}\int_{\Omega_\text{opt}} \rho d\Omega.
\f]

More information regarding objectives and constraints can be found in
[Objectives and constraints](@ref objectives_and_constraints).

The following depicts the optimization history for various volume constraints

![Design field.](/documentation/images/rugby_ball.gif)

![Convergence history.](/documentation/images/rugby_ball_convergence_history.png)
