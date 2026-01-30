# Adjoint sensitivity analysis {#adjoint}

\tableofcontents

## Adjoint fluid
In `neko`, and by extension `neko-top`, we solve the Navier--Stokes equations

\f[
    {\frac {\partial \mathbf {u} }{\partial t}}
    + (\mathbf {u} \cdot \nabla)\mathbf {u}
    =
    -\nabla p
    +{\frac {1}{Re}}\nabla ^{2}\mathbf {u}
    +  \mathbf{f}, \text{ in } \Omega,\\
    \nabla \cdot \mathbf {u} = 0,  \text{ in } \Omega, 
\f]
where \f$\mathbf {u}(\mathbf{x},t)\f$ denotes the velocity field,
\f$p(\mathbf{x},t)\f$
the pressure field, \f$\mathbf{f}\f$ a forcing term and where \f$Re\f$ denotes
the Reynolds number. 
The governing equations are subjected  boundary conditions on
\f[
    \frac{1}{Re} \nabla \mathbf{u} \cdot \mathbf{n} -p \mathbf{n} = 0, 
    \text{ on } \Gamma_O, \\
    \mathbf{u} = \mathbf{u}_\text{in}, \text{ on } \Gamma_D, \\
\f]
where \f$\Gamma_N\f$ and \f$\Gamma_D\f$ denote outflow and dirichlet boundaries
respectively.

A common formulation of the adjoint Navier-Stokes
equations reads

\f[
    {-\frac {\partial \mathbf {u}^\dagger }{\partial t}}
    + (\nabla \mathbf {u})^T \mathbf {u}^\dagger
    - (\mathbf {u} \cdot \nabla) \mathbf {u}^\dagger
    =
    -\nabla p^\dagger
    +{\frac {1}{Re}}\nabla ^{2}\mathbf {u}^\dagger
    +  \mathbf{f}^\dagger, \\
    \nabla \cdot \mathbf {u} ^\dagger= 0,
\f]

where
where \f$\mathbf {u}^\dagger(\mathbf{x},t)\f$ denotes the adjoint velocity field,
\f$p^\dagger(\mathbf{x},t)\f$ the adjoint pressure field and \f$\mathbf{f}^\dagger\f$
denoting a forcing term applied to the adjoint system which generally arises
as a consequence of objective functions being evaluated.

The corresponding boundary conditions generally read

\f[
    \frac{1}{Re} \nabla \mathbf{u}^\dagger \cdot \mathbf{n} 
    -p^\dagger \mathbf{n}
    + (\mathbf{u} \cdot \mathbf{n}) \mathbf{u}^\dagger = \mathbf{0}
    \text{ on } \Gamma_N, \\
    \mathbf{u} = \mathbf{0}, \text{ on } \Gamma_D. \\
\f]

The spectral element method which underpins `neko` solves
the above system of equations using the weak formulation, which has important
implications when solving the adjoint system. Primarily, when deriving the
adjoint system in strong form, one introduces additional boundary terms on
Neumann boundaries, and more importantly, one applies the divergence free
condition strongly.

An alternative to the above adjoint system is to remain in weak form and not
integrate the convective term by parts, resulting in no additional boundary
terms, and more importantly, no pointwise application of the divergence free
condition.

This alternative system is what is implemented in `neko-top`, which reads (in
weak form),

\f[
    -\int_\Omega \mathbf{v}\cdot {\frac {\partial \mathbf {u}^\dagger }
    {\partial t}}
    + \int_\Omega \mathbf{v}\cdot (\nabla \mathbf {u})^T \mathbf {u}^\dagger
    + \int_\Omega \nabla \mathbf{v}\cdot  (\mathbf {u} \otimes \mathbf {u}^\dagger )
    =
    -\int_\Omega \mathbf{v}\cdot \nabla p^\dagger
    +{\frac {1}{Re}}\int_\Omega \nabla \mathbf{v}\cdot \nabla \mathbf {u}^\dagger
    + \int_\Omega \mathbf{v}\cdot  \mathbf{f}^\dagger, \\
    \int_\Omega q  \nabla \cdot \mathbf {u} ^\dagger= 0,
\f]
where \f$\mathbf {v}\f$ is a test function.

## Adjoint scalar
The adjoint scalar equation takes the following form
\f[
    {\frac {\partial \phi }{\partial t}}
    + (\mathbf {u} \cdot \nabla)\phi
    =
    {\frac {1}{Pe}}\nabla ^{2}\phi
    , \text{ in } \Omega,\\
\f]
where \f$\phi(\mathbf{x},t)\f$ denotes the scalar field,
and where \f$Pe\f$ denotes the Peclet number. In the context of conjugate
heat transfer, the velocity equation is often coupled to the scalar equation
through the Boussinesq approximation for instance, however, in lieu of this
coupling the scalar is often referred to as a "passive scalar" to imply the
one way coupling.

In `neko-top` the adjoint scalar reads (in weak form),
\f[
    -\int_\Omega  \psi\cdot {\frac {\partial \phi^\dagger }
    {\partial t}}
    + \int_\Omega (\nabla \psi \cdot  \mathbf {u})  \phi^\dagger 
    =
    {\frac {1}{Pe}}\int_\Omega \nabla \psi \cdot \nabla \phi^\dagger,\\
\f]
where \f$\phi(\mathbf{x},t)^\dagger\f$ denotes the adjoint scalar field and 
where \f$\psi\f$ is a test function. In addition, the perturbation of the term
\f$(\mathbf {u} \cdot \nabla)\phi\f$ results in an additional term in adjoint
velocity momentum equation, which now reads
\f[
    -\int_\Omega \mathbf{v}\cdot {\frac {\partial \mathbf {u}^\dagger }
    {\partial t}}
    + \int_\Omega \mathbf{v}\cdot (\nabla \mathbf {u})^T \mathbf {u}^\dagger
    + \int_\Omega \nabla \mathbf{v}\cdot  (\mathbf {u} \otimes \mathbf {u}^\dagger )
    \underline{+ \int_\Omega (\mathbf{v}\cdot \nabla \phi ) \phi^\dagger}
    =
    -\int_\Omega \mathbf{v}\cdot \nabla p^\dagger
    +{\frac {1}{Re}}\int_\Omega \nabla \mathbf{v}\cdot \nabla \mathbf {u}^\dagger
    + \int_\Omega \mathbf{v}\cdot  \mathbf{f}^\dagger.
\f]

It can be seen from the above equations that due to the one way coupling between
\f$\mathbf{u}\f$ and \f$\phi\f$, when solving the forward problem one must first
solve for \f$\mathbf{u}\f$ and then solve for \f$\phi\f$. However, in the
adjoint the opposite occurs as the equation for \f$\mathbf{u}^\dagger\f$ depends
on \f$\phi^\dagger\f$, hence in `neko-top` we first solve for the adjoint scalar
and then solve for the adjoint velocity.

### Multiple scalars
In the case of using multiple scalar transport equations involving multiple
scalar fields \f$\phi_i(\mathbf{x},t)\f$, the assumption is made that
all scalars are passive, that is, coupled in one way to the velocity, and,
in no way coupled to one another. If this is not the case, additional care
must be taken by the user when deriving their adjoint equations, and will likely
introduce additional coupling terms in the adjoint equations. In lieu of this,
each adjoint equation will take an identical form to that stated above, and
the adjoint velocity momentum equation will contain the summation over all the
additional coupling terms

\f[
    - \int_\Omega \mathbf{v} \cdot \frac{\partial \mathbf{u}^\dagger}{\partial t}
    + \int_\Omega \mathbf{v} \cdot (\nabla \mathbf{u})^T \mathbf{u}^\dagger
    + \int_\Omega \nabla \mathbf{v} \cdot (\mathbf{u} \otimes \mathbf{u}^\dagger)
    + \underline{\sum_i \int_\Omega (\mathbf{v} \cdot \nabla \phi_i) \phi_i^\dagger}
    =
    - \int_\Omega \mathbf{v} \cdot \nabla p^\dagger
    + \frac{1}{Re} \int_\Omega \nabla \mathbf{v} \cdot \nabla \mathbf{u}^\dagger
    + \int_\Omega \mathbf{v} \cdot \mathbf{f}^\dagger.
\f]

\note If a given \f$\phi_i(\mathbf{x},t)\f$ is not contained in the
objective, there will be no resulting forcing to drive the adjoint scalar
equation, effectively yielding \f$\phi_i(\mathbf{x},t)^\dagger = \mathbf{0}\f$,
and hence, will in theory have no contribution to the adjoint momentum equation.

## Immersed Boundary Methods
Following a Brinkman style immersed boundary method, the presence of an
immersed object is imposed by the Brinkman forcing term \f$\mathbf{f} = - \chi
    \mathbf{u}\f$, 
    where \f$\chi\f$ is the spatially dependent Brinkman coefficient, satisfying 
\f[
    \chi =
    \begin{cases} 
    0 & \text{in the fluid region,} \\
    \overline{\chi} & \text{in the solid region.}
    \end{cases}
\f]
Considering \f$\overline{\chi}\f$ to be a large value, this discontinuous forcing 
term models a momentum loss in solid region, and thereby simulating porous 
media with very low permeability. It is worth noting that while the Brinkman 
penalization method is rooted in the idea of modelling solid regions as porous
media with vanishing permeability, the  approach proposed by 
[Goldstein](https://doi.org/10.1006/jcph.1993.1081) frames the interaction of
the solid on the fluid as a control problem, where the feedback force is tuned 
to drive the velocity to zero in the solid region.
Regardless of their different motivation, both methods result in a similar 
mathematical structure. More information regarding the mapping of \f$\chi\f$ can
be found in [Mapping cascade](@ref mapping_cascade).

\note In the future we will provide a full adjoint derivation in this section
of the theory guide. This will tie together all aspects from the mapping, to
how these terms arise etc. For now we are simply documenting the equations
being solved.

## Time integration
This section summarizes the discrete adjoint for the pressure-correction
time integration used in `Neko-TOP`. We start from a right-handed Riemann sum
for the objective,
\f[
  F \approx \Delta t \sum_{n=1}^{N} f(\mathbf{u}_n), \qquad
  \mathrm{d}F = \Delta t \sum_{n=1}^{N}
  \int_\Omega \nabla_{\mathbf{u}_n} f(\mathbf{u}_n)\cdot \delta \mathbf{u}_n \, \mathrm{d}\Omega.
\f]

The velocity-pressure splitting currently implemented `Neko` follows that of
[Karniadakis et al. 1991](https://doi.org/10.1016/0021-9991(91)90007-8) and is commonly
referred to as the  \f$ \mathbb{P}_N-\mathbb{P}_N\f$ formulation, whereby
we introduce an intermediate velocity \f$\mathbf{u}_n^{\star}\f$ and write the
forward time step as three sub-problems:
\f[
  \int_\Omega \mathbf{v}_n \cdot \frac{\mathbf{u}_n^{\star}}{\Delta t} \, \mathrm{d}\Omega
  =
  \int_\Omega \mathbf{v}_n \cdot \frac{\mathbf{u}_{n-1}}{\Delta t} \, \mathrm{d}\Omega
  - \int_\Omega \mathbf{v}_n \cdot \mathcal{C}_{\mathbf{u}_{n-1}}\,\mathbf{u}_{n-1} \, \mathrm{d}\Omega
  - \int_\Omega \mathbf{v}_n \cdot \chi\,\mathbf{u}_{n-1} \, \mathrm{d}\Omega,
\f]
\f[
  \int_\Omega \nabla q_n \cdot \nabla p_n \, \mathrm{d}\Omega
  =
  \int_\Omega \nabla q_n \cdot \frac{\mathbf{u}_n^{\star}}{\Delta t} \, \mathrm{d}\Omega
  - \int_\Omega \nabla q_n \cdot \nabla \times (\nabla \times \mathbf{u}_{n-1}) \, \mathrm{d}\Omega,
\f]
\f[
  \int_\Omega \mathbf{w}_n \cdot \frac{\mathbf{u}_n}{\Delta t} \, \mathrm{d}\Omega
  + \int_\Omega \mathbf{w}_n \cdot \mathcal{L}\mathbf{u}_n \, \mathrm{d}\Omega
  =
  \int_\Omega \mathbf{w}_n \cdot \frac{\mathbf{u}_n^{\star}}{\Delta t} \, \mathrm{d}\Omega
  - \int_\Omega \mathbf{w}_n \cdot \nabla p_n \, \mathrm{d}\Omega.
\f]
Here \f$\mathcal{C}_{\mathbf{u}_{n-1}}\f$ denotes the linearized convection
operator about \f$\mathbf{u}_{n-1}\f$, \f$\chi\f$ denotes any explicit linear
operator (e.g. Brinkman), and \f$\mathcal{L}\f$ is the implicit
diffusion/Helmholtz operator.

Perturbing these equations and collecting terms yields a discrete adjoint that
is solved in three steps (backward in time). First, solve for the adjoint
velocity associated with the implicit step,
\f[
  -\int_\Omega \frac{\mathbf{u}_n^{*\dagger}}{\Delta t} \, \mathrm{d}\Omega
  + \int_\Omega \mathcal{L}^\dagger \mathbf{u}_n^{*\dagger} \, \mathrm{d}\Omega
  =
  \int_\Omega \nabla_{\mathbf{u}_n} f(\mathbf{u}_n)\,\Delta t \, \mathrm{d}\Omega
  + \int_\Omega \frac{\mathbf{u}_{n+1}^{\dagger}}{\Delta t} \, \mathrm{d}\Omega
  - \int_\Omega \mathcal{C}^\dagger_{\mathbf{u}_n}\,\mathbf{u}_{n+1}^{\dagger} \, \mathrm{d}\Omega
  - \int_\Omega \chi^\dagger\, \mathbf{u}_{n+1}^{\dagger} \, \mathrm{d}\Omega.
\f]
Second, solve the adjoint pressure equation,
\f[
  \int_\Omega \nabla q_n \cdot \nabla p_n^{\dagger} \, \mathrm{d}\Omega
  =
  \int_\Omega -\mathbf{u}_n^{*\dagger} \cdot \nabla p_n^{\dagger} \, \mathrm{d}\Omega.
\f]
Finally, recover the adjoint velocity associated with the explicit step,
\f[
  \mathbf{u}_n^{\dagger}
  = \mathbf{u}_n^{*\dagger} - \nabla p_n^{\dagger}.
\f]

With these adjoint variables, the discrete gradient contribution associated
with an explicit linear operator \f$\chi\f$ is
\f[
  \mathrm{d}F(\chi,\delta \chi)
  = -\sum_{n=1}^{N} \int_\Omega \delta \chi\,
  \mathbf{u}_{n+1}^{\dagger} \mathbf{u}_n \, \mathrm{d}\Omega.
\f]
The order above matches the implementation: solve the implicit adjoint step,
then the adjoint pressure correction, then form the projection to obtain the
adjoint velocity used in the sensitivity accumulation.

\warning The curl--curl contribution in the pressure equation is implemented
as the boundary integral
\f$-\int_{\partial\Omega}\big(\mathbf{n}\times\nabla p_{n+1}^{\dagger}\big)\cdot
\big(\nabla\times \mathbf{w}\big)\,\mathrm{d}S\f$, where \f$\mathbf{w}\f$
is the test function of adjoint velocity equation, however
it has not yet been tested extensively.
