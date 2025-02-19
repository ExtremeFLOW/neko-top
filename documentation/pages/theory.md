# Governing equations {#governing-equations}

\tableofcontents

## Adjoint fluid
In `neko`, and by extension `neko-top`, we solve the Navier--Stokes equations

$${\frac {\partial \mathbf {u} }{\partial t}}
+ (\mathbf {u} \cdot \nabla)\mathbf {u}
=
-\nabla p
+{\frac {1}{Re}}\nabla ^{2}\mathbf {u}
+  \mathbf{f}, $$

$$ \nabla \cdot \mathbf {u} = 0, $$
where $\mathbf {u}(\mathbf{x},t)$ denotes the velocity field, $p(\mathbf{x},t)$
the pressure field, $\mathbf{f}$ a forcing term and where $\frac{1}{Re}$ denotes
the Reynolds number.


A common formulation of the adjoint Navier--Stokes
equations reads

$${-\frac {\partial \mathbf {u}^\dagger }{\partial t}}
+ (\nabla \mathbf {u})^T \mathbf {u}^\dagger
- (\mathbf {u} \cdot \nabla) \mathbf {u}^\dagger
=
-\nabla p^\dagger
+{\frac {1}{Re}}\nabla ^{2}\mathbf {u}^\dagger
+  \mathbf{f}^\dagger, $$

$$ \nabla \cdot \mathbf {u} ^\dagger= 0, $$
where
where $\mathbf {u}^\dagger(\mathbf{x},t)$ denotes the adjoint velocity field,
$p^\dagger(\mathbf{x},t)$ the adjoint pressure field and $\mathbf{f}^\dagger$
denoting a forcing term applied to the adjoint system which generally arises
as a consequence of objective functions being evaluated.

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

$$-\int_\Omega \mathbf{v}\cdot {\frac {\partial \mathbf {u}^\dagger }
{\partial t}}
+ \int_\Omega \mathbf{v}\cdot (\nabla \mathbf {u})^T \mathbf {u}^\dagger
+ \int_\Omega \nabla \mathbf{v}\cdot  (\mathbf {u} \otimes \mathbf {u}^\dagger )
=
-\int_\Omega \mathbf{v}\cdot \nabla p^\dagger
+{\frac {1}{Re}}\int_\Omega \nabla \mathbf{v}\cdot \nabla \mathbf {u}^\dagger
+ \int_\Omega \mathbf{v}\cdot  \mathbf{f}^\dagger, $$

$$ \int_\Omega q  \nabla \cdot \mathbf {u} ^\dagger= 0. $$
