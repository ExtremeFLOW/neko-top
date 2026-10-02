# The Method of Moving Asymptotes (MMA) {#MMA}

\tableofcontents

## Overview

The Method of Moving Asymptotes (MMA) is a gradient-based optimization algorithm
widely used in topology optimization. It is particularly effective for
large-scale, constrained, non-linear optimization problems.

The method transforms the original non-convex optimization problem into a
sequence of strictly convex subproblems that are easier to solve.

## Original Optimization Problem

The MMA implementation in Neko-TOP solves problems of the form

\f[
\begin{aligned}
\min_{x,z,y} \quad & f_0(x) + a_0 z + \sum_{i=1}^m \left( c_i y_i + \frac{1}{2} d_i y_i^2 \right) \\
\text{s.t.} \quad & f_i(x) - a_i z - y_i \le 0, \quad i = 1, \dots, m, \\
& x_j^{\min} \le x_j \le x_j^{\max}, \quad j = 1, \dots, n, \\
& z \ge 0, \quad y_i \ge 0.
\end{aligned}
\f]

Here:
- \f$x\f$: design variables
- \f$f_0\f$: objective function
- \f$f_i\f$: constraint functions
- \f$y_i\f$, \f$z\f$: auxiliary variables for constraint relaxation

---
## Convex Approximation

At each iteration, MMA constructs a **separable convex approximation** of the original non-linear problem by replacing each function \f$ f_i(x) \f$ with an **asymptotic approximation** built from (note that \f$f_0\f$ is also approximated in the same way):

- the **current function value** \f$ f_i(x^k) \f$
- the **first-order sensitivities** \f$ \nabla f_i(x^k) \f$

---

### Asymptotic Approximation

Instead of a standard Taylor expansion, MMA uses a **asymptotic model** of the form:

\f[
f_i(x) \;\approx\; \tilde{f}_i(x)
= \sum_{j=1}^n \left(
\frac{p_{ij}}{u_j - x_j} + \frac{q_{ij}}{x_j - l_j}
\right) - b_i
\f]

where:
- \f$ l_j \f$, \f$ u_j \f$ are the **moving asymptotes**
- \f$ p_{ij}, q_{ij} \f$ are **non-negative coefficients constructed from the gradients**
- \f$ b_i \f$ is chosen such that the approximation is **exact at the current iteration**: \f$ \tilde{f}_i(x^k) = f_i(x^k) \f$

---
### Construction from Function Value and Gradient

The coefficients \f$ p_{ij} \f$, \f$ q_{ij} \f$ are derived from the sensitivities \f$ \frac{\partial f_i}{\partial x_j} \f$ to ensure:

- **First-order consistency**:
  \f[
  \nabla \tilde{f}_i(x^k) \approx \nabla f_i(x^k)
  \f]

- **Convexity** (by enforcing \f$ p_{ij}, q_{ij} \ge 0 \f$)

In practice (as implemented in `mma_gensub`):

\f[
\begin{aligned}
p_{ij} &\sim \max\!\left(\frac{\partial f_i}{\partial x_j}, 0\right) \cdot (u_j - x_j)^2 \\
q_{ij} &\sim \max\!\left(-\frac{\partial f_i}{\partial x_j}, 0\right) \cdot (x_j - l_j)^2
\end{aligned}
\f]

with small regularization terms added for numerical stability. The exact
form depends on the subsolver.

For the `dip` subsolver on the CPU backend it follows `MMA::GenSub` in the
topopt_in_petsc code by Aage (with `constraintModification = false`): only the objective
(\f$ i = 0 \f$) is regularized,

\f[
\begin{aligned}
p_{0j} &=
\left(
\max\!\left(\frac{\partial f_0}{\partial x_j}, 0\right)
+ 0.001 \left|\frac{\partial f_0}{\partial x_j}\right|
+ \frac{0.5 \cdot 10^{-6}}{u_j - l_j}
\right)
(u_j - x_j)^2,
&
p_{ij} &= \max\!\left(\frac{\partial f_i}{\partial x_j}, 0\right) (u_j - x_j)^2,
\\[8pt]
q_{0j} &=
\left(
\max\!\left(-\frac{\partial f_0}{\partial x_j}, 0\right)
+ 0.001 \left|\frac{\partial f_0}{\partial x_j}\right|
+ \frac{0.5 \cdot 10^{-6}}{u_j - l_j}
\right)
(x_j - l_j)^2,
&
q_{ij} &= \max\!\left(-\frac{\partial f_i}{\partial x_j}, 0\right) (x_j - l_j)^2,
\end{aligned}
\f]

for \f$ i = 1, \dots, m \f$.

For the `pdip` subsolver, and for both subsolvers on the device backend, it
follows `mmasub` by Svanberg, where the objective
and all constraints (\f$ i = 0, \dots, m \f$) are regularized:

\f[
\begin{aligned}
p_{ij} &=
\left(
1.001 \max\!\left(\frac{\partial f_i}{\partial x_j}, 0\right)
+ 0.001 \max\!\left(-\frac{\partial f_i}{\partial x_j}, 0\right)
+ \frac{10^{-5}}{\max(x_{\text{diff},j}, 10^{-5})}
\right)
(u_j - x_j)^2
\\[8pt]
q_{ij} &=
\left(
0.001 \max\!\left(\frac{\partial f_i}{\partial x_j}, 0\right)
+ 1.001 \max\!\left(-\frac{\partial f_i}{\partial x_j}, 0\right)
+ \frac{10^{-5}}{\max(x_{\text{diff},j}, 10^{-5})}
\right)
(x_j - l_j)^2
\end{aligned}
\f]

Here \f$ x_{\text{diff},j} = x_j^{\max} - x_j^{\min} \f$, restricted by the move
limit if one is set.

![An example of how the upper and lower asymptotes are used to construct a local approximation of the function](mmaGensubExample.png)


---

### Resulting Convex Subproblem

Using these approximations, MMA solves the following **convex separable problem**:

\f[
\min_x \sum_{j=1}^n \left(
\frac{p_{0j}}{u_j - x_j} + \frac{q_{0j}}{x_j - l_j}
\right)
+ a_0 z + \sum_{i=1}^m \left(c_i y_i + \frac{1}{2} d_i y_i^2\right)
\f]

subject to:

\f[
\sum_{j=1}^n \left(
\frac{p_{ij}}{u_j - x_j} + \frac{q_{ij}}{x_j - l_j}
\right)
+ a_i z + y_i \le b_i
\f]

---

### Key Properties of This Approximation

- Uses only **local information**:
  - function values \f$ f_i(x^k) \f$
  - gradients \f$ \nabla f_i(x^k) \f$

- **Separable in** \f$ x_j \f$ → efficient for large-scale problems

- **Convex by construction** → guarantees well-posed subproblems

- **Asymptotic behavior near bounds**:
  - singularities at \f$ x_j \to l_j \f$ and \f$ x_j \to u_j \f$.
  - more conservative bounds (`"alpha"`, `"beta"`) for the updated design variables are chosen to prevent them from hitting bounds caused by too aggressive updates.

---

### Intuition

Unlike a Taylor expansion, MMA builds a **curvature-aware approximation** using **moving asymptotes** where

\f[
\frac{1}{u_j - x_j}, \quad \frac{1}{x_j - l_j}
\f]

act as barrier-like functions that:

- match the local gradient,
- enforce convexity,
- and guide the design smoothly within admissible bounds.

---

## Moving Asymptotes

A key feature of MMA is the adaptive update of asymptotes:

- They **control step size and stability** for each design variable
- They **prevent oscillations**
- They **improve convergence robustness**

The update is governed by:
- `asyinit` (initial spacing)
- `asyincr` (expansion factor to push the asymptotes apart to help with faster convergence using larger steps)
- `asydecr` (contraction factor to pull the asymptotes together to take more conservative steps)

---

## Subproblem Solution

The convex subproblem is solved using a **primal-dual interior point method** `"pdip"` and a pure **dual interior point method** `"dip"`.

In this implementation both subsolvers support both **CPU and device (GPU)** execution.
The alignment of the `dip` subsolver with the topopt_in_petsc code (see the
sections below) is so far only done on the CPU backend. The device backend
still uses the approximation of `mmasub`, the term \f$ \frac{1}{2} z^2 \f$ and
the earlier convergence measure.

---

## KKT Convergence

Convergence is evaluated using Karush-Kuhn-Tucker (KKT) conditions:

- Infinity norm: `residumax`
- Euclidean norm: `residunorm`

The optimizer stops when `residumax` is below `optimization.solver.tolerance`.
The residual is evaluated at the updated design, with the multipliers of the
last subproblem.

For the `dip` subsolver on the CPU backend the residual follows
`MMA::KKTresidual` in the topopt_in_petsc code by Aage. It consists of
- the gradient of the Lagrangian,
  \f$ \partial f_0 / \partial x_j + \sum_{i} \lambda_i \partial f_i / \partial x_j
  - \xi_j + \eta_j \f$,
  where the bound multipliers \f$ \xi_j, \eta_j \f$ are estimated for design
  variables within \f$ 10^{-5} \f$ of \f$ x_j^{\min} \f$ or \f$ x_j^{\max} \f$,
- the complementarity of the bounds, \f$ \xi_j (x_j - x_j^{\min}) \f$ and
  \f$ \eta_j (x_j^{\max} - x_j) \f$,
- the complementarity of the constraints,
  \f$ \sum_{i} \lambda_i (a_i z + y_i - f_i(x)) \f$.

The true bounds \f$ x^{\min}, x^{\max} \f$ are used, not the move limited ones.
Since each subproblem is only solved down to the barrier parameter `epsimin`,
the constraint term does not drop much below \f$ m \f$ times the last barrier
level, which is the smallest power of ten above `epsimin` (between `epsimin`
and \f$ 10 \f$ `epsimin`). The tolerance should be set above that.

For the `pdip` subsolver the residual is the full KKT system of the primal-dual
method (Svanberg).

---

## Implementation Details

The implementation is encapsulated in the `mma_t` type and includes the following parameters that can be set in the case file:

| Name | Description | Default |
|------|-------------|---------|
| `mma.max_iter` | Max iterations for subproblem | `100` |
| `mma.epsimin` | Smallest barrier parameter of the subsolvers, must be positive | \f$ 10^{-9} \sqrt{m + n} \f$ |
| `mma.asyinit` | Initial asymptote distance | `0.2` |
| `mma.asyincr` | Asymptote expansion | `1.05` |
| `mma.asydecr` | Asymptote contraction | `0.65` |
| `mma.move_limit` | Move limit for updating the design variables | `0.2` |
| `mma.a0` | MMA param | `1.0` |
| `mma.a` | MMA param | `0.0` |
| `mma.c` | MMA param | `1000.0` |
| `mma.d` | MMA param | `0.0` |
| `mma.backend` | `cpu` or `device` | auto based on Neko backend |
| `mma.subsolver` | Subsolver type | `dip` |
| `mma.scale` | Scaling factor applied to constraint functions \( f_i \) and their sensitivities. This does **not** affect the objective function \( f_0 \). It is used to improve numerical conditioning when constraint magnitudes and sensitivities differ significantly from those of the objective function. | `1.0` |
| `mma.auto_scale` | If `true`, sensitivity and function values for the constraint \f$ f_i \f$ are scaled at each iteration with a different value such that we get \f$ f_1(x_k)= \f$`mma.scale`. This would be an adaptive scaling based on the value of the first constraint. | `false` |

---

## Notes

- The constant term in the approximation (see Eq. 3.5 in Svanberg) is omitted,
  as it does not affect the minimization (\f$ b_0 \f$ in the approximation of \f$ f_0 \f$).

### Note on the DIP subsolver formulation

The `dip` (Dual Interior Point) subsolver solves a **slightly different but equivalent reformulation** of the MMA subproblem.

Instead of directly minimizing the primal convex approximation in \f$(x, y, z)\f$, the method forms the **Lagrangian dual problem**:

\f[
\Psi(\lambda) =
\sum_{j=1}^{n} \min_{x_j}
\left\{ L_x(x_j, \lambda) \;\middle|\; \alpha_j \le x_j \le \beta_j \right\}
+ \min_{z \ge 0} L_z(z, \lambda)
+ \sum_{i=1}^{m} \min_{y_i \ge 0} L_y(y_i, \lambda)
\f]

and then solves:

\f[
\max_{\lambda \ge 0} \; \Psi(\lambda)
\f]


### Key difference from the primal MMA subproblem

The DIP formulation uses the Lagrangian:

\f[
\begin{aligned}
L(x,y,z,\lambda) =
&\sum_{j=1}^{n}
\left(
\frac{p_{0j} + \sum_{i=1}^{m} \lambda_i p_{ij}}{u_j - x_j}
+
\frac{q_{0j} + \sum_{i=1}^{m} \lambda_i q_{ij}}{x_j - l_j}
\right)
- \sum_{i=1}^{m} \lambda_i b_i \\
&+ \sum_{i=1}^{m}
\left[
(c_i - \lambda_i) y_i + \frac{1}{2} y_i^2
\right]
+ \left(a_0 - \sum_{i=1}^{m} \lambda_i a_i\right) z
+ \frac{1}{20} z^2
\end{aligned}
\f]


so the quadratic terms for \f$y_i\f$ and \f$z\f$ are enforced to make sure that we can solve the minimization problems, analytically.
On the CPU backend the weight \f$ 1/20 \f$ of \f$ z^2 \f$ follows the topopt_in_petsc code by Aage, which gives
\f$ z = \max\left(0, 10 \left(\sum_{i} \lambda_i a_i - a_0\right)\right) \f$.


### Practical implication

Compared to a standard primal MMA subsolve:
- the quadratic terms for \f$y_i\f$ and \f$z\f$ are considered in the solver regardless of what parameters user set in the case file.
- the primal subproblem is **not solved directly**
- instead, each iteration computes:
  - analytical minimizers in \f$x_j, y_i, z\f$
  - and then performs a **dual ascent on** \f$\lambda\f$
- DIP is cheaper per iteration

### Relation to the topopt_in_petsc implementation

The `dip` subsolver on the CPU backend follows `MMA.cc` in the topopt_in_petsc
code by Aage (subproblem, dual solver and its stopping rules, KKT residual), and
reproduces its iterates to round-off on the same input. The remaining
differences are
- the default asymptote parameters: Neko-TOP uses `0.2`, `1.05`, `0.65` for
  `asyinit`, `asyincr`, `asydecr`, while `MMA.cc` uses `0.5`, `1.2`, `0.7`,
- the initial multipliers are \f$ \lambda_i = \max(1, c_i/2) \f$ instead of
  \f$ c_i/2 \f$, which only differs for \f$ c_i < 2 \f$,
- the span \f$ x_j^{\max} - x_j^{\min} \f$ is bounded below by \f$ 10^{-5} \f$
  also for the initial asymptotes, so fixed variables do not break the
  approximation,
- variables with \f$ x_j \f$ within machine precision of \f$ \alpha_j \f$ or
  \f$ \beta_j \f$ are treated as active in the dual Hessian, where `MMA.cc` tests
  whether the unclamped minimizer lies outside \f$ [\alpha_j, \beta_j] \f$,
- the dual Newton system is solved with a pivoted LU factorization (LAPACK).

The topopt_in_petsc driver stops when the largest design change drops below
\f$ 0.01 \f$ and does not evaluate the KKT residual. In Neko-TOP the design
change criterion is `optimization.solver.stop_design_change`, next to the
KKT tolerance.

The topopt_in_petsc driver also normalizes the objective to
\f$ f_0(x^0) = 10 \f$ before every update. Neko-TOP passes the objective
unscaled, while constants such as the regularization \f$ 0.5 \cdot 10^{-6} \f$
and the default \f$ c_i = 1000 \f$ are absolute. For objectives of a very
different magnitude, the weight of the objective can be used to bring it to a
similar scale.

---

## References

- Svanberg, K. (1987). The method of moving asymptotes—a new method for structural optimization. *International Journal for Numerical Methods in Engineering*, 24(2), 359-373.
- Svanberg, K. (1993). The Method of Moving Asymptotes (MMA) with Some Extensions. *Optimization of Large Structural Systems*, 191–207.
- Svanberg, K. (2002). A class of globally convergent optimization methods based on conservative convex separable approximations. *SIAM Journal on Optimization*, 12(2), 555-573.
- Aage, N., & Lazarov, B. S. (2013). Parallel framework for topology optimization using the method of moving asymptotes. *Structural and Multidisciplinary Optimization*, 47(4), 493-505.
