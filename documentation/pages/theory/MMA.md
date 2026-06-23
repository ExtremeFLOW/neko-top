# The Method of Moving Asymptotes (MMA) {#MMA}

\tableofcontents

## Overview

The Method of Moving Asymptotes (MMA) is a gradient-based optimization algorithm
widely used in topology optimization. It is particularly effective for
large-scale, constrained, nonlinear optimization problems.

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

At each iteration, MMA constructs a **separable convex approximation** of the original nonlinear problem by replacing each function \f$ f_i(x) \f$ with an **asymptotic approximation** built from (note that \f$f_0\f$ is also approximated in the same way):

- the **current function value** \f$ f_i(x^k) \f$
- the **first-order sensitivities** \f$ \nabla f_i(x^k) \f$

---

### Asymptotic Approximation

Instead of a standard Taylor expantion, MMA uses a **asymptotic model** of the form:

\f[
f_i(x) \;\approx\; \tilde{f}_i(x)
= \sum_{j=1}^n \left(
\frac{p_{ij}}{u_j - x_j} + \frac{q_{ij}}{x_j - l_j}
\right) - b_i
\f]

where:
- \f$ l_j \f$, \f$ u_j \f$ are the **moving asymptotes**
- \f$ p_{ij}, q_{ij} \f$ are **positive coefficients constructed from the gradients**
- \f$ b_i \f$ is chosen such that the approximation is **exact at the current iteration**: \f$ \tilde{f}_i(x^k) = f_i(x^k) \f$

---
### Construction from Function Value and Gradient

The coefficients \f$ p_{ij} \f$, \f$ q_{ij} \f$ are derived from the sensitivities \f$ \frac{\partial f_i}{\partial x_j} \f$ to ensure:

- **First-order consistency**:
  \f[
  \nabla \tilde{f}_i(x^k) \approx \nabla f_i(x^k)
  \f]

- **Convexity** (by enforcing \f$ p_{ij}, q_{ij} > 0 \f$)

In practice (as implemented in `mma_gensub`):

\f[
\begin{aligned}
p_{ij} &\sim \max\!\left(\frac{\partial f_i}{\partial x_j}, 0\right) \cdot (u_j - x_j)^2 \\
q_{ij} &\sim \max\!\left(-\frac{\partial f_i}{\partial x_j}, 0\right) \cdot (x_j - l_j)^2
\end{aligned}
\f]

with small regularization terms added for numerical stability.
Thus, the exact implementation is:

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
  - more conservative bounds (`"alpha"`, `"beta"`) for the updated design variables are choosen to prevent them from hitting bounds caused by too aggressive updates.

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
- `asydecr` (contraction factor to pull the asymtotes together to take more conservative steps)

---

## Subproblem Solution

The convex subproblem is solved using a **primal-dual interior point method** `"pdip"`  and a pure **dual interior point method** `"dip"`.

In this implementation both subsolvers support both **CPU and device (GPU)** execution

---

## KKT Convergence

Convergence is evaluated using Karush-Kuhn-Tucker (KKT) conditions:

- Infinity norm: `residumax`
- Euclidean norm: `residunorm`

These provide a quantitative measure of optimality.

---

## Implementation Details

The implementation is encapsulated in the `mma_t` type and includes the following parameters that can be set in the case file:

| Name | Description | Default |
|------|-------------|---------|
| `mma.max_iter` | Max iterations for subproblem | `100` |
| `mma.epsimin` | KKT tolerance scaling | \f$ 10^{-9} \sqrt{m + n} \f$ |
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
| `mma.auto_scale` | If `true`, sensitivity and function values for the constriant \f$ f_i \f$ are scaled at each iteration with a differnt value such that we get \f$ f_1(x_k)= \f$`mma.scale`. This would be an adaptive scaling based on the value of the first constriant. | `false` |

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
+ \frac{1}{2} z^2
\end{aligned}
\f]


so the quadratic terms for \f$y_i\f$ and \f$z\f$ are enforced to make sure that we can solve the minimization problems, analytically.


### Practical implication

Compared to a standard primal MMA subsolve:
- the quadratic terms for \f$y_i\f$ and \f$z\f$ are considered in the solver regardless of what parameters user set in the case file.
- the primal subproblem is **not solved directly**
- instead, each iteration computes:
  - analytical minimizers in \f$x_j, y_i, z\f$
  - and then performs a **dual ascent on** \f$\lambda\f$
- DIP is cheaper per iteration

---

## References

- Krister Svanberg, *The Method of Moving Asymptotes*
- Lazarov & Sigmund (2011)
