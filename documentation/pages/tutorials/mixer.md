# Static Mixer {#mixer}

\tableofcontents

## Optimization problem

This tutorial aims to replicate the passive mixer problem by 
[C. S. Andreasen et al. 2009](https://doi.org/10.1002/fld.1964).

The goal of this optimization problem is to design a structure capable of
mixing fluids in a static mannor, that is to say, without the use of moving
parts.

We begin by considering a domain \f$\Omega\f$ which is a square duct seen in
the figure below.

PUT A BLANK FIG

Fluid enters the domain upstream and is ejected downstream, with all other
boundary being considered solid walls. 
In addition, a scalar field, representative of the two species one desires to
have mixed, enters upstream, separated into two distinct species. The goal of
the optimization problem is to design the interior structure such that the
two species are as homogeneous as possible in the downstream location.

In this tutorial you will learn how to
- Formally define this optimization problem in the context of a topology
optimization problem.
- Mesh the domain \f$\Omega\f$
- Set up a fluid, scalar and their respective adjoint simulations
- Define the objective functions being solved

Their objective function was essentially

min \f$\int_{\Gamma_{out}} (\phi - \bar{\phi})^2 d\Gamma\f$ 
(with some normalization)

s.t. \f$\Delta P \leq \beta \Delta P_{ref}\f$

They set \f$P_{out} = 0\f$ implying \f$\Delta P = \int_{\Gamma_{in}}p d \Gamma\f$