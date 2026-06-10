# The Method of Moving Asymptotes (MMA) {#MMA}

\tableofcontents

The Method of Moving Asymptotes (MMA) is the optimization algorithm currently
used in `neko-top`. In a case file, MMA is configured in
`optimization.solver`, with MMA-specific settings under
`optimization.solver.mma`.

## Options overview

| Json entry | Description | Admissible Values | Default value |
|------------|-------------|-------------------|---------------|
| `optimization.solver.type` | Select optimizer implementation. | `"mma"` | `"mma"` |
| `optimization.solver.max_iterations` | Maximum number of outer optimization iterations. | Integer, `> 0` | `100` |
| `optimization.solver.tolerance` | Convergence tolerance used by the optimizer loop. | Real, `> 0` | `1.0e-3` |
| `optimization.solver.enable_output` | Enable writing optimization logs/CSV output. | `.true.` or `.false.` | `.true.` |
| `optimization.solver.mma.max_iter` | Maximum number of iterations in the MMA subsolver. | Integer, `> 0` | `100` |
| `optimization.solver.mma.epsimin` | Minimum barrier parameter used by MMA. | Real, `> 0` | `1.0e-9 * sqrt(m + n_global)` |
| `optimization.solver.mma.asyinit` | Initial asymptote distance factor. | Real, `> 0` | `0.2` |
| `optimization.solver.mma.asyincr` | Asymptote increase factor. | Real, `> 0` | `1.05` |
| `optimization.solver.mma.asydecr` | Asymptote decrease factor. | Real, `> 0` | `0.65` |
| `optimization.solver.mma.backend` | MMA execution backend. | `"cpu"` or `"device"` | `"cpu"` (or `"device"` when `NEKO_BCKND_DEVICE=1`) |
| `optimization.solver.mma.subsolver` | Interior-point variant used in MMA subproblems. | `"dip"` or `"dpip"` | `"dip"` |
| `optimization.solver.mma.xmin` | Global lower bound applied to all design variables. | Real | `0.0` |
| `optimization.solver.mma.xmax` | Global upper bound applied to all design variables. | Real | `1.0` |
| `optimization.solver.mma.a0` | MMA scalar coefficient \f$a_0\f$. | Real | `1.0` |
| `optimization.solver.mma.a` | MMA vector coefficient \f$a_i\f$ (applied uniformly). | Real | `0.0` |
| `optimization.solver.mma.c` | MMA vector coefficient \f$c_i\f$ (applied uniformly). | Real | `100.0` |
| `optimization.solver.mma.d` | MMA vector coefficient \f$d_i\f$ (applied uniformly). | Real | `0.0` |
| `optimization.solver.mma.scale` | Constraint scaling target value. | Real, `> 0` | `10.0` |
| `optimization.solver.mma.auto_scale` | Enable adaptive scaling of constraints each iteration. | `.true.` or `.false.` | `.false.` |

## Description of related options

### Optimizer-loop options
- `optimization.solver.type` selects the optimizer. MMA is currently the
  implemented choice.
- `optimization.solver.max_iterations` and `optimization.solver.tolerance`
  control stopping of the outer optimization loop.
- `optimization.solver.enable_output` controls whether optimizer progress is
  written to output files.

### MMA subproblem and asymptote options
- `optimization.solver.mma.max_iter` and `optimization.solver.mma.epsimin`
  define limits for the MMA subsolver iterations and barrier parameter.
- `optimization.solver.mma.asyinit`, `optimization.solver.mma.asyincr`, and
  `optimization.solver.mma.asydecr` control moving-asymptote updates.
- `optimization.solver.mma.subsolver` selects the interior-point variant used
  inside MMA (`dip`/`dpip`).

### Design bounds and MMA coefficients
- `optimization.solver.mma.xmin` and `optimization.solver.mma.xmax` set global
  lower/upper bounds for all design variables.
- `optimization.solver.mma.a0`, `optimization.solver.mma.a`,
  `optimization.solver.mma.c`, and `optimization.solver.mma.d` are coefficients
  in the standard MMA subproblem formulation.

### Backend and scaling options
- `optimization.solver.mma.backend` selects whether MMA update/KKT operations
  run on CPU or device routines.
- `optimization.solver.mma.scale` sets the constraint scaling magnitude.
- `optimization.solver.mma.auto_scale` enables adaptive scaling from the
  current constraint value at each iteration.
