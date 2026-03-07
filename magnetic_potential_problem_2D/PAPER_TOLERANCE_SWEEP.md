# Paper Tolerance Sweep for `magnetic_potential_problem_2D`

This report tracks Newton and quasi-Newton performance with DCG + HYPRE on the magnetic benchmark.

- tolerances: `2e-01`, `1e-01`, `1e-02`, `1e-03`, `1e-04`, `1e-05`
- completed records: `24`
- current sweep plan size: `2`
- generated from:
  - [run_paper_case.m](scripts/run_paper_case.m)
  - [run_paper_tolerance_sweep.m](scripts/run_paper_tolerance_sweep.m)
  - [magnetic_problem.m](magnetic_problem.m)
  - [paper_tolerance_sweep_results.mat](results/paper_tolerance_sweep_results.mat)

## Inclusion contrast δ=1

- case id: `case_delta_1`
- geometry: `1 x 1`
- delta: `1`, alpha: `0.0003`, beta: `16000`, rho: `40`, kappa: `25`

### Level 3 mesh

| linear tolerance | Newton | qN1 | qN2 |
|---|---|---|---|
| `2e-01` | `~100 / 44 / 0.848` | `~100 / 13 / 0.481` | `~100 / 99000 / 27.5` |
| `1e-01` | `~100 / 50 / 0.849` | `~100 / 23 / 0.607` | `~100 / 99000 / 28.2` |
| `1e-02` | `~100 / 115 / 0.952` | `~100 / 74 / 0.697` | `~100 / 99000 / 28.5` |
| `1e-03` | `~100 / 104 / 0.898` | `~100 / 38 / 0.58` | `~100 / 99000 / 28.2` |
| `1e-04` | `~100 / 182 / 1.2` | `~100 / 1120 / 1.25` | `~100 / 99000 / 28.2` |
| `1e-05` | `~100 / 1238 / 1.84` | `~100 / 3209 / 2.27` | `~100 / 99000 / 27.9` |

### Level 4 mesh

| linear tolerance | Newton | qN1 | qN2 |
|---|---|---|---|
| `2e-01` | `~100 / 20 / 1.92` | `~100 / 25 / 1.89` | `~100 / 99000 / 62.2` |
| `1e-01` | `~13 / 0 / 0.833` | `~56 / 0 / 3.29` | `~100 / 0 / 5.88` |
| `1e-02` | `~100 / 1086 / 3.69` | `~100 / 75 / 2.87` | `~100 / 99000 / 64` |
| `1e-03` | `~100 / 101 / 2.86` | `~100 / 41 / 2.11` | `~100 / 99000 / 63.9` |
| `1e-04` | `~13 / 0 / 0.791` | `~56 / 0 / 3.28` | `~100 / 0 / 5.89` |
| `1e-05` | `~100 / 323 / 7.22` | `~100 / 235 / 5.76` | `~100 / 99000 / 63.9` |

## Inclusion contrast δ=2

- case id: `case_delta_2`
- geometry: `1 x 1`
- delta: `2`, alpha: `0.0003`, beta: `16000`, rho: `40`, kappa: `25`

### Level 3 mesh

| linear tolerance | Newton | qN1 | qN2 |
|---|---|---|---|
| `2e-01` | `~100 / 46 / 0.735` | `-` | `-` |
| `1e-01` | `~100 / 22 / 0.593` | `-` | `-` |
| `1e-02` | `~100 / 83 / 0.749` | `-` | `-` |
| `1e-03` | `~100 / 166 / 1.04` | `-` | `-` |
| `1e-04` | `~100 / 147 / 1.05` | `-` | `-` |
| `1e-05` | `~100 / 222 / 1.29` | `-` | `-` |

### Level 4 mesh

| linear tolerance | Newton | qN1 | qN2 |
|---|---|---|---|
| `2e-01` | `~100 / 24 / 1.97` | `-` | `-` |
| `1e-01` | `~100 / 42 / 2.7` | `-` | `-` |
| `1e-02` | `~100 / 90 / 2.75` | `-` | `-` |
| `1e-03` | `~100 / 1167 / 4.39` | `-` | `-` |
| `1e-04` | `~100 / 226 / 6.19` | `-` | `-` |
| `1e-05` | `~100 / 321 / 7` | `-` | `-` |

Cells are formatted as `nonlinear_iterations / cumulative_linear_iterations / runtime [s]`.
A row prefixed with `~` did not satisfy the stored convergence criterion when stopping.
