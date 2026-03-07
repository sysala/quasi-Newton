# Paper Experiments for `magnetic_potential_problem_2D`

Source: [magnetic_problem.m](magnetic_problem.m) and case wrappers in [scripts/](scripts/).

This note summarizes the magnetic benchmark setup and the replicated
results in this repository from local checkpoints.

## Problem Setup (2D scalar magnetic potential)

- Unknown: scalar potential (u) in 2D
- Constitutive map: S = a(|nabla u|) * grad u
- Scalar conductivity law on the heterogeneous subregion:
  - `a(r) = alpha + (1-alpha) * r^8 / (r^8 + beta)`
- Global source loading strength is set by `rho`.
- Outer methods: `newton`, `newton_quasi1`, `newton_quasi2`.
- Mesh family in this benchmark: `P1` and `P2` are supported; all paper-mimicking
  sweeps in this folder use `P1`.
- Linear inner solver policy used in replication outputs:
  - `DCG_HYPRE_BOOMERAMG` (`DCG + HYPRE`), with scalar BoomerAMG setup
    (`num_functions = 1`, no elasticity near-nullspace).

Note: the parameter set below is the benchmark default used by this
repository workflow; no paper source table is checked into this snapshot.

## Replication Tables from Tolerance Sweep (Full)

- Tolerances used for sweep: `2e-01`, `1e-01`, `1e-02`, `1e-03`, `1e-04`, `1e-05`
- Each table cell: `nonlinear_iterations / cumulative_linear_iterations / runtime [s]`
- A `~` prefix marks a nonconverged/non-fully-validated outer solve entry.

### Inclusion contrast δ=1

- case id: `case_delta_1`
- base geometry: `1 x 1`
- material parameters: `alpha=0.0003`, `beta=16000`, `rho=40`, `kappa=25`, `delta=1`

#### Level 3 mesh

| linear tolerance | Newton | qN1 | qN2 |
|---|---|---|---|
| `2e-01` | `~100 / 44 / 0.848` | `~100 / 13 / 0.481` | `~100 / 99000 / 27.5` |
| `1e-01` | `~100 / 50 / 0.849` | `~100 / 23 / 0.607` | `~100 / 99000 / 28.2` |
| `1e-02` | `~100 / 115 / 0.952` | `~100 / 74 / 0.697` | `~100 / 99000 / 28.5` |
| `1e-03` | `~100 / 104 / 0.898` | `~100 / 38 / 0.58` | `~100 / 99000 / 28.2` |
| `1e-04` | `~100 / 182 / 1.2` | `~100 / 1120 / 1.25` | `~100 / 99000 / 28.2` |
| `1e-05` | `~100 / 1238 / 1.84` | `~100 / 3209 / 2.27` | `~100 / 99000 / 27.9` |

#### Level 4 mesh

| linear tolerance | Newton | qN1 | qN2 |
|---|---|---|---|
| `2e-01` | `~100 / 20 / 1.92` | `~100 / 25 / 1.89` | `~100 / 99000 / 62.2` |
| `1e-01` | `~13 / 0 / 0.833` | `~56 / 0 / 3.29` | `~100 / 0 / 5.88` |
| `1e-02` | `~100 / 1086 / 3.69` | `~100 / 75 / 2.87` | `~100 / 99000 / 64` |
| `1e-03` | `~100 / 101 / 2.86` | `~100 / 41 / 2.11` | `~100 / 99000 / 63.9` |
| `1e-04` | `~13 / 0 / 0.791` | `~56 / 0 / 3.28` | `~100 / 0 / 5.89` |
| `1e-05` | `~100 / 323 / 7.22` | `~100 / 235 / 5.76` | `~100 / 99000 / 63.9` |

### Inclusion contrast δ=2

- case id: `case_delta_2`
- base geometry: `1 x 1`
- material parameters: `alpha=0.0003`, `beta=16000`, `rho=40`, `kappa=25`, `delta=2`

#### Level 3 mesh

| linear tolerance | Newton | qN1 | qN2 |
|---|---|---|---|
| `2e-01` | `~100 / 46 / 0.735` | `-` | `-` |
| `1e-01` | `~100 / 22 / 0.593` | `-` | `-` |
| `1e-02` | `~100 / 83 / 0.749` | `-` | `-` |
| `1e-03` | `~100 / 166 / 1.04` | `-` | `-` |
| `1e-04` | `~100 / 147 / 1.05` | `-` | `-` |
| `1e-05` | `~100 / 222 / 1.29` | `-` | `-` |

#### Level 4 mesh

| linear tolerance | Newton | qN1 | qN2 |
|---|---|---|---|
| `2e-01` | `~100 / 24 / 1.97` | `-` | `-` |
| `1e-01` | `~100 / 42 / 2.7` | `-` | `-` |
| `1e-02` | `~100 / 90 / 2.75` | `-` | `-` |
| `1e-03` | `~100 / 1167 / 4.39` | `-` | `-` |
| `1e-04` | `~100 / 226 / 6.19` | `-` | `-` |
| `1e-05` | `~100 / 321 / 7` | `-` | `-` |

## Selected-Tolerance Replications (DCG + HYPRE)

These tables use the policy: `N` at `1e-4`, `qN1/qN2` at `1e-1`.
Only `DCG + HYPRE` is tracked here to match the selected replication cut.
Cell format remains `nonlinear_iterations / cumulative_linear_iterations / runtime [s]`.
A `~` prefix means the reported nonlinear solve did not satisfy the stored convergence criterion.

### Inclusion contrast δ=1

#### Level 4 mesh

- mesh level: `4`
- measured elements: `51200`
- measured unknowns: `25281`

| method family | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |
|---|---|---|---|
| DCG | `~13 / 0 / 0.843` | `~56 / 0 / 3.82` | `~100 / 0 / 5.95` |

### Inclusion contrast δ=2

#### Level 4 mesh

- mesh level: `4`
- measured elements: `51200`
- measured unknowns: `25281`

| method family | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |
|---|---|---|---|
| DCG | `~100 / 226 / 6.9` | `-` | `-` |

## Mapping to This Repository

- Benchmark entry point: [magnetic_problem.m](magnetic_problem.m)
- Benchmark wrapper: [run_paper_case.m](scripts/run_paper_case.m)
- Sweep controller: [run_paper_tolerance_sweep.m](scripts/run_paper_tolerance_sweep.m)
- Selected replication runner:
  - [run_selected_tolerance_mesh_extension.m](scripts/run_selected_tolerance_mesh_extension.m)
  - [run_selected_tolerance_mesh_method.m](scripts/run_selected_tolerance_mesh_method.m)
- Replication checkpoints:
  - [paper_tolerance_sweep_results.mat](results/paper_tolerance_sweep_results.mat)
  - [selected_tolerance_mesh_results_dcg_hypre.mat](results/selected_tolerance_mesh_results_dcg_hypre.mat)
- Plot generator: [plot_paper_tolerance_sweep.m](scripts/plot_paper_tolerance_sweep.m)
- Current report source: [generate_paper_experiments_md.m](scripts/generate_paper_experiments_md.m)
