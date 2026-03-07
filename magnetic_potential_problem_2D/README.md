# magnetic_potential_problem_2D

This MATLAB benchmark is organized into package folders for the core solver and
plain directories for scripts, figures, and cached experiment outputs.

## Layout

- `magnetic_problem.m`: main benchmark entry point.
- `PAPER_EXPERIMENTS.md`: markdown extraction of the benchmark setup and local
  replication tables.
- `PAPER_TOLERANCE_SWEEP.md`: generated tolerance sweep report.
- `images/`: generated figures used by report scripts.
- `results/`: local `.mat` checkpoints and caches. The folder keeps a local
  `.gitignore` to avoid committing generated data files.
- `scripts/`: experiment drivers and report generators.
- `+ASSEMBLY`: quadrature, basis functions, and assembly helpers.
- `+CONSTITUTIVE_PROBLEM`: constitutive operators and quasi-Newton operators.
- `+LINEAR_SOLVERS`: Krylov/iterative solvers and preconditioners, including HYPRE
  BoomerAMG.
- `+MESH`: mesh generation routines.
- `+NEWTON`: nonlinear and quasi-Newton outer solvers.
- `+VIZ`: lightweight plotting and visualization helpers.

## What the Scripts Produce

- `scripts/run_paper_case.m`
  Runs one fixed case with a unified workspace contract and stores a compact
  result record.
- `scripts/run_paper_tolerance_sweep.m`
  Performs sweep over linear-tolerance values and writes:
  - `results/paper_tolerance_sweep_results.mat`
  - `PAPER_TOLERANCE_SWEEP.md`
  - runtime plots in `images/`.
- `scripts/run_selected_tolerance_mesh_extension.m`
  Runs the selected-mesh method-replication policy for:
  - `newton` at `1e-4`
  - `qn1`, `qn2` at `1e-1`
  using `DCG + HYPRE` on the selected mesh level and writes:
  - `results/selected_tolerance_mesh_results_dcg_hypre.mat`.
- `scripts/run_selected_tolerance_mesh_method.m`
  Runs or marks one selected-mesh nonlinear method case (`newton`, `qn1`, `qn2`) to support
  resumable per-method orchestration.
- `scripts/run_damping_level4_sweep.m`
  Runs the level-4 damping sweep over `δ`, `ρ`, and `κ` with fixed outer tolerances
  (`N = 1e-4`, `qN1/qN2 = 1e-1`) and writes `DAMPING_LEVEL4_EXPERIMENTS.md`.
- `scripts/generate_paper_experiments_md.m`
  Rebuilds `PAPER_EXPERIMENTS.md` from checkpoint data.
- `scripts/plot_paper_tolerance_sweep.m`
  Regenerates tolerance-vs-runtime plots.
- `scripts/setup_hypre_mex.sh`
  Builds/clones HYPRE and builds the BoomerAMG MEX extension for this folder.

## How To Run

Run the benchmark:

```bash
matlab -batch "cd('magnetic_potential_problem_2D'); magnetic_problem"
```

Run the tolerance sweep and regenerate all markdown reports and plots:

```bash
matlab -batch "cd('magnetic_potential_problem_2D'); addpath('scripts'); run_paper_tolerance_sweep"
```

Regenerate the paper-style markdown note only:

```bash
matlab -batch "cd('magnetic_potential_problem_2D'); addpath('scripts'); generate_paper_experiments_md"
```

Run the selected-mesh replication case by case (method by method):

```bash
matlab -batch "cd('magnetic_potential_problem_2D'); addpath('scripts'); run_selected_tolerance_mesh_method('case_delta_1', 'newton')"
```

To run the damped variant for one method, pass `use_damping`:

```bash
matlab -batch "cd('magnetic_potential_problem_2D'); addpath('scripts'); run_selected_tolerance_mesh_method('case_delta_1', 'newton', struct('use_damping', true))"
```

Build HYPRE binaries for this folder:

```bash
bash magnetic_potential_problem_2D/scripts/setup_hypre_mex.sh
```

## HYPRE default in replication scripts

The magnetic benchmark is scalar in 2D; `magnetic_potential_problem_2D/scripts/run_paper_case.m`
and the sweep scripts pass `num_functions = 1` to BoomerAMG, i.e. no 3D
elasticity nullspace basis is injected.
