# nonlinear_Stokes_3D

This MATLAB benchmark is organized into package folders for the core solver and
plain directories for reports, scripts, figures, and cached experiment outputs.

## Layout

- `nonlinear_Stokes_3D.m`: main benchmark entry point.
- `paper.tex`: paper source used as the reference for the extracted experiment notes.
- `PAPER_EXPERIMENTS.md`: compact markdown extraction of the paper setup and replication tables.
- `PAPER_TOLERANCE_SWEEP.md`: tolerance-sweep report.
- `LINEAR_SOLVER_FINDINGS.md`: notes from the solver/preconditioner port.
- `images/`: generated tracked figures used by the markdown reports.
- `results/`: local `.mat` checkpoints and caches. This directory keeps a local `.gitignore`, so result files stay out of git.
- `scripts/`: experiment drivers, report generators, and build helpers.
- `agmg/`: optional local AGMG dependency. If present, the main script adds it to the MATLAB path automatically.

## Packages

- `+ASSEMBLY`: quadrature, basis functions, sparse assembly, and load-vector construction.
- `+CONSTITUTIVE_PROBLEM`: constitutive operators for Newton and quasi-Newton variants.
- `+LINEAR_SOLVERS`: Krylov methods, preconditioners, and solver helpers.
- `+MESH`: mesh generation and boundary-condition setup.
- `+NEWTON`: nonlinear outer solvers.
- `+VIZ`: plotting and postprocessing utilities.

## What The Scripts Produce

- `scripts/run_paper_tolerance_sweep.m`
  Produces `results/paper_tolerance_sweep_results.mat`, refreshes `PAPER_TOLERANCE_SWEEP.md`, and regenerates the PNG plots in `images/`.
- `scripts/generate_paper_experiments_md.m`
  Rebuilds `PAPER_EXPERIMENTS.md` from the saved sweep checkpoints.
- `scripts/run_selected_tolerance_mesh_extension.m`
  Runs the selected-tolerance extra mesh cases and writes `results/selected_tolerance_mesh_results_dcg_hypre.mat`.
- `scripts/run_selected_tolerance_mesh_method.m`
  Runs or marks one selected-mesh nonlinear-method case. This is mainly a helper for resumable automation.
- `scripts/run_selected_tolerance_mesh_hardcap.sh`
  Drives the `~200k` selected-mesh runs with one MATLAB process per nonlinear method and an outer hard timeout.
- `scripts/setup_hypre_mex.sh`
  Clones/builds HYPRE and then builds the BoomerAMG MATLAB MEX binding.

## How To Run

Run the benchmark itself:

```bash
matlab -batch "cd('nonlinear_Stokes_3D'); nonlinear_Stokes_3D"
```

Run the tolerance sweep and regenerate the report:

```bash
matlab -batch "cd('nonlinear_Stokes_3D'); addpath('scripts'); run_paper_tolerance_sweep"
```

Regenerate the paper-experiments markdown from existing checkpoints only:

```bash
matlab -batch "cd('nonlinear_Stokes_3D'); addpath('scripts'); generate_paper_experiments_md"
```

Run the selected `DCG + HYPRE` `~200k` mesh extension with hard time caps:

```bash
bash nonlinear_Stokes_3D/scripts/run_selected_tolerance_mesh_hardcap.sh
```

Build HYPRE and the MATLAB MEX binding:

```bash
bash nonlinear_Stokes_3D/scripts/setup_hypre_mex.sh --jobs 8
```

## Placement Guidelines

- Put element-level integration and sparse assembly in `+ASSEMBLY`.
- Put material-law and tangent logic in `+CONSTITUTIVE_PROBLEM`.
- Put reusable `Ax = b` solvers and preconditioners in `+LINEAR_SOLVERS`.
- Put geometry, surfaces, and boundary masks in `+MESH`.
- Put outer nonlinear iteration logic in `+NEWTON`.
- Put plotting or report visualization in `+VIZ`.
- Put reproducibility drivers and generators in `scripts/`.
- Put generated figures in `images/`.
- Put local cached `.mat` outputs in `results/`.
