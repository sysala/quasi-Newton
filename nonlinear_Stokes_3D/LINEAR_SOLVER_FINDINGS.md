# Linear Solver Port Findings

This note records the `slope_stability` linear-solver port into `nonlinear_Stokes_3D`, the new runtime API, and the test results from this environment.

## What Was Ported

- `LINEAR_SOLVERS.DFGMRES` and its dependencies:
  - `dfgmres_solver.m`
  - `IterationCollector.m`
  - `A_orthogonalize.m`
  - `A_orthogonalize_to.m`
- `LINEAR_SOLVERS.DCG` as the object-based API companion to `DFGMRES`
  - `dcg_solver.m`
- Preconditioner factory and builders:
  - `set_linear_solver.m`
  - `diag_prec_AGMG.m`
  - `diag_prec_ICHOL.m`
  - `diag_prec_HYPRE_BOOMERAMG.m`
  - `check_agmg_present.m`
- HYPRE/BoomerAMG support:
  - `hypre_boomeramg_setup.m`
  - `hypre_boomeramg_apply.m`
  - `hypre_boomeramg_clear.m`
  - `build_hypre_boomeramg_mex.m`
  - `near_null_space_elasticity_3D.m`
  - `mex/hypre_boomeramg_mex.cpp`

## New Nonlinear Stokes API

`nonlinear_Stokes_3D.m` now accepts script-level overrides before execution.

Relevant knobs:

- `solver_type`
- `linear_solver_tolerance`
- `linear_solver_maxit`
- `deflation_basis_tolerance`
- `linear_solver_printing`
- `boomeramg_opts`
- `density`
- `run_postprocess`

Supported solver types:

- `DFGMRES_AGMG`
- `DFGMRES_ICHOL`
- `DFGMRES_HYPRE_BOOMERAMG`
- `DCG_AGMG`
- `DCG_ICHOL`
- `DCG_HYPRE_BOOMERAMG`
- `DIRECT`

Example:

```matlab
density = 4;
run_postprocess = false;
solver_type = "DFGMRES_AGMG";
nonlinear_Stokes_3D
```

## HYPRE Notes

BoomerAMG uses the same near-null-space construction pattern as `slope_stability`.

The project now also contains a one-call setup script in `scripts/`, copied
from the `slope_stability` workflow and adapted for `nonlinear_Stokes_3D`:

```bash
bash nonlinear_Stokes_3D/scripts/setup_hypre_mex.sh --jobs 8
```

One important detail on this machine: MATLAB MEX linking did not work with the
default static `libHYPRE.a` build. The setup script therefore configures HYPRE
with:

- `-DCMAKE_POSITION_INDEPENDENT_CODE=ON`
- `-DBUILD_SHARED_LIBS=ON`

This produces `third_party/hypre-openmp/lib/libHYPRE.so`, which the MEX file
links against successfully.

The port also adds a custom root override:

```matlab
boomeramg_opts = struct( ...
    'hypre_root', '/path/to/hypre-openmp', ...
    'threads', 16, ...
    'print_level', 0, ...
    'use_as_preconditioner', true);
solver_type = "DFGMRES_HYPRE_BOOMERAMG";
```

If `hypre_boomeramg_mex` is not already built, the solver factory tries to build it from:

- default: `third_party/hypre-openmp`
- override: `boomeramg_opts.hypre_root`

On this MATLAB installation, `mex` reports a post-build ENOTMEX warning even
when the shared object is valid. `build_hypre_boomeramg_mex.m` already handles
that case by loading the produced file after link.

## Test Matrix

Environment used:

- MATLAB `R2024b Update 6`
- reduced smoke-test mesh: `density = 4`
- `run_postprocess = false`

### `DCG_AGMG`

Command:

```bash
matlab -batch "cd('nonlinear_Stokes_3D'); density=4; run_postprocess=false; solver_type='DCG_AGMG'; nonlinear_Stokes_3D"
```

Result:

- Newton: `8` iterations, `59` cumulative linear iterations
- Quasi-Newton 1: `63` iterations, `240` cumulative linear iterations
- Quasi-Newton 2: `147` iterations, `203` cumulative linear iterations
- Status: passed

### `DFGMRES_AGMG`

Command:

```bash
matlab -batch "cd('nonlinear_Stokes_3D'); density=4; run_postprocess=false; solver_type='DFGMRES_AGMG'; nonlinear_Stokes_3D"
```

Result:

- Newton: `8` iterations, `57` cumulative linear iterations
- Quasi-Newton 1: `63` iterations, `238` cumulative linear iterations
- Quasi-Newton 2: `147` iterations, `201` cumulative linear iterations
- Status: passed

### `DCG_HYPRE_BOOMERAMG`

Command:

```bash
matlab -batch "cd('nonlinear_Stokes_3D'); density=4; run_postprocess=false; solver_type='DCG_HYPRE_BOOMERAMG'; nonlinear_Stokes_3D"
```

Result:

- Newton: `8` iterations, `33` cumulative linear iterations
- Quasi-Newton 1: `63` iterations, `155` cumulative linear iterations
- Quasi-Newton 2: `147` iterations, `140` cumulative linear iterations
- Status: passed

### `DFGMRES_HYPRE_BOOMERAMG`

Command:

```bash
matlab -batch "cd('nonlinear_Stokes_3D'); density=4; run_postprocess=false; solver_type='DFGMRES_HYPRE_BOOMERAMG'; nonlinear_Stokes_3D"
```

Result:

- Newton: `8` iterations, `33` cumulative linear iterations
- Quasi-Newton 1: `63` iterations, `158` cumulative linear iterations
- Quasi-Newton 2: `147` iterations, `138` cumulative linear iterations
- Status: passed

## Conclusion

- The object-based `DCG` and `DFGMRES` solver API is now available in `nonlinear_Stokes_3D`.
- AGMG-backed runs work for both `DCG` and `DFGMRES`.
- HYPRE/BoomerAMG code is ported and wired, including near-null-space support.
- The copied `setup_hypre_mex.sh` workflow now builds a working local HYPRE install and MEX binding for this project.
