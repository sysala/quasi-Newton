# nonlinear_Stokes_3D Structure

This MATLAB benchmark is organized into package folders by responsibility.

## Entry Point

- `nonlinear_Stokes_3D.m`: benchmark driver, configuration, assembly, solver calls, and postprocessing.
- `agmg/`: optional external AGMG dependency. If present locally, the entry script adds it to the MATLAB path automatically.

## Packages

- `+ASSEMBLY`: quadrature, basis functions, sparse assembly, and load-vector construction.
- `+CONSTITUTIVE_PROBLEM`: constitutive operators for Newton and quasi-Newton variants.
- `+LINEAR_SOLVERS`: Krylov methods, preconditioners, and linear-solver helpers.
- `+MESH`: mesh generation and boundary-condition setup tied to geometry.
- `+NEWTON`: nonlinear outer solvers.
- `+VIZ`: plotting and postprocessing utilities.

## Placement Guidelines

- Put element-level integration and matrix assembly in `+ASSEMBLY`.
- Put stress or tangent-law logic in `+CONSTITUTIVE_PROBLEM`.
- Put reusable `Ax=b` solvers and preconditioners in `+LINEAR_SOLVERS`.
- Put geometry, surfaces, and boundary masks in `+MESH`.
- Put Newton or quasi-Newton iteration logic in `+NEWTON`.
- Put rendering and plotting helpers in `+VIZ`.

## Run

```bash
matlab -batch "cd('nonlinear_Stokes_3D'); nonlinear_Stokes_3D"
```
