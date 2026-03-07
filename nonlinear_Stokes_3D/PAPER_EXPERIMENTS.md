# Paper Experiments for `nonlinear_Stokes_3D`

Source: [paper.tex](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/paper.tex)

This note pulls out the numerical experiment settings and the paper tables from
`paper.tex` into a compact markdown reference.

## Problem Setup Used in the Paper

- Problem type: simplified nonlinear Stokes problem in 3D.
- Simplification: the pore-pressure field is given and the divergence-free
  constraint is not enforced, so the structure is close to nonlinear
  elasticity.
- Constitutive law: Carreau law

```text
a(r) = mu_infty + (mu_0 - mu_infty) * (1 + lambda * r^2)^((p - 2) / 2)
```

- Material parameter range in theory:
  - `p in (1, 2)`
  - `lambda > 0`
  - `mu_0 > mu_infty > 0`
- Nonlinear solvers compared:
  - `N`: Newton
  - `qN1`: quasi-Newton 1, corresponding to the variable preconditioner from
    Section 3.1
  - `qN2`: quasi-Newton 2, corresponding to the fixed preconditioner from
    Section 3.2
- Inner linear solver setup in the paper:
  - deflated and preconditioned conjugate gradient
  - preconditioner choices: incomplete Cholesky (`IC`) and AGMG
- Current repository defaults for these benchmark comparisons:
  - `N`: linear solver tolerance `1e-4`
  - `qN1`: linear solver tolerance `1e-1`
  - `qN2`: linear solver tolerance `1e-1`
- Finite elements: `P1`
- Reported table format:
  - `outer nonlinear iterations / cumulative linear iterations / runtime [s]`

## Example 1: Prismatic Bar

### Geometry and Boundary Conditions

- Domain:

```text
Omega = (0, rho) x (0, rho) x (0, L)
```

- Bar aligned with the `x_3` axis.
- Shell boundary condition: zero velocity in the normal direction
  (paper labels this as slip conditions).
- Left face `x_3 = 0`: prescribed constant axial velocity
  - `u_3(x_1, x_2, 0) = u_3,0`
- Right face `x_3 = L`: homogeneous Neumann condition.
- Volume force:

```text
g = (0, 0, g_3),  g_3 < 0
```

### Parameters

- `L = 50`
- `rho = 5`
- `u_3,0 = 5`
- `mu_0 = 1`
- `mu_infty = 0.001`
- `lambda = 10`
- `p = 1.1`
- varied load values:
  - `g_3 in {-3e-3, -4e-3, -5e-3, -6e-3, -7e-3}`

### Discretization

- `P1` elements
- uniform mesh
- `103,680` elements
- `54,886` unknowns

### Table from the Paper

| `g_3` | `IC-N` | `IC-qN1` | `IC-qN2` | `AGMG-N` | `AGMG-qN1` | `AGMG-qN2` |
|---|---|---|---|---|---|---|
| `-3e-3` | `6/112/1.86` | `14/84/0.95` | `15/58/0.61` | `6/51/2.12` | `14/35/0.99` | `15/44/0.85` |
| `-4e-3` | `6/111/1.84` | `17/90/1.17` | `22/70/0.79` | `6/53/2.01` | `17/40/1.25` | `22/57/1.11` |
| `-5e-3` | `7/130/2.17` | `27/97/1.72` | `35/93/1.12` | `7/63/2.44` | `25/50/1.85` | `36/86/1.68` |
| `-6e-3` | `7/133/2.17` | `41/111/2.69` | `74/167/2.11` | `7/63/2.39` | `41/70/3.05` | `75/170/3.33` |
| `-7e-3` | `8/152/2.53` | `73/143/5.17` | `>200/-/5.53` | `8/76/2.84` | `72/108/5.63` | `>200/-/6.81` |

Paper note: for smaller `|g_3|`, `qN2` is fastest; for `|g_3| = 0.007` or
larger, `qN2` becomes the slowest because the material response is strongly
nonlinear.

### Replications on This Machine

- nonlinear tolerance policy used here:
  - `N`: `1e-4`
  - `qN1`: `1e-1`
  - `qN2`: `1e-1`
- reported cell format remains `outer nonlinear iterations / cumulative linear iterations / runtime [s]`
- this note keeps only the `DCG + HYPRE` cut of the selected-tolerance replications
- `<ins><strong>...</strong></ins>` marks the fastest converged runtime in the row within this `DCG + HYPRE` cut
- `~~>300s~~` marks a nonlinear method stopped by the 300 s cap
- replication note: BC1 is used in the repository replication and reproduces the paper unknown count.

#### Paper-size mesh

- density: `12`
- measured elements: `103680`
- measured unknowns: `54886`

##### DCG with HYPRE

| `g_3` | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |
|---|---|---|---|
| `-3e-03` | `6/33/1.23` | `14/28/0.45` | <ins><strong>15/28/0.41</strong></ins> |
| `-4e-03` | `6/34/1.23` | <ins><strong>17/30/0.58</strong></ins> | `23/32/0.60` |
| `-5e-03` | `7/38/1.43` | <ins><strong>25/37/0.90</strong></ins> | `36/43/1.11` |
| `-6e-03` | <ins><strong>7/39/1.44</strong></ins> | `41/57/1.77` | `74/69/3.10` |
| `-7e-03` | <ins><strong>8/45/1.68</strong></ins> | `73/119/4.61` | `>200/141/16.02` |

#### ~2x unknowns mesh

- density: `15`
- measured elements: `202500`
- measured unknowns: `106048`

##### DCG with HYPRE

| `g_3` | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |
|---|---|---|---|
| `-3e-03` | `6/34/2.58` | `14/30/1.01` | <ins><strong>15/31/0.87</strong></ins> |
| `-4e-03` | `6/34/2.56` | <ins><strong>17/31/1.26</strong></ins> | `23/33/1.35` |
| `-5e-03` | `7/41/3.06` | <ins><strong>25/39/1.91</strong></ins> | `37/44/2.38` |
| `-6e-03` | <ins><strong>7/40/3.06</strong></ins> | `41/60/3.88` | `75/72/6.95` |
| `-7e-03` | <ins><strong>8/46/3.51</strong></ins> | `73/123/11.48` | `>200/144/37.53` |

#### ~200k unknowns mesh

- density: `19`
- measured elements: `411540`
- measured unknowns: `213520`

##### DCG with HYPRE

| `g_3` | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |
|---|---|---|---|
| `-3e-03` | `6/37/5.62` | `14/32/2.34` | <ins><strong>15/36/2.06</strong></ins> |
| `-4e-03` | `6/37/5.60` | `17/37/3.07` | <ins><strong>22/42/2.96</strong></ins> |
| `-5e-03` | `7/44/6.62` | <ins><strong>25/44/4.58</strong></ins> | `37/48/5.74` |
| `-6e-03` | <ins><strong>7/45/6.84</strong></ins> | `41/66/9.15` | `76/77/18.63` |
| `-7e-03` | <ins><strong>8/51/7.86</strong></ins> | `73/130/28.82` | `>200/148/95.83` |

## Example 2: Prismatic Bar with Non-Constant Thickness

### Geometry and Boundary Conditions

- Cross-section size varies linearly along `x_3`:

```text
rho(x_3) = rho_1 + (rho_2 - rho_1) * x_3 / L
```

- Domain:

```text
Omega = {
  (x_1, x_2, x_3) in R^3 |
  x_3 in (0, L),
  x_1, x_2 in (-rho(x_3)/2, rho(x_3)/2)
}
```

- Shell boundary condition in the paper text: `u_3 = 0`.
- Left face `x_3 = 0`: prescribed constant axial velocity
  - `u_3(x_1, x_2, 0) = u_3,0`
- Right face `x_3 = L`: homogeneous Neumann condition.
- Volume force:

```text
g = (0, 0, g_3),  g_3 < 0
```

### Parameters

- `L = 50`
- `rho_1 = 6`
- `rho_2 = 5`
- `u_3,0 = 5`
- `mu_0 = 1`
- `mu_infty = 0.001`
- `lambda = 10`
- `p = 1.1`
- varied load values:
  - `g_3 in {-3e-3, -4e-3, -5e-3, -6e-3, -7e-3, -8e-3}`

### Discretization

- `P1` elements
- uniform mesh
- `86,400` elements
- `50,986` unknowns

### Table from the Paper

| `g_3` | `IC-N` | `IC-qN1` | `IC-qN2` | `AGMG-N` | `AGMG-qN1` | `AGMG-qN2` |
|---|---|---|---|---|---|---|
| `-3e-3` | `5/86/1.28` | `12/72/0.76` | `14/56/0.52` | `5/39/1.52` | `12/33/0.90` | `13/39/0.74` |
| `-4e-3` | `6/100/1.57` | `16/77/0.93` | `17/59/0.59` | `6/50/1.80` | `16/37/1.06` | `17/51/0.93` |
| `-5e-3` | `6/101/1.56` | `19/83/1.15` | `23/72/0.73` | `6/50/1.88` | `19/43/1.35` | `24/70/1.25` |
| `-6e-3` | `7/119/1.85` | `25/90/1.56` | `36/97/1.07` | `7/61/2.21` | `25/53/1.87` | `37/95/1.80` |
| `-7e-3` | `7/120/1.87` | `38/102/2.32` | `68/162/1.77` | `7/61/2.17` | `37/68/2.71` | `70/167/3.09` |
| `-8e-3` | `8/139/2.20` | `62/126/3.87` | `168/347/4.15` | `8/71/2.52` | `61/93/4.44` | `198/359/7.28` |

Paper note: for `|g_3| = 0.008` or larger, `qN2` becomes the slowest; the
paper also states that `IC` is faster than `AGMG` in this example.

### Replications on This Machine

- nonlinear tolerance policy used here:
  - `N`: `1e-4`
  - `qN1`: `1e-1`
  - `qN2`: `1e-1`
- reported cell format remains `outer nonlinear iterations / cumulative linear iterations / runtime [s]`
- this note keeps only the `DCG + HYPRE` cut of the selected-tolerance replications
- `<ins><strong>...</strong></ins>` marks the fastest converged runtime in the row within this `DCG + HYPRE` cut
- `~~>300s~~` marks a nonlinear method stopped by the 300 s cap
- replication note: BC1 is used in the repository replication because it matches the paper table materially better than the literal BC3 interpretation.

#### Paper-size mesh

- density: `12`
- measured elements: `86400`
- measured unknowns: `45786`

##### DCG with HYPRE

| `g_3` | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |
|---|---|---|---|
| `-3e-03` | `5/26/0.82` | `11/22/0.34` | <ins><strong>12/25/0.28</strong></ins> |
| `-4e-03` | `6/31/1.02` | `16/25/0.44` | <ins><strong>18/25/0.40</strong></ins> |
| `-5e-03` | `6/31/1.01` | <ins><strong>19/28/0.54</strong></ins> | `24/31/0.57` |
| `-6e-03` | `7/36/1.17` | <ins><strong>25/35/0.84</strong></ins> | `37/41/0.98` |
| `-7e-03` | <ins><strong>7/37/1.19</strong></ins> | `38/50/1.38` | `67/61/2.33` |
| `-8e-03` | <ins><strong>8/43/1.38</strong></ins> | `62/93/3.31` | `161/123/8.83` |

#### ~2x unknowns mesh

- density: `15`
- measured elements: `168750`
- measured unknowns: `88448`

##### DCG with HYPRE

| `g_3` | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |
|---|---|---|---|
| `-3e-03` | `5/27/1.80` | `11/24/0.77` | <ins><strong>12/26/0.65</strong></ins> |
| `-4e-03` | `6/33/2.14` | `16/28/0.98` | <ins><strong>17/31/0.86</strong></ins> |
| `-5e-03` | `6/32/2.31` | <ins><strong>19/31/1.24</strong></ins> | `25/33/1.26` |
| `-6e-03` | `7/38/2.55` | <ins><strong>25/38/1.74</strong></ins> | `37/42/2.23` |
| `-7e-03` | <ins><strong>7/40/2.58</strong></ins> | `38/57/3.24` | `68/64/5.24` |
| `-8e-03` | <ins><strong>8/45/2.92</strong></ins> | `62/98/7.34` | `163/126/21.10` |

#### ~200k unknowns mesh

- density: `19`
- measured elements: `342228`
- measured unknowns: `177680`

##### DCG with HYPRE

| `g_3` | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |
|---|---|---|---|
| `-3e-03` | `5/29/3.95` | `11/25/1.84` | <ins><strong>13/29/1.70</strong></ins> |
| `-4e-03` | `6/36/4.83` | `16/32/2.43` | <ins><strong>17/34/2.09</strong></ins> |
| `-5e-03` | `6/35/4.67` | `19/34/3.09` | <ins><strong>25/35/2.85</strong></ins> |
| `-6e-03` | `7/41/5.47` | <ins><strong>25/42/4.16</strong></ins> | `38/45/5.45` |
| `-7e-03` | <ins><strong>7/40/5.38</strong></ins> | `38/61/7.34` | `69/67/13.28` |
| `-8e-03` | <ins><strong>8/48/6.36</strong></ins> | `62/97/16.79` | `166/132/52.34` |

## Mapping to This Repository

- `qN1` corresponds to [newton_quasi1.m](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/+NEWTON/newton_quasi1.m)
- `qN2` corresponds to [newton_quasi2.m](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/+NEWTON/newton_quasi2.m)
- The current benchmark entry point is [nonlinear_Stokes_3D.m](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/nonlinear_Stokes_3D.m)
- The replication tables in this note are assembled from [paper_tolerance_sweep_results.mat](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/results/paper_tolerance_sweep_results.mat)

Practical note: the current script defaults are still closer to the paper's Example 2
geometry (`size_xy_0 = 6`, `size_xy_L = 5`, `size_z = 50`, `gamma = 0.008`)
than to Example 1, but the linear solver defaults now use `1e-4` for `N` and
`1e-1` for `qN1/qN2`.
