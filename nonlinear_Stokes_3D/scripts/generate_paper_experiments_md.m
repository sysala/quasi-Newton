function generate_paper_experiments_md()
%GENERATE_PAPER_EXPERIMENTS_MD Refresh PAPER_EXPERIMENTS.md from the sweep checkpoint.

paths = ns3d_paths();
results_path = fullfile(paths.results, 'paper_tolerance_sweep_results.mat');
selected_results_path = fullfile(paths.results, 'selected_tolerance_mesh_results_dcg_hypre.mat');
markdown_path = fullfile(paths.root, 'PAPER_EXPERIMENTS.md');

loaded = load(results_path, 'results');
if ~isfield(loaded, 'results')
    error('generate_paper_experiments_md:missingResults', ...
        'Expected a results variable in %s.', results_path);
end
results = loaded.results;
include_extended_mesh = false;
if isfile(selected_results_path)
    loaded_selected = load(selected_results_path, 'results');
    if isfield(loaded_selected, 'results') && ~isempty(loaded_selected.results)
        results = [results, loaded_selected.results];
        include_extended_mesh = true;
    end
end

tol_newton = 1e-4;
tol_quasi = 1e-1;
examples = local_examples();
mesh_variants = local_mesh_variants(include_extended_mesh);
families = local_linear_families();
methods = local_method_specs(tol_newton, tol_quasi);

fid = fopen(markdown_path, 'w');
if fid < 0
    error('generate_paper_experiments_md:openFailed', ...
        'Unable to open %s for writing.', markdown_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Paper Experiments for `nonlinear_Stokes_3D`\n\n');
fprintf(fid, 'Source: [paper.tex](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/paper.tex)\n\n');
fprintf(fid, 'This note pulls out the numerical experiment settings and the paper tables from\n');
fprintf(fid, '`paper.tex` into a compact markdown reference.\n\n');

fprintf(fid, '## Problem Setup Used in the Paper\n\n');
fprintf(fid, '- Problem type: simplified nonlinear Stokes problem in 3D.\n');
fprintf(fid, '- Simplification: the pore-pressure field is given and the divergence-free\n');
fprintf(fid, '  constraint is not enforced, so the structure is close to nonlinear\n');
fprintf(fid, '  elasticity.\n');
fprintf(fid, '- Constitutive law: Carreau law\n\n');
fprintf(fid, '```text\n');
fprintf(fid, 'a(r) = mu_infty + (mu_0 - mu_infty) * (1 + lambda * r^2)^((p - 2) / 2)\n');
fprintf(fid, '```\n\n');
fprintf(fid, '- Material parameter range in theory:\n');
fprintf(fid, '  - `p in (1, 2)`\n');
fprintf(fid, '  - `lambda > 0`\n');
fprintf(fid, '  - `mu_0 > mu_infty > 0`\n');
fprintf(fid, '- Nonlinear solvers compared:\n');
fprintf(fid, '  - `N`: Newton\n');
fprintf(fid, '  - `qN1`: quasi-Newton 1, corresponding to the variable preconditioner from\n');
fprintf(fid, '    Section 3.1\n');
fprintf(fid, '  - `qN2`: quasi-Newton 2, corresponding to the fixed preconditioner from\n');
fprintf(fid, '    Section 3.2\n');
fprintf(fid, '- Inner linear solver setup in the paper:\n');
fprintf(fid, '  - deflated and preconditioned conjugate gradient\n');
fprintf(fid, '  - preconditioner choices: incomplete Cholesky (`IC`) and AGMG\n');
fprintf(fid, '- Current repository defaults for these benchmark comparisons:\n');
fprintf(fid, '  - `N`: linear solver tolerance `1e-4`\n');
fprintf(fid, '  - `qN1`: linear solver tolerance `1e-1`\n');
fprintf(fid, '  - `qN2`: linear solver tolerance `1e-1`\n');
fprintf(fid, '- Finite elements: `P1`\n');
fprintf(fid, '- Reported table format:\n');
fprintf(fid, '  - `outer nonlinear iterations / cumulative linear iterations / runtime [s]`\n\n');

local_write_example_1(fid, results, examples(1), mesh_variants, families, methods);
local_write_example_2(fid, results, examples(2), mesh_variants, families, methods);

fprintf(fid, '## Mapping to This Repository\n\n');
fprintf(fid, '- `qN1` corresponds to [newton_quasi1.m](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/+NEWTON/newton_quasi1.m)\n');
fprintf(fid, '- `qN2` corresponds to [newton_quasi2.m](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/+NEWTON/newton_quasi2.m)\n');
fprintf(fid, '- The current benchmark entry point is [nonlinear_Stokes_3D.m](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/nonlinear_Stokes_3D.m)\n');
fprintf(fid, '- The replication tables in this note are assembled from [paper_tolerance_sweep_results.mat](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/results/paper_tolerance_sweep_results.mat)\n\n');

fprintf(fid, 'Practical note: the current script defaults are still closer to the paper''s Example 2\n');
fprintf(fid, 'geometry (`size_xy_0 = 6`, `size_xy_L = 5`, `size_z = 50`, `gamma = 0.008`)\n');
fprintf(fid, 'than to Example 1, but the linear solver defaults now use `1e-4` for `N` and\n');
fprintf(fid, '`1e-1` for `qN1/qN2`.\n');
end

function local_write_example_1(fid, results, example, mesh_variants, families, methods)
fprintf(fid, '## Example 1: Prismatic Bar\n\n');
fprintf(fid, '### Geometry and Boundary Conditions\n\n');
fprintf(fid, '- Domain:\n\n');
fprintf(fid, '```text\n');
fprintf(fid, 'Omega = (0, rho) x (0, rho) x (0, L)\n');
fprintf(fid, '```\n\n');
fprintf(fid, '- Bar aligned with the `x_3` axis.\n');
fprintf(fid, '- Shell boundary condition: zero velocity in the normal direction\n');
fprintf(fid, '  (paper labels this as slip conditions).\n');
fprintf(fid, '- Left face `x_3 = 0`: prescribed constant axial velocity\n');
fprintf(fid, '  - `u_3(x_1, x_2, 0) = u_3,0`\n');
fprintf(fid, '- Right face `x_3 = L`: homogeneous Neumann condition.\n');
fprintf(fid, '- Volume force:\n\n');
fprintf(fid, '```text\n');
fprintf(fid, 'g = (0, 0, g_3),  g_3 < 0\n');
fprintf(fid, '```\n\n');

fprintf(fid, '### Parameters\n\n');
fprintf(fid, '- `L = 50`\n');
fprintf(fid, '- `rho = 5`\n');
fprintf(fid, '- `u_3,0 = 5`\n');
fprintf(fid, '- `mu_0 = 1`\n');
fprintf(fid, '- `mu_infty = 0.001`\n');
fprintf(fid, '- `lambda = 10`\n');
fprintf(fid, '- `p = 1.1`\n');
fprintf(fid, '- varied load values:\n');
fprintf(fid, '  - `g_3 in {-3e-3, -4e-3, -5e-3, -6e-3, -7e-3}`\n\n');

fprintf(fid, '### Discretization\n\n');
fprintf(fid, '- `P1` elements\n');
fprintf(fid, '- uniform mesh\n');
fprintf(fid, '- `103,680` elements\n');
fprintf(fid, '- `54,886` unknowns\n\n');

fprintf(fid, '### Table from the Paper\n\n');
local_write_original_table(fid, "example1");
fprintf(fid, 'Paper note: for smaller `|g_3|`, `qN2` is fastest; for `|g_3| = 0.007` or\n');
fprintf(fid, 'larger, `qN2` becomes the slowest because the material response is strongly\n');
fprintf(fid, 'nonlinear.\n\n');

local_write_replication_section(fid, results, example, mesh_variants, families, methods);
end

function local_write_example_2(fid, results, example, mesh_variants, families, methods)
fprintf(fid, '## Example 2: Prismatic Bar with Non-Constant Thickness\n\n');
fprintf(fid, '### Geometry and Boundary Conditions\n\n');
fprintf(fid, '- Cross-section size varies linearly along `x_3`:\n\n');
fprintf(fid, '```text\n');
fprintf(fid, 'rho(x_3) = rho_1 + (rho_2 - rho_1) * x_3 / L\n');
fprintf(fid, '```\n\n');
fprintf(fid, '- Domain:\n\n');
fprintf(fid, '```text\n');
fprintf(fid, 'Omega = {\n');
fprintf(fid, '  (x_1, x_2, x_3) in R^3 |\n');
fprintf(fid, '  x_3 in (0, L),\n');
fprintf(fid, '  x_1, x_2 in (-rho(x_3)/2, rho(x_3)/2)\n');
fprintf(fid, '}\n');
fprintf(fid, '```\n\n');
fprintf(fid, '- Shell boundary condition in the paper text: `u_3 = 0`.\n');
fprintf(fid, '- Left face `x_3 = 0`: prescribed constant axial velocity\n');
fprintf(fid, '  - `u_3(x_1, x_2, 0) = u_3,0`\n');
fprintf(fid, '- Right face `x_3 = L`: homogeneous Neumann condition.\n');
fprintf(fid, '- Volume force:\n\n');
fprintf(fid, '```text\n');
fprintf(fid, 'g = (0, 0, g_3),  g_3 < 0\n');
fprintf(fid, '```\n\n');

fprintf(fid, '### Parameters\n\n');
fprintf(fid, '- `L = 50`\n');
fprintf(fid, '- `rho_1 = 6`\n');
fprintf(fid, '- `rho_2 = 5`\n');
fprintf(fid, '- `u_3,0 = 5`\n');
fprintf(fid, '- `mu_0 = 1`\n');
fprintf(fid, '- `mu_infty = 0.001`\n');
fprintf(fid, '- `lambda = 10`\n');
fprintf(fid, '- `p = 1.1`\n');
fprintf(fid, '- varied load values:\n');
fprintf(fid, '  - `g_3 in {-3e-3, -4e-3, -5e-3, -6e-3, -7e-3, -8e-3}`\n\n');

fprintf(fid, '### Discretization\n\n');
fprintf(fid, '- `P1` elements\n');
fprintf(fid, '- uniform mesh\n');
fprintf(fid, '- `86,400` elements\n');
fprintf(fid, '- `50,986` unknowns\n\n');

fprintf(fid, '### Table from the Paper\n\n');
local_write_original_table(fid, "example2");
fprintf(fid, 'Paper note: for `|g_3| = 0.008` or larger, `qN2` becomes the slowest; the\n');
fprintf(fid, 'paper also states that `IC` is faster than `AGMG` in this example.\n\n');

local_write_replication_section(fid, results, example, mesh_variants, families, methods);
end

function local_write_replication_section(fid, results, example, mesh_variants, families, methods)
fprintf(fid, '### Replications on This Machine\n\n');
fprintf(fid, '- nonlinear tolerance policy used here:\n');
fprintf(fid, '  - `N`: `1e-4`\n');
fprintf(fid, '  - `qN1`: `1e-1`\n');
fprintf(fid, '  - `qN2`: `1e-1`\n');
fprintf(fid, '- reported cell format remains `outer nonlinear iterations / cumulative linear iterations / runtime [s]`\n');
fprintf(fid, '- this note keeps only the `DCG + HYPRE` cut of the selected-tolerance replications\n');
fprintf(fid, '- `<ins><strong>...</strong></ins>` marks the fastest converged runtime in the row within this `DCG + HYPRE` cut\n');
fprintf(fid, '- `~~>300s~~` marks a nonlinear method stopped by the 300 s cap\n');
fprintf(fid, '- replication note: %s\n\n', example.replication_note);

for i_mesh = 1:numel(mesh_variants)
    mesh_variant = mesh_variants(i_mesh);
    measured = local_find_first_result(results, example.id, mesh_variant.label);

    fprintf(fid, '#### %s\n\n', mesh_variant.display_name);
    fprintf(fid, '- density: `%d`\n', mesh_variant.density);
    if ~isempty(measured)
        fprintf(fid, '- measured elements: `%d`\n', measured.n_elements);
        fprintf(fid, '- measured unknowns: `%d`\n', measured.n_unknowns);
    end
    fprintf(fid, '\n');

    for i_family = 1:numel(families)
        family = families(i_family);
        fprintf(fid, '##### %s\n\n', family.display_name);
        fprintf(fid, '| `g_3` | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |\n');
        fprintf(fid, '|---|---|---|---|\n');

        for i_gamma = 1:numel(example.gammas)
            gamma = example.gammas(i_gamma);
            row_cells = local_collect_row_cells(results, example.id, mesh_variant.label, gamma, families, methods);
            fprintf(fid, '| `-%s` ', local_number_label(gamma));
            for i_prec = 1:numel(family.preconditioner_labels)
                for i_method = 1:numel(methods)
                    cell_data = row_cells(i_family, i_prec, i_method);
                    fprintf(fid, '| %s ', local_markdown_cell(cell_data));
                end
            end
            fprintf(fid, '|\n');
        end
        fprintf(fid, '\n');
    end
end
end

function examples = local_examples()
examples = struct([]);

examples(1).id = "example1";
examples(1).gammas = [0.003, 0.004, 0.005, 0.006, 0.007];
examples(1).replication_note = ...
    "BC1 is used in the repository replication and reproduces the paper unknown count.";

examples(2).id = "example2";
examples(2).gammas = [0.003, 0.004, 0.005, 0.006, 0.007, 0.008];
examples(2).replication_note = ...
    "BC1 is used in the repository replication because it matches the paper table materially better than the literal BC3 interpretation.";
end

function mesh_variants = local_mesh_variants(include_extended_mesh)
mesh_variants = struct([]);

mesh_variants(1).label = "paper_mesh";
mesh_variants(1).display_name = "Paper-size mesh";
mesh_variants(1).density = 12;

mesh_variants(2).label = "double_unknowns_mesh";
mesh_variants(2).display_name = "~2x unknowns mesh";
mesh_variants(2).density = 15;

if include_extended_mesh
    mesh_variants(3).label = "roughly_200k_unknowns_mesh";
    mesh_variants(3).display_name = "~200k unknowns mesh";
    mesh_variants(3).density = 19;
end
end

function families = local_linear_families()
families = struct([]);

families(1).display_name = "DCG with HYPRE";
families(1).solver_types = ["DCG_HYPRE_BOOMERAMG"];
families(1).preconditioner_labels = ["HYPRE"];
end

function methods = local_method_specs(tol_newton, tol_quasi)
methods = struct([]);

methods(1).field = "newton";
methods(1).label = "N";
methods(1).tol = tol_newton;
methods(1).iteration_limit = 100;

methods(2).field = "qn1";
methods(2).label = "qN1";
methods(2).tol = tol_quasi;
methods(2).iteration_limit = 100;

methods(3).field = "qn2";
methods(3).label = "qN2";
methods(3).tol = tol_quasi;
methods(3).iteration_limit = 200;
end

function row_cells = local_collect_row_cells(results, example_id, mesh_label, gamma, families, methods)
template = struct('label', '-', 'runtime', inf, 'converged', false, 'timed_out', false, 'highlight', false);
row_cells = repmat(template, numel(families), numel(families(1).solver_types), numel(methods));

best_runtime = inf;
for i_family = 1:numel(families)
    for i_solver = 1:numel(families(i_family).solver_types)
        solver_type = families(i_family).solver_types(i_solver);
        for i_method = 1:numel(methods)
            entry = local_find_result(results, example_id, mesh_label, solver_type, methods(i_method).tol, gamma);
            if isempty(entry) || ~entry.succeeded
                row_cells(i_family, i_solver, i_method).label = 'ERROR';
                continue
            end

            solver_result = entry.(methods(i_method).field);
            row_cells(i_family, i_solver, i_method).label = ...
                local_solver_cell(solver_result, methods(i_method).iteration_limit);
            row_cells(i_family, i_solver, i_method).runtime = solver_result.runtime_seconds;
            row_cells(i_family, i_solver, i_method).converged = logical(solver_result.converged);
            row_cells(i_family, i_solver, i_method).timed_out = ...
                isfield(solver_result, 'timed_out') && logical(solver_result.timed_out);

            if solver_result.converged
                best_runtime = min(best_runtime, solver_result.runtime_seconds);
            end
        end
    end
end

if ~isfinite(best_runtime)
    return
end

tolerance = max(1e-9, 1e-9 * abs(best_runtime));
for i_family = 1:numel(families)
    for i_solver = 1:numel(families(i_family).solver_types)
        for i_method = 1:numel(methods)
            if row_cells(i_family, i_solver, i_method).converged && ...
                    abs(row_cells(i_family, i_solver, i_method).runtime - best_runtime) <= tolerance
                row_cells(i_family, i_solver, i_method).highlight = true;
            end
        end
    end
end
end

function entry = local_find_result(results, example_id, mesh_label, solver_type, tol, gamma)
entry = [];
for i = 1:numel(results)
    if string(results(i).example_id) ~= string(example_id)
        continue
    end
    if string(results(i).mesh_label) ~= string(mesh_label)
        continue
    end
    if string(results(i).solver_type) ~= string(solver_type)
        continue
    end
    if abs(results(i).linear_solver_tolerance - tol) > eps(max(1, abs(tol)))
        continue
    end
    if abs(results(i).gamma - gamma) > eps(max(1, abs(gamma)))
        continue
    end
    entry = results(i);
    return
end
end

function entry = local_find_first_result(results, example_id, mesh_label)
entry = [];
for i = 1:numel(results)
    if string(results(i).example_id) == string(example_id) && ...
            string(results(i).mesh_label) == string(mesh_label) && ...
            results(i).succeeded
        entry = results(i);
        return
    end
end
end

function label = local_solver_cell(solver_result, iteration_limit)
if isfield(solver_result, 'timed_out') && solver_result.timed_out
    label = '>300s';
    return
end
if solver_result.converged
    label = sprintf('%d/%d/%.2f', ...
        solver_result.iterations, ...
        solver_result.linear_iterations, ...
        solver_result.runtime_seconds);
else
    label = sprintf('>%d/%d/%.2f', ...
        iteration_limit, ...
        solver_result.linear_iterations, ...
        solver_result.runtime_seconds);
end
end

function rendered = local_markdown_cell(cell_data)
if cell_data.highlight
    rendered = sprintf('<ins><strong>%s</strong></ins>', cell_data.label);
elseif cell_data.timed_out
    rendered = '~~>300s~~';
elseif strcmp(cell_data.label, 'ERROR') || strcmp(cell_data.label, '-')
    rendered = cell_data.label;
else
    rendered = sprintf('`%s`', cell_data.label);
end
end

function label = local_number_label(value)
label = sprintf('%.0e', value);
end

function local_write_original_table(fid, example_id)
switch string(example_id)
    case "example1"
        fprintf(fid, '| `g_3` | `IC-N` | `IC-qN1` | `IC-qN2` | `AGMG-N` | `AGMG-qN1` | `AGMG-qN2` |\n');
        fprintf(fid, '|---|---|---|---|---|---|---|\n');
        fprintf(fid, '| `-3e-3` | `6/112/1.86` | `14/84/0.95` | `15/58/0.61` | `6/51/2.12` | `14/35/0.99` | `15/44/0.85` |\n');
        fprintf(fid, '| `-4e-3` | `6/111/1.84` | `17/90/1.17` | `22/70/0.79` | `6/53/2.01` | `17/40/1.25` | `22/57/1.11` |\n');
        fprintf(fid, '| `-5e-3` | `7/130/2.17` | `27/97/1.72` | `35/93/1.12` | `7/63/2.44` | `25/50/1.85` | `36/86/1.68` |\n');
        fprintf(fid, '| `-6e-3` | `7/133/2.17` | `41/111/2.69` | `74/167/2.11` | `7/63/2.39` | `41/70/3.05` | `75/170/3.33` |\n');
        fprintf(fid, '| `-7e-3` | `8/152/2.53` | `73/143/5.17` | `>200/-/5.53` | `8/76/2.84` | `72/108/5.63` | `>200/-/6.81` |\n\n');
    case "example2"
        fprintf(fid, '| `g_3` | `IC-N` | `IC-qN1` | `IC-qN2` | `AGMG-N` | `AGMG-qN1` | `AGMG-qN2` |\n');
        fprintf(fid, '|---|---|---|---|---|---|---|\n');
        fprintf(fid, '| `-3e-3` | `5/86/1.28` | `12/72/0.76` | `14/56/0.52` | `5/39/1.52` | `12/33/0.90` | `13/39/0.74` |\n');
        fprintf(fid, '| `-4e-3` | `6/100/1.57` | `16/77/0.93` | `17/59/0.59` | `6/50/1.80` | `16/37/1.06` | `17/51/0.93` |\n');
        fprintf(fid, '| `-5e-3` | `6/101/1.56` | `19/83/1.15` | `23/72/0.73` | `6/50/1.88` | `19/43/1.35` | `24/70/1.25` |\n');
        fprintf(fid, '| `-6e-3` | `7/119/1.85` | `25/90/1.56` | `36/97/1.07` | `7/61/2.21` | `25/53/1.87` | `37/95/1.80` |\n');
        fprintf(fid, '| `-7e-3` | `7/120/1.87` | `38/102/2.32` | `68/162/1.77` | `7/61/2.17` | `37/68/2.71` | `70/167/3.09` |\n');
        fprintf(fid, '| `-8e-3` | `8/139/2.20` | `62/126/3.87` | `168/347/4.15` | `8/71/2.52` | `61/93/4.44` | `198/359/7.28` |\n\n');
    otherwise
        error('generate_paper_experiments_md:unknownExample', ...
            'Unsupported example id "%s".', example_id);
end
end
