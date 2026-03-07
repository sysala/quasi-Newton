function generate_paper_experiments_md()
%GENERATE_PAPER_EXPERIMENTS_MD Rebuild PAPER_EXPERIMENTS.md from local sweeps.

paths = mp2d_paths();
results_path = fullfile(paths.results, 'paper_tolerance_sweep_results.mat');
selected_results_path = fullfile(paths.results, 'selected_tolerance_mesh_results_dcg_hypre.mat');
markdown_path = fullfile(paths.root, 'PAPER_EXPERIMENTS.md');

loaded = load(results_path, 'results');
if ~isfield(loaded, 'results')
    error('generate_paper_experiments_md:missingResults', ...
        'Expected a results variable in %s.', results_path);
end
results = loaded.results;

selected_results = struct([]);
has_selected_results = false;
if isfile(selected_results_path)
    loaded_selected = load(selected_results_path, 'results');
    if isfield(loaded_selected, 'results') && ~isempty(loaded_selected.results)
        selected_results = loaded_selected.results;
        has_selected_results = true;
    end
end

tolerances = [2e-1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5];
cases = local_cases();
mesh_variants = local_mesh_variants();

fid = fopen(markdown_path, 'w');
if fid < 0
    error('generate_paper_experiments_md:openFailed', ...
        'Unable to open %s for writing.', markdown_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Paper Experiments for `magnetic_potential_problem_2D`\n\n');
fprintf(fid, 'Source: [magnetic_problem.m](magnetic_problem.m) and case wrappers in [scripts/](scripts/).\n\n');
fprintf(fid, 'This note summarizes the magnetic benchmark setup and the replicated\n');
fprintf(fid, 'results in this repository from local checkpoints.\n\n');

fprintf(fid, '## Problem Setup (2D scalar magnetic potential)\n\n');
fprintf(fid, '- Unknown: scalar potential (u) in 2D\n');
fprintf(fid, '- Constitutive map: S = a(|nabla u|) * grad u\n');
fprintf(fid, '- Scalar conductivity law on the heterogeneous subregion:\n');
fprintf(fid, '  - `a(r) = alpha + (1-alpha) * r^8 / (r^8 + beta)`\n');
fprintf(fid, '- Global source loading strength is set by `rho`.\n');
fprintf(fid, '- Outer methods: `newton`, `newton_quasi1`, `newton_quasi2`.\n');
fprintf(fid, '- Mesh family in this benchmark: `P1` and `P2` are supported; all paper-mimicking\n');
fprintf(fid, '  sweeps in this folder use `P1`.\n');
fprintf(fid, '- Linear inner solver policy used in replication outputs:\n');
fprintf(fid, '  - `DCG_HYPRE_BOOMERAMG` (`DCG + HYPRE`), with scalar BoomerAMG setup\n');
fprintf(fid, '    (`num_functions = 1`, no elasticity near-nullspace).\n\n');

fprintf(fid, 'Note: the parameter set below is the benchmark default used by this\n');
fprintf(fid, 'repository workflow; no paper source table is checked into this snapshot.\n\n');

fprintf(fid, '## Replication Tables from Tolerance Sweep (Full)\n\n');
tol_labels = local_format_tolerances(tolerances);
fprintf(fid, '- Tolerances used for sweep: `%s`\n', strjoin(tol_labels, '`, `'));
fprintf(fid, '- Each table cell: `nonlinear_iterations / cumulative_linear_iterations / runtime [s]`\n');
fprintf(fid, '- A `~` prefix marks a nonconverged/non-fully-validated outer solve entry.\n\n');

for i_case = 1:numel(cases)
    case_entry = cases(i_case);
    fprintf(fid, '### %s\n\n', case_entry.name);
    fprintf(fid, '- case id: `%s`\n', case_entry.id);
    fprintf(fid, '- base geometry: `%g x %g`\n', case_entry.size_x, case_entry.size_y);
    fprintf(fid, '- material parameters: `alpha=%g`, `beta=%g`, `rho=%g`, `kappa=%g`, `delta=%g`\n\n', ...
        case_entry.alpha, case_entry.beta, case_entry.rho, case_entry.kappa, case_entry.delta);

    for i_mesh = 1:numel(mesh_variants)
        mesh_variant = mesh_variants(i_mesh);
        fprintf(fid, '#### %s\n\n', mesh_variant.display_name);
        fprintf(fid, '| linear tolerance | Newton | qN1 | qN2 |\n');
        fprintf(fid, '|---|---|---|---|\n');
        for i_tol = 1:numel(tolerances)
            tol = tolerances(i_tol);
            entry = local_find_entry(results, case_entry.id, mesh_variant.label, tol);
            if isempty(entry)
                fprintf(fid, '| `%s` | `-` | `-` | `-` |\n', local_format_tolerance(tol));
                continue
            end
            newton_entry = local_format_solver_cell(entry.newton);
            qn1_entry = local_format_solver_cell(entry.qn1);
            qn2_entry = local_format_solver_cell(entry.qn2);
            fprintf(fid, '| `%s` | %s | %s | %s |\n', ...
                local_format_tolerance(tol), newton_entry, qn1_entry, qn2_entry);
        end
        fprintf(fid, '\n');
    end
end

fprintf(fid, '## Selected-Tolerance Replications (DCG + HYPRE)\n\n');
if has_selected_results
    local_write_selected_tolerance_replications(fid, selected_results, cases, local_selected_mesh_variants());
else
    fprintf(fid, 'No `selected_tolerance_mesh_results_dcg_hypre.mat` checkpoint was found, so selected\n');
    fprintf(fid, 'policy replications were not yet generated.\n\n');
end

fprintf(fid, '## Mapping to This Repository\n\n');
fprintf(fid, '- Benchmark entry point: [magnetic_problem.m](magnetic_problem.m)\n');
fprintf(fid, '- Benchmark wrapper: [run_paper_case.m](scripts/run_paper_case.m)\n');
fprintf(fid, '- Sweep controller: [run_paper_tolerance_sweep.m](scripts/run_paper_tolerance_sweep.m)\n');
fprintf(fid, '- Selected replication runner:\n');
fprintf(fid, '  - [run_selected_tolerance_mesh_extension.m](scripts/run_selected_tolerance_mesh_extension.m)\n');
fprintf(fid, '  - [run_selected_tolerance_mesh_method.m](scripts/run_selected_tolerance_mesh_method.m)\n');
fprintf(fid, '- Replication checkpoints:\n');
fprintf(fid, '  - [paper_tolerance_sweep_results.mat](results/paper_tolerance_sweep_results.mat)\n');
fprintf(fid, '  - [selected_tolerance_mesh_results_dcg_hypre.mat](results/selected_tolerance_mesh_results_dcg_hypre.mat)\n');
fprintf(fid, '- Plot generator: [plot_paper_tolerance_sweep.m](scripts/plot_paper_tolerance_sweep.m)\n');
fprintf(fid, '- Current report source: [generate_paper_experiments_md.m](scripts/generate_paper_experiments_md.m)\n');
end

function cases = local_cases()
common = struct( ...
    'elem_type', "P1", ...
    'level', 3, ...
    'size_x', 1, ...
    'size_y', 1, ...
    'alpha', 3e-4, ...
    'beta', 16000, ...
    'rho', 40, ...
    'kappa', 25, ...
    'delta', 1, ...
    'run_newton', true, ...
    'run_qn1', true, ...
    'run_qn2', true, ...
    'run_postprocess', false, ...
    'linear_solver_maxit', 1000, ...
    'deflation_basis_tolerance', 1e-10, ...
    'boomeramg_opts', struct('threads', 16, 'print_level', 0, 'use_as_preconditioner', true, 'num_functions', 1), ...
    'nonlinear_method_max_runtime_seconds', 120);

cases = repmat(common, 1, 2);
cases(1).id = "case_delta_1";
cases(1).name = "Inclusion contrast δ=1";
cases(1).delta = 1;

cases(2).id = "case_delta_2";
cases(2).name = "Inclusion contrast δ=2";
cases(2).delta = 2;
end

function mesh_variants = local_mesh_variants()
mesh_variants(1).label = "level_3";
mesh_variants(1).display_name = "Level 3 mesh";
mesh_variants(1).level = 3;
mesh_variants(2).label = "level_4";
mesh_variants(2).display_name = "Level 4 mesh";
mesh_variants(2).level = 4;
end

function mesh_variants = local_selected_mesh_variants()
mesh_variants(1).label = "level_4";
mesh_variants(1).display_name = "Level 4 mesh";
mesh_variants(1).level = 4;
end

function entry = local_find_entry(results, case_id, mesh_label, tol)
entry = [];
for i = 1:numel(results)
    if string(results(i).case_id) ~= string(case_id)
        continue
    end
    if string(results(i).mesh_label) ~= string(mesh_label)
        continue
    end
    if abs(results(i).linear_solver_tolerance - tol) > eps(max(1, abs(tol)))
        continue
    end
    entry = results(i);
    return
end
end

function entry = local_find_selected_entry(results, case_id, mesh_label, tol)
entry = [];
for i = 1:numel(results)
    if string(results(i).case_id) ~= string(case_id)
        continue
    end
    if string(results(i).mesh_label) ~= string(mesh_label)
        continue
    end
    if abs(results(i).linear_solver_tolerance - tol) > eps(max(1, abs(tol)))
        continue
    end
    entry = results(i);
    return
end
end

function measured = local_find_first_entry(results, case_id, mesh_label)
measured = [];
for i = 1:numel(results)
    if string(results(i).case_id) == string(case_id) && string(results(i).mesh_label) == string(mesh_label) && results(i).succeeded
        measured = results(i);
        return
    end
end
end

function local_write_selected_tolerance_replications(fid, selected_results, cases, mesh_variants)
tol_newton = 1e-4;
tol_quasi = 1e-1;
fprintf(fid, 'These tables use the policy: `N` at `1e-4`, `qN1/qN2` at `1e-1`.\n');
fprintf(fid, 'Only `DCG + HYPRE` is tracked here to match the selected replication cut.\n');
fprintf(fid, 'Cell format remains `nonlinear_iterations / cumulative_linear_iterations / runtime [s]`.\n');
fprintf(fid, 'A `~` prefix means the reported nonlinear solve did not satisfy the stored convergence criterion.\n\n');

for i_case = 1:numel(cases)
    case_entry = cases(i_case);
    fprintf(fid, '### %s\n\n', case_entry.name);

    for i_mesh = 1:numel(mesh_variants)
        mesh_variant = mesh_variants(i_mesh);
        measured = local_find_first_entry(selected_results, case_entry.id, mesh_variant.label);

        fprintf(fid, '#### %s\n\n', mesh_variant.display_name);
        fprintf(fid, '- mesh level: `%d`\n', mesh_variant.level);
        if ~isempty(measured)
            fprintf(fid, '- measured elements: `%d`\n', measured.n_elements);
            fprintf(fid, '- measured unknowns: `%d`\n\n', measured.n_unknowns);
        else
            fprintf(fid, '\n');
        end

        fprintf(fid, '| method family | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |\n');
        fprintf(fid, '|---|---|---|---|\n');
        newton_entry = local_find_selected_entry(selected_results, case_entry.id, mesh_variant.label, tol_newton);
        quasi_entry = local_find_selected_entry(selected_results, case_entry.id, mesh_variant.label, tol_quasi);

        if isempty(newton_entry)
            n_cell = local_format_solver_cell([]);
        else
            n_cell = local_format_solver_cell(newton_entry.newton);
        end
        if isempty(quasi_entry)
            q1_cell = local_format_solver_cell([]);
            q2_cell = local_format_solver_cell([]);
        else
            q1_cell = local_format_solver_cell(quasi_entry.qn1);
            q2_cell = local_format_solver_cell(quasi_entry.qn2);
        end
        fprintf(fid, '| DCG | %s | %s | %s |\n\n', n_cell, q1_cell, q2_cell);
    end
end

end

function text = local_format_solver_cell(solver_result)
if isempty(solver_result) || ~isfield(solver_result, 'iterations') || ~isfinite(solver_result.iterations)
    text = '`-`';
    return
end

if ~isfinite(solver_result.runtime_seconds)
    text = '`-`';
    return
end

if isfield(solver_result, 'timed_out') && solver_result.timed_out
    text = sprintf('`%d / %d / >%g`', ...
        solver_result.iterations, round(solver_result.linear_iterations), solver_result.runtime_seconds);
    return
end

base = sprintf('%d / %d / %.3g', ...
    solver_result.iterations, round(solver_result.linear_iterations), solver_result.runtime_seconds);
if isfield(solver_result, 'converged') && solver_result.converged
    text = sprintf('`%s`', base);
else
    text = sprintf('`~%s`', base);
end
end

function tol_label = local_format_tolerance(tol)
tol_label = sprintf('%.0e', tol);
end

function tol_labels = local_format_tolerances(tolerances)
tol_labels = arrayfun(@local_format_tolerance, tolerances, 'UniformOutput', false);
end
