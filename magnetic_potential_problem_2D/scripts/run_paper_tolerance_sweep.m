function run_paper_tolerance_sweep(options)
%RUN_PAPER_TOLERANCE_SWEEP Sweep the magnetic benchmark over inner linear tolerance.
%
% The sweep is resumable. Results are stored in:
%   results/paper_tolerance_sweep_results.mat
% and the markdown report is regenerated at:
%   PAPER_TOLERANCE_SWEEP.md
% Optional input:
%   options - struct with:
%     run_newton, run_qn1, run_qn2 : logical scalar overrides
%     use_damping : logical scalar

if nargin < 1 || isempty(options)
    options = struct();
end

paths = mp2d_paths();
results_path = fullfile(paths.results, 'paper_tolerance_sweep_results.mat');
markdown_path = fullfile(paths.root, 'PAPER_TOLERANCE_SWEEP.md');

cases = local_cases();
mesh_variants = local_mesh_variants();
tolerances = [2e-1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5];

plan = local_build_plan(cases, mesh_variants, tolerances, options);

results = struct([]);
if isfile(results_path)
    loaded = load(results_path, 'results');
    if isfield(loaded, 'results')
        results = loaded.results;
    end
end

force_rerun = isfield(options, 'force_rerun') && logical(options.force_rerun);
regenerate_only = isfield(options, 'regenerate_only') && logical(options.regenerate_only);
skip_plots = isfield(options, 'skip_plots') && logical(options.skip_plots);

if ~regenerate_only
    for i = 1:numel(plan)
        plan_entry = plan(i);
        key = plan_entry.key;
        if ~force_rerun && local_find_result_index(results, key) > 0
            fprintf('[%d/%d] skip %s\n', i, numel(plan), key);
            continue
        end

        local_clear_hypre_instances();
        fprintf('[%d/%d] run %s\n', i, numel(plan), key);
        try
            case_result = run_paper_case(plan_entry.cfg);
            entry = local_result_entry(plan_entry, case_result, "", true);
        catch ME
            entry = local_result_entry(plan_entry, local_failure_result(), ...
                string(getReport(ME, 'extended', 'hyperlinks', 'off')), false);
        end
        local_clear_hypre_instances();

        if isempty(results)
            results = entry;
        else
            idx = local_find_result_index(results, key);
            if idx > 0
                results(idx) = entry;
            else
                results(end + 1) = entry; %#ok<AGROW>
            end
        end

        save(results_path, 'results', 'plan', '-v7.3');
        local_write_markdown(markdown_path, results, plan, cases, tolerances, mesh_variants);
end
elseif isempty(results)
    error('run_paper_tolerance_sweep:noCheckpoint', ...
        'regenerate_only=true requires an existing results checkpoint at %s.', results_path);
end

local_write_markdown(markdown_path, results, plan, cases, tolerances, mesh_variants);
if ~skip_plots
    try
        plot_paper_tolerance_sweep(results_path);
    catch ME
        warning('run_paper_tolerance_sweep:plotFailed', ...
            'Could not regenerate plots: %s', string(ME.message));
    end
end
end

function local_clear_hypre_instances()
try
    LINEAR_SOLVERS.hypre_boomeramg_clear();
catch
end
end

function cases = local_cases()
base = struct( ...
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

cases = repmat(base, 1, 2);

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

function plan = local_build_plan(cases, mesh_variants, tolerances, options)
plan = struct([]);
counter = 0;

force_newton = local_option_bool(options, 'run_newton', []);
force_qn1 = local_option_bool(options, 'run_qn1', []);
force_qn2 = local_option_bool(options, 'run_qn2', []);
force_damping = local_option_bool(options, 'use_damping', false);

for i_case = 1:numel(cases)
    case_def = cases(i_case);
    if ~local_keep_value(case_def.id, options, 'case_ids')
        continue
    end

    for i_mesh = 1:numel(mesh_variants)
        mesh_variant = mesh_variants(i_mesh);
        if ~local_keep_value(mesh_variant.level, options, 'levels')
            continue
        end

        for i_tol = 1:numel(tolerances)
            tol = tolerances(i_tol);
            if ~local_keep_value(tol, options, 'tolerances')
                continue
            end

            cfg = case_def;
            cfg.level = mesh_variant.level;
            cfg.solver_type = "DCG_HYPRE_BOOMERAMG";
            cfg.linear_solver_tolerance = tol;
            cfg.linear_solver_tolerance_newton = tol;
            cfg.linear_solver_tolerance_qn1 = tol;
            cfg.linear_solver_tolerance_qn2 = tol;
            cfg.use_damping = force_damping;

            if ~isempty(force_newton)
                cfg.run_newton = force_newton;
            end
            if ~isempty(force_qn1)
                cfg.run_qn1 = force_qn1;
            end
            if ~isempty(force_qn2)
                cfg.run_qn2 = force_qn2;
            end

            counter = counter + 1;
            plan(counter).key = sprintf('%s__%s__tol_%s', ...
                char(case_def.id), char(mesh_variant.label), local_number_token(tol));
            plan(counter).case_id = case_def.id;
            plan(counter).case_name = case_def.name;
            plan(counter).mesh_label = mesh_variant.label;
            plan(counter).mesh_display_name = mesh_variant.display_name;
            plan(counter).level = mesh_variant.level;
            plan(counter).solver_type = "DCG_HYPRE_BOOMERAMG";
            plan(counter).solver_label = "DCG-HYPRE";
            plan(counter).solver_method = "DCG";
            plan(counter).preconditioner = "HYPRE";
            plan(counter).linear_solver_tolerance = tol;
            plan(counter).cfg = cfg;
        end
    end
end

function value = local_option_bool(options, field_name, default_value)
if ~isfield(options, field_name)
    value = default_value;
    return
end

value = options.(field_name);
if isempty(value)
    value = default_value;
    return
end

value = logical(value);
end

if isfield(options, 'max_cases') && options.max_cases < numel(plan)
    plan = plan(1:options.max_cases);
end
end

function keep = local_keep_value(value, options, field_name)
if ~isfield(options, field_name)
    keep = true;
    return
end

requested = options.(field_name);
if isstring(value) || ischar(value)
    requested = string(requested);
    keep = any(string(value) == requested);
else
    keep = any(abs(value - requested) < eps(max(1, abs(value))));
end
end

function token = local_number_token(value)
token = regexprep(sprintf('%.0e', value), '\+', '');
token = strrep(token, '.', 'p');
end

function idx = local_find_result_index(results, key)
idx = 0;
for i = 1:numel(results)
    if strcmp(string(results(i).key), string(key))
        idx = i;
        return
    end
end
end

function entry = local_result_entry(plan_entry, case_result, error_text, succeeded)
entry = struct();
entry.key = plan_entry.key;
entry.case_id = plan_entry.case_id;
entry.case_name = plan_entry.case_name;
entry.mesh_label = plan_entry.mesh_label;
entry.mesh_display_name = plan_entry.mesh_display_name;
entry.level = plan_entry.level;
entry.solver_type = plan_entry.solver_type;
entry.solver_label = plan_entry.solver_label;
entry.solver_method = plan_entry.solver_method;
entry.preconditioner = plan_entry.preconditioner;
entry.linear_solver_tolerance = plan_entry.linear_solver_tolerance;
entry.succeeded = succeeded;
entry.error_text = error_text;
entry.n_unknowns = case_result.n_unknowns;
entry.n_elements = case_result.n_elements;
entry.newton = case_result.newton;
entry.qn1 = case_result.qn1;
entry.qn2 = case_result.qn2;
end

function result = local_failure_result()
result = struct( ...
    'n_unknowns', NaN, ...
    'n_elements', NaN, ...
    'newton', local_unused_solver_result(), ...
    'qn1', local_unused_solver_result(), ...
    'qn2', local_unused_solver_result());
end

function result = local_unused_solver_result()
result = struct( ...
    'iterations', NaN, ...
    'linear_iterations', NaN, ...
    'runtime_seconds', NaN, ...
    'final_criterion', NaN, ...
    'converged', false, ...
    'timed_out', false);
end

function local_write_markdown(markdown_path, results, plan, cases, tolerances, mesh_variants)
fid = fopen(markdown_path, 'w');
if fid < 0
    error('Unable to open markdown output: %s', markdown_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Paper Tolerance Sweep for `magnetic_potential_problem_2D`\n\n');
fprintf(fid, 'This report tracks Newton and quasi-Newton performance with DCG + HYPRE on the magnetic benchmark.\n\n');
tol_labels = local_format_tolerances(tolerances);
fprintf(fid, '- tolerances: `%s`\n', strjoin(tol_labels, '`, `'));
fprintf(fid, '- completed records: `%d`\n', numel(results));
fprintf(fid, '- current sweep plan size: `%d`\n', numel(plan));
fprintf(fid, '- generated from:\n');
fprintf(fid, '  - [run_paper_case.m](scripts/run_paper_case.m)\n');
fprintf(fid, '  - [run_paper_tolerance_sweep.m](scripts/run_paper_tolerance_sweep.m)\n');
fprintf(fid, '  - [magnetic_problem.m](magnetic_problem.m)\n');
fprintf(fid, '  - [paper_tolerance_sweep_results.mat](results/paper_tolerance_sweep_results.mat)\n\n');

for i_case = 1:numel(cases)
    case_entry = cases(i_case);
    fprintf(fid, '## %s\n\n', case_entry.name);
    fprintf(fid, '- case id: `%s`\n', case_entry.id);
    fprintf(fid, '- geometry: `%g x %g`\n', case_entry.size_x, case_entry.size_y);
    fprintf(fid, '- delta: `%g`, alpha: `%g`, beta: `%g`, rho: `%g`, kappa: `%g`\n\n', ...
        case_entry.delta, case_entry.alpha, case_entry.beta, case_entry.rho, case_entry.kappa);

    for i_mesh = 1:numel(mesh_variants)
        mesh = mesh_variants(i_mesh);
        fprintf(fid, '### %s\n\n', mesh.display_name);
        fprintf(fid, '| linear tolerance | Newton | qN1 | qN2 |\n');
        fprintf(fid, '|---|---|---|---|\n');
        for i_tol = 1:numel(tolerances)
            tol = tolerances(i_tol);
            entry = local_find_entry(results, case_entry.id, mesh.label, tol);
            if isempty(entry)
                fprintf(fid, '| `%s` | `-` | `-` | `-` |\n', local_format_tolerance(tol));
                continue
            end
            n_entry = local_format_solver_cell(entry.newton);
            q1_entry = local_format_solver_cell(entry.qn1);
            q2_entry = local_format_solver_cell(entry.qn2);
            fprintf(fid, '| `%s` | %s | %s | %s |\n', local_format_tolerance(tol), n_entry, q1_entry, q2_entry);
        end
        fprintf(fid, '\n');
    end
end

fprintf(fid, 'Cells are formatted as `nonlinear_iterations / cumulative_linear_iterations / runtime [s]`.\n');
fprintf(fid, 'A row prefixed with `~` did not satisfy the stored convergence criterion when stopping.\n');
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

function text = local_format_solver_cell(solver_result)
if isempty(solver_result) || ~isfield(solver_result, 'iterations') || ~isfinite(solver_result.iterations)
    text = '`-`';
    return
end

if ~isfinite(solver_result.runtime_seconds)
    text = '`-`';
    return
end

if solver_result.timed_out
    text = sprintf('`%d / %d / >%g`', ...
        solver_result.iterations, round(solver_result.linear_iterations), solver_result.runtime_seconds);
    return
end

base = sprintf('%d / %d / %.3g', ...
    solver_result.iterations, round(solver_result.linear_iterations), solver_result.runtime_seconds);
if solver_result.converged
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
