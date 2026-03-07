function run_damping_level4_sweep(options)
%RUN_DAMPING_LEVEL4_SWEEP Run damped magnetic benchmark sweeps on level-4 mesh.
%
% The sweep is built around paper-relevant parameters:
% - delta: heterogeneity region half-width ratio (0 means no heterogeneity)
% - rho:   right-hand-side scaling, stronger forcing = larger rho
% - kappa: quasi-Newton 2 criticality parameter
%
% Fixed settings for this experiment:
% - level 4 mesh -> 51200 elements / 25281 unknowns (for size_x=size_y=1)
% - inner linear tolerances: Newton 1e-4, qN1/qN2 1e-1
% - solver: DCG + HYPRE on scalar BoomerAMG (num_functions=1)
% - damping enabled for all nonlinear methods

if nargin < 1 || isempty(options)
    options = struct();
end

paths = mp2d_paths();
results_path = fullfile(paths.results, 'damping_level4_results.mat');
markdown_path = fullfile(paths.root, 'DAMPING_LEVEL4_EXPERIMENTS.md');

deltas = local_option_numeric(options, 'deltas', [0, 1, 2]);
rhos = local_option_numeric(options, 'rhos', [40, 80, 160]);
kappas = local_option_numeric(options, 'kappas', [10, 25, 50]);
force_rerun = isfield(options, 'force_rerun') && logical(options.force_rerun);
use_damping = local_option_bool(options, 'use_damping', true);

cases = local_build_cases(deltas, rhos, kappas, use_damping);
results = struct([]);

for i_case = 1:numel(cases)
    cfg = cases(i_case).cfg;
    case_key = cases(i_case).key;

    if ~force_rerun && local_existing_index(results, case_key) > 0
        fprintf('[%d/%d] skip %s\n', i_case, numel(cases), case_key);
        continue
    end

    fprintf('[%d/%d] run %s\n', i_case, numel(cases), case_key);
    try
        case_result = run_paper_case(cfg);
        entry = local_build_entry(case_key, cfg, case_result, true, "");
    catch ME
        entry = local_build_entry(case_key, cfg, local_failure_result(), false, ...
            string(getReport(ME, 'extended', 'hyperlinks', 'off')));
    end
    results = local_store_result(results, entry);
end

save(results_path, 'results', '-v7.3');
local_write_markdown(markdown_path, results);

end

function cases = local_build_cases(deltas, rhos, kappas, use_damping)
cases = struct([]);
counter = 0;
common = struct();
common.level = 4;
common.elem_type = "P1";
common.size_x = 1;
common.size_y = 1;
common.alpha = 3e-4;
common.beta = 16000;
common.solver_type = "DCG_HYPRE_BOOMERAMG";
common.linear_solver_maxit = 1000;
common.deflation_basis_tolerance = 1e-10;
common.nonlinear_method_max_runtime_seconds = 120;
common.run_newton = true;
common.run_qn1 = true;
common.run_qn2 = true;
common.run_postprocess = false;
common.boomeramg_opts = struct('threads', 16, 'print_level', 0, 'use_as_preconditioner', true, 'num_functions', 1);
common.linear_solver_tolerance_newton = 1e-4;
common.linear_solver_tolerance_qn1 = 1e-1;
common.linear_solver_tolerance_qn2 = 1e-1;
common.use_damping = use_damping;

for i_delta = 1:numel(deltas)
    for i_rho = 1:numel(rhos)
        for i_kappa = 1:numel(kappas)
            counter = counter + 1;
            cfg = common;
            cfg.delta = deltas(i_delta);
            cfg.rho = rhos(i_rho);
            cfg.kappa = kappas(i_kappa);
            cfg.case_id = sprintf('case_d_%d_r_%g_k_%g', cfg.delta, cfg.rho, cfg.kappa);
            cfg.case_name = sprintf('delta=%d, rho=%g, kappa=%g', cfg.delta, cfg.rho, cfg.kappa);

            cases(counter).key = sprintf('%s__d_%d__rho_%g__k_%g', ...
                cfg.case_id, cfg.delta, cfg.rho, cfg.kappa);
            cases(counter).cfg = cfg;
            cases(counter).delta = cfg.delta;
            cases(counter).rho = cfg.rho;
            cases(counter).kappa = cfg.kappa;
            cases(counter).case_id = cfg.case_id;
            cases(counter).case_name = cfg.case_name;
        end
    end
end

end

function entry = local_build_entry(key, cfg, case_result, succeeded, error_text)
entry = struct();
entry.key = key;
entry.case_id = string(cfg.case_id);
entry.case_name = string(cfg.case_name);
entry.mesh_level = cfg.level;
entry.mesh_label = "level_4";
entry.solver_type = string(cfg.solver_type);
entry.solver_label = "DCG + HYPRE";
entry.preconditioner = "HYPRE";
entry.linear_solver_tolerance_newton = cfg.linear_solver_tolerance_newton;
entry.linear_solver_tolerance_quasi = cfg.linear_solver_tolerance_qn1;
entry.nonlinear_method_max_runtime_seconds = cfg.nonlinear_method_max_runtime_seconds;
entry.use_damping = cfg.use_damping;
entry.succeeded = succeeded;
entry.error_text = error_text;
entry.delta = cfg.delta;
entry.rho = cfg.rho;
entry.kappa = cfg.kappa;
entry.beta = cfg.beta;
entry.alpha = cfg.alpha;
entry.n_elements = case_result.n_elements;
entry.n_unknowns = case_result.n_unknowns;
entry.newton = case_result.newton;
entry.qn1 = case_result.qn1;
entry.qn2 = case_result.qn2;
entry.boomeramg_opts = cfg.boomeramg_opts;
end

function result = local_failure_result()
stub = local_unused_solver_result();
result = struct();
result.n_elements = NaN;
result.n_unknowns = NaN;
result.newton = stub;
result.qn1 = stub;
result.qn2 = stub;
end

function solver = local_unused_solver_result()
solver = struct( ...
    'iterations', NaN, ...
    'linear_iterations', NaN, ...
    'runtime_seconds', NaN, ...
    'final_criterion', NaN, ...
    'converged', false, ...
    'timed_out', false);
end

function idx = local_existing_index(results, key)
idx = 0;
for i = 1:numel(results)
    if string(results(i).key) == string(key)
        idx = i;
        return
    end
end
end

function results = local_store_result(results, entry)
if isempty(results)
    results = entry;
    return
end

idx = local_existing_index(results, entry.key);
if idx > 0
    results(idx) = entry;
else
    results(end + 1) = entry;
end
end

function value = local_option_numeric(options, field_name, default_value)
if ~isfield(options, field_name)
    value = default_value;
    return
end

value = options.(field_name);
if isempty(value)
    value = default_value;
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

function local_write_markdown(markdown_path, results)
fid = fopen(markdown_path, 'w');
if fid < 0
    error('Unable to open markdown output file %s', markdown_path);
end
cleanup = onCleanup(@() fclose(fid));

cases_count = numel(results);
if cases_count == 0
    fprintf(fid, '# Damping Sweep (Level 4) for `magnetic_potential_problem_2D`\n\n');
    fprintf(fid, 'No cases were executed.\n');
    return
end

first = results(1);
if first.use_damping
    damping_status = 'enabled';
else
    damping_status = 'disabled';
end

fprintf(fid, '# Damping Sweep (Level 4) for `magnetic_potential_problem_2D`\n\n');
fprintf(fid, 'Scope:\n');
fprintf(fid, '- mesh level: `4`\n');
fprintf(fid, '- measured elements: `%d`\n', first.n_elements);
fprintf(fid, '- measured unknowns: `%d`\n', first.n_unknowns);
fprintf(fid, '- linear tolerance policy: `N = 1e-4`, `qN1/qN2 = 1e-1`\n');
fprintf(fid, '- nonlinear mode: damping %s for Newton, qN1, qN2\n', damping_status);
fprintf(fid, '- solver stack: `DCG + HYPRE BoomerAMG` with scalar BoomerAMG setting (`num_functions = 1`).\n\n');
fprintf(fid, 'Cell format: `iterations / runtime [s]`, prefixed by `~` for non-converged results.\n\n');

fprintf(fid, '| case | `δ` | `ρ` | `κ` | `N` | `qN1` | `qN2` |\n');
fprintf(fid, '|---|---|---|---|---|---|---|\n');
for i = 1:cases_count
    entry = results(i);
    n_cell = local_format_solver_cell(entry.newton);
    q1_cell = local_format_solver_cell(entry.qn1);
    q2_cell = local_format_solver_cell(entry.qn2);
    fprintf(fid, '| %s | %g | %g | %g | %s | %s | %s |\n', ...
        entry.case_id, entry.delta, entry.rho, entry.kappa, ...
        n_cell, q1_cell, q2_cell);
end

fprintf(fid, '\n');

local_write_best_hint(fid, results);
end

function local_write_best_hint(fid, results)
best_runtime = inf;
best_idx = 0;

for i = 1:numel(results)
    if ~isfinite(results(i).newton.runtime_seconds)
        continue
    end
    if results(i).newton.runtime_seconds < best_runtime
        best_runtime = results(i).newton.runtime_seconds;
        best_idx = i;
    end
end

if best_idx > 0
    best = results(best_idx);
    fprintf(fid, 'Best `N` runtime in this sweep: case `%s` (`%g / %g / %g`) in `%gs`.\n', ...
        best.case_id, best.newton.iterations, best.newton.linear_iterations, best.newton.final_criterion, best.newton.runtime_seconds);
end

any_qn_better = false;
for i = 1:numel(results)
    candidate = results(i);
    runtimes = [candidate.newton.runtime_seconds, candidate.qn1.runtime_seconds, candidate.qn2.runtime_seconds];
    runtimes = runtimes(isfinite(runtimes));
    if isempty(runtimes)
        continue
    end
    best_runtime = min(runtimes);
    if isfinite(candidate.qn1.runtime_seconds) && candidate.qn1.runtime_seconds == best_runtime
        any_qn_better = true;
        break
    end
    if isfinite(candidate.qn2.runtime_seconds) && candidate.qn2.runtime_seconds == best_runtime
        any_qn_better = true;
        break
    end
end

if any_qn_better
    fprintf(fid, 'There are settings in this sweep where at least one quasi-Newton runtime is the fastest.\n');
else
    fprintf(fid, 'Within this sweep, quasi-Newton methods did not beat Newton in runtime.\n');
end
end

function text = local_format_solver_cell(solver_result)
if isempty(solver_result) || ~isfield(solver_result, 'iterations') || ~isfinite(solver_result.iterations)
    text = '`-`';
    return
end
if isfield(solver_result, 'timed_out') && solver_result.timed_out
    text = sprintf('`~%d / >%g`', solver_result.iterations, solver_result.runtime_seconds);
    return
end

text = sprintf('`%d / %.3g`', solver_result.iterations, solver_result.runtime_seconds);
if isfield(solver_result, 'converged') && solver_result.converged
    text = sprintf('%s', text);
else
    text = sprintf('`~%d / %.3g`', solver_result.iterations, solver_result.runtime_seconds);
end
end
