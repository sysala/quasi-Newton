function run_selected_tolerance_mesh_method(example_id, gamma, method_name, options)
%RUN_SELECTED_TOLERANCE_MESH_METHOD Run or mark one ~200k selected-mesh method case.
%
% This helper is intended for shell-level orchestration with an outer
% `timeout` wrapper so HYPRE setup stalls can still be converted into
% `>300s` report entries.

if nargin < 4 || isempty(options)
    options = struct();
end

paths = ns3d_paths();
results_path = fullfile(paths.results, 'selected_tolerance_mesh_results_dcg_hypre.mat');
if isfield(options, 'results_path') && ~isempty(options.results_path)
    results_path = char(string(options.results_path));
end

force_rerun = isfield(options, 'force_rerun') && logical(options.force_rerun);
mark_timeout = isfield(options, 'mark_timeout') && logical(options.mark_timeout);

plan_entry = local_plan_entry(example_id, gamma);
[entry_key, method_field, tol] = local_entry_key(plan_entry.base_key, method_name);

results = struct([]);
if isfile(results_path)
    loaded = load(results_path, 'results');
    if isfield(loaded, 'results')
        results = loaded.results;
    end
end

existing_idx = local_find_result_index(results, entry_key);
if ~force_rerun && existing_idx > 0 && ...
        local_method_recorded(results(existing_idx).(method_field))
    fprintf('skip %s\n', entry_key);
    return
end

if mark_timeout
    method_result = local_timeout_result();
    case_result = struct( ...
        'n_unknowns', NaN, ...
        'n_elements', NaN, ...
        'newton', local_unused_solver_result(), ...
        'qn1', local_unused_solver_result(), ...
        'qn2', local_unused_solver_result());
else
    cfg = plan_entry.cfg;
    cfg.run_newton = strcmp(method_field, 'newton');
    cfg.run_qn1 = strcmp(method_field, 'qn1');
    cfg.run_qn2 = strcmp(method_field, 'qn2');
    case_result = run_paper_case(cfg);
    method_result = case_result.(method_field);
end

entry = local_make_entry(plan_entry, entry_key, tol);
if existing_idx > 0
    entry = local_merge_entry(entry, results(existing_idx));
end

entry.succeeded = true;
entry.error_text = "";
entry.n_unknowns = local_pick_value(entry.n_unknowns, case_result.n_unknowns);
entry.n_elements = local_pick_value(entry.n_elements, case_result.n_elements);
entry.(method_field) = method_result;

if existing_idx > 0
    results(existing_idx) = entry;
elseif isempty(results)
    results = entry;
else
    results(end + 1) = entry; %#ok<AGROW>
end

save(results_path, 'results', '-v7.3');
fprintf('saved %s\n', entry_key);
end

function plan_entry = local_plan_entry(example_id, gamma)
mesh_variant = struct( ...
    'label', "roughly_200k_unknowns_mesh", ...
    'display_name', "~200k unknowns mesh", ...
    'density', 19);
solver_spec = struct( ...
    'solver_type', "DCG_HYPRE_BOOMERAMG", ...
    'label', "DCG + HYPRE", ...
    'method', "DCG", ...
    'preconditioner', "HYPRE");

example = local_example(example_id);
plan_entry = struct();
plan_entry.base_key = sprintf('%s__%s__%s__g_%s', ...
    char(example.id), char(mesh_variant.label), char(solver_spec.solver_type), local_number_token(gamma));
plan_entry.example_id = example.id;
plan_entry.example_name = example.name;
plan_entry.mesh_label = mesh_variant.label;
plan_entry.mesh_display_name = mesh_variant.display_name;
plan_entry.density = mesh_variant.density;
plan_entry.solver_type = solver_spec.solver_type;
plan_entry.solver_label = solver_spec.label;
plan_entry.solver_method = solver_spec.method;
plan_entry.preconditioner = solver_spec.preconditioner;
plan_entry.gamma = gamma;
plan_entry.boundary = example.boundary;
plan_entry.cfg = example;
plan_entry.cfg.density = mesh_variant.density;
plan_entry.cfg.solver_type = solver_spec.solver_type;
plan_entry.cfg.gamma = gamma;
end

function example = local_example(example_id)
common = struct( ...
    'size_z', 50, ...
    'mu_0', 1, ...
    'mu_infty', 1e-3, ...
    'lambda', 10, ...
    'p', 1.1, ...
    'u_z', 5, ...
    'linear_solver_maxit', 1000, ...
    'deflation_basis_tolerance', 1e-10, ...
    'boomeramg_opts', struct('threads', 16, 'print_level', 0, 'use_as_preconditioner', true), ...
    'linear_solver_tolerance_newton', 1e-4, ...
    'linear_solver_tolerance_quasi', 1e-1, ...
    'nonlinear_method_max_runtime_seconds', 300, ...
    'run_postprocess', false);

switch lower(char(string(example_id)))
    case 'example1'
        example = common;
        example.id = "example1";
        example.name = "Example 1: Prismatic Bar";
        example.size_xy_0 = 5;
        example.size_xy_L = 5;
        example.boundary = "BC1";
    case 'example2'
        example = common;
        example.id = "example2";
        example.name = "Example 2: Tapered Prismatic Bar";
        example.size_xy_0 = 6;
        example.size_xy_L = 5;
        example.boundary = "BC1";
    otherwise
        error('Unknown example_id: %s', char(string(example_id)));
end
end

function [entry_key, method_field, tol] = local_entry_key(base_key, method_name)
method_name = lower(char(string(method_name)));
switch method_name
    case {'n', 'newton'}
        method_field = 'newton';
        tol = 1e-4;
    case {'qn1'}
        method_field = 'qn1';
        tol = 1e-1;
    case {'qn2'}
        method_field = 'qn2';
        tol = 1e-1;
    otherwise
        error('Unknown method_name: %s', method_name);
end

entry_key = sprintf('%s__tol_%s', base_key, local_number_token(tol));
end

function entry = local_make_entry(plan_entry, entry_key, tol)
entry = struct( ...
    'key', string(entry_key), ...
    'example_id', plan_entry.example_id, ...
    'example_name', plan_entry.example_name, ...
    'mesh_label', plan_entry.mesh_label, ...
    'mesh_display_name', plan_entry.mesh_display_name, ...
    'density', plan_entry.density, ...
    'solver_type', plan_entry.solver_type, ...
    'solver_label', plan_entry.solver_label, ...
    'solver_method', plan_entry.solver_method, ...
    'preconditioner', plan_entry.preconditioner, ...
    'linear_solver_tolerance', tol, ...
    'gamma', plan_entry.gamma, ...
    'boundary', plan_entry.boundary, ...
    'succeeded', false, ...
    'error_text', "", ...
    'n_unknowns', NaN, ...
    'n_elements', NaN, ...
    'newton', local_unused_solver_result(), ...
    'qn1', local_unused_solver_result(), ...
    'qn2', local_unused_solver_result());
end

function merged = local_merge_entry(base_entry, existing_entry)
merged = base_entry;
fields = fieldnames(base_entry);
for i = 1:numel(fields)
    field_name = fields{i};
    if isstruct(base_entry.(field_name))
        if isfield(existing_entry, field_name)
            merged.(field_name) = existing_entry.(field_name);
        end
    elseif isfield(existing_entry, field_name)
        merged.(field_name) = existing_entry.(field_name);
    end
end
end

function recorded = local_method_recorded(method_result)
recorded = (isfield(method_result, 'timed_out') && logical(method_result.timed_out)) || ...
    (isfield(method_result, 'runtime_seconds') && isfinite(method_result.runtime_seconds));
end

function value = local_pick_value(current_value, new_value)
if isfinite(current_value)
    value = current_value;
elseif isfinite(new_value)
    value = new_value;
else
    value = current_value;
end
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

function solver_result = local_unused_solver_result()
solver_result = struct( ...
    'iterations', NaN, ...
    'linear_iterations', NaN, ...
    'runtime_seconds', NaN, ...
    'final_criterion', NaN, ...
    'converged', false, ...
    'timed_out', false);
end

function solver_result = local_timeout_result()
solver_result = struct( ...
    'iterations', NaN, ...
    'linear_iterations', NaN, ...
    'runtime_seconds', 300, ...
    'final_criterion', NaN, ...
    'converged', false, ...
    'timed_out', true);
end

function token = local_number_token(value)
token = regexprep(sprintf('%.0e', value), '\+', '');
token = strrep(token, '.', 'p');
end
