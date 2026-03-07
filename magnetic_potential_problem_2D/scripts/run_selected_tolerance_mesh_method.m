function run_selected_tolerance_mesh_method(case_id, method_name, options)
%RUN_SELECTED_TOLERANCE_MESH_METHOD Run one selected-mesh method case.
%
% This helper supports shell-level orchestration with per-method timeouts.
%
% Optional input:
%   options - struct with:
%       force_rerun : logical scalar
%       mark_timeout : logical scalar
%       use_damping : logical scalar

if nargin < 3 || isempty(options)
    options = struct();
end

paths = mp2d_paths();
results_path = fullfile(paths.results, 'selected_tolerance_mesh_results_dcg_hypre.mat');
if isfield(options, 'results_path') && ~isempty(options.results_path)
    results_path = char(string(options.results_path));
end

force_rerun = isfield(options, 'force_rerun') && logical(options.force_rerun);
mark_timeout = isfield(options, 'mark_timeout') && logical(options.mark_timeout);
use_damping = local_option_bool(options, 'use_damping', false);

local_clear_hypre_instances();
plan_entry = local_plan_entry(case_id);
[entry_key, method_field, tol] = local_entry_key(plan_entry.base_key, method_name);

results = struct([]);
if isfile(results_path)
    loaded = load(results_path, 'results');
    if isfield(loaded, 'results')
        results = loaded.results;
    end
end

existing_idx = local_find_result_index(results, entry_key);
if ~force_rerun && existing_idx > 0 && local_method_recorded(results(existing_idx).(method_field))
    fprintf('skip %s\n', entry_key);
    return
end

if mark_timeout
    method_result = local_timeout_result();
    case_result = local_unused_case_result();
else
    cfg = plan_entry.cfg;
    cfg.run_newton = strcmp(method_field, 'newton');
    cfg.run_qn1 = strcmp(method_field, 'qn1');
    cfg.run_qn2 = strcmp(method_field, 'qn2');
    cfg.use_damping = use_damping;
    case_result = run_paper_case(cfg);
    method_result = case_result.(method_field);
end
local_clear_hypre_instances();

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

function plan_entry = local_plan_entry(case_id)
mesh_variant = struct( ...
    'label', "level_4", ...
    'display_name', "Level 4 mesh", ...
    'level', 4);
solver_spec = struct( ...
    'solver_type', "DCG_HYPRE_BOOMERAMG", ...
    'label', "DCG + HYPRE", ...
    'method', "DCG", ...
    'preconditioner', "HYPRE");

case_entry = local_case(case_id);

plan_entry = struct();
plan_entry.base_key = sprintf('%s__%s__%s', ...
    char(case_entry.id), char(mesh_variant.label), char(solver_spec.solver_type));
plan_entry.case_id = case_entry.id;
plan_entry.case_name = case_entry.name;
plan_entry.mesh_label = mesh_variant.label;
plan_entry.mesh_display_name = mesh_variant.display_name;
plan_entry.level = mesh_variant.level;
plan_entry.solver_type = solver_spec.solver_type;
plan_entry.solver_label = solver_spec.label;
plan_entry.solver_method = solver_spec.method;
plan_entry.preconditioner = solver_spec.preconditioner;
plan_entry.linear_solver_tolerance = NaN;
plan_entry.boundary = case_entry.boundary;
plan_entry.cfg = case_entry;
plan_entry.cfg.level = mesh_variant.level;
plan_entry.cfg.solver_type = solver_spec.solver_type;
plan_entry.cfg.boundary = case_entry.boundary;
end

function case_entry = local_case(case_id)
common = struct( ...
    'elem_type', "P1", ...
    'level', 4, ...
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
    'linear_solver_tolerance_newton', 1e-4, ...
    'linear_solver_tolerance_quasi', 1e-1, ...
    'nonlinear_method_max_runtime_seconds', 120);
    common.use_damping = false;

case_id = string(case_id);

switch lower(char(case_id))
    case 'case_delta_1'
        case_entry = common;
        case_entry.id = "case_delta_1";
        case_entry.name = "Inclusion contrast δ=1";
        case_entry.delta = 1;
        case_entry.boundary = "homogeneous_dirichlet";
    case 'case_delta_2'
        case_entry = common;
        case_entry.id = "case_delta_2";
        case_entry.name = "Inclusion contrast δ=2";
        case_entry.delta = 2;
        case_entry.boundary = "homogeneous_dirichlet";
    otherwise
        error('Unknown case_id: %s', char(case_id));
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
    'case_id', plan_entry.case_id, ...
    'case_name', plan_entry.case_name, ...
    'mesh_label', plan_entry.mesh_label, ...
    'mesh_display_name', plan_entry.mesh_display_name, ...
    'level', plan_entry.level, ...
    'solver_type', plan_entry.solver_type, ...
    'solver_label', plan_entry.solver_label, ...
    'solver_method', plan_entry.solver_method, ...
    'preconditioner', plan_entry.preconditioner, ...
    'linear_solver_tolerance', tol, ...
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
    'runtime_seconds', 120, ...
    'final_criterion', NaN, ...
    'converged', false, ...
    'timed_out', true);
end

function case_result = local_unused_case_result()
case_result = struct( ...
    'n_unknowns', NaN, ...
    'n_elements', NaN, ...
    'newton', local_unused_solver_result(), ...
    'qn1', local_unused_solver_result(), ...
    'qn2', local_unused_solver_result());
end

function token = local_number_token(value)
token = regexprep(sprintf('%.0e', value), '\+', '');
token = strrep(token, '.', 'p');
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

function local_clear_hypre_instances()
try
    LINEAR_SOLVERS.hypre_boomeramg_clear();
catch
end
end
