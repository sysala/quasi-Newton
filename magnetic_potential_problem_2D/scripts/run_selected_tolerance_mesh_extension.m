function run_selected_tolerance_mesh_extension(options)
%RUN_SELECTED_TOLERANCE_MESH_EXTENSION Run the selected-mesh method replication.
%
% This runner executes the selected-tolerance policy:
%   N   -> 1e-4
%   qN1 -> 1e-1
%   qN2 -> 1e-1
% on the selected mesh level.
%
% Optional input:
%   options - struct with optional filters:
%       case_ids  : cellstr
%       force_rerun : logical scalar
%       use_damping : logical scalar
%       max_cases : positive integer
%       run_single_pending : logical scalar

if nargin < 1 || isempty(options)
    options = struct();
end

paths = mp2d_paths();
results_path = fullfile(paths.results, 'selected_tolerance_mesh_results_dcg_hypre.mat');

selected_mesh = struct( ...
    'label', "level_4", ...
    'display_name', "Level 4 mesh", ...
    'level', 4);
solver_spec = struct( ...
    'solver_type', "DCG_HYPRE_BOOMERAMG", ...
    'label', "DCG + HYPRE", ...
    'method', "DCG", ...
    'preconditioner', "HYPRE");

cases = local_cases();
plan = local_build_plan(cases, selected_mesh, solver_spec, options);

results = struct([]);
if isfile(results_path)
    loaded = load(results_path, 'results');
    if isfield(loaded, 'results')
        results = loaded.results;
    end
end

force_rerun = isfield(options, 'force_rerun') && logical(options.force_rerun);
run_single_pending = isfield(options, 'run_single_pending') && logical(options.run_single_pending);
force_damping = local_option_bool(options, 'use_damping', false);

for i = 1:numel(plan)
    plan_entry = plan(i);
    base_key = plan_entry.base_key;
    if ~force_rerun && local_has_complete_entries(results, base_key)
        fprintf('[%d/%d] skip %s\n', i, numel(plan), base_key);
        continue
    end

    fprintf('[%d/%d] run %s\n', i, numel(plan), base_key);
    try
        local_clear_hypre_instances();
        case_result = run_paper_case(plan_entry.cfg);
        entries = local_result_entries(plan_entry, case_result, "", true);
    catch ME
        entries = local_failure_entries(plan_entry, ...
            string(getReport(ME, 'extended', 'hyperlinks', 'off')));
    end
    local_clear_hypre_instances();

    for i_entry = 1:numel(entries)
        idx = local_find_result_index(results, entries(i_entry).key);
        if idx > 0
            results(idx) = entries(i_entry);
        elseif isempty(results)
            results = entries(i_entry);
        else
            results(end + 1) = entries(i_entry); %#ok<AGROW>
        end
    end
    save(results_path, 'results', 'plan', '-v7.3');

    if run_single_pending
        return
    end

end
end

function cases = local_cases()
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

cases = repmat(common, 1, 2);
cases(1).id = "case_delta_1";
cases(1).name = "Inclusion contrast δ=1";
cases(1).delta = 1;
cases(1).boundary = "homogeneous_dirichlet";

cases(2).id = "case_delta_2";
cases(2).name = "Inclusion contrast δ=2";
cases(2).delta = 2;
cases(2).boundary = "homogeneous_dirichlet";
end

function plan = local_build_plan(cases, selected_mesh, solver_spec, options)
plan = struct([]);
counter = 0;

for i_case = 1:numel(cases)
    case_entry = cases(i_case);
    if ~local_keep_value(case_entry.id, options, 'case_ids')
        continue
    end

    counter = counter + 1;
    plan(counter).base_key = sprintf('%s__%s__%s', ...
        char(case_entry.id), char(selected_mesh.label), char(solver_spec.solver_type)); %#ok<AGROW>
    plan(counter).case_id = case_entry.id;
    plan(counter).case_name = case_entry.name;
    plan(counter).mesh_label = selected_mesh.label;
    plan(counter).mesh_display_name = selected_mesh.display_name;
    plan(counter).level = selected_mesh.level;
    plan(counter).solver_type = solver_spec.solver_type;
    plan(counter).solver_label = solver_spec.label;
    plan(counter).solver_method = solver_spec.method;
    plan(counter).preconditioner = solver_spec.preconditioner;
    plan(counter).boundary = case_entry.boundary;
    plan(counter).cfg = case_entry;
    plan(counter).cfg.level = selected_mesh.level;
    plan(counter).cfg.solver_type = solver_spec.solver_type;
    plan(counter).cfg.use_damping = force_damping;
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

function complete = local_has_complete_entries(results, base_key)
complete = local_find_result_index(results, base_key + "__tol_1e-04") > 0 && ...
    local_find_result_index(results, base_key + "__tol_1e-01") > 0;
end

function idx = local_find_result_index(results, key)
idx = 0;
for i = 1:numel(results)
    if strcmp(results(i).key, key)
        idx = i;
        return
    end
end
end

function entries = local_result_entries(plan_entry, case_result, error_text, succeeded)
tol_newton = 1e-4;
tol_quasi = 1e-1;
template = local_result_entry_template();

entries = repmat(template, 1, 2);
entries(1) = local_fill_entry(plan_entry, case_result, error_text, succeeded, tol_newton, ...
    case_result.newton, local_unused_solver_result(), local_unused_solver_result());
entries(2) = local_fill_entry(plan_entry, case_result, error_text, succeeded, tol_quasi, ...
    local_unused_solver_result(), case_result.qn1, case_result.qn2);
end

function entry = local_fill_entry(plan_entry, case_result, error_text, succeeded, tol, newton_result, qn1_result, qn2_result)
entry = struct();
entry.key = sprintf('%s__tol_%s', char(plan_entry.base_key), local_number_token(tol));
entry.case_id = plan_entry.case_id;
entry.case_name = plan_entry.case_name;
entry.mesh_label = plan_entry.mesh_label;
entry.mesh_display_name = plan_entry.mesh_display_name;
entry.level = plan_entry.level;
entry.solver_type = plan_entry.solver_type;
entry.solver_label = plan_entry.solver_label;
entry.solver_method = plan_entry.solver_method;
entry.preconditioner = plan_entry.preconditioner;
entry.linear_solver_tolerance = tol;
entry.boundary = plan_entry.boundary;
entry.succeeded = succeeded;
entry.error_text = error_text;
entry.n_unknowns = case_result.n_unknowns;
entry.n_elements = case_result.n_elements;
entry.newton = newton_result;
entry.qn1 = qn1_result;
entry.qn2 = qn2_result;
end

function entries = local_failure_entries(plan_entry, error_text)
failure_case = struct('n_unknowns', NaN, 'n_elements', NaN, ...
    'newton', local_unused_solver_result(), ...
    'qn1', local_unused_solver_result(), ...
    'qn2', local_unused_solver_result());
entries = local_result_entries(plan_entry, failure_case, error_text, false);
end

function case_result = local_unused_case_result()
case_result = struct( ...
    'n_unknowns', NaN, ...
    'n_elements', NaN, ...
    'newton', local_unused_solver_result(), ...
    'qn1', local_unused_solver_result(), ...
    'qn2', local_unused_solver_result());
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

function entry = local_result_entry_template()
entry = struct( ...
    'key', "", ...
    'case_id', "", ...
    'case_name', "", ...
    'mesh_label', "", ...
    'mesh_display_name', "", ...
    'level', NaN, ...
    'solver_type', "", ...
    'solver_label', "", ...
    'solver_method', "", ...
    'preconditioner', "", ...
    'linear_solver_tolerance', NaN, ...
    'boundary', "", ...
    'succeeded', false, ...
    'error_text', "", ...
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

function local_clear_hypre_instances()
try
    LINEAR_SOLVERS.hypre_boomeramg_clear();
catch
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
