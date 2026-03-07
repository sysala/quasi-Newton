function run_selected_tolerance_mesh_extension(options)
%RUN_SELECTED_TOLERANCE_MESH_EXTENSION Run the ~200k selected-tolerance mesh cases.
%
% This runner executes only the chosen benchmark policy:
%   N   -> 1e-4
%   qN1 -> 1e-1
%   qN2 -> 1e-1
% on the added ~200k-unknown mesh, with a 300 s cap per nonlinear method.
%
% Optional input:
%   options - struct with optional filters:
%       example_ids   : cellstr, e.g. {'example1'}
%       solver_types  : cellstr
%       gammas        : numeric vector
%       max_cases     : positive integer
%       force_rerun   : logical scalar
%       run_single_pending : logical scalar

if nargin < 1 || isempty(options)
    options = struct();
end

paths = ns3d_paths();
results_path = fullfile(paths.results, 'selected_tolerance_mesh_results_dcg_hypre.mat');

mesh_variant = struct( ...
    'label', "roughly_200k_unknowns_mesh", ...
    'display_name', "~200k unknowns mesh", ...
    'density', 19);
solver_specs = local_solver_specs();
examples = local_examples();
plan = local_build_plan(examples, solver_specs, mesh_variant, options);

results = struct([]);
if isfile(results_path)
    loaded = load(results_path, 'results');
    if isfield(loaded, 'results')
        results = loaded.results;
    end
end

force_rerun = isfield(options, 'force_rerun') && logical(options.force_rerun);
run_single_pending = isfield(options, 'run_single_pending') && logical(options.run_single_pending);

for i = 1:numel(plan)
    base_key = plan(i).base_key;
    if ~force_rerun && local_has_complete_entries(results, base_key)
        fprintf('[%d/%d] skip %s\n', i, numel(plan), base_key);
        continue
    end

    fprintf('[%d/%d] run %s\n', i, numel(plan), base_key);
    try
        case_result = run_paper_case(plan(i).cfg);
        entries = local_result_entries(plan(i), case_result, "", true);
    catch ME
        entries = local_failure_entries(plan(i), ...
            string(getReport(ME, 'extended', 'hyperlinks', 'off')));
    end

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

function examples = local_examples()
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
    'nonlinear_method_max_runtime_seconds', 300);

examples = repmat(common, 1, 2);

examples(1).id = "example1";
examples(1).name = "Example 1: Prismatic Bar";
examples(1).size_xy_0 = 5;
examples(1).size_xy_L = 5;
examples(1).boundary = "BC1";
examples(1).gammas = [0.003, 0.004, 0.005, 0.006, 0.007];

examples(2).id = "example2";
examples(2).name = "Example 2: Tapered Prismatic Bar";
examples(2).size_xy_0 = 6;
examples(2).size_xy_L = 5;
examples(2).boundary = "BC1";
examples(2).gammas = [0.003, 0.004, 0.005, 0.006, 0.007, 0.008];
end

function solver_specs = local_solver_specs()
solver_specs = struct( ...
    'solver_type', {"DCG_HYPRE_BOOMERAMG"}, ...
    'label', {"DCG + HYPRE"}, ...
    'method', {"DCG"}, ...
    'preconditioner', {"HYPRE"});
end

function plan = local_build_plan(examples, solver_specs, mesh_variant, options)
plan = struct([]);
counter = 0;

for i_example = 1:numel(examples)
    ex = examples(i_example);
    if ~local_keep_value(ex.id, options, 'example_ids')
        continue
    end

    for i_solver = 1:numel(solver_specs)
        solver_spec = solver_specs(i_solver);
        if ~local_keep_value(solver_spec.solver_type, options, 'solver_types')
            continue
        end

        for i_gamma = 1:numel(ex.gammas)
            gamma = ex.gammas(i_gamma);
            if ~local_keep_value(gamma, options, 'gammas')
                continue
            end

            counter = counter + 1;
            plan(counter).base_key = sprintf('%s__%s__%s__g_%s', ...
                char(ex.id), char(mesh_variant.label), char(solver_spec.solver_type), local_number_token(gamma)); %#ok<AGROW>
            plan(counter).example_id = ex.id;
            plan(counter).example_name = ex.name;
            plan(counter).mesh_label = mesh_variant.label;
            plan(counter).mesh_display_name = mesh_variant.display_name;
            plan(counter).density = mesh_variant.density;
            plan(counter).solver_type = solver_spec.solver_type;
            plan(counter).solver_label = solver_spec.label;
            plan(counter).solver_method = solver_spec.method;
            plan(counter).preconditioner = solver_spec.preconditioner;
            plan(counter).gamma = gamma;
            plan(counter).boundary = ex.boundary;
            plan(counter).cfg = ex;
            plan(counter).cfg.density = mesh_variant.density;
            plan(counter).cfg.solver_type = solver_spec.solver_type;
            plan(counter).cfg.gamma = gamma;
        end
    end
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

entries = repmat(local_result_entry_template(), 1, 2);

entries(1) = local_fill_entry(plan_entry, case_result, error_text, succeeded, tol_newton, ...
    case_result.newton, local_unused_solver_result(), local_unused_solver_result());
entries(2) = local_fill_entry(plan_entry, case_result, error_text, succeeded, tol_quasi, ...
    local_unused_solver_result(), case_result.qn1, case_result.qn2);
end

function entry = local_fill_entry(plan_entry, case_result, error_text, succeeded, tol, newton_result, qn1_result, qn2_result)
entry = struct();
entry.key = sprintf('%s__tol_%s', char(plan_entry.base_key), local_number_token(tol));
entry.example_id = plan_entry.example_id;
entry.example_name = plan_entry.example_name;
entry.mesh_label = plan_entry.mesh_label;
entry.mesh_display_name = plan_entry.mesh_display_name;
entry.density = plan_entry.density;
entry.solver_type = plan_entry.solver_type;
entry.solver_label = plan_entry.solver_label;
entry.solver_method = plan_entry.solver_method;
entry.preconditioner = plan_entry.preconditioner;
entry.linear_solver_tolerance = tol;
entry.gamma = plan_entry.gamma;
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
failure_result = struct('n_unknowns', NaN, 'n_elements', NaN, ...
    'newton', local_unused_solver_result(), ...
    'qn1', local_unused_solver_result(), ...
    'qn2', local_unused_solver_result());
entries = local_result_entries(plan_entry, failure_result, error_text, false);
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
    'example_id', "", ...
    'example_name', "", ...
    'mesh_label', "", ...
    'mesh_display_name', "", ...
    'density', NaN, ...
    'solver_type', "", ...
    'solver_label', "", ...
    'solver_method', "", ...
    'preconditioner', "", ...
    'linear_solver_tolerance', NaN, ...
    'gamma', NaN, ...
    'boundary', "", ...
    'succeeded', false, ...
    'error_text', "", ...
    'n_unknowns', NaN, ...
    'n_elements', NaN, ...
    'newton', local_unused_solver_result(), ...
    'qn1', local_unused_solver_result(), ...
    'qn2', local_unused_solver_result());
end
