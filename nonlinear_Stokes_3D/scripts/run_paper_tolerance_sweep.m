function run_paper_tolerance_sweep(options)
%RUN_PAPER_TOLERANCE_SWEEP Sweep paper benchmark settings over solver tolerances.
%
% The sweep is resumable. Results are saved after every case to:
%   paper_tolerance_sweep_results.mat
% and the markdown report is regenerated at the same cadence.
%
% Optional input:
%   options - struct with optional filters:
%       example_ids  : cellstr, e.g. {'example1'}
%       densities    : numeric vector, e.g. [12 15]
%       solver_types : cellstr
%       tolerances   : numeric vector
%       max_cases    : positive integer
%       force_rerun  : logical scalar
%       regenerate_only : logical scalar
%       skip_plots   : logical scalar

if nargin < 1 || isempty(options)
    options = struct();
end

paths = ns3d_paths();
results_path = fullfile(paths.results, 'paper_tolerance_sweep_results.mat');
markdown_path = fullfile(paths.root, 'PAPER_TOLERANCE_SWEEP.md');

examples = local_examples();
solver_specs = local_solver_specs();
mesh_variants = local_mesh_variants();
tolerances = [2e-1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5];

plan = local_build_plan(examples, solver_specs, mesh_variants, tolerances, options);

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
        key = plan(i).key;
        if ~force_rerun && local_find_result_index(results, key) > 0
            fprintf('[%d/%d] skip %s\n', i, numel(plan), key);
            continue
        end

        fprintf('[%d/%d] run %s\n', i, numel(plan), key);
        try
            case_result = run_paper_case(plan(i).cfg);
            entry = local_result_entry(plan(i), case_result, "", true);
        catch ME
            entry = local_result_entry(plan(i), local_failure_result(), ...
                string(getReport(ME, 'extended', 'hyperlinks', 'off')), false);
        end

        if isempty(results)
            results = entry;
        elseif ~isequal(sort(string(fieldnames(results))), sort(string(fieldnames(entry))))
            warning('run_paper_tolerance_sweep:resetCheckpoint', ...
                'Discarding incompatible checkpoint structure and starting a fresh results array.');
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
        local_write_markdown(markdown_path, results, plan, examples, solver_specs, mesh_variants, tolerances);
    end
elseif isempty(results)
    error('run_paper_tolerance_sweep:noCheckpoint', ...
        'regenerate_only=true requires an existing results checkpoint at %s.', results_path);
end

local_write_markdown(markdown_path, results, plan, examples, solver_specs, mesh_variants, tolerances);
if ~skip_plots
    plot_paper_tolerance_sweep(results_path);
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
    'boomeramg_opts', struct('threads', 16, 'print_level', 0, 'use_as_preconditioner', true));

examples = repmat(common, 1, 2);

examples(1).id = "example1";
examples(1).name = "Example 1: Prismatic Bar";
examples(1).size_xy_0 = 5;
examples(1).size_xy_L = 5;
examples(1).boundary = "BC1";
examples(1).boundary_note = "Matches the paper text and reproduces the reported unknown count 54,886.";
examples(1).paper_unknowns = 54886;
examples(1).paper_elements = 103680;
examples(1).gammas = [0.003, 0.004, 0.005, 0.006, 0.007];

examples(2).id = "example2";
examples(2).name = "Example 2: Tapered Prismatic Bar";
examples(2).size_xy_0 = 6;
examples(2).size_xy_L = 5;
examples(2).boundary = "BC1";
examples(2).boundary_note = ...
    "paper.tex states u_3 = 0 on the shell, but that is inconsistent with a constant inlet profile. " ...
    + "A calibration run showed BC1 (u·n = 0 on the shell) is much closer to the paper table than BC3.";
examples(2).paper_unknowns = 50986;
examples(2).paper_elements = 86400;
examples(2).gammas = [0.003, 0.004, 0.005, 0.006, 0.007, 0.008];
end

function mesh_variants = local_mesh_variants()
mesh_variants(1).label = "paper_mesh";
mesh_variants(1).display_name = "Paper-size mesh";
mesh_variants(1).density = 12;

mesh_variants(2).label = "double_unknowns_mesh";
mesh_variants(2).display_name = "~2x unknowns mesh";
mesh_variants(2).density = 15;
end

function solver_specs = local_solver_specs()
solver_specs = struct( ...
    'solver_type', { ...
        "DCG_ICHOL", ...
        "DCG_AGMG", ...
        "DCG_HYPRE_BOOMERAMG", ...
        "DFGMRES_ICHOL", ...
        "DFGMRES_AGMG", ...
        "DFGMRES_HYPRE_BOOMERAMG"}, ...
    'label', { ...
        "DCG + IC", ...
        "DCG + AGMG", ...
        "DCG + HYPRE", ...
        "DFGMRES + IC", ...
        "DFGMRES + AGMG", ...
        "DFGMRES + HYPRE"}, ...
    'method', { ...
        "DCG", "DCG", "DCG", "DFGMRES", "DFGMRES", "DFGMRES"}, ...
    'preconditioner', { ...
        "IC", "AGMG", "HYPRE", "IC", "AGMG", "HYPRE"});
end

function plan = local_build_plan(examples, solver_specs, mesh_variants, tolerances, options)
plan = struct([]);
case_counter = 0;

for i_example = 1:numel(examples)
    ex = examples(i_example);
    if ~local_keep_value(ex.id, options, 'example_ids')
        continue
    end

    for i_mesh = 1:numel(mesh_variants)
        mesh_variant = mesh_variants(i_mesh);
        if ~local_keep_value(mesh_variant.density, options, 'densities')
            continue
        end

        for i_solver = 1:numel(solver_specs)
            solver_spec = solver_specs(i_solver);
            if ~local_keep_value(solver_spec.solver_type, options, 'solver_types')
                continue
            end

            for i_tol = 1:numel(tolerances)
                tol = tolerances(i_tol);
                if ~local_keep_value(tol, options, 'tolerances')
                    continue
                end

                for i_gamma = 1:numel(ex.gammas)
                    gamma = ex.gammas(i_gamma);
                    case_counter = case_counter + 1;
                    plan(case_counter).key = local_case_key(ex.id, mesh_variant.label, solver_spec.solver_type, tol, gamma); %#ok<AGROW>
                    plan(case_counter).example_id = ex.id;
                    plan(case_counter).example_name = ex.name;
                    plan(case_counter).mesh_label = mesh_variant.label;
                    plan(case_counter).mesh_display_name = mesh_variant.display_name;
                    plan(case_counter).density = mesh_variant.density;
                    plan(case_counter).solver_type = solver_spec.solver_type;
                    plan(case_counter).solver_label = solver_spec.label;
                    plan(case_counter).solver_method = solver_spec.method;
                    plan(case_counter).preconditioner = solver_spec.preconditioner;
                    plan(case_counter).linear_solver_tolerance = tol;
                    plan(case_counter).gamma = gamma;
                    plan(case_counter).boundary = ex.boundary;
                    plan(case_counter).cfg = ex;
                    plan(case_counter).cfg.density = mesh_variant.density;
                    plan(case_counter).cfg.solver_type = solver_spec.solver_type;
                    plan(case_counter).cfg.linear_solver_tolerance = tol;
                    plan(case_counter).cfg.gamma = gamma;
                end
            end
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

function key = local_case_key(example_id, mesh_label, solver_type, tol, gamma)
key = sprintf('%s__%s__%s__tol_%s__g_%s', ...
    char(example_id), ...
    char(mesh_label), ...
    char(solver_type), ...
    local_number_token(tol), ...
    local_number_token(gamma));
end

function token = local_number_token(value)
token = regexprep(sprintf('%.0e', value), '\+', '');
token = strrep(token, '.', 'p');
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

function entry = local_result_entry(plan_entry, case_result, error_text, succeeded)
entry = struct();
entry.key = plan_entry.key;
entry.example_id = plan_entry.example_id;
entry.example_name = plan_entry.example_name;
entry.mesh_label = plan_entry.mesh_label;
entry.mesh_display_name = plan_entry.mesh_display_name;
entry.density = plan_entry.density;
entry.solver_type = plan_entry.solver_type;
entry.solver_label = plan_entry.solver_label;
entry.solver_method = plan_entry.solver_method;
entry.preconditioner = plan_entry.preconditioner;
entry.linear_solver_tolerance = plan_entry.linear_solver_tolerance;
entry.gamma = plan_entry.gamma;
entry.boundary = plan_entry.boundary;
entry.succeeded = succeeded;
entry.error_text = error_text;
entry.n_unknowns = case_result.n_unknowns;
entry.n_elements = case_result.n_elements;
entry.newton = case_result.newton;
entry.qn1 = case_result.qn1;
entry.qn2 = case_result.qn2;
end

function result = local_failure_result()
solver_stub = struct( ...
    'iterations', NaN, ...
    'linear_iterations', NaN, ...
    'runtime_seconds', NaN, ...
    'final_criterion', NaN, ...
    'converged', false, ...
    'timed_out', false);

result = struct( ...
    'n_unknowns', NaN, ...
    'n_elements', NaN, ...
    'newton', solver_stub, ...
    'qn1', solver_stub, ...
    'qn2', solver_stub);
end

function local_write_markdown(markdown_path, results, plan, examples, solver_specs, mesh_variants, tolerances)
fid = fopen(markdown_path, 'w');
if fid < 0
    error('Unable to open markdown output: %s', markdown_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Paper Tolerance Sweep for `nonlinear_Stokes_3D`\n\n');
fprintf(fid, 'This report keeps the original paper results and appends a tolerance sweep run on this machine.\n\n');
fprintf(fid, '## Runtime Plots\n\n');
fprintf(fid, 'The figures below show runtime versus linear solver tolerance with uniform spacing between tolerance labels and a log-scale runtime axis.\n');
fprintf(fid, 'Each figure contains four data panels: Example 1/2 crossed with the paper mesh and the `~2x unknowns` mesh.\n\n');
fprintf(fid, '- color: preconditioner\n');
fprintf(fid, '- line style: Krylov method (`DCG` solid, `DFGMRES` dashed)\n');
fprintf(fid, '- marker: `|g_3|`\n');
fprintf(fid, '- hollow markers: the outer nonlinear method hit the iteration limit\n\n');
fprintf(fid, '### N\n\n');
fprintf(fid, '![Runtime vs tolerance for Newton](images/paper_tolerance_sweep_n.png)\n\n');
fprintf(fid, '### qN1\n\n');
fprintf(fid, '![Runtime vs tolerance for qN1](images/paper_tolerance_sweep_qn1.png)\n\n');
fprintf(fid, '### qN2\n\n');
fprintf(fid, '![Runtime vs tolerance for qN2](images/paper_tolerance_sweep_qn2.png)\n\n');
fprintf(fid, 'Generated from:\n');
fprintf(fid, '- [paper.tex](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/paper.tex)\n');
fprintf(fid, '- [PAPER_EXPERIMENTS.md](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/PAPER_EXPERIMENTS.md)\n');
fprintf(fid, '- [run_paper_tolerance_sweep.m](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/scripts/run_paper_tolerance_sweep.m)\n');
fprintf(fid, '- [run_paper_case.m](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/scripts/run_paper_case.m)\n');
fprintf(fid, '- [plot_paper_tolerance_sweep.m](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/scripts/plot_paper_tolerance_sweep.m)\n');
fprintf(fid, '- [paper_tolerance_sweep_results.mat](/home/beremi/repos/quasi-Newton/nonlinear_Stokes_3D/results/paper_tolerance_sweep_results.mat)\n\n');

fprintf(fid, 'Sweep parameters:\n');
fprintf(fid, '- tolerances: `%s`\n', strjoin(compose('%.0e', tolerances), '`, `'));
fprintf(fid, '- solver variants: `%s`\n', strjoin(string({solver_specs.label}), '`, `'));
fprintf(fid, '- mesh variants: `%s`\n', strjoin(string({mesh_variants.display_name}), '`, `'));
fprintf(fid, '- reported cell format in sweep tables: `outer iterations / cumulative linear iterations / runtime [s]`\n');
fprintf(fid, '- completed cases: `%d / %d`\n\n', numel(results), numel(plan));

local_write_fastest_tolerance_summary(fid, results, tolerances);

fprintf(fid, '## Original Paper Results\n\n');
local_write_original_tables(fid);

local_write_selected_tolerance_replications(fid, results, examples, mesh_variants);

fprintf(fid, '\n## Sweep Results on This Machine\n\n');
fprintf(fid, 'Notes:\n');
fprintf(fid, '- Example 1 uses `BC1`, which reproduces the paper unknown count exactly.\n');
fprintf(fid, '- Example 2 is ambiguous in the paper text. The sweep uses `BC1` because it is materially closer to the reported table than the literal `u_3 = 0` shell variant.\n');
fprintf(fid, '- The `~2x unknowns mesh` corresponds to density `15`, which yields approximately double the unknown count compared with the paper-size density `12` mesh.\n\n');

for i_example = 1:numel(examples)
    ex = examples(i_example);
    fprintf(fid, '### %s\n\n', ex.name);
    fprintf(fid, '- boundary used in sweep: `%s`\n', ex.boundary);
    fprintf(fid, '- note: %s\n', ex.boundary_note);
    fprintf(fid, '- paper reported mesh: `%d` elements, `%d` unknowns\n\n', ex.paper_elements, ex.paper_unknowns);

    for i_mesh = 1:numel(mesh_variants)
        mesh_variant = mesh_variants(i_mesh);
        fprintf(fid, '#### %s\n\n', mesh_variant.display_name);
        mesh_entry = local_find_first_result(results, ex.id, mesh_variant.label);
        if ~isempty(mesh_entry)
            fprintf(fid, '- density: `%d`\n', mesh_variant.density);
            fprintf(fid, '- measured elements: `%d`\n', mesh_entry.n_elements);
            fprintf(fid, '- measured unknowns: `%d`\n\n', mesh_entry.n_unknowns);
        else
            fprintf(fid, '- density: `%d`\n\n', mesh_variant.density);
        end

        for i_solver = 1:numel(solver_specs)
            solver_spec = solver_specs(i_solver);
            fprintf(fid, '##### %s\n\n', solver_spec.label);
            fprintf(fid, '| tol | g3 | unknowns | N | qN1 | qN2 | status |\n');
            fprintf(fid, '|---|---|---|---|---|---|---|\n');

            for i_tol = 1:numel(tolerances)
                tol = tolerances(i_tol);
                for i_gamma = 1:numel(ex.gammas)
                    gamma = ex.gammas(i_gamma);
                    key = local_case_key(ex.id, mesh_variant.label, solver_spec.solver_type, tol, gamma);
                    idx = local_find_result_index(results, key);
                    if idx == 0
                        fprintf(fid, '| `%s` | `-%s` | - | PENDING | PENDING | PENDING | pending |\n', ...
                            local_tol_label(tol), local_gamma_label(gamma));
                        continue
                    end

                    entry = results(idx);
                    if entry.succeeded
                        fprintf(fid, '| `%s` | `-%s` | `%d` | `%s` | `%s` | `%s` | ok |\n', ...
                            local_tol_label(tol), ...
                            local_gamma_label(gamma), ...
                            entry.n_unknowns, ...
                            local_solver_cell(entry.newton, 100), ...
                            local_solver_cell(entry.qn1, 100), ...
                            local_solver_cell(entry.qn2, 200));
                    else
                        short_error = local_shorten_error(entry.error_text);
                        fprintf(fid, '| `%s` | `-%s` | - | ERROR | ERROR | ERROR | %s |\n', ...
                            local_tol_label(tol), local_gamma_label(gamma), short_error);
                    end
                end
            end

            fprintf(fid, '\n');
        end
    end
end
end

function entry = local_find_first_result(results, example_id, mesh_label)
entry = [];
for i = 1:numel(results)
    if results(i).example_id == example_id && results(i).mesh_label == mesh_label && results(i).succeeded
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

function label = local_tol_label(tol)
label = sprintf('%.0e', tol);
end

function label = local_gamma_label(gamma)
label = sprintf('%.0e', gamma);
end

function short_error = local_shorten_error(error_text)
short_error = char(error_text);
short_error = regexprep(short_error, '\s+', ' ');
short_error = strrep(short_error, '|', '/');
if numel(short_error) > 80
    short_error = [short_error(1:77) '...'];
end
end

function local_write_fastest_tolerance_summary(fid, results, tolerances)
summary = local_fastest_tolerance_summary(results, tolerances);

fprintf(fid, '## Fastest Tolerance Summary\n\n');
fprintf(fid, ['Counts are taken over base experiment settings `example x mesh x solver variant x |g_3|` ' ...
    'with a complete tolerance sweep and successful solves.\n']);
fprintf(fid, ['Cell format: `within-method wins (overall wins)` where `overall wins` means the same ' ...
    '`method + tolerance` was the fastest across `N`, `qN1`, and `qN2` for that base experiment.\n\n']);
fprintf(fid, '- complete base experiments included: `%d`\n', summary.complete_base_count);
fprintf(fid, '- incomplete base experiments excluded from the count: `%d`\n', summary.incomplete_base_count);
if summary.method_ties > 0 || summary.overall_ties > 0
    fprintf(fid, '- ties counted for every tied winner: method ties `%d`, overall ties `%d`\n', ...
        summary.method_ties, summary.overall_ties);
end
fprintf(fid, '\n');

fprintf(fid, '| tol | N | qN1 | qN2 |\n');
fprintf(fid, '|---|---:|---:|---:|\n');
for i_tol = 1:numel(tolerances)
    fprintf(fid, '| `%s` | `%d (%d)` | `%d (%d)` | `%d (%d)` |\n', ...
        local_tol_label(tolerances(i_tol)), ...
        summary.within_counts(1, i_tol), summary.overall_counts(1, i_tol), ...
        summary.within_counts(2, i_tol), summary.overall_counts(2, i_tol), ...
        summary.within_counts(3, i_tol), summary.overall_counts(3, i_tol));
end
fprintf(fid, '\n');
end

function local_write_selected_tolerance_replications(fid, results, examples, mesh_variants)
selected_results = local_selected_results_merged(results);
selected_mesh_variants = local_selected_mesh_variants(mesh_variants);
families = local_selected_linear_families();
methods = local_selected_method_specs();

fprintf(fid, '## Selected-Tolerance Replications\n\n');
fprintf(fid, 'These tables use the current benchmark policy: `N` at `1e-4` and `qN1/qN2` at `1e-1`.\n');
fprintf(fid, 'Cell format remains `outer nonlinear iterations / cumulative linear iterations / runtime [s]`.\n');
fprintf(fid, '- this section keeps only the `DCG + HYPRE` cut; legacy solver variants stay in the full sweep below\n');
fprintf(fid, ['`<ins><strong>...</strong></ins>` marks the fastest converged runtime in the row ' ...
    'within this `DCG + HYPRE` cut.\n']);
fprintf(fid, '- `~~>300s~~` marks a nonlinear method stopped by the 300 s cap\n\n');

for i_example = 1:numel(examples)
    ex = examples(i_example);
    fprintf(fid, '### %s\n\n', ex.name);
    fprintf(fid, '- boundary used in replication: `%s`\n', ex.boundary);
    fprintf(fid, '- note: %s\n\n', ex.boundary_note);

    for i_mesh = 1:numel(selected_mesh_variants)
        mesh_variant = selected_mesh_variants(i_mesh);
        fprintf(fid, '#### %s\n\n', mesh_variant.display_name);
        mesh_entry = local_find_first_result(selected_results, ex.id, mesh_variant.label);
        if ~isempty(mesh_entry)
            fprintf(fid, '- density: `%d`\n', mesh_variant.density);
            fprintf(fid, '- measured elements: `%d`\n', mesh_entry.n_elements);
            fprintf(fid, '- measured unknowns: `%d`\n\n', mesh_entry.n_unknowns);
        else
            fprintf(fid, '- density: `%d`\n\n', mesh_variant.density);
        end

        for i_family = 1:numel(families)
            family = families(i_family);
            fprintf(fid, '##### %s\n\n', family.display_name);
            fprintf(fid, '| `g_3` | `HYPRE-N` | `HYPRE-qN1` | `HYPRE-qN2` |\n');
            fprintf(fid, '|---|---|---|---|\n');

            for i_gamma = 1:numel(ex.gammas)
                gamma = ex.gammas(i_gamma);
                row_cells = local_selected_row_cells(selected_results, ex.id, mesh_variant.label, gamma, families, methods);
                fprintf(fid, '| `-%s` ', local_gamma_label(gamma));
                for i_prec = 1:numel(family.preconditioner_labels)
                    for i_method = 1:numel(methods)
                        cell_data = row_cells(i_family, i_prec, i_method);
                        fprintf(fid, '| %s ', local_selected_markdown_cell(cell_data));
                    end
                end
                fprintf(fid, '|\n');
            end
            fprintf(fid, '\n');
        end
    end
end
end

function merged_results = local_selected_results_merged(base_results)
merged_results = base_results;
selected_results_path = fullfile(ns3d_paths().results, 'selected_tolerance_mesh_results_dcg_hypre.mat');
if isfile(selected_results_path)
    loaded = load(selected_results_path, 'results');
    if isfield(loaded, 'results') && ~isempty(loaded.results)
        merged_results = [merged_results, loaded.results];
    end
end
end

function selected_mesh_variants = local_selected_mesh_variants(base_mesh_variants)
selected_mesh_variants = base_mesh_variants;
selected_results_path = fullfile(ns3d_paths().results, 'selected_tolerance_mesh_results_dcg_hypre.mat');
if isfile(selected_results_path)
    selected_mesh_variants(end + 1).label = "roughly_200k_unknowns_mesh"; %#ok<AGROW>
    selected_mesh_variants(end).display_name = "~200k unknowns mesh";
    selected_mesh_variants(end).density = 19;
end
end

function summary = local_fastest_tolerance_summary(results, tolerances)
method_fields = ["newton", "qn1", "qn2"];
within_counts = zeros(numel(method_fields), numel(tolerances));
overall_counts = zeros(numel(method_fields), numel(tolerances));
method_ties = 0;
overall_ties = 0;

base_keys = strings(0, 1);
for i = 1:numel(results)
    base_keys(end + 1, 1) = local_base_case_key(results(i)); %#ok<AGROW>
end
base_keys = unique(base_keys);

complete_base_count = 0;
incomplete_base_count = 0;

for i_base = 1:numel(base_keys)
    base_key = base_keys(i_base);
    base_entries = local_collect_base_entries(results, base_key, tolerances);
    if isempty(base_entries)
        incomplete_base_count = incomplete_base_count + 1;
        continue
    end

    complete_base_count = complete_base_count + 1;
    runtime_matrix = nan(numel(method_fields), numel(tolerances));
    for i_tol = 1:numel(tolerances)
        for i_method = 1:numel(method_fields)
            runtime_matrix(i_method, i_tol) = base_entries(i_tol).(method_fields(i_method)).runtime_seconds;
        end
    end

    for i_method = 1:numel(method_fields)
        method_runtimes = runtime_matrix(i_method, :);
        fastest_tol = local_winner_indices(method_runtimes);
        if numel(fastest_tol) > 1
            method_ties = method_ties + 1;
        end
        within_counts(i_method, fastest_tol) = within_counts(i_method, fastest_tol) + 1;
    end

    overall_linear = runtime_matrix(:).';
    overall_winners = local_winner_indices(overall_linear);
    if numel(overall_winners) > 1
        overall_ties = overall_ties + 1;
    end
    for i_winner = 1:numel(overall_winners)
        linear_idx = overall_winners(i_winner);
        method_idx = mod(linear_idx - 1, numel(method_fields)) + 1;
        tol_idx = floor((linear_idx - 1) / numel(method_fields)) + 1;
        overall_counts(method_idx, tol_idx) = overall_counts(method_idx, tol_idx) + 1;
    end
end

summary = struct( ...
    'within_counts', within_counts, ...
    'overall_counts', overall_counts, ...
    'complete_base_count', complete_base_count, ...
    'incomplete_base_count', incomplete_base_count, ...
    'method_ties', method_ties, ...
    'overall_ties', overall_ties);
end

function families = local_selected_linear_families()
families = struct([]);

families(1).display_name = "DCG with HYPRE";
families(1).solver_types = ["DCG_HYPRE_BOOMERAMG"];
families(1).preconditioner_labels = ["HYPRE"];
end

function methods = local_selected_method_specs()
methods = struct([]);

methods(1).field = "newton";
methods(1).tol = 1e-4;
methods(1).iteration_limit = 100;

methods(2).field = "qn1";
methods(2).tol = 1e-1;
methods(2).iteration_limit = 100;

methods(3).field = "qn2";
methods(3).tol = 1e-1;
methods(3).iteration_limit = 200;
end

function row_cells = local_selected_row_cells(results, example_id, mesh_label, gamma, families, methods)
template = struct('label', '-', 'runtime', inf, 'converged', false, 'timed_out', false, 'highlight', false);
row_cells = repmat(template, numel(families), numel(families(1).solver_types), numel(methods));

best_runtime = inf;
for i_family = 1:numel(families)
    for i_solver = 1:numel(families(i_family).solver_types)
        solver_type = families(i_family).solver_types(i_solver);
        for i_method = 1:numel(methods)
            entry = local_selected_find_result(results, example_id, mesh_label, solver_type, methods(i_method).tol, gamma);
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

function entry = local_selected_find_result(results, example_id, mesh_label, solver_type, tol, gamma)
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

function rendered = local_selected_markdown_cell(cell_data)
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

function base_entries = local_collect_base_entries(results, base_key, tolerances)
base_entries = struct([]);

for i_tol = 1:numel(tolerances)
    matched = [];
    for i = 1:numel(results)
        if local_base_case_key(results(i)) ~= base_key
            continue
        end
        if abs(results(i).linear_solver_tolerance - tolerances(i_tol)) > eps(max(1, abs(tolerances(i_tol))))
            continue
        end
        if ~results(i).succeeded
            base_entries = struct([]);
            return
        end
        matched = results(i);
        break
    end

    if isempty(matched)
        base_entries = struct([]);
        return
    end

    if any(~isfinite([matched.newton.runtime_seconds, matched.qn1.runtime_seconds, matched.qn2.runtime_seconds]))
        base_entries = struct([]);
        return
    end

    if isempty(base_entries)
        base_entries = matched;
    else
        base_entries(end + 1) = matched; %#ok<AGROW>
    end
end
end

function key = local_base_case_key(entry)
key = entry.example_id + "|" + entry.mesh_label + "|" + entry.solver_type + "|" + local_number_token(entry.gamma);
end

function winner_idx = local_winner_indices(values)
min_value = min(values);
tolerance = max(1e-9, 1e-9 * abs(min_value));
winner_idx = find(abs(values - min_value) <= tolerance);
end

function local_write_original_tables(fid)
fprintf(fid, '### Example 1: Paper Table\n\n');
fprintf(fid, '| `g_3` | `IC-N` | `IC-qN1` | `IC-qN2` | `AGMG-N` | `AGMG-qN1` | `AGMG-qN2` |\n');
fprintf(fid, '|---|---|---|---|---|---|---|\n');
fprintf(fid, '| `-3e-3` | `6/112/1.86` | `14/84/0.95` | `15/58/0.61` | `6/51/2.12` | `14/35/0.99` | `15/44/0.85` |\n');
fprintf(fid, '| `-4e-3` | `6/111/1.84` | `17/90/1.17` | `22/70/0.79` | `6/53/2.01` | `17/40/1.25` | `22/57/1.11` |\n');
fprintf(fid, '| `-5e-3` | `7/130/2.17` | `27/97/1.72` | `35/93/1.12` | `7/63/2.44` | `25/50/1.85` | `36/86/1.68` |\n');
fprintf(fid, '| `-6e-3` | `7/133/2.17` | `41/111/2.69` | `74/167/2.11` | `7/63/2.39` | `41/70/3.05` | `75/170/3.33` |\n');
fprintf(fid, '| `-7e-3` | `8/152/2.53` | `73/143/5.17` | `>200/-/5.53` | `8/76/2.84` | `72/108/5.63` | `>200/-/6.81` |\n\n');

fprintf(fid, '### Example 2: Paper Table\n\n');
fprintf(fid, '| `g_3` | `IC-N` | `IC-qN1` | `IC-qN2` | `AGMG-N` | `AGMG-qN1` | `AGMG-qN2` |\n');
fprintf(fid, '|---|---|---|---|---|---|---|\n');
fprintf(fid, '| `-3e-3` | `5/86/1.28` | `12/72/0.76` | `14/56/0.52` | `5/39/1.52` | `12/33/0.90` | `13/39/0.74` |\n');
fprintf(fid, '| `-4e-3` | `6/100/1.57` | `16/77/0.93` | `17/59/0.59` | `6/50/1.80` | `16/37/1.06` | `17/51/0.93` |\n');
fprintf(fid, '| `-5e-3` | `6/101/1.56` | `19/83/1.15` | `23/72/0.73` | `6/50/1.88` | `19/43/1.35` | `24/70/1.25` |\n');
fprintf(fid, '| `-6e-3` | `7/119/1.85` | `25/90/1.56` | `36/97/1.07` | `7/61/2.21` | `25/53/1.87` | `37/95/1.80` |\n');
fprintf(fid, '| `-7e-3` | `7/120/1.87` | `38/102/2.32` | `68/162/1.77` | `7/61/2.17` | `37/68/2.71` | `70/167/3.09` |\n');
fprintf(fid, '| `-8e-3` | `8/139/2.20` | `62/126/3.87` | `168/347/4.15` | `8/71/2.52` | `61/93/4.44` | `198/359/7.28` |\n\n');
end
