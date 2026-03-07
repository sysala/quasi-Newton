function output_paths = plot_paper_tolerance_sweep(results_input)
%PLOT_PAPER_TOLERANCE_SWEEP Plot runtime-vs-tolerance figures from magnetic sweep.

if nargin < 1 || isempty(results_input)
    results_input = fullfile(mp2d_paths().results, 'paper_tolerance_sweep_results.mat');
end

results = local_load_results(results_input);
paths = mp2d_paths();

tolerances = [2e-1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5];
cases = local_cases();
mesh_variants = local_mesh_variants();
method_specs = local_method_specs();
styles = local_style_specs();

output_paths = strings(1, numel(method_specs));
for i_method = 1:numel(method_specs)
    out_path = fullfile(paths.images, method_specs(i_method).filename);
    local_make_figure(out_path, results, tolerances, cases, mesh_variants, method_specs(i_method), styles);
    output_paths(i_method) = string(out_path);
end
end

function results = local_load_results(results_input)
if isstruct(results_input)
    results = results_input;
    return
end

if isstring(results_input) || ischar(results_input)
    loaded = load(char(results_input), 'results');
    if ~isfield(loaded, 'results')
        error('plot_paper_tolerance_sweep:missingResults', ...
            'The checkpoint does not contain a ''results'' variable.');
    end
    results = loaded.results;
    return
end

error('plot_paper_tolerance_sweep:badInput', ...
    'Expected a results struct array or a path to paper_tolerance_sweep_results.mat.');
end

function cases = local_cases()
base = struct( ...
    'id', "case_delta_1", ...
    'name', "Inclusion contrast δ=1", ...
    'delta', 1, ...
    'display', "δ=1");
cases(1) = base;

base.id = "case_delta_2";
base.name = "Inclusion contrast δ=2";
base.delta = 2;
base.display = "δ=2";
cases(2) = base;
end

function mesh_variants = local_mesh_variants()
mesh_variants(1).label = "level_3";
mesh_variants(1).display_name = "Level 3";
mesh_variants(1).level = 3;

mesh_variants(2).label = "level_4";
mesh_variants(2).display_name = "Level 4";
mesh_variants(2).level = 4;
end

function method_specs = local_method_specs()
method_specs = struct([]);

method_specs(1).field = "newton";
method_specs(1).label = "N";
method_specs(1).title = "Newton";
method_specs(1).filename = "paper_tolerance_sweep_n.png";

method_specs(2).field = "qn1";
method_specs(2).label = "qN1";
method_specs(2).title = "Quasi-Newton 1";
method_specs(2).filename = "paper_tolerance_sweep_qn1.png";

method_specs(3).field = "qn2";
method_specs(3).label = "qN2";
method_specs(3).title = "Quasi-Newton 2";
method_specs(3).filename = "paper_tolerance_sweep_qn2.png";
end

function style_specs = local_style_specs()
style_specs = struct();
style_specs.palette = [ ...
    0.00, 0.45, 0.74; ...
    0.85, 0.33, 0.10; ...
    0.47, 0.67, 0.19; ...
    0.55, 0.34, 0.67];
style_specs.markers = {'o', 's', 'd', '^'};
end

function local_make_figure(out_path, results, tolerances, cases, mesh_variants, method_spec, styles)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1500, 820]);
cleanup = onCleanup(@() close(fig));

method_field = char(method_spec.field);
hold on;
box on;
grid on;
ax = gca;
ax.FontSize = 11;
ax.LineWidth = 0.9;
ax.GridAlpha = 0.18;
ax.MinorGridAlpha = 0.08;
ax.YScale = 'log';

legend_labels = strings(1, 0);
legend_handles = gobjects(1, 0);

for i_case = 1:numel(cases)
    for i_mesh = 1:numel(mesh_variants)
        color = styles.palette(mod((i_case - 1) * numel(mesh_variants) + i_mesh - 1, ...
            size(styles.palette, 1)) + 1, :);
        marker = styles.markers{mod((i_case - 1) * numel(mesh_variants) + i_mesh - 1, ...
            numel(styles.markers)) + 1};

        y_values = nan(1, numel(tolerances));
        for i_tol = 1:numel(tolerances)
            entry = local_find_entry(results, cases(i_case).id, mesh_variants(i_mesh).label, tolerances(i_tol));
            if isempty(entry) || ~entry.succeeded
                continue
            end
            solver_result = entry.(method_field);
            y_values(i_tol) = solver_result.runtime_seconds;
        end

        valid = isfinite(y_values);
        if ~any(valid)
            continue
        end

        h = plot(1:numel(tolerances), y_values, ...
            'Color', color, ...
            'Marker', marker, ...
            'MarkerSize', 7, ...
            'LineStyle', '-', ...
            'LineWidth', 1.3);

        legend_handles(end + 1) = h;
        legend_labels(end + 1) = sprintf('%s / %s', cases(i_case).name, mesh_variants(i_mesh).display_name);
    end
end

if ~isempty(legend_handles)
    legend(legend_handles, cellstr(legend_labels), ...
        'Location', 'bestoutside', 'Box', 'off');
end

xticks(1:numel(tolerances));
xticklabels(arrayfun(@(x) sprintf('%.0e', x), tolerances, 'UniformOutput', false));
xlabel('Linear solver tolerance');
ylabel('Runtime [s]');
title(sprintf('%s runtime sweep (DCG + HYPRE)', method_spec.title), 'FontWeight', 'bold');
grid on;

method_field = char(method_spec.field);
[runtime_min, runtime_max, has_runtime] = local_runtime_limits(results, method_field);
if has_runtime
    if runtime_max <= runtime_min
        runtime_max = runtime_min * 10;
    end
    ylim([runtime_min, runtime_max * 1.2]);
end

exportgraphics(fig, out_path, 'Resolution', 180);
end

function [runtime_min, runtime_max, has_data] = local_runtime_limits(results, method_field)
runtime_min = inf;
runtime_max = 0;
has_data = false;

for i = 1:numel(results)
    if ~isfield(results(i), method_field)
        continue
    end
    if ~isfield(results(i).(method_field), 'runtime_seconds')
        continue
    end
    solver_result = results(i).(method_field);
    if ~isfield(results(i), 'succeeded') || ~logical(results(i).succeeded)
        continue
    end
    if ~isfinite(solver_result.runtime_seconds)
        continue
    end
    runtime_min = min(runtime_min, solver_result.runtime_seconds);
    runtime_max = max(runtime_max, solver_result.runtime_seconds);
    has_data = true;
end
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
