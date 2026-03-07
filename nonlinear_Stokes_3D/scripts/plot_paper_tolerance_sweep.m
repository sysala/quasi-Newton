function output_paths = plot_paper_tolerance_sweep(results_input)
%PLOT_PAPER_TOLERANCE_SWEEP Generate runtime-vs-tolerance figures from the sweep checkpoint.

if nargin < 1 || isempty(results_input)
    results_input = fullfile(ns3d_paths().results, 'paper_tolerance_sweep_results.mat');
end

results = local_load_results(results_input);
paths = ns3d_paths();

tolerances = [2e-1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5];
examples = local_examples();
mesh_variants = local_mesh_variants();
method_specs = local_method_specs();
style_specs = local_style_specs();

output_paths = strings(1, numel(method_specs));
for i_method = 1:numel(method_specs)
    out_path = fullfile(paths.images, method_specs(i_method).filename);
    local_make_figure(out_path, results, tolerances, examples, mesh_variants, method_specs(i_method), style_specs);
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

function examples = local_examples()
examples = struct([]);

examples(1).id = "example1";
examples(1).short_label = "Example 1";
examples(1).name = "Prismatic Bar";
examples(1).gammas = [0.003, 0.004, 0.005, 0.006, 0.007];

examples(2).id = "example2";
examples(2).short_label = "Example 2";
examples(2).name = "Tapered Bar";
examples(2).gammas = [0.003, 0.004, 0.005, 0.006, 0.007, 0.008];
end

function mesh_variants = local_mesh_variants()
mesh_variants = struct([]);

mesh_variants(1).label = "paper_mesh";
mesh_variants(1).short_label = "Paper mesh";

mesh_variants(2).label = "double_unknowns_mesh";
mesh_variants(2).short_label = "~2x unknowns";
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

style_specs.preconditioners = struct( ...
    'name', {"IC", "AGMG", "HYPRE"}, ...
    'color', { ...
        [0.00, 0.45, 0.74], ...
        [0.85, 0.33, 0.10], ...
        [0.47, 0.67, 0.19]});

style_specs.krylov = struct( ...
    'name', {"DCG", "DFGMRES"}, ...
    'line_style', {"-", "--"});

style_specs.gamma_markers = struct( ...
    'value', {0.003, 0.004, 0.005, 0.006, 0.007, 0.008}, ...
    'marker', {'o', 's', 'd', '^', 'v', 'p'});
end

function local_make_figure(out_path, results, tolerances, examples, mesh_variants, method_spec, style_specs)
fig = figure( ...
    'Visible', 'off', ...
    'Color', 'w', ...
    'Position', [50, 50, 1800, 980]);

cleanup = onCleanup(@() close(fig));
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

method_field = char(method_spec.field);
[runtime_min, runtime_max] = local_runtime_limits(results, method_field);

tile_indices = [1, 2, 4, 5];
panel_counter = 0;
for i_example = 1:numel(examples)
    for i_mesh = 1:numel(mesh_variants)
        panel_counter = panel_counter + 1;
        ax = nexttile(layout, tile_indices(panel_counter));
        hold(ax, 'on');
        box(ax, 'on');
        grid(ax, 'on');
        ax.GridAlpha = 0.12;
        ax.MinorGridAlpha = 0.08;
        ax.FontSize = 10;
        ax.LineWidth = 0.9;

        example = examples(i_example);
        mesh_variant = mesh_variants(i_mesh);
        local_plot_panel(ax, results, tolerances, example, mesh_variant, method_field, style_specs);

        unknowns = local_unknown_count(results, example.id, mesh_variant.label);
        if isnan(unknowns)
            title(ax, sprintf('%s / %s', example.short_label, mesh_variant.short_label), ...
                'FontWeight', 'bold');
        else
            title(ax, { ...
                sprintf('%s / %s', example.short_label, mesh_variant.short_label), ...
                sprintf('%d unknowns', unknowns)}, ...
                'FontWeight', 'bold');
        end

        xticks(ax, 1:numel(tolerances));
        xticklabels(ax, compose('%.0e', tolerances));
        xlim(ax, [0.8, numel(tolerances) + 0.2]);
        ax.YScale = 'log';
        ax.YMinorGrid = 'on';
        ylim(ax, [runtime_min, runtime_max]);

        if panel_counter > 2
            xlabel(ax, 'Linear solver tolerance');
        end
        if any(panel_counter == [1, 3])
            ylabel(ax, 'Runtime [s]');
        end
    end
end

ax_info = nexttile(layout, 3, [2, 1]);
local_draw_info_panel(ax_info, style_specs);

sgtitle(layout, sprintf('%s: runtime vs tolerance', method_spec.label), ...
    'FontSize', 18, ...
    'FontWeight', 'bold');

exportgraphics(fig, out_path, 'Resolution', 180);
end

function local_plot_panel(ax, results, tolerances, example, mesh_variant, method_field, style_specs)
x_values = 1:numel(tolerances);

for i_prec = 1:numel(style_specs.preconditioners)
    prec = style_specs.preconditioners(i_prec);
    for i_krylov = 1:numel(style_specs.krylov)
        krylov = style_specs.krylov(i_krylov);
        solver_type = string(sprintf('%s_%s', krylov.name, local_solver_preconditioner_name(prec.name)));

        for i_gamma = 1:numel(example.gammas)
            gamma = example.gammas(i_gamma);
            marker = local_gamma_marker(gamma, style_specs.gamma_markers);

            y_values = nan(1, numel(tolerances));
            converged = false(1, numel(tolerances));

            for i_tol = 1:numel(tolerances)
                entry = local_find_entry(results, example.id, mesh_variant.label, solver_type, tolerances(i_tol), gamma);
                if isempty(entry) || ~entry.succeeded
                    continue
                end

                solver_result = entry.(method_field);
                y_values(i_tol) = solver_result.runtime_seconds;
                converged(i_tol) = logical(solver_result.converged);
            end

            valid = isfinite(y_values);
            if ~any(valid)
                continue
            end

            plot(ax, x_values, y_values, ...
                'Color', prec.color, ...
                'LineStyle', krylov.line_style, ...
                'LineWidth', 1.2);

            converged_mask = valid & converged;
            if any(converged_mask)
                plot(ax, x_values(converged_mask), y_values(converged_mask), ...
                    'LineStyle', 'none', ...
                    'Marker', marker, ...
                    'MarkerSize', 6.5, ...
                    'MarkerEdgeColor', prec.color, ...
                    'MarkerFaceColor', prec.color);
            end

            stalled_mask = valid & ~converged;
            if any(stalled_mask)
                plot(ax, x_values(stalled_mask), y_values(stalled_mask), ...
                    'LineStyle', 'none', ...
                    'Marker', marker, ...
                    'MarkerSize', 7.5, ...
                    'LineWidth', 1.2, ...
                    'MarkerEdgeColor', prec.color, ...
                    'MarkerFaceColor', 'w');
            end
        end
    end
end
end

function ax = local_draw_info_panel(ax, style_specs)
hold(ax, 'on');
axis(ax, [0, 1, 0, 1]);
axis(ax, 'off');

text(ax, 0.02, 0.97, 'Encoding', 'FontSize', 13, 'FontWeight', 'bold', 'VerticalAlignment', 'top');

text(ax, 0.02, 0.88, 'Color = preconditioner', 'FontSize', 11, 'FontWeight', 'bold');
for i = 1:numel(style_specs.preconditioners)
    y = 0.82 - 0.07 * (i - 1);
    plot(ax, [0.05, 0.18], [y, y], 'Color', style_specs.preconditioners(i).color, 'LineWidth', 2.2);
    text(ax, 0.22, y, style_specs.preconditioners(i).name, 'FontSize', 10, 'VerticalAlignment', 'middle');
end

text(ax, 0.02, 0.57, 'Line style = Krylov method', 'FontSize', 11, 'FontWeight', 'bold');
for i = 1:numel(style_specs.krylov)
    y = 0.51 - 0.07 * (i - 1);
    plot(ax, [0.05, 0.18], [y, y], 'Color', [0.15, 0.15, 0.15], ...
        'LineStyle', style_specs.krylov(i).line_style, 'LineWidth', 2.2);
    text(ax, 0.22, y, style_specs.krylov(i).name, 'FontSize', 10, 'VerticalAlignment', 'middle');
end

text(ax, 0.02, 0.33, 'Marker = |g_3|', 'FontSize', 11, 'FontWeight', 'bold');
for i = 1:numel(style_specs.gamma_markers)
    y = 0.27 - 0.045 * (i - 1);
    plot(ax, 0.10, y, ...
        'LineStyle', 'none', ...
        'Marker', style_specs.gamma_markers(i).marker, ...
        'MarkerSize', 7, ...
        'MarkerEdgeColor', [0.15, 0.15, 0.15], ...
        'MarkerFaceColor', [0.15, 0.15, 0.15]);
    text(ax, 0.22, y, sprintf('%.0e', style_specs.gamma_markers(i).value), ...
        'FontSize', 10, 'VerticalAlignment', 'middle');
end

plot(ax, 0.10, 0.02, ...
    'LineStyle', 'none', ...
    'Marker', 'o', ...
    'MarkerSize', 8, ...
    'LineWidth', 1.2, ...
    'MarkerEdgeColor', [0.15, 0.15, 0.15], ...
    'MarkerFaceColor', 'w');
text(ax, 0.22, 0.02, 'Hollow marker: outer solve hit the iteration limit', ...
    'FontSize', 10, 'VerticalAlignment', 'middle');
end

function solver_preconditioner = local_solver_preconditioner_name(preconditioner)
switch string(preconditioner)
    case "IC"
        solver_preconditioner = "ICHOL";
    case "AGMG"
        solver_preconditioner = "AGMG";
    case "HYPRE"
        solver_preconditioner = "HYPRE_BOOMERAMG";
    otherwise
        error('plot_paper_tolerance_sweep:unknownPreconditioner', ...
            'Unsupported preconditioner "%s".', preconditioner);
end
end

function marker = local_gamma_marker(gamma, gamma_markers)
marker = 'o';
for i = 1:numel(gamma_markers)
    if abs(gamma - gamma_markers(i).value) < 1e-12
        marker = gamma_markers(i).marker;
        return
    end
end
end

function [runtime_min, runtime_max] = local_runtime_limits(results, method_field)
runtime_min = inf;
runtime_max = 0;

for i = 1:numel(results)
    if ~results(i).succeeded
        continue
    end

    solver_result = results(i).(method_field);
    if ~isfield(solver_result, 'runtime_seconds')
        continue
    end

    if solver_result.runtime_seconds > 0
        runtime_min = min(runtime_min, solver_result.runtime_seconds);
    end
    runtime_max = max(runtime_max, solver_result.runtime_seconds);
end

if ~isfinite(runtime_min) || runtime_min <= 0
    runtime_min = 1e-1;
else
    runtime_min = runtime_min / 1.15;
end

if ~(isfinite(runtime_max) && runtime_max > 0)
    runtime_max = 1;
else
    runtime_max = runtime_max * 1.08;
end
end

function unknowns = local_unknown_count(results, example_id, mesh_label)
unknowns = NaN;

for i = 1:numel(results)
    if ~results(i).succeeded
        continue
    end

    if string(results(i).example_id) ~= string(example_id)
        continue
    end

    if string(results(i).mesh_label) ~= string(mesh_label)
        continue
    end

    unknowns = results(i).n_unknowns;
    return
end
end

function entry = local_find_entry(results, example_id, mesh_label, solver_type, tol, gamma)
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
