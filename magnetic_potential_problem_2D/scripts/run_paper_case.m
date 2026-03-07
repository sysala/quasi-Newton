function result = run_paper_case(cfg)
%RUN_PAPER_CASE Execute one magnetic_potential_problem_2D benchmark case.
%
% The wrapper lifts local fields from a case struct into the legacy
% "magnetic_problem.m" workspace and collects a normalized result record.

paths = mp2d_paths();

% Geometry and material settings.
elem_type = local_cfg_field(cfg, 'elem_type', 'P1');
level = local_cfg_field(cfg, 'level', 3);
size_x = local_cfg_field(cfg, 'size_x', 1);
size_y = local_cfg_field(cfg, 'size_y', 1);
delta = local_cfg_field(cfg, 'delta', 1);
alpha = local_cfg_field(cfg, 'alpha', 0.0003);
beta = local_cfg_field(cfg, 'beta', 16000);
rho = local_cfg_field(cfg, 'rho', 40);
kappa = local_cfg_field(cfg, 'kappa', 25);

% Solver settings.
solver_type = string(local_cfg_field(cfg, 'solver_type', "DCG_HYPRE_BOOMERAMG"));
linear_solver_tolerance_newton = local_cfg_field(cfg, 'linear_solver_tolerance_newton', 1e-4);
linear_solver_tolerance_qn1 = local_cfg_field(cfg, 'linear_solver_tolerance_qn1', local_cfg_field(cfg, 'linear_solver_tolerance', 1e-1));
linear_solver_tolerance_qn2 = local_cfg_field(cfg, 'linear_solver_tolerance_qn2', local_cfg_field(cfg, 'linear_solver_tolerance', 1e-1));
linear_solver_maxit = local_cfg_field(cfg, 'linear_solver_maxit', 1000);
deflation_basis_tolerance = local_cfg_field(cfg, 'deflation_basis_tolerance', 1e-10);
linear_solver_printing = local_cfg_field(cfg, 'linear_solver_printing', 0);

boomeramg_opts = local_cfg_field(cfg, 'boomeramg_opts', struct());
boomeramg_opts = local_set_default_boomeramg_opts(boomeramg_opts);

nonlinear_method_max_runtime_seconds = local_cfg_field(cfg, 'nonlinear_method_max_runtime_seconds', 300);

run_newton = local_cfg_field(cfg, 'run_newton', true);
run_qn1 = local_cfg_field(cfg, 'run_qn1', true);
run_qn2 = local_cfg_field(cfg, 'run_qn2', true);
run_postprocess = local_cfg_field(cfg, 'run_postprocess', false);
use_damping = local_cfg_field(cfg, 'use_damping', false);

agmg_dir = '';

% Run the packaged benchmark in a controlled workspace.
evalc("run(fullfile(paths.root, 'magnetic_problem.m'));");

result = struct();
result.elem_type = elem_type;
result.level = level;
result.size_x = size_x;
result.size_y = size_y;
result.delta = delta;
result.alpha = alpha;
result.beta = beta;
result.rho = rho;
result.kappa = kappa;
result.solver_type = char(solver_type);
result.linear_solver_tolerance_newton = linear_solver_tolerance_newton;
result.linear_solver_tolerance_qn1 = linear_solver_tolerance_qn1;
result.linear_solver_tolerance_qn2 = linear_solver_tolerance_qn2;
result.linear_solver_maxit = linear_solver_maxit;
result.boomeramg_opts = boomeramg_opts;
result.use_damping = use_damping;
result.n_nodes = n_n;
result.n_unknowns = n_unknown;
result.n_elements = n_e;
result.n_integration_points = n_int;

result.newton = local_collect_solver_result(it0, itcg0, time_Newton, crit_hist0, timed_out0);
result.qn1 = local_collect_solver_result(it1, itcg1, time_qn1, crit_hist1, timed_out1);
result.qn2 = local_collect_solver_result(it2, itcg2, time_qn2, crit_hist2, timed_out2);
end

function value = local_cfg_field(cfg, field_name, default_value)
if isstruct(cfg) && isfield(cfg, field_name)
    value = cfg.(field_name);
else
    value = default_value;
end
end

function solver_opts = local_set_default_boomeramg_opts(solver_opts)
if ~isstruct(solver_opts)
    solver_opts = struct();
end

if ~isfield(solver_opts, 'threads')
    solver_opts.threads = 16;
end
if ~isfield(solver_opts, 'print_level')
    solver_opts.print_level = 0;
end
if ~isfield(solver_opts, 'use_as_preconditioner')
    solver_opts.use_as_preconditioner = true;
end
if ~isfield(solver_opts, 'num_functions')
    % Magnetic problem is a scalar PDE (no elasticity-like vector nullspace).
    solver_opts.num_functions = 1;
end
end

function solver_result = local_collect_solver_result(iterations, linear_iterations, runtime_seconds, crit_hist, timed_out)
solver_result = struct();
solver_result.iterations = iterations;
solver_result.linear_iterations = linear_iterations;
solver_result.runtime_seconds = runtime_seconds;
if isempty(crit_hist)
    solver_result.final_criterion = NaN;
    solver_result.converged = false;
else
    solver_result.final_criterion = crit_hist(end);
    solver_result.converged = crit_hist(end) < 1e-10;
end
solver_result.timed_out = logical(timed_out);
end
