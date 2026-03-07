function result = run_paper_case(cfg)
%RUN_PAPER_CASE Execute one nonlinear_Stokes_3D benchmark case.
%
% The benchmark entry point is a script. This helper wraps it in a function
% workspace so repeated sweeps do not leak variables between cases.

paths = ns3d_paths();

% Geometry and benchmark inputs.
elem_type = 'P1';
density = cfg.density;
size_xy_0 = cfg.size_xy_0;
size_xy_L = cfg.size_xy_L;
size_z = cfg.size_z;
mu_0 = cfg.mu_0;
mu_infty = cfg.mu_infty;
lambda = cfg.lambda;
p = cfg.p;
gamma = cfg.gamma;
boundary = cfg.boundary;
u_z = cfg.u_z;

% Linear solver controls.
solver_type = string(cfg.solver_type);
if isfield(cfg, 'linear_solver_tolerance')
    linear_solver_tolerance = cfg.linear_solver_tolerance;
end
if isfield(cfg, 'linear_solver_tolerance_newton')
    linear_solver_tolerance_newton = cfg.linear_solver_tolerance_newton;
end
if isfield(cfg, 'linear_solver_tolerance_quasi')
    linear_solver_tolerance_quasi = cfg.linear_solver_tolerance_quasi;
end
if isfield(cfg, 'linear_solver_tolerance_qn1')
    linear_solver_tolerance_qn1 = cfg.linear_solver_tolerance_qn1;
end
if isfield(cfg, 'linear_solver_tolerance_qn2')
    linear_solver_tolerance_qn2 = cfg.linear_solver_tolerance_qn2;
end
linear_solver_maxit = cfg.linear_solver_maxit;
deflation_basis_tolerance = cfg.deflation_basis_tolerance;
linear_solver_printing = 0;
boomeramg_opts = cfg.boomeramg_opts;
run_postprocess = false;
if isfield(cfg, 'nonlinear_method_max_runtime_seconds')
    nonlinear_method_max_runtime_seconds = cfg.nonlinear_method_max_runtime_seconds;
end
if isfield(cfg, 'run_newton')
    run_newton = logical(cfg.run_newton);
end
if isfield(cfg, 'run_qn1')
    run_qn1 = logical(cfg.run_qn1);
end
if isfield(cfg, 'run_qn2')
    run_qn2 = logical(cfg.run_qn2);
end

evalc("run(fullfile(paths.root, 'nonlinear_Stokes_3D.m'));");

result = struct();
result.elem_type = elem_type;
result.density = density;
result.size_xy_0 = size_xy_0;
result.size_xy_L = size_xy_L;
result.size_z = size_z;
result.boundary = boundary;
result.gamma = gamma;
result.solver_type = char(solver_type);
if exist('linear_solver_tolerance', 'var')
    result.linear_solver_tolerance = linear_solver_tolerance;
end
if exist('linear_solver_tolerance_newton', 'var')
    result.linear_solver_tolerance_newton = linear_solver_tolerance_newton;
end
if exist('linear_solver_tolerance_quasi', 'var')
    result.linear_solver_tolerance_quasi = linear_solver_tolerance_quasi;
end
if exist('linear_solver_tolerance_qn1', 'var')
    result.linear_solver_tolerance_qn1 = linear_solver_tolerance_qn1;
end
if exist('linear_solver_tolerance_qn2', 'var')
    result.linear_solver_tolerance_qn2 = linear_solver_tolerance_qn2;
end
if exist('nonlinear_method_max_runtime_seconds', 'var')
    result.nonlinear_method_max_runtime_seconds = nonlinear_method_max_runtime_seconds;
end
result.linear_solver_maxit = linear_solver_maxit;
result.deflation_basis_tolerance = deflation_basis_tolerance;
result.n_nodes = n_n;
result.n_unknowns = n_unknown;
result.n_elements = n_e;
result.n_integration_points = n_int;

result.newton = local_collect_solver_result(it_N, itcg_N, time_Newton, crit_hist_N, timed_out_N);
result.qn1 = local_collect_solver_result(it_qN1, itcg_qN1, time_qNewton1, crit_hist_qN1, timed_out_qN1);
result.qn2 = local_collect_solver_result(it_qN2, itcg_qN2, time_qNewton2, crit_hist_qN2, timed_out_qN2);

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
