% =========================================================================
%
%  Magnetic potential problem in 2D (packaged version).
%
%  This program compares Newton and quasi-Newton methods for the 2D magnetic
%  potential problem. The structure mirrors the other benchmark packages in this
%  repository with packaged assembly, mesh and nonlinear solvers.
%
% =========================================================================

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

% ------------------------------------------------------------------------
% Inputs (can be overridden from caller workspace in wrapper scripts)
% ------------------------------------------------------------------------
if ~exist('elem_type', 'var')
    elem_type = 'P1';
end
if ~exist('level', 'var')
    level = 3;
end
if ~exist('size_x', 'var')
    size_x = 1;
end
if ~exist('size_y', 'var')
    size_y = 1;
end
if ~exist('delta', 'var')
    delta = 1;
end
if ~exist('alpha', 'var')
    alpha = 0.0003;
end
if ~exist('beta', 'var')
    beta = 16000;
end
if ~exist('rho', 'var')
    rho = 40;
end
if ~exist('kappa', 'var')
    kappa = 25;
end
if ~exist('solver_type', 'var')
    solver_type = "DCG_HYPRE_BOOMERAMG";
end
if ~exist('linear_solver_tolerance_newton', 'var')
    linear_solver_tolerance_newton = 1e-4;
end
if ~exist('linear_solver_tolerance_qn1', 'var')
    if exist('linear_solver_tolerance', 'var')
        linear_solver_tolerance_qn1 = linear_solver_tolerance;
    else
        linear_solver_tolerance_qn1 = 1e-1;
    end
end
if ~exist('linear_solver_tolerance_qn2', 'var')
    if exist('linear_solver_tolerance', 'var')
        linear_solver_tolerance_qn2 = linear_solver_tolerance;
    else
        linear_solver_tolerance_qn2 = 1e-1;
    end
end
if ~exist('linear_solver_maxit', 'var')
    linear_solver_maxit = 1000;
end
if ~exist('deflation_basis_tolerance', 'var')
    deflation_basis_tolerance = 1e-10;
end
if ~exist('linear_solver_printing', 'var')
    linear_solver_printing = 0;
end
if ~exist('boomeramg_opts', 'var')
    boomeramg_opts = struct('threads', 16, 'print_level', 0, ...
        'use_as_preconditioner', true);
end
if ~exist('nonlinear_method_max_runtime_seconds', 'var')
    nonlinear_method_max_runtime_seconds = 300;
end
if ~exist('run_newton', 'var')
    run_newton = true;
end
if ~exist('run_qn1', 'var')
    run_qn1 = true;
end
if ~exist('run_qn2', 'var')
    run_qn2 = true;
end
if ~exist('run_postprocess', 'var')
    run_postprocess = false;
end
if ~exist('use_damping', 'var')
    use_damping = false;
end

agmg_dir = '';

% ------------------------------------------------------------------------
% Material data
% ------------------------------------------------------------------------
r_crit = a_function(kappa, alpha, beta);

% ------------------------------------------------------------------------
% Mesh generation
% ------------------------------------------------------------------------
switch(elem_type)
    case 'P1'
        [COORD, ELEM, heter, Q] = MESH.mesh_P1(level, size_x, size_y, delta);
        fprintf('P1 elements:\n');
    case 'P2'
        [COORD, ELEM, heter, Q] = MESH.mesh_P2(level, size_x, size_y, delta);
        fprintf('P2 elements:\n');
    otherwise
        error('Bad choice of element type.');
end

% Sizes and quick printout
n_n = size(COORD, 2);
n_unknown = sum(Q(:));
n_e = size(ELEM, 2);

[Xi, WF] = ASSEMBLY.quadrature_volume_2D(elem_type);
[Xi_s, WF_s] = ASSEMBLY.quadrature_surface(elem_type); %#ok<NASGU>
[HatP, DHatP1, DHatP2] = ASSEMBLY.local_basis_volume_2D(elem_type, Xi);
[HatP_s, DHatP1_s] = ASSEMBLY.local_basis_surface(elem_type, Xi_s); %#ok<NASGU>

n_q = length(WF);
n_int = n_e * n_q;
fprintf('\nThe coarsest mesh:\n');
fprintf('  number of nodes = %d\n', n_n);
fprintf('  number of unknowns = %d\n', n_unknown);
fprintf('  number of elements = %d\n', n_e);
fprintf('  number of integration points = %d\n', n_int);

% ------------------------------------------------------------------------
% Assembly and forcing vector
% ------------------------------------------------------------------------
heter_int = repmat(heter, n_q, 1);
heter_int = heter_int(:);

[K_fix, B, WEIGHT] = ASSEMBLY.stiffness_matrix(ELEM, COORD, heter_int, Q, DHatP1, DHatP2, WF, alpha);
f_V = ASSEMBLY.vector_volume(ELEM, COORD, HatP, WEIGHT);
f = rho * f_V(Q);

% ------------------------------------------------------------------------
% Prepare linear solvers (copy before using, so each nonlinear run can keep own state)
% ------------------------------------------------------------------------
linear_system_solver_newton = LINEAR_SOLVERS.set_linear_solver(agmg_dir, solver_type, ...
    linear_solver_tolerance_newton, linear_solver_maxit, deflation_basis_tolerance, ...
    linear_solver_printing, Q, COORD, boomeramg_opts);
linear_system_solver_qn1 = LINEAR_SOLVERS.set_linear_solver(agmg_dir, solver_type, ...
    linear_solver_tolerance_qn1, linear_solver_maxit, deflation_basis_tolerance, ...
    linear_solver_printing, Q, COORD, boomeramg_opts);
linear_system_solver_qn2 = LINEAR_SOLVERS.set_linear_solver(agmg_dir, solver_type, ...
    linear_solver_tolerance_qn2, linear_solver_maxit, deflation_basis_tolerance, ...
    linear_solver_printing, Q, COORD, boomeramg_opts);

% ------------------------------------------------------------------------
% Newton-like solvers
% ------------------------------------------------------------------------
U0_n = nan;
it0 = NaN;
itcg0 = NaN;
crit_hist0 = [];
time_Newton = NaN;
timed_out0 = false;

U1 = nan;
it1 = NaN;
itcg1 = NaN;
crit_hist1 = [];
omega_hist1 = [];
time_qn1 = NaN;
timed_out1 = false;

U2 = nan;
it2 = NaN;
itcg2 = NaN;
crit_hist2 = [];
omega_hist2 = [];
time_qn2 = NaN;
timed_out2 = false;

U_ini = zeros(n_unknown, 1);

if run_newton
    tic;
    if use_damping
        [U0_n, it0, crit_hist0, ~] = NEWTON.newton_damped(U_ini, WEIGHT, K_fix, B, f, ...
            heter_int, alpha, beta);
        itcg0 = 0;
        timed_out0 = false;
    else
        [U0_n, it0, crit_hist0, itcg0, timed_out0] = NEWTON.newton(U_ini, WEIGHT, K_fix, B, f, ...
            heter_int, alpha, beta, linear_system_solver_newton.copy(), nonlinear_method_max_runtime_seconds);
    end
    time_Newton = toc;
end

if run_qn1
    tic;
    if use_damping
        [U1, it1, crit_hist1, omega_hist1] = NEWTON.newton_quasi1_damped(U_ini, WEIGHT, K_fix, B, f, ...
            heter_int, alpha, beta);
        itcg1 = 0;
        timed_out1 = false;
    else
        [U1, it1, crit_hist1, omega_hist1, itcg1, timed_out1] = NEWTON.newton_quasi1(U_ini, WEIGHT, K_fix, B, f, ...
            heter_int, alpha, beta, linear_system_solver_qn1.copy(), nonlinear_method_max_runtime_seconds);
    end
    time_qn1 = toc;
end

if run_qn2
    tic;
    if use_damping
        [U2, it2, crit_hist2, omega_hist2] = NEWTON.newton_quasi2_damped(U_ini, WEIGHT, K_fix, B, f, ...
            heter_int, r_crit, alpha, beta);
        itcg2 = 0;
        timed_out2 = false;
    else
        [U2, it2, crit_hist2, omega_hist2, itcg2, timed_out2] = NEWTON.newton_quasi2(U_ini, WEIGHT, K_fix, B, f, ...
            heter_int, r_crit, alpha, beta, linear_system_solver_qn2.copy(), nonlinear_method_max_runtime_seconds);
    end
    time_qn2 = toc;
end

if run_postprocess
    U0_full = zeros(n_n, 1);
    if ~isnan(U0_n(1))
        U0_full(Q) = U0_n;
    end

    E = B * U0_n;
    grad_mag = zeros(1, n_int);
    if ~isempty(E)
        E = reshape(E, 2, n_int);
        grad_mag = sqrt(E(1,:).^2 + E(2,:).^2);
        grad_node = VIZ.transformation(grad_mag, ELEM, WEIGHT);
    else
        grad_node = zeros(1, n_n);
    end

    VIZ.draw_quantity(COORD, ELEM, U0_full', size_x, size_y);
    VIZ.draw_quantity(COORD, ELEM, grad_node, size_x, size_y);
end

% ------------------------------------------------------------------------
% Result record (consumed by wrapper scripts)
% ------------------------------------------------------------------------
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
result.use_damping = use_damping;
result.linear_solver_tolerance_newton = linear_solver_tolerance_newton;
result.linear_solver_tolerance_qn1 = linear_solver_tolerance_qn1;
result.linear_solver_tolerance_qn2 = linear_solver_tolerance_qn2;
result.linear_solver_maxit = linear_solver_maxit;
result.n_nodes = n_n;
result.n_unknowns = n_unknown;
result.n_elements = n_e;
result.n_integration_points = n_int;
result.newton = struct('iterations', it0, 'linear_iterations', itcg0, ...
    'runtime_seconds', time_Newton, 'final_criterion', local_final_crit(crit_hist0), ...
    'timed_out', timed_out0);
result.qn1 = struct('iterations', it1, 'linear_iterations', itcg1, ...
    'runtime_seconds', time_qn1, 'final_criterion', local_final_crit(crit_hist1), ...
    'timed_out', timed_out1);
result.qn2 = struct('iterations', it2, 'linear_iterations', itcg2, ...
    'runtime_seconds', time_qn2, 'final_criterion', local_final_crit(crit_hist2), ...
    'timed_out', timed_out2);

% Keep old variable naming compatibility where scripts expect these names.
it0d = it0;
it1d = it1;
it2d = it2;

% ------------------------------------------------------------------------
% Helpers
% ------------------------------------------------------------------------
function value = local_final_crit(crit_hist)
    if isempty(crit_hist)
        value = NaN;
    else
        value = crit_hist(end);
    end
end
