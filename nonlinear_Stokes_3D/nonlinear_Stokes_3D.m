% =========================================================================
%
%  This program compares Newton and quasi-Newton methods for a static
%  and simplified nonlinear Stokes problem in 3D.
%  It is considered a benchmark with a prismatic bar.
%  One can set optionally 2 types of finite elements,
%  different levels of mesh density and many other parameters and
%  inner solvers. More details can be found in the paper
%  J. Karatson, S. Sysala. M. Beres:  Quasi-Newton variable preconditioning
%  for nonlinear elasticity systems in 3D, CAMWA, 2025.
%
% ======================================================================
%

% Resolve local dependencies relative to this script.
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
agmg_dir = fullfile(script_dir, 'agmg');

if exist(fullfile(agmg_dir, 'agmg.m'), 'file') == 2
    addpath(agmg_dir);
end

%
% Main input data
%

% Element, mesh and geometrical data
if ~exist('elem_type', 'var')
    elem_type='P1'; % available choices of finite elements: 'P1', 'P2'
end
% for P2 elements, it is necessary to change
% the regularization parameter within inner solvers
if ~exist('density', 'var')
    density=12;      % a positive integer defining mesh density in x-direction
end
if ~exist('size_xy_0', 'var')
    size_xy_0 = 6;  % size of the body in directions x and y on the left
end
if ~exist('size_xy_L', 'var')
    size_xy_L = 5;  % size of the body in directions x and y on the right
end
if ~exist('size_z', 'var')
    size_z = 50;      % size of the body in z-direction
end
% It is assumed that size_xy_L < size_xy_0 < size_z.

% Material parameters
if ~exist('mu_0', 'var')
    mu_0 = 1;
end
if ~exist('mu_infty', 'var')
    mu_infty = 1e-3;
end
if ~exist('lambda', 'var')
    lambda = 10;
end
if ~exist('p', 'var')
    p = 1.1;
end

% Constant volume force representing pore pressure gradient
if ~exist('gamma', 'var')
    gamma = 0.008;
end

% Boundary conditions - two available choices:
% 'BC1' - zero velocity field in the normal direction on the shell
% 'BC2' - zero velocity field in all directions on the shell
% 'BC3' - paper-style variant: u_z fixed on the shell, constant inlet in u_z
if ~exist('boundary', 'var')
    boundary='BC1';
end
if ~exist('u_z', 'var')
    u_z = 5;  % prescribed velocity in the direction z and at the prismatic centre
end

% Linear solver configuration
if ~exist('solver_type', 'var')
    solver_type = "DCG_HYPRE_BOOMERAMG";
end
if ~exist('linear_solver_tolerance_newton', 'var')
    if exist('linear_solver_tolerance', 'var')
        linear_solver_tolerance_newton = linear_solver_tolerance;
    else
        linear_solver_tolerance_newton = 1e-4;
    end
end
if ~exist('linear_solver_tolerance_qn1', 'var')
    if exist('linear_solver_tolerance_quasi', 'var')
        linear_solver_tolerance_qn1 = linear_solver_tolerance_quasi;
    elseif exist('linear_solver_tolerance', 'var')
        linear_solver_tolerance_qn1 = linear_solver_tolerance;
    else
        linear_solver_tolerance_qn1 = 1e-1;
    end
end
if ~exist('linear_solver_tolerance_qn2', 'var')
    if exist('linear_solver_tolerance_quasi', 'var')
        linear_solver_tolerance_qn2 = linear_solver_tolerance_quasi;
    elseif exist('linear_solver_tolerance', 'var')
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
if ~exist('run_postprocess', 'var')
    run_postprocess = true;
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

%
% Mesh generation
%
switch(elem_type)
    case 'P1'
        [COORD,ELEM,SURF1,SURF2,SURF3,SURF4,SURF5,SURF6]=...
            MESH.mesh_P1(density,size_xy_0,size_xy_L,size_z);
        fprintf('P1 elements: \n')
    case 'P2'
        [COORD,ELEM,SURF1,SURF2,SURF3,SURF4,SURF5,SURF6]=...
            MESH.mesh_P2(density,size_xy_0,size_xy_L,size_z);
        fprintf('P2 elements: \n')
    otherwise
        disp('bad choice of element type');
end
SURF=[SURF1,SURF2,SURF3,SURF4,SURF5,SURF6];

%
% Arrays representing the prescribed boundary conditions
%
switch(boundary)
    case 'BC1'
        [U_3,Q]=MESH.boundary_BC1(COORD,SURF1,SURF3,SURF4,SURF5,SURF6);
    case 'BC2'
        [U_3,Q]=MESH.boundary_BC2(COORD,SURF1,SURF3,SURF4,SURF5,SURF6,size_xy_0);
    case 'BC3'
        [U_3,Q]=MESH.boundary_BC3(COORD,SURF1,SURF3,SURF4,SURF5,SURF6);
    otherwise
        disp('bad choice of element type');
end
U_3=u_z*U_3;

%
% Data from the reference element
%

% quadrature points and weights for volume and surface integration
[Xi, WF] = ASSEMBLY.quadrature_volume(elem_type);
[Xi_s, WF_s] = ASSEMBLY.quadrature_surface(elem_type);

% local basis functions and their derivatives for volume and surface
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume(elem_type, Xi);
[HatP_s, DHatP1_s, DHatP2_s] = ASSEMBLY.local_basis_surface(elem_type, Xi_s);

%
% Number of nodes, elements and integration points + print
%
n_n = size(COORD, 2); % number of nodes
n_unknown = length(COORD(Q)); % number of unknowns
n_e = size(ELEM, 2); % number of elements
n_q = length(WF); % number of quadratic points
n_int = n_e * n_q; % total number of integrations points
fprintf('number of nodes =%d ',n_n);
fprintf('\n');
fprintf('number of unknowns =%d ', n_unknown);
fprintf('\n');
fprintf('number of elements =%d ',n_e);
fprintf('\n');
fprintf('number of integration points =%d ',n_int);
fprintf('\n');

%
% Values of material parameters at integration points
%
mu_0 = mu_0*ones(1,n_int) ;
mu_infty = mu_infty*ones(1,n_int) ;
lambda =lambda*ones(1,n_int);
gamma=gamma*ones(1,n_int);

%
% Assembling of auxiliary arrays for Newton's and quasi-Newton's methods
%
[B,K_elast,WEIGHT]=ASSEMBLY.auxiliary_matrices(ELEM,COORD,mu_0,...
    DHatP1,DHatP2,DHatP3,WF);

%
% Assembling of the vector of volume forces.
%
% Volume forces at integration points: f_V_int has size (3, n_int).
f_V_int = [zeros(1, n_int); zeros(1, n_int); -gamma];
% Compute the volume force vector.
f = ASSEMBLY.vector_volume_3D(ELEM, COORD, f_V_int, HatP, WEIGHT);

%
% Newton's and quasi-Newton's solvers
%

linear_system_solver_newton = LINEAR_SOLVERS.set_linear_solver(agmg_dir, solver_type, ...
    linear_solver_tolerance_newton, linear_solver_maxit, deflation_basis_tolerance, ...
    linear_solver_printing, Q, COORD, boomeramg_opts);
linear_system_solver_qn1 = LINEAR_SOLVERS.set_linear_solver(agmg_dir, solver_type, ...
    linear_solver_tolerance_qn1, linear_solver_maxit, deflation_basis_tolerance, ...
    linear_solver_printing, Q, COORD, boomeramg_opts);
linear_system_solver_qn2 = LINEAR_SOLVERS.set_linear_solver(agmg_dir, solver_type, ...
    linear_solver_tolerance_qn2, linear_solver_maxit, deflation_basis_tolerance, ...
    linear_solver_printing, Q, COORD, boomeramg_opts);

% initialization displacement
U_it = [zeros(1,n_n); zeros(1,n_n); U_3] ;

% standard Newton method
U_N = nan(size(U_it));
it_N = NaN;
crit_hist_N = [];
itcg_N = NaN;
timed_out_N = false;
time_Newton = NaN;
if run_newton
    tic;
    [U_N, it_N, crit_hist_N, itcg_N, timed_out_N]=NEWTON.newton(U_it,WEIGHT,B,f,Q,mu_0,mu_infty,lambda,p, ...
        linear_system_solver_newton.copy(), nonlinear_method_max_runtime_seconds);
    time_Newton=toc;
    fprintf("     solver's runtime:  " + time_Newton + "-----\n");
end

% quasi-Newton method - preconditioner 1 (simplified stiffness matrix with variable coefficients)
U_qN1 = nan(size(U_it));
it_qN1 = NaN;
crit_hist_qN1 = [];
omega_hist_qN1 = [];
itcg_qN1 = NaN;
timed_out_qN1 = false;
time_qNewton1 = NaN;
if run_qn1
    tic;
    [U_qN1, it_qN1, crit_hist_qN1, omega_hist_qN1, itcg_qN1, timed_out_qN1]=NEWTON.newton_quasi1(U_it,WEIGHT,K_elast,B,f,Q,mu_0,mu_infty,lambda,p, ...
        linear_system_solver_qn1.copy(), nonlinear_method_max_runtime_seconds);
    time_qNewton1=toc;
    fprintf("     solver's runtime:  " + time_qNewton1 + "-----\n");
end

% quasi-Newton method - preconditioner 2 (simplified stiffness matrix with fixed coefficients)
U_qN2 = nan(size(U_it));
it_qN2 = NaN;
crit_hist_qN2 = [];
omega_hist_qN2 = [];
itcg_qN2 = NaN;
timed_out_qN2 = false;
time_qNewton2 = NaN;
if run_qn2
    tic;
    [U_qN2, it_qN2, crit_hist_qN2, omega_hist_qN2, itcg_qN2, timed_out_qN2]=NEWTON.newton_quasi2(U_it,WEIGHT,K_elast,B,f,Q,mu_0,mu_infty,lambda,p, ...
        linear_system_solver_qn2.copy(), nonlinear_method_max_runtime_seconds);
    time_qNewton2=toc;
    fprintf("     solver's runtime:  " + time_qNewton2 + "-----\n");
end

if contains(upper(string(solver_type)), "BOOMERAMG")
    LINEAR_SOLVERS.hypre_boomeramg_clear();
end

%
% Visualization of selected results
%
if run_postprocess && run_newton
    % mesh
    if density<5
        VIZ.draw_mesh(COORD,SURF,elem_type)
    end

    % total displacements + deformed shape - newton
    U_total = sqrt(U_N(1,:).^2 + U_N(2,:).^2 + U_N(3,:).^2);
    VIZ.draw_quantity(COORD,SURF,1*U_N,U_total,elem_type,size_xy_0,size_xy_L,size_z)

    % U_total = sqrt(U_qN1(1,:).^2 + U_qN1(2,:).^2 + U_qN1(3,:).^2);
    % VIZ.draw_quantity(COORD,SURF,1*U_qN1,U_total,elem_type,size_xy_0,size_xy_L,size_z)
    % U_total = sqrt(U_qN2(1,:).^2 + U_qN2(2,:).^2 + U_qN2(3,:).^2);
    % VIZ.draw_quantity(COORD,SURF,1*U_qN2,U_total,elem_type,size_xy_0,size_xy_L,size_z)

    % values of the function a
    E= reshape( B*U_N(:) , 6,[] ) ;
    IDENT=diag([1,1,1,1/2,1/2,1/2]);
    tilde_E=IDENT*E;                 % deviatoric part of E
    z=max(0,sum(E.*tilde_E));        % scalar product of the deviatoric strain
    r=sqrt(z);                       % norm of the deviatoric strain
    a=mu_infty+(mu_0-mu_infty).*(1+lambda.*(r.^2)).^(p/2-1);
    a_node=VIZ.transformation(a,ELEM,WEIGHT);
    VIZ.draw_quantity(COORD,SURF,0*U_N,a_node,elem_type,size_xy_0,size_xy_L,size_z);
    colorbar off; colorbar('location','eastoutside')

    % values of div u
    div=E(1,:)+E(2,:)+E(3,:);
    div_node=VIZ.transformation(div,ELEM,WEIGHT);
    VIZ.draw_quantity(COORD,SURF,0*U_N,div_node,elem_type,size_xy_0,size_xy_L,size_z);
    colorbar off; colorbar('location','eastoutside')

    % convergence of the Newton-like solvers
    VIZ.figure_convergence(1:it_N, crit_hist_N, ...
        1:it_qN1, crit_hist_qN1, ...
        1:it_qN2, crit_hist_qN2)

    % line search coefficients
    VIZ.figure_omega(1:it_qN1-1, omega_hist_qN1,...
        1:it_qN2-1, omega_hist_qN2)
end
