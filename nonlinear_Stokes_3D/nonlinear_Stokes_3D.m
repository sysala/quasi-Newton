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

% if there is AGMG solver present in path (function agmg(...))
% "Notay, Y. AGMG software and documentation. Available at http://agmg.eu"
% Note: academic licence is free upon an email request
global AGMG_present
AGMG_present = 1;

%
% Main input data 
% 

  % Element, mesh and geometrical data 
  elem_type='P1'; % available choices of finite elements: 'P1', 'P2'
                  % for P2 elements, it is necessary to change
                  % the regularization parameter within inner solvers
  density=12;      % a positive integer defining mesh density in x-direction
  size_xy_0 = 6;  % size of the body in directions x and y on the left
  size_xy_L = 5;  % size of the body in directions x and y on the right
  size_z = 50;      % size of the body in z-direction  
  % It is assumed that size_xy_L < size_xy_0 < size_z.

  % Material parameters
  mu_0 = 1;                         
  mu_infty = 1e-3;   
  lambda = 10;
  p = 1.1;

  % Constant volume force representing pore pressure gradient
  gamma = 0.008 ; 

  % Boundary conditions - two available choices:
  % 'BC1' - zero velocity field in the normal direction on the shell
  % 'BC1' - zero velocity field in all directions on the shell
  boundary='BC1';
  u_z = 5 ;  % prescribed velocity in the direction z and at the prismatic centre 

%
% Mesh generation
%
  switch(elem_type)
    case 'P1'
        [COORD,ELEM,SURF1,SURF2,SURF3,SURF4,SURF5,SURF6]=...
                               mesh_P1(density,size_xy_0,size_xy_L,size_z);
        fprintf('P1 elements: \n')
    case 'P2'
        [COORD,ELEM,SURF1,SURF2,SURF3,SURF4,SURF5,SURF6]=...
                               mesh_P2(density,size_xy_0,size_xy_L,size_z);
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
        [U_3,Q]=boundary_BC1(COORD,SURF1,SURF3,SURF4,SURF5,SURF6);
    case 'BC2'
        [U_3,Q]=boundary_BC2(COORD,SURF1,SURF3,SURF4,SURF5,SURF6,size_xy_0);
    otherwise
        disp('bad choice of element type');
  end        
  U_3=u_z*U_3;

%
% Data from the reference element
%

% quadrature points and weights for volume and surface integration
[Xi, WF] = quadrature_volume(elem_type);
[Xi_s, WF_s] = quadrature_surface(elem_type);

% local basis functions and their derivatives for volume and surface
[HatP, DHatP1, DHatP2, DHatP3] = local_basis_volume(elem_type, Xi);
[HatP_s, DHatP1_s, DHatP2_s] = local_basis_surface(elem_type, Xi_s);

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
  [B,K_elast,WEIGHT]=auxiliary_matrices(ELEM,COORD,mu_0,...
                                                  DHatP1,DHatP2,DHatP3,WF);  

%
% Assembling of the vector of volume forces.
% 
  % Volume forces at integration points: f_V_int has size (3, n_int).
  f_V_int = [zeros(1, n_int); zeros(1, n_int); -gamma];
  % Compute the volume force vector.
  f = vector_volume_3D(ELEM, COORD, f_V_int, HatP, WEIGHT);

%
% Newton's and quasi-Newton's solvers
%    

  % initialization displacement
  U_it = [zeros(1,n_n); zeros(1,n_n); U_3] ;     

  % standard Newton method
  tic; 
  [U_N, it_N, crit_hist_N]=newton(U_it,WEIGHT,B,f,Q,mu_0,mu_infty,lambda,p);
  time_Newton=toc;
  fprintf("     solver's runtime:  " + time_Newton + "-----\n");

  % quasi-Newton method - preconditioner 1 (simplified stiffness matrix with variable coefficients)
  tic; 
  [U_qN1, it_qN1, crit_hist_qN1, omega_hist_qN1]=newton_quasi1(U_it,WEIGHT,K_elast,B,f,Q,mu_0,mu_infty,lambda,p);
  time_qNewton1=toc;
  fprintf("     solver's runtime:  " + time_qNewton1 + "-----\n");

  % quasi-Newton method - preconditioner 2 (simplified stiffness matrix with fixed coefficients)
  tic; 
  [U_qN2, it_qN2, crit_hist_qN2, omega_hist_qN2]=newton_quasi2(U_it,WEIGHT,K_elast,B,f,Q,mu_0,mu_infty,lambda,p);
  time_qNewton2=toc;
  fprintf("     solver's runtime:  " + time_qNewton2 + "-----\n");

%  
% Visualization of selected results
%

  % mesh  
  if density<5
   draw_mesh(COORD,SURF,elem_type) 
  end

  % total displacements + deformed shape - newton
  U_total = sqrt(U_N(1,:).^2 + U_N(2,:).^2 + U_N(3,:).^2);
  draw_quantity(COORD,SURF,1*U_N,U_total,elem_type,size_xy_0,size_xy_L,size_z) 

  % U_total = sqrt(U_qN1(1,:).^2 + U_qN1(2,:).^2 + U_qN1(3,:).^2);
  % draw_quantity(COORD,SURF,1*U_qN1,U_total,elem_type,size_xy_0,size_xy_L,size_z) 
  % U_total = sqrt(U_qN2(1,:).^2 + U_qN2(2,:).^2 + U_qN2(3,:).^2);
  % draw_quantity(COORD,SURF,1*U_qN2,U_total,elem_type,size_xy_0,size_xy_L,size_z) 

  % values of the function a
  E= reshape( B*U_N(:) , 6,[] ) ;
  IDENT=diag([1,1,1,1/2,1/2,1/2]); 
  tilde_E=IDENT*E;                 % deviatoric part of E
  z=max(0,sum(E.*tilde_E));        % scalar product of the deviatoric strain
  r=sqrt(z);                       % norm of the deviatoric strain
  a=mu_infty+(mu_0-mu_infty).*(1+lambda.*(r.^2)).^(p/2-1);
  a_node=transformation(a,ELEM,WEIGHT);
  draw_quantity(COORD,SURF,0*U_N,a_node,elem_type,size_xy_0,size_xy_L,size_z);
  colorbar off; colorbar('location','eastoutside')  

  % values of div u
  div=E(1,:)+E(2,:)+E(3,:);
  div_node=transformation(div,ELEM,WEIGHT);
  draw_quantity(COORD,SURF,0*U_N,div_node,elem_type,size_xy_0,size_xy_L,size_z);
  colorbar off; colorbar('location','eastoutside')  

  % convergence of the Newton-like solvers 
  figure_convergence(1:it_N, crit_hist_N, ...
                     1:it_qN1, crit_hist_qN1, ...
                     1:it_qN2, crit_hist_qN2)

  % line search coefficients
  figure_omega(1:it_qN1-1, omega_hist_qN1,...
               1:it_qN2-1, omega_hist_qN2)