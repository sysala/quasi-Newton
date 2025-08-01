
% =========================================================================
%
%  Magnetic potential problem in 2D 
%
%  This program solves the magnetic potential problem in 2D. The problem is
%  discretized by P1 or P2 elements. In the latter case, the 7-point Gauss 
%  quadrature is used. The aim is to compare the Newton method with the
%  quasi-Newton/variable preconditioning method. All Newton-like solvers are 
%  considered either without line search or with line search. More details
%  can be found in the report "magnetic_problem_remark.pdf".
% 
% ======================================================================
%

%
% The main input data 
%

  % elem_type - the type of finite elements; available choices: 'P1', 'P2'
  elem_type='P1';
  
  % mesh data
  level=3;       % density of the initial mesh

  % geometrical parameters (choose only integers)
  size_x = 1;  % size of the rectangular body in x-direction
  size_y = 1;  % size of the rectangular body in y-direction
  delta = 1;   % parameter defining splitting of the computational domain
               % on subdomains with different material parameters
               % delta=0 --> no heterogeneity Ω1=Ω
               % delta=1 --> Ω1 = [0.1, 0.9] × [0.1, 0.9]
               % delta=2 --> Ω1 = [0.2, 0.8] × [0.2, 0.8]

  % material parameters
  alpha=0.0003;     
  beta=16000;  
  rho=40;
%   rho=4*10^(6);    % electric density (right hand side)
  kappa=25;
  r_crit=a_function(kappa,alpha,beta);
 
%
% Data from the reference element
%
  
  % quadrature points and weights for volume and surface integration
  [Xi, WF] = quadrature_volume_2D(elem_type);
  [Xi_s, WF_s] = quadrature_surface(elem_type);
  
  % local basis functions and their derivatives for volume and surface
  [HatP,DHatP1,DHatP2] = local_basis_volume_2D(elem_type, Xi);
  [HatP_s,DHatP1_s] = local_basis_surface(elem_type, Xi_s);    
  
%
% Finite element mesh
%  
    
  % mesh generation
  switch(elem_type)
    case 'P1'
        [COORD,ELEM,heter,Q]=mesh_P1(level,size_x,size_y,delta);
        fprintf('P1 elements: \n')
    case 'P2'
        [COORD,ELEM,heter,Q]=mesh_P2(level,size_x,size_y,delta);
        fprintf('P2 elements: \n')
    otherwise
        disp('bad choice of element type');
  end     
%   draw_zones(COORD,ELEM,heter)

  % number of nodes, elements and integration points + print
  n_n=size(COORD,2);            % number of nodes
  n_unknown=length(COORD(1,Q)); % number of unknowns
  n_e=size(ELEM,2);             % number of elements
  n_q=length(WF);               % number of quadratic points
  n_int = n_e*n_q ;             % total number of integrations points 
  % 
  fprintf('\n');   
  fprintf('The coarsest mesh:'); 
  fprintf('  number of nodes =%d ',n_n);  
  fprintf('  number of unknowns =%d ',n_unknown);
  fprintf('  number of elements =%d ',n_e);
  fprintf('  number of integration points =%d ',n_int);
  fprintf('\n'); 
  
%
% Assembling of the stiffness matrix and the load vector
%
  heter1=repmat(heter,n_q,1);
  heter_int=heter1(:);   

  % assembling of auxiliary matrices
  [K_fix,B,WEIGHT]=stiffness_matrix(ELEM,COORD,heter_int,Q,DHatP1,DHatP2,WF,alpha);  
  
  % assembling of the right hand side vector (only volume forces)
  f_V=vector_volume(ELEM,COORD,HatP,WEIGHT);
  f=rho*f_V(Q);

% %
% % Linear solver of the problem on the last line on page 6 in the report
% %
%   AUX=reshape(1:2*n_int,2,n_int);
%   AUX1=AUX(:,heter_int);
%   WEIGHT1=WEIGHT(heter_int);
%   vD = repmat(WEIGHT1,2,1) ;   
%   D_p=sparse(AUX1(:),AUX1(:),vD(:),2*n_int, 2*n_int);
%   K = K_fix+B'*D_p*B;  
%   U0d = K\(f); 

%
% Newton-like solvers
%

  % the Newton method
  tic;
  U_ini=zeros(n_unknown,1);
  [U0, it0, crit_hist0]=newton(U_ini,WEIGHT,K_fix,B,f,heter_int,alpha,beta);
  assembly_time=toc; 
  fprintf('time spent on K: %e seconds, ',assembly_time(end));
  fprintf('\n');  

  % the Newton method with line search
  tic;
  U_ini=zeros(n_unknown,1);
  [U0d, it0d, crit_hist0d, omega_hist0d]=newton_damped(U_ini,WEIGHT,K_fix,B,f,heter_int,alpha,beta);
  assembly_time=toc; 
  fprintf('time spent on K: %e seconds, ',assembly_time(end));
  fprintf('\n');  

  % the Quasi-Newton method 1
  tic;
  U_ini=zeros(n_unknown,1);
  [U1, it1, crit_hist1, omega_hist1]=newton_quasi1(U_ini,WEIGHT,K_fix,B,f,heter_int,alpha,beta);
  assembly_time=toc; 
  fprintf('time spent on K: %e seconds, ',assembly_time(end));
  fprintf('\n');  
 
  % the Quasi-Newton method 1 with line search
  tic;
  U_ini=zeros(n_unknown,1);
  [U1d, it1d, crit_hist1d, omega_hist1d]=newton_quasi1_damped(U_ini,WEIGHT,K_fix,B,f,heter_int,alpha,beta);
  assembly_time=toc; 
  fprintf('time spent on K: %e seconds, ',assembly_time(end));
  fprintf('\n');  

  % the Quasi-Newton method 2
  tic;
  U_ini=zeros(n_unknown,1);
  [U2, it2, crit_hist2, omega_hist2]=newton_quasi2(U_ini,WEIGHT,K_fix,B,f,heter_int,r_crit,alpha,beta);
  assembly_time=toc; 
  fprintf('time spent on K: %e seconds, ',assembly_time(end));
  fprintf('\n');  

  % the Quasi-Newton method 2 with line search
  tic;
  U_ini=zeros(n_unknown,1);
  [U2d, it2d, crit_hist2d, omega_hist2d]=newton_quasi2_damped(U_ini,WEIGHT,K_fix,B,f,heter_int,r_crit,alpha,beta);
  assembly_time=toc; 
  fprintf('time spent on K: %e seconds, ',assembly_time(end));
  fprintf('\n');  


%  
% Postprocessing - visualization of selected results
%      
  
%   % mesh visualization
%   draw_mesh(COORD,ELEM,size_x,size_y) 

  % visualization of the solution
  U0_full=zeros(n_n,1);
  U0_full(Q)=U0d;
  draw_quantity(COORD,ELEM,0*U0_full,U0_full',size_x,size_y)
  %
  figure;
  hold on;
  trisurf(ELEM(1:3,:)',COORD(1,:)',COORD(2,:)' ,U0_full,'EdgeColor','none');
  hold off;

  % visualization of the gradient to U
  E = zeros(2,n_int);  % values of the gradient at integration points
  E(:) = B*U0d(:) ;   % strain at integration points
  r=sqrt(E(1,:).^2+E(2,:).^2);
  r_node=transformation(r,ELEM,WEIGHT);
  draw_quantity(COORD,ELEM,0*U0_full,r_node,size_x,size_y);

  % convergence evolution for algorithm without line search
  figure;
  hold on;
  plot(1:it0, crit_hist0)
  plot(1:it1, crit_hist1)
  plot(1:it2, crit_hist2)
  hold off;

  % convergence evolution for algorithm with line search
  figure;
  hold on;
  plot(1:it0d, crit_hist0d)
  plot(1:it1d, crit_hist1d)
  plot(1:it2d, crit_hist2d)
  hold off;

  % convergence evolution for algorithm with line search
  figure;
  hold on;
  plot(1:it1-1, omega_hist1)
  plot(1:it2-1, omega_hist2)
  hold off;

  % convergence evolution for algorithm with line search
  figure;
  hold on;
  plot(1:it0d-1, omega_hist0d)
  plot(1:it1d-1, omega_hist1d)
  plot(1:it2d-1, omega_hist2d)
  hold off;

 
 