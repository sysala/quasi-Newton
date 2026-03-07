
function [coord,elem,heter,Q]=mesh_P1(level,size_x,size_y,delta)

% =========================================================================
%
%  This function creates a triangular mesh for P1 elements and specify
%  material heterogeneity
%
%  input data:
%    level   - an integer defining a density of a uniform mesh
%    size_x  - size of the body in x-direction (integer)
%    size_y  - size of the body in y-direction (integer)
%    body=(0,size_x)  x (0,size_y)
%    delta   - parameter defining splitting of the computational domain
%              on subdomains with different material parameters
%              delta=0 --> no heterogeneity Ω1=Ω
%              delta=1 --> Ω1 = [0.1, 0.9] × [0.1, 0.9]
%              delta=2 --> Ω1 = [0.2, 0.8] × [0.2, 0.8]
%
%  output data:
%    coord    - coordinates of the nodes, size(coord)=(2,n_n) where n_n is
%               a number of nodes
%    elem     - array containing numbers of nodes defining each element, 
%               size(elem)=(3,n_e), n_e = number of elements (triangles)
%    heter    - logical array describing material heterogeneity, size(heter)=(1,n_e)
%    Q        - logical array indicating the nodes where the homogeneous
%               Dirichlet boundary condition is considered, size(Q)=(n_n,1)
%
% ======================================================================
%

%
% Numbers of segments, nodes and elements
%

  N_x = 10*round(size_x)*2^level;      % number of segments in x direction
  N_y = 10*round(size_y)*2^level;      % number of segments in y direction
  n_e = 2*N_x*N_y    ;       % total number of elements

%
% Coordinates of nodes - the array coord
%
  
  % coordinates in directions x and y
  coord_x=linspace(0,size_x,N_x+1);
  coord_y=linspace(0,size_y,N_y+1);
  
  % long 1D arrays containing coordinates of all nodes in x,y directions
  c_x=repmat(coord_x,1,N_y+1);     
  c_y=repmat(kron(coord_y,ones(1,N_x+1)),1);  
         
  % the required array of coordinates, size(coord)=(2,n_n)    
  coord=[c_x; c_y] ;  

% 
% The array elem
%
  
  % V - 2D auxilliary array containing node numbers used for 
  %     the construction of the array elem
  V=reshape(1:(N_x+1)*(N_y+1),N_x+1,N_y+1);

  %  V1,...,V4 are 1 x (Nx*Ny) arrays containing numbers of nodes for each
  %  square cell of the triangulation. Location of the nodes is explained
  %  on the unit square:
  %  V1: [0 0], V2: [1 0], V3: [1 1], V4: [0,1]
  V1=V(1: N_x   ,1: N_y   ); V1=V1(:)';
  V2=V(2:(N_x+1),1: N_y   ); V2=V2(:)';
  V3=V(2:(N_x+1),2:(N_y+1)); V3=V3(:)';
  V4=V(1: N_x   ,2:(N_y+1)); V4=V4(:)';
  
  % used division of a square cell into 2 triangles:   
  %   V2 V3 V1
  %   V4 V1 V3 
  aux_elem=[V2; V3; V1;
            V4; V1; V3 ];
  elem=reshape(aux_elem,3,n_e);     

%
% The array heter
%  
  N1_x = delta*round(size_x)*2^level;     
  N1_y = delta*round(size_y)*2^level;    
  N2_x = N_x-N1_x;     
  N2_y = N_y-N1_y;     
  R=false(N_x,N_y);
  R(N1_x+1:N2_x,N1_y+1:N2_y)=true;
  R1=R(:)'; R2=[R1; R1];
  heter=R2(:)';

%
% Boundary conditions
%   

  % Dirichlet boundary conditions are represented by logical array Q which
  % restricts the corresponding degrees of freedom; size(Q)=(n_n,1)
  n_n=size(coord,2);
  Q=true(n_n,1);
  Q((coord(1,:)==size_x)) = 0;   
  Q((coord(1,:)==0)) = 0;   
  Q((coord(2,:)==size_y)) = 0;   
  Q((coord(2,:)==0)) = 0;   
  
end
