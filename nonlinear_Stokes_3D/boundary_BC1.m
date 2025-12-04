
function [U_3,Q]=boundary_BC1(coord,surf1,surf3,surf4,surf5,surf6)

% =========================================================================
%
%  This function creates tetrahedral mesh for P1 elements
%
%  input data:
%    coord   - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%              number of nodes
%    surf1    - array containing numbers of nodes defining each surface 
%               element on the bottom face;
%               size(surf)=(3,n_s1), n_s1 = number of surface elements
%    surf3    - array containing numbers of nodes defining each surface 
%               element on the front face;
%               size(surf)=(3,n_s3), n_s3 = number of surface elements
%    surf4    - array containing numbers of nodes defining each surface 
%               element on the right face;
%               size(surf)=(3,n_s4), n_s4 = number of surface elements
%    surf5    - array containing numbers of nodes defining each surface 
%               element on the back face;
%               size(surf)=(3,n_s5), n_s5 = number of surface elements
%    surf6    - array containing numbers of nodes defining each surface 
%               element on the left face;
%               size(surf)=(3,n_s6), n_s6 = number of surface elements
%
%  output data:
%    U_3     - initial displacements in the z-direction that catches
%              nonhomogeneous Dirichlet boundary conditions, size(U_z)=(1,n_n)
%    Q       - logical array indicating the nodes where the 
%              Dirichlet boundary condition is considered, size(Q)=(3,n_n)
%
% ======================================================================
%
 

%% Logical array indicating the nodes with the Dirichlet boundary cond.    
  Q = true(size(coord));
  Q(3,surf1(:)) = 0;  
  Q(2,surf3(:)) = 0;     
  Q(1,surf4(:)) = 0;   
  Q(2,surf5(:)) = 0;     
  Q(1,surf6(:)) = 0;   

%% Nonhomogeneous part of the velocity in the direction z
U_3=ones(1,size(coord,2));

end
