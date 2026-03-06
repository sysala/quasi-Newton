
function [coord,elem,surf1,surf2,surf3,surf4,surf5,surf6]=...
                              mesh_P1(density,size_xy_0,size_xy_L,size_z)

% =========================================================================
%
%  This function creates tetrahedral mesh for P1 elements
%
%  input data:
%    density   - an integer defining a density of a uniform mesh
%    size_xy_0 - size of the body in directions x and y on the left
%    size_xy_L - size of the body in directions x and y on the right
%    size_z  - size of the body in z-direction (integer) 
%    body=(0,size_xy)x(0,size_xy)x(0,size_z)
%
%  output data:
%    coord   - coordinates of the nodes, size(coord)=(3,n_n) where n_n is a
%              number of nodes
%    elem    - array containing numbers of nodes defining each element,
%              size(elem)=(4,n_e), n_e = number of elements
%    surf1    - array containing numbers of nodes defining each surface 
%               element on the bottom face;
%               size(surf)=(3,n_s1), n_s1 = number of surface elements
%    surf2    - array containing numbers of nodes defining each surface 
%               element on the top face;
%               size(surf)=(3,n_s2), n_s2 = number of surface elements
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
% ======================================================================
%

%
% numbers of segments, nodes and elements
%
  
  N_x = density;              % number of segments in x direction
  N_y = density;              % number of segments in y direction
  N_z = round((size_z/size_xy_0)*density); 
                              % number of segments in z direction
 
  % 
  n_cell_xy = N_x*N_y;          % number of cells in xy plane
  n_e = n_cell_xy*N_z*6;        % total number of elements

%
% C - 3D auxilliary array that contains node numbers and that is important 
% for the mesh construction. 
%
  C=reshape(1:(N_x+1)*(N_y+1)*(N_z+1),N_x+1,N_y+1,N_z+1); 
  
%
% coordinates of nodes
%
  % coordinates in directions x, y and z (before transformation)
  coord_x=linspace(-size_xy_0/2,size_xy_0/2,N_x+1);
  coord_y=linspace(-size_xy_0/2,size_xy_0/2,N_y+1);
  coord_z=linspace(0,size_z,N_z+1);
  
  % long 1D arrays containing coordinates of all nodes in x,y,z directions
  c_x=repmat(coord_x,1,(N_y+1)*(N_z+1));     
  c_y=repmat(kron(coord_y,ones(1,N_x+1)),1,N_z+1);     
  c_z=kron(coord_z,ones(1,(N_x+1)*(N_y+1)));    
         
  % the required array of coordinates, size(coord)=(3,n_n)
  coord=[c_x.*(1-c_z*(1-size_xy_L/size_xy_0)/size_z) 
         c_y.*(1-c_z*(1-size_xy_L/size_xy_0)/size_z)
         c_z] ;  

% 
% construction of the array elem
%
  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0 0], V2 -> [1 0 0], V3 -> [1 1 0], V4 -> [0 1 0]
  %  V5 -> [0 0 1], V6 -> [1 0 1], V7 -> [1 1 1], V8 -> [0 1 1]
  %  V1,...,V8 are logical 3D arrays which enable to select appropriate
  %  nodes from the array C.

  V1=false(N_x+1,N_y+1,N_z+1);
  V1(1:N_x,1:N_y,1:N_z)=1;
  %
  V2=false(N_x+1,N_y+1,N_z+1);
  V2(2:(N_x+1),1:N_y,1:N_z)=1;
  %
  V3=false(N_x+1,N_y+1,N_z+1);
  V3(2:(N_x+1),2:(N_y+1),1:N_z)=1;
  %
  V4=false(N_x+1,N_y+1,N_z+1);
  V4(1:N_x,2:(N_y+1),1:N_z)=1;
  %
  V5=false(N_x+1,N_y+1,N_z+1);
  V5(1:N_x,1:N_y,2:(N_z+1))=1;
  %
  V6=false(N_x+1,N_y+1,N_z+1);
  V6(2:(N_x+1),1:N_y,2:(N_z+1))=1;
  %
  V7=false(N_x+1,N_y+1,N_z+1);
  V7(2:(N_x+1),2:(N_y+1),2:(N_z+1))=1;
  %
  V8=false(N_x+1,N_y+1,N_z+1);
  V8(1:N_x,2:(N_y+1),2:(N_z+1))=1;
 
  % used division of a prism into 6 tetrahedrons:   
  %   V1 V2 V4 V6
  %   V1 V4 V5 V6
  %   V4 V5 V6 V8
  %   V2 V3 V4 V6
  %   V3 V6 V7 V4
  %   V4 V6 V7 V8
  % size(aux_elem)=(6*4,n_e/6)
  aux_elem=[C(V1)'; C(V2)'; C(V4)'; C(V6)';
            C(V1)'; C(V4)'; C(V5)'; C(V6)';
            C(V4)'; C(V5)'; C(V6)'; C(V8)';
            C(V2)'; C(V3)'; C(V4)'; C(V6)';
            C(V3)'; C(V6)'; C(V7)'; C(V4)';
            C(V4)'; C(V6)'; C(V7)'; C(V8)' ];
        
  % the array elem, size(elem)=(4,n_e)
  elem=reshape(aux_elem,4,n_e);     
 
%
% Surfaces of the body - the array "surf1,...,surf6"
  
  % For each face of the body, we define the restriction C_s of the array C 
  % and logical 2D arrays V1_s,...,V4_s which enable to select appropriate
  % nodes from the array C_s. We consider the following ordering of the 
  % nodes within the unit square:
  %   V1_s -> [0 0], V2_s -> [1 0], V3_s -> [1 1], V4_s -> [0 1]
  % Finally, we use the division of a rectangle into 2 triangles which is
  % in accordance to the division of a prism into 6 tetrahedrons, see above.

  % Face 1: z=0 (the bottom of the body)
  C_s=zeros(N_x+1,N_y+1);
  C_s(:,:)=C(:,:,1);  
  V1_s=false(N_x+1,N_y+1);  V1_s(1:N_x    ,1:N_y    )=1;
  V2_s=false(N_x+1,N_y+1);  V2_s(2:(N_x+1),1:N_y    )=1;
  V3_s=false(N_x+1,N_y+1);  V3_s(2:(N_x+1),2:(N_y+1))=1;
  V4_s=false(N_x+1,N_y+1);  V4_s(1:N_x    ,2:(N_y+1))=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)' ];
  surf1=reshape(aux_surf,3,2*N_x*N_y);       
  
  % Face 2: z=size_z (the top of the body)
  C_s=zeros(N_x+1,N_y+1);
  C_s(:,:)=C(:,:,end);  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)' ];
  surf2=reshape(aux_surf,3,2*N_x*N_y);       
  
  % Face 3: y=-size_xy/2 (the front of the body)
  C_s=zeros(N_x+1,N_z+1);
  C_s(:,:)=C(:,1,:);  
  V1_s=false(N_x+1,N_z+1);  V1_s(1:N_x    ,1:N_z    )=1;
  V2_s=false(N_x+1,N_z+1);  V2_s(2:(N_x+1),1:N_z    )=1;
  V3_s=false(N_x+1,N_z+1);  V3_s(2:(N_x+1),2:(N_z+1))=1;
  V4_s=false(N_x+1,N_z+1);  V4_s(1:N_x    ,2:(N_z+1))=1;  
  aux_surf=[C_s(V4_s)'; C_s(V1_s)'; C_s(V3_s)';
            C_s(V2_s)'; C_s(V3_s)'; C_s(V1_s)' ];
  surf3=reshape(aux_surf,3,2*N_x*N_z);       
  
  % Face 4: x=size_xy/2 (the right hand side of the body)
  C_s=zeros(N_y+1,N_z+1);
  C_s(:,:)=C(end,:,:);  
  V1_s=false(N_y+1,N_z+1);  V1_s(1:N_y    ,1:N_z    )=1;
  V2_s=false(N_y+1,N_z+1);  V2_s(2:(N_y+1),1:N_z    )=1;
  V3_s=false(N_y+1,N_z+1);  V3_s(2:(N_y+1),2:(N_z+1))=1;
  V4_s=false(N_y+1,N_z+1);  V4_s(1:N_y    ,2:(N_z+1))=1;  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)' ];
  surf4=reshape(aux_surf,3,2*N_y*N_z);       
  
  % Face 5: y=size_xy/2 (the back of the body)
  C_s=zeros(N_x+1,N_z+1);
  C_s(:,:)=C(:,end,:);  
  V1_s=false(N_x+1,N_z+1);  V1_s(1:N_x    ,1:N_z)=1;
  V2_s=false(N_x+1,N_z+1);  V2_s(2:(N_x+1),1:N_z)=1;
  V3_s=false(N_x+1,N_z+1);  V3_s(2:(N_x+1),2:(N_z+1))=1;
  V4_s=false(N_x+1,N_z+1);  V4_s(1:N_x    ,2:(N_z+1))=1;  
  aux_surf=[C_s(V4_s)'; C_s(V1_s)'; C_s(V3_s)';
            C_s(V2_s)'; C_s(V3_s)'; C_s(V1_s)' ];
  surf5=reshape(aux_surf,3,2*N_x*N_z);       
  
  % Face 6: x=-size_xy/2 (the left hand side of the body)
  C_s=zeros(N_y+1,N_z+1);
  C_s(:,:)=C(1,:,:);  
  V1_s=false(N_y+1,N_z+1);  V1_s(1:N_y    ,1:N_z    )=1;
  V2_s=false(N_y+1,N_z+1);  V2_s(2:(N_y+1),1:N_z    )=1;
  V3_s=false(N_y+1,N_z+1);  V3_s(2:(N_y+1),2:(N_z+1))=1;
  V4_s=false(N_y+1,N_z+1);  V4_s(1:N_y    ,2:(N_z+1))=1;  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)'; C_s(V4_s)';
            C_s(V3_s)'; C_s(V4_s)'; C_s(V2_s)' ];
  surf6=reshape(aux_surf,3,2*N_y*N_z);      
    
end
