
function  [K,B,WEIGHT]=stiffness_matrix(ELEM,COORD,heter_int,Q,DHatP1,DHatP2,WF,alpha)
 
% =========================================================================
%
%  This function assembles auxilliary arrays which are important for 
%  assembling of linearized matrices.
%
%  input data:
%    ELEM - n_p x n_e array containing numbers of nodes defining each
%           element, n_e = number of elements
%    COORD - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    heter_int - a logical array indicating the heterogeneity
%    Q -     logical array indicating the nodes with the Dirichlet boundary 
%            conditions. size(Q)=(1,n_n)
%    DHatP1 - derivatives of basis functions at the quadrature points 
%             in the direction xi_1, size(DHatP1)=(n_p,n_q)
%    DHatP2 - derivatives of basis functions at the quadrature points 
%             in the direction xi_2, size(DHatP2)=(n_p,n_q)
%        WF - weight factors on the reference element, size(WF)=(1,n_q)
%
%  output data:
%    K_fix - the linearized matrix representing the subdomain Omega_2, size(K)=(n_unkown, n_uknown)
%    B     - the matrix transforming nodal values to values at integration
%            points, size(B)=(2*n_int,n_unknown)
%    WEIGHT - the weight coefficients for each quadrature point, 
%             size(WEIGHT)=(1, n_int)
%
% ======================================================================
%

%
% auxilliary notation
%

  n_n=size(COORD,2);    % number of nodes including midpoints
  n_e=size(ELEM,2);     % number of elements
  n_p=size(ELEM,1);     % number of vertices per element
  n_q=length(WF);       % number of quadrature points
  n_int = n_e*n_q ;     % total number of integrations points
     
%
% Jacobian, its determinant and inverse, 
% derivative of local basis functions
% 

  % extension of the input arrays DHatP1,DHatP2 by replication
  % size(DHatPhi1)=size(DHatPhi2)=(n_p, n_int)
  DHatPhi1=repmat(DHatP1,1,n_e);
  DHatPhi2=repmat(DHatP2,1,n_e);
  
  % coordinates of nodes defining each element
  % size(COORDe1)=size(COORDe2)=(n_p, n_e)
  COORDe1=reshape(COORD(1,ELEM(:)),n_p,n_e);
  COORDe2=reshape(COORD(2,ELEM(:)),n_p,n_e);
  
  % coordinates of nodes around each integration point
  % size(COORDint1)=size(COORDint2)=(n_p, n_int)
  COORDint1=kron(COORDe1,ones(1,n_q)); 
  COORDint2=kron(COORDe2,ones(1,n_q)); 
  
  % components of the Jacobians: size=(1, n_int)
  J11=sum(COORDint1.*DHatPhi1); J12=sum(COORDint2.*DHatPhi1); 
  J21=sum(COORDint1.*DHatPhi2); J22=sum(COORDint2.*DHatPhi2); 
  
  % determinant of the Jacobian: size=(1, n_int)
  DET=J11.*J22 - J12.*J21;
  
  % components of the inverse to the Jacobian: size=(1, n_int)
  Jinv11 =  J22./DET; Jinv12 = -J12./DET;
  Jinv21 = -J21./DET; Jinv22 =  J11./DET; 
  
  DET=abs(DET);
  
  % derivatives of local basis functions w.r.t the coordinates x_1,x_2:
  % size(DPhi1)=size(DPhi2)=(n_p, n_int)
  DPhi1 = repmat(Jinv11,n_p,1).*DHatPhi1 + repmat(Jinv12,n_p,1).*DHatPhi2;
  DPhi2 = repmat(Jinv21,n_p,1).*DHatPhi1 + repmat(Jinv22,n_p,1).*DHatPhi2;  

%
% assembling of the matrix B
% size(B)=(2*n_int, n_unknown)
%   
  % nonzero values of B
  n_b = 2*n_p ;
  vB=zeros(n_b,n_int);        % size(vB)=(2*n_p,n_int)
  vB(1:2:n_b-1,:)=DPhi1; 
  vB(2:2:n_b  ,:)=DPhi2; 
  % i-th and j-th indices of B: size(iB)=size(jB)=(2*n_p,n_int)
  AUX=reshape(1:2*n_int,2,n_int);
  iB=repmat(AUX,n_p,1);  
  jB=kron(ELEM,ones(2,n_q));  
  % the sparse matrix B
  B=sparse(iB(:),jB(:),vB(:), 2*n_int,n_n);  
  B=B(:,Q);

%  
% the weight coefficients, size(WEIGHT)=(1, n_int) 
%
  WEIGHT = DET.*repmat(WF,1,n_e);
  
%  
% assembling of the constitutive matrix D
% size(D)=(2*n_int, 2*n_int)
%
  AUX2=AUX(:,~heter_int);
  vD=repmat(alpha*WEIGHT(~heter_int),2,1); 
  D=sparse(AUX2(:),AUX2(:),vD(:),2*n_int, 2*n_int);

%
% elastic stiffness matrix: size(K)=(n_unkown, n_uknown)
%
  K = B'*D*B ;    
 
end  % end of function