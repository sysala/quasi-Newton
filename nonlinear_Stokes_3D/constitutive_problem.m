
function [S,DS]=constitutive_problem(E,mu_0,mu_infty,lambda,p)

% =========================================================================
%
% The aim of this function is to construct constitutive and consistent 
% tangent operators at integration points 1,2,...,n_int for the semismooth 
% Newton method.
%
% Input data:
%  E       - current strain tensor, size(E)=(6,n_int)
%  mu_0,mu_infty,lambda,p - material parameters
%
% Output data:
%  S      - stress tensors at integration points, size(S)=(6,n_int)
%  DS     - consistent tangent operators at integr. points,
%           size(DS)=(36,n_int)
%
% =========================================================================
%

%
% Norm of the strain tensor ae any integration point 
%
  IDENT=diag([1,1,1,1/2,1/2,1/2]); 
  tilde_E=IDENT*E;                 % deviatoric part of E
  z=max(0,sum(E.*tilde_E));        % scalar product of the deviatoric strain
  r=sqrt(z);                       % norm of the deviatoric strain
  
%  
% Values of the function a its derivative at any integration point
%
  a=mu_infty+(mu_0-mu_infty).*(1+lambda.*(r.^2)).^(p/2-1);
  a_der=(p-2)*(mu_0-mu_infty).*lambda.*r.*(1+lambda.*(r.^2)).^(p/2-2);

%
% The stress tensor
%
  S=repmat(a,6,1).*tilde_E;   

%
% The tangent stiffness matrix
%
  DS=IDENT(:)*a;
  IND=(r>0);
  N_hat=tilde_E(:,IND)./repmat(r(IND),6,1);
  NN_hat=repmat(N_hat,6,1).*kron(N_hat,ones(6,1));
  DS(:,IND)=DS(:,IND)+repmat(r(IND).*a_der(IND),36,1).*NN_hat; 

 end
