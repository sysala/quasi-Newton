
function [S,m,M]=constitutive_problem_quasi2(E,mu_0,mu_infty,lambda,p)

% =========================================================================
%
% The aim of this function is to construct constitutive operator and
% auxiliary arrays for the quasi-Newton method 2 at integration points 
% 1,2,...,n_int. 
% 
% Input data:
%  E       - current strain tensor, size(E)=(3,n_int)
%  mu_0,mu_infty,lambda,p - material parameters
%
% Output data:
%  S      - stress tensors at integration points, size(S)=(4,n_int)
%  M,m    - parameters defining the quasi-Newton step
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
% Values of the functions a and a_bar at any integration point
%
  a=mu_infty+(mu_0-mu_infty).*(1+lambda.*(r.^2)).^(p/2-1);
  a_bar=mu_infty+(mu_0-mu_infty).*(1+(p-1)*lambda.*(r.^2)).*(1+lambda.*(r.^2)).^(p/2-2);

%
% The stress tensor
%
  S=repmat(a,6,1).*tilde_E;   

%
% Parameters m and M defining the quasi-Newton step
%
  
  m=min(a_bar./mu_0);
  M=max(a./mu_0);
 
 end
