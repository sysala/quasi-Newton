% ************************************************************************
%
%  constitutive_problem   -  solution of the constitutive model
%                            to the magnetic problem 
%
% ************************************************************************

function   [S,DS1]= constitutive_problem (E,heter_int,alpha,beta)
% =========================================================================
%
% Input data:
%  E     - current values of the gradient at integration points, size(E)=(2,n_int)
%  heter_int - a logical array indicating the heterogeneity
%  alpha, beta - material parameters
%
% Output data:
%  S     - constitutive operator at integration points, size(S)=(2,n_int)
%  DS    - derivative of the constitutive operator at integr. points
%          corresponding to the domain Omega1
%
% =========================================================================                                   

% constitutive operator
  r=sqrt(E(1,:).^2+E(2,:).^2);
  r1=r(heter_int);
  a1=alpha+(1-alpha)*r1.^8./(r1.^8+beta);
  a=zeros(1,length(r));
  a(heter_int)=a1;
  a(~heter_int)=alpha;
  S=[a; a].*E;
  
% derivative of the constitutive operator   
  if nargout>1
    a_der_r=8*beta*(1-alpha)*(r1.^6)./((r1.^8+beta).^2);
    E1=E(:,heter_int);
    Mat=[E1(1,:).*E1(1,:)
       E1(2,:).*E1(1,:)
       E1(1,:).*E1(2,:)
       E1(2,:).*E1(2,:)];
    DS1=[1;0;0;1]*a1+Mat.*repmat(a_der_r,4,1);
  end

 end
