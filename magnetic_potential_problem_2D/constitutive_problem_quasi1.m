% ************************************************************************
%
%  constitutive_problem_quasi1   -  solution of the constitutive model
%                                    to the magnetic problem 
%
% ************************************************************************

function   [S,DS1,m,M]= constitutive_problem_quasi1 (E,heter_int,alpha,beta)
% =========================================================================
%
% Input data:
%  E     - current values of the gradient at integration points, size(E)=(2,n_int)
%  heter_int - a logical array indicating the heterogeneity
%  alpha, beta - material parameters
%
% Output data:
%  S      - projection of E at integration points, size(S)=(2,n_int)
%  DS1    - an array representing the quasi-Newton operator at integr. points,
%           corresponding to the domain Omega1
%  M,m    - parameters defining the quasi-Newton step
% =========================================================================                                   

% constitutive operator
  r=sqrt(E(1,:).^2+E(2,:).^2);
  r1=r(heter_int);
  a1=alpha+(1-alpha)*r1.^8./(r1.^8+beta);
  a=zeros(1,length(r));
  a(heter_int)=a1;
  a(~heter_int)=alpha;
  S=[a; a].*E;
  
  if nargout>1
% approximation of the derivative of the constitutive operator    
    a1_der=8*beta*(1-alpha)*(r1.^7)./((r1.^8+beta).^2);
    a1_bar=a1+r1.*a1_der;
    gamma=1/2;  
    DS1=(1-gamma)*a1+gamma*a1_bar;

% parameters m and M defining the quasi-Newton step
    rho=max(a1_bar./a1);
    if size(rho,2)>0
      m=1/(1-gamma+rho*gamma);
      M=rho/(1-gamma+rho*gamma);
    else
      m=1; M=1;
    end

  end
  
 end
