% ************************************************************************
%
%  constitutive_problem_quasi2   -  solution of the constitutive model
%                                     to the magnetic problem 
%
% ************************************************************************

function   [S,DS1,m,M]= constitutive_problem_quasi2 (E,heter_int,r_crit,alpha,beta)
% =========================================================================
%
% Input data:
%  E     - current values of the gradient at integration points, size(E)=(2,n_int)
%  heter_int - a logical array indicating the heterogeneity
%  r_crit - critical values of |\nabla U|
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

% approximation of the derivative of the constitutive operator and 
% parameters m and M defining the quasi-Newton step
  if nargout>1

    a1_der=8*beta*(1-alpha)*(r1.^7)./((r1.^8+beta).^2);
    a1_bar=a1+r1.*a1_der;
    gamma=0.5;  
    m=1; M=1;
    DS1=zeros(1,length(r1));
    %
    n_r=length(r_crit);
    if n_r==1
        lambda=min(a1); Lambda=max(a1_bar);
        c=(1-gamma)*lambda+gamma*Lambda;
        DS1(:)=c;
        m=lambda/c;
        M=Lambda/c;
    else
        for j=2:n_r
            Om=(r1>=r_crit(j-1))&(r1<r_crit(j));
            a11=a1(Om); a11_bar=a1_bar(Om);
            if size(a11,2)>0
              lambda=min(a11); Lambda=max(a11_bar);
              c=(1-gamma)*lambda+gamma*Lambda;
              DS1(Om)=c;            
              m=min(m,lambda/c);
              M=max(M,Lambda/c);
            end
        end
        Om=(r1>=r_crit(n_r));
        a11=a1(Om); a11_bar=a1_bar(Om);
        if size(a11,2)>0
              lambda=min(a11); Lambda=max(a11_bar);
              c=(1-gamma)*lambda+gamma*Lambda;
              DS1(Om)=c;            
              m=min(m,lambda/c);
              M=max(M,Lambda/c);
       end
    end  % if n_r

  end % nargout
  
end % function
