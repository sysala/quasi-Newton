function [U, it, crit_hist]=newton(U_ini,WEIGHT,K_fix,B,f,heter_int,alpha,beta)
                               
%--------------------------------------------------------------------------
% The Newton method for solution of the system 
%              find U:   F(U)=f
%
% Input data:
%   U_ini - initial choice of U
%   WEIGHT - weight coefficients of integration points, size(WEIGHT)=(1,n_int)
%   K_fix - the linearized matrix representing the subdomain Omega_2, size(K)=(n_unkown, n_uknown)
%   B - the gradient matrix, size(B)=(2*n_int,n_uknown)
%   f - vector of external forces, size(f)=(n_uknown,1)
%   heter_int - a logical array indicating the heterogeneity
%   alpha, beta - material parameters
%
% Output data:
%   U - approximation of the solution, size(U)=(n_uknown,1)
%   it - number of Newton's iteration
%   crit_hist - evaluation of the stopping criterion
%
%--------------------------------------------------------------------------

%
% Initialization 
%

  it_max=100;
  tol=1e-9;
  crit_hist=zeros(1,it_max);
  n_int=length(WEIGHT);
  E = zeros(2,n_int);  % values of the gradient at integration points               
  U=U_ini;             % initial approximation of the solution

%
% Auxiliary arrays
%
  AUX=reshape(1:2*n_int,2,n_int);
  AUX1=AUX(:,heter_int);
  iD=repmat(AUX1,2,1); 
  jD=kron(AUX1,ones(2,1));
  WEIGHT1=WEIGHT(heter_int);

%  
% Newton's solver (the semismooth Newton method)
%
      
  it=0;               % iteration number
  while true         
      
     it=it+1;  

     % constitutive operator and its derivative
     E(:) = B*U ;   % strain at integration points
     [S,DS1]=constitutive_problem(E,heter_int,alpha,beta);
                          % solution of the constitutive problem
       
     % tangential stiffness matrix
     vD = repmat(WEIGHT1,4,1).*DS1 ;      
     D_p = sparse( iD(:),jD(:),vD(:), 2*n_int,2*n_int ) ;   
     K = K_fix+B'*D_p*B;   
 
     % vector of internal forces
     F = B'*reshape(repmat(WEIGHT,2,1).*S, 2*n_int,1) ;    

     % stopping criterion 
     criterion = norm(f-F)/norm(f);   
     crit_hist(it)=criterion;
     if  criterion < tol
        crit_hist=crit_hist(1:it);
        fprintf(' number of iteration=%d  ',it); 
        fprintf(' stopping criterion=%e  ',criterion); 
        fprintf('\n'); 
        break
     end      

     % test on number of iteration
     if  it == it_max
        fprintf('     Newton solver converges slowly: stopping criterion=%e  ',criterion)
        fprintf('\n'); 
        break
     end     
        
     % Newton's increment and next iteration
     dU = K\(f-F);  
     U= U + dU ;   
     
  end % true     
  
end % function