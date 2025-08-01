function omega=damping(U_it,dU,B,f,WEIGHT,heter_int,alpha,beta)

%--------------------------------------------------------------------------
% Computation of the damping coefficient for Newton's iteration
%--------------------------------------------------------------------------
 
  
%
% initialization 
%
  omega_min=0;    % lower bound of omega
  omega_max=2;    % upper bound of omega
  it_damp_max=10; % maximal number of damping steps

%
% finding omega by the bisection method
%  
  
  omega=(omega_min+omega_max)/2;
  it_damp=0;
  
  while it_damp<it_damp_max    
    it_damp=it_damp+1;
    U_omega = U_it + omega*dU ;
    E_omega = reshape( B*U_omega, 2,[] ) ;
    S_omega=constitutive_problem(E_omega,heter_int,alpha,beta);
    F_omega = B'*reshape(repmat(WEIGHT,2,1).*S_omega, [],1);
    decrease = (F_omega'-f')*dU;
    if decrease<0
        omega_min=omega;
    else
        omega_max=omega;
    end
    omega=(omega_min+omega_max)/2;
  end

end % function
     